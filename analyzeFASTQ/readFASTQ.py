import sys
sys.path.append('../')
from simCRN.ml_nupack import gen_complement
import config
import time
import functools
import concurrent.futures
import numpy as np
import jellyfish as jf

# # Create dictionary of barcodes
# with open('barcodes.txt', 'r') as barcode_file:
#     # Create empty dictionary
#     barcode_dict = {}
#     barcodec_dict = {}
#     # Create a counter 
#     i = 0
#     # Read lines
#     while True:
#         # Read line
#         line = barcode_file.readline()
#         # Break if the line is empty
#         if not line:
#             break
#         # Remove the newline
#         barcode = line.rstrip('\n')
#         # Add the sequence and complement to the dictionary
#         barcode_dict['barcode{}'.format(i)] = barcode
#         barcodec_dict['barcode{}c'.format(i)] = gen_complement(seq)
        
#         # Increase the counter
#         i += 1

def check_match(subseq, ref_dict, tol):
    """ check_match(subseq, ref_dict, tol)
    Checks to see if a subsequence is in the provided dictionary of sequences.
    If there is no perfect match, it calculates the Levenshtein distances 
    between the subsequence and the dictionary, normalizes the distance to a 
    score between [0,1] ([worst, perfect]) to see if there is a match with a 
    score high enough to pass the provided tolerance. The normalized score is: 
        1.0 - levenshtein_distance/ref_len

    Arguments:
        subseq (str): a string specifying the subsequence to analyze
        ref_dict (dict): a dictionary of sequences to compare the subsequence 
            against
        tol (float): a float between [0,1] to specify the accepted tolerance.
            1.0 = perfect match, 0.0 = worst match.

    Returns:
        (lev_score, ID) (tuple):
            lev_score (float): the normalized Levenshtein distance of the best 
                matching key in the dictionary
            ID (int): the value of the best matching key in the dictionary
    """
    # First check if there's a perfect match in the dictionary
    check_subseq = ref_dict.get(subseq)
    # If there is, return the value of the corresponding match
    if check_subseq != None:
        return (1.0, check_subseq)
    # If not, calculate the Levenshtein distance between the subsequence and the 
    # keys in the dictionary
    else:
        # Store the score in a tuple (lev_score, index, key)
        score = (0.0, -1)
        for key in ref_dict.keys():
            # Calculate the Levenshtein distance
            lev_dist = jf.levenshtein_distance(subseq, key)
            # Normalize the Levenshtein distance to [0, 1]
            lev_score = 1.0 - lev_dist / len(key)
            # If the score is greater than the tolerance, return the index of 
            # the corresponding key
            if lev_score >= tol:
                return (lev_score, ref_dict.get(key))
            # If not, check if the score is better than the last saved score.
            else:
                # Save the score and index if the score is better
                if lev_score > score[0]:
                    score = (lev_score, ref_dict.get(key))
                else:
                    pass
        # Getting to this point means a match wasn't found, so return the score
        # tuple
        return score

def pick_best_score(start, start_del, start_ins, seq_len, sequence, ref_len, \
    ref_sticky_dict, tol):
    """
    pick_best_scorestart, start_del, start_ins, seq_len, sequence, ref_len, 
        ref_sticky_dict, tol)

    This function uses 3 start indices to check if subsequences with those start
    indices have a match in the reference sticky dictionary, within the desired
    tolerance. 

    Arguments:
        start (int): The starting index of a subsequence to check.
        start_del (int): The starting index (shifted left, -1) of a subsequence 
            to check.
        start_ins (int): The starting index (shifted right, +1) of a subsequence 
            to check.
        seq_len (int): The length of the sequence being analyzed.
        sequence (str): The sequence being analyzed.
        ref_len (int): The length of the sticky reference sequences (28 for 
            targets, 33 for barcodes).
        ref_sticky_dict (dict): The dictionary of sticky reference sequences.
            Keys = sticky end + barcode/target + sticky end, values = ID.
        tol (float): A float between [0,1] specifying how well a subsequence
            should match the reference sequence to be considered a match. This
            float is equivalent to a normalized Levenshtein distance. Default is
            0.9.
    Returns:
        score (tuple): A score tuple of (lev_score, ID) where
            lev_score (float): the normalized Levenshtein distance of the best 
                matching key in the dictionary
            ID (int): the value of the best matching key in the dictionary
        subseq (str): A string representing the subsequence that best matched a
            reference sequence in the provided dictionary.
        start (int): The start index of the best matching subsequence.
        end (int): The end index of the best matching subsequence.
    """
    # Calculate the subsequence end indices
    end = start + ref_len
    end_del = end - 1
    end_ins = end + 1
    # If the end indices are out of bounds, cap them at the end of the \
    # sequence
    if end_del > seq_len:
        end = seq_len
        end_del = seq_len
        end_ins = seq_len
    elif end_del == seq_len:
        end = seq_len
        end_ins = seq_len
    elif end == seq_len:
        end_ins = seq_len
    else:
        pass
    # Check matches for the subsequence, and the subsequence shifted over 
    # by 1
    subseq = sequence[start:end]
    subseq_del = sequence[start_del:end_del]
    subseq_ins = sequence[start_ins:end_ins]
    # Check if the sequence matches the sticky target
    score = check_match(subseq, ref_sticky_dict, tol)
    score_del = check_match(subseq_del, ref_sticky_dict, tol)
    score_ins = check_match(subseq_ins, ref_sticky_dict, tol)
    # Pick the best match
    # If the +1 (ins) subsequence is best, save its score, subsequence, and
    # start and end indices.
    if score_ins[0] > score[0] and score_ins[0] > score_del[0]:
        score = score_ins
        subseq = subseq_ins
        start = start_ins
        end = end_ins
    # If the -1 (del) subsequence is best, save its score, subsequence, and
    # start and end indices.
    elif score_del[0] > score[0] and score_del[0] > score_ins[0]:
        score = score_del
        subseq = subseq_del
        start = start_del
        end = end_del
    # If the original subsequence is best, save its start index.
    else:
        start = start
    return(score, subseq, start, end)

def list_check_match(sequence, seq_ID, seq_list, fastq_file, barcode_dict, \
    target_dict, lengths, count_vec, error_list, tol):
    """
    slide_check_match(sequence, seq_ID, fastq_file, sticky_barcode_dict, \
    sticky_target_dict, count_vec, error_list, tol)

    Takes in a sequence, and checks subsequences to see if they match 
    sticky end + target + sticky end or sticky end + barcode + sticky end. Once
    a subsequence is analyzed for a match, the function slides to the next 
    region and repeats the analysis. Target counts are stored in the count 
    array and errors are stored in the provided error list. The barcode ID of 
    the sequence is also determined, an a 0/1 error status is set.

    Arguments:
        sequence (str): a string representing a DNA sequence
        seq_ID (str): the sequence ID, starts with @
        seq_list (list): a list of the subsequences split by the sticky end 
            sequence 
        fastq_file (str): the file name of the fastq file that is being analyzed
        barcode_dict (dict): a dictionary containing barcode sequences (keys) 
            and barcode IDs (values)
        target_dict (dict): a dictionary containing target sequences (keys) and
            target IDs (values)
        lengths (tuple): (sticky_len, barcode_len, target_len, \
            sticky_barcode_len, sticky_target_len) 
            A tuple of the significant sequence lengths (int).
            sticky_len: length of sticky end
            barcode_len: length of the barcode sequences
            target_len: length of the target sequences
            sticky_barcode_len: length of 2 * sticky_end len + barcode_len
            sticky_target_len: length of 2 * sticky_end len + target_len
        count_vec (numpy array): a nx1 numpy array to store target sequence 
            counts, where n is the number of target sequences
        error_list (list): a list to store error messages in
        tol (float): A float between [0,1] specifying how well a subsequence
            should match the reference sequence to be considered a match. This
            float is equivalent to a normalized Levenshtein distance. Default is
            0.9.
    Returns:
        (barcode_ID, has_error) (tuple): 
            barcode_ID (int): An integer representing the barcode that the 
                sequence was tagged with. barcode_ID = -2 if a barcode 
                was not found, -1 if a barcode mismatch was found, and 0+ if a
                barcode was found. 
            has_error (bool): True if an error is encountered during analysis,
                False otherwise.
    """
    # Unpack the lengths tuple
    barcode_len = lengths[1]
    target_len = lengths[2]
    # Set barcode ID to -2 (unidentified)
    barcode_ID = -2
    # Set error flag to False
    has_error = False
    # Go through the subsequences and identify them
    for subseq in seq_list:
        # Check if the subsequence is in the target dictionary
        if len(subseq) == target_len:
            target_score = check_match(subseq, target_dict, tol)
            # If there's a match, increase the count in the array
            if target_score[0] >= tol:
                count_vec[target_score[1]] += 1
            # If there's no match, create a string of the fastq file name, 
            # sequence ID, error type, subsequence, and best match 
            # score to add to the error list
            else:
                error_list.append('{}\t{}\tno_match_subseq\t{}\t{}\t{}\n'\
                    .format(fastq_file, seq_ID, subseq, target_score[0], \
                        target_score[1]))
                # Set error status to True
                has_error = True
        # Check if the subsequence is in the barcode dictionary
        elif len(subseq) == barcode_len:
            barcode_score = check_match(subseq, barcode_dict, tol)
            # If there's a match, set the barcode ID
            if barcode_ID == -2:
                barcode_ID = barcode_score[1]
            elif barcode_score[0] >= tol and barcode_ID > -1 and \
                barcode_score[1] != barcode_ID:
                # Record error
                error_list.append('{}\t{}\tbarcode_mismatch\t{}\n'\
                    .format(fastq_file, seq_ID, sequence))
                # Set barcode ID to -1 (mismatch)
                barcode_ID = -1
                # Set error status to True
                has_error = True
            else:
                pass
        else:
            pass
    return (barcode_ID, has_error)

def slide_check_match(sequence, seq_ID, fastq_file, sticky_barcode_dict, \
    sticky_target_dict, lengths, count_vec, error_list, tol):
    """
    slide_check_match(sequence, seq_ID, fastq_file, sticky_barcode_dict, \
    sticky_target_dict, count_vec, error_list, tol)

    Takes in a sequence, and checks subsequences to see if they match 
    sticky end + target + sticky end or sticky end + barcode + sticky end. Once
    a subsequence is analyzed for a match, the function slides to the next 
    region and repeats the analysis. Target counts are stored in the count 
    array and errors are stored in the provided error list. The barcode ID of 
    the sequence is also determined, an a 0/1 error status is set.

    Arguments:
        sequence (str): a string representing a DNA sequence
        seq_ID (str): the sequence ID, starts with @
        fastq_file (str): the file name of the fastq file that is being analyzed
        sticky_barcode_dict (dict): a dictionary containing barcode sequences 
            with sticky ends at both ends (keys) and barcode IDs (values)
        sticky_target_dict (dict): a dictionary containing target sequences
            with sticky ends at both ends (keys) and target IDs (values)
        lengths (tuple): (sticky_len, barcode_len, target_len, \
            sticky_barcode_len, sticky_target_len) 
            A tuple of the significant sequence lengths (int).
            sticky_len: length of sticky end
            barcode_len: length of the barcode sequences
            target_len: length of the target sequences
            sticky_barcode_len: length of 2 * sticky_end len + barcode_len
            sticky_target_len: length of 2 * sticky_end len + target_len
        count_vec (numpy array): a nx1 numpy array to store target sequence 
            counts, where n is the number of target sequences
        error_list (list): a list to store error messages in
        tol (float): A float between [0,1] specifying how well a subsequence
            should match the reference sequence to be considered a match. This
            float is equivalent to a normalized Levenshtein distance. Default is
            0.9.
    Returns:
        (barcode_ID, has_error) (tuple): 
            barcode_ID (int): An integer representing the barcode that the 
                sequence was tagged with. barcode_ID = -2 if a barcode 
                was not found, -1 if a barcode mismatch was found, and 0+ if a
                barcode was found. 
            has_error (bool): True if an error is encountered during analysis,
                False otherwise.
    """
    # Unpack the lengths tuple
    sticky_len, barcode_len, target_len, sticky_barcode_len, \
    sticky_target_len = lengths
    # Save the length of the sequence
    seq_len = len(sequence)
    # Set starting index to 0
    start = 0
    # Set barcode ID to -2 (unidentified)
    barcode_ID = -2
    # Set subsequence identification flag to True
    identified = True
    # Set start point of unidentified subsequence to 0
    unID_start = 0
    # Set sequence error flag to False
    has_error = False
    # Make sure the start point is within the sequence
    while start + sticky_len < seq_len:
        # Set the start points for deletion and insertion subsequences
        start_del = start - 1
        start_ins = start + 1
        # If the deletion start point is out of bounds, set it to 0
        if start_del < 0:
            start_del = start
        else:
            pass
        # Pick the best of the unshifted/shifted subsequence target matches
        target_score, target_subseq, target_start, target_end = \
        pick_best_score(start, start_del, start_ins, seq_len, sequence, \
            sticky_target_len, sticky_target_dict, tol)
        # If there's a match, add it to the count
        if target_score[0] >= tol:
            # Update the count
            count_vec[target_score[1]] += 1
            # If a preceding subsequence was not identified, save the
            # unidentified subsequence as an error
            if identified == False:
                # Set the start point to the end point of the unidentified 
                # sequence and record the unidentified sequence in the error 
                # list
                unID_subseq = sequence[unID_start:target_start]
                error_list.append('{}\t{}\tunID_subseq\t{}\n'\
                .format(fastq_file, seq_ID, unID_subseq))
                # Set status to identified 
                identified = True
                # Set error status to True
                has_error = True
            else:
                pass
            # Set the next start index to the target end index - sticky length
            start = target_end - sticky_len
        else:
            # Pick the best of the unshifted/shifted subsequence target matches
            barcode_score, barcode_subseq, barcode_start, barcode_end = \
            pick_best_score(start, start_del, start_ins, seq_len, sequence, \
                sticky_barcode_len, sticky_barcode_dict, tol)
            # Check to see if the target or barcode match is better
            # If the target match is better
            if target_score[0] >= barcode_score[0]:
                # If the status was previously identified, save the start 
                # point and set the status to unidentified
                if identified == True:
                    unID_start = start
                    identified = False
                # If the status was unidentified, and we've hit the end of the 
                # sequence, save the end subsequence as an unidentified error
                elif identified == False and target_end == seq_len:
                    unID_subseq = sequence[unID_start:target_end]
                    error_list.append('{}\t{}\tunID_subseq\t{}\n'\
                    .format(fastq_file, seq_ID, unID_subseq))
                    # Set error status to True
                    has_error = True
                    break
                else:
                    pass
                # Set start to check next 3 indices
                start = start + 3
            else:
                if barcode_score[0] >= tol: 
                    # If a preceding subsequence was not identified, save the
                    # unidentified subsequence
                    if identified == False:
                        # Set the start point to the end point of the 
                        # unidentified sequence and record the unidentified 
                        # sequence in the error list
                        unID_subseq = sequence[unID_start:barcode_start]
                        error_list.append('{}\t{}\tunID_subseq\t{}\n'\
                        .format(fastq_file, seq_ID, unID_subseq))
                        # Set status to identified 
                        identified = True
                        # Set error status to True
                        has_error = True
                    else:
                        pass
                    # If the barcode was previously not found, save the barcode
                    # ID
                    if barcode_ID == -2:
                        barcode_ID = barcode_score[1]
                    # If the barcode has already been found, and a previous
                    # barcode error has not been found, record the error
                    elif barcode_score[0] >= tol and barcode_ID > -1 and \
                        barcode_score[1] != barcode_ID:
                        # Record error
                        error_list.append('{}\t{}\tbarcode_mismatch\t{}\n'\
                            .format(fastq_file, seq_ID, sequence))
                        # Set barcode ID flag to -1 (mismatch)
                        barcode_ID = -1
                        # Set error status to True
                        has_error = True
                    else:
                        pass
                    # Set the start to the barcode end - sticky length
                    start = barcode_end - sticky_len
                else:
                    # If the status was previously identified, save the start 
                    # point and set the status to unID
                    if identified == True:
                        unID_start = start
                        identified = False
                    elif identified == False and barcode_end == seq_len:
                        unID_subseq = sequence[unID_start:barcode_end]
                        error_list.append('{}\t{}\tunID_subseq\t{}\n'\
                        .format(fastq_file, seq_ID, unID_subseq))
                        # Set error status to True
                        has_error = True
                        break
                    else:
                        pass
                    # Set start to check next 3 indices
                    start = start + 3
    return barcode_ID, has_error

def analyze_sequence(sequence, seq_ID, fastq_file, sticky_end, lengths, 
    barcode_dict, target_dict, sticky_barcode_dict, sticky_target_dict, \
    count_vec, error_list, tol=0.9):
    """ 
    analyze_sequence(sequence, seq_ID, fastq_file, count_vec, error_list, tol) 
    Analyzes a given sequence for barcode and target sequences. It identifies 
    which barcode the sequence is tagged with, and adds the target counts to the
    appropriate index in a provided numpy array. For sequences that don't 
    perfectly match the reference target sequences, it calculates a match score 
    using a normalized Levenshtein distance. The match is accepted if the match
    score is >= the provided tolerance value. If errors arise, it adds a string 
    documenting the error and sequence to the provided error list.

    Arguments:
        sequence (str): a string representing a DNA sequence
        seq_ID (str): the sequence ID, starts with @
        fastq_file (str): the file name of the fastq file that is being analyzed
        sticky_end (str): a string representing a sticky end sequence
        barcode_dict (dict): a dictionary containing barcode sequences (keys) 
            and barcode IDs (values)
        target_dict (dict): a dictionary containing target sequences (keys) and
            target IDs (values)
        sticky_barcode_dict (dict): a dictionary containing barcode sequences 
            with sticky ends at both ends (keys) and barcode IDs (values)
        sticky_target_dict (dict): a dictionary containing target sequences
            with sticky ends at both ends (keys) and target IDs (values)
        lengths (tuple): (sticky_len, barcode_len, target_len, 
            sticky_barcode_len, sticky_target_len) 
            A tuple of the significant sequence lengths (int).
            sticky_len: length of sticky end
            barcode_len: length of the barcode sequences
            target_len: length of the target sequences
            sticky_barcode_len: length of 2 * sticky_end len + barcode_len
            sticky_target_len: length of 2 * sticky_end len + target_len
        count_vec (numpy array): a nx1 numpy array to store target sequence 
            counts, where n is the number of target sequences
        error_list (list): a list to store error messages in
        tol (float): A float between [0,1] specifying how well a subsequence
            should match the reference sequence to be considered a match. This
            float is equivalent to a normalized Levenshtein distance. Default is
            0.9.

    Returns:
        (barcode_ID, has_error) (tuple): 
            barcode_ID (int): An integer representing the barcode that the 
                sequence was tagged with. barcode_ID = -2 if a barcode 
                is not found, -1 if a barcode mismatch was found, and 0+ if a
                barcode was found. 
            has_error (bool): True if an error is encountered during analysis,
                False otherwise.
    Errors:
        no target/barcode match: a subsequence (correct length) not identifiable
            as barcode or target
            "fastq_file [tab] seq_ID [tab] subseq [tab] lev_score [tab] ID"
        no match: a subsequence (correct or incorrect length) not identifiable
            as barcode or target
            "fastq_file [tab] seq_ID [tab] subseq"
        no barcode: no barcode ID successfully determined for the sequence
            "fastq_file [tab] seq_ID [tab] sequence"

    """
    # Assume we will be able to chop the sequence into a list of subsequences
    check_list = True
    # Set the count vector to 0
    count_vec[:] = np.zeros(len(count_vec))
    # Unpack the lengths tuple
    sticky_len, barcode_len, target_len, sticky_barcode_len, \
    sticky_target_len = lengths
    # Check for sticky ends at the beginning and end of the sequence
    head = sequence[:sticky_len]
    tail = sequence[-sticky_len:]
    # If the ends match the sticky ends, remove them
    if head == sticky_end and tail == sticky_end:
        stripped_seq = sequence[sticky_len:-sticky_len]
        # Try splitting the sequence up into 
        seq_list = stripped_seq.split(sticky_end)
        # Set barcode ID flag to -2 (not identified):
        barcode_ID = -2
        # Check to make sure the subsequences are all the correct sizes
        for subseq in seq_list:
            if len(subseq) == target_len or len(subseq) == barcode_len:
                pass
            # If a subsequence isn't the right size, set the check_list variable
            # to False
            else:
                check_list = False
    # Otherwise set check_list to false to parse through the sequence instead of
    # analyzing subseequences 
    else:
        check_list = False
    
    # If check_list is True, identify each subsequence in the list
    if check_list:
        barcode_ID, has_error = list_check_match(sequence, seq_ID, seq_list, \
            fastq_file, barcode_dict, target_dict, lengths, count_vec, \
            error_list, tol)
    # If check list is false, use a sliding check
    else:
        barcode_ID, has_error = slide_check_match(sequence, seq_ID, fastq_file, \
            sticky_barcode_dict, sticky_target_dict, lengths, count_vec, \
            error_list, tol)
        # sliding check

    # If no barcode was found, 
    if barcode_ID == -2:
        error_list.append('{}\t{}\tno_barcode\t{}\n'.format(fastq_file, seq_ID, \
            sequence))
        # Set error status to true
        has_error = True
    else:
        pass  
    return (barcode_ID, has_error)

def test_pick_info_files(barcode_file, target_file):
    """
    pick_info_files(barcode_file, target_file)

    Used to specify the files containing the barcode
    """
    def save_file_decorator(func):
        @functools.wraps(func)
        def wrapper(fastq_file):
            return func(fastq_file, barcode_file, target_file)
        return wrapper
    return save_file_decorator

@test_pick_info_files('barcode teehee', 'target uwu')
def test_func(fastq_file, barcode_file, target_file):
    """
    my docstring oolala
    """
    return fastq_file

# test_func('file owo')

class EmptyLineErr(Exception):
    def __init__(self, arg):
        self.strerror = arg
        self.args = {arg}

def make_dictionary(file_name, sticky_end, is_sticky=True):
    """
    make_dictionary(file_name, sticky_end)

    Takes in the name of a file that contains sequences (each on a different 
    line), and a sticky end sequence. Creates 2 dictionaries of the sequences.
    One dictionary has the sequences and their complement as keys, and the
    sequence ID as values. One dictionary has the sequences and complements 
    sandwiched by the sticky end sequence as keys, and the sequence ID as values.
    This function also identifies the length of the sequences, as well as the 
    length of the sequence sandwiched by sticky ends.

    Arguments:
        file_name (str): A string representing the file of sequences to create
            dictionaries from. Each sequence should be in its own line, with no
            empty lines.
        sticky_end (str): a string representing a sticky end sequence
        is_sticky (bool): True if the sequences in the sequence file are 
            sandwiched by sticky ends. False if not. Defaults to True.

    Returns: 
        (dictionary, sticky_dict, seq_len, sticky_seq_len, num_seq) (tuple):
            dictionary (dict): A dictionary with sequences and complements as 
            keys and the sequence ID as values.
            sticky_dict (dict): A dictionary with sequences and complements 
                sandwiched by the sticky end sequence as keys and the sequence 
                ID as values.
            seq_len (int): The length of the sequences/complements in dictionary.
            sticky_seq_len (int): The length of the sequences/complements in 
                sticky_dict.
            num_seq (int): The number of sequences in the file.

    Errors: 
        ValueError: raises an exception if an empty line is found in the
            sequence file. Exits script.
        AssertionError: raises an exception if a sequence with a different 
            length is found. Exits script.
    """
    # Create empty dictionaries
    dictionary = {}
    sticky_dict = {}
    sticky_len = len(sticky_end)
    # Open the file with the sequences
    with open(file_name, 'r') as file:
        # # Read all the sequences in the file
        lines = file.readlines()
        # Save the sequence count
        num_seq = len(lines)
        # Save the sequence length of the first line as a reference
        # If the file has sequences sandwiched by sticky ends
        if is_sticky:
            ref_len = len(lines[0])-1-2*sticky_len
        # If the file has sequences without sticky ends
        else:
            ref_len = len(lines[0])-1
        try:
            for i in range(len(lines)):
                # Read in the whole sequence
                whole_sequence = lines[i].rstrip('\n')
                # If the line is empty, skip the line
                if whole_sequence != '':
                    # If the sequences have sticky ends, cut the sticky ends off
                    # to get just the sequence of interest
                    if is_sticky:
                        sequence = whole_sequence[sticky_len:-sticky_len]
                        sticky_sequence = whole_sequence
                    # If the sequences don't have sticky ends, add the sticky
                    # ends to get the sticky end sequence
                    else:
                        sequence = whole_sequence
                        sticky_sequence = sticky_end + whole_sequence + sticky_end
                    # Figure out if the sequence is the right length
                    assert len(sequence) == ref_len,\
                    'AssertionError: Sequence {} is a different length than '\
                    .format(i) + 'expected in sequence file {}.'.format(file_name) + \
                    '\nExpected length {}, got {}.'.format(ref_len, len(sequence))
                    # Generate complements
                    sticky_comp = gen_complement(sticky_sequence)
                    comp = gen_complement(sequence)
                    # Add the sequences to the dictionaries
                    sticky_dict[sticky_sequence] = i
                    sticky_dict[sticky_comp] = i
                    dictionary[sequence] = i
                    dictionary[comp] = i
                    # Increase the sequence count
                else:
                    raise ValueError('ValueError: Found an empty line in {}.'\
                        .format(file_name) + 
                        '\nMake sure sequence files do not contain empty lines.')
            # Figure out the sequence length and the sticky sequence length
            seq_len = len(list(dictionary.keys())[0])
            sticky_seq_len = int(seq_len + 2*len(sticky_end))
            return (dictionary, sticky_dict, seq_len, sticky_seq_len, num_seq)

        except (ValueError, AssertionError) as msg:
            print(msg)
            # Stop the script execution
            sys.exit(1)

def write_errors(error_file, error_list):
    '''
    write_errors(error_file, error_list)
    
    Appends the errors in the list of errors to an error log file.

    Arguments:
        error_file (str): A string representing the name of the file to append
            the errors to
        error_list (list): A list of error messages (str)

    Returns:
        nothing
    '''
    with open(error_file, 'a') as file:
        for error in error_list:
            file.write(error)

def pick_info_files(barcode_file, target_file, sticky_end, tol=0.9):
    """
    pick_info_files(barcode_file, target_file)

    Used to specify the files containing the barcode
    """
    def save_file_decorator(func):
        @functools.wraps(func)
        def wrapper(fastq_file):
            return func(fastq_file, barcode_file, target_file, sticky_end, \
                tol=tol)
        return wrapper
    return save_file_decorator

def specify_tol(tol=0.9):
    """
    specify_tol(tol=0.9)

    Used to specify the tolerance that will be used for sequence matching.
    Use as a decorator to decorate analyze_fastq(fastq_file, tol=0.9).

    Arguments:
        tol (float): a float between [0,1] inclusive, specifying the tolerance
            to use for sequence matching. 1 is a perfect match.

    Returns:
        A decorated function with tol set to the specified tolerance.
    """
    def save_file_decorator(func):
        @functools.wraps(func)
        def wrapper(fastq_file):
            return func(fastq_file, tol=tol)
        return wrapper
    return save_file_decorator



def analyze_barcodes_file(fastq_file):
    """
    This function reads in a fastq file, calculates and returns the number of 
    barcode instances in it.
    Arguments:
        fastq_file (str): A string representing the name of the fastq file to
        read.
            ex: test.fastq
    Returns:
        num_barcodes (int): The number of times the barcode showed up in the 
        file
    """
    with open(fastq_file, 'r') as file:
        # Save the number of barcodes
        num_barcodes = 0
        # for i in range(10):
        while True:
            # Read the sequence ID
            seq_ID = file.readline().rstrip('\n')
            # If we've reached the end of the file, break the loop
            if not seq_ID:
                break
            # Read in the sequence and strip the new line
            seq = file.readline().rstrip('\n')
            file.readline()
            Q_score = file.readline().rstrip('\n')

            # Figure out how many barcodes are in the sequence
            # Update the count
            num_barcodes += int((len(seq)-4)/29)        
    return num_barcodes


if __name__ == "__main__":
    sticky_end = 'TGCA'
    # Test values
    file_list = ['tiny.fastq', 'approx.fastq', 'tiny_error.fastq']
    barcode_file = 'test_barcodes.txt'
    target_file = 'test_targets.txt'
    error_file = 'test_errors.txt'
    count_file = 'test_count.txt'

    # Make the barcode and sticky barcode dictionaries
    barcode_dict, sticky_barcode_dict, barcode_len, sticky_barcode_len, \
    num_barcodes = make_dictionary(barcode_file, sticky_end)
    # Make the target and sticky target dictionaries
    target_dict, sticky_target_dict, target_len, sticky_target_len, num_targets\
     = make_dictionary(target_file, sticky_end, is_sticky=False)
    # Determine the sticky end length
    sticky_len = len(sticky_end)
    # Make the lengths tuple
    lengths = (sticky_len, barcode_len, target_len, sticky_barcode_len, \
        sticky_target_len) 
    # Create an array to store counts in
    total_count_array = np.zeros((num_barcodes, num_targets))
    # Create a variable to store the number of sequences analyzed
    total_analyzed_seq = 0
    # Create a variable to store the number of sequences with errors
    total_seq_error = 0
    # Create a variable to store the total number of errors found
    total_error = 0

    
    @specify_tol()            
    def analyze_fastq(fastq_file, tol=0.9):
        print('Tol = {}'.format(tol))
        count_array = np.zeros((num_barcodes, num_targets))
        # Create a vector to store target counts in 
        count_vec = np.zeros(num_targets)
        # Open the target file and analyze the sequences
        with open(fastq_file, 'r') as file:
            # Count the number of sequences
            num_seq = 0
            num_errors = 0
            # Create a list to store errors in
            error_list = []
            while True:
                # Read the sequence ID
                seq_ID = file.readline().rstrip('\n')
                # If we've reached the end of the file, break the loop
                if not seq_ID:
                    break
                # Read in the sequence and strip the new line
                seq = file.readline().rstrip('\n')
                # Update the sequence count
                num_seq += 1
                # Read in the + line
                file.readline()
                # Read in the Q_score
                Q_score = file.readline().rstrip('\n')
                # Analyze the sequence
                barcode_ID, has_error = analyze_sequence(seq, seq_ID, \
                    fastq_file, sticky_end, lengths, barcode_dict, target_dict, \
                    sticky_barcode_dict, sticky_target_dict, count_vec, error_list, \
                    tol=tol)
                # Update the number of sequences with errors
                num_errors += has_error
                # Add the target counts to the count array
                if barcode_ID >= 0:
                    count_array[barcode_ID,:] += count_vec
                else:
                    pass
        return (count_array, num_seq, error_list, num_errors)

    # Parallelize the file analysis process
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(analyze_fastq, file_list)
        for result in results:
            count_arr, num_seq, error_list, num_errors = result
            # Add the counts from each file to the total array
            total_count_array += count_arr
            # Write the errors from each file to the error file
            write_errors(error_file, error_list)
            # Update the number of sequences analyzed
            total_analyzed_seq += num_seq
            # Update the number of sequences with errors
            total_seq_error += num_errors
            # Update the number of errors found
            total_error += len(error_list)

    print('total count array')
    print(total_count_array)
    print('total analyzed sequences: ', total_analyzed_seq)
    print('total sequences with errors: ', total_seq_error)
    print('total errors: ', total_error)
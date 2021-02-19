import config
import time
import functools
import concurrent.futures
import numpy as np
import jellyfish as jf

sticky_end = config.sticky_end
barcode = config.barcode
barcodec = config.barcodec
barcode1 = config.barcode1
target = config.target
targetc = config.targetc
target1 = config.target1

barcode_dict = config.barcode_dict
target_dict = config.target_dict
sticky_barcode_dict = config.sticky_barcode_dict
sticky_target_dict = config.sticky_target_dict

sticky_len = config.sticky_len
target_len = config.target_len
barcode_len = config.barcode_len

file_name = config.file_name

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

def list_check_match(sequence, seq_id, seq_list, file_name, barcode_dict, \
    target_dict, count_array, error_list, tol):
    """
    slide_check_match(sequence, seq_id, file_name, sticky_barcode_dict, \
    sticky_target_dict, count_array, error_list, tol)

    Takes in a sequence, and checks subsequences to see if they match 
    sticky end + target + sticky end or sticky end + barcode + sticky end. Once
    a subsequence is analyzed for a match, the function slides to the next 
    region and repeats the analysis. Target counts are stored in the count 
    array and errors are stored in the provided error list. The barcode ID of 
    the sequence is also determined, an a 0/1 error status is set.

    Arguments:
        sequence (str): a string representing a DNA sequence
        seq_id (str): the sequence ID, starts with @
        seq_list (list): a list of the subsequences split by the sticky end 
            sequence 
        file_name (str): the file name of the fastq file that is being analyzed
        barcode_dict (dict): a dictionary containing barcode sequences (keys) 
            and barcode IDs (values)
        target_dict (dict): a dictionary containing target sequences (keys) and
            target IDs (values)
        count_array (numpy array): a nx1 numpy array to store target sequence 
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
            has_error (int): A 1 if an error is encountered during analysis,
                0 otherwise.
    """
    # Set barcode ID to -2 (unidentified)
    barcode_ID = -2
    # Set error flag to 0 (false)
    has_error = 0
    # Go through the subsequences and identify them
    for subseq in seq_list:
        # Check if the subsequence is in the target dictionary
        if len(subseq) == target_len:
            target_score = check_match(subseq, target_dict, tol)
            # If there's a match, increase the count in the array
            if target_score[0] >= tol:
                count_array[target_score[1]] += 1
            # If there's no match, create a string of the fastq file name, 
            # sequence ID, error type, subsequence, and best match 
            # score to add to the error list
            else:
                error_list.append('{}\t{}\tno_match_subseq\t{}\t{}\t{}'\
                    .format(file_name, seq_id, subseq, target_score[0], \
                        target_score[1]))
                # Set error status to 1
                has_error = 1
        # Check if the subsequence is in the barcode dictionary
        elif len(subseq) == barcode_len:
            barcode_score = check_match(subseq, barcode_dict, tol)
            # If there's a match, set the barcode ID
            if barcode_ID == -2:
                barcode_ID = barcode_score[1]
            elif barcode_score[0] >= tol and barcode_ID > -1 and \
                barcode_score[1] != barcode_ID:
                # Record error
                error_list.append('{}\t{}\tbarcode_mismatch\t{}'\
                    .format(file_name, seq_id, sequence))
                # Set barcode ID to -1 (mismatch)
                barcode_ID = -1
                # Set error status to 1
                has_error = 1
            else:
                pass
        else:
            pass
    return (barcode_ID, has_error)

def slide_check_match(sequence, seq_id, file_name, sticky_barcode_dict, \
    sticky_target_dict, count_array, error_list, tol):
    """
    slide_check_match(sequence, seq_id, file_name, sticky_barcode_dict, \
    sticky_target_dict, count_array, error_list, tol)

    Takes in a sequence, and checks subsequences to see if they match 
    sticky end + target + sticky end or sticky end + barcode + sticky end. Once
    a subsequence is analyzed for a match, the function slides to the next 
    region and repeats the analysis. Target counts are stored in the count 
    array and errors are stored in the provided error list. The barcode ID of 
    the sequence is also determined, an a 0/1 error status is set.

    Arguments:
        sequence (str): a string representing a DNA sequence
        seq_id (str): the sequence ID, starts with @
        file_name (str): the file name of the fastq file that is being analyzed
        sticky_barcode_dict (dict): a dictionary containing barcode sequences 
            with sticky ends at both ends (keys) and barcode IDs (values)
        sticky_target_dict (dict): a dictionary containing target sequences
            with sticky ends at both ends (keys) and target IDs (values)
        count_array (numpy array): a nx1 numpy array to store target sequence 
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
            has_error (int): A 1 if an error is encountered during analysis,
                0 otherwise.
    """
    # Save the lengths of the sticky target/barcode subsequences
    sticky_target_len = int(sticky_len * 2 + target_len)
    sticky_barcode_len = int(sticky_len * 2 + barcode_len)
    seq_len = len(sequence)
    # Set starting index to 0
    start = 0
    # Set barcode ID to -2 (unidentified)
    barcode_ID = -2
    # Set subsequence identification flag to 1 (true)
    identified = 1
    # Set start point of unidentified subsequence to 0
    unID_start = 0
    # Set sequence error flag to 0
    has_error = 0
    # Make sure the start point is within the sequence
    while start + 4 < seq_len:
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
            count_array[target_score[1]] += 1
            # If a preceding subsequence was not identified, save the
            # unidentified subsequence as an error
            if identified == 0:
                # Set the start point to the end point of the unidentified 
                # sequence and record the unidentified sequence in the error 
                # list
                unID_subseq = sequence[unID_start:target_start]
                error_list.append('{}\t{}\tunID_subseq\t{}\t{}'\
                .format(file_name, seq_id, unID_subseq))
                # Set status to identified 
                identified = 1
                # Set error status to 1
                has_error = 1
            else:
                pass
            # Set the next start index to the target end index - 4
            start = target_end - 4
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
                if identified == 1:
                    unID_start = start
                    identified = 0
                # If the status was unidentified, and we've hit the end of the 
                # sequence, save the end subsequence as an unidentified error
                elif identified == 0 and target_end == seq_len:
                    unID_subseq = sequence[unID_start:target_end]
                    error_list.append('{}\t{}\tunID_subseq\t{}\t{}'\
                    .format(file_name, seq_id, unID_subseq))
                    # Set error status to 1
                    has_error = 1
                    break
                else:
                    pass
                # Set start to check next 3 indices
                start = start + 3
            else:
                if barcode_score[0] >= tol: 
                    # If a preceding subsequence was not identified, save the
                    # unidentified subsequence
                    if identified == 0:
                        # Set the start point to the end point of the 
                        # unidentified sequence and record the unidentified 
                        # sequence in the error list
                        unID_subseq = sequence[unID_start:barcode_start]
                        error_list.append('{}\t{}\tunID_subseq\t{}\t{}'\
                        .format(file_name, seq_id, unID_subseq))
                        # Set status to identified 
                        identified = 1
                        # Set error status to 1
                        has_error = 1
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
                        error_list.append('{}\t{}\tbarcode_mismatch\t{}'\
                            .format(file_name, seq_id, sequence))
                        # Set barcode ID flag to -1 (mismatch)
                        barcode_ID = -1
                        # Set error status to 1
                        has_error = 1
                    else:
                        pass
                    # Set the start to the barcode end
                    start = barcode_end - 4
                else:
                    # If the status was previously identified, save the start 
                    # point and set the status to unID
                    if identified == 1:
                        unID_start = start
                        identified = 0
                    elif identified == 0 and barcode_end == seq_len:
                        unID_subseq = sequence[unID_start:barcode_end]
                        error_list.append('{}\t{}\tunID_subseq\t{}\t{}'\
                        .format(file_name, seq_id, unID_subseq))
                        # Set error status to 1
                        has_error = 1
                        break
                    else:
                        pass
                    # Set start to check next 3 indices
                    start = start + 3
    return barcode_ID, has_error

def analyze_sequence(sequence, seq_id, count_array, error_list, tol=0.9):
    """ 
    analyze_sequence(sequence, seq_id, count_array, error_list, tol) 
    Analyzes a given sequence for barcode and target sequences. It identifies 
    which barcode the sequence is tagged with, and adds the target counts to the
    appropriate row in a provided numpy array. For sequences that don't 
    perfectly match the reference target sequences, it calculates a match score 
    using a normalized Levenshtein distance. The match is accepted if it is 
    within the provided tolerance value. If errors arise, it adds a string 
    documenting the error and sequence to the provided error list.

    Arguments:
        sequence (str): a string representing a DNA sequence
        seq_id (str): the sequence ID, starts with @
        count_array (numpy array): a nx1 numpy array to store target sequence 
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
            has_error (int): A 1 if an error is encountered during analysis,
                0 otherwise.
    Errors:
        no target/barcode match: a subsequence (correct length) not identifiable
            as barcode or target
            "file_name [tab] seq_id [tab] subseq [tab] lev_score [tab] ID"
        no match: a subsequence (correct or incorrect length) not identifiable
            as barcode or target
            "file_name [tab] seq_id [tab] subseq"
        no barcode: no barcode ID successfully determined for the sequence
            "file_name [tab] seq_id [tab] sequence"

    """
    # Assume we will be able to chop the sequence into a list of subsequences
    check_list = True
    # Check for sticky ends at the beginning and end of the sequence
    head = sequence[:4]
    tail = sequence[-4:]
    # If the ends match the sticky ends, remove them
    if head == sticky_end and tail == sticky_end:
        stripped_seq = sequence[4:-4]
        # Try splitting the sequence up into 
        seq_list = stripped_seq.split(sticky_end)
        # Set barcode ID flag to -2 (not identified):
        barcode_ID = -2
        # Check to make sure the subsequences are all the correct sizes
        for subseq in seq_list:
            if len(subseq) == target_len or len(subseq) == barcode_len:
                pass
            # If a subsequence isn't the right size, set the check_list variable to 
            # False
            else:
                check_list = False
    # Otherwise set check_list to false to parse through the sequence instead of
    # analyzing subseequences 
    else:
        check_list = False
    
    # If check_list is True, identify each subsequence in the list
    if check_list:
        barcode_ID, has_error = list_check_match(sequence, seq_id, seq_list, \
            file_name, barcode_dict, target_dict, count_array, error_list, tol)
    # If check list is false, use a sliding check
    else:
        barcode_ID, has_error = slide_check_match(sequence, seq_id, file_name, \
            sticky_barcode_dict, sticky_target_dict, count_array, error_list, \
            tol)
        # sliding check

    # If no barcode was found, 
    if barcode_ID == -2:
        error_list.append('{}\t{}\tno_barcode\t{}'.format(file_name, seq_id, \
            sequence))
        # Set error status to 1
        has_error = 1
    else:
        pass  
    return (barcode_ID, has_error)

def analyze_file(file_name):
    """
    This function reads in a fastq file, calculates and returns the number of 
    barcode instances in it.
    Arguments:
        file_name (str): A string representing the name of the fastq file to
        read.
            ex: test.fastq
    Returns:
        num_barcodes (int): The number of times the barcode showed up in the 
        file
    """
    with open(file_name, 'r') as file:
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
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

def analyze_sequence(sequence, seq_id, count_array, error_list, tol=0.9):
    """ analyze_sequence(sequence, seq_id, count_array, error_list, tol) 
    analyzes a given sequence for barcode and target sequences. It identifies 
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

    Errors:
        no target match: subsequence not identifiable as barcode or target
            "(file_name    seq_id, subseq, target_score[0])"

    """
    # Assume we will be able to chop the sequence into a list of subsequences
    check_list = True
    # Check for sticky ends at the beginning and end of the sequence
    head = sequence[:4]
    tail = sequence[-4:]
    # If the ends match the sticky ends, remove them
    if head == sticky_end and tail == sticky_end:
        print('Stripping sticky ends')
        stripped_seq = sequence[4:-4]
        print('stripped sequence: {}'.format(stripped_seq))
    # Otherwise set check_list to false to parse through the sequence instead of
    # analyzing subseequences 
    else:
        print('Sticky ends error, check_list set to false')
        check_list = False
    # Try splitting the sequence up into 
    print('splitting sequence')
    seq_list = stripped_seq.split(sticky_end)
    print('This is seq_list: {}'.format(seq_list))
    # Set barcode ID flag to -1 (not identified):
    barcode_ID = -1
    barcode_error = 0
    # Check to make sure the subsequences are all the correct sizes
    for subseq in seq_list:
        if len(subseq) == target_len or len(subseq) == barcode_len:
            pass
        # If a subsequence isn't the right size, set the check_list variable to 
        # False
        else:
            check_list = False
    # If check_list is True, identify each subsequence
    print('Check_list is {}'.format(check_list))
    if check_list:
        for subseq in seq_list:
            print('Checking sequence: {}'.format(subseq))
            # Check if the subsequence is in the target dictionary
            if len(subseq) == target_len:
                print('Checking target')
                target_score = check_dict(subseq, target_dict, tol)
                # If there's a match, increase the count in the array
                if target_score[0] >= tol:
                    print('Found target match')
                    count_array[target_score[1]] += 1
                # If there's no match, create a string of the fastq file name, 
                # sequence ID, error type, subsequence, and best match 
                # score to add to the error list
                else:
                    print('Didn\'t find target match')
                    error_list.append('{}\t{}\tsubseq\t{}\t{}'\
                        .format(file_name, seq_id, subseq, target_score[0]))
            # Check if the subsequence is in the barcode dictionary
            elif len(subseq) == barcode_len:
                print('Finding barcode')
                barcode_score = check_dict(subseq, barcode_dict, tol)
                # If there's a match, set the barcode ID
                if barcode_score[0] >= tol and barcode_ID == -1:
                    print('Found barcode')
                    barcode_ID = barcode_score[1]
                # If there's a match and the ID doesn't match the previously 
                # found ID, and an error hasn't previously been found, record an
                # error
                elif barcode_score[0] >= tol and barcode_ID != -1 and \
                    barcode_score[1] != barcode_ID and barcode_error == 0:
                    # Record error
                    error_list.append('{}\t{}\tbarcode_mismatch\t{}'\
                        .format(file_name, seq_id, seq))
                    # Set error flag to 1
                    barcode_error = 1
                else:
                    print('Barcode already found, or didn\'t find barcode')
                    pass
            else:
                print('Passing on barcode')
                pass
    # If check list is false, use a sliding check
    else:
        print('sliiiide check time')
        barcode_ID, barcode_error = slide_check_match(sequence, seq_id, \
            sticky_barcode_dict, sticky_target_dict, count_array, error_list, \
            tol)
        # sliding check

    # If no barcode was found, 
    if barcode_ID < 0:
        print('Didn\'t find barcode')
        error_list.append('{}\t{}\tno_barcode\t{}'.format(file_name, seq_id, \
            sequence))
    # If multiple barcodes were found, set the ID to -1
    elif barcode_error == 1:
        barcode_ID = -1
    else:
        pass  
    return barcode_ID


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
        (score, value) (tuple):
            score (float): the normalized Levenshtein distance of the best 
                matching key in the dictionary
            value (int): the value of the best matching key in the dictionary
    """
    # First check if there's a perfect match in the dictionary
    check_subseq = ref_dict.get(subseq)
    # If there is, return the value of the corresponding match
    if check_subseq != None:
        # print('found a perfect match')
        return (1.0, check_subseq)
    # If not, calculate the Levenshtein distance between the subsequence and the 
    # keys in the dictionary
    else:
        # print('calculating Levenshtein distances')
        # Store the score in a tuple (lev_score, index, key)
        score = (0.0, -1)
        for key in ref_dict.keys():
            # print('trying key: {}'.format(ref_dict.get(key)))
            # Calculate the Levenshtein distance
            lev_dist = jf.levenshtein_distance(subseq, key)
            # Normalize the Levenshtein distance to [0, 1]
            lev_score = 1.0 - lev_dist / len(key)
            # If the score is greater than the tolerance, return the index of 
            # the corresponding key
            if lev_score >= tol:
                # print('found match within tolerance')
                return (lev_score, ref_dict.get(key))
            # If not, check if the score is better than the last saved score.
            else:
                # Save the score and index if the score is better
                if lev_score > score[0]:
                    # print('found better match')
                    score = (lev_score, ref_dict.get(key))
                else:
                    # print('found worse match')
                    pass
        # Getting to this point means a match wasn't found, so return the score
        # tuple
        return score

def slide_check_match(sequence, seq_id, sticky_barcode_dict, sticky_target_dict, 
    count_array, error_list, tol):
    print('sliiiide check time')
    sticky_target_len = int(sticky_len * 2 + target_len)
    sticky_barcode_len = int(sticky_len * 2 + barcode_len)
    start = 0
    check_target = 1
    barcode_ID = -1
    barcode_error = 0
    # Make sure the start point is within the sequence
    print('Starting while loop')
    i = 1
    while start + 4 < len(sequence):
        print('@ beginning of loop: {}'.format(i))
        print('start: {}'.format(start))
        # Set the start points for deletion and insertion
        start_del = start - 1
        start_ins = start + 1
        # If the deletion start point is out of bounds, set it to 0
        if start_del < 0:
            start_del = start
        else:
            pass
        # Calculate the target end index
        target_end = start + sticky_target_len
        target_end_del = target_end - 1
        target_end_ins = target_end + 1
        print('target end: {}'.format(target_end))
        print('target end1: {}'.format(target_end_ins))
        # If the end indexes are out of bounds, cap them at the end of the \
        # sequence
        if target_end_del > len(sequence):
            print('Target end del > sequence length, capping all ends')
            target_end = len(sequence)
            target_end_del = len(sequence)
            target_end_ins = len(sequence)
        elif target_end_del == len(sequence):
            print('Target end del = sequence length, capping reg/ins end')
            target_end = len(sequence)
            target_end_ins = len(sequence)
        elif target_end == len(sequence):
            print('Target end = sequence length, capping ins end')
            target_end_ins = len(sequence)
        else:
            pass

        # Check matches for the subsequence, and the subsequence shifted over 
        # by 1
        target_subseq = sequence[start:target_end]
        target_subseq_del = sequence[start_del:target_end_del]
        target_subseq_ins = sequence[start_ins:target_end_ins]
        print('target subseq: {}'.format(target_subseq))
        print('target subseq del: {}'.format(target_subseq_del))
        print('target subseq ins: {}'.format(target_subseq_ins))
        # Check if the sequence matches the sticky target
        target_score = check_match(target_subseq, sticky_target_dict, tol)
        target_score_del = check_match(target_subseq_del, sticky_target_dict, tol)
        target_score_ins = check_match(target_subseq_ins, sticky_target_dict, tol)
        print('target score: {}'.format(target_score))
        print('target score_del: {}'.format(target_score_del))
        print('target score_ins: {}'.format(target_score_ins))
        # Pick the better of the two matches
        if target_score_ins[0] > target_score[0] and \
        target_score_ins[0] > target_score_del[0]:
            print('insert score is best')
            target_score = target_score_ins
            target_subseq = target_subseq_ins
            print('saving insert start/end point')
            # Save the the insert indices as start and end indices 
            target_start = start_ins
            target_end = target_end_ins
        # If the first subsequence is better, save the start and end indices
        elif target_score_del[0] > target_score[0] and \
        target_score_del[0] > target_score_ins[0]:
            print('del score is best')
            target_score = target_score_del
            target_subseq = target_subseq_del
            print('saving del start/end point')
            # Save the the delete indices as start and end indices 
            target_start = start_del
            target_end = target_end_del
        else:
            print('unshifted score is better')
            print('saving unshifted start point')
            target_start = start
        print('saved start point {}'.format(target_start))
        print('saved end point {}'.format(target_end))
        # If there's a match, add it to the count
        if target_score[0] >= tol:
            print('found a good target match')
            # Update the count
            print('updating count array')
            count_array[target_score[1]] += 1
            # Set the start to the end
            print('setting start to target end')
            start = target_end - 4
        else:
            print('didn\'t find a good target match')
            print('Checking barcodes now')
            # Calculate the barcode end index
            barcode_end = start + sticky_barcode_len
            barcode_end_del = barcode_end - 1
            barcode_end_ins = barcode_end + 1
            print('start: {}'.format(start))
            print('barcode end: {}'.format(barcode_end))
            print('barcode end1: {}'.format(barcode_end_ins))
            # If the end index is out of bounds, set it to the end of the 
            # sequence
            if barcode_end_del > len(sequence):
                print('barcode end del > sequence length, capping all ends')
                barcode_end = len(sequence)
                barcode_end_del = len(sequence)
                barcode_end_ins = len(sequence)
            elif barcode_end_del == len(sequence):
                print('barcode end del = sequence length, capping reg/ins end')
                barcode_end = len(sequence)
                barcode_end_ins = len(sequence)
            elif barcode_end == len(sequence):
                print('barcode end = sequence length, capping ins end')
                barcode_end_ins = len(sequence)
            else:
                pass
            # Check matches for the subsequence, and the subsequence shifted over 
            # by 1
            barcode_subseq = sequence[start:barcode_end]
            barcode_subseq_del = sequence[start_del:barcode_end_del]
            barcode_subseq_ins = sequence[start_ins:barcode_end_ins]
            print('barcode subseq: {}'.format(barcode_subseq))
            print('barcode subseq del: {}'.format(barcode_subseq_del))
            print('barcode subseq ins: {}'.format(barcode_subseq_ins))
            # Check if the sequence matches the sticky target
            barcode_score = check_match(barcode_subseq, sticky_barcode_dict, tol)
            barcode_score_del = check_match(barcode_subseq_del, \
                sticky_barcode_dict, tol)
            barcode_score_ins = check_match(barcode_subseq_ins, \
                sticky_barcode_dict, tol)
            print('barcode score: {}'.format(barcode_score))
            print('barcode score_del: {}'.format(barcode_score_del))
            print('barcode score_ins: {}'.format(barcode_score_ins))
            # Pick the better of the two matches
            if barcode_score_ins[0] > barcode_score[0] and \
            barcode_score_ins[0] > barcode_score_del[0]:
                print('insert score is best')
                barcode_score = barcode_score_ins
                barcode_subseq = barcode_subseq_ins
                print('saving insert start/end point')
                # Save the the insert indices as start and end indices 
                barcode_start = start_ins
                barcode_end = barcode_end_ins
            # If the first subsequence is better, save the start and end indices
            elif barcode_score_del[0] > barcode_score[0] and \
            barcode_score_del[0] > barcode_score_ins[0]:
                print('del score is best')
                barcode_score = barcode_score_del
                barcode_subseq = barcode_subseq_del
                print('saving del start/end point')
                # Save the the delete indices as start and end indices 
                barcode_start = start_del
                barcode_end = barcode_end_del
            else:
                print('unshifted score is better')
                print('saving unshifted start point')
                barcode_start = start
            print('saved barcode start: {}'.format(barcode_start))
            print('saved barcode end: {}'.format(barcode_end))
            # Check to see if the target or barcode match is better
            # If the target match is better
            if target_score[0] >= barcode_score[0]:
                print('target is better than barcode score')
                print('recording match error')
                error_list.append('{}\t{}\tsubseq\t{}\t{}'.format(file_name, \
                        seq_id, target_subseq, target_score[0]))
                # Set the start to the target end
                print('setting start to target end - 4')
                start = target_end - 4
            else:
                if barcode_score[0] >= tol: 
                    print('found barcode match')
                    if barcode_ID == -1:
                        print('setting barcode ID to {}'.format(barcode_score[1]))
                        barcode_ID = barcode_score[1]
                    elif barcode_score[0] >= tol and barcode_ID != -1 and \
                        barcode_score[1] != barcode_ID and barcode_error == 0:
                        print('found barcode: {}'.format(barcode_score[1]))
                        print('already found a barcode, doesn\'t match')
                        print('recording barcode mismatch')
                        # Record error
                        error_list.append('{}\t{}\tbarcode_mismatch\t{}'\
                            .format(file_name, seq_id, sequence))
                        # Set error flag to 1
                        print('setting barcode error flag to 1')
                        barcode_error = 1
                    else:
                        print('Barcode already found, or didn\'t find barcode')
                        pass
                else:
                    print('Barcode is better match, but not a good match')
                    print('recording error')
                    error_list.append('{}\t{}\tsubseq\t{}\t{}'\
                        .format(file_name, seq_id, barcode_subseq, \
                            barcode_score[0])) 
                # Set the start to the barcode end
                print('setting start to barcode end - 4')
                start = barcode_end - 4
        i += 1

    return(barcode_ID, barcode_error)


                
             


def replace_check(sequence, barcode_dict, target_dict, tol):
    for key in target_dict.keys():
        sequence = sequence.replace(key, 't{}'.format(target_dict.get(key)))
    for key in barcode_dict.keys():
        sequence = sequence.replace(key, 't{}'.format(barcode_dict.get(key)))
    return sequence



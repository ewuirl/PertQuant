import sys
# sys.path.append('../')
from PertQuant.simCRN.ml_nupack import gen_complement
import concurrent.futures
import numpy as np
import argparse
from pathlib import Path
import os
import time
import matplotlib.pyplot as plt 
import beepy as bp


def read_settings(settings_path):
    """
    read_settings(settings_path)

    This function takes in a path to a settings file and reads it to obtain 
    settings to use for counting subsequence matches.

    Arguments:
        settings_path (str): /path/to/count/settings/file. This is the count 
        settings

    Returns: 
        (save_path, target_file_path, barcode_file_path, handle_repeat_error, \
        repeat_list, n_repeat, target_sticky, barcode_sticky, sticky_end, min_len, \
        record) (tuple)
            save_path (str): /path/to/save/folder/. This is the path to the 
            folder to save the count data files to.
            target_file_path (str): This is the path to a text file containing
            the target sequences. 
            barcode_file_path (str): 
            handle_repeat_error (bool): 
            repeat_list (list): 
            n_repeat (int): 
            target_sticky (bool): 
            barcode_sticky (bool): 
            sticky_end (str): 
            min_len (int): 
            record (bool): 
    """
    # Optional settings, set the default value
    stride = 1
    min_stretch = 5

    with open(settings_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            # Remove the new line character
            line = line.rstrip("\n")
            # Split each line into a setting name and the setting
            line_list = line.split(": ", 1)
            # Handle the mandatory settings.
            if line_list[0] == "target_file_path":
                assert line_list[1] != "", "Target file path (str) must be provided in settings file."
                target_file_path = line_list[1]
            elif line_list[0] == "handle_repeat_error":
                assert line_list[1] != "", "handle_repeat_error (bool) must be provided in settings file."
                if line_list[1] == "True":
                    handle_repeat_error = True
                else:
                    handle_repeat_error = False
            # Handle optional settings. If they are not provided, default settings 
            # are used.
            elif line_list[0] == "# save_path":
                if line_list[1] != "":
                    save_path = line_list[1]
                else:
                    save_path = "default"
            elif line_list[0] == "# barcode_file_path":
                if line_list[1] != "":
                    barcode_file_path = line_list[1]
                else:
                    barcode_file_path = ""
            elif line_list[0] == "# repeat_list":
                if line_list[1] != "":
                    repeat_str = line_list[1]
                    # Make the list
                    repeat_list = repeat_str.split()
                else:
                    repeat_list = ["TG", "ATT"]
            elif line_list[0] == "# n_repeat":
                if line_list[1] != "":
                    n_repeat = line_list[1]
                else:
                    n_repeat = 3
            elif line_list[0] == "# target_sticky":
                if line_list[1] != "" and line_list[1] == "True":
                    target_sticky = True
                else:
                    target_sticky = False
            elif line_list[0] == "# barcode_sticky":
                if line_list[1] != "" and line_list[1] == "True":
                    barcode_sticky = True
                else:
                    barcode_sticky = False
            elif line_list[0] == "# sticky_end":
                if line_list[1] != "":
                    sticky_end = line_list[1].rstrip(" ")
                else:
                    sticky_end = "TGCA"
            elif line_list[0] == "# min_len":
                if line_list[1] != "":
                    min_len = int(line_list[1])
                else:
                    min_len = 5
            elif line_list[0] == "# record":
                if line_list[1] != "" and line_list[1] == "True":
                    record = True
                else:
                    record = False
            elif line_list[0] == "# stride":
                if line_list[1] != "":
                    stride = int(line_list[1])
                else:
                    pass
            elif line_list[0] == "# min_stretch":
                if line_list[1] != "":
                    min_stretch = int(line_list[1])
                else:
                    pass
            else:
                pass
    return(save_path, target_file_path, barcode_file_path, handle_repeat_error, \
        repeat_list, n_repeat, target_sticky, barcode_sticky, sticky_end, min_len, \
        record, stride, min_stretch)


def make_dictionary(file_name, sticky_end="TGCA", is_sticky=True):
    """
    make_dictionary(file_name, sticky_end, is_sticky=True)

    Takes in the name of a file that contains sequences (each on a different 
    line), and a sticky end sequence. Creates 2 dictionaries of the sequences.
    The dictionaries have the sequence ID as keys. One has the sequences as 
    values, and the other has complements as values. This function also 
    identifies the length of the sequences, and counts the number of sequences
    it read (does not include complements).

    Arguments:
        file_name (str): A string representing the file of sequences to create
            dictionaries from. Each sequence should be in its own line, with no
            empty lines.
        sticky_end (str): a string representing a sticky end sequence. Defaults
            to "TGCA".
        is_sticky (bool): True if the sequences in the sequence file are 
            sandwiched by sticky ends. False if not. Defaults to True.

    Returns: 
        (dictionary, comp_dictionary, seq_list, seq_len, num_seq) (tuple):
            dictionary (dict): A dictionary with sequence ID as keys and 
            sequences as values.
            comp_dictionary (dict): A dictionary with sequence ID as keys and 
            complementary sequences as values.
            seq_comp_list (list): A list with all the sequences and their 
                complements.
            seq_len (int): The length of the sequences/complements in dictionary.
            num_seq (int): The number of sequences in the file.

    Errors: 
        ValueError: raises an exception if an empty line is found in the
            sequence file. Exits script.
        AssertionError: raises an exception if a sequence with a different 
            length is found. Exits script.
    """
    # Create empty dictionaries and an empty list
    dictionary = {}
    comp_dictionary = {}
    seq_comp_list = []
    # seq_list = []
    # comp_list = []

    # Determine the length of the sticky sequence
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
            ref_len = len(lines[0].rstrip('\n'))-2*sticky_len
        # If the file has sequences without sticky ends
        else:
            ref_len = len(lines[0].rstrip('\n'))
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
                    # If the sequences don't have sticky ends, add the sticky
                    # ends to get the sticky end sequence
                    else:
                        sequence = whole_sequence
                    # Figure out if the sequence is the right length
                    assert len(sequence) == ref_len,\
                    'AssertionError: Sequence {} is a different length than '\
                    .format(i) + 'expected in sequence file {}.'.format(file_name) + \
                    '\nExpected length {}, got {}.'.format(ref_len, len(sequence))
                    # Generate complements

                    comp = gen_complement(sequence)
                    # Add the sequences to the dictionaries
                    dictionary[i] = sequence
                    comp_dictionary[i] = comp
                    seq_comp_list.append(sequence)
                    seq_comp_list.append(comp)
                    # seq_list.append(sequence)
                    # comp_list.append(comp)
                    # Increase the sequence count
                else:
                    raise ValueError('ValueError: Found an empty line in {}.'\
                        .format(file_name) + 
                        '\nMake sure sequence files do not contain empty lines.')
            # seq_comp_list = seq_list + comp_list
            return (dictionary, comp_dictionary, seq_comp_list, ref_len, num_seq)

        except (ValueError, AssertionError) as msg:
            print(msg)
            # Stop the script execution
            sys.exit(1)

def check_list_min_unique_len(target_list, min_len, barcode_list=[], \
    record=False):
    """
    check_list_min_unique_len(target_list, min_len, barcode_list=[], 
    record=False))

    Takes in a list of target sequences and their complements and determines if 
    all the subsequences of length min_len are unique. If they are, returns True.
    If not, returns False. A list of barcode sequences and their complements can
    also be provided, and the function will check if the minimum length results 
    in unique subsequences of the targets, barcodes, and their complements.

    Arguments:
        target_list (list): a list of target sequences and their complementary 
            sequences (str)
        min_len (int): an integer representing the length of the subsequences
            to check
        barcode_list (list): a list of barcode sequences and their complementary 
            sequences (str)
        record (bool): If True, returns unique (bool) and lists of tuples 
            identifying the sequences with duplicate subsequences. If False, 
            returns unique.

    Returns:
        unique (bool): True if all the subsequences of length min_len of the 
            sequences in the target_list and barcode_list are unique. False 
            otherwise.
        target_target_list (list): A list of tuples of (t1, t2, j, min_len, 
            nonunique_subseq). t1 and t2 are the indices of the targets that 
            have subsequences of length min_len that aren't unique. j is the 
            index of where the subsequence starts in sequence t1. 
            nonunique_subseq is the nonunique subsequence.
        target_barcode_list (list): A list of tuples of (t1, b, j, min_len, 
            nonunique_subseq). t1 and b are the indices of the target and 
            barcode that have subsequences of length min_len that aren't unique. 
            j is the index of where the subsequence starts in sequence t1. 
            nonunique_subseq is the nonunique subsequence.
        barcode_barcode_list (list): A list of tuples of (b1, b2, j, min_len,
            nonunique_subseq). b1 and b2 are the indices of the barcodes that 
            have subsequences of length min_len that aren't unique. j is the 
            index of where the subsequence starts in sequence b1. 
            nonunique_subseq is the nonunique subsequence.
    """
    unique = True
    target_len = len(target_list[0])
    # Create lists to store the identities of the duplicate subsequences in
    target_target_list = []
    target_barcode_list = []
    barcode_barcode_list = []

    # Figure out if there are barcodes
    if len(barcode_list) > 0:
        barcoded = True
        barcode_len = (len(barcode_list[0]))
    else:
        barcoded = False
    # Iterate through the sequences and their complements
    for t1 in range(len(target_list)): 
        # Iterate through the different subsequences
        for j in range(target_len + 1 - min_len):
            # Iterate through the other sequences
            for t2 in range(len(target_list)):
                if t2 == t1:
                    pass
                else:
                    if target_list[t1][j:j + min_len] in target_list[t2]:
                        print(f'Found target {t1}\'s subseq {target_list[t1][j:j + min_len]} of len {min_len} in target {t2}')
                        if record:
                            # Add the duplicate to the target-target list
                            target_target_list.append((t1, t2, j, min_len, target_list[t1][j:j + min_len]))
                        else:
                            pass
                        unique = False
                    else:
                        pass
            if barcoded:
                # Iterate through the barcodes if there are barcodees
                for b in range(len(barcode_list)):
                    if target_list[t1][j:j + min_len] in barcode_list[b]:
                        print(f'Found target {t1}\'s subseq {target_list[t1][j:j + min_len]} of len {min_len} in barcode {b}')
                        if record:
                            # Add the duplicate to the target-barcode list    
                            target_barcode_list.append((t1, b, j, min_len, target_list[t1][j:j + min_len]))
                        else:
                            pass
                        unique = False
                    else:
                        pass
            else:
                pass

    # Iterate through the barcodes and their complements if the sequences are 
    # barcoded. Checks barcodes against barcodes
    if barcoded:
        # Iterate through the barcodes andd their complements
        for b1 in range(len(barcode_list)):
            # Iterate through the different subsequences
            for j in range(barcode_len + 1 - min_len):
                # Iterate through the other barcodes
                for b2 in range(len(barcode_list)):
                    if b2 == b1:
                        pass
                    else:
                        if barcode_list[b1][j:j + min_len] in barcode_list[b2]:
                            print(f'Found barcode {b1}\'s subseq {barcode_list[b1][j:j + min_len]} of len {min_len} in target {b2}')
                            if record:
                                # Add the duplicate to the target-target list
                                barcode_barcode_list.append((b1, b2, j, min_len, barcode_list[b1][j:j + min_len]))
                            else:
                                pass
                            unique = False
                        else:
                            pass
    else:
        pass
    return (unique, target_target_list, target_barcode_list, barcode_barcode_list)
    

def find_min_unique_len(target_list, min_unique_len=5, barcode_list=[], \
    record=False):
    """
    find_min_unique_len(target_list, min_unique_len=5, barcode_list=[], \
    record=False)

    Takes in a list of target sequences and their complements and finds the 
    minimum length min_len such that subsequences of this length of the target and
    complementary sequences are unique. Returns this minimum length. A list of 
    barcode sequences and their complements can also be provided, and the 
    function will find the min_len such that subsequences of the targets, 
    barcodes, and their complementary sequences are unique.

    Arguments:
        target_list (list): a list of target sequences and their complementary 
            sequences (str)
        min_unique_len (int): an integer representing the minimum subsequence 
        length to check. Defaults to 4.
        barcode_list (list): a list of barcode sequences and their complementary 
            sequences (str)
        record (bool): If True, returns unique (bool) and lists of tuples 
            identifying the sequences with duplicate subsequences. If False, 
            returns unique.

    Outputs:
        min_unique_len (int): The minimum length such that subsequences of 
            length min_unique_len of the sequences in the target_list and 
            barcode_list are unique. 
        target_target_list (list): A list of tuples of (t1, t2, j, min_len, 
            nonunique_subseq). t1 and t2 are the indices of the targets that 
            have subsequences of length min_len that aren't unique. j is the 
            index of where the subsequence starts in sequence t1. 
            nonunique_subseq is the nonunique subsequence.
        target_barcode_list (list): A list of tuples of (t1, b, j, min_len, 
            nonunique_subseq). t1 and b are the indices of the target and 
            barcode that have subsequences of length min_len that aren't unique. 
            j is the index of where the subsequence starts in sequence t1. 
            nonunique_subseq is the nonunique subsequence.
        barcode_barcode_list (list): A list of tuples of (b1, b2, j, min_len,
            nonunique_subseq). b1 and b2 are the indices of the barcodes that 
            have subsequences of length min_len that aren't unique. j is the 
            index of where the subsequence starts in sequence b1. 
            nonunique_subseq is the nonunique subsequence.
    """
    found = False
    # Create lists to store the identities of the duplicate subsequences in
    target_target_list = []
    target_barcode_list = []
    barcode_barcode_list = []

    # Set the maximum unique length to either the minimum length of the targets
    # and the barcodees
    if len(barcode_list) > 0:
        max_len = min(len(target_list[0]), len(barcode_list[0]))
    else:
        max_len = len(target_list[0])

    # Subtract 1 from the min length because the loop adds one immediately
    min_unique_len -= 1

    # Check for the minimum subsequence length, starting at min_unique_len
    while found == False and min_unique_len <= max_len:
        min_unique_len += 1
        print(f'checking {min_unique_len}')
        found, tt_list, tb_list, bb_list \
            = check_list_min_unique_len(target_list, min_unique_len, \
                barcode_list=barcode_list, record=record)
        if record:
            target_target_list.append(tt_list)
            target_barcode_list.append(tb_list)
            barcode_barcode_list.append(bb_list)
        else:
            pass
        print(f'found = {found}')
    print(f'found min len: {min_unique_len}')
    
    return (min_unique_len, found, target_target_list, target_barcode_list, \
            barcode_barcode_list)

def record_nonunique_subseq_list(nonunique_subseq_list, file):
    """
    record_nonunique_subseq_list(nonunique_subseq_list, file)

    Write the nonunique subsequence information from a list to the provided file.

    Arguments:
        nonunique_subseq_list (list): A list of tuples of nonunique subsequence
            information. See check_list_min_unique_len or find_min_unique_len 
            for more details.
        file (file): The file to save the nonunique subsequence information to.
    Returns:
        Nothing.
    """
    # Write the nonunique subsequence information
    if len(nonunique_subseq_list) > 0:
        for entry in nonunique_subseq_list:
            for item in entry:
                file.write(f"{item}\t")
            file.write("\n")
    # Write an empty line if the list is empty
    else:
        file.write("\n")

def record_nonunique_subseq(save_file_name, save_folder, target_target_list,\
    target_barcode_list, barcode_barcode_list):
    """
    record_nonunique_subseq(save_file_name, save_folder, target_target_list,\
    target_barcode_list, barcode_barcode_list)

    Write the nonunique subsequence information from a list to the provided file.

    Arguments:
        save_file_name (str): The custom part of the name of the file to save 
            the nonunique subsequence information to.
        save_folder (str): The path to the save folder.
        target_target_list (list): A list of tuples of (t1, t2, j, min_len, 
            nonunique_subseq). t1 and t2 are the indices of the targets that 
            have subsequences of length min_len that aren't unique. j is the 
            index of where the subsequence starts in sequence t1. 
            nonunique_subseq is the nonunique subsequence.
        target_barcode_list (list): A list of tuples of (t1, b, j, min_len, 
            nonunique_subseq). t1 and b are the indices of the target and 
            barcode that have subsequences of length min_len that aren't unique. 
            j is the index of where the subsequence starts in sequence t1. 
            nonunique_subseq is the nonunique subsequence.
        barcode_barcode_list (list): A list of tuples of (b1, b2, j, min_len,
            nonunique_subseq). b1 and b2 are the indices of the barcodes that 
            have subsequences of length min_len that aren't unique. j is the 
            index of where the subsequence starts in sequence b1. 
            nonunique_subseq is the nonunique subsequence.
    Returns:
        Nothing.
    """
    with open(f"{save_folder}/{save_file_name}_nonunique_subseq.txt", 'w') as file:
        # Record the target-target nonunique subsequences
        file.write("Target-Target Nonunique Subsequences\n")
        record_nonunique_subseq_list(target_target_list, file)
        # Record the target-barcode nonunique subsequences
        file.write("Target-Barcode Nonunique Subsequences\n")
        record_nonunique_subseq_list(target_barcode_list, file)
        # Record the barcode-barcode nonunique subsequences
        file.write("Barcode-Barcode Nonunique Subsequences\n")
        record_nonunique_subseq_list(barcode_barcode_list, file)
        

def get_avg_Q_score(Q_score_index, seq_len, lines):
    """
    get_avg_Q_score(Q_score_index, seq_len, lines):

    This function takes in an index for the line containing per-base Phred
    quality scores encoded as 0-93 in ASCII 33-126. It also takes in a sequence
    length, and the lines of a FASTQ file. It converts the per-base Q scores to 
    per-base P(error) probabilities. It uses these to calculate the average 
    P(error) for the sequence. This is used to calculate the average Q score for
    the sequence, which is returned.

    Arguments:
        Q_score_index (int): The index of the line containing the Q score data.
        seq_len (int): The length of the sequence being analyzed.
        lines (list): A list of lines read in from a FASTQ file.

    Outputs:
        avg_Q_score (float): The average Q score for the entire sequence.
    """

    # Read in the quality score
    Q_score = lines[Q_score_index].rstrip("\n")

    # Add up the per base P(error)
    tot_p_error = 0

    # Find the Q score and P(error) for each base
    for char in Q_score:
        # Q scores are encoded from 0-93 using ASCII 33-126.
        Q = ord(char) - 33
        tot_p_error += 10.0 ** (-Q/10.0)

    # Calculate the average P(error)
    avg_p_error = tot_p_error / seq_len
    
    # Convert to an average Q_score
    avg_Q_score = -10.0*np.log10(avg_p_error)
    return avg_Q_score


def get_barcode_ID(sequence, index, lines, barcode_list):
    pass

def count_matches_seq(target, targetc, sequence, target_len, max_num_subseqs, \
    array_len, min_len):
    """
    count_matches(target, targetc, sequence, target_len, max_num_subseqs, \
    array_len, min_len)

    This function counts the subsequence instances of a target sequence and its
    complement in a provided sequence, and returns these counts in separate 
    arrays for the target and its complement.

    Arguments:
        target (str): a string of the target sequence to count subsequences for
        targetc (str): a string of the complementary target sequence to count
            subsequences for
        sequence (str): a string of the sequence to check for subsequences in
        target_len (int): the length of the target sequence
        max_num_subseqs (int): the number of the subsequences of the minimum
            length that will be counted
        array_len (int): the total number of subsequences that will be 
            checked (based on min_len and target_len)
        min_len (int): the minimum subsequence length that will be checked

    Outputs:
        count_array (numpy arr): A numpy array containing the subsequence counts
            of the target sequence in the provided sequence. Starts with smaller 
            subsequence lengths, and increases to the length of the target sequence.
            Let the target sequence be indexed starting at 0 from the left. For a 
            particular subsequence length, the counts are ordered by smallest to 
            largest starting index of the subsequence. For example, for a target
            sequence of ACGT, the order of the counts of subsequences of length 2 is 
            AC, CG, GT.
        comp_count_array (numpy arr): A numpy array containing the subsequence 
            counts of the target complementary sequence in the provided sequence. 
            This array is ordered in a similar manner to count_array.
    """
    # Make arrays to store 
    count_array = np.zeros(array_len)
    comp_count_array = np.zeros(array_len)
    # Create array index 
    array_index = -max_num_subseqs - 1
    for i in range(max_num_subseqs):
        # Set the length of the subsequence to store
        n = i + min_len
        num_subseqs = target_len + 1 - n
        # print(f"subseq len:{n}")
        array_index = array_index + num_subseqs + 1
        for j in range(num_subseqs):
            # print(f"starting index: {j}")
            count_array[array_index+j] = sequence.count(target[j:j + n])
            comp_count_array[array_index+j] = sequence.count(targetc[j:j + n])
    return (count_array, comp_count_array)

def make_target_subseq_array(target_list, target_len, min_len):
    num_seq = len(target_list)
    target_subseq_array = np.empty((num_seq, target_len-min_len+1),dtype='<U6')
    for i in range(num_seq):
        for j in range(target_len-min_len+1):
            target_subseq_array[i,j] = target_list[i][j:j+5]
    return target_subseq_array

def tag_subsequences(read_sequence, target_subseq_array, num_targets, min_len, \
    stride=1):
    num_steps = int((len(read_sequence)-min_len)/stride)
    indices = np.linspace(0,(num_steps)*stride, num_steps+1, dtype=int)
    tag_array = np.full((2,num_steps+1), np.nan)
    tag_index = 0
    for i in indices:
        identity = np.asarray(target_subseq_array == read_sequence[i:i+min_len]).nonzero()
        if len(identity[0])> 0:
            tag_array[0,tag_index] = identity[1][0]
            tag_array[1,tag_index] = identity[0][0]
        else:
            pass
        tag_index += 1
    return tag_array

def ID_sequences(tag_array, num_targets, min_stretch):
    clean_tag_array = tag_array[:,~np.all(np.isnan(tag_array), axis=0)].astype(int)
    # An array to store the counts in
    sequence_counts_arr = np.zeros(num_targets, dtype=int)
    if clean_tag_array.shape[1] >= min_stretch:
        previous_tag = clean_tag_array[0,0]
        previous_sequence = clean_tag_array[1,0]
        fragment_len = 1
        for i in range(np.shape(clean_tag_array)[1]-1):
            current_tag = clean_tag_array[0,i+1]
            current_sequence = clean_tag_array[1,i+1]
            if current_sequence == previous_sequence and previous_tag < current_tag:
                # Add to the fragment
                fragment_len += 1                
            else:
                # Add to the count if the fragment is long enough
                if fragment_len >= min_stretch:
                    sequence_counts_arr[previous_sequence] += 1
                else:
                    pass
                # Reset the fragment length
                fragment_len = 1
                # Reset the sequence
                previous_sequence = current_sequence
            previous_tag = current_tag
    else:
        pass
    return sequence_counts_arr

def check_repeat_errors(sequence, repeat_list, n_repeat):
    """
    check_repeat_errors(sequence, repeat_list, n_repeat)

    Determines whether a sequence has a repeat error.

    Arguments:
        sequence (str): a string representing the sequence to check 
        repeat_list (list): a list of strings of the repeat sequences to check for
        n (int): the minimum repeat number of the repeat sequence to look for, eg
            if the repeat sequence is "TG", the sequence is considered to contain
            a repeat error if it contains a sequence of "TG" repeated n times.

    Returns:
        has_error (int): 1 if a repeat sequence was found, 0 otherwise.
    """
    has_error = 0
    for repeat in repeat_list:
        if n_repeat*repeat in sequence:
            has_error = 1
        else:
            pass
    return has_error

def get_seq_info(index):
    """
    analyze_seq(index)

    Takes in an index of the line containing a sequence ID of a sequence to be
    analyzed. It uses this to identify the sequence ID, barcode ID, and average 
    Q score of the associated sequence. It also counts the number of 
    target/complementary subsequences in the sequence. It only takes in one 
    argument in order to be easily parallelized. It also relies on a set of 
    global arguments.

    Arguments:
        index (int): The index of the line containing the sequence ID data.

    Global Arguments;
        lines (list): A list of lines read in from a FASTQ file.
        min_len (int): the minimum subsequence length that will be checked
        target_list (list):
        barcode_list (list): A list of barcode sequences and their complements. 
            If barcode_list is an empty list, the function will not try to 
            dentify the barcode ID of the sequence, and the barcode ID will be 
            set to -1. If it is not empty, the function will try to identify 
            that barcode ID, and dictionaries of barcode and complementary 
            sequences will be required.
        handle_repeat_error (bool): If True, the function will handle repeat 
            errors. If False, all sequences will be analyzed in the same way.

        Optional:
            Barcode ID - Required if barcode_list is NOT an empty list.
            barcode_dict (dict): A dictionary with barcode ID as keys and
                barcode sequences as values. 
            barcode_comp_dict (dict): A dictionary with barcode ID as keys and
                barcode complementary sequences as values. 

            Handling Repeat Errors - Required if handle_repeat_error is True.
            read_file_name (str): The name of the FASTQ file being analyzed.
            repeat_list (list): A list of strings representing the repeat 
                sequences to check for.
            n_repeat (int) the threshold number of subsequent repeats of the 
                repeat sequences for a sequence to be considered to have a 
                repeat error.

    Outputs:
        seq_ID (str): The Nanopore sequence ID (contains information about the
            read).
        avg_Q_score (float): The base-averaged per base Q score for the read.
        seq_len (int): The length of the sequence.
        barcode_ID (int): Returns -1 if a barcode ID was not found. Returns the
            barcode ID if it was found. Defaults to -1 if the barcode dict
            is empty.
        has_repeat_error (int): If handle_repeat_error is False, returns -1. If
            handle_repeat_error is True, returns 1 if a repeat error was found, 
            0 if not.
        target_count_list (list): A list containing the subsequence count arrays
            for all of the targets. 
        target_comp_list (list): A list containing the subsequence count arrays 
            for all of the target complements. 
    """
    global lines 
    global min_len
    global target_list 
    global barcode_dict
    global handle_repeat_error
    # Read in the sequence ID
    seq_ID = lines[index].rstrip("\n")

    # Read in the sequence
    sequence = lines[index+1].rstrip("\n")
    seq_len = len(sequence)

    # Get the average Q_score
    avg_Q_score = get_avg_Q_score(index+3, seq_len, lines)

    # Find the barcode ID of the sequence 
    if len(barcode_dict) > 0:
        global barcode_comp_dict
        global barcode_list
        barcode_ID = get_barcode_ID(sequence, index, lines, barcode_list)
    # Otherwise set it to -1 (no barcode)
    else:
        barcode_ID = -1

    # If we want to handle repeat errors, check sequences for repeat errors
    if handle_repeat_error:
        # print(f'checking for repeat errors: {index}')
        global read_file_name
        global repeat_list
        global n_repeat
        # Check for repeat error 
        has_repeat_error = check_repeat_errors(sequence, repeat_list, n_repeat)
        if has_repeat_error == 1:
            print(f'file: {read_file_name[-9:]}, index: {index}, seq_len: {len(sequence)}')
        else:
            pass
    else:
        has_repeat_error = -1

    
    return(seq_ID, avg_Q_score, seq_len, barcode_ID, has_repeat_error, sequence)

def count_seq_prob(sequence):
    global target_list
    global min_len
    # Make lists to store the count arrays in
    num_targets = int(len(target_list)/2)
    target_count_list = []
    targetc_count_list = []
    # Count the target subsequences (loop through the targets)
    for i in range(num_targets):
        # Get the target and comp sequence to search for
        target = target_list[int(2*i)]
        targetc = target_list[int(2*i+1)]
        target_len = len(target)
        # Figure out the maximum number of subsequences possible
        max_num_subseqs = target_len - min_len + 1
        # Calculate the count array length
        array_len = int(max_num_subseqs * (max_num_subseqs+1)/2)
        # Count the subsequences and check for repeat errors
        count_array, comp_count_array = count_matches_seq(target, targetc, \
            sequence, target_len, max_num_subseqs, array_len, min_len)
        # Store the counts
        target_count_list.append(count_array)
        targetc_count_list.append(comp_count_array)
    return(target_count_list, targetc_count_list)

def count_seq_slide(sequence):
    global target_subseq_array
    global num_target_comp
    global min_len
    global stride
    global min_stretch
    tag_array = tag_subsequences(sequence, target_subseq_array, num_target_comp, min_len, \
    stride=stride)
    sequence_counts_arr = ID_sequences(tag_array, num_target_comp, min_stretch)
    return sequence_counts_arr

def analyze_seq_prob(index):
    seq_ID, avg_Q_score, seq_len, barcode_ID, has_repeat_error, sequence = get_seq_info(index)
    target_count_list, targetc_count_list = count_seq_prob(sequence)
    return(seq_ID, avg_Q_score, seq_len, barcode_ID, has_repeat_error, \
        target_count_list, targetc_count_list)

def analyze_seq_slide(index):
    seq_ID, avg_Q_score, seq_len, barcode_ID, has_repeat_error, sequence = get_seq_info(index)
    sequence_counts_arr = count_seq_slide(sequence)
    return(seq_ID, avg_Q_score, seq_len, barcode_ID, has_repeat_error, \
        sequence_counts_arr)

# def analyze_seq(index):
#     """
#     analyze_seq(index)

#     Takes in an index of the line containing a sequence ID of a sequence to be
#     analyzed. It uses this to identify the sequence ID, barcode ID, and average 
#     Q score of the associated sequence. It also counts the number of 
#     target/complementary subsequences in the sequence. It only takes in one 
#     argument in order to be easily parallelized. It also relies on a set of 
#     global arguments.

#     Arguments:
#         index (int): The index of the line containing the sequence ID data.

#     Global Arguments;
#         lines (list): A list of lines read in from a FASTQ file.
#         min_len (int): the minimum subsequence length that will be checked
#         target_list (list):
#         barcode_list (list): A list of barcode sequences and their complements. 
#             If barcode_list is an empty list, the function will not try to 
#             dentify the barcode ID of the sequence, and the barcode ID will be 
#             set to -1. If it is not empty, the function will try to identify 
#             that barcode ID, and dictionaries of barcode and complementary 
#             sequences will be required.
#         handle_repeat_error (bool): If True, the function will handle repeat 
#             errors. If False, all sequences will be analyzed in the same way.

#         Optional:
#             Barcode ID - Required if barcode_list is NOT an empty list.
#             barcode_dict (dict): A dictionary with barcode ID as keys and
#                 barcode sequences as values. 
#             barcode_comp_dict (dict): A dictionary with barcode ID as keys and
#                 barcode complementary sequences as values. 

#             Handling Repeat Errors - Required if handle_repeat_error is True.
#             read_file_name (str): The name of the FASTQ file being analyzed.
#             repeat_list (list): A list of strings representing the repeat 
#                 sequences to check for.
#             n_repeat (int) the threshold number of subsequent repeats of the 
#                 repeat sequences for a sequence to be considered to have a 
#                 repeat error.

#     Outputs:
#         seq_ID (str): The Nanopore sequence ID (contains information about the
#             read).
#         avg_Q_score (float): The base-averaged per base Q score for the read.
#         seq_len (int): The length of the sequence.
#         barcode_ID (int): Returns -1 if a barcode ID was not found. Returns the
#             barcode ID if it was found. Defaults to -1 if the barcode dict
#             is empty.
#         has_repeat_error (int): If handle_repeat_error is False, returns -1. If
#             handle_repeat_error is True, returns 1 if a repeat error was found, 
#             0 if not.
#         target_count_list (list): A list containing the subsequence count arrays
#             for all of the targets. 
#         target_comp_list (list): A list containing the subsequence count arrays 
#             for all of the target complements. 
#     """
#     global lines 
#     global min_len
#     global target_list 
#     global barcode_dict
#     global handle_repeat_error
#     # Read in the sequence ID
#     seq_ID = lines[index].rstrip("\n")

#     # Read in the sequence
#     sequence = lines[index+1].rstrip("\n")
#     seq_len = len(sequence)

#     # Get the average Q_score
#     avg_Q_score = get_avg_Q_score(index+3, seq_len, lines)

#     # Find the barcode ID of the sequence 
#     if len(barcode_dict) > 0:
#         global barcode_comp_dict
#         global barcode_list
#         barcode_ID = get_barcode_ID(sequence, index, lines, barcode_list)
#     # Otherwise set it to -1 (no barcode)
#     else:
#         barcode_ID = -1

#     # Make lists to store the count arrays in
#     num_targets = int(len(target_list)/2)
#     target_count_list = []
#     targetc_count_list = []

#     # If we want to handle repeat errors, check sequences for repeat errors
#     if handle_repeat_error:
#         # print(f'checking for repeat errors: {index}')
#         global read_file_name
#         global repeat_list
#         global n_repeat
#         # Check for repeat error 
#         has_repeat_error = check_repeat_errors(sequence, repeat_list, n_repeat)
#         if has_repeat_error == 1:
#             print(f'file: {read_file_name[-9:]}, index: {index}, seq_len: {len(sequence)}')
#         else:
#             pass
#     else:
#         has_repeat_error = -1

#     # Count the target subsequences (loop through the targets)
#     for i in range(num_targets):
#         # Get the target and comp sequence to search for
#         target = target_list[int(2*i)]
#         targetc = target_list[int(2*i+1)]
#         target_len = len(target)
#         # Figure out the maximum number of subsequences possible
#         max_num_subseqs = target_len - min_len + 1
#         # Calculate the count array length
#         array_len = int(max_num_subseqs * (max_num_subseqs+1)/2)
#         # Count the subsequences and check for repeat errors
#         count_array, comp_count_array = count_matches_seq(target, targetc, \
#             sequence, target_len, max_num_subseqs, array_len, min_len, lines)
#         # Store the counts
#         target_count_list.append(count_array)
#         targetc_count_list.append(comp_count_array)
#     return(seq_ID, avg_Q_score, seq_len, barcode_ID, has_repeat_error, target_count_list, targetc_count_list)

def make_save_folder(fastq_folder, save_path="default", count_method="prob"):
    """
    make_save_folder(fastq_folder, save_path="default")

    This function creates a save folder to save the count matches dat files. It 
    takes in a fastq_folder path and a save path. It uses the save path as the
    path of the folder to make. The default save path is to create a "counts" 
    folder in the fastq folder. If the save folder already exists, this function
    will not overwrite it.

    Arguments:
        fastq_folder (str): The path to the folder with fastq files to analyze.
        save_path (str): The path to the desired save folder. Defaults to 
            "counts", which then makes the save folder "counts" inside the fastq
            folder.

    Output:
        save_folder (str): The path to the save folder.
    """
    if save_path == "default" and count_method == "prob":
        save_folder = f"{fastq_folder}/counts"
    elif save_path == "default" and count_method == "slide":
        save_folder = f"{fastq_folder}/slide_counts"
    else:
        save_folder = save_path

    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    return save_folder

def write_header_file(save_folder, save_file_name, target_dict, \
    handle_repeat_error, count_method, barcode_dict={}, min_len=5, \
    repeat_list=["TG", "ATT"], n_repeat=3, note="", stride=1, min_stretch=5):
    """
    write_header_file(save_folder, save_file_name, target_dict, 
    handle_repeat_error, barcode_dict={}, min_len=5, repeat_list=["TG", "ATT"], 
    n_repeat=3, note="")

    Creates a header file called "count_settings_{save_file_name}.txt" in the
    provided save folder. The file contains the settings used while counting 
    subsequence matches. If a file with the same name already exists, it asks
    for confirmation (Y/N, Yes/No) to continue. If the user says no, the program
    will exit.

    Arguments:
        save_folder (str): The path to the save folder.
        save_file_name (str): The custom name to add to header file name.
        target_dict (dict): A dictionary of the target sequences and their IDs.
        handle_repeat_error (bool): Whether or not repeat errors were handled.

    Optional Arguments:
        barcode_dict (dict): A dictionary of barcodes and their IDs. Defaults to
            an empty dictionary.
        min_len (int): The minimum subsequence length to check counts for. 
            Defaults to 5.
        repeat_list (list of strings): The list of repeat sequences to check for
            if repeat errors are handled. Defaults to ["TG", "ATT"].
        n_repeat (int): The number of consecutive repeats to flag as a repeat 
            error. Defaults to 3.
        note (str): A string representing a note to add to the header file. 
            Defaults to an empty string.

    Outputs:
        header_file_name (str): The name of the header file that was written.

    """
    # Check if the header file exists
    header_file_name = f"count_settings_{save_file_name}.txt"
    header_file_path = f"{save_folder}/{header_file_name}"
    if os.path.isfile(header_file_path):
        continue_ans = str(input(f"Header file {header_file_name} already exists. Continuing will overwrite this file and may overwrite previous count data files in the same folder. Continue? (Y/N) \n"))
        yes_list = ["Y", "y", "Yes", "yes"]
        no_list = ["N", "n", "No", "no"]
        if continue_ans in yes_list:
            pass
        elif continue_ans in no_list:
            sys.exit()
    else:
        pass

    # Open header file
    header_file = open(header_file_path, "w")

    # Write the target sequence info to the save file
    header_file.write(f"# Targets = {len(target_dict)} \n")
    for key in target_dict:
        header_file.write(f'{key} {target_dict[key]} \n')

    # Write the barcode sequence info to the save file
    header_file.write(f"# Barcodes = {len(barcode_dict)} \n")
    # If there are barcodes, write the keys and sequences.
    if len(barcode_dict) == 0:
        pass
    else:
        for key in barcode_dict:
            header_file.write(f'{key} {barcode_dict[key]} \n')

    # Write the analysis settings info
    header_file.write(f'# Analysis Settings \n')
    header_file.write(f'count_method = {count_method} \n')
    if count_method == 'slide':
        header_file.write(f'stride = {stride} \n')
        header_file.write(f'min_stretch = {min_stretch} \n')
    else:
        pass
    header_file.write(f'min_len = {min_len} \n')
    header_file.write(f'handle_repeat_error = {handle_repeat_error} \n')
    if handle_repeat_error:
        pass
    else:
        # If repeat errors are not handled, set repeat_list to empty and 
        # n_repeat to 0
        repeat_list = []
        n_repeat = 0
    if repeat_list == []:
        header_file.write(f'repeat_list = {repeat_list} \n')
    else:
        header_file.write('repeat_list =')
        for repeat in repeat_list:
            header_file.write(f" {repeat}")
        header_file.write("\n")
    header_file.write(f'n_repeat = {n_repeat} \n')

    # Add note section
    header_file.write("# Notes \n")
    if note == "":
        pass
    else:
        header_file.write(f"{note} \n")
    header_file.close()

    print(f"Wrote header file: {header_file_name}")
    return header_file_name

def init_save_file(read_file_path, save_folder, save_file_name, header_file_name, \
    count_method = "prob"):
    """
    init_save_file(read_file_path, save_folder, save_file_name, header_file_name)

    This function initiates a .dat save file to store match count data in. It 
    writes the name of the header file containing information about the match 
    counting analysis into the first line of the file. It returns the save file 
    object. Each fastq file gets its own save file.

    Arguments:
        read_file_path (str): the path to the read file that is about to be
            analyzed
        save_folder (str): the path to the save folder to save results in
        save_file_name (str): the custom name provided at the beginning of the 
            run
        header_file_name (str): the name of the header file for the match 
            counting run

    Outputs:
        save_file (file): the save file that match count data for the read file
            will be saved to
    """
    global target_num_seq
    read_file_list = read_file_path.split("/")
    read_file_name = read_file_list[-1].replace(".fastq","")

    if count_method == "prob":
        extension = "counts"
    else:
        extension = "slide_counts"

    # Open the save file
    save_file = open(f"{save_folder}/{read_file_name}_{save_file_name}_{extension}.dat", "w")

    if count_method == "prob":
        # Write where the analysis settings are saved
        save_file.write(f"# Data info and analysis settings in {header_file_name}\n")
    else:
        save_file.write(f"seq_ID\tread\tchannel\tread_time\tavg_q_score\tseq_len\tbarcode_ID\thas_repeat_error")
        for i in range(target_num_seq):
            save_file.write(f"\t{i}\t{i}C")
        save_file.write(f"\n")
    return save_file

def write_dat_file_prob(save_file, seq_ID, avg_Q_score, seq_len, barcode_ID, \
    has_repeat_error, target_count_list, targetc_count_list):
    """
    write_dat_file(save_file, seq_ID, seq_len, avg_Q_score, barcode_ID, 
    has_repeat_error, target_count_list, targetc_count_list)

    This function writes count match data to a save file. Each fastq file gets 
    its own save file.

    Arguments:
        save_file (file):
        seq_ID (int):
        seq_len (int):
        avg_Q_score (float):
        barcode_ID (int):
        has_repeat_error (int):
        target_count_list (list):
        targetc_count_list (list):
    
    Output:
        None. Writes data to a save file.
    """
    save_file.write(f"{seq_ID}\n")
    save_file.write(f"{avg_Q_score} {seq_len} {barcode_ID} {has_repeat_error}\n")
    for i in range(len(target_count_list)):
        # Write target subsequence counts to save file
        target_array = target_count_list[i]
        for count in target_array:
            save_file.write(f"{int(count)} ")
        save_file.write("\n")

        # Write target complement subsequence counts to save file
        targetc_array = targetc_count_list[i]
        for count in targetc_array:
            save_file.write(f"{int(count)} ")
        save_file.write("\n")

def write_dat_file_slide(save_file, seq_ID, avg_Q_score, seq_len, barcode_ID, \
    has_repeat_error, sequence_counts_arr):
    """
    write_dat_file(save_file, seq_ID, seq_len, avg_Q_score, barcode_ID, 
    has_repeat_error, target_count_list)

    This function writes count match data to a save file. Each fastq file gets 
    its own save file.

    Arguments:
        save_file (file):
        seq_ID (int):
        seq_len (int):
        avg_Q_score (float):
        barcode_ID (int):
        has_repeat_error (int):
        target_count_list (list):
        targetc_count_list (list):
    
    Output:
        None. Writes data to a save file.
    """
    seq_ID_list = seq_ID.split()
    short_seq_id = seq_ID_list[0]
    read = seq_ID_list[2].split('=')[1]
    channel = seq_ID_list[3].split('=')[1]
    read_time_str = seq_ID_list[4].split('=')[1]
    save_file.write(f"{short_seq_id}\t{read}\t{channel}\t{read_time_str}\t{avg_Q_score}\t{seq_len}\t{barcode_ID}\t{has_repeat_error}")
    for count in sequence_counts_arr:
        save_file.write(f"\t{int(count)}")
    save_file.write("\n")

# def analyze_fastq_file(read_file_path, save_folder, save_file_name, header_file_name):
#     """ """
#     # Read in the data from a fastq file
#     read_file =  open(read_file_path, 'r')
#     lines = read_file.readlines()
#     # index_list = range(int(len(lines)/4))
#     index_arr = np.arange(0,len(lines),4)
#     # Create a save file
#     save_file = init_save_file(read_file_path, save_folder, save_file_name, \
#         header_file_name)
#     # Parallelize the sunsequence counting
#     with concurrent.futures.ProcessPoolExecutor() as executor:
#         results = executor.map(analyze_seq, index_arr)
#         for result in results:
#             # Unpack the subsequence counts
#             seq_ID, avg_Q_score, barcode_ID, has_repeat_error, \
#             target_count_list, targetc_count_list = result
#             # Write the results to the save file
#             write_dat_file(save_file, seq_ID, avg_Q_score, barcode_ID, \
#                 has_repeat_error, target_count_list, targetc_count_list)
#     read_file.close()
#     print("analyzed file")
#     save_file.close()


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="Analyzes fastq files for \
        perfect sequence matches of different lengths")
    parser.add_argument("fastq_folder", type=str, help="The path to the folder \
        with the fastq files to analyze. Fastq files in subdirectories will also \
        be analyzed.")
    parser.add_argument("save_file", type=str, help="The name extension to add to \
        save file names that results are saved to.")
    parser.add_argument("--settings", type=str, help="The path to the file \
        containing count analysis settings.")
    parser.add_argument("--method", type=str, help="The count method to use. \
        Available options are 'prob' and 'slide'. Defaults to prob.")
    parser.add_argument("--barcoded", type=bool, help="If True, the program will \
        try to identify the barcode ID of each sequence.")
    parser.add_argument("--serial", type=str, help="If True, counts matches serially \
        instead of in parallel. Defaults to False.")
    parser.add_argument("--note", type=str, help="Adds a note to the header file.")
    parser.add_argument("--prog", type=str, help="If True, prints progress \
        messages. Defaults to False.")
    parser.add_argument("--beep", type=int, help="Plays a sound using beepy when \
        the program finishes running. To pick a sound, provide an integer from \
        1-7. To not play a sound, set to 0. Defaults to 1.")
    parser.add_argument("--min_stretch", type=int, help="Overrides min_stretch in \
        settings file.")
    parser.add_argument("--save_path", type=str, help="Overrides save_path in \
        settings file.")
    args = parser.parse_args()

    # Parse the arguments
    fastq_folder = args.fastq_folder
    save_file_name = args.save_file
    if args.settings:
        settings_path = args.settings
    else:
        settings_path_search = Path(f"{fastq_folder}/counts").glob("*_settings.txt")
        settings_path_list = [str(path) for path in settings_path_search]
        # If there are no default settings files 
        if len(settings_path_list) == 0:
            settings_path = str(input(f"No settings file found in counts folder. Please enter the path to the settings file.\n"))
        # If there is 1 default settings file
        elif len(settings_path_list) == 1:
            settings_path = settings_path_list[0]
        # If there are multiple settings files.
        else:
            for i in range(len(settings_path_list)):
                print("-1: Not listed")
                print(f"{i}: {settings_path_list[i]}")
            settings_path_index = str(input(f"Multiple settings file found in counts folder. Please enter the number of the desired settings file.\n"))
            if settings_path_index == -1:
                settings_path = str(input(f"Please enter the path to the settings file.\n"))
            else:
                settings_path = settings_path_list[settings_path_index]
    if args.method:
        count_method = args.method
    else:
        count_method = "prob"

    if args.serial:
        serial = eval(args.serial)
    else:
        serial = False
    if args.note:
        note = args.note
    else:
        note = ""
    if args.prog:
        prog = eval(args.prog)
    else:
        prog = False
    if args.beep:
        which_beep = args.beep 
    else:
        which_beep = 1

    # Read in the settings file
    save_path, target_file_path, barcode_file_path, handle_repeat_error, \
    repeat_list, n_repeat, target_sticky, barcode_sticky, sticky_end, min_len, \
    record, stride, min_stretch = read_settings(settings_path)

    if args.min_stretch:
        min_stretch = int(args.min_stretch)
    else:
        pass

    if args.save_path:
        save_path = f"{fastq_folder}/{args.save_path}"
    else:
        pass

    # Make the target dictionary
    target_dict, target_comp_dict, target_list, target_len, target_num_seq = \
    make_dictionary(target_file_path, sticky_end=sticky_end, \
        is_sticky=target_sticky)
    # Get the number of target and comp sequences
    num_target_comp = 2*target_num_seq

    # Set the count method
    if count_method=="prob":
        analyze_seq = analyze_seq_prob
        write_dat_file = write_dat_file_prob

    else:
        analyze_seq = analyze_seq_slide
        write_dat_file = write_dat_file_slide
        target_subseq_array = make_target_subseq_array(target_list, target_len, min_len)

    # Read in the barcode file if the data is barcoded and make a dictionary
    if args.barcoded:
        barcoded=True
        barcode_dict, barcode_comp_dict, barcode_list, barcode_len, barcode_num_seq = \
        make_dictionary(barcode_file_path, sticky_end=sticky_end, \
            is_sticky=barcode_sticky) 
    else:
        barcoded=False
        barcode_dict={}

    # Find min unique subsequence length
    if args.barcoded:
        subseq_min, found, target_target_list, target_barcode_list, \
            barcode_barcode_list = find_min_unique_len(target_list, \
                min_unique_len=min_len, barcode_list=barcode_list, record=record)
    else:
        subseq_min, found, target_target_list, target_barcode_list, \
            barcode_barcode_list = find_min_unique_len(target_list, \
                min_unique_len=min_len, record=record)
    
    assert found==True,"No unique minimum subsequence length found."

    # If the provided min subsequence length is greater than the minimum 
    # unique subsequence length, use it.
    if subseq_min > min_len:
        min_len = subseq_min
        print(f"Provided min subsequence is non unique. Using found minimum subsequence length: {subseq_min}")
    else:
        print(f"Using default or provided minimum subsequence length: {min_len}")

    # Create the save folder
    save_folder = make_save_folder(fastq_folder, save_path=save_path, count_method=count_method)

    # Record any nonunique subsequences
    if record:
        record_nonunique_subseq(save_file_name, save_folder, target_target_list,\
            target_barcode_list, barcode_barcode_list)
    else:
        pass

    # Create header file
    header_file_name = write_header_file(save_folder, save_file_name, target_dict, \
    handle_repeat_error, count_method, barcode_dict=barcode_dict, min_len=min_len, \
    repeat_list=repeat_list, n_repeat=n_repeat, note="", stride=stride, \
    min_stretch=min_stretch)



    # Create list of fastq files to analyze
    fastq_files = Path(fastq_folder).rglob("*.fastq")
    fastq_file_list = [str(file) for file in fastq_files] 

    # Make a list to contain sequence lengths
    len_list = []

    # Analyze all the fastq files
    index = 0
    start = time.perf_counter()

    # Serial counting
    if serial:
        print("Beginning serial subsequence counting")
        for read_file_path in fastq_file_list:
            # Read in the data from a fastq file
            read_file =  open(read_file_path, 'r')
            lines = read_file.readlines()
            # index_list = range(int(len(lines)/4))
            index_arr = np.arange(0,len(lines),4)
            # Create a save file
            save_file = init_save_file(read_file_path, save_folder, save_file_name, \
                header_file_name, count_method=count_method)
            # Parallelize the subsequence counting
            for index in index_arr:
                # seq_ID, avg_Q_score, seq_len, barcode_ID, has_repeat_error, \
                #     target_count_list, targetc_count_list = analyze_seq(index)
                # # Write the results to the save file
                # write_dat_file(save_file, seq_ID, avg_Q_score, seq_len, \
                #     barcode_ID, has_repeat_error, target_count_list, \
                #     targetc_count_list)
                result = analyze_seq(index)
                # Write the results to the save file
                write_dat_file(save_file, *result)
            read_file.close()
            save_file.close()
        if prog:
            index += 1
            print(f"Analyzed {index}/{len(fastq_file_list)} files.")
        else:
            pass
    # Parallelized counting
    else:
        print("Beginning parallelized subsequence counting")
        for read_file_path in fastq_file_list:
            # Read in the data from a fastq file
            read_file =  open(read_file_path, 'r')
            lines = read_file.readlines()
            # index_list = range(int(len(lines)/4))
            index_arr = np.arange(0,len(lines),4)
            # Create a save file
            save_file = init_save_file(read_file_path, save_folder, save_file_name, \
                header_file_name, count_method=count_method)
            # Parallelize the subsequence counting
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = executor.map(analyze_seq, index_arr)
                for result in results:
                    # Unpack the subsequence counts
                    # seq_ID, avg_Q_score, seq_len, barcode_ID, has_repeat_error, \
                    # target_count_list, targetc_count_list = result
                    # Write the results to the save file
                    # write_dat_file(save_file, seq_ID, avg_Q_score, seq_len, \
                    #     barcode_ID, has_repeat_error, target_count_list, \
                    #     targetc_count_list)

                    # Write the results to the save file
                    write_dat_file(save_file, *result)
            read_file.close()
            save_file.close()
            if prog:
                index += 1
                print(f"Analyzed {index}/{len(fastq_file_list)} files.")
            else:
                pass
    end = time.perf_counter()
    print("Time elapsed: {} second(s)".format(end-start))

    if which_beep > 0:
        bp.beep(sound=which_beep)
    else:
        pass

    # # # Test Functions

    # # Test read_setting function
    # save_path = "/Users/emilywu/Dropbox (MIT)/Emily <-> Ashwin/PerturbativeQuantification/PertQuant/analyzeFASTQ/test_count_matches/test_counts"
    # target_file_path = "/Users/emilywu/Dropbox (MIT)/Emily <-> Ashwin/PerturbativeQuantification/PertQuant/analyzeFASTQ/test_count_matches/test_targets.txt"
    # barcode_file_path = "/Users/emilywu/Dropbox (MIT)/Emily <-> Ashwin/PerturbativeQuantification/PertQuant/analyzeFASTQ/test_count_matches/test_barcodes.txt"
    # repeat_list = ["TG", "ATT"]
    # n_repeat = 3
    # sticky_end = "TGCA"
    # min_len = 5
    # record = False
    # test_settings_dict = { \
    # "test_settings_tar_bar.txt": ("counts", target_file_path, barcode_file_path, \
    #     False, repeat_list, n_repeat, False, False, sticky_end, min_len, record),
    # "test_settings_tar_sticky.txt": (save_path, target_file_path, "", \
    #     False, repeat_list, n_repeat, True, False, sticky_end, min_len, record),
    # "test_settings_bar_sticky.txt": (save_path, target_file_path, barcode_file_path, \
    #     False, repeat_list, n_repeat, False, True, sticky_end, min_len, record),
    # "test_settings_sticky_end.txt": (save_path, target_file_path, "", \
    #     False, repeat_list, n_repeat, False, False, 'sticky', min_len, record),
    # "test_settings_min_len.txt": (save_path, target_file_path, "", \
    #     False, repeat_list, n_repeat, False, False, sticky_end, 6, record),
    # "test_settings_record.txt": (save_path, target_file_path, "", \
    #     False, repeat_list, n_repeat, False, False, sticky_end, min_len, True)}

    # def test_read_settings(test_folder, test_settings_dict, print_test=False):
    #     test_files = Path(test_folder).rglob("*.txt")

    #     test_file_list = [str(file) for file in test_files] 
    #     failed_list = []
    #     for test_file in test_file_list:
    #         fail = False
    #         test_file_name = test_file.split("/")[-1]
    #         try:
    #             settings = read_settings(test_file)
    #         except:
    #             print(test_file_name)
    #         ans = test_settings_dict[test_file_name]
            
    #         # Record the failed test
    #         if ans != settings:
    #             failed_list.append(test_file_name)
    #             if print_test:
    #                 print("Answer")
    #                 print(ans)
    #                 print("Found")
    #                 print(settings)
    #             else:
    #                 pass
    #         else:
    #             pass
    #     print(f"Passed {len(test_file_list) - len(failed_list)}/{len(test_file_list)} read settings tests")
    # test_folder = "/Users/emilywu/Dropbox (MIT)/Emily <-> Ashwin/PerturbativeQuantification/PertQuant/analyzeFASTQ/test_count_matches/read_settings_tests"
    # test_read_settings(test_folder, test_settings_dict,print_test=True)

    # # Test make_dictionary
    # barcode_dict, barcode_comp_dict, barcode_list, barcode_len, barcode_num_seq = \
    # make_dictionary("test_barcodes.txt", is_sticky=True) 
    # test_dict, test_comp_dict , test_seq_comp_list, test_seq_len, test_num_seq = \
    # make_dictionary("test_barcodes_empty_begin.txt", is_sticky=True)
    # test_dict, test_comp_dict , test_seq_comp_list, test_seq_len, test_num_seq = \
    # make_dictionary("test_barcodes_empty_line.txt", is_sticky=True)
    # test_dict, test_comp_dict , test_seq_comp_list, test_seq_len, test_num_seq = \
    # make_dictionary("test_barcodes_wrong_len.txt", is_sticky=True)
    # print(test_dict)
    # print(test_comp_dict)
    # print(test_seq_comp_list)
    # print(test_seq_len)
    # print(test_num_seq)
    # target_dict, target_comp_dict, target_list, target_len, target_num_seq = \
    # make_dictionary("test_targets.txt", is_sticky=False)


    # # Test find unique min
    # print('just barcode')
    # subseq_min, found = find_min_unique_len(barcode_list)
    # print(f'found?: {found}')
    # print('just target')
    # subseq_min, found = find_min_unique_len(target_list)
    # print(f'found?: {found}')
    # print('target and barcode')
    # subseq_min, found = find_min_unique_len(target_list, barcode_list=barcode_list)
    # print(f'found?: {found}')


    # # Test finding unique min and recording the duplicates
    # print('just barcode')
    # subseq_min, found, target_target_list, target_barcode_list, \
    # barcode_barcode_list = find_min_unique_len(barcode_list, record=True)
    # print(f'found?: {found}')
    # print(f'tt_list = {target_target_list}')
    # print(f'tb_list = {target_barcode_list}')
    # print(f'bb_list = {barcode_barcode_list}')
    # print('just target')
    # subseq_min, found, target_target_list, target_barcode_list, \
    # barcode_barcode_list = find_min_unique_len(target_list, record=True)
    # print(f'found?: {found}')
    # print(f'tt_list = {target_target_list}')
    # print(f'tb_list = {target_barcode_list}')
    # print(f'bb_list = {barcode_barcode_list}')
    # print('target and barcode')
    # subseq_min, found, target_target_list, target_barcode_list, \
    # barcode_barcode_list = find_min_unique_len(target_list, \
    #     barcode_list=barcode_list, record=True)
    # print(f'found?: {found}')
    # print(f'tt_list = {target_target_list}')
    # print(f'tb_list = {target_barcode_list}')
    # print(f'bb_list = {barcode_barcode_list}')


    # # Test counting
    # fastq_folder = "/Users/emilywu/OneDrive - Massachusetts Institute of Technology/Minion/test_analyzeFASTQ/fastq_pass"
    # read_file_name = fastq_folder+"/AFO090_pass_45c91af1_0.fastq"
    # read_file_name = "test_count_matches.fastq"
    # with open(read_file_name, 'r') as read_file:
    #     lines = read_file.readlines()
    #     # index_list = range(int(len(lines)/4))
    #     index_arr = np.arange(0,len(lines),4)

    # # Globals for testing
    # min_len = 5
    # handle_repeat_error = False
    # repeat_list = ["TG", "ATT"]
    # n_repeat = 3

    # target = "TATGAGGACGAATCTCCCGCTTATA"
    # targetc = gen_complement(target)
    # index = 0
    # target_len = len(target)    
    # sequence = lines[1].rstrip("\n")


    # # Test the get_avg_Q_score function
    # Q_score_index = 3
    # Q_score_ascii = lines[3].rstrip("\n")
    # seq_len = len(Q_score_ascii)
    # avg_Q_score = get_avg_Q_score(Q_score_index, seq_len, lines)


    # # Test the count_matches_seq function
    # max_num_subseqs = target_len - min_len + 1
    # array_len = int(max_num_subseqs * (max_num_subseqs+1)/2)
    # count_array, comp_count_array = \
    # count_matches_seq(target, targetc, sequence, target_len, max_num_subseqs, \
    # array_len, min_len, lines)
    # print(count_array)
    # print(comp_count_array)


    # # Test make_save_folder function
    # default settings
    # save_folder = make_save_folder(fastq_folder, save_path="counts")
    # save_folder = make_save_folder(fastq_folder, save_path=save_path)


    # # Test make_header_file function
    # # Default settings
    # header_file_name = write_header_file(save_folder, save_file_name, target_dict, \
    # handle_repeat_error, barcode_dict={}, min_len=4, repeat_list=["TG", "ATT"], n_repeat=3, \ 
    # note="")
    # # Default settings 2
    # header_file_name = write_header_file(save_folder, save_file_name, target_dict, \
    # handle_repeat_error=True)
    # Set settings
    # handle_repeat_error=True
    # header_file_name = write_header_file(save_folder, save_file_name, target_dict, \
    # handle_repeat_error, barcode_dict=barcode_dict, min_len=min_len, \
    # repeat_list=['potato'], n_repeat=5, note=":)")


    # # Test init_save_file function
    # print("Initiating save file")
    # read_file_path = "/Users/emilywu/Dropbox (MIT)/Emily <-> Ashwin/PerturbativeQuantification/PertQuant/analyzeFASTQ/test_count_matches.fastq"
    # save_file = init_save_file(read_file_path, save_folder, save_file_name, \
    # header_file_name)
    

    # # # Test the analyze_seq function and write_dat_file
    # # Make target dictionary
    # target_dict, target_comp_dict, target_list, target_len, target_num_seq = \
    # make_dictionary("test_barcodes.txt", is_sticky=True)
    # barcode_dict = {} 
    # # Find min unique length
    # min_len, found = find_min_unique_len(target_list)
    # print(f'min_len: {min_len}')
    # print(f'found?: {found}')
    # # Write header file
    # header_file_name = write_header_file(save_folder, save_file_name, target_dict, \
    # handle_repeat_error, min_len=min_len, repeat_list=['potato'], n_repeat=5, note=":)")
    # # Initialize save file
    # print("Initiating save file")
    # read_file_path = "/Users/emilywu/Dropbox (MIT)/Emily <-> Ashwin/PerturbativeQuantification/PertQuant/analyzeFASTQ/test_count_matches.fastq"
    # save_file = init_save_file(read_file_path, save_folder, save_file_name, \
    # header_file_name)
    # # Count and write to the save file
    # print("Parallelized sequence analysis")
    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     results = executor.map(analyze_seq, index_arr)
    #     print("Writing count data")
    #     for result in results:
    #         seq_ID, avg_Q_score, barcode_ID, has_repeat_error, \
    #         target_count_list, targetc_count_list = result
    #         print("about to write")
    #         write_dat_file(save_file, seq_ID, avg_Q_score, barcode_ID, \
    #             has_repeat_error, target_count_list, targetc_count_list)
    # # Close the save file
    # save_file.close()
    

    # # # Test analyze_fastq_file
    # # Make target dictionary
    # target_dict, target_comp_dict, target_list, target_len, target_num_seq = \
    # make_dictionary("test_barcodes.txt", is_sticky=True)
    # barcode_dict = {} 
    # # Find min unique length
    # min_len, found = find_min_unique_len(target_list)
    # print(f'min_len: {min_len}')
    # print(f'found?: {found}')
    # # Write header file
    # header_file_name = write_header_file(save_folder, save_file_name, target_dict, \
    # handle_repeat_error, min_len=min_len, repeat_list=['potato'], n_repeat=5, note=":)")
    # # Initialize save file
    # print("Initiating save file")
    # read_file_path = "/Users/emilywu/Dropbox (MIT)/Emily <-> Ashwin/PerturbativeQuantification/PertQuant/analyzeFASTQ/test_count_matches.fastq"
    # save_file = init_save_file(read_file_path, save_folder, save_file_name, \
    #     header_file_name)
    # analyze_fastq_file(read_file_path, save_folder, save_file_name, header_file_name)


    # # to do
    # write a function that analyzes per file
    # write functions to read the data
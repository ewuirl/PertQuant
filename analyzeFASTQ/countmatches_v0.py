import sys
sys.path.append('../')
from simCRN.ml_nupack import gen_complement
import concurrent.futures
import numpy as np
import argparse
from pathlib import Path
import time
# import scipy.optimize
import matplotlib.pyplot as plt 
import functools
import os
import datetime as dt


# barcode  = "TATGAGGACGAATCTCCCGCTTATA"
# barcodec = "TATAAGCGGGAGATTCGTCCTCATA"
# barcode1 = "GGTCTTGACAAACGTGTGCTTGTAC"
# barcode1c = "GTACAAGCACACGTTTGTCAAGACC"
repeat_list = ["TG", "ATT"]
# seq_len = len(barcode)
seq_len = 25
min_len = 5
sticky_end = "TGCA"
sticky_len = len(sticky_end)

# seq_list = [barcode, barcodec, barcode1, barcode1c]

def check_min_unique_len(seq_list):
    # Iterate through the sequences and their complements
    for k in range(len(seq_list)):
        # Iterate through the different subsequence lengths
        for i in range(seq_len - min_len + 1):
            # Set the length of the subsequence to examine
            n = i + min_len
            # Iterate through the different subsequences
            for j in range(seq_len + 1 - n):
                # Iterate through the other sequences
                for l in range(len(seq_list)):
                    if l == k:
                        pass
                    else:
                        if seq_list[k][j:j + n] in seq_list[l]:
                            print(f'Found {k}\'s subseq of len {n} in {l}')
                        else:
                            pass

def count_sticky(file_name):
    """
    count_matches(file_name)

    Takes in a file name of a fastq file and counts perfect matches of 
    the sticky end.
    """
    sticky_count = 0
    with open(file_name, "r") as file:        
        while True:
            # Read the sequence ID
            seq_ID = file.readline()
            if not seq_ID:
                break
            # Read in the sequence and strip the new line
            sequence = file.readline().rstrip("\n")
            file.readline()
            Q_score = file.readline()
            sticky_count += sequence.count(sticky_end)
    return sticky_count

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
        Bool: True if a repeat sequence was found, False otherwise.
    """
    has_error = False
    for repeat in repeat_list:
        if n_repeat*repeat in sequence:
            has_error = True
        else:
            pass
    return has_error

# def count_matches(file_name, barcode, handle_repeat_error=False, n_repeat=6):
#     """
#     count_matches(file_name, barcode, handle_repeat_error=False, n_repeat=6

#     Takes in a file name of a fastq file and searches for perfect matches of 
#     subsequences of the barcode and its complement. It adds these into an array.
#     """
#     count_array = np.zeros((seq_len - min_len + 1, seq_len - min_len + 1))
#     comp_count_array = np.zeros((seq_len - min_len + 1, seq_len - min_len + 1))

#     # Generate the complement strand of the barcode
#     barcodec = gen_complement(barcode)
#     print(f'handle repeat error: {handle_repeat_error}')
#     with open(file_name, "r") as file:  
#         index = 0      
#         while True:
#             # Read the sequence ID
#             seq_ID = file.readline()
#             if not seq_ID:
#                 break
#             # Read in the sequence and strip the new line
#             # print(f'new sequence: {index}')
#             sequence = file.readline().rstrip("\n")
#             # If we want to handle repeat errors, check sequences for repeat errors
#             if handle_repeat_error:
#                 # print(f'checking for repeat errors: {index}')
#                 # Check for repeat error 
#                 repeat_error = check_repeat_errors(sequence, repeat_list, n_repeat)
#             else:
#                 # Don't check for repeat errors
#                 repeat_error = False
                
#             # Only analyze sequences without repeat errors
#             if repeat_error == False:
#                 file.readline()
#                 Q_score = file.readline()
#                 for i in range(seq_len - min_len + 1):
#                     # Set the length of the subsequence to store
#                     n = i + min_len
#                     # print(f"subseq len:{n}")
#                     for j in range(seq_len + 1 - n):
#                         # print(f"starting index: {j}")
#                         subseq_count = sequence.count(barcode[j:j + n])
#                         comp_count = sequence.count(barcodec[j:j + n])
#                         count_array[i,j] += subseq_count
#                         comp_count_array[i,j] += comp_count
#             else:
#                 print(f'not counting sequence: {index} in {file_name[-8:]}')
#             index += 1
#     return (count_array, comp_count_array)

def count_matches(file_name):
    """
    count_matches(file_name)

    Takes in a file name of a fastq file and searches for perfect matches of 
    subsequences of the barcode and its complement. It adds these into an array.
    """
    count_array = np.zeros((seq_len - min_len + 1, seq_len - min_len + 1))
    comp_count_array = np.zeros((seq_len - min_len + 1, seq_len - min_len + 1))

    # Generate the complement strand of the barcode
    global barcode
    # print(f'barcode: {barcode}')
    global barcodec 
    # print(f'barcodec: {barcodec}')
    global n_repeat
    global repeat_list
    # print(f'handle repeat error: {handle_repeat_error}')
    n_repeat_errors = 0
    with open(file_name, "r") as file:  
        index = 0      
        while True:
            # Read the sequence ID
            seq_ID = file.readline()
            if not seq_ID:
                break
            # Read in the sequence and strip the new line
            # print(f'new sequence: {index}')
            sequence = file.readline().rstrip("\n")
            file.readline()
            Q_score = file.readline()
            # If we want to handle repeat errors, check sequences for repeat errors
            if handle_repeat_error:
                # print(f'checking for repeat errors: {index}')
                # Check for repeat error 
                has_repeat_error = check_repeat_errors(sequence, repeat_list, n_repeat)
            else:
                # Don't check for repeat errors
                has_repeat_error = False
            # Only analyze sequences without repeat errors
            if has_repeat_error == False:

                for i in range(seq_len - min_len + 1):
                    # Set the length of the subsequence to store
                    n = i + min_len
                    # print(f"subseq len:{n}")
                    for j in range(seq_len + 1 - n):
                        # print(f"starting index: {j}")
                        subseq_count = sequence.count(barcode[j:j + n])
                        comp_count = sequence.count(barcodec[j:j + n])
                        count_array[i,j] += subseq_count
                        comp_count_array[i,j] += comp_count
            else:
                print(f'file: {file_name[-9:]}, index: {index}, seq_len: {len(sequence)}')
                n_repeat_errors += 1
            index += 1
    # print(f'done with file {file_name[-9:]}')
    return (count_array, comp_count_array, n_repeat_errors)

def deco_count_matches(func, barcode, handle_repeat_error, n_repeat=6):
    """
    deco_count_matches(func, barcode, handle_repeat_error):

    Used to specify 
        (1) which barcode to analyze the file with
        (2) whether or not the decorated function should handle repeat
            errors. Use as a decorator to decorate 
            count_matches(file_name, barcode, handle_repeat_error=False, n_repeat=6).

    Arguments:
        barcode (str): the barcode to look for subsequences for
        handle_repeat_error (bool): True or False. If True, count_matches will
            check for repeat errors.
        n (int)

    Returns:
        A decorated function with tol set to the specified tolerance.
    """
    
    @functools.wraps(func)
    def wrapper(fastq_file):
        return func(fastq_file, barcode, handle_repeat_error=handle_repeat_error, n_repeat=n_repeat)
    return wrapper

def test_deco_func(fastq_file, barcode, handle_repeat_error, n_repeat=6):
    """my doc"""
    print(f'fastq_file: {fastq_file}')
    print(f'barcode: {barcode}')
    print(f'handle_repeat_error: {handle_repeat_error}')
    print(f'n_repeat: {n_repeat}')

def record_lengths(file_name, save=False):
    """
    record_lengths(file_name)

    Takes in a file name of a fastq file and records all the lengths of the 
    sequences in the file in a list. It also counts how many sequences are the
    "correct" sequence length, and how many are incorrect.
    """
    len_list = []
    correct_len = 0
    wrong_len = 0
    num_repeat_errors = 0
    global n_repeat
    global repeat_list
    if save != False:
        save_file = open(save, "a")
    else:
        pass
    index = 0
    with open(file_name, "r") as file:        
        while True:
            # Read the sequence ID
            seq_ID = file.readline()
            if not seq_ID:
                break
            # Read in the sequence and strip the new line
            sequence = file.readline().rstrip("\n")
            file.readline()
            Q_score = file.readline()
            len_list.append(len(sequence))
            # # Check if the sequence is longer than the cutoff
            # if len(sequence) > 600 and save == False:
            #     print(f'seq_len = {len(sequence)}\n{sequence}')
            # elif len(sequence) > 600 and save != False:
            #     save_file.write(f'seq_len = {len(sequence)}\n{sequence}\n')
            # else:
            #     pass
            # Check if there's a repeat error
            has_repeat_error = check_repeat_errors(sequence, repeat_list, n_repeat)
            if has_repeat_error == False:
                pass
            else:
                print(f'seq in {file_name[-9:]} has repeat error')
                save_file.write(f'file: {file_name[-9:]}, index: {index}, seq_len: {len(sequence)}\n')
                # save_file.write(f'{sequence}\n')
                num_repeat_errors += 1
            # elif has_repeat_error and save != False:
            #     print(f'seq in {file_name[-9:]} has repeat error')
            #     save_file.write(f'file: {file_name[-9:]}, index: {index}, seq_len: {len(sequence)}\n')
            #     save_file.write(f'{sequence}\n')
            #     num_repeat_errors += 1
            # else:
            #     print(f'seq_len = {len(sequence)}\n{sequence}')
            #     num_repeat_errors += 1
            if len(sequence) % (seq_len + sticky_len) - 4 == 0:
                correct_len += 1
            else:
                wrong_len += 1
            index += 1
    if save != False:
        save_file.close()
    else:
        pass
    return (len_list, correct_len, wrong_len, num_repeat_errors)


def filter_folder_time_step(time_step, run_length, run_folder):
    """
    filter_folder_time_step(time_step, run_length, run_folder)

    Takes in a time step in minutes, and a run length in hours, and run_folder,
    the path to a folder enclosing fastq_pass and fastq_fail folders with fastq 
    files in them. It generates a folder in run_folder called 
    "time_step_{time_step}" and creates fastq_fail and fastq_pass subfolders, 
    and provides the path to the time_step folder as an output. It also figures
    out the start and end times of the run, and generates lists of 
    fastq_pass/fail files to process, as well as a time step range and list of 
    corresponding timedeltas. 

    Arguments:
        time_step (int/float): a time step in minutes to bin the data into.
        run_length (int/float): the length of the sequencing run in hours.
        run_folder (str): a string representing the path to a folder enclosing 
            fastq_pass and fastq_fail folders that have fastq_files in them.

    Returns:
        start_time (datetime): a datetime object of the start time of the run.
        end_time (datetime): a datetime object of the end time of the run.
        time_step_range (range): a range object with bounds [time_step, run_length (min)]
        time_range (list): a list of timedelta objects corresponding to the 
            values in time_step_range
        fastq_pass_file_list (list): a list of the fastq files in the fastq_pass
            folder
        fastq_fail_file_list (list): a list of the fastq files in the fastq_fail
            folder
        save_folder

    """
    run_folder_list = run_folder.split("/")
    folder_name = run_folder_list[-1] # change to -1
    folder_name_list = folder_name.split("_")
    date = folder_name_list[0]
    time = folder_name_list[1]
    start_time = dt.datetime(int(date[:4]), int(date[4:6]), \
        int(date[6:8]), hour=int(time[:2]), minute=int(time[2:]))

    # Prepare lists of time steps
    time_step_range = range(time_step,run_length*60+time_step, time_step)
    time_range = []
    for i in range(len(time_step_range)):
        delta_time = dt.timedelta(minutes=time_step_range[i])
        time_range.append(start_time+delta_time)

    end_time = time_range[-1]

    fastq_pass_files = Path(run_folder+"/fastq_pass").rglob("*.fastq")
    fastq_pass_file_list = [str(file) for file in fastq_pass_files] 
    fastq_fail_files = Path(run_folder+"/fastq_fail").rglob("*.fastq")
    fastq_fail_file_list = [str(file) for file in fastq_fail_files] 

    # Create save folders
    save_folder = run_folder+f"/time_step_{time_step}"
    save_folder_pass = save_folder + "/fastq_pass"
    save_folder_fail = save_folder + "/fastq_fail"
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
        os.makedirs(save_folder_pass)
        os.makedirs(save_folder_fail)

    return(start_time, end_time, time_step_range, time_range, \
        fastq_pass_file_list, fastq_fail_file_list, save_folder)

    # time_range = np.arange(start_time+time_step,start_time+run_length*60, time_step)

def count_matches_time_step(file_name):
    """
    count_matches(file_name)

    Takes in a file name of a fastq file and searches for perfect matches of 
    subsequences of the barcode and its complement. It adds these into an array.
    """
    count_list = []
    comp_count_list = []
    time_list = []
    count_list.append(np.zeros((seq_len - min_len + 1, seq_len - min_len + 1)))
    comp_count_list.append(np.zeros((seq_len - min_len + 1, seq_len - min_len + 1)))

    # Figure out the time step to use
    global time_step
    global start_time
    global end_time
    current_time = start_time
    time_index = 0

    # Generate the complement strand of the barcode
    global barcode
    # print(f'barcode: {barcode}')
    global barcodec 
    # print(f'barcodec: {barcodec}')
    find_current_time = True

    with open(file_name, "r") as file:  
        # index = 0      
        while True:
            # Read the sequence ID
            seq_ID = file.readline()
            if not seq_ID:
                break
            # Determine the time of the run
            seq_ID_list = seq_ID.split(" ")
            run_time_str = seq_ID_list[4]
            run_time = dt.datetime(int(run_time_str[11:15]), int(run_time_str[16:18]), \
                int(run_time_str[19:21]), hour=int(run_time_str[22:24]), \
                minute=int(run_time_str[25:27]), second=int(run_time_str[28:30]))

            # Read in the sequence and strip the new line
            # print(f'new sequence: {index}')
            sequence = file.readline().rstrip("\n")
            file.readline()
            Q_score = file.readline()

            # Find the time bracket of the first sequence
            while find_current_time:
                if run_time < current_time:
                    # Save the time of the first bracket
                    time_list.append(current_time)
                    find_current_time = False
                else:
                    current_time = current_time + dt.timedelta(minutes=time_step)
            
            # Figure out whether to add another time bracket
            # if run_time <= current_time:                
            if run_time > current_time and current_time < end_time:
                current_time = current_time + dt.timedelta(minutes=time_step)
                count_list.append(np.zeros((seq_len - min_len + 1, seq_len - min_len + 1)))
                comp_count_list.append(np.zeros((seq_len - min_len + 1, seq_len - min_len + 1)))
                time_list.append(current_time)
                time_index += 1
            else:
                count_array = count_list[time_index]
                comp_count_array = comp_count_list[time_index]

            for i in range(seq_len - min_len + 1):
                # Set the length of the subsequence to store
                n = i + min_len
                # print(f"subseq len:{n}")
                for j in range(seq_len + 1 - n):
                    # print(f"starting index: {j}")
                    subseq_count = sequence.count(barcode[j:j + n])
                    comp_count = sequence.count(barcodec[j:j + n])
                    count_array[i,j] += subseq_count
                    comp_count_array[i,j] += comp_count
            # else:
            #     print(f'file: {file_name[-9:]}, index: {index}, seq_len: {len(sequence)}')
            #     n_repeat_errors += 1
            # index += 1
    # print(f'done with file {file_name[-9:]}')
    return (count_list, comp_count_list, time_list)

def write_time_step_counts(save_folder, passfail, time_step, count_array, barcode_index, comp=False):
    save_file_base = f"{save_folder}/fastq_{passfail}/fastq_{passfail}_{time_step}_{barcode_index}"
    if comp:
        save_file = save_file_base+"_comp_counts.txt"
    else:
        save_file = save_file_base+"_counts.txt"
    with open(save_file,"w") as file:
        for i in range(seq_len - min_len + 1):
            # Set the length of the subsequence to store
            n = i + min_len
            # print(f"subseq len:{n}")
            for j in range(seq_len + 1 - n):
                file.write(str(count_array[i,j])+"\t")
            file.write("\n")

def count_matches_time_filter_parallel(time_step_range, time_range, \
    fastq_pass_file_list, fastq_fail_file_list, barcode_index):
    total_count_list_pass = []
    total_comp_count_list_pass = []
    total_count_list_fail = []
    total_comp_count_list_fail = []
    
    for i in range(len(time_step_range)):
        total_count_list_pass.append(np.zeros((seq_len - min_len + 1, seq_len - min_len + 1)))
        total_comp_count_list_pass.append(np.zeros((seq_len - min_len + 1, seq_len - min_len + 1)))
        total_count_list_fail.append(np.zeros((seq_len - min_len + 1, seq_len - min_len + 1)))
        total_comp_count_list_fail.append(np.zeros((seq_len - min_len + 1, seq_len - min_len + 1)))

    # Parallelize the counting
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(count_matches_time_step, fastq_pass_file_list)
        for result in results:
            count_array, comp_count_array, time_list = result

            for i in range(len(time_list)):
                step_index = time_range.index(time_list[i])

                # Add the counts from each file to the total count array
                total_count_list_pass[step_index] += count_array[i]
                # Add the comp counts from each file to the total comp count array
                total_comp_count_list_pass[step_index] += comp_count_array[i]

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(count_matches_time_step, fastq_fail_file_list)
        for result in results:
            count_array, comp_count_array, time_list = result

            for i in range(len(time_list)):
                step_index = time_range.index(time_list[i])

                # Add the counts from each file to the total count array
                total_count_list_fail[step_index] += count_array[i]
                # Add the comp counts from each file to the total comp count array
                total_comp_count_list_fail[step_index] += comp_count_array[i]

    # Save the results
    for i_time in range(len(time_step_range)):
        time_step = time_step_range[i_time]
        total_count_array_pass = total_count_list_pass[i_time]
        total_comp_count_array_pass = total_comp_count_list_pass[i_time]
        total_count_array_fail = total_count_list_fail[i_time]
        total_comp_count_array_fail = total_comp_count_list_fail[i_time]
        write_time_step_counts(save_folder, "pass", time_step, total_count_array_pass, barcode_index)
        write_time_step_counts(save_folder, "pass", time_step, total_comp_count_array_pass, barcode_index, comp=True)
        write_time_step_counts(save_folder, "fail", time_step, total_count_array_fail, barcode_index)
        write_time_step_counts(save_folder, "fail", time_step, total_comp_count_array_fail, barcode_index, comp=True)

if __name__ == "__main__":
    # Generate the counts
    parser = argparse.ArgumentParser(description="Analyzes fastq files for \
        perfect sequence matches of different lengths")
    parser.add_argument("fastq_folder", type=str, help="The path to the folder \
        with the fastq files to analyze.")
    parser.add_argument("save_file", type=str, help="The name of to the file \
        to save results to.")
    args = parser.parse_args()

    # Parse the arguments
    fastq_folder = args.fastq_folder
    save_file = args.save_file
    # Create list of fastq files to analyze
    fastq_files = Path(fastq_folder).rglob("*.fastq")
    fastq_file_list = [str(file) for file in fastq_files] 
    # print(fastq_file_list)

    # # Create arrays for the total counts
    total_count_array0 = np.zeros((seq_len - min_len + 1, seq_len - min_len + 1))
    total_comp_count_array0 = np.zeros((seq_len - min_len + 1, seq_len - min_len + 1))
    total_repeat_errors = 0
    total_count_array1 = np.zeros((seq_len - min_len + 1, seq_len - min_len + 1))
    total_comp_count_array1 = np.zeros((seq_len - min_len + 1, seq_len - min_len + 1))

    # Check the min unique sequence length
    # check_min_unique_len(seq_list)

    # start timing
    start = time.perf_counter()

    # Decide whether or not to handle repeat errors
    handle_repeat_error = False
    n_repeat = 5
    # Parallelize the counting 
    # Count barcode 0
    barcode  = "TATGAGGACGAATCTCCCGCTTATA"
    barcodec = gen_complement(barcode)

    print('counting')
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(count_matches, fastq_file_list)
        for result in results:
            count_array, comp_count_array, n_repeat_errors = result
            # Add the counts from each file to the total count array
            total_count_array0 += count_array
            # Add the comp counts from each file to the total comp count array
            total_comp_count_array0 += comp_count_array

    # # Serial counting        
    # for file in fastq_file_list:
    #     count_array, comp_count_array, n_repeat_errors = count_matches(file)
    #     # Add the counts from each file to the total count array
    #     total_count_array0 += count_array
    #     # Add the comp counts from each file to the total comp count array
    #     total_comp_count_array0 += comp_count_array
    #     total_repeat_errors += n_repeat_errors
    # print(total_repeat_errors)
    print('writing')

    with open(fastq_folder+"/"+save_file+"_0_counts.txt","w") as file:
        for i in range(seq_len - min_len + 1):
            # Set the length of the subsequence to store
            n = i + min_len
            # print(f"subseq len:{n}")
            for j in range(seq_len + 1 - n):
                file.write(str(total_count_array0[i,j])+"\t")
            file.write("\n")

    with open(fastq_folder+"/"+save_file+"_0_comp_counts.txt","w") as file:
        for i in range(seq_len - min_len + 1):
            # Set the length of the subsequence to store
            n = i + min_len
            # print(f"subseq len:{n}")
            for j in range(seq_len + 1 - n):
                file.write(str(total_comp_count_array0[i,j])+"\t")
            file.write("\n")

    # # Count barcode 1
    # barcode  = "GGTCTTGACAAACGTGTGCTTGTAC"
    # barcodec = gen_complement(barcode)

    # with concurrent.futures.ProcessPoolExecutor() as executor:
    #     results = executor.map(count_matches, fastq_file_list)
    #     for result in results:
    #         count_array, comp_count_array, n_repeat_errors = result
    #         # Add the counts from each file to the total count array
    #         total_count_array1 += count_array
    #         # Add the comp counts from each file to the total comp count array
    #         total_comp_count_array1 += comp_count_array

    # with open(fastq_folder+"/"+save_file+"_1_counts.txt","w") as file:
    #     for i in range(seq_len - min_len + 1):
    #         # Set the length of the subsequence to store
    #         n = i + min_len
    #         # print(f"subseq len:{n}")
    #         for j in range(seq_len + 1 - n):
    #             file.write(str(total_count_array1[i,j])+"\t")
    #         file.write("\n")

    # with open(fastq_folder+"/"+save_file+"_1_comp_counts.txt","w") as file:
    #     for i in range(seq_len - min_len + 1):
    #         # Set the length of the subsequence to store
    #         n = i + min_len
    #         # print(f"subseq len:{n}")
    #         for j in range(seq_len + 1 - n):
    #             file.write(str(total_comp_count_array1[i,j])+"\t")
    #         file.write("\n")

    # end = time.perf_counter()
    # print("Time elapsed: {} second(s)".format(end-start))

    # # Parallelize counting of time-filtered strands
    # barcode  = "TATGAGGACGAATCTCCCGCTTATA"
    # barcodec = gen_complement(barcode)
    
    # # Read in the arguments
    # parser = argparse.ArgumentParser(description="Analyzes fastq files for \
    #     perfect sequence matches of different lengths")
    # parser.add_argument("run_folder", type=str, help="The path to the folder \
    #     with the fastq files to analyze.")
    # parser.add_argument("run_length", type=str, help="The length of the run in \
    #     hours.")
    # parser.add_argument("time_step", type=str, help="The time step to sort \
    #     reads by in minutes.")

    # args = parser.parse_args()
    #  # Parse the arguments
    # run_folder = args.run_folder
    # run_length = int(args.run_length)
    # time_step = int(args.time_step)
    
    # # Create lists of fastq files to analyze
    # start_time, end_time, time_step_range, time_range, fastq_pass_file_list, \
    # fastq_fail_file_list, save_folder \
    # = filter_folder_time_step(time_step, run_length, run_folder)
    
    # start = time.perf_counter()
    # # Count barcode 0
    # barcode  = "TATGAGGACGAATCTCCCGCTTATA"
    # barcodec = gen_complement(barcode)
    # count_matches_time_filter_parallel(time_step_range, time_range, \
    # fastq_pass_file_list, fastq_fail_file_list, 0)

    # # Count barcode 1
    # barcode  = "GGTCTTGACAAACGTGTGCTTGTAC"
    # barcodec = gen_complement(barcode)
    # count_matches_time_filter_parallel(time_step_range, time_range, \
    # fastq_pass_file_list, fastq_fail_file_list, 1)
    # end = time.perf_counter()
    # print("Time elapsed: {} second(s)".format(end-start))

    # # # Parallelize the counting of sticky ends
    # # total_sticky_ends = 0
    # # with concurrent.futures.ProcessPoolExecutor() as executor:
    # #     results = executor.map(count_sticky, fastq_file_list)
    # #     for result in results:
    # #         # Add the total number of sticky ends together
    # #         total_sticky_ends += result

    # # print(f"Total number of sticky ends: {total_sticky_ends}")
    
    # # # Parallelize the counting of sequence lengths
    # all_length_list = []
    # count_list = []
    # total_correct_len = 0
    # total_wrong_len = 0
    # # with concurrent.futures.ProcessPoolExecutor() as executor:
    # #     results = executor.map(record_lengths, fastq_file_list)
    # #     for result in results:
    # #         len_list, correct_len, wrong_len = result
    # #         # Organize the length data
    # #         total_correct_len += correct_len
    # #         total_wrong_len += wrong_len
    # #         for length in len_list:
    # #             if length < 1750:
    # #                 if length not in all_length_list:
    # #                     all_length_list.append(length)
    # #                     count_list.append(1)
    # #                 else:
    # #                     count_list[all_length_list.index(length)] += 1
    # #             else:
    # #                 pass
    # # Serial recording of excess lengths
    # recorded_num_repeat_errors = 0
    # for file in fastq_file_list:
    #     len_list, correct_len, wrong_len, num_repeat_errors = record_lengths(file,save=fastq_folder+"/"+save_file)
    #     recorded_num_repeat_errors += num_repeat_errors
    # print(recorded_num_repeat_errors)
    # print(f"Correct length: {total_correct_len}")
    # print(f"Incorrect length: {total_wrong_len}")
    # print(f"Total Sequences: {total_correct_len+total_wrong_len}")

    # end = time.perf_counter()
    # print("Time elapsed: {} second(s)".format(end-start))

    # print(all_length_list)
    # plt.bar(all_length_list, count_list)
    # plt.xlabel('Sequence Length')
    # plt.ylabel('Number of Sequences')
    # plt.title('Read Lengths of Barcode Flongle Sequencing')
    # plt.savefig(fastq_folder+"/read_lengths_hist_no_outliers.png")
    # plt.show()


    # print(test_deco_func.__name__)
    # print(test_deco_func.__doc__)
    # test_deco_func('undecorated', 'barcode lol', 'handle err lol')
    # test_deco_func2 = deco_count_matches(test_deco_func, 'barcode yay', 'handle error ftw', n_repeat=10)
    # print(test_deco_func2.__name__)
    # print(test_deco_func2.__doc__)
    # test_deco_func2('decorated')

    # time_step = 30
    # run_length = 24
    # run_folder= "/Users/emilywu/OneDrive - Massachusetts Institute of Technology/Minion/BarcodeRatio/0-1ratio2-1/20210315_2220_MC-110826_0_AFO272_6d8d594b/time_step_test_folder"
    # save_folder = run_folder + "/time_step_30"

    # # filter_folder_time_step(time_step, run_length, run_folder)

    # test_fastq_file = run_folder + "/fastq_pass/AFO272_pass_1d85e24b_0.fastq"
    # start_time = dt.datetime(2021,3,15,hour=22,minute=20)
    # count_list, comp_count_list, time_list = count_matches_time_step(test_fastq_file)
    # print(count_list)
    # print(comp_count_list)
    # print(time_list)

    # time_step_range = range(0,run_length*60, time_step)
    # make_files_filter_fastq_time(save_folder, time_step_range)
    # pass_save_files = open_files_filter_fastq_time(save_folder, "pass", time_step_range)
    

    # time_range = []
    # for i in range(len(time_step_range)):
    #     delta_time = dt.timedelta(minutes=time_step_range[i])
    #     time_range.append(start_time+delta_time)
    
    # filter_fastq_time_step(test_fastq_file, start_time, time_range, time_step_range, pass_save_files)
    # close_files_filter_fastq_time(pass_save_files)
    # test_run_time = dt.datetime(2021,3,16,hour=23,minute=50)
    # test_dt = test_run_time-start_time
    # print(test_dt/dt.timedelta(minutes=1))
    # print(time_range)
    # print(len(time_range))
import concurrent.futures
import numpy as np
import argparse
from pathlib import Path
import os
import time
import datetime as dt
import matplotlib.pyplot as plt 
import beepy as bp
from PertQuant.analyzeFASTQ.countmatches import make_save_folder
from PertQuant.analyzeFASTQ.sum_count_dat import str_2_date_time
from PertQuant.analyzeFASTQ.sum_count_dat import get_time_range

def decide_time_bin(read_time, start_time, time_delta, range_len):
    time_diff = read_time - start_time
    time_floor = np.floor(time_diff/time_delta)
    time_mod = time_diff % time_delta

    # Figure out which bin the time goes into
    if time_mod > dt.timedelta(0):
        time_bin = int(time_floor)
    else:
        time_bin = int(time_floor) - 1

    # Make sure the time bin is in the bin range
    if time_bin < 0:
        time_bin = 0
    elif time_bin >= range_len:
        time_bin = range_len - 1
    return(time_bin)

def get_time(index):
    """
    get_time(index)

    Takes in an index of the line containing a sequence ID of a sequence to be
    analyzed. It uses this to identify the read time of the associated sequence
    and returns this as a datetime object. It only takes in one argument in order
    to be easily parallelized. It also relies on a set of global arguments.

    Arguments:
        index (int): The index of the line containing the sequence ID data.

    Global Arguments;
        lines (list): A list of lines read in from a FASTQ file.

    Output:
        time_bin (datetime): the bin of 
    """
    global lines 
    global start_time 
    global time_step_td
    global range_len
    # Read in the sequence ID
    seq_ID = lines[index].rstrip("\n")
    # Split the sequence ID into parts
    seq_ID_list = seq_ID.split()
    # Get the read time string
    read_time_str = seq_ID_list[4]
    # Convert the string to a date time object
    read_time = str_2_date_time(read_time_str)

    time_bin = decide_time_bin(read_time, start_time, time_step_td, range_len)
    return(time_bin)    

if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="Analyzes fastq files for \
        perfect sequence matches of different lengths")
    parser.add_argument("fastq_folder", type=str, help="The path to the folder \
        with the fastq files to analyze.")
    parser.add_argument("save_file", type=str, help="The name extension to add to \
        save file names that results are saved to.")
    parser.add_argument("--time", type=int, help="Bins counts by the \
        provided time step in minutes. Default is to not bin by time.")
    parser.add_argument("--run", type=int, help="The length of the sequencing run \
        in hours. Defaults to 24 hours.")
    parser.add_argument("--prog", type=bool, help="If True, prints progress \
        messages. Defaults to False.")
    parser.add_argument("--pf", type=bool, help="If True, finds the pass and fail \
        folders and analyzes fastq files in both. Defaults to False.")
    parser.add_argument("--beep", type=int, help="Plays a sound using beepy when \
        the program finishes running. To pick a sound, provide an integer from \
        1-7. To not play a sound, set to 0. Defaults to 1.")
    args = parser.parse_args()

    # Parse the arguments
    fastq_folder = args.fastq_folder
    save_file_name = args.save_file
    if args.prog:
        prog = args.prog
    else:
        prog = False
    if args.pf:
        pf = args.pf 
    else:
        pf = False
    if args.time:
        time_step = int(args.time)
    else:
        time_step = 10
    if args.run:
        run_length = args.run
    else:
        run_length = 24 # Default sequencing run length is 24 hours
    if args.beep:
        which_beep = args.beep 
    else:
        which_beep = 1

    # Create the save folder
    save_folder = make_save_folder(fastq_folder)

    # Get the time range info
    start_time, end_time, time_step_range, time_range = \
    get_time_range(time_step, run_length, save_folder)
    time_step_td = dt.timedelta(minutes=time_step)
    range_len = len(time_range)

    # Create list of fastq files to analyze
    if pf:
        fastq_pass_files = Path(f"{fastq_folder}/fastq_pass").glob("*.fastq")
        fastq_pass_file_list = [str(file) for file in fastq_pass_files]
        fastq_fail_files = Path(f"{fastq_folder}/fastq_fail").glob("*.fastq")
        fastq_fail_file_list = [str(file) for file in fastq_fail_files]
        fastq_file_list = fastq_pass_file_list + fastq_fail_file_list
    else:    
        fastq_files = Path(fastq_folder).rglob("*.fastq")
        fastq_file_list = [str(file) for file in fastq_files] 

    # Make a list to contain sequence lengths
    read_counts_arr = np.zeros(len(time_range))

    # Analyze all the fastq files
    index = 0
    start = time.perf_counter()
    for read_file_path in fastq_file_list:
        # Read in the data from a fastq file
        with open(read_file_path, 'r') as read_file:
            lines = read_file.readlines()
            # index_list = range(int(len(lines)/4))
            index_arr = np.arange(0,len(lines),4)
            # Parallelize the subsequence counting
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = executor.map(get_time, index_arr)
                for result in results:
                    read_counts_arr[result] += 1
        if prog:
            index += 1
            print(f"Analyzed {index}/{len(fastq_file_list)} files.")
        else:
            pass

    # Make save file
    with open(f"{save_folder}/{save_file_name}_counts_tstep_{time_step}.txt", "w") as save_file:
        for i in read_counts_arr:
            save_file.write(f"{i:.0f} ")
    end = time.perf_counter()
    print("Time elapsed: {} second(s)".format(end-start))

    if which_beep > 0:
        bp.beep(sound=which_beep)
    else:
        pass

    # # Test decide_time_bin
    # start_time =  dt.datetime(2021,3,10)
    # time_step_range = range(0,24*60, 30)
    # time_step = dt.timedelta(minutes=30)
    # time_range = []
    # for i in range(len(time_step_range)):
    #     delta_time = dt.timedelta(minutes=time_step_range[i])
    #     time_range.append(start_time+delta_time)
    # range_len = len(time_range)

    # test_time_range = []
    # test_time_step_range = range(-15,25*60, 15)
    # for i in range(len(test_time_step_range)):
    #     delta_time = dt.timedelta(minutes=test_time_step_range[i])
    #     test_time_range.append(start_time+delta_time)

    # for read_time in test_time_range:
    #     time_bin = decide_time_bin(read_time, start_time, time_step, range_len)
    #     print(f"read_time: {read_time}, bin: {time_bin}, bin_lower_bound: {time_range[time_bin]}")
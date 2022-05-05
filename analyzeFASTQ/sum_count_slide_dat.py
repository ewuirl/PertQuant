import argparse
from pathlib import Path
import numpy as np
import concurrent.futures
# import os
import datetime as dt
import time
import beepy as bp
from PertQuant.analyzeFASTQ.sum_count_dat import get_count_settings
from PertQuant.analyzeFASTQ.sum_count_dat import get_time_range
from PertQuant.analyzeFASTQ.sum_count_dat import decide_time_bin

def stripped_str_2_date_time(read_time_str):
    """
    str_2_date_time(read_time_str)

    This function takes a read time string and creates a date time object 
    containing the information.

    Arguments:
        read_time_str (str): A string describing the time the sequence was read.

    Returns: 
        read_time (datetime): A datetime object representing the time the 
            sequence was read..
    """
    read_time = dt.datetime(int(read_time_str[:4]), int(read_time_str[5:7]), \
                int(read_time_str[8:10]), hour=int(read_time_str[11:13]), \
                minute=int(read_time_str[14:16]), second=int(read_time_str[17:19]))
    return read_time 

def read_slide_dat_file_time(dat_file_path):
    """
    This function takes a path to a dat file of counts and bins each read's 
    counts by its read time. It sums the counts in each bin, and returns an array
    of these summed, binned counts. Each column in the array represents the 
    counts for a different subsequence, and each row represents the counts for a
    different bin. Each time bin is upper bound inclusive and lower bound 
    exclusive, except the first and last bins. To ensure all reads are binned, 
    the first bin is inclusive of read times equal to or before the start time, 
    and the last bin is inclusive of read times equal to or after the end of the 
    run.

    Arguments:
        dat_file_path (str): The path to the dat file to get base-averaged Q 
            scores for.

    Global Arguments:
        n_targets (int): The number of target sequences in the run.
        target_lengths (list): A list of the target lengths (int).
        n_barcodes (int): The number of barcode sequences.
        min_len (int): The minimum subsequence length that was counted.
        barcoded (bool): If True, the read data is barcoded, and will be sorted
            according to barcode as well. Specified by command-line argument. 
            Defaults to False.
        prog (bool): If True, progress statements are printed. Specified by     
            command-line argument. Defaults to False.
        start_time (datetime): A datetime object of the start of the sequencing
            run.
        time_range_len (int): The number of time bins.
        time_step_td (timedelta): A timedelta object of the timestep of the bins.

    Returns:
        target_sum_array_list (list): A list of arrays of subsequence counts 
            binned and summed by read time. Each target gets its own array. 
        targetc_sum_array_list (list): A list of arrays of complement subsequence 
            counts binned and summed by read time. Each target gets its own array. 
    """
    global n_targets
    global target_lengths
    global n_barcodes
    global min_len
    global barcoded 
    global prog
    global start_time
    global time_range_len
    global time_step_td

    with open(dat_file_path, 'r') as dat_file:
        header = dat_file.readline()
        # Set up the lists of summed arrays
        
        target_comp_sum_array = np.zeros((len(time_range), int(n_targets*2)),dtype=int)
        while True:
            # Read the sequence ID
            line = dat_file.readline()
            if not line:
                break
            line_list = line.rstrip("\n").split("\t")
            # Get read time
            read_time_str = line_list[3]
            read_time = stripped_str_2_date_time(read_time_str)
            # Decide time bin
            current_bin = decide_time_bin(read_time, start_time, time_step_td, \
                time_range_len)

            target_comp_sum_array[current_bin,:] += np.array(line_list[8:],dtype=int)

        if prog:
            dat_file_name = dat_file_path.split("/")[-1]
            dat_file_name_list = dat_file_name.split("_")
            print(f"Finished counting {dat_file_name_list[1]} file {dat_file_name_list[3]} count data.")
        else:
            pass

        return target_comp_sum_array

def read_slide_dat_file_all(dat_file_path):
    """
    read_dat_file_all(dat_file_path)

    This function reads in all the subsequence counts from a dat file and
    sums the subsequence counts of each read together.

    Arguments:
        dat_file_path (str): The path to the dat file to read subsequence
            count data from.
    
    Global Arguments:
        n_targets (int): The number of target sequences in the run.
        target_lengths (list): A list of the target lengths (int).
        n_barcodes (int): The number of barcode sequences.
        min_len (int): The minimum subsequence length that was counted.
        barcoded (bool): If True, the read data is barcoded, and will be sorted
            according to barcode as well. Specified by command-line argument. 
            Defaults to False.
        prog (bool): If True, progress statements are printed. Specified by     
            command-line argument. Defaults to False.

    Returns:
        target_sum_array_list (list): A list of arrays of summed subsequence 
            counts. Each target gets its own array. 
        targetc_sum_array_list (list): A list of arrays of summed complement 
            subsequence counts. Each target gets its own array. 
    """
    global n_targets
    global n_barcodes
    global min_len
    global barcoded 
    global prog

    with open(dat_file_path, 'r') as dat_file:
        header = dat_file.readline()
        # Set up the lists of summed arrays
    
        target_comp_sum_array = np.zeros(int(n_targets*2),dtype=int)

        while True:
            # Read the sequence ID
            line = dat_file.readline()
            if not line:
                break
            line_list = line.rstrip("\n").split("\t")
            # print(np.array(line_list[7:],dtype=int))
            target_comp_sum_array += np.array(line_list[8:],dtype=int)

        if prog:
            dat_file_name = dat_file_path.split("/")[-1]
            dat_file_name_list = dat_file_name.split("_")
            print(f"Finished counting {dat_file_name_list[1]} file {dat_file_name_list[3]} count data.")
        else:
            pass
        return target_comp_sum_array

def read_slide_dat_file_bin_parallel(bin_range, bin_func, dat_file_list):
    """
    read_dat_file_bin_parallel(bin_range, bin_func, dat_file_list, target_lengths)

    This function takes a list of dat files, and a list of target lengths, and 
    it parallelizes the summing and binning using the provided bin range and 
    read binning function of the subsequence counts of the reads in all of the 
    files. The counts are also summed across files.
    
    Arguments:
        bin_range (arr): An array of the lower bounds of the bins to bin the
            subsequence counts with.
        bin_func (func): The function to use to bin the subsequence counts with,
            such as read_dat_file_Pbin and read_dat_file_time.
        dat_file_list (list): A list of paths to dat file of subsequence counts 
            to bin and sum.
        target_lengths (list): A list of the target lengths (int).
    
    Returns:
        total_target_sum_array_list (list): A list of arrays of subsequence counts 
            binned and summed and summed using bin_func. Results of each file 
            are also summed. Each target gets its own array. 
        total_targetc_sum_array_list (list): A list of arrays of complement 
            subsequence counts binned and summed using bin_func. Results of each
            file are also summed. Each target gets its own array. 
    """
    total_target_comp_sum_array = np.zeros((len(bin_range),int(2*n_targets)), dtype=int)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(bin_func, dat_file_list)
        for result in results:
            # Unpack the subsequence counts
            total_target_comp_sum_array += result
    return total_target_comp_sum_array



def read_slide_dat_file_all_parallel(dat_file_list, n_targets):
    """
    read_dat_file_all_parallel(dat_file_list, n_targets)

    This function takes a list of dat file paths and sums all the subsequence
    counts together by read and file.

    Arguments:
        dat_file_list (list): A list of paths to dat file of subsequence counts 
            to bin and sum.
        target_lengths (list): A list of the target lengths (int).

    Returns:
        total_target_sum_array_list (list): A list of arrays of summed subsequence 
            counts. Results of each file are also summed. Each target gets its 
            own array. 
        total_targetc_sum_array_list (list): A list of arrays of summed 
            complement subsequence counts. Results of each file are also summed. 
            Each target gets its own array. 
    """
    total_target_comp_sum_array = np.zeros(int(2*n_targets), dtype=int)

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(read_slide_dat_file_all, dat_file_list)
        for result in results:
            # Unpack the subsequence counts
            total_target_comp_sum_array += result
    return total_target_comp_sum_array

def write_all_summed_counts_array(save_file_path, summed_counts_array):
    """
    write_all_summed_counts_array(save_file_path, summed_counts_array)

    This function writes the data from a summed count array (1D) into a file. 
    
    Arguments:
        save_file_path (str): The path to the file to save the data to.
        summed_counts_array (arr): A 1D array containing summed subsequence count
            data.

    Returns:
        Nothing. Writes the data in the provided array to the file specified by
            the path.
    """
    with open(save_file_path, 'w') as save_file:
        for i in range(int(len(summed_counts_array)/2)):
            save_file.write(f"{i}\t{i}C\t")
        save_file.write(f"\n")
        for i in range(len(summed_counts_array)):
            save_file.write(f"{summed_counts_array[i]}\t")

def write_2d_summed_counts_array(save_file_path, summed_counts_arrays, max_stretch=0):
    """
    write_all_summed_counts_array(save_file_path, summed_counts_array)

    This function writes the data from a summed count array (1D) into a file. 
    
    Arguments:
        save_file_path (str): The path to the file to save the data to.
        summed_counts_arrays (tuple): 

    Returns:
        Nothing. Writes the data in the provided array to the file specified by
            the path.
    """
    global n_targets
    arr_shape = summed_counts_arrays.shape
    with open(save_file_path, 'w') as save_file:
        if max_stretch > 0:
            num_targets = int(arr_shape[1]/(max_stretch*2))
            for i in range(max_stretch):
                for j in range(num_targets):
                    save_file.write(f"{j}_{i+1}\t{j}C_{i+1}\t")
        else:
            for i in range(int(arr_shape[1]/2)):
                save_file.write(f"{i}\t{i}C\t")
        save_file.write(f"\n")
        for i in range(arr_shape[0]):
            for j in range(arr_shape[1]):
                save_file.write(f"{summed_counts_arrays[i][j]}\t")
            if i < arr_shape[0] - 1:
                save_file.write("\n")
            else:
                pass


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="Analyzes dat files of subsequence \
        matches. Optional arguments are provided to bin the summed data. \
        Otherwise, the counts are summed together without binning.")
    parser.add_argument("dat_folder", type=str, help="The path to the folder \
        with the dat files to analyze.")
    parser.add_argument("save_file", type=str, help="The name extension to add to \
        save file names that results are saved to.")
    parser.add_argument("--settings", type=str, help="The path to the file \
        containing count analysis settings.")
    parser.add_argument("--barcoded", type=str, help="If True, the program will \
        sort results by the barcode.")
    parser.add_argument("--time", type=int, help="Bins counts by the \
        provided time step in minutes. Default is to not bin by time.")
    parser.add_argument("--run", type=int, help="The length of the sequencing run \
        in hours. Defaults to 24 hours.")
    parser.add_argument("--prog", type=str, help="If True, prints progress \
        messages.")
    parser.add_argument("--pf", type=str, help="If True, sums pass and fail files \
        separately.")
    parser.add_argument("--max_stretch", type=str, help="If True, sums pass and fail files \
        separately.")
    # parser.add_argument("--Lhist", type=bool, help="If True, makes a histogram \
    #     of the average Q scores of the sequences.")
    parser.add_argument("--beep", type=int, help="Plays a sound using beepy when \
        the program finishes running. To pick a sound, provide an integer from \
        1-7. To not play a sound, set to 0. Defaults to 1.")
    args = parser.parse_args()

    # Parse the arguments
    dat_folder = args.dat_folder
    save_file_name = args.save_file

    if args.settings:
        count_settings_path = args.settings
    else:
        count_settings_path  = f"{dat_folder}/count_settings_{save_file_name}.txt"

    # if args.sum:
    #     do_sum = args.sum
    # else:
    #     do_sum = False

    if args.barcoded:
        barcoded = eval(args.barcoded)
    else:
        barcoded = False

    if args.time:
        time_step = args.time
    else:
        time_step = 0

    if args.run:
        run_length = args.run
    else:
        run_length = 24 # Default sequencing run length is 24 hours

    if args.prog:
        prog = eval(args.prog)
    else:
        prog = False

    if args.beep:
        which_beep = args.beep 
    else:
        which_beep = 1

    if args.pf:
        pf = eval(args.pf)
    else:
        pf = False

    if args.max_stretch:
        max_stretch = int(args.max_stretch)
    else:
        max_stretch = 0

    # Get all the dat files
    dat_files = Path(dat_folder).rglob("*.dat")
    dat_file_list = [str(file) for file in dat_files]

    # Get the number of targets
    n_targets, target_list, target_lengths, n_barcodes, settings_dict = \
    get_count_settings(count_settings_path)

    if max_stretch > 0:
        n_targets = n_targets*max_stretch
    else:
        pass

    if time_step > 0.0:
        # Get the time range info
        start_time, end_time, time_step_range, time_range = get_time_range(time_step, run_length, dat_folder)
        time_step_td = dt.timedelta(minutes=time_step)
        time_range_len = len(time_range)

        start = time.perf_counter()
        # Sum and bin the counts according to average P(corr)
        binned_target_comp_array = \
        read_slide_dat_file_bin_parallel(time_range, read_slide_dat_file_time, dat_file_list)

        save_file_path = f"{dat_folder}/{save_file_name}_tstep-{time_step}_target_comp_slide_counts.txt"
        # Write the summedd counts to save files
        write_2d_summed_counts_array(save_file_path, binned_target_comp_array, max_stretch=max_stretch)
    else:
        save_file_path = f"{dat_folder}/{save_file_name}_all_target_comp_slide_counts.txt"
        start = time.perf_counter()
        if pf:
            dat_file_arr = np.array(dat_file_list)
            pass_file_arr = dat_file_arr[["pass" in file_name for file_name in dat_file_arr]]
            fail_file_arr = dat_file_arr[["fail" in file_name for file_name in dat_file_arr]]
            pass_target_comp_sum_array = read_slide_dat_file_all_parallel(pass_file_arr, n_targets)
            fail_target_comp_sum_array = read_slide_dat_file_all_parallel(fail_file_arr, n_targets)
            total_target_comp_sum_array = pass_target_comp_sum_array + fail_target_comp_sum_array
            all_summed_counts_array = np.vstack((pass_target_comp_sum_array, fail_target_comp_sum_array, \
            total_target_comp_sum_array))
            write_2d_summed_counts_array(save_file_path, all_summed_counts_array, max_stretch=max_stretch)
        else:
            total_target_comp_sum_array = read_slide_dat_file_all_parallel(dat_file_list, n_targets)
            write_all_summed_counts_array(save_file_path, total_target_comp_sum_array)

    end = time.perf_counter()
    print("Time elapsed: {} second(s)".format(end-start))
    
    if which_beep > 0:
        bp.beep(sound=which_beep)
    else:
        pass
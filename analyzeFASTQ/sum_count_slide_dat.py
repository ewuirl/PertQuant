import argparse
from pathlib import Path
import numpy as np
import concurrent.futures
# import os
import time
import beepy as bp
from PertQuant.analyzeFASTQ.sum_count_dat import get_count_settings

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
            target_comp_sum_array += np.array(line_list[7:],dtype=int)

        if prog:
            dat_file_name = dat_file_path.split("/")[-1]
            dat_file_name_list = dat_file_name.split("_")
            print(f"Finished counting {dat_file_name_list[1]} file {dat_file_name_list[3]} count data.")
        else:
            pass
        return target_comp_sum_array

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

def write_all_summed_counts_array_pf(save_file_path, summed_counts_arrays):
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
        for i in range(int(len(summed_counts_arrays[0])/2)):
            save_file.write(f"{i}\t{i}C\t")
        save_file.write(f"\n")
        for i in range(len(summed_counts_arrays)):
            for j in range(len(summed_counts_arrays[i])):
                save_file.write(f"{summed_counts_arrays[i][j]}\t")
            if i < len(summed_counts_arrays) - 1:
                save_file.write("\n")

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

    # Get all the dat files
    dat_files = Path(dat_folder).rglob("*.dat")
    dat_file_list = [str(file) for file in dat_files]




    # Get the number of targets
    n_targets, target_list, target_lengths, n_barcodes, settings_dict = \
    get_count_settings(count_settings_path)
    save_file_path = f"{dat_folder}/{save_file_name}_all_target_comp_slide_counts.txt"

    start = time.perf_counter()
    if pf:
        dat_file_arr = np.array(dat_file_list)
        pass_file_arr = dat_file_arr[["pass" in file_name for file_name in dat_file_arr]]
        fail_file_arr = dat_file_arr[["fail" in file_name for file_name in dat_file_arr]]
        pass_target_comp_sum_array = read_slide_dat_file_all_parallel(pass_file_arr, n_targets)
        fail_target_comp_sum_array = read_slide_dat_file_all_parallel(fail_file_arr, n_targets)
        total_target_comp_sum_array = pass_target_comp_sum_array + fail_target_comp_sum_array
        summed_counts_arrays = [pass_target_comp_sum_array, fail_target_comp_sum_array, \
        total_target_comp_sum_array]
        write_all_summed_counts_array_pf(save_file_path, summed_counts_arrays)
    else:
        total_target_comp_sum_array = read_slide_dat_file_all_parallel(dat_file_list, n_targets)
        write_all_summed_counts_array(save_file_path, total_target_comp_sum_array)

    end = time.perf_counter()
    print("Time elapsed: {} second(s)".format(end-start))
    
    if which_beep > 0:
        bp.beep(sound=which_beep)
    else:
        pass
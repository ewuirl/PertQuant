# import sys
# sys.path.append('../')
from PertQuant.simCRN.ml_nupack import gen_complement
import argparse
from pathlib import Path
import numpy as np
import concurrent.futures
import datetime as dt
import matplotlib.pyplot as plt 
import os
import time
import beepy as bp

def get_count_settings(count_settings_path):
    # Open the settings file to find the target file path
    with open(count_settings_path, 'r') as settings_file:
        lines = settings_file.readlines()
        # Get number of targets and target lengths
        line = lines[0].rstrip("\n")
        line_list = line.split()
        n_targets = int(line_list[3])
        target_list = []
        target_lengths = []
        for i in range(n_targets):
            target = lines[i+1].rstrip("\n").split()[1]
            target_list.append(target)
            target_lengths.append(len(target))

        # Get number of barcodes
        line = lines[n_targets+1].rstrip("\n")
        line_list = line.split()
        n_barcodes = int(line_list[3])

        # Get min_len
        line = lines[n_targets+3+n_barcodes].rstrip("\n")
        line_list = line.split()
        min_len = int(line_list[2])

        # Get handle_repeat_error
        line = lines[n_targets+4+n_barcodes].rstrip("\n")
        line_list = line.split()
        handle_repeat_error = bool(line_list[2])

        # Get repeat_list
        line = lines[n_targets+5+n_barcodes].rstrip("\n")
        line_list = line.split()
        if line_list[2] == "[]":
            repeat_list = []
        else:
            repeat_list = line_list[2:]

        # Get n_repeat
        line = lines[n_targets+6+n_barcodes].rstrip("\n")
        line_list = line.split()
        n_repeat = int(line_list[2])

    return (n_targets, target_list, target_lengths, n_barcodes, min_len, \
        handle_repeat_error, repeat_list, n_repeat)

def get_seq_info(seq_ID, features):
    seq_ID_list = seq_ID.split()
    read_time_str = seq_ID_list[4]
    feature_list = features.split()
    avg_Q_score = float(feature_list[0])
    read_len = int(feature_list[1])
    barcode_ID = int(feature_list[2])
    has_repeat_error = int(feature_list[3])
    return(read_time_str, avg_Q_score, read_len, barcode_ID, has_repeat_error)

def str_2_date_time(read_time_str):
    read_time = dt.datetime(int(read_time_str[11:15]), int(read_time_str[16:18]), \
                int(read_time_str[19:21]), hour=int(read_time_str[22:24]), \
                minute=int(read_time_str[25:27]), second=int(read_time_str[28:30]))
    return read_time 

def make_count_save_folder(dat_folder, save_file, time_step, Qbin, Pbin):
    if time_step > 0.0:
        time_step_name = f"_tstep-{time_step}"
    else:
        time_step_name = ""

    if Qbin > 0.0:
        Qbin_name = f"_Qbin-{Qbin:.2f}".replace(".", "-")
    else:
        Qbin_name = ""

    if Pbin > 0.0:
        Pbin_name = f"_Pbin-{Pbin:.2f}".replace(".", "-")
    else:
        Pbin_name = ""

    # Put together the save folder
    save_folder = f"{dat_folder}/summed_counts_{save_file_name}{time_step_name}{Qbin_name}{Pbin_name}"

    # Make the folder if it doesn't already exist
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)

    return save_folder

def get_Q_scores(dat_file_path):
    global n_targets

    with open(dat_file_path, 'r') as dat_file:
        # Read in the lines
        lines = dat_file.readlines()
        index_arr = np.arange(1, len(lines)-1, 2+2*n_targets)

        # Make an array to store all the Q scores
        avg_Q_score_arr = np.zeros(len(index_arr))

        # Read the data and save it to the arrays
        for index in index_arr:
            seq_ID = lines[index].rstrip('\n')
            features = lines[index+1].rstrip('\n')
            read_time_str, avg_Q_score, read_len, barcode_ID, has_repeat_error =\
            get_seq_info(seq_ID, features)
            avg_Q_score_arr[int((index-1)/(2+2*n_targets))] = avg_Q_score
        
        # Return the array
        return avg_Q_score_arr

def get_Q_scores_parallel(dat_file_list):
    all_Q_score_arr = np.array([])
    # Parallelize the Q score data collection
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(get_Q_scores, dat_file_list)
        for result in results:
            # Unpack the subsequence counts
            avg_Q_score_arr = result
            all_Q_score_arr = np.concatenate((all_Q_score_arr,avg_Q_score_arr))
    return all_Q_score_arr

def read_count_line(dat_file):
    count_line = dat_file.readline().rstrip('\n')
    count_line_list = [int(i) for i in count_line.split()]
    count_arr = np.array(count_line_list)
    return count_arr

def read_dat_file_Pbin(dat_file_path):
    global n_targets
    global target_lengths
    global n_barcodes
    global min_len
    global barcoded 
    global prog
    global Pbin
    global P_corr_bins

    with open(dat_file_path, 'r') as dat_file:
        header = dat_file.readline()
        # Set up the lists of summed arrays
    
        target_sum_array_list = []
        targetc_sum_array_list = []
        for i in range(n_targets):
            target_len = target_lengths[i]
            # Figure out the maximum nunber of subsequences possible
            max_num_subseqs = target_len - min_len + 1
            # Calculate the count array length
            array_len = int(max_num_subseqs * (max_num_subseqs+1)/2)
            target_sum_array_list.append(np.zeros((len(P_corr_bins),array_len), dtype=int))
            targetc_sum_array_list.append(np.zeros((len(P_corr_bins),array_len), dtype=int))
        else:
            pass

        while True:
            # Read the sequence ID
            seq_ID = dat_file.readline()
            if not seq_ID:
                break
            # Get sequence info
            seq_ID = seq_ID.rstrip('\n')
            features = dat_file.readline().rstrip('\n')
            read_time_str, avg_Q_score, read_len, barcode_ID, has_repeat_error =\
            get_seq_info(seq_ID, features)

            # Calculate P(corr)
            seq_P_corr = 1.0 - 10.0 ** (-avg_Q_score/10.0)

            # Figure out which bin to use
            P_corr_bin = 0
            # Check if the next bin up should be used and make sure the last bin
            # Catches everything greater than the lower bound of the last bin
            while seq_P_corr >= P_corr_bins[P_corr_bin] + Pbin \
            and P_corr_bin + 1 < len(P_corr_bins):
                P_corr_bin += 1

            # Get count info
            for i in range(n_targets):
                # Get target count info
                target_sum_array_list[i][P_corr_bin,:] += read_count_line(dat_file)
                # Get complement count info
                targetc_sum_array_list[i][P_corr_bin,:] += read_count_line(dat_file)
        if prog:
            dat_file_name = dat_file_path.split("/")[-1]
            dat_file_name_list = dat_file_name.split("_")
            print(f"Finished counting {dat_file_name_list[1]} file {dat_file_name_list[3]} count data.")
        else:
            pass
        return(target_sum_array_list, targetc_sum_array_list)

def read_dat_file_bin_parallel(bin_range, bin_func, dat_file_list, target_lengths):
    total_target_sum_array_list = []
    total_targetc_sum_array_list = []
    for i in range(n_targets):
        target_len = target_lengths[i]
        # Figure out the maximum nunber of subsequences possible
        max_num_subseqs = target_len - min_len + 1
        # Calculate the count array length
        array_len = int(max_num_subseqs * (max_num_subseqs+1)/2)
        total_target_sum_array_list.append(np.zeros((len(bin_range),array_len), dtype=int))
        total_targetc_sum_array_list.append(np.zeros((len(bin_range),array_len), dtype=int))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(bin_func, dat_file_list)
        for result in results:
            # Unpack the subsequence counts
            target_sum_array_list, targetc_sum_array_list = result
            
            for i in range(len(total_target_sum_array_list)):
                total_target_sum_array_list[i] += target_sum_array_list[i]
                total_targetc_sum_array_list[i] += targetc_sum_array_list[i]

    return (total_target_sum_array_list, total_targetc_sum_array_list)

def read_dat_file_time(dat_file_path):
    global n_targets
    global target_lengths
    global n_barcodes
    global min_len
    global barcoded 
    global prog
    global time_range
    global time_step_td

    with open(dat_file_path, 'r') as dat_file:
        header = dat_file.readline()
        # Set up the lists of summed arrays
    
        target_sum_array_list = []
        targetc_sum_array_list = []
        for i in range(n_targets):
            target_len = target_lengths[i]
            # Figure out the maximum nunber of subsequences possible
            max_num_subseqs = target_len - min_len + 1
            # Calculate the count array length
            array_len = int(max_num_subseqs * (max_num_subseqs+1)/2)
            target_sum_array_list.append(np.zeros((len(time_range),array_len), dtype=int))
            targetc_sum_array_list.append(np.zeros((len(time_range),array_len), dtype=int))
        else:
            pass

        while True:
            # Read the sequence ID
            seq_ID = dat_file.readline()
            if not seq_ID:
                break
            # Get sequence info
            seq_ID = seq_ID.rstrip('\n')
            features = dat_file.readline().rstrip('\n')
            read_time_str, avg_Q_score, read_len, barcode_ID, has_repeat_error =\
            get_seq_info(seq_ID, features)
            read_time = str_2_date_time(read_time_str)

            # Figure out which bin to use
            current_bin = 0

            # Check if the next bin up should be used and make sure the last bin
            # Catches everything greater than the lower bound of the last bin
            while read_time > time_range[current_bin] + time_step_td \
            and current_bin + 1 < len(time_range):
                current_bin += 1

            # Get count info
            for i in range(n_targets):
                # Get target count info
                target_sum_array_list[i][current_bin,:] += read_count_line(dat_file)
                # Get complement count info
                targetc_sum_array_list[i][current_bin,:] += read_count_line(dat_file)
        
        if prog:
            dat_file_name = dat_file_path.split("/")[-1]
            dat_file_name_list = dat_file_name.split("_")
            print(f"Finished counting {dat_file_name_list[1]} file {dat_file_name_list[3]} count data.")
        else:
            pass
        return(target_sum_array_list, targetc_sum_array_list)

def write_summed_counts_array(save_file_path, summed_counts_array):
    with open(save_file_path, 'w') as save_file:
        array_dim = summed_counts_array.shape 
        for i in range(array_dim[0]):
            for j in range(array_dim[1]):
                save_file.write(f"{summed_counts_array[i,j]}\t")
            save_file.write("\n")


def write_summed_counts(save_folder, save_file_name, Pbin, Qbin, time_step, \
    total_target_sum_array_list, total_targetc_sum_array_list, prog):
    for i in range(len(total_targetc_sum_array_list)):
        # Make the paths for the save files
        if Pbin > 0.0:
            Pbin_name = f"_Pbin-{Pbin:.2f}".replace(".", "-")
        else:
            Pbin_name = ""
        if Qbin > 0.0:
            Qbin_name = f"_Qbin-{Qbin:.2f}".replace(".", "-")
        else:
            Qbin_name = ""
        if time_step > 0.0:
            tstep_name = f"_tstep-{time_step}"
        else:
            tstep_name = ""

        target_save_file_path = f"{save_folder}/{save_file_name}{tstep_name}{Qbin_name}{Pbin_name}_{i}_counts.txt"
        targetc_save_file_path = f"{save_folder}/{save_file_name}{tstep_name}{Qbin_name}{Pbin_name}_{i}_comp_counts.txt"
        # Write the counts to the save files
        write_summed_counts_array(target_save_file_path, total_target_sum_array_list[i])
        write_summed_counts_array(targetc_save_file_path, total_targetc_sum_array_list[i])
        if prog:
            print(f"Finished writing target {i} count data.")
        else:
            pass

def get_sum_settings(sum_count_folder, has_time_step=False, has_Pbin=False, has_Qbin=False):
    sum_count_list = sum_count_folder.split("_")
    # Set default values
    time_step = 0
    Pbin = 0.0
    Qbin = 0.0

    assert has_Pbin==False or has_Qbin==False, "Pbin and Qbin cannot both be > 0"

    # Figure out what the sum settings were
    # Figure out the time step if it was used
    if has_time_step and has_Pbin:
        time_step_str = sum_count_list[-2]
        time_step = int(time_step_str.replace("tstep-",""))
    elif has_time_step and has_Qbin:
        time_step_str = sum_count_list[-2]
        time_step = int(time_step_str.replace("tstep-",""))    
    elif has_time_step and has_Pbin == False and has_Qbin == False:
        time_step_str = sum_count_list[-1]
        time_step = int(time_step_str.replace("tstep-",""))
    else:
        pass
    # Handle Q score or P(correct) binning if they were used
    if has_Pbin:
        Pbin_str = sum_count_list[-1]
        Pbin = float(Pbin_str.replace("Pbin-","").replace("-","."))
    elif has_Qbin:
        Qbin_str = sum_count_list[-1]
        Qbin = float(Qbin_str.replace("Qbin-","").replace("-","."))
    else:
        pass

    return (time_step, Pbin, Qbin)

def get_time_range(time_step, run_length, sum_count_folder):
    """
    get_time_range(time_step, run_length, sum_count_folder)

    Takes in a time step in minutes, and a run length in hours, and sum_count_folder,
    the path to a folder containing dat files. It figures out the start and end 
    times of the run, and generates a time step range and list of corresponding 
    timedeltas. 

    Arguments:
        time_step (int/float): a time step in minutes to bin the data into.
        run_length (int/float): the length of the sequencing run in hours.
        sum_count_folder (str): a string representing the path to a folder 
            containing dat files.

    Returns:
        start_time (datetime): a datetime object of the start time of the run.
        end_time (datetime): a datetime object of the end time of the run.
        time_step_range (range): a range object with bounds [time_step, run_length (min)]
        time_range (list): a list of timedelta objects corresponding to the 
            values in time_step_range
        

    """
    sum_count_folder_list = sum_count_folder.split("/")
    folder_name = sum_count_folder_list[-2] # change to -2
    folder_name_list = folder_name.split("_")
    date = folder_name_list[0]
    time = folder_name_list[1]
    start_time = dt.datetime(int(date[:4]), int(date[4:6]), \
        int(date[6:8]), hour=int(time[:2]), minute=int(time[2:]))

    # Prepare lists of time steps
    time_step_range = range(0,run_length*60, time_step)
    time_range = []
    for i in range(len(time_step_range)):
        delta_time = dt.timedelta(minutes=time_step_range[i])
        time_range.append(start_time+delta_time)

    end_time = time_range[-1]
    return(start_time, end_time, time_step_range, time_range)

def read_dat_file_all(dat_file_path):
    global n_targets
    global target_lengths
    global n_barcodes
    global min_len
    global barcoded 
    global prog

    with open(dat_file_path, 'r') as dat_file:
        header = dat_file.readline()
        # Set up the lists of summed arrays
    
        target_sum_array_list = []
        targetc_sum_array_list = []
        for i in range(n_targets):
            target_len = target_lengths[i]
            # Figure out the maximum nunber of subsequences possible
            max_num_subseqs = target_len - min_len + 1
            # Calculate the count array length
            array_len = int(max_num_subseqs * (max_num_subseqs+1)/2)
            target_sum_array_list.append(np.zeros(array_len, dtype=int))
            targetc_sum_array_list.append(np.zeros(array_len, dtype=int))
        else:
            pass

        while True:
            # Read the sequence ID
            seq_ID = dat_file.readline()
            if not seq_ID:
                break
            # Get sequence info
            seq_ID = seq_ID.rstrip('\n')
            features = dat_file.readline().rstrip('\n')

            # Get count info
            for i in range(n_targets):
                # Get target count info
                target_sum_array_list[i] += read_count_line(dat_file)
                # Get complement count info
                targetc_sum_array_list[i] += read_count_line(dat_file)
        if prog:
            dat_file_name = dat_file_path.split("/")[-1]
            dat_file_name_list = dat_file_name.split("_")
            print(f"Finished counting {dat_file_name_list[1]} file {dat_file_name_list[3]} count data.")
        else:
            pass
        return(target_sum_array_list, targetc_sum_array_list)

def read_dat_file_all_parallel(dat_file_list, target_lengths):
    total_target_sum_array_list = []
    total_targetc_sum_array_list = []
    for i in range(n_targets):
        target_len = target_lengths[i]
        # Figure out the maximum nunber of subsequences possible
        max_num_subseqs = target_len - min_len + 1
        # Calculate the count array length
        array_len = int(max_num_subseqs * (max_num_subseqs+1)/2)
        total_target_sum_array_list.append(np.zeros(array_len, dtype=int))
        total_targetc_sum_array_list.append(np.zeros(array_len, dtype=int))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(read_dat_file_all, dat_file_list)
        for result in results:
            # Unpack the subsequence counts
            target_sum_array_list, targetc_sum_array_list = result
            for i in range(len(target_sum_array_list)):
                total_target_sum_array_list[i] += target_sum_array_list[i]
                total_targetc_sum_array_list[i] += targetc_sum_array_list[i]
    return (total_target_sum_array_list, total_targetc_sum_array_list)

def write_all_summed_counts_array(save_file_path, summed_counts_array):
    with open(save_file_path, 'w') as save_file:
        for i in range(len(summed_counts_array)):
            save_file.write(f"{summed_counts_array[i]}\t")

def create_time_test_dat_file(bin_range, save_file_name):
    with open(save_file_name, "w") as save_file:
        save_file.write("# Data info and analysis settings in count_settings_0-1ratio1-1_0.txt\n")
        for i in range(len(bin_range)):
            time = bin_range[i]
            start_time = f"{time.year}-{time.month:02d}-{time.day:02d}T{time.hour:02d}:{time.minute:02d}:{time.second:02d}Z"
            save_file.write(f"@seq_id{i} runid={i} read=16 ch=122 start_time={start_time} flow_cell_id=AGP710 protocol_group_id=ExactBarcode sample_id=0-1_ratio_1-1_0\n")
            save_file.write(f"{(i+1)/20} -1 -1\n")
            for j in range(4):
                for k in range(231):
                    save_file.write(f"{j+1} ")
                save_file.write("\n")


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="Analyzes dat files of subsequence \
        matches.")
    parser.add_argument("dat_folder", type=str, help="The path to the folder \
        with the dat files to analyze.")
    parser.add_argument("count_settings_path", type=str, help="The path to the file \
        containing count analysis settings.")
    parser.add_argument("save_file", type=str, help="The name extension to add to \
        save file names that results are saved to.")
    parser.add_argument("--sum", type=bool, help="If True, sums subsequence counts \
        together and saves the sums into ")
    parser.add_argument("--barcoded", type=bool, help="If True, the program will \
        sort results by the barcode.")
    parser.add_argument("--time", type=int, help="Bins counts by the \
        provided time step in minutes. Default is to not bin by time.")
    parser.add_argument("--run", type=int, help="The length of the sequencing run \
        in hours. Defaults to 24 hours.")
    parser.add_argument("--prog", type=bool, help="If True, prints progress \
        messages.")
    parser.add_argument("--Pbin", type=float, help="Bins counts by the ")
    parser.add_argument("--Phist", type=bool, help="If True, makes a histogram \
        of the average P(corr) of the sequences.")
    parser.add_argument("--Qbin", type=float, help="Bins counts by the ")
    parser.add_argument("--Qhist", type=bool, help="If True, makes a histogram \
        of the average Q scores of the sequences.")
    parser.add_argument("--Lhist", type=bool, help="If True, makes a histogram \
        of the average Q scores of the sequences.")
    parser.add_argument("--beep", type=int, help="Plays a sound using beepy when \
        the program finishes running. To pick a sound, provide an integer from \
        1-7. To not play a sound, set to 0. Defaults to 1.")
    args = parser.parse_args()

    # Parse the arguments
    dat_folder = args.dat_folder
    count_settings_path = args.count_settings_path
    save_file_name = args.save_file

    if args.sum:
        do_sum = args.sum
    else:
        do_sum = False

    if args.barcoded:
        barcoded = args.barcoded
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

    if args.Qbin:
        Qbin = args.Qbin
    else:
        Qbin = 0.0

    if args.prog:
        prog = args.prog
    else:
        prog = False

    if args.Qhist:
        Qhist = args.Qhist
    else:
        Qhist = False

    if args.Pbin:
        Pbin = args.Pbin
    else:
        Pbin = 0.0

    if args.Phist:
        Phist = args.Phist
    else:
        Phist = False

    if args.Lhist:
        Lhist = args.Lhist
    else:
        Lhist = False

    if args.beep:
        which_beep = args.beep 
    else:
        which_beep = 1


    # Get all the dat files
    dat_files = Path(dat_folder).rglob("*.dat")
    dat_file_list = [str(file) for file in dat_files]

    # Get the number of targets
    n_targets, target_list, target_lengths, n_barcodes, min_len, \
    handle_repeat_error, repeat_list, n_repeat = \
    get_count_settings(count_settings_path)

    
    # Get avg Q scores for histograms
    if Qhist or Phist:
        all_Q_score_arr = get_Q_scores_parallel(dat_file_list)
    else:
        pass

    # Make avg Q score histogram
    if Qhist:
        plt.hist(all_Q_score_arr, bins=np.arange(0,23,1))
        plt.title('Average Q Scores')
        plt.xlabel('Average Q Scores')
        plt.ylabel('Number of Sequences')
        plt.savefig(f'{dat_folder}/average_Q_score_hist_{save_file_name}.png')
    else:
        pass

    # Make avg P(corr) histogram
    if Phist:
        all_P_corr_arr = 1.0 - 10 ** (-all_Q_score_arr/10.0)
        plt.hist(all_P_corr_arr,bins=np.arange(0,1.1,0.05))
        plt.title('Average P(correct)')
        plt.xlabel('Average P(correct)')
        plt.ylabel('Number of Sequences')
        plt.savefig(f'{dat_folder}/average_P_corr_hist_{save_file_name}.png')

    else:
        pass

    # Make length histogram
    if Lhist:
        pass
    else:
        pass

    # Make save folder if counts are being summed
    if time_step > 0 or Pbin > 0 or Qbin > 0:
        save_folder = make_count_save_folder(dat_folder, save_file_name, time_step, Qbin, Pbin)
        start = time.perf_counter()
    else:
        pass

    if Pbin > 0.0 and time_step == 0.0:
        # Set some globals
        P_corr_bins = np.arange(0, 1, Pbin)
        
        # Sum and bin the counts according to average P(corr)
        total_target_sum_array_list, total_targetc_sum_array_list = \
        read_dat_file_bin_parallel(P_corr_bins, read_dat_file_Pbin, dat_file_list, target_lengths)

        # Write the summedd counts to save files
        write_summed_counts(save_folder, save_file_name, Pbin, Qbin, time_step, \
            total_target_sum_array_list, total_targetc_sum_array_list, prog)
    else:
        pass

    if Pbin == 0.0 and Qbin == 0.0 and time_step > 0.0:
        # Get the time range info
        start_time, end_time, time_step_range, time_range = get_time_range(time_step, run_length, dat_folder)
        time_step_td = dt.timedelta(minutes=time_step)

        # Sum and bin the counts according to average P(corr)
        total_target_sum_array_list, total_targetc_sum_array_list = \
        read_dat_file_bin_parallel(time_range, read_dat_file_time, dat_file_list, target_lengths)

        # Write the summedd counts to save files
        write_summed_counts(save_folder, save_file_name, Pbin, Qbin, time_step, \
            total_target_sum_array_list, total_targetc_sum_array_list, prog)
    else:
        pass

    # start = time.perf_counter()
    # total_target_sum_array_list, total_targetc_sum_array_list = read_dat_file_all_parallel(dat_file_list, target_lengths)
    # time_step = 1
    # for i in range(len(total_target_sum_array_list)):
    #     save_file_path = f"{dat_folder}/all_target_{i}_counts.txt"
    #     write_all_summed_counts_array(save_file_path, total_target_sum_array_list[i])
    #     save_file_path = f"{dat_folder}/all_target_comp_{i}_counts.txt"
    #     write_all_summed_counts_array(save_file_path, total_targetc_sum_array_list[i])

    if time_step > 0 or Pbin > 0 or Qbin > 0:
        end = time.perf_counter()
        print("Time elapsed: {} second(s)".format(end-start))
    else:
        pass
    
    if which_beep > 0:
        bp.beep(sound=which_beep)
    else:
        pass

    # # # Test functions
    # dat_file_path = "/Users/emilywu/OneDrive - Massachusetts Institute of Technology/Minion/BarcodeRatio/0-1ratio2-1/20210315_2220_MC-110826_0_AFO272_6d8d594b/counts/AFO272_fail_1d85e24b_0_0-1ratio2-1_counts.dat"
    # dat_file_path = "/Users/emilywu/OneDrive - Massachusetts Institute of Technology/Minion/ExactBarcode/0-1_ratio_1-1_0/20210619_2330_MC-110826_0_AGP710_fa5cb0a9/counts/AGP710_pass_bb5f66e6_293_0-1ratio1-1_0_counts.dat"
    # dat_file_path = "/Users/emilywu/OneDrive - Massachusetts Institute of Technology/Minion/ExactBarcode/0-1_ratio_1-1_0/20210619_2330_MC-110826_0_AGP710_fa5cb0a9/counts/AGP710_fail_bb5f66e6_0_0-1ratio1-1_0_counts.dat"


    # # Test get_Q_scores
    # avg_Q_score_arr = get_Q_scores(dat_file_path)
    # print(avg_Q_score_arr)
    # print(len(avg_Q_score_arr))


    # # Test binning decision
    # # Test P(corr) binning
    # P_corr_bins = np.arange(0, 1, 0.05)
    # seq_P_corr_arr = np.arange(0,1.025,0.025)
    # for seq_P_corr in seq_P_corr_arr:
    #     P_corr_bin = 0
    #     while seq_P_corr >= P_corr_bins[P_corr_bin] + Pbin and P_corr_bin + 1 < len(P_corr_bins):
    #         P_corr_bin += 1
    #     print(f"P_corr: {seq_P_corr:.2f}, Bin lower bound: {P_corr_bins[P_corr_bin]:.2f}")

    # # Test time binning
    # # Prepare lists of time steps
    # time_step = 30
    # start_time, end_time, time_step_range, time_range = get_time_range(time_step, run_length, dat_folder)
    # time_step_td = dt.timedelta(minutes=time_step)
    # print(len(time_range))
    # time_step_range = range(-30,run_length*60+30, 15)
    # test_time_range = []
    # for i in range(len(time_step_range)):
    #     delta_time = dt.timedelta(minutes=time_step_range[i])
    #     test_time_range.append(start_time+delta_time)
    
    # for read_time in test_time_range:
    #     current_bin = 0
    #     while read_time > time_range[current_bin] + time_step_td \
    #     and current_bin + 1 < len(time_range):
    #         current_bin += 1
    #     print(f"time: {read_time}, Bin lower bound: {time_range[current_bin]}, Bin # {current_bin}")


    # # Test read_dat_file_Pbin
    # target_sum_array_list, targetc_sum_array_list = read_dat_file_Pbin(dat_file_path)
    # print(target_sum_array_list[0][int(0.6/0.05)+1,:])
    # print(targetc_sum_array_list[0][int(0.6/0.05)+1,:])


    # Test make_count_save_folder
    # save_folder = make_count_save_folder(dat_folder, save_file_name, time_step, Qbin, Pbin)

    # Test read_dat_file_time
    # time_step = 30
    # start_time, end_time, time_step_range, time_range = get_time_range(time_step, run_length, dat_folder)
    # time_step_td = dt.timedelta(minutes=time_step)
    # print(start_time, end_time, time_step_range, time_range)
    # print(time_step_td)
    # target_sum_array_list, targetc_sum_array_list = read_dat_file_time(dat_file_path)
    # print(target_sum_array_list[0])


    # # Make test time binning file
    # time_step = 15
    # run_length = 24.5
    # time_step_td = dt.timedelta(minutes=time_step)
    # time_step_range = np.arange(0, run_length*60, time_step)
    # start_time = dt.datetime(2021, 6, 19, hour=23, minute=15)
    # time_range = []
    # save_file = "/Users/emilywu/OneDrive - Massachusetts Institute of Technology/Minion/ExactBarcode/0-1_ratio_1-1_0/20210619_2330_MC-110826_0_AGP710_fa5cb0a9/testcounts/test_time_counts.dat"
    # for i in range(len(time_step_range)):
    #     delta_time = dt.timedelta(minutes=time_step_range[i])
    #     time_range.append(start_time+delta_time)
    # create_time_test_dat_file(time_range, save_file)
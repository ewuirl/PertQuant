import sys
sys.path.append('../')
from simCRN.ml_nupack import gen_complement
from analyzeFASTQ.countmatches import read_settings
import numpy as np
import argparse
from pathlib import Path

def get_seq_info(seq_ID, features):
    seq_ID_list = seq_ID.split()
    run_time_str = seq_ID_list[4]
    start_time = dt.datetime(int(run_time_str[11:15]), int(run_time_str[16:18]), \
                int(run_time_str[19:21]), hour=int(run_time_str[22:24]), \
                minute=int(run_time_str[25:27]), second=int(run_time_str[28:30]))
    feature_list = features.split()
    avg_Q_score = feature_list[0]
    barcode_ID = feature_list[1]
    has_repeat_error = feature_list[2]
    return(start_time, avg_Q_score, barcode_ID, has_repeat_error)

def get_Q_score_hist_data(dat_file_path):
    global n_targets

    with open(dat_file_path, 'r') as dat_file:
        # Read in the lines
        lines = dat_file.readlines()
        index_arr = np.arange(0, len(lines)-2, 2+2*n_targets)

        # Make an array to store all the Q scores
        avg_Q_score_arr = np.array(len(index_arr))

        # Read the data and save it to the arrays
        for index in index_arr:
            seq_ID = lines[index].rstrip('\n')
            features = lines[index+1].rstrip('\n')
            start_time, avg_Q_score, barcode_ID, has_repeat_error = get_seq_info(seq_ID, features)
            avg_Q_score_arr[int(index_arr/(2+2*n_targets))] = avg_Q_score
        
        # Return the array
        return avg_Q_score_arr

# def read_dat_file(dat_file_path):
#     global n_targets
#     global barcoded 
#     global time_step
#     global prog

#     with open(dat_file_path, 'r') as dat_file:
#         header = dat_file.readline()
#         while True:
#             # Read the sequence ID
#             seq_ID = file.readline()
#             if not seq_ID:
#                 break
#             # Get sequence info
#             seq_ID = seq_ID.rstrip('\n')
#             features = dat_file.readline().rstrip('\n')
#             start_time, avg_Q_score, barcode_ID, has_repeat_error = get_seq_info(seq_ID, features)

#             # Get count info


#             if barcoded=False and time_step=0 and Qbi:




if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="Analyzes dat files of subsequence \
        matches.")
    parser.add_argument("dat_folder", type=str, help="The path to the folder \
        with the dat files to analyze.")
    parser.add_argument("settings_path", type=str, help="The path to the file \
        containing count analysis settings.")
    parser.add_argument("save_file", type=str, help="The name extension to add to \
        save file names that results are saved to.")
    parser.add_argument("--sum", type=bool, help="If True, sums subsequence counts \
        together and saves the sums into ")
    parser.add_argument("--barcoded", type=bool, help="If True, the program will \
        sort results by the barcode.")
    parser.add_argument("--time", type=int, help="Bins counts by the \
        provided time step in minutes. Default is to not bin by time.")
    parser.add_argument("--Qbin", type=float, help="Bins counts by the ")
    parser.add_argument("--prog", type=bool, help="If True, prints progress \
        messages.")
    parser.add_argument("--Qhist", type=bool, help="Produces and saves a histogram \
        of the average quality scores.")
    args = parser.parse_args()

    # Parse the arguments
    dat_folder = args.dat_folder
    settings_path = args.settings_path
    save_file_name = args.save_file

    if args.sum:
        do_sum = args.sum
    else:
        do_sum = False

    if args.barcoded:
        prog = args.barcoded
    else:
        prog = False

    if args.time:
        time_step = args.time
    else:
        time_step = 0

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

    # Get all the dat files

    if Qhist:
        pass
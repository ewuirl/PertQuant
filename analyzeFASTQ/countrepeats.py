from PertQuant.analyzeFASTQ.countmatches import check_repeat_errors
import argparse
from pathlib import Path
import os
import concurrent.futures
import time
import beepy as bp


def analyze_sequence_repeats(index):
    global lines 
    global repeat_list
    global n_repeat

    # Read in the sequence
    sequence = lines[index+1].rstrip("\n")
    has_repeat_error = check_repeat_errors(sequence, repeat_list, n_repeat)
    return has_repeat_error

def cut_Qscore(fastq_file_name, edited_file_name):
    with open(edited_file_name, 'w') as edited_file:
        with open(fastq_file_name, 'r') as fastq_file:
            index = 0
            while True:
                line = fastq_file.readline()
                if not line:
                    break
                index += 1
                if index % 4 == 1 or index % 4 == 2:
                    edited_file.write(line)
                else:
                    pass

def edit_fastq_Qscore(read_file_path):
    global save_folder
    global fastq_file_list
    # Make the edit file name
    read_file_list = read_file_path.split("/")
    read_file_name = read_file_list[-1].replace(".fastq","")
    edited_file_name = f"{save_folder}/{read_file_name}_no_Qscore.txt"
    cut_Qscore(read_file_path, edited_file_name)
    if prog:
        index += 1
        print(f"Analyzed file #{read_file_list[-1][-1]} (total # files: {len(fastq_file_list)}).")
    else:
        pass


if __name__ == "__main__":
    # Argument parser
    parser = argparse.ArgumentParser(description="Analyzes fastq files for \
        perfect sequence matches of different lengths")
    parser.add_argument("fastq_folder", type=str, help="The path to the folder \
        with the fastq files to analyze. Fastq files in subdirectories will also \
        be analyzed.")
    parser.add_argument("--cutQ", type=str, help="If True, creates copies of the \
        fastq files without the '+' and Q score lines. Defaults to False.")
    parser.add_argument("--serial", type=str, help="If True, counts matches serially \
        instead of in parallel. Defaults to False.")
    parser.add_argument("--note", type=str, help="Adds a note to the header file.")
    parser.add_argument("--prog", type=bool, help="If True, prints progress \
        messages. Defaults to False.")
    parser.add_argument("--beep", type=int, help="Plays a sound using beepy when \
        the program finishes running. To pick a sound, provide an integer from \
        1-7. To not play a sound, set to 0. Defaults to 1.")
    args = parser.parse_args()

    # Parse the arguments
    fastq_folder = args.fastq_folder
    if args.serial:
        serial = args.serial
    else:
        serial = False
    if args.cutQ:
        cutQ = args.cutQ
        # Make a save folder
        save_folder = f"{fastq_folder}/fastq_no_Qscore"
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)
    else:
        cutQ = False
    if args.note:
        note = args.note
    else:
        note = ""
    if args.prog:
        prog = args.prog
    else:
        prog = False
    if args.beep:
        which_beep = args.beep 
    else:
        which_beep = 1

    # Create list of fastq files to analyze
    fastq_files = Path(fastq_folder).rglob("*.fastq")
    fastq_file_list = [str(file) for file in fastq_files] 

    # Analyze all the fastq files
    file_index = 0
    start = time.perf_counter()

    # Serial counting
    if cutQ:
        if serial:
            for read_file_path in fastq_file_list:
                edit_fastq_Qscore(read_file_path)
        else:
            with concurrent.futures.ProcessPoolExecutor() as executor:
            # with concurrent.futures.ThreadPoolExecutor() as executor:
                results = executor.map(edit_fastq_Qscore, fastq_file_list)
    else:
        pass

    end = time.perf_counter()
    print("Time elapsed: {} second(s)".format(end-start))

    if which_beep > 0:
        bp.beep(sound=which_beep)
    else:
        pass
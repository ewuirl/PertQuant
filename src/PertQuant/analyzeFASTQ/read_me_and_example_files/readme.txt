Last Updated 09/01/2021

# Connecting to Minion with SSH

    To connect to the Minion, use the following command.

        ssh minit@IP.address

    The password is minit. The terminal should prompt you for the password.

    Notes:
        1.) You need to be connected to the same network as the Minion, which may
        require VPN.


# Copying Remote Files->Local Using SCP

    To copy files from the Minion to a local directory using scp, use the following 
    command. This is good for Unix-like operating systems. For Windows, you can use
    pscp.exe

        scp -r username@remoteHost:/remote/dir/file.txt /local/dir/

    To close the connection, use the following command.

        exit

    Notes:
        1.) The -r flag is recursive and will copy subfolders to the specified 
        directory as well.
        2.) You need to be connected to the same network as the Minion, which may
        require VPN.
        3.) If either directory has spaces in the path, enclose it in quotes.
        4.) To find the path of the data folder on the Minion, connect via SSH
        and go up directories until you find a subdirectory called data. The 
        data for a particular experiment will be in the subfolder
        /data/name_of_experiment/name_of_sample/

    Example:

        scp -r minit@IP.address:/remote/dir/file.txt "/my local/dir/"


# Updating the Minion Remotely
    
    If the Minion has issues updating locally, it can be updated remotely. To 
    remotely update the Minion, connect via SSH and run the following command.

        sudo apt update

    Reboot the Minion.

    Notes:
        1.) To connect via SSH, see "Connecting to Minion with SSH" above. 


# Unzipping .gz files

    As of 2021, Nanopore fastq files are saved as .fastq.gz files. Before analysis
    these files must be unzipped. To unzip all of the .gz files in a folder, use
    the following command.

        gunzip /path/to/.gz/files/*.gz


# Making Read Length Histograms with NanoPlot

    To make Read Length Histograms with NanoPlot, place a copy of 
    NanoPlot_hist.sh in /usr/local/bin/. Then move to the folder you want the 
    plots to be saved in and run the following command:

        NanoPlot_hist.sh --fastq fastq_file [fastq_file]

    Example:

        To look in multiple specific folders for fastq files:
        NanoPlot_hist.sh --fastq "/path/to/folder/with/fastq/files1/"*.fastq "/path/to/folder/with/fastq/files2/"*.fastq

        To look in all the subdirectories of a folder for fastq files:
            NanoPlot_hist.sh --fastq "/path/to/folder/with/fastq/files/"*/*.fastq



# Nanopore Folder Structure
    Experiment_Folder
        |__ Sample_Folder (a particular run)
            |__ Dated_Sample_Folder
                |__ fast5_fail
                |__ fast5_pass
                |__ fastq_fail
                |__ fastq_pass
                |__ other_reports


# Counting Subsequence Matches
    
    Counting subsequence matches is done by running countmatches.py. 

    Required files:
        (1) target sequence file
        (2) countmatches settings file
        (3) fastq file(s)

    Optional files:
        (4) barcode sequence file

    (1) Create a target sequence file (see ex_target_file.txt)
    (2) Create a countmatches.py settings file (see ex_countmatches_settings.txt) 
        with the desired settings. 
    (3) Nanopore places fastq files in the fastq_fail and fastq_pass folders
        (see Nanopore Folder Structure above).
    (4) See ex_barcode_file.txt

    Generated files:
        (5) .dat file(s)
        (6) count_settings file

    (5) subsequence count data files created by countmatches.py. Each fastq file 
        analyzed gets its own .dat file. See ex_count_data.dat for more details.
    (6) A file summarizing the sequences and settings used by countmatches.py.
        See ex_sum_count_data_settings.txt.
 
    Default Use
    For default use, place the countmatches settings file in the 
    Dated_Sample_Folder (see Nanopore Folder Structure above), and make sure the
    name ends in "_settings.txt". Then run the following command.

        python3 countmatches.py /path/to/Dated_Sample_Folder custom_name 

    Notes:
        (1) fastq files in /path/to/Dated_Sample_Folder and its subdirectories 
        will be analyzed.

    Other Use
      - The countmatches settings file can be placed elsewhere, but the path to the
        file must be provide with the --settings argument.

# Summing Count Data

    Summing count data for further analysis is done by running sum_count_dat.py.
    Files are placed in the .dat folder if the data is not binned. The data can 
    be binned by P (per base probability of calling a base correctly), Q score, 
    or time. If the data is binned, produced files are placed in a subfolder.

    Required files:
        (1) .dat file(s)
        (2) count_settings file

    (1) subsequence count data files created by countmatches.py. Each fastq file 
        analyzed gets its own .dat file. See ex_count_data.dat for more details.
    (2) A file summarizing the sequences and settings used by countmatches.py.
        This is created by countmatches.py See ex_sum_count_data_settings.txt.

    Generated files:
        (3) counts.txt files

    (3) counts text files contain summed subsequence count data. Separate files
        are created for each target and target complement. See 

    Default Use (No Binning)
        python3 sum_count_dat.py /path/to/dat_folder custom_name 

    Time-Binned Use
        python3 sum_count_dat.py /path/to/dat_folder custom_name --time TIME

    Q Score-Binned Use
        python3 sum_count_dat.py /path/to/dat_folder custom_name --Qbin TIME
        (not actually written)

    P-Binned Use
        python3 sum_count_dat.py /path/to/dat_folder custom_name --Pbin PBIN
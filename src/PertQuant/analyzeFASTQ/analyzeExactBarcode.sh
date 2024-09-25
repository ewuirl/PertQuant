#!/bin/bash

# place a copy in usr/local/bin

# define help function for usage
helpFunction()
{
    echo ""
    echo "Usage: analyzeExact_barcode.sh -r read_folder \
    -s settings_file -t time - p progress"
    echo "Run this in the Exact Barcode directory containing the read folders."
    echo -e "\t-r The name of the read folder. "
    echo -e "\t-s Optional. The name of the settings file to use."
    echo -e "\t-t Optional. The time bin to use for summing counts. If provided,\
    will also sum counts binned by time."
    echo -e "\t -p Optional. If True, prints the progress of counting and \
    summing the counts. Defaults to False."
    exit 1 # Exit script after printing help
}

# Get arguments
while getopts "r:s:t:p:" opt
do
    case "$opt" in
        r ) read_folder="$OPTARG" ;;
        s ) settings_file="$OPTARG" ;;
        t ) time="$OPTARG" ;;
        p ) progress="$OPTARG" ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
    esac
done

# Print helpFunction in case parameters are empty
if [ -z "$read_folder" ] 
then
    echo "Some or all of the required arguments were not provided.";
    helpFunction
fi

# Get the current working directory
cwd=$(pwd)

# Figure out if the read folder is in the working directory
if [ -d "$read_folder" ]
then
    :
else
    echo "$read_folder directory does not exist in $cwd."
    exit 1
fi

# Figure out if the settings file was provided
if [ -z "$settings_file" ];
then
    # Assume the settings file is in the count folder
    :
else
   # Check if the settings file is in the working directory
    if [ -d "$read_folder" ]
    then
        settings_file="$cwd/$settings_file"
        echo "updated settings file"
    else
        # Assume the full settings path was provided
        :
    fi
fi

# Figure out if progress should be printed
if [ -z "$progress" ]
then
    progress="False"
else
    :
fi

# Figure out where the subfolder is
cd $read_folder
subfolder=$(ls)

# Change directory to the python package directory
python_wd="/usr/local/lib/python3.7/site-packages/PertQuant/analyzeFASTQ"
cd $python_wd

# Count the matches
if [ -z "$settings_file" ]
then
    echo "Settings file not provided. Counting matches."
    python3 countmatches.py "$cwd/$read_folder/$subfolder" $read_folder --prog $progress
else
    echo "Settings file provided. Counting matches."
    python3 countmatches.py "$cwd/$read_folder/$subfolder" "$read_folder" --settings "$settings_file" --prog $progress
fi

# Sum the counts
# Define the counts folder
counts_folder="$cwd/$read_folder/$subfolder/counts"
# Go back to the python package directory
cd $python_wd

# Sum the counts
echo "Summing counts."
python3 sum_count_dat.py "$counts_folder" "$read_folder" --prog $progress
if [ -z "$time" ]
then
    :
else
    echo "Summing counts binned by provided time bin."
    # Sum the counts binned by time
    python3 sum_count_dat.py "$counts_folder" "$read_folder" --time 30 --prog $progress
fi
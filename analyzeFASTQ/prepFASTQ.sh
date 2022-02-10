#!/bin/zsh

# place a copy in usr/local/bin
# then make it executable

# define help function for usage
helpFunction()
{
    echo ""
    echo "Usage: analyze_Exact_barcode.sh -r read_folder "
    echo "Run this in the Exact Barcode directory containing the read folders."
    echo -e "\t-r The name of the read folder containing the sequencing experiment \
    data. "
    exit 1 # Exit script after printing help
}

# Get arguments
while getopts "r:" opt
do
    case "$opt" in
        r ) read_folder="$OPTARG" ;;
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
cd "$cwd/$read_folder"

# Unzip the fastq files
echo "Unzipping *.gz files"
gunzip */*.gz

# Run Nanoplot_hist.sh
subfolder=$(ls)
cd "$cwd/$read_folder/$subfolder"
echo "Running NanoPlot_hist"
NanoPlot_hist.sh --fastq */*.fastq -o "NanoPlot_Summary"
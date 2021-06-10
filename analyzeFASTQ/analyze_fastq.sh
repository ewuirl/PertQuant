#!/bin/bash

# define help function for usage
helpFunction()
{
   echo ""
   echo "Usage: $0 -s save_file \
   -f path/to/fastq/folder [path/to/other/fastq/folders ...] \
   -r fastq_file [other_fastq_files ...] \
   -b path/to/barcode/file \
   -t path/to/target/file \
   -o number_of_non-target_strands_in_target_file"
   echo -e "\t-s The name of the save file (no extension)."
   echo -e "\t-f The path to the folder(s) of fastq files to analyze, separated
   by commas (no spaces). If -r is not specified, all fastq files in folder(s) \
   are analyzed."
   echo -e "\t-r The name of fastq file(s) to analyze, separated by commas (no \
   spaces).  Optional, must be used with -f. If not specified, all fastq files \
   in specified folders are analyzed."
   echo -e "\t-b The path to the barcode file to use. Optional, will look for a \
   barcode file ending in 'barcodes.txt' in the folder of the save file first."
   echo -e "\t-t The path to the target file to use. Optional, will look for a \
   target file ending in 'targets.txt' in the folder of the save file first."
   echo -e "\t-o (int) The number of non-target strands in the target file. \
   Optional. The non-target strands should be placed at the end of the target \
   file."
   exit 1 # Exit script after printing help
}

# Get arguments
while getopts "s:f:r:b:t:o:" opt
do
   case "$opt" in
      s ) save_file="$OPTARG" ;;
      f ) folder_paths="$OPTARG" ;;
      r ) read_files="$OPTARG" ;;
      b ) barcode_path="$OPTARG" ;;
      t ) target_path="$OPTARG" ;;
      o ) other="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$save_file" ] || [ -z "$folder_paths" ] 
then
   echo "Some or all of the required parameters are empty";
   helpFunction
fi

# Save working directory of analyzeFASTQ
cwd=$(pwd)
cd ..
echo "Move to Nupack"
cd ..
echo "Move to $save_file"
cd $save_file
# Save path of save file
echo "Saving save file path"
save_file_path="$(pwd)/$save_file.txt"
echo "$save_file_path"

# Find barcode file
if [ -z "$barcode_path" ] 
then
   echo "No barcode path provided"
   # look for barcode file
   barcode_file=$(ls | grep barcodes.txt)
   echo "barcode_file = $barcode_file"
   # Check the number of barcode files found
   num_barcode_files=$(ls | grep -c barcodes.txt)
   echo "num barcode files: $num_barcode_files"
   # If 0 barcode files are found
   if [ $num_barcode_files == 0 ]
   then
      echo "Error: No barcode file found ending in \"barcodes.txt\""
      exit 1
   # If multiple barcode files are found
   elif [ $num_barcode_files -gt 1 ]
      then
      echo "Error: Multiple barcode files found ending in \"barcodes.txt\""
      echo "Found $num_barcode_files files:"
      echo "$barcode_file"
      echo "Please enter the name of desired \"barcodes.txt\""
      read barcode_file
      {
         echo "Saving barcode path"
         barcode_path="$(pwd)/$save_file.txt"
         echo "$barcode_path"
      } || {
         echo "Error: $barcode_file not found"
         cd "$cwd"
         exit 1
      }   
   # If 1 barcode file is found
   else
      echo "Saving barcode path"
      barcode_path="$(pwd)/$barcode_file"
      echo "barcode path: $(pwd)/$barcode_file"
   fi
else
   echo "Barcode path provided"
   pass
fi

# Find target file
if [ -z "$target_path" ] 
then
   echo "No target path provided"
   # look for target file
   target_file=$(ls | grep targets.txt)
   # Check the number of target files found
   num_target_files=$(ls | grep -c targets.txt)
   echo "num target files: $num_target_files"
   # If 0 target files are found
   if [ $num_target_files == 0 ]
   then
      echo "Error: No target file found ending in \"targets.txt\""
      exit 1
   # If multiple target files are found
   elif [ $num_target_files -gt 1 ]
      then
      echo "Error: Multiple target files found ending in \"targets.txt\""
      echo "Found $num_target_files files:"
      echo "$target_file"
      echo "Please enter the name of desired \"targets.txt\""
      read target_file
      {
         echo "Saving target path"
         target_path="$(pwd)/$target_file"
         echo "target path: $target_path"
      } || {
         echo "Error: $target_file not found"
         cd "$cwd"
         exit 1
      }   
   # If 1 target file is found
   else
      echo "Saving target path"
      target_path="$(pwd)/$target_file"
      echo "target path: $(pwd)/$target_file"
   fi
else
   :
   echo "target path provided"
   pass
fi

# Set internal field separator to "," instead of " "
IFS=","
# If the read file names are specified, create an array 
if [ -n read_files ]
then
   # Iterate through the read file names to get paths to the fastq files
   read_file_array=($read_files)
   # Make an array to keep track of whether or not a specified file name was 
   # found
   num_read_files=${#read_file_array[@]}
   echo "# read files: $num_read_files"
   found_array=()
   while [ $num_read_files -gt 0 ]
   do
      found_array+=(0)
      let num_read_files-=1
   done
   echo "found_array: ${found_array[@]}"
else
   :
fi

echo "folder paths: $folder_paths"
# Iterate through the folders to get paths to fastq files
folder_array=($folder_paths)
for folder in ${folder_array[@]}
do 
   echo "folder: $folder"
   # Move to folder directory
   cd "$folder"
   if [ -z read_files ]
   then 
      # Read all the fastq files in the folder
      folder_files=$(ls | grep ".fastq")
      # Save the paths of the fastq files to a text file
      for file in $(echo "$folder_files" | tr "\n" ",")
      do
         echo "$(pwd)/$file" >> "$cwd/fastq_files.txt"
      done
   else
      echo "read files: $read_files"
      i=0
      for file_name in ${read_file_array[@]}
      do 
         # Look for the file(s)
         files=$(ls | grep -w $file_name)
         is_file=$(ls | grep -c -w $file_name)
         # If there's one match, record it in the fastq file name text file.
         if [ $is_file == 1 ]
         then
            echo "$(pwd)/$files" >> "$cwd/fastq_files.txt"
            # Record finding the filename
            found_array[i]=1
         # If there's more than one match, record them in the fastq file name
         # text file
         elif [ $is_file -gt 1 ]
         then
            for file in $(echo "$files" | tr "\n" ",")
            do
               echo "$(pwd)/$file" >> "$cwd/fastq_files.txt"
            done
            # Record finding the filename
            found_array[i]=1
         else
            :
         fi
         let i+=1
      done
   fi
done

# Check if files were found for all the supplied read file names
if [ -n read_files ]
then
   # Add up how many of the names were successfully found
   total=0
   for count in ${found_array[@]}
   do
      let total+=$count
   done
   # If the total doesn't match the number of supplied names, print an error
   # statement and ask if the user wants to continue
   if [ $total != ${#read_file_array[@]} ]
   then
      echo "Error: One or more of the supplied fastq file names were not found."
      echo "Continue analyzing? [Y/N]"
      read continue
      # If yes, continue
      if [ $continue == "y" ] || [ $continue == "Y" ]
      then
         :
      # Otherwise exit the script.
      else 
         exit 1
      fi
   else
      :
   fi
else
   :
fi


# IFS=" "
# echo "Moving to analyzeFASTQ"
# cd "$cwd"


# figure out start time
start_time=`date`

# echo "Move $save_file files to save folder"
# mv $save_file.* "$swd"
# echo "Move $barcode_file to barcode folder"
# mv $barcode_file "$bwd"
# echo "Move $target_file to target folder"
# mv $target_file "$twd"

echo "I'm done"
echo "start time: $start_time"
end_time=`date`
echo "end time: $end_time"
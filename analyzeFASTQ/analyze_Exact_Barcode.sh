#!/bin/bash

# place a copy in usr/local/bin

# define help function for usage
helpFunction()
{
   echo ""
   echo "Usage: $0 -r read_folder \
   -t settings_file"
   echo -e "\t-r The name of the read folder."
   echo -e "\t-s Optional. The name of the settings file to use."
   exit 1 # Exit script after printing help
}

# Get arguments
while getopts "s:t:" opt
do
   case "$opt" in
      s ) save_file="$OPTARG" ;;
      t ) settings_file="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$read_folder" ] 
then
   echo "The name of the read folder must be specified";
   helpFunction
fi

# Get the current working directory
cwd=$(pwd)

# Figure out if the read folder is in the working directory
if [ -d "$read_folder" ];
then
    pass
    
else
   echo "$DIR directory does not exist."
fi

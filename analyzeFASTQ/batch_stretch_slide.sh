#!/bin/bash

# Save working directory of analyzeFASTQ
cwd=$(pwd)

exact_dir="/Users/emilywu/OneDrive - Massachusetts Institute of Technology/Minion/ExactBarcode"
# settings_file="/Users/emilywu/OneDrive - Massachusetts Institute of Technology/Minion/ExactBarcode/A-Bratio_settings.txt"

echo $exact_dir

# folder_array=( "0-1ratio1-1_z" "A-Bratio1-1_control" "A-Bratio1-1_z" \
#     "A-Bratio1-1_zCC" "A-Bratio1-2d" "A-Bratio1-2s" "A-Bratio1-4d" "A-Bratio1-4s" \
#     "A-Bratio1-4sv1" "A-Bratio1-4sv1_1" "A-Bratio2-1sv1" "A-Bratio2-1sv1_1" \
#     "A-Bratio4-1sv1" "A-Bratio4-1sv1_1" "A-GC30" "A-GC50" "A-GC70" )

folder_array=( "A-GC30" "A-GC50" "A-GC70" )

for folder in ${folder_array[@]}
do
    # Customize settings_file
    settings_file="/Users/emilywu/OneDrive - Massachusetts Institute of Technology/Minion/ExactBarcode/$folder""_settings.txt"
    # Go to the folder
    folder_path=$exact_dir"/"$folder
    cd "$folder_path"
    # Get the dat folder
    dat_folder=$folder_path/$(ls)
    echo $dat_folder
    echo $settings_file
    # # Go back to the current working directory
    # cd "$cwd"
    # # Slide count
    # echo "Slide counts"
    # python3 countmatches.py "$dat_folder" "slide_counts_stretch" --settings \
    # "$settings_file" --method "slide" --max_stretch 22 --save_path "slide_counts_stretch"
    # # Sum the counts
    # echo "Sum counts"
    # python3 sum_count_slide_dat.py "$dat_folder/slide_counts_stretch" "slide_counts_stretch" --pf True --max_stretch 22
    # # Time bin the counts
    # echo "Time bin sum counts"
    # python3 sum_count_slide_dat.py "$dat_folder/slide_counts_stretch" "slide_counts_stretch" --max_stretch 22 --time 30


done
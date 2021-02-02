#!/bin/bash

# define help function for usage
helpFunction()
{
   echo ""
   echo "Usage: $0 -f file_name -r n_runs -i Cmin -a Cmax -n N -m M -l L"
   echo -e "\t-f The name of the input file (no extension)."
   echo -e "\t-r The number of runs to calculate."
   echo -e "\t-i The minimum value of Ci"
   echo -e "\t-a The maximum value of Ci"
   echo -e "\t-n The number of A strands"
   echo -e "\t-m The number of B strands"
   echo -e "\t-l The number of C strands"
   echo -e "\t-s The skew fraction. Optional, must be used with -c"
   echo -e "\t-c Cskew, the minimum value to skew Ci towards. Optional, 
            must be used with -s"
   exit 1 # Exit script after printing help
}

# Get arguments
while getopts "f:i:a:r:n:m:l:s:c:" opt
do
   case "$opt" in
      f ) file_name="$OPTARG" ;;
      r ) n_runs="$OPTARG" ;;
      i ) Cmin="$OPTARG" ;;
      a ) Cmax="$OPTARG" ;;
      n ) N="$OPTARG" ;;
      m ) M="$OPTARG" ;;
      l ) L="$OPTARG" ;;
      s ) fskew="$OPTARG" ;;
      c ) Cskew="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$file_name" ] || [ -z "$n_runs" ] || [ -z "$N" ] || [ -z "$M" ] \
   || [ -z "$L" ] || [ -z "$Cmin" ] || [ -z "$Cmax" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi
# figure out start time
start_time=`date`

echo "file name: $file_name"
echo "n_runs: $n_runs" 
echo "N: $N" 
echo "M: $M" 
echo "L: $L" 
echo "Cmin = $Cmin"
echo "Cmax = $Cmax"
echo "fskew = $fskew"
echo "Cskew = $Cskew"

cd ..
echo "Move to Nupack"
cd ..
echo "Move to $file_name"
cd $file_name
echo "Move .in file to simCRN"
mv $file_name.in ../
cd ../
mv $file_name.in pertquant
cd pertquant
mv $file_name.in simCRN
cd simCRN

echo "Running complexes -material dna -pairs -mfe -degenerate"
complexes -material dna -pairs -mfe -degenerate $file_name

echo "Making data file"
Ai=1e-6

if [ -z "$fskew" ] && [ -z "$Cskew" ]
then
   python3 init_dat.py $file_name $n_runs $Cmin $Cmax $N $M $L $Ai
elif [ -n "$fskew" ] && [ -z "$Cskew" ]
then 
   echo "Cskew is empty";
   helpFunction
elif [ -z "$fskew" ] && [ -n "$Cskew" ]
then
   echo "fskew is empty";
   helpFunction
else
   python3 init_dat.py $file_name $n_runs $Cmin $Cmax $N $M $L $Ai -s $fskew -c $Cskew
fi

# Make while loop to complete runs
while [ $n_runs -gt 0 ]
do
	echo "current run: $n_runs"
   if [ -z "$fskew" ] && [ -z "$Cskew" ]
   then
      python3 gen_con.py $file_name $Cmin $Cmax $N $M $L
   else
      python3 gen_con.py $file_name $Cmin $Cmax $N $M $L -s $fskew -c $Cskew
   fi
   concentrations -pairs -sort 0 -quiet $file_name
   python3 write_dat.py $file_name $N $M $L
	n_runs=`expr $n_runs - 1`
done

echo "Move $file_name files to pertquant"
mv $file_name.* ../
cd ../
echo "Move $file_name files to Nupack"
mv $file_name.* ../
cd ../
echo "Move $file_name files to $file_name"
mv $file_name.* $file_name

echo "I'm done"
echo "start time: $start_time"
end_time=`date`
echo "end time: $end_time"
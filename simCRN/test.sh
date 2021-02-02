#!/bin/bash

echo "move up"
cd ..
echo "move to hoho"
cd hoho
ls

mv 4-4-4-asym-AB-AC.in ../
cd ../
mv 4-4-4-asym-AB-AC.in simCRN
cd simCRN

printf "Running complexes -material dna -pairs -mfe -degenerate \n"
complexes -material dna -pairs -mfe -degenerate 4-4-4-asym-AB-AC

echo "make init file"
python3 init_dat2.py 4-4-4-asym-AB-AC 5 0 1e-6 4 4 4 1e-6

n_runs=2
echo "Nruns = $n_runs"

while [ $n_runs -gt 0 ]
do
	echo "current run: $n_runs"
	echo "gen a con file"
	python3 gen_con2.py 4-4-4-asym-AB-AC 0 1e-6 4 4 4
	concentrations -pairs -sort 0 -quiet 4-4-4-asym-AB-AC
	python3 write_dat2.py 4-4-4-asym-AB-AC 4 4 4
	n_runs=`expr $n_runs - 1`
done

echo "move filed to hoho"
mv 4-4-4-asym-AB-AC.* ..
cd ..
mv 4-4-4-asym-AB-AC.* hoho
import random as rand 
import numpy as np
import math
# import argparse


# parser = argparse.ArgumentParser(description='Generates a .txt file.')
# parser.add_argument('file_name', type=str, help='The name of the .txt file to be made.')
# parser.add_argument('--Cmin', help='Cmin is minimum Ci value (real)')
# parser.add_argument('Cmax', help='Cmax is maximum Ci value (real)')
# parser.add_argument('d', help='d is minimum difference between Ci values (real)')
# parser.add_argument('L', type=int, help='L is the number of C strands (int)')
# parser.add_argument('Nruns', type=int, help='L is the number of C strands (int)')
# parser.add_argument('-s', '--skew', help='the fraction of the data set to skew')
# parser.add_argument('-c', '--Cskew', help='the minimum value of Ci to skew the data towards')
# args = parser.parse_args()

# # Unpack the args
# file_name = args.file_name
# if args.Cmin:
# 	Cmin = float(args.Cmin)
# else: 
# 	Cmin = 0.0
# Cmax = float(args.Cmax)
# d = float(args.d)
# L = args.L
# Nruns = args.Nruns
# if args.skew:
# 	fskew = float(args.skew)
# 	Cskew = float(args.Cskew)
# else:
# 	fskew = 0.0
# 	is_skew = 1.0

Cmin = 0.0
Cmax = 1.0
d = 0.1
fskew = 0.5
Cskew = 0.8
is_skew = 1.0
file_name = 'test_gen_con_exp'
Nruns = 10
L = 4

# Seed the random number generator
rand.seed()

# Create an array of allowed concentrations
conc_array = np.arange(Cmin, Cmax+d,d)

# If the data set is skewed, create an array of skewed concentrations
skew_conc_array = np.arange(Cskew, Cmax+d, d)
skew_conc_array = skew_conc_array[skew_conc_array <= Cmax]

# Check the concentration lists
print(skew_conc_array)
print(conc_array)

with open(file_name + '.txt', 'w') as file:
	# Write the header
	for i in range(L):
		file.write('C{}'.format(i))
		file.write('	')
	file.write('\n')

	for n in range(Nruns):
		# If the dataset is skewed, determine if this run is skewed.
		if fskew != 0:
			is_skew = rand.uniform(0,1)
		# If the dataset is not skewed, pass.
		else:
			pass

		# Check if is_skew is less than f_skew to see if the run is skewed.
		if is_skew <= fskew:
			# Figure out which strand to skew the concentration for
			which_skew = rand.choice(range(L))
			for i in range(L):
				# Generate the skewed value
				if i == which_skew:
					Ci = rand.choice(skew_conc_array)
				# Generate the non-skewed values
				else:
					Ci = rand.choice(conc_array)
				file.write('{:.2f}'.format(Ci))
				file.write('	')
		# Generate concentrations for a non-skewed run.
		else:
			print('not skewed')
			for i in range(L):
				Ci = rand.choice(conc_array)
				file.write('{:.2f}'.format(Ci))
				file.write('	')

		file.write('\n')

print('Done writing file')








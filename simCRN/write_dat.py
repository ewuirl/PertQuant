# import sys
# sys.path.append('../')
from ml_nupack import calc_Am_NP
import argparse

parser = argparse.ArgumentParser(description='Initializes a simCRN data file')
parser.add_argument('file_name', type=str, help='The name of the .dat file to \
	edit.')
parser.add_argument('N', type=int, help='N is the number of A strands (int)')
parser.add_argument('M', type=int, help='M is the number of B strands (int)')
parser.add_argument('L', type=int, help='L is the number of C strands (int)')
args = parser.parse_args()

# Unpack the arguments
file_name = args.file_name
N = int(args.N)
M = int(args.M)
L = int(args.L)

# Open the data file
with open(file_name + '.txt', 'a') as dat:

	# Figure out the Ci
	with open(file_name + '.con','r') as con:
		con_lines = con.readlines()
		for i in range(len(con_lines) - N - M):
			Ci_line = con_lines[N + M + i]
			Ci = Ci_line.rstrip('\n')
			# Write the Ci to the file
			dat.write(str(Ci) + '\t')

	# Calculate Am_array
	Am_array = calc_Am_NP(file_name, N, M, L)
	for i in range(len(Am_array)):
		Am = Am_array[i]
		# Write the Am to the file
		dat.write(str(Am) + '\t')

	# Add a newline character so the next sample gets written to the next line.
	dat.write('\n')


	



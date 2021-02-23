import argparse

parser = argparse.ArgumentParser(description='Initializes a simCRN data file')
parser.add_argument('file_name', type=str, help='The name of the .con file to \
	be made.')
parser.add_argument('N_runs', type=int, help='N_runs is the number of samples \
	for the data set')
parser.add_argument('Cmin', help='Cmin is minimum Ci value (real)')
parser.add_argument('Cmax', help='Cmax is maximum Ci value (real)')
parser.add_argument('N', type=int, help='N is the number of A strands (int)')
parser.add_argument('M', type=int, help='M is the number of B strands (int)')
parser.add_argument('L', type=int, help='L is the number of C strands (int)')
parser.add_argument('Ai', help='Ai is the initial A concentration \
	(real)')
parser.add_argument('-s', '--skew', help='the fraction of the data set to skew')
parser.add_argument('-c', '--Cskew', help='the minimum value of Ci to skew the data towards')
parser.add_argument('-n', '--note', help='notes about the data')
args = parser.parse_args()

# Unpack the args
file_name = args.file_name
N_runs = int(args.N_runs)
Cmin = float(args.Cmin)
Cmax = float(args.Cmax)
N = args.N
M = args.M
L = args.L
Ai = float(args.Ai)

if args.skew:
	fskew = float(args.skew)
	Cskew = float(args.Cskew)
else:
	fskew = "None"
	Cskew = "None"

if args.note:
	note = args.note
else:
	note = "None"

with open(file_name + '.txt', 'w') as file:
	file.write('# N = {:d} \n'.format(N))
	file.write('# M = {:d} \n'.format(M))
	file.write('# L = {:d} \n'.format(L))
	file.write('# N_runs = {:d} \n'.format(N_runs))
	file.write('# Cmin = {} \n'.format(Cmin))
	file.write('# Cmax = {} \n'.format(Cmax))
	file.write('# Ai = {} \n'.format(Ai))
	file.write('# fskew = {} \n'.format(fskew))
	file.write('# Cskew = {} \n'.format(Cskew))
	file.write('# note = {} \n'.format(note))
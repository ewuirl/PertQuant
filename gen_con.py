# import sys
# sys.path.append('../')
from ml_nupack import gen_con_file as gen_con_file
from ml_nupack import gen_skewed_con_file as gen_skewed_con_file
import argparse

parser = argparse.ArgumentParser(description='Generates a .con file.')
parser.add_argument('file_name', type=str, help='The name of the .con file to be made.')
parser.add_argument('Cmin', help='Cmin is minimum Ci value (real)')
parser.add_argument('Cmax', help='Cmax is maximum Ci value (real)')
parser.add_argument('N', type=int, help='N is the number of A strands (int)')
parser.add_argument('M', type=int, help='M is the number of B strands (int)')
parser.add_argument('L', type=int, help='L is the number of C strands (int)')
parser.add_argument('-s', '--skew', help='the fraction of the data set to skew')
parser.add_argument('-c', '--Cskew', help='the minimum value of Ci to skew the data towards')
args = parser.parse_args()

# Unpack the args
file_name = args.file_name
Cmin = float(args.Cmin)
Cmax = float(args.Cmax)
N = args.N
M = args.M
L = args.L


# Generate the .con file
if args.skew:
	fskew = float(args.skew)
	Cskew = float(args.Cskew)
	gen_skewed_con_file(file_name, Cmin, Cmax, Cskew, fskew, N, M, L)
else:
	gen_con_file(file_name, Cmin, Cmax, N, M, L)
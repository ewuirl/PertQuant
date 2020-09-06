import argparse

parser = argparse.ArgumentParser(description='Generates a .fasta file of DNA \
	strands with specified GC content and Tm range. Tm is computed using the \
	nearest-neighbor model. Sequences with 5 or more repeated bases are \
	excluded. Sequences with 4 or more double nucleotide sequence repeats are \
	also excluded.')
parser.add_argument('file_name', type=str, help='The name of the .fasta file \
	to be made.')
parser.add_argument('N', type=int, help='N is the number of strands to make.')
parser.add_argument('L', type=int, help='L is the strand length.')
parser.add_argument('--GCmin', help='GCmin is the minimum GC content (0-1.0). \
	Defaults to 0.40.')
parser.add_argument('--GCmax', help='GCmax is the maximum GC content (0-1.0). \
	Defaults to 0.60.')
parser.add_argument('-c', '--conc', help='conc is the per strand concentration \
	(uM). Defaults to 1 uM.')
parser.add_argument('-s', '--salt', help='s is the salt concentration (mM). \
	Defaults to 150 mM NaCl.')
parser.add_argument('--Tmin', help='Tmin is the minimum Tm (C). Defaults to \
	58C.')
parser.add_argument('--Tmax', help='Tmax is the maximum Tm (C). Defaults to \
	68C.')
args = parser.parse_args()

# Define constant values
R = 1.987 		# cal / K mol
K = 273.15 		# To convert from Celsius to Kelvin

# Unpack the args
file_name = args.file_name
N = args.N
L = args.L

if args.GCmin:
	GCmin = float(args.GCmin)
else:
	GCmin = 0.4
if args.GCmax:
	GCmax = float(args.GCmax)
else:
	GCmax = 0.6 
if args.conc:
	conc = float(args.conc) * 1e-6	# Convert to molar
else:
	conc = 1e-6
if args.salt:
	salt = float(args.salt) * 1e-3	# Convert to molar
else:
	salt = 150e-3	# Convert to molar
if args.Tmin:
	Tmin = args.Tmin + K	# Convert to Kelvin
else:
	Tmin = 58 + K 	# Convert to Kelvin
if args.Tmax:
	Tmax = args.Tmax + K 	# Convert to Kelvin
else:
	Tmax = 68 + K 	# Convert to Kelvin


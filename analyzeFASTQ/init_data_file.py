import argparse

parser = argparse.ArgumentParser(description='Initializes a data file to save \
    results from fastq file analysis')
parser.add_argument('file_name', type=str, help='The name of the data file to \
    be made.')
parser.add_argument('-N_runs', type=int, help='N_runs is the number of samples \
    for the data set')
parser.add_argument('-Pmin', help='Pmin is minimum i value (real)')
parser.add_argument('-Pmax', help='Pmax is maximum i value (real)')
parser.add_argument('-N', type=int, help='N is the number of A strands (int)')
parser.add_argument('-M', type=int, help='M is the number of B strands (int)')
parser.add_argument('-L', type=int, help='L is the number of perturbations (int)')
parser.add_argument('-Ai', help='Ai is the initial A concentration \
    (real)')
parser.add_argument('-s', '--skew', help='the fraction of the data set to skew')
parser.add_argument('-c', '--Cskew', help='the minimum value of Ci to skew the data towards')
parser.add_argument('-n', '--note', help='notes about the data')
args = parser.parse_args()

# Unpack the args
file_name = args.file_name

if args.N_runs:
    N_runs = int(args.N_runs)
else:
    N_runs = ""

if args.Pmin:
    Pmin = float(args.Pmin)
else:
    Pmin = ""

if args.Pmax:
    Pmax = float(args.Pmax)
else:
    Pmax = ""

if args.N:
    N = args.N
else:
    N = ""

if args.M:
    M = args.M
else:
    M = ""

if args.L:
    L = args.L
else:
    L = ""

if args.Ai:
    Ai = float(args.Ai)
else:
    Ai = ""

if args.note:
    note = args.note
else:
    note = ""

with open(file_name + '.txt', 'w') as file:
    file.write(f'# N = {N} \n'.format(N))
    file.write(f'# M = {M} \n'.format(M))
    file.write(f'# L = {L} \n'.format(L))
    file.write(f'# N_runs = {N_runs} \n'.format(N_runs))
    file.write(f'# Pmin = {Pmin} \n'.format(Pmin))
    file.write(f'# Pmax = {Pmax} \n'.format(Pmax))
    file.write(f'# Ai = {Ai} \n'.format(Ai))
    file.write(f'# note = {note} \n'.format(note))
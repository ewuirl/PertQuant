import matplotlib.pyplot as plt 
import random as rand
import numpy as np

# Initialize the random number generator
rand.seed()

def generate_strand(length):
	'''This function takes desired strand length (integer) as an input, randomly 
	generates a DNA sequence of the desired length, and outputs the sequence
	as a string.'''

	# Make sure length is an integer
	assert length == int(length), 'length must be an integer.'
	# If the length is an integer but provided as a float, convert it to an int.
	length = int(length)

	# Initialize the strand sequence
	strand = ''

	# Randomly choose bases and add them to the sequence
	for i in range(length):
		base = rand.choice(['A','T','C','G'])
		strand = strand + base
	return strand

def generate_strands(N_strands, length=0, min_l=0, max_l=0):
	'''This function takes an integer number of strands N_strands. It outputs 
	a list of unique DNA sequences (str) that is N_strands long. It also takes
	in optional integer arguments of length, min_l, and max_l. If length > 0, 
	all the DNA sequences will be of length=length. If max_l > min_l > 0, the 
	DNA sequences will vary randomly in length between min_l and max_l.'''

	# Create a list to store the sequences in.
	strand_list = []

	# Error Handling
	# Make sure N_strands is an integer
	assert N_strands == int(N_strands), 'N_strands must be an integer.'
	# If N_strands is an integer but provided as a float, convert it to an int.
	N_strands = int(N_strands)

	# Make sure either length or min_l and max_l are being used for strand #
	# length 
	assert (length > 0 and min_l == 0 and max_l == 0) or \
	(length == 0 and min_l > 0 and max_l > 0), \
	'length or both min_l and max_l must be nonzero and positive.'
	
	# More error handling
	# Varying lengths
	if length == 0 and min_l > 0 and max_l > 0:
		# More error handling
		# Make sure max_l > min_l > 0.
		assert 0 < min_l < max_l, \
		'max_l must be greater than min_l, and both must be greater than zero.'

		# Make sure min_l is an integer
		assert min_l == int(min_l), 'min_l must be an integer.'
		# If the length is an integer but provided as a float, convert it to an int.
		min_l = int(min_l)

		# Make sure max_l is an integer
		assert max_l == int(max_l), 'max_l must be an integer.'
		# If the length is an integer but provided as a float, convert it to an int.
		max_l = int(max_l)

	else:
		# Make sure length is an integer
		assert length == int(length), 'length must be an integer.'
		# If the length is an integer but provided as a float, convert it to an int.
		length = int(length)
		

	# Generate the sequences and add them to strand_list
	for i in range(N_strands):
		
		# Varying strand length
		if min_l > 0 and max_l > 0:
			# Randomly decide the length
			length = rand.randint(min_l, max_l)
		else:
			pass
		
		# Generate the strand
		strand = generate_strand(length)
		
		# Add the sequence if it is not in strand_list already.
		if strand not in strand_list:
			strand_list.append(strand)

	return strand_list


def pick_base(base):
	'''This function takes a string representing a DNA base (A C G T) and 
	returns a base that is not its complement.'''
	base_list = ['A','T','C','G']
	comp_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
	# Remove the complementary base from the list of bases
	comp_base = comp_dict.get(base)
	base_list.remove(comp_base)
	# Get a new base
	new_base = rand.choice(base_list)
	return(new_base)


def gen_complement(strand, comp=1):
	'''This function takes a string representing a DNA strand, and generates a 
	semi-complementary strand. The degree of complementarity is specified by 
	comp, which should be between [0,1]. The default value of comp = 1.'''
	assert comp >= 0 and comp <= 1, 'comp must be between [0,1]'
	# Figure out how long the strand is
	N = len(strand)

	# Figure out how much of the strand should be complementary
	N_comp = round(comp*N)

	# A dictionary specifying base pairs
	comp_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

	# Start the new strand
	new_strand = ''

	# Start making the complementary strand
	for i in range(N):
		# Pick out the base
		base = strand[i]
		# Find its complement
		if i < N_comp:
			new_base = comp_dict.get(base)
		# Or pick a base that isn't its complement
		else:
			new_base = pick_base(base)
		new_strand = new_base + new_strand

	return(new_strand)


def generate_input_file(file_name, strand_list, Lmax):
	'''This function takes a file name, a list of strands, and Lmax, the maximum
	complex size. The output is a file_name.in file with the following entries
	on  separate lines:
	The number of distinct strand species
	The sequences for each distinct strand species (each on a separate line)
	The maximum complex size.'''

	assert Lmax == int(Lmax), 'Lmax must be an integer'

	# Open a file of the right name
	with open(file_name + '.in', 'w') as input_f:
		# Write the number of distinct strands in
		N_strands = len(strand_list)
		input_f.write('{}\n'.format(N_strands))
		# Write the strands in
		for i in range(len(strand_list)):
			input_f.write(strand_list[i] + '\n')
		# Write Lmax
		input_f.write('{}'.format(Lmax))
	print('{}.in file generated.'.format(file_name))
	

def gen_con_file(file_name, Cmin, Cmax, N, M, L):
	'''This function takes:
	file_name = a string specifying the name of the .eq file to read.
	Cmin = the minimum value of Ci
	Cmax = the maximum value of Ci
	N = the number of A strands in the system
	M = the number of B strands in the system
	L = the number of C strands in the system
	It generates a .con file to be used with concentrations. The concentrations
	of A and B are 1 uM, and the concentrations of C are randomly generated 
	between [Cmin, Cmax].
	'''
	# Open a file of the right name
	with open(file_name + '.con', 'w') as con:
		for i in range(N+M):
			con.write('1e-6\n')
		for l in range(L):
			Ci = rand.uniform(Cmin, Cmax)
			con.write('{}\n'.format(Ci))


def calc_Am_NP(file_name, N, M, L):
	'''This function takes:
	file_name = a string specifying the name of the .eq file to read.
	N = the number of A strands in the system
	M = the number of B strands in the system
	L = the number of C strands in the system
	It calculates {Am} from the .eq file and returns it as a Numpy
	array.'''
	with open(file_name + '.eq', 'r') as f:
		# Read in the lines
		lines = f.readlines()

		# Figure out how big the header is
		is_header = True
		header = 0
		while is_header == True:
			if '%' in lines[header]:
				header = header + 1
			else:
				is_header = False

		# Check the number of strands
		n_strand_line = lines[header - 3 - N - M - L]
		n_strand_list = n_strand_line.split(' ')
		n_strands = int(n_strand_list[-1])
		assert n_strands == int(N + M + L), 'n_strands must be equal to N+M+L'

		# Figure out the init Ai concentrations
		Ai_lines = lines[8:8+N]
		Ai_array = np.empty(N)
		for i in range(len(Ai_lines)):
			Ai_line = Ai_lines[i]
			Ai_line_list = Ai_line.split(' ')
			Ai = Ai_line_list[4]
			Ai = float(Ai)
			Ai_array[i] = Ai

		# Initialize the AB_array
		AB_array = np.zeros(N)
		# Figure out the AB concentrations
		start = header + N 
		for i in range(N):
			start = start + M + L + N - i
			for j in range(M):
				AB_line = lines[start + j]
				AB_line = AB_line.rstrip('	\n')
				AB_line_list = AB_line.split('	')
				AiBj = float(AB_line_list[-1])
				AB_array[i] = AB_array[i] + AiBj

		# Calculate Am_array
		Am_array = Ai_array - AB_array	
		return Am_array


def check_conc(file_name):
	'''Check all complexes in a file_name.eq file with an equilibrium 
	concentration > 1e-8. This function writes these complex .eq lines to a file
	 called check.txt'''
	with open(file_name + '.eq', 'r') as f:
		# Read in the lines
		lines = f.readlines()
		# Figure out how big the header is
		header = 0
		for i in range(int(len(lines)/2)):
			split = lines[i].split('%')
			if len(split) > 1:
				header = header + 1
		with open('check.txt', 'w') as check:
			for i in range(len(lines)-header):
				line = lines[header + i].rstrip('	\n')
				line_list = line.split('	')
				conc = float(line_list[-1])
				if conc >= float(1e-08):
					check.write(line + '\n')

def get_Kvals(file_name, N, M, L):
	'''This function takes
	file_name = a string specifying the name of the .eq file to read.
	N = the number of A strands in the system
	M = the number of B strands in the system
	L = the number of C strands in the system
	It reads file_name.eq and generates K_AB, KBC, K_AC, K_AA, K_BB, and K_CC
	numpy arrays. These arrays are returned as a tuple:
	(K_AB, KBC, K_AC, K_AA, K_BB, K_CC)
	'''
	with open(file_name + '.ocx', 'r') as f:
		# Read in the lines
		lines = f.readlines()

		# Figure out how big the header is
		header = 0
		for i in range(int(len(lines)/2)):
			split = lines[i].split('%')
			if len(split) > 1:
				header = header + 1

		# Check the number of strands
		n_strand_line = lines[header - 3 - N - M - L]
		n_strand_list = n_strand_line.split(' ')
		n_strands = int(n_strand_list[-1])
		assert n_strands == int(N + M + L), 'n_strands must be equal to N+M+L'

		# Initialize arrays
		KAB_array = np.zeros((N,M))
		KAC_array = np.zeros((N,L))
		KBC_array = np.zeros((M,L))
		KAA_array = np.zeros((N,N))
		KBB_array = np.zeros((M,M))
		KCC_array = np.zeros((L,L))

def gen_skewed_con_file(file_name, Cmin, Cmax, Cskew, f_skew, N, M, L):
	'''This function takes:
	file_name = a string specifying the name of the .eq file to read.
	Cmin = the minimum value of Ci
	Cmax = the maximum value of Ci
	Cskew = the minimum value that should be sampled more
	f_skew = fraction of data set that should be sampled between [Cskew, Cmax]
	N = the number of A strands in the system
	M = the number of B strands in the system
	L = the number of C strands in the system
	It generates a .con file to be used with concentrations. The concentrations
	of A and B are 1 uM, and the concentrations of C are randomly generated 
	between [Cmin, Cmax], with f_skew of them having at least one Ci between
	[Cskew, Cmax].
	'''
	# Open a file of the right name
	with open(file_name + '.con', 'w') as con:
		for i in range(N+M):
			con.write('1e-6\n')

		# Determine whether to skew the concentrations this time
		is_skew = rand.uniform(0,1)
		if is_skew <= f_skew:
			# If the concentrations are skewed, pick which one has a value
			# between (Cskew, Cmax)
			which_skew = rand.choice(range(L))
			for l in range(L):
				# Generate the skewed value
				if l == which_skew:
					Ci = rand.uniform(Cskew, Cmax)
				# Generate the non-skewed values
				else:
					Ci = rand.uniform(Cmin, Cmax)
				con.write('{}\n'.format(Ci))
		else:
			# Generate the rest of the data set
			for l in range(L):
				Ci = rand.uniform(Cmin, Cskew)
				con.write('{}\n'.format(Ci))


def wombo_combo(strand_1, strand_2, L1, L2, spacer_len=0):
	# Create complementary strands for strand 1 and strand 2
	comp_1 = gen_complement(strand_1)
	comp_2 = gen_complement(strand_2)

	# Create the lead and tail sequences
	lead = comp_1[:L1]
	tail = comp_2[-L2:]

	# A list of bases
	bases = ['A', 'T', 'G', 'C']

	# Figure out what the complementary edge bases are
	edge_1 = comp_1[L1]
	edge_2 = comp_2[-L2-1]
	edges = [edge_1, edge_2]

	# Choose spacer edge bases that are not complementary
	spacers = []
	for i in range(2):
		base_copy = bases.copy()
		base_copy.remove(edges[i])
		spacer_edge = rand.choice(base_copy)
		spacers.append(spacer_edge)

	# Make the default spacer nothing
	spacer = ''

	# If 1 spacer is specified, make sure it isn't a complementary base to 
	# either of the strands
	if spacer_len == 1:
		bases.remove(edge_1)
		if edge_2 in bases:
			bases.remove(edge_2)
		spacer = rand.choice(bases)
	
	# If 2 spacers are specified, make sure they aren't complementary bases
	elif spacer_len == 2:
		spacer = spacers[0] + spacers[1]
		# spacer = 'hoho'
	
	# If more than 2 spacers are specified, make sure the edges aren't
	# complementary bases, and randomly choose the middle
	elif spacer_len > 2:
		spacer = spacers[0]
		for i in range(spacer_len-2):
			base = rand.choice(bases)
			spacer = spacer + base
		spacer = spacer + spacer[1]

	# Make the Combo
	combo = lead + spacer + tail

	return(combo)


# Nearest neighbor parameeters for generating orthogonal strands for 
# SantaLucia (1998).

# Nearest neighbor entropy parameters for 1M NaCl in cal / K mol
NN_S_dict = {'AA': -22.2, 'TT': -22.2, 'AT': -20.4, 'TA': -21.3, 'CA': -22.7, \
'TG': -22.7, 'GT': -22.4, 'AC': -22.4, 'CT': -21.0, 'AG': -21.0, 'GA': -22.2, \
'TC': -22.2, 'CG': -27.2, 'GC': -24.4, 'GG': -19.9, 'CC': -19.9}

NN_corr_S_dict = {'G': -2.8, 'C': -2.8, 'A': 4.1, 'T': 4.1, 'sym': -1.4}

# Nearest neighbor enthalpy parameters for 1M NaCl in cal / mol
NN_H_dict = {'AA': -7900.0, 'TT': -7900.0, 'AT': -7200.0, 'TA': -7200.0, \
'CA': -8500.0, 'TG': -8500.0, 'GT': -8400.0, 'AC': -8400.0, 'CT':-7800.0, \
'AG': -7800.0, 'GA': -8200.0, 'TC': -8200.0, 'CG': -10600.0, 'GC': -9800.0, \
'GG': -8000.0, 'CC': -8000.0}

NN_corr_H_dict = {'G': 100.0, 'C': 100.0, 'A': 2300.0, 'T': 2300.0, 'sym': 0.0}

# Functions for generating orthogonal strands for SantaLucia (1998) Tm
# calculation.

def check_self_complement(strand):
	'''This function takes a DNA sequence as a string and checks if it is 
	self-complementary. It returns True if it is, or False if it isn't.'''
	complement = gen_complement(strand)
	check = strand == complement
	return check

def calc_NN_S(strand):
	'''This function takes a DNA sequence as a string and calculates and returns
	the nearest neighbor entropy in 1M NaCl using the method described by 
	SantaLucia (1998).'''
	
	# Initiate the value of S
	S = 0
	
	# Add the terminal corrections
	S = S + NN_corr_S_dict.get(strand[0]) + NN_corr_S_dict.get(strand[-1])
	
	# Check if the self-complementary correction is needed
	if check_self_complement(strand):
		S = S + NN_corr_S_dict.get('sym')
	else:
		pass
	
	# Loop through the neighbor pairs to add their entropy
	for i in range(len(strand)-1):
		pair = strand[i:i+2]
		S = S + NN_S_dict.get(pair)
	return S

def calc_NN_H(strand):
	'''This function takes a DNA sequence as a string and calculates and returns
	the nearest neighbor enthalpy in 1M NaCl using the method described by 
	SantaLucia (1998).'''
	
	# Initiate the value of H
	H = 0
	
	# Add the terminal corrections
	H = H + NN_corr_H_dict.get(strand[0]) + NN_corr_H_dict.get(strand[-1])
	
	# Loop through the neighbor pairs to add their entropy
	for i in range(len(strand)-1):
		pair = strand[i:i+2]
		H = H + NN_H_dict.get(pair)
	return H

def calc_Tm(strand, conc, salt):
	'''This function takes a DNA sequence as a string (>13 bases), the strand 
	concentration (molar), and the salt concentration (Na+, molar). It computes 
	and returns the melting temperature Tm using the method described by 
	SantaLucia (1998).'''

	# Calculate the enthalpy (assumed to be salt independent)
	H = calc_NN_H(strand)
	
	# Calculate the entropy
	S = calc_NN_S(strand)

	# Figure out the number of phosphates to calculate the salt correction
	# This assumes that each base has a phosphate.
	n = len(strand)/2.0
	# Add the salt correction to the entropy
	S = S + 0.368 * n * np.log(salt)
	print(S)
	
	# Check if the strand is self-complementary to figure out what concentration
	# to use
	if check_self_complement:
		pass
	# If the strand is not self-complementary, assume it is at an equal 
	# concentration with its complement.
	else:
		conc = conc / 4.0

	# Calculate the melting temperature
	R = 1.987	# gas constant, cal / (K mol)
	Tm = H / (S + R * np.log(conc))
	return Tm

# Tm calculation via Sigma Aldrich method: 
# https://www.sigmaaldrich.com/technical-documents/articles/biology/oligos-melting-temp.html

# Nearest neighbor entropy parameters for Sigma Aldrich in cal / K mol
NN_S_dict_SA = {'AA': -24.0, 'TT': -24.0, 'AT': -23.9, 'TA': -16.9, 'CA': -12.9, \
'TG': -12.9, 'GT': -17.3, 'AC': -17.3, 'CT': -20.8, 'AG': -20.8, 'GA': -13.5, \
'TC': -13.5, 'CG': -27.8, 'GC': -26.7, 'GG': -26.6, 'CC': -26.6}

# Nearest neighbor enthalpy parameters for Sigma Aldrich in cal / mol
NN_H_dict_SA = {'AA': -9100.0, 'TT': -9100.0, 'AT': -8600.0, 'TA': -6000.0, \
'CA': -5800.0, 'TG': -5800.0, 'GT': -6500.0, 'AC': -6500.0, 'CT':-7800.0, \
'AG': -7800.0, 'GA': -5600.0, 'TC': -5600.0, 'CG': -11900.0, 'GC': -11100.0, \
'GG': -11000.0, 'CC': -11000.0}

def calc_NN_S_SA(strand):
	'''This function takes a DNA sequence as a string and calculates and returns
	the nearest neighbor entropy in 1M NaCl using the method described by
	Sigma Aldrich.
	https://www.sigmaaldrich.com/technical-documents/articles/biology/oligos-melting-temp.html'''
	
	# Initiate the value of S
	S = 0
	
	# Loop through the neighbor pairs to add their entropy
	for i in range(len(strand)-1):
		pair = strand[i:i+2]
		S = S + NN_S_dict_SA.get(pair)
	return S

def calc_NN_H_SA(strand):
	'''This function takes a DNA sequence as a string and calculates and returns
	the nearest neighbor enthalpy in 1M NaCl using the method described by 
	Sigma Aldrich.
	https://www.sigmaaldrich.com/technical-documents/articles/biology/oligos-melting-temp.html'''
	
	# Initiate the value of H
	H = 0
	
	# Loop through the neighbor pairs to add their entropy
	for i in range(len(strand)-1):
		pair = strand[i:i+2]
		H = H + NN_H_dict_SA.get(pair)
	return H

def calc_Tm_sigma_aldrich(strand, conc, salt):
	'''This function takes a DNA sequence as a string (>13 bases), the strand 
	concentration (molar), and the salt concentration (Na+, molar). It computes 
	and returns the melting temperature Tm using the method described by Sigma 
	Aldrich.
	https://www.sigmaaldrich.com/technical-documents/articles/biology/oligos-melting-temp.html'''

	# Calculate the enthalpy (assumed to be salt independent)
	H = calc_NN_H_SA(strand)
	
	# Calculate the entropy
	S = calc_NN_S_SA(strand)

	# Calculate the melting temperature
	R = 1.987	# gas constant, cal / (K mol)
	A = -10.8	# for helix initiation, cal / (K mol)
	Tm = H / (A + S + R * np.log(conc/4.0)) + 16.6 * np.log10(salt)
	return Tm

# 
def calc_Tm_short(strand, salt):
	'''This function takes a DNA sequence as a string (<14 bases), and the salt
	concentration (Na+, molar). The concentration of the DNA is assumed to be 
	0.5 uM. It computes and returns the melting tmperature Tm using the Marmur
	Doty method, with a salt correction. 
	Kibbe WA. 'OligoCalc: an online oligonucleotide properties calculator'. 
	(2007) Nucleic Acids Res. 35(web server issue)'''
	# Count A's and T's
	AT = strand.count('A') + strand.count('T')
	# Count G's and C's
	GC = strand.count('G') + strand.count('C')

	# Calculate Tm
	# Subtract 16.6 log(50e-3) to adjust for the assumed 50 mM
	Tm = 2.0 * AT + 4.0 * GC - 16.6 * np.log10(50e-3) + 16.6 * np.log10(salt)
	return Tm





# def generate_rand_input_file(file_name, N_strands, max_complex, length=0, \
# 	min=0, max=0):
# 	'''This function takes a file name, number of strands, and optional 
# 	arguments of length (strand length), min (minimum strand length), and max
# 	(maximum strand length). The output is a .in file with '''

# 	input = open(file_name + '.in', 'w')
# 	input.write('%d\n'.format(N_strands))

# 	input.close()

# 	print('{}.in file generated.' % file_name)
# 	pass


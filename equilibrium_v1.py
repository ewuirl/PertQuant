import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.special import comb
from timeit import default_timer as timer
	

def AA_interactions(AA_array, F, n):
	'''This function takes in an array of in-set double stranded concentrations, 
	such as AA, BB, or CC. It also takes in an value F = Ai - Af - {AB} - {AC}, 
	and an index n. This function computes the in-set interactions that 
	contribute to the final concentration of strand An. In other words for the 
	nth strand of A, this function computes the contribution of {AA} to the 
	relation 0 = Ai - Af - {AB} - {AC} - {AA} and returns 
	F = Ai - Af - {AB} - {AC} - {AA} as a float.'''
	if len(AA_array) > 0:
			F = F - AA_array[0]

	elif len(AA_array) > 1:
			
		if n < N-1:
			AnA_row_interactions = np.sum(AA_array[n])

		else:
			AnA_row_interactions = 0

		# Account for column interactions
		AnA_col_interactions = 0.0

		if n > 0:
			i = N - n
			for j in range(n):
				AnA_col_interactions = AnA_col_interactions + AA_array[j][-i]
		else:
			pass

		F = F - AnA_row_interactions - AnA_col_interactions

	else:
		pass
	return F

def calc_Af_equation(Ai_array, Af_array, N, M, L, AB_array, AC_array, AA_array):
	'''This function takes in:
	Ai_array, a 1xN array of initial concentrations for {A}
	Af_array, a 1xN array of final concentrations for {A}
	N, M, and L, the number of strands in {A}, {B}, and {C}
	AB_array, a NxM array of AB concentrations
	AC_array, a NxL array of AC concentrations
	AA_array, a 1x(N-1) array of 1x(N-1) to 1x1 arrays of AA concentrations.

	It computes and returns F_array, an 1xN array whose nth element is
	F = Ai - Af - {AB} - {AC} - {AA} for the nth strand of {A}.'''

	# Create an array to store F values in 
	F_array = np.empty(N)

	# Cycle through the A strands
	for n in range(N):

		# Calculate the contributions of {AB}
		# Handle 1 B case:
		if M == 1:
			AnB_interactions = AB_array[n]

		# Handle 1 A and multiple B case
		elif N == 1 and M > 1:
			AnB_interactions = np.sum(AB_array)

		# Handle multiple A and B case
		else:
			AnB_interactions = np.sum(AB_array[n,:])

		# Calculate the contributions of {AC}
		# Handle 1 C case:
		if L == 1:
			AnC_interactions = AC_array[n]

		# Handle 0 C case
		elif L == 0:
			AnC_interactions = 0

		# Handle 1 A and multiple C case
		elif N == 1 and L > 1:
			AnC_interactions = np.sum(AC_array)

		# Handle multiple A and C case
		else:
			AnC_interactions = np.sum(AC_array[n,:])
		F = Ai_array[n] - Af_array[n] - AnB_interactions - AnC_interactions

		# Account for AA interactions
		F = AA_interactions(AA_array, F, n)

		F_array[n] = F

	return F_array

def calc_Bf_equation(Bi_array, Bf_array, N, M, L, AB_array, BC_array, BB_array):
	'''This function takes in:
	Bi_array, a 1xM array of initial concentrations for {B}
	Bf_array, a 1xM array of final concentrations for {B}
	N, M, and L, the number of strands in {A}, {B}, and {C}
	AB_array, a NxM array of AB concentrations
	BC_array, a MxL array of BC concentrations
	BB_array, a 1x(M-1) array of 1x(M-1) to 1x1 arrays of BB concentrations.

	It computes and returns F_array, an 1xM array whose mth element is
	F = Bi - Bf - {AB} - {BC} - {BB} for the mth strand of {B}.'''
	# Create an array to store F values in 
	F_array = np.empty(M)

	# Cycle throug the B strands
	for m in range(M):

		# Calculate the contribution of {AB}
		# Handle 1 A case:
		if N == 1:
			ABm_interactions = AB_array[m]

		# Handle multiple A case
		else:
			ABm_interactions = np.sum(AB_array[:,m])

		# Calculate the contribution of {BC}
		# Handle 1 C case:
		if L == 1:
			BmC_interactions = BC_array[m]

		# Handle no C case
		elif L == 0:
			BmC_interactions = 0

		# Handle 1 B and multiple C case
		elif M == 1 and L > 1:
			BmC_interactions = np.sum(BC_array)

		# Handle multiple B and C case
		else:
			BmC_interactions = np.sum(BC_array[m,:])

		F = Bi_array[m] - Bf_array[m] - ABm_interactions - BmC_interactions

		# Account for BB interactions
		F = AA_interactions(BB_array, F, m)

		F_array[m] = F

	return F_array

def calc_Cf_equation(Ci_array, Cf_array, N, M, L, AC_array, BC_array, CC_array):
	'''This function takes in:
	Ci_array, a 1xL array of initial concentrations for {C}
	Cf_array, a 1xL array of final concentrations for {C}
	N, M, and L, the number of strands in {A}, {B}, and {C}
	AC_array, a NxL array of AC concentrations
	BC_array, a MxL array of BC concentrations
	CC_array, a 1x(L-1) array of 1x(L-1) to 1x1 arrays of CC concentrations.

	It computes and returns F_array, an 1xL array whose nth element is
	F = Ci - Cf - {AC} - {BC} - {CC} for the lth strand of {C}.'''
	# Create an array to store F values in 
	F_array = np.empty(L)

	# Cycle throug the B strands
	for l in range(L):

		# Calculate the contribution of {AC}
		# Handle 1 A case:
		if N == 1:
			ACl_interactions = AC_array[l]

		else:
			ACl_interactions = np.sum(AC_array[:,l])
			# print('AC %d = %f' % l, ACl_interactions)

		# Calculate the contribution of {BC}
		# Handle 1 B case:
		if M == 1:
			BCl_interactions = BC_array[l]

		# Handle no B case
		elif M == 0:
			BCl_interactions = 0

		# Handle multiple B and C case
		else:
			BCl_interactions = np.sum(BC_array[:,l])
			# print('BC %d = %f' % l, BCl_interactions)

		F = Ci_array[l] - Cf_array[l] - ACl_interactions - BCl_interactions

		# Account for CC interactions
		F = AA_interactions(CC_array, F, l)

		F_array[l] = F

	return F_array


def calc_AB_equations(K_array, AB_array, Af_array, Bf_array):
	'''This function takes in:
	K_array, an NxM array of KAB values
	AB_array, an NxM array of AB concentrations
	Af_array, a 1xN array of Af concentrations
	Bf_array, a 1xM array of Bf concentrations.

	It calculates and returns F_array, whose (n,m)th element is
	F = KAnBm * Afn * Bfm - AnBm.'''
	N = len(Af_array)
	M = len(Bf_array)

	# Handle N = M = 1 case
	if N == 1 and M == 1:
		# Make an array to hold all of the equations
		F_array = np.empty(1)

		# Calculate the F = 0 equation
		F = K_array[0]*Af_array[0]*Bf_array[0] - AB_array[0]
		F_array = np.array([F])

	# Handle the N = 1, M > 1 case
	elif N == 1 and M > 1:
		F_array = np.empty(K_array.shape)
		Af = Af_array[0]

		for m in range(M):
				Bf = Bf_array[m]
				F = K_array[m]*Af*Bf - AB_array[m]
				F_array[m] = F

	# Handle all other cases
	else:
		# Make an array to hold all of the equations
		F_array = np.empty(K_array.shape)

		# print('this is Karray')
		# print(K_array)
		# print('this is ABarray')
		# print(AB_array)
		# Cycle through {A}
		for n in range(N):
			Af = Af_array[n]

			# Cycle through {B} for each strand in {A}
			for m in range(M):
				Bf = Bf_array[m]
				F = K_array[n][m]*Af*Bf - AB_array[n][m]
				F_array[n][m] = F

	return F_array


def calc_AA_equations(K_array, AA_array, Af_array):
	'''This function takes in:
	K_array, a 1x(N-1) array of arrays (size 1x(N-1) to 1x1) of KAA values
	AA_array, a 1x(N-1) array of arrays (size 1x(N-1) to 1x1) of AA values
	Af_array, a 1xN array of Af concentrations

	It calculates and returns F_array, a 1x(N-1) array of arrays 
	(size 1x(N-1) to 1x1). The nth row of F_array is a matrix that contains 
	{Fnm} reactions not already contained in a previous row.
	Fnm = KAnAm * Afn * Afm - AnAm.'''

	N = len(Af_array)

	# Create an array to store arrays of F equations in
	F_array = np.empty(N-1,dtype=object)

	# Handle 1 AA case:
	if N == 2:
		F_array[0] = K_array[0]*Af_array[0]*Af_array[1] - AA_array[0]

	else: 
		# Iterate through the rows
		for n in range(N-1):
			# Figure out how many elements are in each row
			M = len(K_array[n])
			
			# Create an array to store Fn equations in
			F_row = np.empty(M)
			for m in range(M):
				Fnm = K_array[n][m] * Af_array[n] * Af_array[n+1+m] - AA_array[n][m]
				F_row[m] = Fnm

			F_array[n] = F_row

	return F_array


def reshape_AA_array(N, AA_flat_array, KAA_array, AA_array):
	print("reshape is being called")
	'''This function takes in:
	N = # of strands
	AA_flat_array = a 1xcomb(N,2) array of AA concentrations
	KAA_array = a 1x(N-1) array of KAA values whose rows go from 1x(N-1) to 1x1
	AA_array = an array with the same dimensions as KAA_array to store AA 
				concentrations in.

	It returns AA_array with the values of AA_flat_array stored inside.'''

	row_start = 0
	row_stop = 0
	for i in range(len(KAA_array)):
		# Figure out how long the row is
		row_length = len(KAA_array[i])
		row_stop = row_stop + row_length
		# Find the row
		AA_row = AA_flat_array[row_start:row_stop]
		# Insert the row
		AA_array[i] = AA_row

		row_start = row_stop

	return(AA_array)

def fill_AA_array(x_array, start, stop, N, KAA_array, AA_array):
	'''This function takes in:
	x_array = a 1 x (N + M + NxM + L + MxL + NxL + comb(N,2) + comb(M,2) 
	+ comb(L,2)) array consisting of the following: 
	[Af1, ...AfN - first N values, = A final concentrations
	 Bf1, ...BfM - next M values, = B final concentrations
	 AB_11, ...AB_NM - next NxM values, = AB concentrations
	 Cf1, ...CfL - next L values, = C final concentrations
	 BC_11, ...BC_ML - next MxL values, = BC final concentrations
	 AC_11, ...AC_NL - next NxL values, = AC final concentrations
	 AA_11, ...AA_NN - next comb(N,2) values, = AA final concentrations 
	 BB_11, ...BB_MM - next comb(M,2) values, = BB final concentrations 
	 CC_11, ...CC_LL - last comb(L,2) values, = CC final concentrations.]
	stop = an indice used to slice x_array
	N = # of strands
	AA_flat_array = a 1xcomb(N,2) array of AA concentrations
	KAA_array = a 1x(N-1) array of KAA values whose rows go from 1x(N-1) to 1x1
	AA_array = an array with the same dimensions as KAA_array to store AA 
				concentrations in.

	It returns a tuple of (start, stop, AA_array) in which start and stop have 
	been updated and AA_array has the values of AA_flat_array stored inside.'''
	
	# Handle the case where there is only one KAA value
	if len(KAA_array) == 1:
		start = stop 
		stop = int(stop + 1)
		AA_array[0] = x_array[start:stop]

	# Handle the case where there is more than one KAA value
	elif len(KAA_array) > 1:
		start = stop 
		stop = int(stop + comb(N,2))
		AA_flat_array = x_array[start:stop]
		AA_array = reshape_AA_array(N, AA_flat_array, KAA_array, AA_array)

	else: 
		pass

	return(start, stop, AA_array)

def Keq(x_array, *args):
	'''This function creates a system of equations to be used with 
	scipy.optimize.fsolve. This represents a DNA network with 3 sets of strands
	(A, B, C), which may have pairwise interactions. For example, the ith A 
	strand and the jth B strand may interact as Ai + Bj <-> ABij, with 
	equilibrium constant KABij. At equilibrium this interaction satisfies
	KABij[Af][Bf] - [AB] = 0. 

	x_array is a 1 x (N + M + NxM + L + MxL + NxL + comb(N,2) + comb(M,2) 
	+ comb(L,2)) array consisting of the following: 
	[Af1, ...AfN - first N values, = A final concentrations
	 Bf1, ...BfM - next M values, = B final concentrations
	 AB_11, ...AB_NM - next NxM values, = AB concentrations
	 Cf1, ...CfL - next L values, = C final concentrations
	 BC_11, ...BC_ML - next MxL values, = BC final concentrations
	 AC_11, ...AC_NL - next NxL values, = AC final concentrations
	 AA_11, ...AA_NN - next comb(N,2) values, = AA final concentrations 
	 BB_11, ...BB_MM - next comb(M,2) values, = BB final concentrations 
	 CC_11, ...CC_LL - last comb(L,2) values, = CC final concentrations.]

	The function also takes an array of "optional arguments".
	arguments is a numpy array:
	[[Ai1, ...AiN] (1xN), = A initial concentrations
	 [Bi1, ...BiM] (1xM), = B initial concentrations
	 [AB_11, ...AB_NM] (NxM), = KAB values
	 [Ci1, ...CiL] (1xL), = C initial concentrations
	 [KBC_11, ...KBC_ML] (MxL), = KBC values
	 [KAC_11, ...KAC_NL] (NxL), = KAC values
	 [KAA_11, ...KAA_NN] (1x(N-1) array of whose rows are arrays 1x(N-1) to 1x1), 
	 					= KAA values
	 [KBB_11, ...KBB_MM] (1x(M-1) array of whose rows are arrays 1x(M-1) to 1x1), 
	 					= KBB values
	 [KCC_11, ...KCC_LL] (1x(L-1) array of whose rows are arrays 1x(L-1) to 1x1) 
	 					= KCC values]

	Output: It returns F_array, (same dimensions as x_array) of a system of 
	equations = 0 to be	used with scipy.optimize.fsolve.'''

	# args = arguments[0]
	# print('These are the arguments')
	# print(arguments)

	# print('These are args')
	# print(args)

	# Unpack the arguments
	Ai_array = args[0].astype(np.float64)
	Bi_array = args[1].astype(np.float64)
	KAB_array = args[2].astype(np.float64)

	# Figure out dimensions of system
	N = len(Ai_array)
	M = len(Bi_array)

	# Unpack the guess
	# start and stop are slicing indices for x-array
	start = 0
	stop = N
	Af_array = x_array[start:stop].astype(np.float64)
	# print('I am Af_array')
	# print(Af_array)
	
	start = stop
	stop = stop + M
	Bf_array = x_array[start:stop].astype(np.float64)
	# print('I Am Bf_array')
	# print(Bf_array)

	start = stop 
	stop = int(stop + N*M)
	
	if N == 1 and M > 1:
		AB_array = x_array[start:stop].astype(np.float64)
	else:
		AB_array = x_array[start:stop].reshape(N, M).astype(np.float64)

	# print('I am AB_array')
	# print(AB_array)

	# Handle the {A} + {B} case where there are no in-set interactions	
	if len(args) == 3:
		AA_array = np.array([]).astype(np.float64)
		BB_array = np.array([]).astype(np.float64)

	# Handle the {A} + {B} case where there are {A} and/or {B} in-set interactions 
	if len(args) == 5:
		# Unpack the KAA and KBB arrays
		KAA_array = args[3]
		KBB_array = args[4]

		# Make arrays to store in-set interactions in
		AA_array = KAA_array.astype(np.float64)
		BB_array = KBB_array.astype(np.float64)

		# Fill in the AA array
		start, stop, AA_array = fill_AA_array(x_array, start, stop, N, \
			KAA_array, AA_array)

		# Fill in the BB array
		start, stop, BB_array = fill_AA_array(x_array, start, stop, M, \
			KBB_array, BB_array)

		# print('I am AA_array')
		# print(AA_array)
		# print('I am BB_array')
		# print(BB_array)

	# Handle the {A} + {B} + {C} case with no in-set interactions
	if len(args) > 5:
		Ci_array = args[3].astype(np.float64)
		KBC_array = args[4].astype(np.float64) 
		KAC_array = args[5].astype(np.float64)
		# Figure out the number of strands in {C}
		L = len(Ci_array)

		start = stop 
		stop = int(stop + L)
		Cf_array = x_array[start:stop].astype(np.float64)
		
		start = stop 
		stop = int(stop + M*L)
		if M == 1 and L > 1:
			BC_array = x_array[start:stop].astype(np.float64)
		else:
			BC_array = x_array[start:stop].reshape(M,L).astype(np.float64)
		
		start = stop 
		stop = int(stop + N*L)
		if N == 1 and L > 1:
			AC_array = x_array[start:stop].astype(np.float64)
		else:
			AC_array = x_array[start:stop].reshape(N,L).astype(np.float64)

		# Assign empty arrays for in-set interaction arrays
		AA_array = np.array([]).astype(np.float64)
		BB_array = np.array([]).astype(np.float64)
		CC_array = np.array([]).astype(np.float64)
		
		# print('I am Cf_array')
		# print(Cf_array)
		# print('I am BC_array')
		# print(BC_array)
		# print('I am AC_array')
		# print(AC_array)
		# print('I am AA_array')
		# print(AA_array)
		# print('I am BB_array')
		# print(BB_array)
		# print('I am CC_array')
		# print(CC_array)
	
	else:
		L = 0
		BC_array = np.array([]).astype(np.float64)
		AC_array = np.array([]).astype(np.float64)

	# Handle the {A} + {B} + {C} case with self-interactions
	if len(args) == 9:
		# Unpack the K_arrays
		KAA_array = args[6].astype(np.float64)
		KBB_array = args[7].astype(np.float64)
		KCC_array = args[8].astype(np.float64)

		# Make arrays to store in-set interaction concentrations in 
		AA_array = KAA_array.astype(np.float64)
		BB_array = KBB_array.astype(np.float64)
		CC_array = KCC_array.astype(np.float64)

		# Fill in the AA array
		start, stop, AA_array = fill_AA_array(x_array, start, stop, N, \
			KAA_array, AA_array)

		# Fill in the BB array
		start, stop, BB_array = fill_AA_array(x_array, start, stop, M, \
			KBB_array, BB_array)

		# Fill in the CC array
		start, stop, CC_array = fill_AA_array(x_array, start, stop, L, \
			KCC_array, CC_array)

		# print('I am AA_array')
		# print(AA_array)
		# print('I am BB_array')
		# print(BB_array)
		# print('I am CC_array')
		# print(CC_array)
	else:
		pass

	
	
	# Generate F_array
	F_array = np.empty(len(args),dtype=object)

	# Calculate the Af equations
	F_Af_array = calc_Af_equation(Ai_array, Af_array, N, M, L, AB_array, \
		AC_array, AA_array)

	# Calculate the Bf equations
	F_Bf_array = calc_Bf_equation(Bi_array, Bf_array, N, M, L, AB_array, \
		BC_array, BB_array)

	# Calculate the AB equations
	F_AB_array = calc_AB_equations(KAB_array, AB_array, Af_array, Bf_array)

	# Enter them into F_array
	# f is an index for F to determine where in F to put an array
	f = 0
	F_array[f] = F_Af_array
	f = f+1
	F_array[f] = F_Bf_array	
	f = f+1
	F_array[f] = F_AB_array

	# Handle the A + B + C cases
	if len(args) >= 6:
		# Calculate the Cf equations
		F_Cf_array = calc_Cf_equation(Ci_array, Cf_array, N, M, L, AC_array, \
			BC_array, CC_array)

		# Calculate the BC equations
		F_BC_array = calc_AB_equations(KBC_array, BC_array, Bf_array, Cf_array)

		# Calculate the AC equations
		F_AC_array = calc_AB_equations(KAC_array, AC_array, Af_array, Cf_array)

		# Enter them into F_array
		f = f+1
		F_array[f] = F_Cf_array
		f = f+1
		F_array[f] = F_BC_array
		f = f+1
		F_array[f] = F_AC_array


	# Handle in-set interactions
	if len(args) == 5 or len(args) == 9:
		# Calculate the AA equations if there are any interactions
		if len(AA_array) != 0:
			F_AA_array = calc_AA_equations(KAA_array, AA_array, Af_array)
			# Enter the equations into F_array
			f = f+1
			F_array[f] = F_AA_array

		else:
			pass

		# Calculate the BB equations if there are any interactions
		if len(BB_array) != 0:
			F_BB_array = calc_AA_equations(KBB_array, BB_array, Bf_array)
			f = f+1
			F_array[f] = F_BB_array

		else:
			pass
	else:
		pass

	if len(args) == 9:
		# Calculate the CC equations if there are any interactions
		if len(CC_array) != 0:
			F_CC_array = calc_AA_equations(KCC_array, CC_array, Cf_array)
			f = f+1
			F_array[f] = F_CC_array
	else:
		pass

	# print('this is F_array')
	# print(F_array)

	# Flatten F_array
	F_flat_array = np.array([])
	for i in range(len(F_array)):
		F_row = F_array[i]
		F_flat_array = np.concatenate((F_flat_array, F_row),axis=None).astype(np.float64)

	# Remove 'none' elements
	F_flat_array = F_flat_array[F_flat_array != np.array(None)]

	return(F_flat_array)


def calc_Am_array(Ai_array, N, M, AB_array):
	'''This function takes in:
	Ai_array, a 1xN array of initial concentrations for {A}
	N and M, the number of strands in {A} and {B}
	AB_array, a NxM array of AB concentrations

	It computes and returns Am_array, an 1xN array whose nth element is
	Am = Ai - {AB} for the nth strand of {A}.'''

	# Create an array to store Am values in 
	Am_array = np.empty(N)

	# Cycle through the A strands
	for n in range(N):

		# Calculate the contributions of {AB}
		# Handle 1 B case:
		if M == 1:
			AnB_interactions = AB_array[n]

		# Handle 1 A and multiple B case
		elif N == 1 and M > 1:
			AnB_interactions = np.sum(AB_array)

		# Handle multiple A and B case
		else:
			AnB_interactions = np.sum(AB_array[n,:])

		Am = Ai_array[n] - AnB_interactions

		Am_array[n] = Am

	return Am_array


def calc_A_measured(x_array, *args):
	'''This function is to be used with x_array, the output of 
	scipy.optimize.fsolve used on Keq. We are interested in a system where we 
	remove all {B} strands and anything attached to {B} and measure the
	concentrations of {A}. This function calculates the concentrations of {A} 
	that should be measured.

	x_array is a 1 x (N + M + NxM + L + MxL + NxL + comb(N,2) + comb(M,2) 
	+ comb(L,2)) array consisting of the following: 
	[Af1, ...AfN - first N values, = A final concentrations
	 Bf1, ...BfM - next M values, = B final concentrations
	 AB_11, ...AB_NM - next NxM values, = AB concentrations
	 Cf1, ...CfL - next L values, = C final concentrations
	 BC_11, ...BC_ML - next MxL values, = BC final concentrations
	 AC_11, ...AC_NL - next NxL values, = AC final concentrations
	 AA_11, ...AA_NN - next comb(N,2) values, = AA final concentrations 
	 BB_11, ...BB_MM - next comb(M,2) values, = BB final concentrations 
	 CC_11, ...CC_LL - last comb(L,2) values, = CC final concentrations.]

	The function also takes an array of "optional arguments".
	arguments is a numpy array:
	[[Ai1, ...AiN] (1xN), = A initial concentrations
	 [Bi1, ...BiM] (1xM), = B initial concentrations
	 [AB_11, ...AB_NM] (NxM), = KAB values
	 [Ci1, ...CiL] (1xL), = C initial concentrations
	 [KBC_11, ...KBC_ML] (MxL), = KBC values
	 [KAC_11, ...KAC_NL] (NxL), = KAC values
	 [KAA_11, ...KAA_NN] (1x(N-1) array of whose rows are arrays 1x(N-1) to 1x1), 
	 					= KAA values
	 [KBB_11, ...KBB_MM] (1x(M-1) array of whose rows are arrays 1x(M-1) to 1x1), 
	 					= KBB values
	 [KCC_11, ...KCC_LL] (1x(L-1) array of whose rows are arrays 1x(L-1) to 1x1) 
	 					= KCC values]

	Output: It returns Am_array, a 1xN array of whose nth element is 
	Am = Ai - {AB} for the nth strand in {A}.'''

	# Unpack the arguments
	Ai_array = args[0].astype(np.float64)
	Bi_array = args[1].astype(np.float64)

	# Figure out dimensions of system
	N = len(Ai_array)
	M = len(Bi_array)

	# Unpack the guess
	# start and stop are slicing indices for x-array
	start = 0
	stop = N
	Af_array = x_array[start:stop].astype(np.float64)
	# print('I am Af_array')
	# print(Af_array)
	
	start = stop
	stop = stop + M
	Bf_array = x_array[start:stop].astype(np.float64)
	# print('I Am Bf_array')
	# print(Bf_array)

	start = stop 
	stop = int(stop + N*M)
	
	if N == 1 and M > 1:
		AB_array = x_array[start:stop].astype(np.float64)
	else:
		AB_array = x_array[start:stop].reshape(N, M).astype(np.float64)

	# print('I am AB_array')
	# print(AB_array)

	# Calculate the Af equations
	Am_array = calc_Am_array(Ai_array, N, M, AB_array)

	return(Am_array)

# # Test A + B case for Keq
# Af_array = np.array([1])
# Bf_array = np.array([5])
# AB_array = np.array([3])
# x_array = np.concatenate((Af_array, Bf_array, AB_array),axis=None)
# print('this is x_array')
# print(x_array)

# Ai_array = np.array([4])
# Bi_array = np.array([7])
# KAB_array = np.array([10000])
# args = (Ai_array, Bi_array, KAB_array)

# trial = Keq(x_array, *args)
# print(trial)
# print(float(Ai_array[0]-Af_array[0]-AB_array[0]))
# print(float(Bi_array[0]-Bf_array[0]-AB_array[0]))
# print(float(KAB_array[0]*Af_array[0]*Bf_array[0]-AB_array[0]))

# # Test A + B + C case for Keq
# Af_array = np.array([1])
# Bf_array = np.array([5])
# Cf_array = np.array([7])
# AB_array = np.array([3])
# BC_array = np.array([2])
# AC_array = np.array([8])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, Cf_array, BC_array, AC_array),axis=None)
# print('this is x_array')
# print(x_array)

# Ai_array = np.array([4])
# Bi_array = np.array([7])
# Ci_array = np.array([9])
# KAB_array = np.array([10000])
# KBC_array = np.array([241])
# KAC_array = np.array([35.78])
# args = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(float(Ai_array[0]-Af_array[0]-AB_array[0]-AC_array[0]))
# print(float(Bi_array[0]-Bf_array[0]-AB_array[0]-BC_array[0]))
# print(float(KAB_array[0]*Af_array[0]*Bf_array[0]-AB_array[0]))
# print(float(Ci_array[0]-Cf_array[0]-BC_array[0]-AC_array[0]))
# print(float(KBC_array[0]*Cf_array[0]*Bf_array[0]-BC_array[0]))
# print(float(KAC_array[0]*Af_array[0]*Cf_array[0]-AC_array[0]))


# # Test {A} + B case with no interactions
# Af_array = np.array([1,2])
# Bf_array = np.array([5])
# AB_array = np.array([[3],[8]])
# x_array = np.concatenate((Af_array, Bf_array, AB_array),axis=None)

# Ai_array = np.array([4,6])
# Bi_array = np.array([7])
# KAB_array = np.array([[10000],[45]])
# args = (Ai_array, Bi_array, KAB_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(Ai_array[0]-Af_array[0]-AB_array[0])
# print(Ai_array[1]-Af_array[1]-AB_array[1])
# print(Bi_array[0]-Bf_array[0]-AB_array[0][0]-AB_array[1])
# print(KAB_array[0]*Af_array[0]*Bf_array[0]-AB_array[0])
# print(KAB_array[1]*Af_array[1]*Bf_array[0]-AB_array[1])


# # Test {A} + B case with interactions
# Af_array = np.array([1,2])
# Bf_array = np.array([5])
# AB_array = np.array([[3],[8]])
# AA_array = np.array([3.2])
# BB_array = np.array([])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, AA_array, BB_array),axis=None)
# print('this is x_array')
# print(x_array)

# Ai_array = np.array([4,6])
# Bi_array = np.array([7])
# KAB_array = np.array([[10000],[45]])
# KAA_array = np.array([35])
# KBB_array = np.array([])
# args = (Ai_array, Bi_array, KAB_array, KAA_array, KBB_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(float(Ai_array[0]-Af_array[0]-AB_array[0]-AA_array[0]))
# print(float(Ai_array[1]-Af_array[1]-AB_array[1]-AA_array[0]))
# print(float(Bi_array[0]-Bf_array[0]-AB_array[0]-AB_array[1]))
# print(float(KAB_array[0]*Af_array[0]*Bf_array[0]-AB_array[0]))
# print(float(KAB_array[1]*Af_array[1]*Bf_array[0]-AB_array[1]))
# print(float(KAA_array[0]*Af_array[0]*Af_array[1]-AA_array[0]))


# # Test A + {B} + case with no interactions
# Bf_array = np.array([1,2])
# Af_array = np.array([5])
# AB_array = np.array([3,8])
# x_array = np.concatenate((Af_array, Bf_array, AB_array),axis=None)
# print(x_array)

# Bi_array = np.array([4,6])
# Ai_array = np.array([7])
# KAB_array = np.array([10000,45])
# args = (Ai_array, Bi_array, KAB_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(Ai_array[0]-Af_array[0]-AB_array[0]-AB_array[1])
# print(Bi_array[0]-Bf_array[0]-AB_array[0])
# print(Bi_array[1]-Bf_array[1]-AB_array[1])
# print(KAB_array[0]*Bf_array[0]*Af_array[0]-AB_array[0])
# print(KAB_array[1]*Bf_array[1]*Af_array[0]-AB_array[1])


# # Test A + {B} case with interactions
# Bf_array = np.array([1,2])
# Af_array = np.array([5])
# AB_array = np.array([3,8])
# AA_array = np.array([])
# BB_array = np.array([1.6])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, AA_array, BB_array),axis=None)
# print(x_array)

# Bi_array = np.array([4,6])
# Ai_array = np.array([7])
# KAB_array = np.array([10000,45])
# KAA_array = np.array([])
# KBB_array = np.array([43])
# args = (Ai_array, Bi_array, KAB_array, KAA_array, KBB_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(Ai_array[0]-Af_array[0]-AB_array[0]-AB_array[1])
# print(Bi_array[0]-Bf_array[0]-AB_array[0]-BB_array[0])
# print(Bi_array[1]-Bf_array[1]-AB_array[1]-BB_array[0])
# print(KAB_array[0]*Bf_array[0]*Af_array[0]-AB_array[0])
# print(KAB_array[1]*Bf_array[1]*Af_array[0]-AB_array[1])
# print(KBB_array[0]*Bf_array[1]*Bf_array[0]-BB_array[0])


# # Test {A} + {B} case without interactions
# Bf_array = np.array([1,2])
# Af_array = np.array([5,2.1])
# AB_array = np.array([[3,8],[12,30]])
# x_array = np.concatenate((Af_array, Bf_array, AB_array),axis=None)
# print(x_array)

# Bi_array = np.array([4,6])
# Ai_array = np.array([7,3.5])
# KAB_array = np.array([[10000,45],[61,85]])
# args = (Ai_array, Bi_array, KAB_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(Ai_array[0]-Af_array[0]-AB_array[0][0]-AB_array[0][1])
# print(Ai_array[1]-Af_array[1]-AB_array[1][0]-AB_array[1][1])
# print(Bi_array[0]-Bf_array[0]-AB_array[0][0]-AB_array[1][0])
# print(Bi_array[1]-Bf_array[1]-AB_array[0][1]-AB_array[1][1])
# print(KAB_array[0][0]*Bf_array[0]*Af_array[0]-AB_array[0][0])
# print(KAB_array[0][1]*Bf_array[1]*Af_array[0]-AB_array[0][1])
# print(KAB_array[1][0]*Bf_array[0]*Af_array[1]-AB_array[1][0])
# print(KAB_array[1][1]*Bf_array[1]*Af_array[1]-AB_array[1][1])


# # Test {A} + {B} case with interactions
# Bf_array = np.array([1,2])
# Af_array = np.array([5,2.1])
# AB_array = np.array([[3,8],[12,30]])
# AA_array = np.array([3.2])
# BB_array = np.array([1.6])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, AA_array, BB_array),axis=None)
# print(x_array)

# Bi_array = np.array([4,6])
# Ai_array = np.array([7,3.5])
# KAB_array = np.array([[10000,45],[61,85]])
# KAA_array = np.array([35])
# KBB_array = np.array([43])
# args = (Ai_array, Bi_array, KAB_array, KAA_array, KBB_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(Ai_array[0]-Af_array[0]-AB_array[0][0]-AB_array[0][1]-AA_array[0])
# print(Ai_array[1]-Af_array[1]-AB_array[1][0]-AB_array[1][1]-AA_array[0])
# print(Bi_array[0]-Bf_array[0]-AB_array[0][0]-AB_array[1][0]-BB_array[0])
# print(Bi_array[1]-Bf_array[1]-AB_array[0][1]-AB_array[1][1]-BB_array[0])
# print(KAB_array[0][0]*Bf_array[0]*Af_array[0]-AB_array[0][0])
# print(KAB_array[0][1]*Bf_array[1]*Af_array[0]-AB_array[0][1])
# print(KAB_array[1][0]*Bf_array[0]*Af_array[1]-AB_array[1][0])
# print(KAB_array[1][1]*Bf_array[1]*Af_array[1]-AB_array[1][1])
# print(KAA_array[0]*Af_array[1]*Af_array[0]-AA_array[0])
# print(KBB_array[0]*Bf_array[1]*Bf_array[0]-BB_array[0])


# # Test {A} + {B} + {C} case without interactions
# Bf_array = np.array([1,2])
# Af_array = np.array([5,2.1])
# AB_array = np.array([[3,8],[12,30]])
# Cf_array = np.array([16, 5.9])
# BC_array = np.array([[4,5.1],[11,0.9]])
# AC_array = np.array([[20,10],[2,0.43]])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, Cf_array, BC_array, AC_array),axis=None)
# print(x_array)

# Bi_array = np.array([4,6])
# Ai_array = np.array([7,3.5])
# KAB_array = np.array([[10000,45],[61,85]])
# Ci_array = np.array([16, 5.9])
# KBC_array = np.array([[3,2.7],[13,14]])
# KAC_array = np.array([[2.3,1.8],[8.5,7.2]])
# args = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(Ai_array[0]-Af_array[0]-AB_array[0][0]-AB_array[0][1]-AC_array[0][0]-AC_array[0][1])
# print(Ai_array[1]-Af_array[1]-AB_array[1][0]-AB_array[1][1]-AC_array[1][0]-AC_array[1][1])
# print(Bi_array[0]-Bf_array[0]-AB_array[0][0]-AB_array[1][0]-BC_array[0][0]-BC_array[0][1])
# print(Bi_array[1]-Bf_array[1]-AB_array[0][1]-AB_array[1][1]-BC_array[1][0]-BC_array[1][1])
# print(KAB_array[0][0]*Bf_array[0]*Af_array[0]-AB_array[0][0])
# print(KAB_array[0][1]*Bf_array[1]*Af_array[0]-AB_array[0][1])
# print(KAB_array[1][0]*Bf_array[0]*Af_array[1]-AB_array[1][0])
# print(KAB_array[1][1]*Bf_array[1]*Af_array[1]-AB_array[1][1])
# print(Ci_array[0]-Cf_array[0]-AC_array[0][0]-AC_array[1][0]-BC_array[0][0]-BC_array[1][0])
# print(Ci_array[1]-Cf_array[1]-AC_array[0][1]-AC_array[1][1]-BC_array[0][1]-BC_array[1][1])
# print(KBC_array[0][0]*Cf_array[0]*Bf_array[0]-BC_array[0][0])
# print(KBC_array[0][1]*Cf_array[1]*Bf_array[0]-BC_array[0][1])
# print(KBC_array[1][0]*Cf_array[0]*Bf_array[1]-BC_array[1][0])
# print(KBC_array[1][1]*Cf_array[1]*Bf_array[1]-BC_array[1][1])
# print(KAC_array[0][0]*Cf_array[0]*Af_array[0]-AC_array[0][0])
# print(KAC_array[0][1]*Cf_array[1]*Af_array[0]-AC_array[0][1])
# print(KAC_array[1][0]*Cf_array[0]*Af_array[1]-AC_array[1][0])
# print(KAC_array[1][1]*Cf_array[1]*Af_array[1]-AC_array[1][1])


# # Test {A} + {B} + {C} case with interactions
# Bf_array = np.array([1,2])
# Af_array = np.array([5,2.1])
# AB_array = np.array([[3,8],[12,30]])
# Cf_array = np.array([16, 5.9])
# BC_array = np.array([[4,5.1],[11,0.9]])
# AC_array = np.array([[20,10],[2,0.43]])
# AA_array = np.array([3.2])
# BB_array = np.array([1.6])
# CC_array = np.array([4.1])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, Cf_array, BC_array, \
# 	AC_array, AA_array, BB_array, CC_array),axis=None)
# print(x_array)

# Bi_array = np.array([4,6])
# Ai_array = np.array([7,3.5])
# KAB_array = np.array([[10000,45],[61,85]])
# Ci_array = np.array([16, 5.9])
# KBC_array = np.array([[3,2.7],[13,14]])
# KAC_array = np.array([[2.3,1.8],[8.5,7.2]])
# KAA_array = np.array([35])
# KBB_array = np.array([43])
# KCC_array = np.array([33])
# args = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array, \
# 	KAA_array, KBB_array, KCC_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(Ai_array[0]-Af_array[0]-AB_array[0][0]-AB_array[0][1]-AC_array[0][0]-AC_array[0][1]-AA_array[0])
# print(Ai_array[1]-Af_array[1]-AB_array[1][0]-AB_array[1][1]-AC_array[1][0]-AC_array[1][1]-AA_array[0])
# print(Bi_array[0]-Bf_array[0]-AB_array[0][0]-AB_array[1][0]-BC_array[0][0]-BC_array[0][1]-BB_array[0])
# print(Bi_array[1]-Bf_array[1]-AB_array[0][1]-AB_array[1][1]-BC_array[1][0]-BC_array[1][1]-BB_array[0])
# print(KAB_array[0][0]*Bf_array[0]*Af_array[0]-AB_array[0][0])
# print(KAB_array[0][1]*Bf_array[1]*Af_array[0]-AB_array[0][1])
# print(KAB_array[1][0]*Bf_array[0]*Af_array[1]-AB_array[1][0])
# print(KAB_array[1][1]*Bf_array[1]*Af_array[1]-AB_array[1][1])
# print(Ci_array[0]-Cf_array[0]-AC_array[0][0]-AC_array[1][0]-BC_array[0][0]-BC_array[1][0]-CC_array)
# print(Ci_array[1]-Cf_array[1]-AC_array[0][1]-AC_array[1][1]-BC_array[0][1]-BC_array[1][1]-CC_array)
# print(KBC_array[0][0]*Cf_array[0]*Bf_array[0]-BC_array[0][0])
# print(KBC_array[0][1]*Cf_array[1]*Bf_array[0]-BC_array[0][1])
# print(KBC_array[1][0]*Cf_array[0]*Bf_array[1]-BC_array[1][0])
# print(KBC_array[1][1]*Cf_array[1]*Bf_array[1]-BC_array[1][1])
# print(KAC_array[0][0]*Cf_array[0]*Af_array[0]-AC_array[0][0])
# print(KAC_array[0][1]*Cf_array[1]*Af_array[0]-AC_array[0][1])
# print(KAC_array[1][0]*Cf_array[0]*Af_array[1]-AC_array[1][0])
# print(KAC_array[1][1]*Cf_array[1]*Af_array[1]-AC_array[1][1])
# print(KAA_array[0]*Af_array[1]*Af_array[0]-AA_array[0])
# print(KBB_array[0]*Bf_array[1]*Bf_array[0]-BB_array[0])
# print(KCC_array[0]*Cf_array[1]*Cf_array[0]-CC_array[0])


# # Test {A} + {B} + C case without interactions
# Bf_array = np.array([1,2])
# Af_array = np.array([5,2.1])
# AB_array = np.array([[3,8],[12,30]])
# Cf_array = np.array([16])
# BC_array = np.array([[4],[11]])
# AC_array = np.array([[20],[2]])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, Cf_array, BC_array, \
# 	AC_array),axis=None)
# print(x_array)

# Bi_array = np.array([4,6])
# Ai_array = np.array([7,3.5])
# KAB_array = np.array([[10000,45],[61,85]])
# Ci_array = np.array([16])
# KBC_array = np.array([[3],[13]])
# KAC_array = np.array([[2.3],[8.5]])
# args = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(Ai_array[0]-Af_array[0]-AB_array[0][0]-AB_array[0][1]-AC_array[0])
# print(Ai_array[1]-Af_array[1]-AB_array[1][0]-AB_array[1][1]-AC_array[1])
# print(Bi_array[0]-Bf_array[0]-AB_array[0][0]-AB_array[1][0]-BC_array[0])
# print(Bi_array[1]-Bf_array[1]-AB_array[0][1]-AB_array[1][1]-BC_array[1])
# print(KAB_array[0][0]*Bf_array[0]*Af_array[0]-AB_array[0][0])
# print(KAB_array[0][1]*Bf_array[1]*Af_array[0]-AB_array[0][1])
# print(KAB_array[1][0]*Bf_array[0]*Af_array[1]-AB_array[1][0])
# print(KAB_array[1][1]*Bf_array[1]*Af_array[1]-AB_array[1][1])
# print(Ci_array[0]-Cf_array[0]-AC_array[0][0]-AC_array[1][0]-BC_array[0][0]-BC_array[1][0])
# print(KBC_array[0][0]*Cf_array[0]*Bf_array[0]-BC_array[0][0])
# print(KBC_array[1][0]*Cf_array[0]*Bf_array[1]-BC_array[1][0])
# print(KAC_array[0][0]*Cf_array[0]*Af_array[0]-AC_array[0][0])
# print(KAC_array[1][0]*Cf_array[0]*Af_array[1]-AC_array[1][0])


# # Test {A} + {B} + C case with interactions
# Bf_array = np.array([1,2])
# Af_array = np.array([5,2.1])
# AB_array = np.array([[3,8],[12,30]])
# Cf_array = np.array([16])
# BC_array = np.array([[4],[11]])
# AC_array = np.array([[20],[2]])
# AA_array = np.array([3.2])
# BB_array = np.array([1.6])
# CC_array = np.array([])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, Cf_array, BC_array, \
# 	AC_array, AA_array, BB_array, CC_array),axis=None)
# print(x_array)

# Bi_array = np.array([4,6])
# Ai_array = np.array([7,3.5])
# KAB_array = np.array([[10000,45],[61,85]])
# Ci_array = np.array([16])
# KBC_array = np.array([[3],[13]])
# KAC_array = np.array([[2.3],[8.5]])
# KAA_array = np.array([35])
# KBB_array = np.array([43])
# KCC_array = np.array([])
# args = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array, \
# 	KAA_array, KBB_array, KCC_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(Ai_array[0]-Af_array[0]-AB_array[0][0]-AB_array[0][1]-AC_array[0]-AA_array[0])
# print(Ai_array[1]-Af_array[1]-AB_array[1][0]-AB_array[1][1]-AC_array[1]-AA_array[0])
# print(Bi_array[0]-Bf_array[0]-AB_array[0][0]-AB_array[1][0]-BC_array[0]-BB_array[0])
# print(Bi_array[1]-Bf_array[1]-AB_array[0][1]-AB_array[1][1]-BC_array[1]-BB_array[0])
# print(KAB_array[0][0]*Bf_array[0]*Af_array[0]-AB_array[0][0])
# print(KAB_array[0][1]*Bf_array[1]*Af_array[0]-AB_array[0][1])
# print(KAB_array[1][0]*Bf_array[0]*Af_array[1]-AB_array[1][0])
# print(KAB_array[1][1]*Bf_array[1]*Af_array[1]-AB_array[1][1])
# print(Ci_array[0]-Cf_array[0]-AC_array[0][0]-AC_array[1][0]-BC_array[0][0]-BC_array[1][0])
# print(KBC_array[0][0]*Cf_array[0]*Bf_array[0]-BC_array[0][0])
# print(KBC_array[1][0]*Cf_array[0]*Bf_array[1]-BC_array[1][0])
# print(KAC_array[0][0]*Cf_array[0]*Af_array[0]-AC_array[0][0])
# print(KAC_array[1][0]*Cf_array[0]*Af_array[1]-AC_array[1][0])
# print(KAA_array[0]*Af_array[1]*Af_array[0]-AA_array[0])
# print(KBB_array[0]*Bf_array[1]*Bf_array[0]-BB_array[0])


# # Test {A} + B + {C} case without interactions
# Af_array = np.array([5,2.1])
# Bf_array = np.array([16])
# AB_array = np.array([[20],[2]])
# Cf_array = np.array([1,2])
# BC_array = np.array([4,11])
# AC_array = np.array([[3,8],[12,30]])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, Cf_array, BC_array, \
# 	AC_array),axis=None)
# print(x_array)


# Ai_array = np.array([7,3.5])
# Bi_array = np.array([16])
# KAB_array = np.array([[2.3],[8.5]])
# Ci_array = np.array([4,6])
# KBC_array = np.array([3,13])
# KAC_array = np.array([[10000,45],[61,85]])
# args = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(Ai_array[0]-Af_array[0]-AC_array[0][0]-AC_array[0][1]-AB_array[0][0])
# print(Ai_array[1]-Af_array[1]-AC_array[1][0]-AC_array[1][1]-AB_array[1][0])
# print(Bi_array[0]-Bf_array[0]-AB_array[0][0]-AB_array[1][0]-BC_array[0]-BC_array[1])
# print(KAB_array[0][0]*Bf_array[0]*Af_array[0]-AB_array[0][0])
# print(KAB_array[1][0]*Bf_array[0]*Af_array[1]-AB_array[1][0])
# print(Ci_array[0]-Cf_array[0]-AC_array[0][0]-AC_array[1][0]-BC_array[0])
# print(Ci_array[1]-Cf_array[1]-AC_array[0][1]-AC_array[1][1]-BC_array[1])
# print(KBC_array[0]*Bf_array[0]*Cf_array[0]-BC_array[0])
# print(KBC_array[1]*Bf_array[0]*Cf_array[1]-BC_array[1])
# print(KAC_array[0][0]*Af_array[0]*Cf_array[0]-AC_array[0][0])
# print(KAC_array[0][1]*Af_array[0]*Cf_array[1]-AC_array[0][1])
# print(KAC_array[1][0]*Af_array[1]*Cf_array[0]-AC_array[1][0])
# print(KAC_array[1][1]*Af_array[1]*Cf_array[1]-AC_array[1][1])


# # Test {A} + B + {C} case with interactions
# Af_array = np.array([5,2.1])
# Bf_array = np.array([16])
# AB_array = np.array([[20],[2]])
# Cf_array = np.array([1,2])
# BC_array = np.array([4,11])
# AC_array = np.array([[3,8],[12,30]])
# AA_array = np.array([3.2])
# BB_array = np.array([])
# CC_array = np.array([1.6])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, Cf_array, BC_array, \
# 	AC_array, AA_array, BB_array, CC_array),axis=None)
# print(x_array)


# Ai_array = np.array([7,3.5])
# Bi_array = np.array([16])
# KAB_array = np.array([[2.3],[8.5]])
# Ci_array = np.array([4,6])
# KBC_array = np.array([3,13])
# KAC_array = np.array([[10000,45],[61,85]])
# KAA_array = np.array([35])
# KBB_array = np.array([])
# KCC_array = np.array([43])
# args = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array, \
# 	KAA_array, KBB_array, KCC_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(Ai_array[0]-Af_array[0]-AC_array[0][0]-AC_array[0][1]-AB_array[0][0]-AA_array[0])
# print(Ai_array[1]-Af_array[1]-AC_array[1][0]-AC_array[1][1]-AB_array[1][0]-AA_array[0])
# print(Bi_array[0]-Bf_array[0]-AB_array[0][0]-AB_array[1][0]-BC_array[0]-BC_array[1])
# print(KAB_array[0][0]*Bf_array[0]*Af_array[0]-AB_array[0][0])
# print(KAB_array[1][0]*Bf_array[0]*Af_array[1]-AB_array[1][0])
# print(Ci_array[0]-Cf_array[0]-AC_array[0][0]-AC_array[1][0]-BC_array[0]-CC_array[0])
# print(Ci_array[1]-Cf_array[1]-AC_array[0][1]-AC_array[1][1]-BC_array[1]-CC_array[0])
# print(KBC_array[0]*Bf_array[0]*Cf_array[0]-BC_array[0])
# print(KBC_array[1]*Bf_array[0]*Cf_array[1]-BC_array[1])
# print(KAC_array[0][0]*Af_array[0]*Cf_array[0]-AC_array[0][0])
# print(KAC_array[0][1]*Af_array[0]*Cf_array[1]-AC_array[0][1])
# print(KAC_array[1][0]*Af_array[1]*Cf_array[0]-AC_array[1][0])
# print(KAC_array[1][1]*Af_array[1]*Cf_array[1]-AC_array[1][1])
# print(KAA_array[0]*Af_array[1]*Af_array[0]-AA_array[0])
# print(KCC_array[0]*Cf_array[1]*Cf_array[0]-CC_array[0])


# # Test A + {B} + {C} case without interactions
# Af_array = np.array([16])
# Bf_array = np.array([5,2.1])
# AB_array = np.array([20,2])
# Cf_array = np.array([1,2])
# BC_array = np.array([[3,8],[12,30]])
# AC_array = np.array([4,11])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, Cf_array, BC_array, \
# 	AC_array),axis=None)
# print(x_array)


# Ai_array = np.array([16])
# Bi_array = np.array([7,3.5])
# KAB_array = np.array([3,13])
# Ci_array = np.array([4,6])
# KBC_array = np.array([[10000,45],[61,85]])
# KAC_array = np.array([2.3,8.5])
# args = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(Ai_array[0]-Af_array[0]-AC_array[0]-AC_array[1]-AB_array[0]-AB_array[1])
# print(Bi_array[0]-Bf_array[0]-AB_array[0]-BC_array[0][0]-BC_array[0][1])
# print(Bi_array[1]-Bf_array[1]-AB_array[1]-BC_array[1][0]-BC_array[1][1])
# print(KAB_array[0]*Bf_array[0]*Af_array[0]-AB_array[0])
# print(KAB_array[1]*Bf_array[1]*Af_array[0]-AB_array[1])
# print(Ci_array[0]-Cf_array[0]-BC_array[0][0]-BC_array[1][0]-AC_array[0])
# print(Ci_array[1]-Cf_array[1]-BC_array[0][1]-BC_array[1][1]-AC_array[1])
# print(KBC_array[0][0]*Bf_array[0]*Cf_array[0]-BC_array[0][0])
# print(KBC_array[0][1]*Bf_array[0]*Cf_array[1]-BC_array[0][1])
# print(KBC_array[1][0]*Bf_array[1]*Cf_array[0]-BC_array[1][0])
# print(KBC_array[1][1]*Bf_array[1]*Cf_array[1]-BC_array[1][1])
# print(KAC_array[0]*Af_array[0]*Cf_array[0]-AC_array[0])
# print(KAC_array[1]*Af_array[0]*Cf_array[1]-AC_array[1])


# # Test A + {B} + {C} case without interactions
# Af_array = np.array([16])
# Bf_array = np.array([5,2.1])
# AB_array = np.array([20,2])
# Cf_array = np.array([1,2])
# BC_array = np.array([[3,8],[12,30]])
# AC_array = np.array([4,11])
# AA_array = np.array([])
# BB_array = np.array([3.2])
# CC_array = np.array([1.6])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, Cf_array, BC_array, \
# 	AC_array, AA_array, BB_array, CC_array),axis=None)
# print(x_array)


# Ai_array = np.array([16])
# Bi_array = np.array([7,3.5])
# KAB_array = np.array([3,13])
# Ci_array = np.array([4,6])
# KBC_array = np.array([[10000,45],[61,85]])
# KAC_array = np.array([2.3,8.5])
# KAA_array = np.array([])
# KBB_array = np.array([35])
# KCC_array = np.array([43])
# args = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array, \
# 	KAA_array, KBB_array, KCC_array)

# trial = Keq(x_array, *args)
# print(trial)

# print(Ai_array[0]-Af_array[0]-AC_array[0]-AC_array[1]-AB_array[0]-AB_array[1])
# print(Bi_array[0]-Bf_array[0]-AB_array[0]-BC_array[0][0]-BC_array[0][1]-BB_array[0])
# print(Bi_array[1]-Bf_array[1]-AB_array[1]-BC_array[1][0]-BC_array[1][1]-BB_array[0])
# print(KAB_array[0]*Bf_array[0]*Af_array[0]-AB_array[0])
# print(KAB_array[1]*Bf_array[1]*Af_array[0]-AB_array[1])
# print(Ci_array[0]-Cf_array[0]-BC_array[0][0]-BC_array[1][0]-AC_array[0]-CC_array[0])
# print(Ci_array[1]-Cf_array[1]-BC_array[0][1]-BC_array[1][1]-AC_array[1]-CC_array[0])
# print(KBC_array[0][0]*Bf_array[0]*Cf_array[0]-BC_array[0][0])
# print(KBC_array[0][1]*Bf_array[0]*Cf_array[1]-BC_array[0][1])
# print(KBC_array[1][0]*Bf_array[1]*Cf_array[0]-BC_array[1][0])
# print(KBC_array[1][1]*Bf_array[1]*Cf_array[1]-BC_array[1][1])
# print(KAC_array[0]*Af_array[0]*Cf_array[0]-AC_array[0])
# print(KAC_array[1]*Af_array[0]*Cf_array[1]-AC_array[1])
# print(KBB_array[0]*Bf_array[1]*Bf_array[0]-BB_array[0])
# print(KCC_array[0]*Cf_array[1]*Cf_array[0]-CC_array[0])


# Test A + B case for fsolve
# Af_array = np.array([1])
# Bf_array = np.array([1])
# AB_array = np.array([0])
# x_array = np.concatenate((Af_array, Bf_array, AB_array),axis=None)
# print('this is x_array')
# print(x_array)

# Ai_array = np.array([1])
# Bi_array = np.array([1])
# KAB_array = np.array([1])
# args = (Ai_array, Bi_array, KAB_array)


# # # # Test A + B + C case for fsolve
# Af_array = np.array([0])
# Bf_array = np.array([0])
# Cf_array = np.array([0])
# AB_array = np.array([0])
# BC_array = np.array([0])
# AC_array = np.array([0])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, Cf_array, BC_array, AC_array),axis=None)
# print('this is x_array')
# print(x_array)

# Ai_array = np.array([1])
# Bi_array = np.array([1])
# Ci_array = np.array([1])
# KAB_array = np.array([10000])
# KBC_array = np.array([10000])
# KAC_array = np.array([10000])
# args = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array)


# # Test {A} + {B} + {C} case without interactions for fsolve
# Bf_array = np.array([0,0])
# Af_array = np.array([0,0])
# AB_array = np.array([[0,0],[0,0]])
# Cf_array = np.array([0, 0])
# BC_array = np.array([[0,0],[0,0]])
# AC_array = np.array([[0,0],[0,0]])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, Cf_array, BC_array, AC_array),axis=None)
# print(x_array)

# Bi_array = np.array([1,1])
# Ai_array = np.array([1,1])
# KAB_array = np.array([[1000,1000],[1000,1000]])
# Ci_array = np.array([1,1])
# KBC_array = np.array([[1000,1000],[1000,1000]])
# KAC_array = np.array([[1000,1000],[1000,1000]])
# args = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array)


# Test {A} + {B} + {C} case with interactions fsolve
# Af_array = np.array([0,0])
# Bf_array = np.array([0,0])
# AB_array = np.array([[0,0],[0,0]])
# Cf_array = np.array([0,0])
# BC_array = np.array([[0,0],[0,0]])
# AC_array = np.array([[0,0],[0,0]])
# AA_array = np.array([0])
# BB_array = np.array([0])
# CC_array = np.array([0])
# x_array = np.concatenate((Af_array, Bf_array, AB_array, Cf_array, BC_array, \
# 	AC_array, AA_array, BB_array, CC_array),axis=None)
# print(x_array)

# Bi_array = np.array([1,1])
# Ai_array = np.array([1,1])
# KAB_array = np.array([[1000,1000],[1000,1000]])
# Ci_array = np.array([1, 1])
# KBC_array = np.array([[1000,1000],[1000,1000]])
# KAC_array = np.array([[1000,1000],[1000,1000]])
# KAA_array = np.array([1000])
# KBB_array = np.array([1000])
# KCC_array = np.array([1000])
# args = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array, \
# 	KAA_array, KBB_array, KCC_array)

# trial = Keq(x_array, *args)
# print(trial)

# start = timer()
# z = fsolve(Keq,x_array, args)
# end = timer()
# t = end-start
# print('time = %f' % t)
# print(z)

# Am_array = calc_A_measured(z, *args)
# print(Am_array)


# Test a lot of stuff
Af_array = np.array([1,1])
Bf_array = np.array([1,1])
AB_array = np.array([[0,0],[0,0]])
Cf_array = np.array([1,1])
BC_array = np.array([[0,0],[0,0]])
AC_array = np.array([[0,0],[0,0]])
AA_array = np.array([0])
BB_array = np.array([0])
CC_array = np.array([0])

x_array = np.concatenate((Af_array, Bf_array, AB_array, AA_array, BB_array),axis=None)
xc_array = np.concatenate((Af_array, Bf_array, AB_array, Cf_array, BC_array, \
	AC_array, AA_array, BB_array, CC_array),axis=None)
# print(x_array)

Bi_array = np.array([1,1])
Ai_array = np.array([1,1])
# KAB_array = np.array([[1000,1000],[1000,1000]])
Ci_array = np.array([1, 1])
KBC_array = np.array([[1,1],[1,1]])
KAC_array = np.array([[1,1],[1,1]])
KAA_array = np.array([1])
KBB_array = np.array([1])
KCC_array = np.array([1])

KAB_vals = np.logspace(-5, 5, 20)

Am_no_c_array = np.empty((20,2))
Am_with_c_array = np.empty((20,2))
Am_diff_array = np.empty((20,2))

start = timer()
for i in range(len(KAB_vals)):
	KAB = KAB_vals[i]
	KAB_array = np.array([[KAB,KAB],[KAB,KAB]])
	args = (Ai_array, Bi_array, KAB_array, KAA_array, KBB_array)
	args_c = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array, \
	KAA_array, KBB_array, KCC_array)

	# Calculate the equilibrium result
	no_c = fsolve(Keq, x_array, args)
	with_c = fsolve(Keq, xc_array, args_c)
	print('this is KAB %f' % KAB)
	print('this is no c')
	print(no_c)
	print('this is with c')
	print(with_c)

	Am_no_c = calc_A_measured(no_c, *args)
	Am_with_c = calc_A_measured(with_c, *args_c)

	Am_diff = Am_no_c - Am_with_c

	Am_no_c_array[i] = Am_no_c
	Am_with_c_array[i] = Am_with_c
	Am_diff_array[i] = Am_diff

end = timer()
t = end-start
print('time = %f' % t)

print(Am_no_c_array)
print(Am_with_c_array)
print(Am_diff_array)

plt.semilogx(KAB_vals, Am_no_c_array, 'r')
plt.semilogx(KAB_vals, Am_with_c_array, 'g')
plt.show()


# 
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.optimize import root
from scipy.special import comb
from timeit import default_timer as timer

def AA_interactions(AA_array, F, n, N):
    '''This function takes in an array of in-set double stranded concentrations, 
    such as AA, BB, or CC. It also takes in an value F = Ai - Af - {AB} - {AC}, 
    and an index n. This function computes the in-set interactions that 
    contribute to the final concentration of strand An. In other words for the 
    nth strand of A, this function computes the contribution of {AA} to the 
    relation 0 = Ai - Af - {AB} - {AC} - {AA} and returns 
    F = Ai - Af - {AB} - {AC} - {AA} as a float.'''

    # Handle case with only 1 strand and 1 in-set interaction
    if len(AA_array) == 1:
        F = F - AA_array[0]

    # Handle case with > 2 strands and > 1 in-set interactions
    elif len(AA_array) > 1 and n == 0:
        F = F - np.sum(AA_array[:N])
    elif len(AA_array) > 1 and n > 0:
        start = 0
        for i in range(n):
            start += N-i
        AnA_row_interactions = np.sum(AA_array[start:start+N-n])

        start = 0
        AnA_col_interactions = 0
        for i in range(n):
            AnA_col_interactions+=AA_array[start+n]
            start += N-i-1
        F = F - AnA_row_interactions - AnA_col_interactions
    # No in-set interactions
    else:
        pass
    return F

def AA_interactions_no_self(AA_array, F, n, N):
    '''This function takes in an array of in-set double stranded concentrations, 
    such as AA, BB, or CC. It also takes in an value F = Ai - Af - {AB} - {AC}, 
    and an index n. This function computes the in-set interactions that 
    contribute to the final concentration of strand An. In other words for the 
    nth strand of A, this function computes the contribution of {AA} to the 
    relation 0 = Ai - Af - {AB} - {AC} - {AA} and returns 
    F = Ai - Af - {AB} - {AC} - {AA} as a float.'''

    # Handle case with only 2 strands and 1 in-set interaction
    if len(AA_array) == 1:
        F = F - AA_array[0]

    # Handle case with > 2 strands and > 1 in-set interactions
    elif len(AA_array) > 1 and n == 0:
        F = F - np.sum(AA_array[:N-1])
    elif len(AA_array) > 1 and n > 0:
        start = 0
        for i in range(n):
            start += N-i-1
        print(AA_array[start:start+N-n-1])
        AnA_row_interactions = np.sum(AA_array[start:start+N-n-1])

        start = 0
        AnA_col_interactions = 0
        for i in range(n):
            print(AA_array[start+n-1])
            AnA_col_interactions+=AA_array[start+n-1]
            start += N-i-2
        F = F - AnA_row_interactions - AnA_col_interactions
    # No in-set interactions
    else:
        pass
    return F

def calc_Af_equation(Ai_array, Af_array, N, M, L, AB_array, AC_array, AA_array):
    '''This function takes in:
    Ai_array (arr): a 1xN array of initial concentrations for {A}
    Af_array (arr): a 1xN array of final concentrations for {A}
    N (int): the number of strands in {A}
    M (int): the number of strands in {B}
    L (int): the number of strands in {C}
    AB_array (arr): a NxM array of AB concentrations
    AC_array (arr): a NxL array of AC concentrations
    AA_array (arr): an array of AA concentrations. It is 1xcomb(N+1,2) if 
        self-interactions are permitted, and 1xcomb(N,2) if not.

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
        F = AA_interactions(AA_array, F, n, N)
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
        F = AA_interactions(BB_array, F, m, M)

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
        F = AA_interactions(CC_array, F, l, L)

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

def calc_AA_equations_no_self(K_array, AA_array, Af_array):
    '''This function takes in:
    K_array (arr): a 1xcomb(N,2) array of KAA values
    AA_array (arr): a 1xcomb(N,2) array of AA values
    Af_array (arr): a 1xN array of Af concentrations

    This is compatible with a network allowing no self-interactions, and just
    reduces the number of equations.

    It calculates and returns F_array, a 1xcomb(N,2) array. The interactions are ordered
    as [F_0.1, F_0.2, ..., F_0.N-1, F_1.2, F_1.3, ..., F_1.N-1, ..., F_N-2.N-1]
    
    F_n.m = KAnAm * Afn * Afm - AnAm.'''

    N = len(Af_array)

    # Create an array to store F equations in
    F_array = np.zeros(len(AA_array))

    # Handle 1 AA case:
    if N == 1:
        F_array[0] = K_array[0]*Af_array[0]*Af_array[1] - AA_array[0]
    else: 
        i = 0
        # Iterate through the nth strands
        for n in range(N-1):
            M = N-n
            for m in range(M):
                F_array[i] = K_array[i] * Af_array[n] * Af_array[n+m+1] - AA_array[i]
                i+=1
    return F_array

def calc_AA_equations(K_array, AA_array, Af_array):
    '''This function takes in:
    K_array (arr): a 1xcomb(N+1,2) array of KAA values
    AA_array (arr): a 1xcomb(N+1,2) array of AA values
    Af_array (arr): a 1xN array of Af concentrations

    It calculates and returns F_array, a 1xcomb(N+1,2) array. The interactions are ordered
    as [F_0.0, F_0.1, ..., F_0.N-1, F_1.1, F_1.2, ..., F_1.N-1, ..., F_N-1.N-1]
    
    F_n.m = KAnAm * Afn * Afm - AnAm.'''

    N = len(Af_array)

    # Create an array to store F equations in
    F_array = np.zeros(len(AA_array))

    # Handle 1 AA case:
    if N == 1:
        F_array[0] = K_array[0]*Af_array[0]**2 - AA_array[0]
    else: 
        i = 0
        # Iterate through the nth strands
        for n in range(N):
            M = N-n
            for m in range(M):
                F_array[i] = K_array[i] * Af_array[n] * Af_array[n+m] - AA_array[i]
                i+=1
    return F_array


# def reshape_AA_array(N, AA_flat_array, KAA_array, AA_array):
#     '''This function takes in:
#     N = # of strands
#     AA_flat_array = a 1xcomb(N,2) array of AA concentrations
#     KAA_array = a 1x(N-1) array of KAA values whose rows go from 1x(N-1) to 1x1
#     AA_array = an array with the same dimensions as KAA_array to store AA 
#                 concentrations in.

#     It returns AA_array with the values of AA_flat_array stored inside.'''

#     row_start = 0
#     row_stop = 0
#     for i in range(len(KAA_array)):
#         # Figure out how long the row is
#         row_length = len(KAA_array[i])
#         row_stop = row_stop + row_length
#         # Find the row
#         AA_row = AA_flat_array[row_start:row_stop]
#         # Insert the row
#         AA_array[i] = AA_row

#         row_start = row_stop

#     return(AA_array)

# def fill_AA_array(x_array, start, stop, N, KAA_array, AA_array):
#     '''This function takes in:
#     x_array = a 1 x (N + M + NxM + L + MxL + NxL + comb(N,2) + comb(M,2) 
#     + comb(L,2)) array consisting of the following: 
#     [Af1, ...AfN - first N values, = A final concentrations
#      Bf1, ...BfM - next M values, = B final concentrations
#      AB_11, ...AB_NM - next NxM values, = AB concentrations
#      Cf1, ...CfL - next L values, = C final concentrations
#      BC_11, ...BC_ML - next MxL values, = BC final concentrations
#      AC_11, ...AC_NL - next NxL values, = AC final concentrations
#      AA_11, ...AA_NN - next comb(N,2) values, = AA final concentrations 
#      BB_11, ...BB_MM - next comb(M,2) values, = BB final concentrations 
#      CC_11, ...CC_LL - last comb(L,2) values, = CC final concentrations.]
#     stop = an indice used to slice x_array
#     N = # of strands
#     AA_flat_array = a 1xcomb(N,2) array of AA concentrations
#     KAA_array = a 1x(N-1) array of KAA values whose rows go from 1x(N-1) to 1x1
#     AA_array = an array with the same dimensions as KAA_array to store AA 
#                 concentrations in.

#     It returns a tuple of (start, stop, AA_array) in which start and stop have 
#     been updated and AA_array has the values of AA_flat_array stored inside.'''
    
#     # Handle the case where there is only one KAA value
#     if len(KAA_array) == 1:
#         start = stop 
#         stop = int(stop + 1)
#         AA_array[0] = x_array[start:stop]

#     # Handle the case where there is more than one KAA value
#     elif len(KAA_array) > 1:
#         start = stop 
#         stop = int(stop + comb(N,2))
#         AA_flat_array = x_array[start:stop]
#         AA_array = reshape_AA_array(N, AA_flat_array, KAA_array, AA_array)

#     else: 
#         pass

#     return(start, stop, AA_array)

def flatten_AA_list(AA_list):
    '''This function flattens an in-set interaction array ({AA}, {BB}, {CC}) so
    that it can be concatenated with other arrays into x_array.'''
    flattened = np.array([])
    for i in range(len(AA_list)):
        row = np.asarray(AA_list[i])
        flattened = np.concatenate((flattened, row),axis=None)
    return(flattened)

def gen_x_array(*array_args):
    '''This function concatenates arrays into args and returns x_array, a 
    flattened version'''
    # Unpack the A and B arrays
    Af_array = array_args[0]
    Bf_array = array_args[1]
    AB_array = array_args[2]
    x_array = np.concatenate((Af_array, Bf_array, AB_array),axis=None)
    
    # Unpack the AA and BB in-set interactions for a system with no {C}
    if len(array_args) == 5:
        AA_array = array_args[3]
        BB_array = array_args[4]
        flat_AA_array = flatten_AA_list(AA_array)
        flat_BB_array = flatten_AA_list(BB_array)
        x_array = np.concatenate((x_array, flat_AA_array, flat_BB_array),axis=None)
    else:
        pass
    
    # Unpack the C arrays if they are included
    if len(array_args) >= 6:
        Cf_array = array_args[3]
        BC_array = array_args[4]
        AC_array = array_args[5]
        x_array = np.concatenate((x_array, Cf_array, BC_array, AC_array),axis=None)
    else:
        pass
    
    # Unpack the self-interaction arrays in a system with {C}
    if len(array_args) == 9:
        AA_array = array_args[6]
        BB_array = array_args[7]
        CC_array = array_args[8]
        # Flatten the self-interaction arrays
        flat_AA_array = flatten_AA_list(AA_array)
        flat_BB_array = flatten_AA_list(BB_array)
        flat_CC_array = flatten_AA_list(CC_array)
        x_array = np.concatenate((x_array, flat_AA_array, flat_BB_array, flat_CC_array),axis=None)
    else:
        pass

    return x_array

def unpack_x_array(x_array, *args):
    '''This function takes in x_array and uses args to unpack it into its 
    constituent arrays. It also unpacks args and determines the dimensions of 
    the system. If certain arrays are not necessary, they are rep

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

    This function takes a tuple of "optional arguments".
    args is a tuple consisting of:
    ([Ai1, ...AiN] (arr): a 1xN array of A initial concentrations
     [Bi1, ...BiM] (arr): a 1xM array of B initial concentrations
     [KAB_11, ...KAB_NM] (NxM), = KAB values
     [Ci1, ...CiL] (arr): a 1xN array of  initial concentrations
     [KBC_11, ...KBC_ML] (MxL), = KBC values
     [KAC_11, ...KAC_NL] (NxL), = KAC values
     [KAA_11, ...KAA_NN] (array of whose rows are arrays 1x(N-1) to 1x1), 
                         = KAA values
     [KBB_11, ...KBB_MM] (1x(M-1) array of whose rows are arrays 1x(M-1) to 1x1), 
                         = KBB values
     [KCC_11, ...KCC_LL] (1x(L-1) array of whose rows are arrays 1x(L-1) to 1x1) 
                         = KCC values))

    Output: It returns (x_tuple, args_tuple, size_tuple) where
    x_tuple is a tuple of (Af_array, Bf_array, AB_array, Cf_array, BC_array, 
    AC_array, AA_array, BB_array, CC_array) with empty arrays taking the place
    of arrays that aren't provided in x_array.

    args_tuple is a tuple of (Ai_array, Bi_array, KAB_array, Ci_array, 
    KBC_array, KAC_array, KAA_array, KBB_array, KCC_array) with empty arrays 
    taking the place of arrays that aren't provided in args.

    size_tuple is a tuple of (N, M, L), which are the number of strands in {A},
    {B}, and {C}.'''

    # Unpack the arguments
    Ai_array = args[0].astype(np.float64)
    Bi_array = args[1].astype(np.float64)
    KAB_array = args[2].astype(np.float64)

    # Figure out dimensions of system
    N = len(Ai_array)
    M = len(Bi_array)
    # Set L = 0 for now
    L = 0

    # Unpack the guess
    # start and stop are slicing indices for x-array
    # Unpack Af_array
    start = 0
    stop = N
    Af_array = x_array[start:stop].astype(np.float64)
    # print('I am Af_array')
    # print(Af_array)

    # Unpack Bf_array
    start = stop
    stop = stop + M
    Bf_array = x_array[start:stop].astype(np.float64)
    # print('I Am Bf_array')
    # print(Bf_array)

    # Unpack AB_array
    start = stop 
    stop = int(stop + N*M)
    
    # Handle the case where there is only 1 A strand and multiple B strands
    if N == 1 and M > 1:
        AB_array = x_array[start:stop].astype(np.float64)

    # Handle all other cases
    else:
        AB_array = x_array[start:stop].reshape(N, M).astype(np.float64)

    # print('I am AB_array')
    # print(AB_array)

    # Make a bunch of empty arrays for x tuple
    Cf_array = np.array([]).astype(np.float64)
    BC_array = np.array([]).astype(np.float64)
    AC_array = np.array([]).astype(np.float64)
    AA_array = np.array([]).astype(np.float64)
    BB_array = np.array([]).astype(np.float64)
    CC_array = np.array([]).astype(np.float64)

    # Make a bunch of empty arrays for args
    Ci_array = np.array([]).astype(np.float64)
    KBC_array = np.array([]).astype(np.float64)
    KAC_array = np.array([]).astype(np.float64)
    KAA_array = np.array([]).astype(np.float64)
    KBB_array = np.array([]).astype(np.float64)
    KCC_array = np.array([]).astype(np.float64)

    # Handle the {A} + {B} case where there are {A} and/or {B} in-set interactions 
    if len(args) == 5:
        # Unpack the KAA and KBB arrays
        KAA_array = args[3]
        KBB_array = args[4]

        # Figure out the shape of KAA and KBB arrays
        AA_len = len(KAA_array)
        BB_len = len(KBB_array)

        # Make the AA array
        start = stop 
        stop = stop + AA_len
        AA_array = x_array[start:stop].astype(np.float64)

        # Make the AA array
        start = stop 
        stop = stop + BB_len
        BB_array = x_array[start:stop].astype(np.float64)

        # print('I am AA_array')
        # print(AA_array)
        # print('I am BB_array')
        # print(BB_array)

    # Handle the {A} + {B} + {C} case with no in-set interactions
    if len(args) > 5:
        # Unpack Ci, KBC, and KAC
        Ci_array = args[3].astype(np.float64)
        KBC_array = args[4].astype(np.float64) 
        KAC_array = args[5].astype(np.float64)

        # Figure out the number of strands in {C}
        L = len(Ci_array)

        # Unpack Cf_array
        start = stop 
        stop = int(stop + L)
        Cf_array = x_array[start:stop].astype(np.float64)
        
        # Unpack BC_array
        start = stop 
        stop = int(stop + M*L)

        # Handle the case where there is only 1 B strand and multiple C strands
        if M == 1 and L > 1:
            BC_array = x_array[start:stop].astype(np.float64)

        # Handle all other cases
        else:
            BC_array = x_array[start:stop].reshape(M,L).astype(np.float64)
        
        # Unpack AC_array
        start = stop 
        stop = int(stop + N*L)

        # Handle the case where there is only 1 A strand and multiple C strands
        if N == 1 and L > 1:
            AC_array = x_array[start:stop].astype(np.float64)

        # Handle all other cases
        else:
            AC_array = x_array[start:stop].reshape(N,L).astype(np.float64)
        
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

    # Handle the {A} + {B} + {C} case with self-interactions
    if len(args) == 9:
        # Unpack the K_arrays
        KAA_array = args[6]
        KBB_array = args[7]
        KCC_array = args[8]

        # Figure out the shape of KAA and KBB arrays
        AA_len = len(KAA_array)
        BB_len = len(KBB_array)
        CC_len = len(KCC_array)

        # Make the AA array
        start = stop 
        stop = stop + AA_len
        AA_array = x_array[start:stop].astype(np.float64)

        # Make the BB array
        start = stop 
        stop = stop + BB_len
        BB_array = x_array[start:stop].astype(np.float64)

        # Make the CC array
        start = stop 
        stop = stop + CC_len
        CC_array = x_array[start:stop].astype(np.float64)

        # print('I am AA_array')
        # print(AA_array)
        # print('I am BB_array')
        # print(BB_array)
        # print('I am CC_array')
        # print(CC_array)

    # Make the x tuple
    x_tuple = (Af_array, Bf_array, AB_array, Cf_array, BC_array, AC_array, \
            AA_array, BB_array, CC_array)

    # Make the args tuple
    args_tuple = (Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, \
        KAC_array, KAA_array, KBB_array, KCC_array)

    # Make the size tuple
    size_tuple = (N,M,L)

    return (x_tuple, args_tuple, size_tuple)

def Keq_v0(x_array, *args):
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
    equations = 0 to be    used with scipy.optimize.fsolve.'''

    # args = arguments[0]
    # print('These are the arguments')
    # print(arguments)

    # print('These are args')
    # print(args)

    # Unpack x_array, args, and dimensions of the system
    x_tuple, args_tuple, size_tuple = unpack_x_array(x_array, *args)
    # Unpack x_tuple
    Af_array, Bf_array, AB_array, Cf_array, BC_array, AC_array, AA_array, \
    BB_array, CC_array = x_tuple
    # Unpack the args
    Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array, KAA_array, \
    KBB_array, KCC_array = args_tuple
    # Unpack the dimensions of the system
    N, M, L = size_tuple

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
            # Flatten F_AA_array so it isn't an array of arrays
            flat_F_AA_array = flatten_AA_array(F_AA_array)
            # Enter the equations into F_array
            f = f+1
            F_array[f] = flat_F_AA_array

        else:
            pass

        # Calculate the BB equations if there are any interactions
        if len(BB_array) != 0:
            F_BB_array = calc_AA_equations(KBB_array, BB_array, Bf_array)
            # Flatten F_BB_array so it isn't an array of arrays
            flat_F_BB_array = flatten_AA_array(F_BB_array)
            # Enter the equations into F_array
            f = f+1
            F_array[f] = flat_F_BB_array

        else:
            pass
    else:
        pass

    if len(args) == 9:
        # Calculate the CC equations if there are any interactions
        if len(CC_array) != 0:
            F_CC_array = calc_AA_equations(KCC_array, CC_array, Cf_array)
            # Flatten F_CC_array so it isn't an array of arrays
            flat_F_CC_array = flatten_AA_array(F_CC_array)
            # Ebter the equations into F_array
            f = f+1
            F_array[f] = flat_F_CC_array
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


def Keq(x_array, *args):
    '''This function creates a system of equations to be used with 
    scipy.optimize.minimize. This represents a DNA network with 3 sets of strands
    (A, B, C), which may have pairwise interactions. For example, the ith A 
    strand and the jth B strand may interact as Ai + Bj <-> ABij, with 
    equilibrium constant KABij. At equilibrium this interaction satisfies
    KABij[Af][Bf] - [AB] = 0. All strands must also obey this equality:
    Ai - Af - {AB} - {AC} - {AA} = 0. They must also obey this inequality:
    0 <= Af <= Ai. These inequalities are accounted for by using bounds with 
    scipy.optimize.minimize.

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

    This function takes a tuple of "optional arguments".
    args is a tuple consisting of:
    ([Ai1, ...AiN] (1xN), = A initial concentrations
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
                         = KCC values)

    Output: It returns F_array, (same dimensions as x_array) of a system of 
    equations = 0 to be    used with scipy.optimize.minimize. Since 
    scipy.optimize.minimize looks for the minimum of the function, and our 
    system of equations = [0], we minimize the equations 
    (KABij[Af][Bf] - [AB])^2 = 0 and (Ai - Af - {AB} - {AC} - {AA})^2 = 0'''

    # args = arguments[0]
    # print('These are the arguments')
    # print(arguments)

    # print('These are args')
    # print(args)

    # Unpack x_array, args, and dimensions of the system
    x_tuple, args_tuple, size_tuple = unpack_x_array(x_array, *args)
    # Unpack x_tuple
    Af_array, Bf_array, AB_array, Cf_array, BC_array, AC_array, AA_array, \
    BB_array, CC_array = x_tuple
    # Unpack the args
    Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array, KAA_array, \
    KBB_array, KCC_array = args_tuple
    # Unpack the dimensions of the system
    N, M, L = size_tuple
    
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
            F_array[f] = flat_F_AA_array

        else:
            pass

        # Calculate the BB equations if there are any interactions
        if len(BB_array) != 0:
            F_BB_array = calc_AA_equations(KBB_array, BB_array, Bf_array)
            # Enter the equations into F_array
            f = f+1
            F_array[f] = flat_F_BB_array

        else:
            pass
    else:
        pass

    if len(args) == 9:
        # Calculate the CC equations if there are any interactions
        if len(CC_array) != 0:
            F_CC_array = calc_AA_equations(KCC_array, CC_array, Cf_array)
            # Enter the equations into F_array
            f = f+1
            F_array[f] = flat_F_CC_array
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

    # print('This is F_flat_array')
    # print(F_flat_array)

    # Calculate F dot F
    F_out = np.dot(F_flat_array, F_flat_array)

    return(F_out)

def gen_bounds_A(Ai_array, bound_list):
    '''This function takes in an array of initial strand concentrations and a 
    list of bounds to be used with scipy.optimize.minimize. For each strand it
    generates a tuple of bounds (0, Ai) and adds it to the list of bounds. This
    function returns the updated list of bounds.'''
    # print('gen A is being called for')
    # print(Ai_array)

    for n in range(len(Ai_array)):
        An = Ai_array[n]
        bound = (0, An)
        bound_list.append(bound)

    # print('this is the bound_list now')
    # print(bound_list)
    return(bound_list)

def gen_bounds_AB(Ai_array, Bi_array, bound_list):
    '''This function takes in two arrays of initial strand concentrations and a 
    list of bounds to be used with scipy.optimize.minimize. For each strand pair
    AnBm it generates a tuple of bounds (0, min(An,Bm)) and adds it to the list 
    of bounds. This function returns the updated list of bounds.'''

    # print('gen AB is being called for')
    # print(Ai_array)
    # print(Bi_array)

    # Iterate through the A strands
    for n in range(len(Ai_array)):
        An = Ai_array[n]

        # Iterate through the B strands
        for m in range(len(Bi_array)):
            Bm = Bi_array[m]

            # Pick the minimum of An and Bm as the upper bound for AnBm
            bound = (0, min(An,Bm))
            bound_list.append(bound)

    # print('this is the bound_list now')
    # print(bound_list)
    return(bound_list)

def gen_bounds_AA_no_self(Ai_array, bound_list):
    '''This function takes in an array of initial strand concentrations and a 
    list of bounds to be used with scipy.optimize.minimize. For each strand pair
    AnAm it generates a tuple of bounds (0, min(An,Am)) and adds it to the list 
    of bounds. This function returns the updated list of bounds.'''
    # print('gen AA is being called for')
    # print(Ai_array)

    N = len(Ai_array)
    # Iterate through the A strands
    for n in range(N-1):
        An = Ai_array[n]

        # Iterate through the other A strands
        for m in range(N-1-n):
            Am = Ai_array[n + m + 1]

            # Pick the minimum of An and Bm as the upper bound for AnBm
            bound = (0, min(An,Am))
            bound_list.append(bound)

    # print('this is the bound_list now')
    # print(bound_list)
    return(bound_list)

def gen_bounds_AA(Ai_array, bound_list):
    '''This function takes in an array of initial strand concentrations and a 
    list of bounds to be used with scipy.optimize.minimize. For each strand pair
    AnAm it generates a tuple of bounds (0, min(An,Am)) and adds it to the list 
    of bounds. This function returns the updated list of bounds.'''
    # print('gen AA is being called for')
    # print(Ai_array)

    N = len(Ai_array)
    # Iterate through the A strands
    for n in range(N):
        An = Ai_array[n]

        # Iterate through the other A strands
        for m in range(N-n):
            Am = Ai_array[n + m]

            # Pick the minimum of An and Bm as the upper bound for AnBm
            bound = (0, min(An,Am))
            bound_list.append(bound)

    # print('this is the bound_list now')
    # print(bound_list)
    return(bound_list)

def gen_bounds(*args):
    '''This function creates bounds for a system of equations to be used with 
    scipy.optimize.minimize. This system is a DNA network with 3 sets of strands
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

    This function takes a tuple of "optional arguments".
    args is a tuple consisting of:
    [[Ai1, ...AiN] (1xN), = A initial concentrations
     [Bi1, ...BiM] (1xM), = B initial concentrations
     [KAB_11, ...KAB_NM] (NxM), = KAB values
     [Ci1, ...CiL] (1xL), = C initial concentrations (optional)
     [KBC_11, ...KBC_ML] (MxL), = KBC values (optional)
     [KAC_11, ...KAC_NL] (NxL), = KAC values (optional)
     [KAA_11, ...KAA_NN] (1x(N-1) array of whose rows are arrays 1x(N-1) to 1x1), 
                         = KAA values (optional)
     [KBB_11, ...KBB_MM] (1x(M-1) array of whose rows are arrays 1x(M-1) to 1x1), 
                         = KBB values (optional)
     [KCC_11, ...KCC_LL] (1x(L-1) array of whose rows are arrays 1x(L-1) to 1x1) 
                         = KCC values (optional)]

    Output: It returns a tuple of bounds for the elements of x_array to be used 
    with scipy.optimize.minimize.'''

    # Create a list to store bounds in 
    bound_list = []

    # Unpack the arguments
    Ai_array = args[0].astype(np.float64)
    Bi_array = args[1].astype(np.float64)

    # Figure out the bounds for Af, Bf, and AB
    bound_list = gen_bounds_A(Ai_array, bound_list)
    bound_list = gen_bounds_A(Bi_array, bound_list)
    bound_list = gen_bounds_AB(Ai_array, Bi_array, bound_list)

    # Handle the {A} + {B} case where there are {A} and/or {B} in-set interactions 
    if len(args) == 5:
        # Figure out the bounds for AA and BB
        bound_list = gen_bounds_AA(Ai_array, bound_list)
        bound_list = gen_bounds_AA(Bi_array, bound_list)
    else:
        pass

    # Handle the {A} + {B} + {C} case with no in-set interactions
    if len(args) > 5:
        # Unpack Ci
        Ci_array = args[3].astype(np.float64)

        # Figure out the bounds for Cf, BC, and AC.
        bound_list = gen_bounds_A(Ci_array, bound_list)
        bound_list = gen_bounds_AB(Bi_array, Ci_array, bound_list)
        bound_list = gen_bounds_AB(Ai_array, Ci_array, bound_list)

    else:
        pass

    # Handle the {A} + {B} + {C} case with self-interactions
    if len(args) == 9:
        # Figure out the bounds for AA, BB, and CC
        bound_list = gen_bounds_AA(Ai_array, bound_list)
        bound_list = gen_bounds_AA(Bi_array, bound_list)
        bound_list = gen_bounds_AA(Ci_array, bound_list)

    else:
        pass

    bound_tuple = tuple(bound_list)
    return(bound_tuple)

# def create_empty_AA_array(N):
#     '''This function creates an array filled with zeros that is shaped like 
#     AA_array for N strands.'''
#     empty_AA_array = np.empty(N-1, dtype='object')
#     for i in range(N-1):
#         row = np.zeros(N-1-i)
#         empty_AA_array[i] = row
#     return empty_AA_array

# def calc_AA_jacobian_matrix(n, N):
#     '''This function takes in an strand index n and a number of strands N and 
#     returns AA_array filled with zeros except wherever there is an interaction 
#     with the nth strand. In these spots it returns -1.'''
#     AA_array = create_empty_AA_array(N)
#     for i in range(len(AA_array)):
#         row = AA_array[i]
#         for j in range(len(row)):
#             y = i + j + 1
#             if i == n or y == n:
#                 AA_array[i][j] = -1
#     return AA_array

# def calc_jacobian(x_array, *args):
#     '''This function takes in x array and "optional" args that define a DNA
#     network system. It calculates the jacobian for F_array, the output of Keq. 
#     The jacobian is a numpy matrix to be used with scipy.optimize.'''
#     # Unpack x_array, args, and dimensions of the system
#     x_tuple, args_tuple, size_tuple = unpack_x_array(x_array, *args)
#     # Unpack x_tuple
#     Af_array, Bf_array, AB_array, Cf_array, BC_array, AC_array, AA_array, \
#     BB_array, CC_array = x_tuple
#     # Unpack the args
#     Ai_array, Bi_array, KAB_array, Ci_array, KBC_array, KAC_array, KAA_array, \
#     KBB_array, KCC_array = args_tuple
#     # Unpack the dimensions of the system
#     N, M, L = size_tuple

#     # Make an array to store the Jacobian in
#     jac = np.zeros((len(x_array),len(x_array)),dtype='int')

#     # Calculate partial derivatives for {Af} equations
#     for n in range(N):
#         # derivative for Anf
#         jac[n][n] = -1
#         # derivatives for AnB
#         for m in range(M):
#             jac[n][N + M + n*M + m] = -1
#         # derivatives for AnC
#         for l in range(L):
#             jac[n][N + M + N*M + L + M*L + n*L + l] = -1
#         # derivatives for AnA if they're included
#         if len(KAA_array) > 0:
#             AA_jacobian = calc_AA_jacobian_matrix(n,N)
#             flat_AA_jacobian = flatten_AA_array(AA_jacobian)
#             for i in range(len(flat_AA_jacobian)):
#                 jac[n][N + M + N*M + L + M*L + N*L + i] = flat_AA_jacobian[i]

#     # Calculate partial derivatives for {Bf} equations
#     for m in range(M):
#         # derivative for Bnf
#         jac[N + m][N + m]= -1
#         # derivatives for ABm
#         for n in range(N):
#             jac[N + m][N + M + n*M + m] = -1
#         # derivatives for BmC
#         for l in range(L):
#             jac[N + m][N + M + N * M + L + m*L + l] = -1
#         # derivatives for BnB if they're included
#         if len(KBB_array) > 0:
#             BB_jacobian = calc_AA_jacobian_matrix(m,M)
#             flat_BB_jacobian = flatten_AA_array(BB_jacobian)
#             for i in range(len(flat_BB_jacobian)):
#                 jac[N + m][N + M + N*M + L + M*L + N*L + int(comb(N,2)) + i] = \
#                 flat_BB_jacobian[i]

#     # Calculate partial derivatives for {AB} equations
#     for n in range(N):
#         An = Af_array[n]
#         for m in range(M):
#             KAB = KAB_array[n][m]
#             Bm = Bf_array[m]
#             row = N + M + n*M + m
#             # Add Bm derivative
#             jac[row][n] = KAB*Bm
#             # Add An derivative
#             jac[row][N + m] = KAB*An
#             # Add AnBm derivative
#             jac[row][row] = -1

#     # Calculate partial derivatives if the system includes {C}
#     if L > 0:
#         # Calculate partial derivatives for {Cf} equations
#         for l in range(L):
#             row = N + M + N*M + l
#             # derivative for Clf
#             jac[row][row] = -1
#             # derivatives for BCl
#             for m in range(M):
#                 jac[row][N + M + N*M + L + m*L + l] = -1
#             # derivatives for ACl
#             for n in range(N):
#                 jac[row][N + M + N*M + L + M*L +n*L + l] = -1
#             # derivatives for ClC if they're included
#             if len(KCC_array) > 0:
#                 CC_jacobian = calc_AA_jacobian_matrix(l,L)
#                 flat_CC_jacobian = flatten_AA_array(CC_jacobian)
#                 for i in range(len(flat_CC_jacobian)):
#                     jac[row][N + M + N*M + L + M*L + N*L + int(comb(N,2)) + i \
#                     + int(comb(M,2))] = flat_CC_jacobian[i]

#         # Calculate partial derivatives for {BC} equations
#         for m in range(M):
#             Bm = Bf_array[m]
#             for l in range(L):
#                 KBC = KBC_array[m][l]
#                 Cl = Cf_array[l]
#                 row = N + M + N*M + L + m*L + l
#                 # Add Cl derivative
#                 jac[row][N + m] = KBC*Cl
#                 # Add Bm derivative
#                 jac[row][N + M + N*M + l] = KBC*Bm
#                 # Add AnBm derivative
#                 jac[row][row] = -1

#         # Calculate partial derivatives for {AC} equations
#         for n in range(N):
#             An = Af_array[n]
#             for l in range(L):
#                 KAC = KAC_array[n][l]
#                 Cl = Cf_array[l]
#                 row = N + M + N*M + L + M*L + n*L + l
#                 # Add Cl derivative
#                 jac[row][n] = KAC*Cl
#                 # Add Bm derivative
#                 jac[row][N + M + N*M + l] = KAC*An
#                 # Add AnBm derivative
#                 jac[row][row] = -1
#     else:
#         pass

#     # Calculate partial derivatives for {AA} equations
#     # Handle 1 AA interaction case
#     if len(KAA_array) == 1:
#         A0 = Af_array[0]
#         A1 = Af_array[1]
#         KAA = KAA_array[0]
#         row = N + M + N*M + L + M*L + N*L 
#         # A0 derivative
#         jac[row][0] = KAA * A1
#         # A1 derivative
#         jac[row][1] = KAA * A0
#         # A0A1 derivative
#         jac[row][row] = -1

#     # Handle multiple AA interaction case
#     elif len(KAA_array) > 1:
#         row = N + M + N*M + L + M*L + N*L 
#         for n in range(len(KAA_array)):
#             An = Af_array[n]
#             for i in range(len(KAA_array[n])):
#                 KAA = KAA_array[n][i]
#                 Ai = Af_array[n + i + 1]
#                 # An derivative
#                 jac[row][n] = KAA * Ai
#                 # Ai derivative
#                 jac[row][n + i + 1] = KAA * An
#                 # AnAi derivative
#                 jac[row][row] = -1
#                 row = row + 1
#     else:
#         pass

#     # Calculate partial derivatives for {BB} equations
#     # Handle 1 BB interaction case
#     if len(KBB_array) == 1:
#         B0 = Bf_array[0]
#         B1 = Bf_array[1]
#         KBB = KBB_array[0]
#         row = N + M + N*M + L + M*L + N*L + int(comb(N,2)) 
#         # B0 derivative
#         jac[row][N] = KBB * B1
#         # B1 derivative
#         jac[row][N + 1] = KBB * B0
#         # B0B1 derivative
#         jac[row][row] = -1

#     # Handle multiple BB interaction case
#     elif len(KBB_array) > 1:
#         row = N + M + N*M + L + M*L + N*L + int(comb(N,2))
#         for m in range(len(KBB_array)):
#             Bm = Bf_array[m]
#             for i in range(len(KBB_array[m])):
#                 KBB = KBB_array[m][i]
#                 Bi = Bf_array[m + i + 1]
#                 # Bm derivative
#                 jac[row][N + m] = KBB * Bi
#                 # Bi derivative
#                 jac[row][N + m + i + 1] = KBB * Bm
#                 # BmBi derivative
#                 jac[row][row] = -1
#                 row = row + 1
#     else:
#         pass

#     # Calculate partial derivatives for {CC} equations
#     # Handle 1 CC interaction case
#     if len(KCC_array) == 1:
#         C0 = Cf_array[0]
#         C1 = Cf_array[1]
#         KCC = KCC_array[0]
#         row = N + M + N*M + L + M*L + N*L + int(comb(N,2)) + int(comb(M,2))
#         # C0 derivative
#         jac[row][N + M + N*M] = KCC * C1
#         # C1 derivative
#         jac[row][N + M + N*M + 1] = KCC * C0
#         # C0C1 derivative
#         jac[row][row] = -1

#     # Handle multiple CC interaction case
#     elif len(KCC_array) > 1:
#         row = N + M + N*M + L + M*L + N*L + int(comb(N,2))+ int(comb(M,2))
#         for l in range(len(KCC_array)):
#             Cl = Cf_array[l]
#             for i in range(len(KCC_array[l])):
#                 KCC = KCC_array[l][i]
#                 Ci = Cf_array[l + i + 1]
#                 # Cl derivative
#                 jac[row][N + M + N*M + l] = KCC * Ci
#                 # Ci derivative
#                 jac[row][N + M + N*M + l + i + 1] = KCC * Cl
#                 # ClCi derivative
#                 jac[row][row] = -1
#                 row = row + 1
#     else:
#         pass

#     return jac


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
    args is a tuple of arrays:
        Ai_array (arr): a 1xN array [Ai1, ...AiN] of A initial concentrations
        Bi_array (arr): a 1xN array [Bi1, ...BiN] of B initial concentrations
        KAB_array (arr): a NxM array [AB_11, ...AB_NM] of KAB values
        Ci_array (arr): a 1xN array [Ci1, ...CiN] of C initial concentrations
        KBC_array (arr): a MxL array [BC_11, ...BC_ML] of KBC values
        KAC_array (arr): a NxL array [AC_11, ...AC_NL] of KAC values
        KAA_array (arr): a 1xcomb(N+1,2) array [KAA_11, ...KAA_NN] of KAA values
        KBB_array (arr): a 1xcomb(M+1,2) array [KBB_11, ...KBB_NN] of KBB values
        KCC_array (arr): a 1xcomb(L+1,2) array [KCC_11, ...KCC_NN] of KCC values
            Note: if self-interactions are not allowed, KAA_array is 1xcomb(N,2), 
            KBB_array is 1xcomb(M,2), and KCC_array is 1xcomb(L,2)

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


# def calc_Am_KAB_range(KAB_vals, x_array, args, x_array_c, args_c):
#     '''This function takes in numpy array of a range of KAB values as well as
#     arrays that define a DNA network system of {A}, {B}, and {C} strands. It 
#     computes the equilibrium concentrations of {A}, {B}, {C}, {AB}, {BC}, {AC},
#     {AA}, {BB}, and {CC}. It then computes the measured {A} concentrations and
#     returns 4 arrays: {A} measured with no {C} in the system, 
#     {A} measured with {C} in the system, 
#     The difference in {A} measured between the two systems, and
#     The percent difference in {A} measured between the two systems using 
#     {A} measured with no {C} in the system as the divisor.
#     '''
    
#     # Figure out the number of KAB values
#     num_KAB = len(KAB_vals)

#     # Figure out the number of {A} strands
#     N = len(args[0])

#     # Create arrays to store the {A} measured values for the system with {C} and
#     # without {C}, the difference in {A} measured between the two, and the 
#     # percent difference between the two.
#     Am_no_c_array = np.empty((num_KAB,N))
#     Am_with_c_array = np.empty((num_KAB,N))
#     Am_diff_array = np.empty((num_KAB,N))
#     Am_percent_diff_array = np.empty((num_KAB,N))

#     # Iterate through the KAB values and calculate A measured for the system 
#     # with {C} and without {C}.
#     for i in range(len(KAB_vals)):
#         # Figure out KAB
#         KAB = KAB_vals[i]

#         # Make KAB array
#         KAB_array = np.array([[KAB,KAB],[KAB,KAB]])

#         # Make the bounds to use with scipy.optimize.minimize
#         bounds_tuple = gen_bounds(*args)
#         bounds_tuple_c = gen_bounds(*args_c)

#         # Calculate the equilibrium result
#         no_c_result = minimize(Keq, x_array, args, method='SLSQP', \
#             bounds=bounds_tuple)
#         with_c_result = minimize(Keq, x_array_c, args_c, method='SLSQP', \
#             bounds=bounds_tuple_c)

#         # Pick out the result
#         no_c = no_c_result.x
#         with_c = with_c_result.x

#         # Calculate {A} measured
#         Am_no_c = calc_A_measured(no_c, *args)
#         Am_with_c = calc_A_measured(with_c, *args_c)

#         # Calculate the difference between the system with {C} and without {C}.
#         Am_diff = Am_no_c - Am_with_c

#         # Calculat the percentage change of {A} measured between the system with
#         # {C} and without {C}.
#         Am_percent_diff = np.divide(Am_diff, Am_no_c)

#         # Add the results to the arrays
#         Am_no_c_array[i] = Am_no_c
#         Am_with_c_array[i] = Am_with_c
#         Am_diff_array[i] = Am_diff
#         Am_percent_diff_array[i] = Am_percent_diff 

#     return(Am_no_c_array, Am_with_c_array, Am_diff_array, Am_percent_diff_array)
import numpy as np
from scipy.optimize import minimize

# from PertQuant.simCRN.gen_eq_data import compose_association_array

def calc_XY_equilibrium(Xf_array, Yf_array, XY_array, KXY_array):
    ''' This function calculates (KXY[X][Y] - [XY])^2 in order to minimize the
    equation (KXY[X][Y] - [XY])^2 = 0. 
    Parameters:
    Xf_array (arr): a 1xI array of X final concentrations.
    Yf_array (arr): a 1xJ array of Y final concentrations.
    XY_array (arr): a IxJ array of XY concentrations.
    KXY_array (arr): a IxJ array of KXY association constants.

    Returns:
    XY_F_array (arr): a 1xF (F <= I*J) array of equilibrium equations.

    '''
    # Get [X][Y]
    XY_outer = np.outer(Xf_array, Yf_array)
    # Calculated XY equilibrium expressions
    XY_eq = np.multiply(KXY_array, XY_outer) - XY_array
    # Get rid of the elements that have K = 0
    XY_nonzero_array = XY_eq[KXY_array>0]
    # Square the elements, since we are trying to minimize the equations
    XY_F_array = np.power(XY_nonzero_array,2)
    # Sum the elements
    XY_F = np.sum(XY_F_array)

    return XY_F

def calc_XX_equilibrium(Xf_array, XX_array, KXX_array):
    ''' This function calculates (KXX[X][X] - [XX])^2 in order to minimize the
    equation (KXX[X][X] - [XX])^2 = 0. 
    Parameters:
    Xf_array (arr): a 1xI array of X final concentrations.
    XX_array (arr): a IxI array of XX concentrations.
    KXX_array (arr): a IxI array of KXX association constants.

    Returns:
    XY_F_array (arr): a 1xF (F <= I*I) array of equilibrium equations.

    '''
    # Get [A][B]
    XX_outer = np.outer(Xf_array, Xf_array)
    # Calculated AB equilibrium expressions
    XX_eq = np.multiply(KXX_array, XX_outer) - XX_array
    # Get rid of the elements that have K = 0
    # Use the upper triangle of KXX_array to avoid duplicate reactions
    XX_nonzero_array = XX_eq[np.triu(KXX_array)>0]
    # Square the elements, since we are trying to minimize the equations
    XX_F_array = np.power(XX_nonzero_array,2)
    # Sum the elements
    XX_F = np.sum(XX_F_array)

    return XX_F

def calc_X_conservation(Xi_array, Xf_array, XY_array, XZ_array, XX_array, \
    X_connected):
    ''' This function calculates (Xi - Xf - {XY} - {XZ} - {X})^2  in order to 
    minimize the equation (Xi - Xf - {XY} - {XZ} - {X})^2 = 0.
    Parameters:
    Xi_array (arr): a 1xI array of X initial concentrations
    Xf_array (arr): a 1xI array of X final concentrations
    XY_array (arr): a IxJ array of XY concentrations
    XZ_array (arr): a IxK array of XZ concentrations
    XX_array (arr): a IxI array of XX concentrations
    X_connected (arr): a 1xI boolean array of whether Xi is connected to the 
        chemical reaction network

    Returns:
    X_conservation_F_array (arr): a 1xF (F <= I) array of conservation equations.
    '''
    # Add the {XiY}, {XiZ}, and {XiX} contributions
    XY_sum = np.sum(XY_array,axis=1)
    XZ_sum = np.sum(XZ_array,axis=1)
    XX_sum = np.sum(XX_array,axis=1)
    # Compute the conservation equations. 
    X_conservation = Xi_array - Xf_array - XY_sum - XZ_sum - XX_sum
    # Return the X concentrations that participate in the network
    X_conservation_connected_array = X_conservation[X_connected]
    # Square the elements since we are trying to minize the equations
    X_conservation_F_array = X_conservation_connected_array ** 2
    # Sum the elements
    X_conservation_F = np.sum(X_conservation_F_array)

    return X_conservation_F

def get_X_connected(KXY_array, KXZ_array, KXX_array):
    ''' This function takes in the K arrays relevant to {X} and figures out
    which elements of X participate in the chemical reactio network.
    Parameters:
    KXY_array (arr): a IxJ array of KXY association constants
    KXZ_array (arr): a IxK array of KXZ association constants
    KXX_array (arr): a IxI array of KXX association constants

    Returns:
    X_connected (arr): a 1xI boolean array of whether Xi is connected to the 
        chemical reaction network
    '''
    # Add the {XiY}, {XiZ}, and {XiX} contributions
    KXY_sum = np.sum(KXY_array,axis=1)
    KXZ_sum = np.sum(KXZ_array,axis=1)
    KXX_sum = np.sum(KXX_array,axis=1)
    KX_sum = KXY_sum + KXZ_sum + KXX_sum
    return KX_sum > 0

def ID_network_species_reactions(KAB_array, KBC_array, KAC_array, KAA_array, \
    KBB_array, KCC_array):
    '''This function takes in association constant arrays and identifies which
    species react in the network and which reactions are nonzero. For each
    species sets {A}, {B}, and {C} the reacting species are returned as a 
    boolean array. The indices of the reacting species are
    also returned as arrays. For the reactions, tuples of the arrays of the 
    nonzero values in the K_arrays are returned.
    
    The total number of reacting species in {A}, {B}, {C}, {AB}, {BC}, {AC}, 
    {AA}, {BB}, and {CC} is also returned as N_reactants.

    Parameters:
    KAB_array (arr): a NxM array of KAB association constants
    KBC_array (arr): a MxL array of KBC association constants
    KAC_array (arr): a NxL array of KAC association constants
    KAA_array (arr): a NxN array of KAA association constants
    KBB_array (arr): a MxM array of KBB association constants
    KCC_array (arr): a LxL array of KCC association constants

    Returns:
    ABC_connected (list): a list of [A_connected, B_connected, C_connected] where
        A_connected (arr): a 1xN boolean array of whether Ai is connected to 
                the chemical reaction network
        B_connected (arr): a 1xM boolean array of whether Bi is connected to 
                the chemical reaction network
        C_connected (arr): a 1xL boolean array of whether Ci is connected to 
                the chemical reaction network
    ABC_nonzero (list): a list of [A_nonzero, B_nonzero, C_nonzero] where              
        X_nonzero (arr): a 1D array of the indices of the connected X species
    duplex_nonzero (list): a list of [KAB_nonzero, KBC_nonzero, KAC_nonzero, 
        KAA_nonzero, KBB_nonzero, KCC_nonzero] where
        KXY_nonzero (tuple): a tuple of arrays of the indices of the nonzero 
            KAB_array values 
    
    N_reactants (int): the total number of reacting species in {A}, {B},
        {C}, {AB}, {BC}, {AC}, {AA}, {BB}, and {CC}.
    '''

    # Get the reacting species in {A}
    A_connected = get_X_connected(KAB_array, KAC_array, KAA_array)
    A_nonzero = np.nonzero(A_connected)[0]
    # Get the reacting species in {B}
    B_connected = get_X_connected(np.transpose(KAB_array), KBC_array, KBB_array)
    B_nonzero = np.nonzero(B_connected)[0]
    # Get the reacting species in {C}
    C_connected = get_X_connected(np.transpose(KAC_array), \
        np.transpose(KBC_array), KCC_array)
    C_nonzero = np.nonzero(C_connected)[0]
    ABC_connected = [A_connected, B_connected, C_connected]
    ABC_nonzero = [A_nonzero, B_nonzero, C_nonzero]
    
    # Get the number of reacting species in {A}, {B}, and {C}
    N_reactants = np.sum(A_connected) + np.sum(B_connected) + np.sum(C_connected)

    # ID the nonzero reactions and get the number of duplex species
    KAB_nonzero = np.nonzero(KAB_array)
    N_reactants += len(KAB_nonzero[0])
    KBC_nonzero = np.nonzero(KBC_array)
    N_reactants += len(KBC_nonzero[0])
    KAC_nonzero = np.nonzero(KAC_array)
    N_reactants += len(KAC_nonzero[0])
    
    # Handle the in-set interactions
    # Use the upper triangle of each association matrix to avoid double-counting
    KAA_nonzero = np.nonzero(np.triu(KAA_array))
    N_reactants += len(KAA_nonzero[0])
    KBB_nonzero = np.nonzero(np.triu(KBB_array))
    N_reactants += len(KBB_nonzero[0])
    KCC_nonzero = np.nonzero(np.triu(KCC_array))
    N_reactants += len(KCC_nonzero[0])
    
    duplex_nonzero = [KAB_nonzero, KBC_nonzero, KAC_nonzero, KAA_nonzero, \
    KBB_nonzero, KCC_nonzero]
    
    return (ABC_connected, ABC_nonzero, duplex_nonzero, N_reactants)

def gen_XY_bounds(Xi_array, Yi_array, KXY_nonzero):
    '''This function takes in arrays of initial concentrations for X and Y, and
    an array of the indices of the nonzero values of flattened KXY. It computes
    the bounds for the reactions that take place between Xi and Yj and returns 
    them as a list. 
    Parameters:
    Xi_array (arr): a 1xI array of initial concentrations of X
    Yi_array (arr): a 1xJ array of initial concentrations of Y
    KXY_nonzero (arr): a tuple of arrays of the indices of the nonzero KXY_array
     values.

    Returns:
    bounds_list (list): a list of (min,max) bound tuples. The list is the same 
        length as KXY_nonzero.
    '''
    bounds_list = []

    for index in range(len(KXY_nonzero[0])):
        i = KXY_nonzero[0][index]
        j = KXY_nonzero[1][index]
        bounds_list.append((0,min(Xi_array[i],Yi_array[j])))

    return bounds_list

def gen_AB_bounds(Ai_array, Bi_array, ABC_nonzero, duplex_nonzero):
    '''This function takes in an array of C initial concentrations and a 
    dictionary of settings, and creates bounds to use with Keq_fast for 
    scipy.optimize.minimize.

    Arguments:
    Ai_array (arr): a 1xN array of initial A concentrations
    Bi_array (arr): a 1xM array of initial B concentrations
    ABC_nonzero (list): [A_nonzero, B_nonzero, C_nonzero] where
        A_nonzero (arr): a 1D array of the indices of the connected A species
        B_nonzero (arr): a 1D array of the indices of the connected B species
        C_nonzero (arr): a 1D array of the indices of the connected C species
    duplex_nonzero (list): [KAB_nonzero, KBC_nonzero, KAC_nonzero, KAA_nonzero, 
        KBB_nonzero, KCC_nonzero] where
        KXY_nonzero (tuple): a tuple of arrays of the indices of the nonzero 
            KXY_array values. 

    Returns:
    fixed_bounds (list): [A_B_bounds, AB_bounds, AA_bounds, BB_bounds] where
        A_B_bounds (list): a list of (min,max) concentration bound tuples for 
            {A} and {B} to use with Keq_fast in scipy.optimize.minimize.
        XY_bounds (list): a list of (min,max) concentration bound tuples for 
            {XY} to use with Keq_fast in scipy.optimize.minimize. 
    '''
    
    # Create list to store {A} and {B} bounds in
    A_B_bounds = []

    # Collect initial concentrations
    init_concentrations = [Ai_array, Bi_array]

    # Add bounds for {A}, {B}, and {C}
    for i, init_conc_array in enumerate(init_concentrations):
        reactant_indices = ABC_nonzero[i]
        for index in reactant_indices:
            A_B_bounds.append((0,init_conc_array[index]))

    # Unpack the duplex nonzero indices
    KAB_nonzero, KBC_nonzero, KAC_nonzero, KAA_nonzero, KBB_nonzero, KCC_nonzero \
    = duplex_nonzero

    # Add bounds for {AB}, (AA), and {BB}
    AB_bounds = gen_XY_bounds(Ai_array, Bi_array, KAB_nonzero)
    AA_bounds = gen_XY_bounds(Ai_array, Ai_array, KAA_nonzero)
    BB_bounds = gen_XY_bounds(Bi_array, Bi_array, KBB_nonzero)

    fixed_bounds = [A_B_bounds, AB_bounds, AA_bounds, BB_bounds]
    return fixed_bounds

def add_C_bounds(Ai_array, Bi_array, Ci_array, fixed_bounds, ABC_nonzero, \
    duplex_nonzero):
    '''This function takes in an array of C initial concentrations and a 
    dictionary of settings, and creates bounds to use with Keq_fast for 
    scipy.optimize.minimize.

    Arguments:
    Ai_array (arr): a 1xN array of initial A concentrations
    Bi_array (arr): a 1xM array of initial B concentrations
    Ci_array (arr): a 1xL array of initial C concentrations
    fixed_bounds (list): [A_B_bounds, AB_bounds, AA_bounds, BB_bounds] where
        A_B_bounds (list): a list of (min,max) concentration bound tuples for 
            {A} and {B} to use with Keq_fast in scipy.optimize.minimize.
        XY_bounds (list): a list of (min,max) concentration bound tuples for 
            {XY} to use with Keq_fast in scipy.optimize.minimize. 
    ABC_nonzero (list): [A_nonzero, B_nonzero, C_nonzero] where
        A_nonzero (arr): a 1D array of the indices of the connected A species
        B_nonzero (arr): a 1D array of the indices of the connected B species
        C_nonzero (arr): a 1D array of the indices of the connected C species
    duplex_nonzero (list): [KAB_nonzero, KBC_nonzero, KAC_nonzero, KAA_nonzero, 
        KBB_nonzero, KCC_nonzero] where
        KXY_nonzero (tuple): a tuple of arrays of the indices of the nonzero 
            KXY_array values. 

    Returns:
    bounds_list (list): a list of (min,max) bound tuples to use with Keq_fast in
        scipy.optimize.minimize. The list is the same length as x_array.
    '''

    # Create list to store bounds in
    C_bounds = []

    # Unpack the duplex nonzero indices
    KAB_nonzero, KBC_nonzero, KAC_nonzero, KAA_nonzero, KBB_nonzero, KCC_nonzero \
    = duplex_nonzero

    # Get bounds for {C}
    for index in ABC_nonzero[2]:
        C_bounds.append((0,Ci_array[index]))

    # Get bounds for {BC}, {AC}, and {CC}
    BC_bounds = gen_XY_bounds(Bi_array, Ci_array, KBC_nonzero)
    AC_bounds = gen_XY_bounds(Ai_array, Ci_array, KAC_nonzero)
    CC_bounds = gen_XY_bounds(Ci_array, Ci_array, KCC_nonzero)

    bounds_list = fixed_bounds[0] + C_bounds + fixed_bounds[1] + BC_bounds \
    + AC_bounds + fixed_bounds[2] + fixed_bounds[3] + CC_bounds

    return bounds_list

def unpack_x_array(x_array, N, M, L, ABC_nonzero, duplex_nonzero):
    '''This function takes in x_array, set_sizes, and indices of nonzero species, 
    and unpacks x_array into concentration arrays Af_array, Bf_array, Cf_array, 
    AB_array, BC_array, AC_array, AA_array, BB_array, and CC_array.

    Parameters:
    N (int): the number of A species
    M (int): the number of B species
    L (int): the number of C species
    ABC_nonzero (list): [A_nonzero, B_nonzero, C_nonzero] where
        A_nonzero (arr): a 1D array of the indices of the connected A species
        B_nonzero (arr): a 1D array of the indices of the connected B species
        C_nonzero (arr): a 1D array of the indices of the connected C species
    duplex_nonzero (list): [KAB_nonzero, KBC_nonzero, KAC_nonzero, KAA_nonzero, 
        KBB_nonzero, KCC_nonzero] where
        KXY_nonzero (tuple): a tuple of arrays of the indices of the nonzero 
            KXY_array values. 

    Returns:
    Af_array (arr): a 1xN array of A concentrations
    Bf_array (arr): a 1xM array of B concentrations
    Cf_array (arr): a 1xL array of C concentrations
    AB_array (arr): a NxM array of AB concentrations
    BC_array (arr): a MxL array of BC concentrations
    AC_array (arr): a NxL array of AC concentrations
    AA_array (arr): a NxN array of AA concentrations
    BB_array (arr): a MxM array of BB concentrations
    CC_array (arr): a LxL array of CC concentrations
    '''
    
    # Create arrays to store the concentrations in.
    Af_array = np.zeros(N)
    Bf_array = np.zeros(M)
    Cf_array = np.zeros(L)
    AB_array = np.zeros((N,M))
    BC_array = np.zeros((M,L))
    AC_array = np.zeros((N,L))
    AA_array = np.zeros((N,N))
    BB_array = np.zeros((M,M))
    CC_array = np.zeros((L,L))

    # Organize the concentration arrays
    ABC_arrays = [Af_array, Bf_array, Cf_array]
    duplex_arrays = [AB_array, BC_array, AC_array, AA_array, BB_array, CC_array]

    # Unpack the x_array {A}, {B}, and {C} concentrations
    index = 0
    for i, reactant_indices in enumerate(ABC_nonzero):
        species_array = ABC_arrays[i]
        for reactant_index in reactant_indices:
            species_array[reactant_index] = x_array[index]
            index+=1
    
    # Unpack the x_array duplex concentrations
    for i, reactant_indices in enumerate(duplex_nonzero):
        species_array = duplex_arrays[i]
        for j in range(len(reactant_indices[0])):
            species_array[reactant_indices[0][j],reactant_indices[1][j]] = x_array[index]
            index+=1

    # Fill in the lower triangle of the XX_arrays
    AA_array = np.transpose(np.triu(AA_array,k=1))+AA_array
    BB_array = np.transpose(np.triu(BB_array,k=1))+BB_array
    CC_array = np.transpose(np.triu(CC_array,k=1))+CC_array

    return (Af_array, Bf_array, Cf_array, AB_array, BC_array, AC_array, \
        AA_array, BB_array, CC_array)

def Keq(x_array, Ai_array, Bi_array, Ci_array, set_sizes, K_array_list, \
    ABC_connected, ABC_nonzero, duplex_nonzero):
    '''Keq(x_array, Ai_array, Bi_array, Ci_array, set_sizes, K_array_list, 
           ABC_connected, ABC_nonzero, duplex_nonzero)
    The output of this function is intended to be used with 
    scipy.optimize.minimize. It represents a system of equations for a DNA 
    chemical reaction network with 3 sets of strand species ({A}, {B}, {C}), 
    which may have pairwise interactions to produce duplexes ({AB}, {BC, {AC}, 
    {AA}, {BB}, {CC}). For example, the ith X strand and the jth Y strand may 
    interact as Xi + Yj <-> XYij, with equilibrium constant KXYij. At 
    equilibrium this interaction satisfies KXYij[Xf][Yf] - [XY] = 0. This 
    conservation equality must be maintained: Xi - Xf - {XY} - {XZ} - {XX} = 0. 
    Concentrations must also obey the following two inequalitys: 0 <= Xf <= Xi
    and 0 <= XYij <= min(Xi, Yj). These inequalities are accounted for by 
    generating bounds with gen_bounds to be used with scipy.optimize.minimize.

    For every equilibrium equation with nonzero K, the function computes 
    (KXYij[Xf][Yf] - [XY])^2. For every conservation equation the function 
    computes (Xi - Xf - {XY} - {XZ} - {XX})^2. These are summed and output as 
    F_out.

    Parameters:
    x_array (arr): a 1D array of concentrations for species that react in the 
        CRN. [*A, *B, *C, *AB, *BC, *AC, *AA, *BB, *CC]
    Ai_array (arr): a 1xL array of initial A concentrations
    Bi_array (arr): a 1xL array of initial B concentrations
    Ci_array (arr): a 1xL array of initial C concentrations
    set_sizes (list): [N, M, L] where
        N (int): the number of A species
        M (int): the number of B species
        L (int): the number of C species
    K_array_list: [KAB_array, KBC_array, KAC_array, KAA_array, KBB_array, 
        KCC_array] where
        KXY_array (arr): a IxJ array of KXY association constants
    ABC_connected (list): [A_connected, B_connected, C_connected] where
        X_connected (arr): a 1D array of the indices of the connected A species
    ABC_nonzero (list): [A_nonzero, B_nonzero, C_nonzero] where
        A_nonzero (arr): a 1D array of the indices of the connected A species
        B_nonzero (arr): a 1D array of the indices of the connected B species
        C_nonzero (arr): a 1D array of the indices of the connected C species
    duplex_nonzero (list): [KAB_nonzero, KBC_nonzero, KAC_nonzero, KAA_nonzero, 
        KBB_nonzero, KCC_nonzero] where
        KXY_nonzero (tuple): a tuple of arrays of the indices of the nonzero 
            KXY_array values. 

    Returns: 
    F_out (float): the sum of the squares of the equilibrium and conservation
        equations. It should equal 0 at equilibrium.
    '''
    # Unpack set sizes
    N, M, L = set_sizes
    # Unpack K_arrays
    KAB_array, KBC_array, KAC_array, KAA_array, KBB_array, KCC_array = K_array_list
    # Unpack ABC_connected
    A_connected, B_connected, C_connected = ABC_connected

    # Unpack x_array
    Af_array, Bf_array, Cf_array, AB_array, BC_array, AC_array, AA_array, \
    BB_array, CC_array = unpack_x_array(x_array, N, M, L, ABC_nonzero, \
        duplex_nonzero)

    # Calculate the equilibrium equation contributions
    # Calculate the AB equations
    AB_F = calc_XY_equilibrium(Af_array, Bf_array, AB_array, KAB_array)
    # Calculate the BC equations
    BC_F = calc_XY_equilibrium(Bf_array, Cf_array, BC_array, KBC_array)
    # Calculate the AC equations
    AC_F = calc_XY_equilibrium(Af_array, Cf_array, AC_array, KAC_array)
    # Calculate the AA equations
    AA_F = calc_XX_equilibrium(Af_array, AA_array, KAA_array)
    # Calculate the AA equations
    BB_F = calc_XX_equilibrium(Bf_array, BB_array, KBB_array)
    # Calculate the AA equations
    CC_F = calc_XX_equilibrium(Cf_array, CC_array, KCC_array)

    # Calculate the conservation equation contributions
    # Calculate the A conservation equations
    A_conservation_F = calc_X_conservation(Ai_array, Af_array, AB_array, \
        AC_array, AA_array, A_connected)
    # Calculate the B conservation equations
    B_conservation_F = calc_X_conservation(Bi_array, Bf_array, \
        np.transpose(AB_array), BC_array, BB_array, B_connected)
    # Calculate the C conservation equations
    C_conservation_F = calc_X_conservation(Ci_array, Cf_array, \
        np.transpose(AC_array), np.transpose(BC_array), CC_array, \
        C_connected)

    # Sum all the equation contributions
    F_out = AB_F + BC_F + AC_F + AA_F + BB_F + CC_F + A_conservation_F \
    + A_conservation_F + B_conservation_F + C_conservation_F
    return F_out

def calc_A_measured(xf_array, Ai_array, N, M, ABC_nonzero, duplex_nonzero):
    '''
    This function takes in an array of initial A concentrations and an array of 
    AB concentrations, and computes Ai_meas = Ai_array - {AiB}. 

    Parameters:
    xf_array (arr): a 1D array of final concentrations for species that react in 
        the CRN. [*A, *B, *C, *AB, *BC, *AC, *AA, *BB, *CC]
    Ai_array (arr): a 1xL array of initial A concentrations
    N (int): the number of A species
    M (int): the number of B species
    ABC_nonzero (list): [A_nonzero, B_nonzero, C_nonzero] where
        A_nonzero (arr): a 1D array of the indices of the connected A species
        B_nonzero (arr): a 1D array of the indices of the connected B species
        C_nonzero (arr): a 1D array of the indices of the connected C species
    duplex_nonzero (list): [KAB_nonzero, KBC_nonzero, KAC_nonzero, KAA_nonzero, 
        KBB_nonzero, KCC_nonzero] where
        KXY_nonzero (tuple): a tuple of arrays of the indices of the nonzero 
            KXY_array values. 

    Returns:
    A_measured (arr): a 1xN of measured A concentrations
    '''

    # Get AB_array from xf_array
    # Create arrays to store the concentrations in.
    AB_array = np.zeros((N,M))
    
    # Get starting index of AB concentrations
    index = len(ABC_nonzero[0]) + len(ABC_nonzero[1]) + len(ABC_nonzero[2])
    
    # Unpack the x_array duplex concentrations
    for i in range(len(duplex_nonzero[0][0])):
        AB_array[duplex_nonzero[0][0][i],duplex_nonzero[0][1][i]] = xf_array[index]
        index+=1

    # Sum the AiB concentrations
    AiB_array = np.sum(AB_array, axis=1)
    # Compute A_measured
    A_measured = Ai_array-AiB_array

    return A_measured


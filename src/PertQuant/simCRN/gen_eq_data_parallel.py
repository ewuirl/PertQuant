import numpy as np
import time
import concurrent.futures
import random as rand
import math
import os
import argparse
import platform
from scipy.optimize import minimize
from PertQuant.simCRN.read_gen_eq_data_settings import read_array
from PertQuant.simCRN.read_gen_eq_data_settings import read_eq_data_settings
from PertQuant.simCRN.gen_eq_data import gen_Ci_array
from PertQuant.simCRN.gen_eq_data import compose_association_array
from PertQuant.simCRN.gen_eq_data import write_array
from PertQuant.simCRN.gen_eq_data import gen_detailed_eq_data_file
from PertQuant.simCRN.equilibrium_v3 import ID_network_species_reactions
from PertQuant.simCRN.equilibrium_v3 import gen_XY_bounds
from PertQuant.simCRN.equilibrium_v3 import gen_AB_bounds
from PertQuant.simCRN.equilibrium_v3 import add_C_bounds
from PertQuant.simCRN.equilibrium_v3 import Keq
from PertQuant.simCRN.equilibrium_v3 import calc_A_measured

def gen_A_measured(Cmin, Cmax, x_array, Ai_array, Bi_array, set_sizes, \
    K_array_list, ABC_connected, ABC_nonzero, duplex_nonzero, fixed_bounds):
    '''This function generates random Ci concentrations between [Cmin, Cmax] and
    uses the given args with the function Keq and scipy.optimize.minimize to
    find the equilibrium concentrations of the reactant species. A_measured for 
    the nth A strand is computed as An_measured = Ani - {AnB}. The Ci 
    concentrations are stored in Ci_array. The results are returned as 
    (Ci_array, A_measured).

    Parameters:
    Cmin (int): Minimum value of Ci
    Cmax (int): Maximum value of Ci
    x_array (arr): a 1D array of initial guesses for the equilibrium 
        concentrations of species that react in the CRN. 
        [*A, *B, *C, *AB, *BC, *AC, *AA, *BB, *CC]
    Ai_array (arr): a 1xL array of initial A concentrations
    Bi_array (arr): a 1xL array of initial B concentrations
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
    fixed_bounds (list): [A_B_bounds, AB_bounds, AA_bounds, BB_bounds] where
        A_B_bounds (list): a list of (min,max) concentration bound tuples for 
            {A} and {B} to use with Keq_fast in scipy.optimize.minimize.
        XY_bounds (list): a list of (min,max) concentration bound tuples for 
            {XY} to use with Keq_fast in scipy.optimize.minimize. 

    Returns:
    Ci_array: a 1xL array of initial C concentrations
    A_measured (arr): a 1xN of measured A concentrations
    '''

    # Figure out how many A and C strands there are
    N, M, L = set_sizes
    
    # Generate new Ci_array
    Ci_array = gen_Ci_array(L, Cmin, Cmax)

    # Generate the bounds
    bounds_list = add_C_bounds(Ai_array, Bi_array, Ci_array, fixed_bounds, \
        ABC_nonzero, duplex_nonzero)

    # Calculate x_array that solves Keq
    result = minimize(Keq, x_array, (Ai_array, Bi_array, Ci_array, set_sizes, \
            K_array_list, ABC_connected, ABC_nonzero, duplex_nonzero), \
        method='SLSQP', bounds=bounds_list)
    x_final = result.x

    # Calculate the measured A concentration
    A_measured = calc_A_measured(x_final, Ai_array, N, M, ABC_nonzero, \
        duplex_nonzero)

    return(Ci_array, A_measured)    

def gen_eq_data_parallel_main():
    rand.seed()
    start = time.perf_counter()

    # Argument parser
    parser = argparse.ArgumentParser(description='Generates eq_data ')    
    parser.add_argument('settings_file', type=str, help='The settings input file.')
    parser.add_argument('--save_folder', type=str, help='The path to the folder \
        to save data to. Defaults to current working directory')
    parser.add_argument('--quiet', type=int, help='Verbosity level. 0 for settings \
        and timing. 1 for no commandline output. Defaults to 0.')
    parser.add_argument('--time_file', type=str, help='File to append simulation time \
        data to.')
    args = parser.parse_args()

    # Check the platform
    if platform.system() == 'Windows':
        sep = '\\'
    else:
        sep = '/'
    
    # Parse the arguments
    settings_file = args.settings_file 
    settings_file_path, settings_file_name = os.path.split(settings_file)
    # Get save file name
    save_file_name = settings_file_name.removesuffix('_settings.txt')
    # Use provided save folder
    if args.save_folder:
        save_folder = args.save_folder
    # Use the folder of the input file
    else:
        save_folder = os.getcwd()
    if args.quiet:
        quiet = args.quiet
    else:
        quiet = 0

    if quiet==0:
        print(f'Reading in settings file {settings_file_name}')
    
    # Read in settings
    settings_dict = read_eq_data_settings(settings_file, quiet=quiet)
    
    # Unpack settings
    settings = settings_dict.keys()
    N = settings_dict['N']
    M = settings_dict['M']
    L = settings_dict['L']
    set_sizes = [N, M, L]

    N_runs = settings_dict['N_runs']
    Cmin = settings_dict['Cmin']
    Cmax = settings_dict['Cmax']

    # Make array to store results
    Ci_all_array = np.zeros((N_runs, L))
    Am_array = np.zeros((N_runs, N))

    # Get initial concentration arrays for {A} and {B}
    Ai_array = settings_dict['Ai_array']
    Bi_array = settings_dict['Bi_array']   

    # Unpack K arrays
    KAB_array = settings_dict['KAB_array']
    KBC_array = settings_dict['KBC_array']
    KAC_array = settings_dict['KAC_array']
    KAA_array = settings_dict['KAA_array']
    KBB_array = settings_dict['KBB_array']
    KCC_array = settings_dict['KCC_array']

    K_array_list = [KAB_array, KBC_array, KAC_array, KAA_array, KBB_array, \
    KCC_array]
    
    # Get connected {A}, {B}, {C} species, and nonzero (reacting) species
    ABC_connected, ABC_nonzero, duplex_nonzero, N_reactants = \
    ID_network_species_reactions(KAB_array, KBC_array, KAC_array, KAA_array, \
        KBB_array, KCC_array)

    # Generate {A}, {B}, {AB}, {AA}, and {BB} bounds
    fixed_bounds = gen_AB_bounds(Ai_array, Bi_array, ABC_nonzero, \
        duplex_nonzero)

    # Generate x_array initial guess
    x_array = np.zeros(N_reactants)

    # Compute Am
    if quiet==0:
        print('Computing Am.')
    N_range = range(N_runs)
    index = 0
    start = time.perf_counter()
    lap_time = start
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = {executor.submit(gen_A_measured, Cmin, Cmax, x_array, Ai_array, \
            Bi_array, set_sizes, K_array_list, ABC_connected, ABC_nonzero, \
            duplex_nonzero, fixed_bounds): i for i in N_range}
        for result in concurrent.futures.as_completed(results):
            Ci_all_array[index,:], Am_array[index,:] = result.result()
            index+=1
            if index % 100 == 0 and quiet == 0:
                lap_time = time.perf_counter() 
                print(f'Finished with {index}/{N_runs}\tTime elapsed: \
                    {lap_time-start} second(s)')

    # Save results
    if quiet==0:
        print('Saving results.')
    gen_detailed_eq_data_file(f'{save_folder}{sep}{save_file_name}_data.txt', Ci_all_array, \
        Am_array, Cmin, Cmax, Bi_array, Ai_array, KAB_array, KBC_array, \
        KAC_array, KAA_array=KAA_array, KBB_array=KBB_array, KCC_array=KCC_array)

    end=time.perf_counter()
    
    if args.time_file: 
        time_file = args.time_file.strip('"')
        with open(time_file,'a') as file:
            file.write(f'\n{settings_file_name}\t{end-start}')
    else:
        pass
        
    if quiet==0:
        print(f"Time elapsed: {end-start} second(s)")

if __name__ == '__main__':
    gen_eq_data_parallel_main()
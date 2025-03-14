import numpy as np
import time
import concurrent.futures
import random as rand
import math
import os
import argparse
import platform
import datetime
from scipy.optimize import minimize
from scipy.optimize import NonlinearConstraint
from PertQuant.simCRN.read_gen_eq_data_settings import read_array
from PertQuant.simCRN.read_gen_eq_data_settings import read_eq_data_settings
from PertQuant.simCRN.gen_eq_data import gen_Ci_array
from PertQuant.simCRN.gen_eq_data import compose_association_array
from PertQuant.simCRN.gen_eq_data import write_array
from PertQuant.simCRN.gen_eq_data import initiate_detailed_eq_data_file
from PertQuant.simCRN.gen_eq_data import append_data
from PertQuant.simCRN.equilibrium_v3 import ID_network_species_reactions
from PertQuant.simCRN.equilibrium_v3 import gen_base_x_array
from PertQuant.simCRN.equilibrium_v3 import fill_x_array
from PertQuant.simCRN.equilibrium_v3 import gen_XY_bounds
from PertQuant.simCRN.equilibrium_v3 import gen_AB_bounds
from PertQuant.simCRN.equilibrium_v3 import add_C_bounds
from PertQuant.simCRN.equilibrium_v3 import Keq
from PertQuant.simCRN.equilibrium_v3 import calc_A_measured
from PertQuant.simCRN.equilibrium_v3 import Keq_nonlinear_constraint
from PertQuant.simCRN.equilibrium_v3 import calculate_conservation_constraints

def gen_A_measured(Cmin, Cmax, x_array, Ai_array, Bi_array, set_sizes, \
    K_array_list, ABC_connected, ABC_nonzero, duplex_nonzero, ABC_guess,\
    fixed_bounds, method, guess, use_constraints):
    '''This function generates random Ci concentrations between [Cmin, Cmax] and
    uses the given args with the function Keq and scipy.optimize.minimize to
    find the equilibrium concentrations of the reactant species. A_measured for 
    the nth A strand is computed as An_measured = Ani - {AnB}. The Ci 
    concentrations are stored in Ci_array. The results are returned as 
    (Ci_array, A_measured).

    Parameters:
    Cmin (int): Minimum value of Ci
    Cmax (int): Maximum value of Ci
    base_x_array (arr): a 1D array of initial guesses for the equilibrium 
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
        X_connected (arr): a 1xI boolean array of whether Xi is connected to the 
            chemical reaction network
    ABC_nonzero (list): [A_nonzero, B_nonzero, C_nonzero] where
        X_nonzero (arr): a 1D array of the indices of the connected X species
    duplex_nonzero (list): [KAB_nonzero, KBC_nonzero, KAC_nonzero, KAA_nonzero, 
        KBB_nonzero, KCC_nonzero] where
        KXY_nonzero (tuple): a tuple of arrays of the indices of the nonzero 
            KXY_array values. 
    ABC_guess (list): [A_guess, B_guess, C_reactions] where
        X_guess (arr) = initial guesses for the final concentrations of X species.
        C_reactions (arr) = the number of reactions Ci participates in.
    fixed_bounds (list): [A_B_bounds, AB_bounds, AA_bounds, BB_bounds] where
        A_B_bounds (list): a list of (min,max) concentration bound tuples for 
            {A} and {B} to use with Keq_fast in scipy.optimize.minimize.
        XY_bounds (list): a list of (min,max) concentration bound tuples for 
            {XY} to use with Keq_fast in scipy.optimize.minimize. 
    method (str): The solver method to use with scipy.optimize.minimize.

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

    # Generate the initial guess:
    if guess:
        x_array = fill_x_array(N, M, L, x_array, Ci_array, ABC_guess, \
            K_array_list, duplex_nonzero)
    else:
        pass

    Keq_args = (Ai_array, Bi_array, Ci_array, set_sizes, K_array_list, \
        ABC_connected, ABC_nonzero, duplex_nonzero)

    # Calculate x_array that solves Keq
    if use_constraints:
        constraints={'type':'eq','fun': calculate_conservation_constraints, \
        'args': Keq_args}
        result = minimize(Keq_nonlinear_constraint, x_array, \
            Keq_args, method=method, bounds=bounds_list, constraints=constraints)
    else:
        
        result = minimize(Keq, x_array, Keq_args, method=method, \
            bounds=bounds_list)
    x_final = result.x

    # Calculate the measured A concentration
    A_measured = calc_A_measured(x_final, Ai_array, N, M, ABC_nonzero, \
        duplex_nonzero)

    return(Ci_array, A_measured)


def gen_eq_data_parallel_main():
    start = time.perf_counter()
    rand.seed()
    # Argument parser
    parser = argparse.ArgumentParser(description='Generates eq_data ')    
    parser.add_argument('settings_file', type=str, help='The settings input file.')
    parser.add_argument('--save_folder', type=str, help='The path to the folder \
        to save data to. Defaults to current working directory')
    parser.add_argument('--quiet', type=int, help='Verbosity level. 0 for settings \
        and timing. 1 for no commandline output. Defaults to 0.')
    parser.add_argument('--time_file', type=str, help='File to append simulation time \
        data to.')
    parser.add_argument('--method', type=str, help='The solver method to use with \
        scipy.optimize.minimize. Defaults to SLSQP.')
    parser.add_argument('--guess', type=bool, help='Whether to guess equally \
        distributed concentrations or zero concentrations. Defaults to False (zero).')
    parser.add_argument('--constraints', type=bool, help='Whether to use nonlinear \
        constraints. Defaults to False.')
    parser.add_argument('--overwrite', type=bool, help='Whether to overwrite \
        existing data file. Defaults to False.')
    parser.add_argument('--lap', type=int, help='Frequency to save data. Defaults \
        to 100 samples.')
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

    if args.method:
        method=args.method
    else:
        method='SLSQP'

    if args.guess:
        guess = args.guess 
        guess_name = '_guess'
    else:
        guess = False
        guess_name = ''

    if args.constraints:
        use_constraints = args.constraints 
        constraints_name = '_constraints'
    else:
        use_constraints = False
        constraints_name = ''

    if args.overwrite:
        overwrite = args.overwrite 
    else:
        overwrite = False

    if args.lap:
        lap=args.lap 
    else:
        lap=100
    
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

    # Generate the fixed portion of the x_array initial guess
    if guess:
        x_array, ABC_guess = gen_base_x_array(N, M, L, N_reactants, Ai_array, \
            Bi_array, K_array_list, duplex_nonzero)
        
    else:
        x_array = np.zeros(N_reactants)
        ABC_guess = [np.zeros(N), np.zeros(M), np.zeros(L)]

    # Prepare the data file
    data_file_name = f'{save_folder}{sep}{save_file_name}_data{guess_name}{constraints_name}.txt'
    try: 
        # Check if the data file exists
        with open(data_file_name, 'r') as data_file:
            header = True
            lines = data_file.readlines()
            header_size = 0
            index = 0
            # Check how much data has been generated
            while header:
                if '#' in lines[index]:
                    header_size += 1
                    index+=1
                else:
                    header=False
            dataset_size = len(lines)-header_size
            N_runs = N_runs - dataset_size
        if quiet==0:
            print(f'Found existing data file with {dataset_size} samples')
        # Overwrite the file if desired
        if overwrite:
            if quiet==0:
                print('Overwriting existing data file')
            initiate_detailed_eq_data_file(data_file_name, Cmin, Cmax, N_runs, \
                Ai_array, Bi_array, L, Ai_array, KAB_array, KBC_array, KAC_array, \
             KAA_array=KAA_array, KBB_array=KBB_array, KCC_array=KCC_array)
    # Inititate the file if it doesn't exist
    except:
        if quiet==0:
            print('Initiating data file')
        initiate_detailed_eq_data_file(data_file_name, Cmin, Cmax, N_runs, \
            Ai_array, Bi_array, L, KAB_array, KBC_array, KAC_array, \
            KAA_array=KAA_array, KBB_array=KBB_array, KCC_array=KCC_array)

    # Compute Am
    if quiet==0:
        print('Computing Am.')
    laps = [lap]*int(np.floor(N_runs/lap)) + [N_runs % lap]
    simulation_start = time.perf_counter()
    lap_time = simulation_start
    dataset_size = 0
    for lap in laps:
        N_range = range(lap)
        # Make array to store results
        Ci_all_array = np.zeros((lap, L))
        Am_array = np.zeros((lap, N))
        array_index=0
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = {executor.submit(gen_A_measured, Cmin, Cmax, x_array, Ai_array, \
                Bi_array, set_sizes, K_array_list, ABC_connected, ABC_nonzero, \
                duplex_nonzero, ABC_guess, fixed_bounds, method, guess, \
                use_constraints): i for i in N_range}
            for result in concurrent.futures.as_completed(results):
                Ci_all_array[array_index,:], Am_array[array_index,:] = result.result()
                array_index+=1
                dataset_size += 1
        # Append data to file
        append_data(data_file_name, Ci_all_array, Am_array)
        if quiet == 0:
            lap_time = time.perf_counter() 
            print(f'Finished with {dataset_size}/{N_runs}\tTime elapsed: \
                {str(datetime.timedelta(seconds=lap_time-simulation_start))}')
    
    time_file_name = f'{save_file_name}_data{guess_name}{constraints_name}.txt'
    end=time.perf_counter()
    
    if args.time_file: 
        time_file = args.time_file.strip('"')
        with open(time_file,'a') as file:
            file.write(f'\n{time_file_name}\t{N_runs} runs\t{end-simulation_start}')
    else:
        pass
        
    if quiet==0:
        print(f"Total time: {str(datetime.timedelta(seconds=end-start))}")

if __name__ == '__main__':
    gen_eq_data_parallel_main()
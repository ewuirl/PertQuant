from PertQuant.simCRN.equilibrium_v2 import calc_A_measured
from PertQuant.simCRN.equilibrium_v2 import gen_bounds
from PertQuant.simCRN.equilibrium_v2 import gen_bounds_A
from PertQuant.simCRN.equilibrium_v2 import gen_bounds_AB
from PertQuant.simCRN.equilibrium_v2 import gen_bounds_AA
from PertQuant.simCRN.equilibrium_v2 import Keq
from scipy.optimize import minimize
import numpy as np
import pandas as pd
import random as rand 

def gen_Ci_array(L, Cmin, Cmax):
    '''This function takes L, an integer number of strands in {C} and randomly 
    generates initial Ci concentrations betwen [Cmin,Cmax]. It returns these 
    concentrations as Ci_array, a numpy array.'''

    # Create an array to store the concentrations in 
    Ci_array = np.zeros(L)

    for i in range(L):
        Ci = rand.uniform(Cmin, Cmax)
        Ci_array[i] = Ci

    return(Ci_array)

def gen_Am(x_array, Cmin, Cmax, N_runs, *args):
    '''This function generates random Ci concentrations and uses this with given
    args to calculate the resulting Am concentrations. The Ci concentrations and
    Am concentrations are stored in Ci_all_array and Am_array, respectively. 
    This is repeated N_runs times, and the concentrations are returned in a 
    tuple (Ci_all_array, Am_array).'''
    # Figure out how many A and C strands there are
    N = len(args[0])
    L = len(args[3])

    # Create arrays to store values in
    Ci_all_array = np.zeros((N_runs, L))
    Am_array = np.zeros((N_runs, N))

    for i in range(N_runs):
        # Generate new Ci_array
        Ci_array = gen_Ci_array(L, Cmin, Cmax)

        # Make new args with the
        run_args = args[:3] + (Ci_array,) + args[4:]
    
        bounds_tuple = gen_bounds(*run_args)

        # Calculate x_array that solves Keq
        result = minimize(Keq, x_array, run_args, method='SLSQP', bounds=bounds_tuple)
        x_final = result.x

        # Calculate the measured A concentration
        Am = calc_A_measured(x_final, *run_args)

        Ci_all_array[i,:] = Ci_array
        Am_array[i] = Am

    return(Ci_all_array, Am_array)

def gen_eq_data_file(file_name, Ci_all_array, Am_array, Cmin, Cmax, Ai):
    '''This function takes data consisting of:
    Ci_all_array = N_runs x L array of initial C concentrations
    Am_array = N_runs x N array of A measured concentrations calculated based on
               the initial C conentrations.
    Cmin = the minimum concentration a strand C can take
    Cmax = the maximum concentration a strand C can take
    Ai = the initial concentration of strand A, equivalent to the max Am can be

    It writes the data to a file named file_name and prints a message when it is
    done. The first 6 rows give L (the number of C strands), N (the number of A
    strands), and N_runs (the number of computed runs), Cmin, Cmax, and Ai. Each
    row after that consists of the Ci values and Am values for a run, separated 
    by tabs.'''
    L = len(Ci_all_array[0])
    N = len(Am_array[0])
    N_runs = len(Ci_all_array)
    with open(file_name, 'w') as save_file:
        save_file.write('L = {:d} \n'.format(L))
        save_file.write('N = {:d} \n'.format(N))
        save_file.write('N_runs = {:d} \n'.format(N_runs))
        save_file.write('Cmin = {} \n'.format(Cmin))
        save_file.write('Cmax = {} \n'.format(Cmax))
        save_file.write('Ai = {} \n'.format(Ai))

        # Write the data to thefile
        for i in range(len(Ci_all_array)):
            # Pick out rows to write 
            Ci_row = Ci_all_array[i]
            Am_row = Am_array[i]

            for l in range(L):
                save_file.write(str(Ci_row[l]))
                save_file.write('\t')

            for n in range(N):
                save_file.write(str(Am_row[n]))
                save_file.write('\t')

            save_file.write('\n')

    print('Done writing file')

def write_array(save_file, array_name, array):
    save_file.write(f'\n# {array_name}')
    array_shape = np.shape(array)
    if np.size(array) > 0 and len(array_shape)==1:
        save_file.write('\n# '+'\t'.join(array.astype(str)))
    elif np.size(array) > 0 and len(array_shape)==2:
        for i in range(array_shape[0]):
            save_file.write('\n# '+'\t'.join(array[i,:].astype(str)))
    else:
        save_file.write('\n# ')

def compose_association_array(N, M, L, KAB_array, KBC_array, KAC_array, KAA_array, \
    KBB_array, KCC_array):
    total_strands = N+M+L
    association_array = np.zeros((total_strands,total_strands))
    # Set KAB
    association_array[:N,N:N+M] = KAB_array
    # Set KBC
    association_array[N:N+M,N+M:] = KBC_array
    # Set KAC
    association_array[:N,N+M:] = KAC_array
    
    # Set AA in-set interactions
    if len(KAA_array)> 0:
        for i in range(N):
            association_array[i,i:N] = KAA_array[i]
    else:
        pass
    # Set BB in-set interactions
    if len(KBB_array)> 0:
        for i in range(M):
            association_array[N+i,N+i:N+M] = KBB_array[i]
    else:
        pass        
    # Set CC in-set interactions
    if len(KCC_array)> 0:
        for i in range(L):
            association_array[N+M+i,N+M+i:N+M+L] = KCC_array[i]
    else:
        pass
    # Make array symmetric
    association_array += np.triu(association_array,k=1).transpose()
    return association_array

def gen_detailed_eq_data_file(file_name, Ci_all_array, Am_array, Cmin, Cmax, \
    Bi_array, Ai_array, KAB_array, KBC_array, KAC_array, KAA_array=[], \
    KBB_array=[], KCC_array=[], measured=True):
    '''This function takes data consisting of:
    Ci_all_array = N_runs x L array of initial C concentrations
    Am_array = N_runs x N array of A measured concentrations calculated based on
               the initial C conentrations.
    Cmin = the minimum concentration a strand C can take
    Cmax = the maximum concentration a strand C can take
    Ai = the initial concentration of strand A, equivalent to the max Am can be

    It writes the data to a file named file_name and prints a message when it is
    done. The first 6 rows give L (the number of C strands), N (the number of A
    strands), and N_runs (the number of computed runs), Cmin, Cmax, and Ai. Each
    row after that consists of the Ci values and Am values for a run, separated 
    by tabs.'''
    N = len(Am_array[0])
    M = len(Bi_array)
    L = len(Ci_all_array[0])
    N_runs = len(Ci_all_array)
    with open(file_name, 'w') as save_file:
        save_file.write(f'# N = {N:d} \n')
        save_file.write(f'# M = {M:d} \n')
        save_file.write(f'# L = {L:d} \n')
        save_file.write(f'# N_runs = {N_runs:d} \n')
        save_file.write(f'# Cmin = {Cmin} \n')
        save_file.write(f'# Cmax = {Cmax} \n')
        save_file.write(f'# Measured = {measured} \n# ')
        # Save the association constants
        write_array(save_file, 'Ai array', Ai_array)
        write_array(save_file, 'Bi array', Bi_array)
        # Save the association constants
        associaton_array = compose_association_array(N, M, L, KAB_array, \
            KBC_array, KAC_array, KAA_array, KBB_array, KCC_array)
        write_array(save_file, 'Association array', associaton_array)

        # Write the data to the file
        save_file.write('\n# \n# Data')
        for i in range(len(Ci_all_array)):
            # write the Ci and generated Am to a line
            save_file.write('\n'+'\t'.join(Ci_all_array[i,:].astype(str)))
            save_file.write('\t')
            save_file.write('\t'.join(Am_array[i,:].astype(str)))

    print('Done writing file')

def gen_Ap(x_array, x_array_C, Cmin, Cmax, N_runs, *args, verbose=True, 
    measured=True):
    '''This function takes in:
    x_array = a 1 x [N + M + NxM  + comb(N,2) + comb(M,2)] array of an initial 
        guess
    x_array_C = a 1 x [N + M + NxM + L + MxL + NxL + comb(N,2) + comb(M,2) 
        + comb(L,2)] of an initial guess
    Cmin = minimum Ci value
    Cmax = maximum Ci value
    N_run = number of equilibrium concentrations to calculate
    args = args for x_array_C to be used with Keq

    This function generates random Ci concentrations between [Cmin, Cmax] and
    uses the given args with the function Keq to calculate the resulting Am 
    concentrations for both x_array and x_array_C. It uses the given args to 
    create a corresponding set of args to be used with x_array. It calculates
    Ap = Am_C - Am, and stores these values in Ap_array. The corresponding 
    Ci concentrations are stored in Ci_all_array. This is repeated N_runs times,
    and the concentrations are returned in a tuple (Ci_all_array, Ap_array).'''\
    
    # Figure out how many A and C strands there are
    N = len(args[0])
    L = len(args[3])

    # Create arrays to store values in
    Ci_all_array = np.zeros((N_runs, L))
    Ap_array = np.zeros((N_runs, N))

    # Make args for the case with no C
    # Handle case with in-set interactions
    if len(args) == 9:
        run_args = args[:3] + args[6:8]
    # Handle case without in-set interactions
    else:
        run_args = args[:3]

    # Calculate x_array that solves Keq for the case with no C
    bounds_tuple = gen_bounds(*run_args)
    result = minimize(Keq, x_array, run_args, method='SLSQP', bounds=bounds_tuple)
    x_final = result.x
    Am = calc_A_measured(x_final, *run_args)
    if verbose:
        print(f'Am: {Am}')
    else:
        pass

    for i in range(N_runs):
        # Generate new Ci_array
        Ci_array = gen_Ci_array(L, Cmin, Cmax)

        # Make new args with the
        run_args_C = args[:3] + (Ci_array,) + args[4:]
    
        bounds_tuple_C = gen_bounds(*run_args_C)

        # Calculate x_array that solves Keq
        result_C = minimize(Keq, x_array_C, run_args_C, method='SLSQP', bounds=bounds_tuple_C)
        x_final_C = result_C.x

        # Calculate the measured A concentration
        Am_C = calc_A_measured(x_final_C, *run_args_C)
        if verbose:
            print(f'run: {i}/{N_runs}\tAm_C: {Am_C}')
        else:
            pass

        # Add the concentrations to the storage arrays
        Ci_all_array[i,:] = Ci_array
        if measured:
            Ap_array[i] = Am_C
        else:
            # Calculate the perturbation
            Ap = Am_C - Am
            Ap_array[i] = Ap
    return(Ci_all_array, Ap_array)
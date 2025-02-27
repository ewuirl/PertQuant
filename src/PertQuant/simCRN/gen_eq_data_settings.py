import numpy as np
import os
import argparse
import platform

def gen_KAB_array_base(N, M, KAB):
    assert N==M, 'N must be equal to M'
    return np.identity(N) * KAB

def gen_KBC_array_base(M, L):
    assert M>=L, 'M must be >= L'
    return np.zeros((M, L))

def gen_KBC_array_case_a(M, distinct, L, KBC):
    assert M>=L, 'M must be >= L'
    assert distinct>0, 'distinct must be > 0'
    KBC_array = np.zeros((M, L))
    KBC_array[:distinct,:distinct] = np.identity(distinct)*KBC
    return KBC_array

def gen_KBC_array_case_b(M, distinct, L, KBC):
    assert M>=L, 'M must be >= L'
    assert distinct>0, 'distinct must be > 0'
    KBC_array = np.zeros((M, L))
    KBC_array[1:1+distinct,:distinct] = np.identity(distinct)*KBC
    return KBC_array

def gen_KAC_array_case_1_all(N, L, KAC_vals):
    assert N>=L, 'N must be >= L'
    assert len(KAC_vals)==L, 'Length of KAC_vals must be equal to L'
    KAC_array = np.zeros((N,L))
    for i, KAC in enumerate(KAC_vals):
        KAC_array[i,:]=KAC*np.ones(L)
    return KAC_array

def gen_KAC_array_case_1_subset(N, L, KAC_scale, KAC_vals):
    assert N>=L, 'N must be >= L'
    assert 2<=len(KAC_vals)<L, 'Length of KAC_vals must be less than L, >=2'
    num_case=len(KAC_vals)
    case_1_array = np.ones((num_case,num_case))
    for i, KAC in enumerate(KAC_vals):
        case_1_array[i,:]=KAC*case_1_array[i,:]
    case_1_fill = np.identity(L-num_case)*KAC_scale
    KAC_array = np.zeros((N,L))
    KAC_array[:num_case,:num_case] = case_1_array
    KAC_array[num_case:L,num_case:L] = case_1_fill
    return KAC_array

def gen_KAC_array_case_2_all(N, L, column, scale):
    assert N>=L, 'N must be >= L'
    assert len(column)==N, 'Length of column must be equal to N'
    assert len(scale)==L, 'Length of scale must be equal to L'
    KAC_array = np.zeros((N,L))
    for i, KAC in enumerate(scale):
        KAC_array[:,i]=KAC*column
    return KAC_array

def gen_KAC_array_case_3_all(N, L, row):
    assert N>=L, 'N must be >= L'
    assert len(row)==L, 'Length of row must be equal to L'
    KAC_array = np.zeros((N,L))
    KAC_array[0,:] = row
    return KAC_array

def gen_KAC_array_case_4(N, L, KAC):
    assert N>=L, 'N must be >= L'
    if N == L:
        KAC_array = np.identity(N) * KAC
    else:
        KAC_array = np.zeros((N,L)) * KAC
        KAC_array[:L,:L] = np.identity(L)
    return KAC_array

def write_settings_array(array, file):
    shape = array.shape
    if len(shape)==1:
        file.write('\n')
        for i in array:
            file.write(f'{i} ')
    elif len(shape)==2:
        for i in range(shape[0]):
            file.write('\n')
            for j in range(shape [1]):
                file.write(f'{array[i,j]} ')
    else:
        pass

def write_settings_file(file_name, case_title, N, M, L, N_runs, Cmin, 
    Cmax, Ai_array, Bi_array, K_array_dict):
    '''This function takes in settings to use for generating data with 
    gen_eq_data_parallel.py and generates a settings file. 

    Parameters:
    file_name (str): the name to save the settings file with. Saves as 
        file_name_settings.txt
    case (str): the kind of case
    case_title (str): the title to use for the kind of case
    N (int): the number of {A} species
    M (int): the number of {B} species
    L (int): the number of {C} species
    N_runs (int): the number of data points to generate
    Cmin (float): the minimum concentration value of Ci
    Cmax (float): the minimum concentration value of Ci
    Ai_array (arr): a 1xN array of initial concentrations of {A}
    Bi_array (arr): a 1xM array of initial concentrations of {B}
    K_array_dict (dict): a dictionary with K_arrays as values and the K_array
        names as keys.
    '''
    with open(f'{file_name}_settings.txt', 'w') as file:
        if len(case_title) > 0:
            file.write(f'case_title = {case_title}')
        file.write(f'\nN = {N}')
        file.write(f'\nM = {M}')
        file.write(f'\nL = {L}')
        file.write(f'\nN_runs = {N_runs}')
        file.write(f'\nCmin = {Cmin}')
        file.write(f'\nCmax = {Cmax}')
        file.write(f'\nAi_array')
        write_settings_array(Ai_array, file)
        file.write(f'\nBi_array')
        write_settings_array(Bi_array, file)
        for K_array in K_array_dict:
            if np.sum(K_array_dict[K_array]) > 0:
                file.write(f'\n{K_array}')
                write_settings_array(K_array_dict[K_array], file)

if __name__ == '__main__':
    # Argument parser
    parser = argparse.ArgumentParser(description='Generates eq_data')    
    parser.add_argument('save_folder', type=str, help='The path to the save_folder')
    args = parser.parse_args()
    
    # Check the platform
    if platform.system() == 'Windows':
        sep = '\\'
    else:
        sep = '/'

    save_folder = args.save_folder

    N = 20
    M = 20
    L = 20
    N_runs = 5000
    Cmin = 0
    Cmax = 2
    KAB = 100
    KAC = 1000
    # Default initial concentrations
    DB_init=1
    Ai_array = np.ones(N)*DB_init
    Bi_array = np.ones(M)*DB_init
    # Default KAB and KBC
    K_array_dict = {}
    K_array_dict['KAB_array'] = gen_KAB_array_base(N, M, KAB)

    KBC_list = []
    case_title_list = []
    folder_list = []

    # # Case 1
    # KAC_vals = np.linspace(1000, 2000, L)
    # K_array_dict['KAC_array'] = gen_KAC_array_case_1_all(N, L, KAC_vals)
    # # base case
    # KBC_list.append(np.zeros((M,L)))
    # case_title_list.append('Case 1')
    # folder_list.append(f'{N}-{M}-{L}_case-1')
    # folder_list.append(f'{N}-{M}-{L}_case-1_DB-{DB_init}')
    # # KBC = 1, 1 distinct
    # KBC = 1
    # distinct = 1
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append(f'Case 1a, KBC = {KBC}, 1 distinct')
    # # folder_list.append(f'{N}-{M}-{L}_case-1a_KBC-{KBC}_1d')
    # folder_list.append(f'{N}-{M}-{L}_case-1a_KBC-{KBC}_DB-{DB_init}_1d')
    # # KBC = 1, all distinct
    # distinct = L-1
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append(f'Case 1a, KBC = {KBC}, all distinct')
    # folder_list.append(f'{N}-{M}-{L}_case-1a_KBC-{KBC}_all')
    # # folder_list.append(f'{N}-{M}-{L}_case-1a_KBC-{KBC}_DB-{DB_init}_all')
    # # KBC = 1, L distinct
    # distinct = L
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append('Case 1a, KBC = 1, all KBC')
    # # folder_list.append(f'{N}-{M}-{L}_case-1a_KBC-{KBC}_{L}d')
    # folder_list.append(f'{N}-{M}-{L}_case-1a_KBC-{KBC}_DB-{DB_init}_{L}d')
    # # KBC = 100, 1 distinct
    # KBC = 100
    # distinct = 1
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append('Case 1a, KBC = 100, 1 distinct')
    # # folder_list.append(f'{N}-{M}-{L}_case-1a_KBC-{KBC}_1d')
    # folder_list.append(f'{N}-{M}-{L}_case-1a_KBC-{KBC}_DB-{DB_init}_1d')
    # KBC = 100, all distinct
    # distinct = L-1
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append('Case 1a, KBC = 1, all distinct')
    # # folder_list.append(f'{N}-{M}-{L}_case-1a_KBC-{KBC}_all')
    # folder_list.append(f'{N}-{M}-{L}_case-1a_KBC-{KBC}_DB-{DB_init}_all')
    # # KBC = 100, 5 distinct
    # distinct = L
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append('Case 1a, KBC = 100, all KBC')
    # # folder_list.append(f'{N}-{M}-{L}_case-1a_KBC-{KBC}_{L}d')
    # folder_list.append(f'{N}-{M}-{L}_case-1a_KBC-{KBC}_DB-{DB_init}_{L}d')

    # # # Case 2 
    # KAC_vals = np.linspace(1,2,L)
    # column = np.linspace(100,1000,N)
    # K_array_dict['KAC_array'] = gen_KAC_array_case_2_all(N, L, column, KAC_vals)
    # # base case
    # KBC_list.append(np.zeros((M,L)))
    # case_title_list.append('Case 2')
    # folder_list.append(f'{N}-{M}-{L}_case-2')
    # # KBC = 1, 1 distinct
    # KBC = 100
    # # KBC = 10
    # distinct = 1
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append(f'Case 2a, KBC = {KBC}, 1 distinct')
    # folder_list.append(f'{N}-{M}-{L}_case-2a_KBC-{KBC}_1d')
    # # KBC = 1, all distinct
    # distinct = L-1
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append(f'Case 2a, KBC = {KBC}, all distinct')
    # folder_list.append(f'{N}-{M}-{L}_case-2a_KBC-{KBC}_all')
    # # KBC = 1, all L distinct
    # distinct = L
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append(f'Case 2a, KBC = {KBC}, all KBC')
    # folder_list.append(f'{N}-{M}-{L}_case-2a_KBC-{KBC}_{L}d')
    # KBC = 100, 1 distinct
    # KBC = 100
    # distinct = 1
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append('Case 2a, KBC = 100, 1 distinct')
    # folder_list.append(f'{N}-{M}-{L}_case-2a_KBC-100_1d')
    # # KBC = 100, all distinct
    # distinct = L-1
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append('Case 2a, KBC = 100, all distinct')
    # folder_list.append(f'{N}-{M}-{L}_case-2a_KBC-100_all')
    # # KBC = 100, all L distinct
    # distinct = L
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append('Case 2a, KBC = 100, all KBC')
    # folder_list.append(f'{N}-{M}-{L}_case-2a_KBC-100_{L}d')

    # # # Case 3
    # row = np.linspace(1000,2000,L)
    # K_array_dict['KAC_array'] = gen_KAC_array_case_3_all(N, L, row)
    # base case
    # KBC_list.append(np.zeros((M,L)))
    # case_title_list.append('Case 3')
    # folder_list.append(f'{N}-{M}-{L}_case-3')
    # KBC = 1
    # KBC = 10
    # distinct = 1
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append('Case 3a, KBC = 1')
    # folder_list.append(f'{N}-{M}-{L}_case-3a_KBC-{KBC}')
    # # KBC = 1, 1 distinct
    # distinct = 1
    # KBC_list.append(gen_KBC_array_case_b(M, distinct, L, KBC))
    # case_title_list.append('Case 3b, KBC = 1, 1 distinct')
    # folder_list.append(f'{N}-{M}-{L}_case-3b_KBC-{KBC}_1d')
    # # KBC = 1, all distinct
    # distinct = L-1
    # KBC_list.append(gen_KBC_array_case_b(M, distinct, L, KBC))
    # case_title_list.append('Case 3b, KBC = 1, all distinct')
    # folder_list.append(f'{N}-{M}-{L}_case-3b_KBC-{KBC}_all')
    # # KBC = 1, all KBC
    # distinct = L-1
    # KBC_array = gen_KBC_array_case_b(M, distinct, L, KBC)
    # KBC_array[0,-1] = KBC
    # KBC_list.append(KBC_array)
    # case_title_list.append('Case 3b, KBC = 1, all KBC')
    # folder_list.append(f'{N}-{M}-{L}_case-3b_KBC-{KBC}_{L}d')
    # # KBC = 100
    # KBC = 1
    # distinct = 1
    # KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    # case_title_list.append('Case 3a, KBC = 100')
    # folder_list.append(f'{N}-{M}-{L}_case-3a_KBC-100')
    # # KBC = 100, 1 distinct
    # distinct = 1
    # KBC_list.append(gen_KBC_array_case_b(M, distinct, L, KBC))
    # case_title_list.append('Case 3b, KBC = 100, 1 distinct')
    # folder_list.append(f'{N}-{M}-{L}_case-3b_KBC-100_1d')
    # # KBC = 100, all distinct
    # distinct = L-1
    # KBC_list.append(gen_KBC_array_case_b(M, distinct, L, KBC))
    # case_title_list.append('Case 3b, KBC = 100, all distinct')
    # folder_list.append(f'{N}-{M}-{L}_case-3b_KBC-100_all')
    # # KBC = 100, all KBC
    # distinct = L-1
    # KBC_array = gen_KBC_array_case_b(M, distinct, L, KBC)
    # KBC_array[0,-1] = KBC
    # KBC_list.append(KBC_array)
    # case_title_list.append('Case 3b, KBC = 100, all KBC')
    # folder_list.append(f'{N}-{M}-{L}_case-3b_KBC-100_{L}d')

    # # # Case 4
    # K_array_dict['KAC_array'] = gen_KAC_array_case_4(N, L, KAC)
    # case_title = 'Case 4'
    # folder = f'{N}-{M}-{L}_case-4'

    # # # Case 5
    # K_array_dict['KAC_array'] = np.ones((N,L))*1000
    # scale_arr = np.array([1,2.5,5])
    # order_arr = np.logspace(start=-2,stop=3, num=6)
    # # KBC_vals = np.array([0, 0.75])
    # KBC_vals = np.array([0.75])
    # for order in order_arr:
    #     KBC_vals = np.concatenate((KBC_vals, order*scale_arr))
    # KBC_vals = np.concatenate((KBC_vals,np.array([10000])))
    # distinct=L-1
    # for KBC in KBC_vals:
    #     KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    #     if np.floor(KBC)==KBC:
    #         case_title = f'Case 5 (K={int(KBC)})'
    #         folder =f'{N}-{M}-{L}_case-5_K-{int(KBC)}'
    #     else:
    #         case_title = f'Case 5 (K={KBC})'
    #         folder=f'{N}-{M}-{L}_case-5_K-{KBC}'.replace('.','-')
    #     case_title_list.append(case_title)
    #     folder_list.append(folder)
    # distinct=L
    # for KBC in KBC_vals:
    #     KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    #     if np.floor(KBC)==KBC:
    #         case_title = f'Case 5c (K={int(KBC)})'
    #         folder =f'{N}-{M}-{L}_case-5c_K-{int(KBC)}'
    #     else:
    #         case_title = f'Case 5 (K={KBC})'
    #         folder=f'{N}-{M}-{L}_case-5c_K-{KBC}'.replace('.','-')
    #     case_title_list.append(case_title)
    #     folder_list.append(folder)

    # # # Case 5b
    # Ai_array = np.ones(N)*2
    # Ai_array[-1] = 1
    # Bi_array = np.ones(M)*2
    # Bi_array[-1] = 1
    # K_array_dict['KAC_array'] = np.ones((N,L))*1000
    # scale_arr = np.array([1,2.5,5])
    # order_arr = np.array([100, 1000])
    # KBC_vals = np.array([])
    # for order in order_arr:
    #     KBC_vals = np.concatenate((KBC_vals, order*scale_arr))
    # KBC_vals = np.concatenate((KBC_vals,np.array([10000])))
    # distinct=4
    # for KBC in KBC_vals:
    #     KBC_list.append(gen_KBC_array_case_a(M, distinct, L, KBC))
    #     if np.floor(KBC)==KBC:
    #         case_title = f'Case 5b (K={int(KBC)})'
    #         folder =f'{N}-{M}-{L}_case-5b_K-{int(KBC)}'
    #     else:
    #         case_title = f'Case 5b (K={KBC})'
    #         folder=f'{N}-{M}-{L}_case-5b_K-{KBC}'.replace('.','-')
    #     case_title_list.append(case_title)
    #     folder_list.append(folder)


    # Case 6
    K_array_dict['KAC_array'] = np.zeros((N, L))
    scale_arr = np.array([1,2.5,5])
    order_arr = np.logspace(start=-2,stop=3, num=6)
    # KBC_vals = np.array([0, 0.75])
    KBC_vals = np.array([0.75])
    for order in order_arr:
        KBC_vals = np.concatenate((KBC_vals, order*scale_arr))
    KBC_vals = np.concatenate((KBC_vals,np.array([10000])))
    for KBC in KBC_vals:
        KBC_list.append(gen_KAC_array_case_4(M, L, KBC))
        if np.floor(KBC)==KBC:
            case_title = f'Case 6 (K={int(KBC)})'
            folder =f'{N}-{M}-{L}_case-6_K-{int(KBC)}'
        else:
            case_title = f'Case 6 (K={KBC})'
            folder=f'{N}-{M}-{L}_case-6_K-{KBC}'.replace('.','-')
        case_title_list.append(case_title)
        folder_list.append(folder)

    # # For single case
    # if folder not in os.listdir(save_folder):
    #     os.mkdir(f'{save_folder}{sep}{folder}')
    
    # file_name = f'{save_folder}{sep}{folder}{sep}{folder}'
    # write_settings_file(file_name, case_title, N, M, L, N_runs, Cmin, 
    #     Cmax, Ai_array, Bi_array, K_array_dict)

    # For case variations
    for i in range(len(KBC_list)):
        folder = folder_list[i]
        if folder not in os.listdir(save_folder):
            os.mkdir(f'{save_folder}{sep}{folder}')
        K_array_dict['KBC_array'] = KBC_list[i]
        case_title = case_title_list[i]
        file_name = f'{save_folder}{sep}{folder}{sep}{folder}'
        write_settings_file(file_name, case_title, N, M, L, N_runs, Cmin, 
            Cmax, Ai_array, Bi_array, K_array_dict)
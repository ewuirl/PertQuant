import numpy as np
import os
import math
import time
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import argparse
import pickle
from PertQuant.simCRN.multivariate_reg_v2 import read_detailed_eq_data_file
from PertQuant.simCRN.multivariate_reg_v2 import get_partitioned_data
from PertQuant.simCRN.multivariate_reg_v2 import plot_raw_data
from PertQuant.simCRN.multivariate_reg_v2 import predict_model

import sklearn
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import r2_score
from sklearn.metrics import mean_absolute_error

def gen_Ti_specific_metrics_main():
    # Argument parser
    parser = argparse.ArgumentParser(description='Generates eq_data ')    
    parser.add_argument('data_file_name', type=str, help='The name of the data_file \
        to use')
    parser.add_argument('model_type', type=str, help='The ML regularization model.')
    parser.add_argument('--MAE', type=str, help='Whether to generate MAE metrics. Defaults \
        to False.')
    parser.add_argument('--r2', type=str, help='Whether to generate r2 metrics. Defaults \
        to False.')
    
    args = parser.parse_args()

    # Parse the arguments
    data_file_name = args.data_file_name
    case_name, remainder = data_file_name.split('_data')
    suffix = remainder.removesuffix('.txt')

    model_type = args.model_type
    parent_folder = os.getcwd()

    if args.MAE:
        MAE=args.MAE
    else:
        MAE=False

    if args.r2:
        r2=args.r2
    else:
        r2=False

    if MAE or r2:
        # Read in data
        settings_dict = read_detailed_eq_data_file(f'{case_name}_data{suffix}.txt')
        N = settings_dict['N']
        M = settings_dict['M']
        L = settings_dict['L']
        Cmin = settings_dict['Cmin']
        Cmax = settings_dict['Cmax']
        case = case_name.removeprefix(f'{N}-{M}-{L}_')
        folder = f'{case}_{model_type}{suffix}'

        # Get partitioned data
        A_out_array = settings_dict['A_out_array']
        Ci_array = settings_dict['Ci_array']
        dataset_csv = pd.read_csv(f'{N}-{M}-{L}_{case}_dataset_partitions{suffix}.csv', index_col=0)
        data_dict = get_partitioned_data(A_out_array, Ci_array, dataset_csv)

        # Get regression results
        results_df = pd.read_csv(f'{N}-{M}-{L}_{case}_results_summary_{model_type}{suffix}.csv', index_col=0)
        dataset_size_list = results_df['dataset_size'].to_numpy()

        # Get trained models
        model_list = []
        for i, dataset_size in enumerate(dataset_size_list):
            model_file_name = f'{folder}/{N}-{M}-{L}_{case}_{model_type}_{dataset_size}{suffix}.pkl'
            
            with open(model_file_name, 'rb') as model_file:
                model = pickle.load(model_file)
                model_list.append(model)

        # Predict with trained models
        pred_train_list = []
        pred_test_list = []
        for i, dataset_size in enumerate(dataset_size_list):
            pred_train_list.append(model_list[i].predict(data_dict[dataset_size][0]))
            pred_test_list.append(model_list[i].predict(data_dict['test_set'][0]))
    else:
        print('No metric provided')

    Ti_columns = ['dataset_size'] + [f'Train T{i+1}' for i in range(L)] + [f'Test T{i+1}' for i in range(L)]
    if MAE:
        # Get Ti-specific MAE
        print('Getting Ti-specific MAE')
        
        MAE_df = pd.DataFrame(columns=Ti_columns)
        MAE_df['dataset_size'] = dataset_size_list
        for i, dataset_size in enumerate(dataset_size_list):
            MAE_df.loc[i,'Train T1':f'Train T{L}'] = mean_absolute_error(data_dict[dataset_size][1], pred_train_list[i], \
                                                          multioutput='raw_values')
            MAE_df.loc[i,'Test T1':f'Test T{L}'] = mean_absolute_error(data_dict['test_set'][1], pred_test_list[i], \
                                                          multioutput='raw_values')
        MAE_df.to_csv(f'{N}-{M}-{L}_{case}_MAE_Ti_{model_type}{suffix}.csv')
    else:
        pass 

    if r2:
        # Get Ti-specific r2
        print('Getting Ti-specific r2')
        r2_df = pd.DataFrame(columns=Ti_columns)
        r2_df['dataset_size'] = dataset_size_list
        for i, dataset_size in enumerate(dataset_size_list):
            r2_df.loc[i,'Train T1':f'Train T{L}'] = r2_score(data_dict[dataset_size][1], pred_train_list[i], \
                                                          multioutput='raw_values')
            r2_df.loc[i,'Test T1':f'Test T{L}'] = r2_score(data_dict['test_set'][1], pred_test_list[i], \
                                                          multioutput='raw_values')
        r2_df.to_csv(f'{N}-{M}-{L}_{case}_r2_Ti_{model_type}{suffix}.csv')
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
from PertQuant.simCRN.multivariate_reg_v2 import plot_true_vs_pred
from PertQuant.simCRN.multivariate_reg_v2 import plot_residuals
from PertQuant.simCRN.multivariate_reg_v2 import plot_metric_summary
from PertQuant.simCRN.multivariate_reg_v2 import plot_metric_best_vs_alt

import sklearn
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import r2_score
from sklearn.metrics import mean_absolute_error

def gen_regularization_plots_main():
    # Argument parser
    parser = argparse.ArgumentParser(description='Generates plots for regression \
        results.')    
    parser.add_argument('data_file_name', type=str, help='The name of the data_file \
        to use')
    parser.add_argument('case_title', type=str, help='The name to describe the case \
        in plot titles')
    parser.add_argument('model_type', type=str, help='The ML regularization model.')
    parser.add_argument('--max_cols', type=int, help='The maximum number of columns \
        to plot subplots with. Defaults to 4.')
    parser.add_argument('--show_plot', type=bool, help='Whether to not show the plots. \
        Defaults to False.')
    parser.add_argument('--no_save', type=bool, help='Whether to not save the plots. \
        Defaults to False (save).')
    args = parser.parse_args()

    # Parse the arguments
    data_file_name = args.data_file_name
    case_name, remainder = data_file_name.split('_data')
    suffix = remainder.removesuffix('.txt')

    case_title = args.case_title
    model_type = args.model_type
    parent_folder = os.getcwd()

    if args.max_cols:
        max_cols=args.max_cols
    else:
        max_cols=4

    if args.show_plot:
        show_plot = args.show_plot
    else:
        show_plot = False

    if args.no_save:
        save_plot = not args.no_save
    else:
        save_plot = True


    # plotting settings
    mpl.rcParams['font.size']=18

    # Read in data
    settings_dict = read_detailed_eq_data_file(f'{case_name}_data{suffix}.txt')
    N = settings_dict['N']
    M = settings_dict['M']
    L = settings_dict['L']
    Cmin = settings_dict['Cmin']
    Cmax = settings_dict['Cmax']
    case = case_name.removeprefix(f'{N}-{M}-{L}_')
    folder = f'{case}_{model_type}{suffix}'

    # Take a look at the raw data
    print('Plotting raw data')
    title=f'D=B=T={N} {case_title} Raw Data'
    save_file = f'{N}-{M}-{L}_{case}_raw-data{suffix}'
    if N==5:
        y_raw_title=0.92
    elif N==10:
        y_raw_title=0.89
    else:
        y_raw_title=0.9
    raw_data_plot = plot_raw_data(settings_dict, f'D=B=T={N} {case_title} Raw Data', y_raw_title, save_file=save_file, \
                                  save=save_plot)
    if show_plot:
        plt.show(block=True)

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

    # Plot true vs pred and residuals
    print('Plotting true vs predicted Ti and residuals')
    for i, dataset_size in enumerate(dataset_size_list):
        dataset_size = dataset_size_list[i]
        train_title = f'{case_title}: True vs Predicted Train $T_i$, Best {model_type} with $N_s$={dataset_size}' 
        train_save = f'{folder}/{N}-{M}-{L}_{case}_train_true_vs_pred_best_{model_type}_{dataset_size}{suffix}'
        plot_true_vs_pred(data_dict[dataset_size][1], pred_train_list[i], Cmin, Cmax, \
            train_title, train_save, max_cols=max_cols, save=save_plot)
        test_title = f'{case_title}: True vs Predicted Test $T_i$, Best {model_type} with $N_s$={dataset_size}'
        test_save = f'{folder}/{N}-{M}-{L}_{case}_test_true_vs_pred_best_{model_type}_{dataset_size}{suffix}'
        plot_true_vs_pred(data_dict['test_set'][1], pred_test_list[i], Cmin, Cmax, \
            test_title, test_save, max_cols=max_cols, save=save_plot, color='orange')
        train_title = f'{case_title}: Predicted Train $T_i$ Residuals, Best {model_type} with $N_s$={dataset_size}' 
        train_save = f'{folder}/{N}-{M}-{L}_{case}_train_residuals_best_{model_type}_{dataset_size}{suffix}'
        plot_residuals(data_dict[dataset_size][1], pred_train_list[i], Cmin, Cmax, \
            train_title, train_save, max_cols=max_cols, save=save_plot)
        test_title = f'{case_title}: Predicted Test $T_i$ Residuals, Best {model_type} with $N_s$={dataset_size}' 
        test_save = f'{folder}/{N}-{M}-{L}_{case}_test_residuals_best_{model_type}_{dataset_size}{suffix}'
        plot_residuals(data_dict['test_set'][1], pred_test_list[i], Cmin, Cmax, \
            test_title, test_save, max_cols=max_cols, save=save_plot, color='orange')
        if show_plot:
            plt.show(block=True)
        plt.close('all')

    # Plot the overall R^2 results
    print('Plotting R^2 and MAE summary plots')
    r2_save_file=f'{N}-{M}-{L}_{case}_{model_type}_summary_plot_all_r2{suffix}'
    plot_metric_summary(results_df, case_title, model_type, 'r2', r2_save_file, save=save_plot)

    # Plot the overall MAE results
    MAE_save_file=f'{N}-{M}-{L}_{case}_{model_type}_summary_plot_all_mae{suffix}'
    plot_metric_summary(results_df, case_title, model_type, 'neg_mean_absolute_error', MAE_save_file, save=save_plot)

    if show_plot:
        plt.show(block=True)

    # Get Ti-specific MAE
    print('Getting Ti-specific MAE')
    MAE_columns = ['dataset_size'] + [f'Train T{i+1}' for i in range(L)] + [f'Test T{i+1}' for i in range(L)]
    MAE_df = pd.DataFrame(columns=MAE_columns)
    MAE_df['dataset_size'] = dataset_size_list
    for i, dataset_size in enumerate(dataset_size_list):
        MAE_df.loc[i,'Train T1':f'Train T{L}'] = mean_absolute_error(data_dict[dataset_size][1], pred_train_list[i], \
                                                      multioutput='raw_values')
        MAE_df.loc[i,'Test T1':f'Test T{L}'] = mean_absolute_error(data_dict['test_set'][1], pred_test_list[i], \
                                                      multioutput='raw_values')
    MAE_df.to_csv(f'{N}-{M}-{L}_{case}_MAE_Ti_{model_type}{suffix}.csv')
    MAE_df
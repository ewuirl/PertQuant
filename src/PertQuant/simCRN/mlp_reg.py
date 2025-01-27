import numpy as np
import os
import math
import time
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import pickle
import argparse
import platform
from contextlib import redirect_stdout
from PertQuant.simCRN.multivariate_reg_v2 import read_detailed_eq_data_file
from PertQuant.simCRN.multivariate_reg_v2 import plot_raw_data
from PertQuant.simCRN.multivariate_reg_v2 import partition_data
from PertQuant.simCRN.multivariate_reg_v2 import save_grid_search_results
from PertQuant.simCRN.multivariate_reg_v2 import best_model_predict
from PertQuant.simCRN.multivariate_reg_v2 import best_model_metrics
from PertQuant.simCRN.multivariate_reg_v2 import train_model
from PertQuant.simCRN.multivariate_reg_v2 import predict_model
from PertQuant.simCRN.multivariate_reg_v2 import get_best_metric_model
from PertQuant.simCRN.multivariate_reg_v2 import add_predictions_to_df
from PertQuant.simCRN.multivariate_reg_v2 import get_alt_model_results
from PertQuant.simCRN.multivariate_reg_v2 import compare_alt_model
from PertQuant.simCRN.multivariate_reg_v2 import plot_true_vs_pred
from PertQuant.simCRN.multivariate_reg_v2 import plot_residuals
from PertQuant.simCRN.multivariate_reg_v2 import plot_metric_summary
from PertQuant.simCRN.multivariate_reg_v2 import plot_metric_best_vs_alt

import sklearn
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import r2_score
from sklearn.metrics import mean_absolute_error

def mlp_reg_main():
    # Argument parser
    parser = argparse.ArgumentParser(description='Runs MLP regression')    
    parser.add_argument('data_file', type=str, help='The settings input file.')
    parser.add_argument('--max_iter', type=int, help='The value of max_iter to use \
        with MLPRegressor. Defaults to 1000.')
    parser.add_argument('--save_folder', type=str, help='The path to the parent \
        folder to save data to. Defaults to the folder of the input file.')
    parser.add_argument('--raw_data', type=bool, help='If True, plots and shows \
    the raw data and exits. Defaults to False.')
    # parser.add_argument('--quiet', type=int, help='Verbosity level. 0 for settings \
        # and timing. 1 for no commandline output. Defaults to 0.')
    args = parser.parse_args()

    # Check the platform
    if platform.system() == 'Windows':
        sep = '\\'
    else:
        sep = '/'

    # Parse the arguments
    data_file = args.data_file
    data_file_path, data_file_name = os.path.split(data_file)
    save_name = data_file_name.removesuffix('_data.txt')
    if args.max_iter:
        max_iter = args.max_iter
    else:
        max_iter = 1000
    # Get current working directory
    if args.save_folder:
        parent_folder = args.save_folder
    else:
        parent_folder = os.getcwd()

    # Read in the settings and data
    print(f'max_iter = {max_iter}')
    print('Reading in data')
    settings_dict = read_detailed_eq_data_file(data_file)
    N = settings_dict['N']
    M = settings_dict['M']
    L = settings_dict['L']
    A_out_array = settings_dict['A_out_array']
    Ci_array = settings_dict['Ci_array']
    case = save_name.removeprefix(f'{N}-{M}-{L}_')

    # # Make a log file
    # log_file_name = f'{save_name}_MLP_reg.log'
    # log_file = open(log_file_name, 'w')
    # log_file.close()

    # Make a folder for results
    print('Making folder for results')
    model_type = 'MLP'
    folder = f'{case}_{model_type}'
    save_folder = f'{parent_folder}{sep}{folder}'
    if folder not in os.listdir(parent_folder):
        os.mkdir(save_folder)

    # # # Data preparation
    # Take a look at the raw data
    if args.raw_data:
        # plotting settings
        mpl.rcParams['font.size']=18
        title=f"D=B=T={N} {case} Raw Data"
        raw_data_plot = plot_raw_data(settings_dict, f"{N}-{M}-{L} {case} Raw Data", \
            0.92, save=False)
        plt.show()
    else:
        pass

    print('Partitioning data')
    # # Partition the data into subsets
    # Save the dataset indices
    full_set_size = np.shape(A_out_array)[0]
    train_set_size = int(0.8*full_set_size)
    dataset_size_list = [train_set_size, 3000, 2000, 1000, 500, 250, 100, 50, 25]
    data_set_df = pd.DataFrame(data=np.zeros((np.shape(A_out_array)[0], \
        len(dataset_size_list)+1), dtype='int'), \
    columns=["test_set"]+dataset_size_list)
    # Split data into 80/20 train/test
    D_train, T_train, D_final_test, T_final_test, train_array, test_list = \
    partition_data(A_out_array, Ci_array, train_set_size, data_set_df, \
        np.arange(full_set_size), test=True)
    D_train_list = [D_train]
    T_train_list = [T_train]
    # Sequentially partition the data sets so the smaller subsets don't have 
    # access to data larger don't have access to
    for i in range(len(dataset_size_list)-1):
        D_train, T_train, D_test, T_test, train_array, out_list = \
        partition_data(A_out_array, Ci_array, dataset_size_list[i+1], data_set_df, \
            train_array, test=False)
        D_train_list.append(D_train)
        T_train_list.append(T_train)
    data_set_df.to_csv(f'{parent_folder}{sep}{N}-{M}-{L}_{case}_dataset_partitions.csv')

    # # # MLP regression
    # # Hyperparameter optimization
    pipeline_MLP = Pipeline([('scaler', StandardScaler()), ('model', MLPRegressor(max_iter=max_iter))])
    parameters_MLP = {'model__activation': ['relu','tanh','logistic'], \
    # 2-2-2
    # 'model__hidden_layer_sizes': [(5,),(10,),(15,)], \
    # 10-10-10
    'model__hidden_layer_sizes': [(10,),(15,),(20,),(25,)], \
    'model__solver': ['lbfgs', 'adam'], \
    'model__alpha': [0.0001, 0.0005, 0.001, 0.005]}
    scoring_MLP = ['r2', 'neg_mean_absolute_error']
    MLP_metrics = [r2_score, mean_absolute_error]
    refit_MLP='r2'
    combo_number_MLP = 1
    for value_list in parameters_MLP.values():
        combo_number_MLP *= len(value_list)

    # Create dataframe to store results
    grid_search_list = []
    pred_train_list = []
    pred_test_list = []
    results_columns = ["dataset_size"] + [key.lstrip("model__") for key in parameters_MLP.keys()] + \
    ['mean_test_'+ score.lstrip("neg_") for score in scoring_MLP] + \
    ['std_test_'+ score.lstrip("neg_") for score in scoring_MLP] + \
    ['rank_test_'+ score.lstrip("neg_") for score in scoring_MLP] + ['refit_time'] + \
    ['train_'+ score.lstrip("neg_") for score in scoring_MLP] + ['train_pred_time'] + \
    ['test_'+ score.lstrip("neg_") for score in scoring_MLP] + ['test_pred_time']
    results_df = pd.DataFrame(columns=results_columns)
    results_df['dataset_size'] = dataset_size_list
    # alt_model_df = pd.DataFrame(columns=results_columns)

    # Grid search 5 fold cross validation
    print('Performing grid search with 5 fold cross validation')
    for i, dataset_size in enumerate(dataset_size_list):
        print(f'Dataset size: {dataset_size}')
        # with open(log_file_name, 'a') as log_file:
        #     with redirect_stdout(log_file):
        grid_search_MLP = GridSearchCV(pipeline_MLP, param_grid=parameters_MLP, \
            verbose=3, scoring=scoring_MLP, refit=refit_MLP)
        grid_search_MLP.fit(D_train_list[i], T_train_list[i])
        grid_search_list.append(grid_search_MLP)

        print('Saving grid search results')
        save_grid_search_results(results_df, i, grid_search_list[i], \
            dataset_size_list[i], refit_MLP, N, M, L, case, model_type, \
            scoring_MLP, save_folder=save_folder)

        # Predict with the best model
        print('Predicting with best fit model (R^2)')
        pred_train, pred_test, time_array = best_model_predict(D_train_list[i], \
            D_final_test, grid_search_list[i].best_estimator_)
        pred_train_list.append(pred_train)
        pred_test_list.append(pred_test)
        train_metrics = best_model_metrics(T_train_list[i], pred_train, MLP_metrics)
        test_metrics = best_model_metrics(T_final_test, pred_test, MLP_metrics)

        # Store the results
        add_predictions_to_df(results_df, i, train_metrics, time_array, scoring_MLP, \
            dataset_type='train')
        add_predictions_to_df(results_df, i, test_metrics, time_array, scoring_MLP, \
            dataset_type='test')

    # Save the results dataframe
    print('Saving results to csv')
    results_df.to_csv(f'{parent_folder}{sep}{N}-{M}-{L}_{case}_results_summary_{model_type}.csv')

    # # Look at the MAE rank 1 model
    # if results_df.loc[i,'rank_test_mean_absolute_error'] != 1:
    #     alt_model_params = get_best_metric_model(alt_model_df, i, \
    #         dataset_size_list, grid_search_list, \
    #         'rank_test_neg_mean_absolute_error', scoring_MLP)
    #     alt_model = MLPRegressor(activation=alt_model_params['model__activation'], \
    #                                alpha=alt_model_params['model__alpha'], \
    #                                hidden_layer_sizes=alt_model_params['model__hidden_layer_sizes'], \
    #                                solver=alt_model_params['model__solver'], \
    #                                max_iter=1000)
    #     # Save the alternate
    #     alt_model_save_file = f'{save_folder}/{N}-{M}-{L}_{case}_alt_{model_type}_{dataset_size_list[i]}.pkl'
    #     with open(alt_model_save_file,'wb') as model_file:
    #             pickle.dump(alt_model, model_file)
    #     alt_model_train_pred, alt_model_test_pred, alt_model_time_array, \
    #     alt_model_train_metrics, alt_model_test_metrics = \
    #     get_alt_model_results(alt_model_df, alt_model, i, D_train_list, \
    #         T_train_list, D_final_test, T_final_test, StandardScaler(), \
    #         MLP_metrics, scoring_MLP)

if __name__ == '__main__':
    mlp_reg_main()
    
        
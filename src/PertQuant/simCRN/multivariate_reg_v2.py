import numpy as np
import matplotlib.pyplot as plt
import math
import time
import pandas as pd
import pickle

def read_detailed_eq_data_file(file_name):
    with open(file_name, 'r') as file:
        # Read all the lines in the file
        lines = file.readlines()
        index = 0
        settings_dict = {}
        settings = True
        # Read in the settings
        while settings:
            line = lines[index]
            line = line.lstrip('# ').rstrip(' \n')
            if '=' in line:
                line_list = line.split(' = ')
                if line_list[0] == 'Measured':
                    settings_dict[line_list[0]] = bool(line_list[1])
                elif 'Cm' in line_list[0]:
                    settings_dict[line_list[0]] = float(line_list[1])
                else:
                    settings_dict[line_list[0]] = int(line_list[1])
                index += 1
            else:
                settings = False
        # Read in the header arrays
        header = True
        while header:
            line = lines[index]
            if '# ' in line:
                if 'i array' in line:
                    index += 1
                    array_name = line.lstrip('# ').rstrip('\n')
                    array = np.array(lines[index].lstrip('# ').rstrip('\n').split('\t'),dtype='float')
                    settings_dict[array_name] = array
                elif 'Association array' in line:
                    total_strands = settings_dict['N']+settings_dict['M']+settings_dict['L']
                    association_array = np.zeros((total_strands,total_strands))
                    for i in range(total_strands):
                        index += 1
                        association_array[i,:] = np.array(lines[index].lstrip('# ').rstrip('\n').split('\t'),dtype='float')
                    settings_dict['Association array'] = association_array
                index += 1
            else:
                header = False
            
        Ci_array = np.zeros((settings_dict['N_runs'],settings_dict['L']))
        A_out_array = np.zeros((settings_dict['N_runs'],settings_dict['N']))
        
        settings_dict["Ci_array"] = Ci_array
        settings_dict["A_out_array"] = A_out_array
        for i in np.arange(index,len(lines)):
            line_array = np.array(lines[i].lstrip('# ').rstrip('\n').split('\t'),dtype='float')
            Ci_array[i-index,:] = line_array[:settings_dict['L']]
            A_out_array[i-index,:] = line_array[settings_dict['L']:]
            
    return settings_dict

def plot_raw_data(settings_dict, title, y_title, save_file='', save=True):
    L = settings_dict['L']
    N = settings_dict['N']
    
    # Figure out the dimensions of the plot
    width = N*6.5
    height = L*5
    fig, ax = plt.subplots(L,N,figsize=(width, height))
    
    for i in range(L):
        for j in range(N):
            ax[i,j].scatter(settings_dict['A_out_array'][:,j], settings_dict['Ci_array'][:,i], marker='.', \
                           color='skyblue', alpha=0.5)
            ax[i,j].set_ylabel(f"$T_{i+1}$"+ " initial conc.")
            ax[i,j].set_xlabel(f"$D_{j+1}$"+ " measured conc.")
            ax[i,j].set_xlim(left=0,right=1.05)
            ax[i,j].set_ylim(bottom=0,top=2.1)
    fig.suptitle(title, y=y_title)
    if save and len(save_file)>0:
        fig.savefig(f"{save_file}.png", bbox_inches="tight")
    elif save and len(save_file)==0:
        fig.savefig(f"{title}.png", bbox_inches="tight")
    else:
        pass
    return (fig, ax)

def partition_data(D_data, T_data, train_size, data_set_df, index_array, test=False):
    # make arrays for the new data
    D_train_data = np.zeros((train_size, np.shape(D_data)[1]))
    T_train_data = np.zeros((train_size, np.shape(T_data)[1]))
    
    # Get random number generator
    rng = np.random.default_rng()
    # Keep track of which data is in which set
    out_list = list(index_array)    
    in_array = np.zeros(train_size, dtype='int')
    
    # Generate training set
    for i in range(train_size):
        index = rng.choice(out_list)
        in_array[i] = index
        D_train_data[i,:] = D_data[index,:]
        T_train_data[i,:] = T_data[index,:]
        out_list.remove(index)
        data_set_df.loc[index, train_size] = 1
    if test:
        D_test_data = np.zeros((np.shape(D_data)[0]-train_size, np.shape(D_data)[1]))
        T_test_data = np.zeros((np.shape(T_data)[0]-train_size, np.shape(T_data)[1]))
        for i, index in enumerate(out_list):
            D_test_data[i,:] = D_data[index,:]
            T_test_data[i,:] = T_data[index,:]
            data_set_df.loc[index, "test_set"] = 1
    else:
        D_test_data = np.array([])
        T_test_data = np.array([])
    
    return D_train_data, T_train_data, D_test_data, T_test_data, in_array, out_list

def get_partitioned_data(D_data, T_data, dataset_csv):
    columns = dataset_csv.columns.values
    data_dict = {}
    for i, column in enumerate(columns):
        try:
            column_key = int(column)
        except:
            column_key = column
        D_subset = D_data[dataset_csv[column]==1]
        T_subset = T_data[dataset_csv[column]==1]
        data_dict[column_key]=(D_subset, T_subset)
    return data_dict

def save_grid_search_results(results_df, index, grid_search, dataset_size, refit, \
case, model_type, scoring_list, save_folder=''):
    if len(save_folder) > 0:
        save_folder=f'{save_folder}/'
    # Save info about the best model to the dataframe
    print('Best 5 fold cross validated model parameters')
    print(grid_search.best_params_)
    print(f'Mean cross-validated {refit} score of best model: {grid_search.best_score_}')
    for param in grid_search.best_params_:
        column_name = param.split('model__')[1]
        results_df.loc[index,column_name]=grid_search.best_params_[param]
    
    # Save refitted model with best parameters
    with open(f'{save_folder}2-2-2_{case}_{model_type}_{dataset_size}.pkl','wb') as model_file:
        pickle.dump(grid_search.best_estimator_, model_file)
    print(f"Model refit time: {grid_search.refit_time_}")
    results_df.loc[index,'refit_time']= grid_search.refit_time_
    
    # Save grid search cross validation results
    cv_results_df = pd.DataFrame(grid_search.cv_results_)
    cv_results_df.to_csv(f'{save_folder}2-2-2_{case}_cv_results_{model_type}_{dataset_size}.csv')
    best_model = cv_results_df[cv_results_df[f'rank_test_{refit}']==1]
    
    grid_metric_list = ['mean','std','rank']
    # Add the grid search score metrics to the dataframe
    for score in scoring_list:
        for grid_metric in grid_metric_list:
            if 'neg' in score and 'mean' in grid_metric:
                value=-best_model[f'{grid_metric}_test_{score}'].values[0]
            else:
                value=best_model[f'{grid_metric}_test_{score}'].values[0]
            column_name= f'{grid_metric}_test_{score.lstrip('neg_')}'
            results_df.loc[index,column_name]=value
            print(f"{column_name}: {value}")

def best_model_predict(D_train, D_test, estimator):
    # Predict with the trained model and time it
    time_array = np.zeros(2)
    time_0 = time.perf_counter()
    pred_train = estimator.predict(D_train)
    time_1 = time.perf_counter()
    pred_test = estimator.predict(D_test)
    time_2 = time.perf_counter()
    time_array[0]=time_1-time_0
    time_array[1]=time_2-time_1
    return pred_train, pred_test, time_array

def best_model_metrics(T_true, T_pred, metrics):
    metric_array = np.zeros(len(metrics))
    for i, metric in enumerate(metrics):
        metric_array[i] = metric(T_true, T_pred)
    return metric_array

def add_predictions_to_df(results_df, index, metrics_results, time_array, \
    scoring_list, dataset_type='train'):
    if dataset_type=='train':
        pred_time = time_array[0]
    else:
        pred_time = time_array[1]
    results_df.loc[index, f'{dataset_type}_pred_time'] = pred_time
    print(f'Prediction time: {pred_time}')
    
    for i, score in enumerate(scoring_list):
        column_name= f'{dataset_type}_{score.lstrip('neg_')}'
        results_df.loc[index,column_name]=metrics_results[i]
        print(f"{column_name}: {metrics_results[i]}")

def train_model(D_train, T_train, scaler, model):
    # Train scaler on training data
    data_scaler = scaler.fit(D_train)
    D_train_scaled = data_scaler.transform(D_train)
    # Fit the model with scaled training data
    time_0 = time.perf_counter()
    model.fit(D_train_scaled, T_train)
    time_1 = time.perf_counter()
    return data_scaler, time_1-time_0

def predict_model(D_data, fitted_scaler, fitted_model):
    # Scale data
    D_data_scaled = fitted_scaler.transform(D_data)
    time_0 = time.perf_counter()
    pred = fitted_model.predict(D_data_scaled)
    time_1 = time.perf_counter()
    return pred, time_1-time_0

def get_best_metric_model(alt_model_df, dataset_index, dataset_size_list, \
    grid_search_list, rank_metric, scoring_list):
    # Add row to alt model results dataframe
    alt_model_df.loc[dataset_index, 'dataset_size'] = dataset_size_list[dataset_index]
    # Pick out model with best MAE
    cv_results_df = pd.DataFrame(grid_search_list[dataset_index].cv_results_)
    alt_model_results = cv_results_df[cv_results_df["rank_test_neg_mean_absolute_error"]==1]
    alt_model_params = alt_model_results['params'].values[0]
    # Add parameters to alt model results dataframe
    for param in alt_model_params:
        column_name = param.split('model__')[1]
        alt_model_df.loc[dataset_index,column_name]=alt_model_params[param]
    # Add cross validated model results to alt model results dataframe
    grid_metric_list = ['mean','std','rank']
    # Add the grid search score metrics to the dataframe
    for score in scoring_list:
        for grid_metric in grid_metric_list:
            if 'neg' in score and 'mean' in grid_metric:
                value=-alt_model_results[f'{grid_metric}_test_{score}'].values[0]
            else:
                value=alt_model_results[f'{grid_metric}_test_{score}'].values[0]
            column_name= f'{grid_metric}_test_{score.lstrip('neg_')}'
            alt_model_df.loc[dataset_index, column_name]=value
    return alt_model_params

def get_alt_model_results(alt_model_df, alt_model, dataset_index, D_train, \
    T_train, D_test, T_test, scaler, metrics, scoring_list):
    # Train model
    data_scaler, refit_time = train_model(D_train, T_train, scaler, alt_model)
    alt_model_df.loc[dataset_index,'refit_time'] = refit_time
    # Predict with model
    train_pred, train_pred_time = predict_model(D_train, data_scaler, alt_model)
    test_pred, test_pred_time = predict_model(D_test, data_scaler, alt_model)
    time_array = np.array([train_pred_time, test_pred_time])
    train_metrics = best_model_metrics(T_train, train_pred, metrics)
    test_metrics = best_model_metrics(T_test, test_pred, metrics)
    # Save prediction results
    add_predictions_to_df(alt_model_df, dataset_index, train_metrics, time_array, \
        scoring_list, dataset_type='train')
    add_predictions_to_df(alt_model_df, dataset_index, test_metrics, time_array, \
        scoring_list, dataset_type='test')
    return train_pred, test_pred, time_array, train_metrics, test_metrics

def compare_alt_model(results_df, alt_model_df):
    for index in alt_model_df.index.values:
        for set_type in ['train', 'test']:
            alt_score = alt_model_df.loc[index,f'{set_type}_r2']
            main_score = results_df.loc[index,f'{set_type}_r2']
            alt_model_df.loc[index,f'delta_{set_type}_r2'] = main_score-alt_score
            if alt_score > main_score:
                alt_model_df.loc[index,f'better_{set_type}_r2'] = True
            else:
                alt_model_df.loc[index,f'better_{set_type}_r2'] = False
            alt_score = alt_model_df.loc[index,f'{set_type}_mean_absolute_error']
            main_score = results_df.loc[index,f'{set_type}_mean_absolute_error']
            alt_model_df.loc[index,f'delta_{set_type}_mean_absolute_error'] = main_score-alt_score
            if alt_score < main_score:
                alt_model_df.loc[index,f'better_{set_type}_mean_absolute_error'] = True
            else:
                alt_model_df.loc[index,f'better_{set_type}_mean_absolute_error'] = False

def plot_true_vs_pred(T_true, T_pred, Cmin, Cmax, title, save_file='', max_cols=4, \
    save=True, col_scaler=5, row_scaler=5, color='cadetblue', alpha=0.5):
    L = np.shape(T_true)[1]
    # Figure size of plot
    if L <= max_cols:
        cols = L
    else:
        cols = max_cols
    rows = math.ceil(L/max_cols)
    width = cols*col_scaler
    height = rows*row_scaler
    fig, ax = plt.subplots(rows,cols,figsize=(width, height),layout="constrained")
    ylim_array = np.zeros((L,2))

    # Handle 1 row case
    if rows < 2:
        for i in range(L):
            ax[i].scatter(T_true[:,i],T_pred[:,i], alpha=alpha, marker='.', color=color)
            ax[i].plot([Cmin, Cmax], [Cmin, Cmax], color='k', linestyle=(0,(10,5)))
            ax[i].set_xlabel(f'True $T_{i+1}$')
            ax[i].set_ylabel(f'Predicted $T_{i+1}$')
            ylim_array[i,:]=ax[i].get_ylim()
        for i in range(L):
            ax[i].set_ylim(bottom=min(ylim_array[:,0]),top=max(ylim_array[:,1]))
    else:
        for i in range(L):
            row = math.floor(i)
            col = i % max_cols
            ax[row,col].scatter(T_true[:,i],T_pred[:,i], alpha=alpha, marker='.', color=color)
            ax[row,col].plot([Cmin, Cmax], [Cmin, Cmax], color='k', linestyle=(0,(10,5)))
            ax[row,col].set_xlabel(f'True $T_{i+1}$')
            ax[row,col].set_ylabel(f'Predicted $T_{i+1}$')
            ylim_array[i,:]=ax[row,col].get_ylim()
        for i in range(L):
            row = math.floor(i)
            col = i % max_cols
            ax[row,col].set_ylim(bottom=min(ylim_array[:,0]),top=max(ylim_array[:,1]))
    fig.suptitle(title)
    if save and len(save_file)>0:
        fig.savefig(f"{save_file}.png", bbox_inches="tight")
    elif save and len(save_file)==0:
        fig.savefig(f"{title}.png", bbox_inches="tight")
    else:
        pass
    return (fig, ax)

def plot_residuals(T_true, T_pred, Cmin, Cmax, title, save_file, \
    max_cols=4, save=True, col_scaler=6.5, row_scaler=5, color='cadetblue'):
    L = np.shape(T_true)[1]
    # Figure size of plot
    if L <= max_cols:
        cols = L
    else:
        cols = max_cols
    rows = math.ceil(L/max_cols)
    width = cols*col_scaler
    height = rows*row_scaler
    fig, ax = plt.subplots(rows,cols,figsize=(width, height),layout='constrained')
    ylim_array = np.zeros((L,2))
    # Handle 1 row case
    if rows < 2:
        for i in range(L):
            ax[i].scatter(T_true[:,i],T_true[:,i]-T_pred[:,i], alpha=0.5, \
                marker='.', color=color)
            ax[i].plot([Cmin, Cmax], [0, 0], color='k', linestyle=(0,(10,5)))
            ax[i].set_xlabel(f'True $T_{i}$')
            ax[i].set_ylabel(f'Residuals (True-Predicted $T_{i}$)')
            ylim_array[i,:]=ax[i].get_ylim()
        for i in range(L):
            ax[i].set_ylim(bottom=min(ylim_array[:,0]),top=max(ylim_array[:,1]))
    else:
        for i in range(L):
            row = math.floor(i)
            col = i % max_cols
            ax[row,col].scatter(T_true[:,i],T_true[:,i]-T_pred[:,i], alpha=0.5, \
                marker='.', color=color)
            ax[row,col].plot([Cmin, Cmax], [0, 0], color='k', linestyle=(0,(10,5)))
            ax[row,col].set_xlabel(f'True $T_{i}$')
            ax[row,col].set_ylabel(f'Residuals (True-Predicted $T_{i}$)')
            ylim_array[i,:]=ax[row,col].get_ylim()
        for i in range(L):
            row = math.floor(i)
            col = i % max_cols
            ax[row,col].set_ylim(bottom=min(ylim_array[:,0]),top=max(ylim_array[:,1]))
    fig.suptitle(title)
    if save and len(save_file)>0:
        fig.savefig(f"{save_file}.png", bbox_inches="tight")
    elif save and len(save_file)==0:
        fig.savefig(f"{title}.png", bbox_inches="tight")
    else:
        pass
    return (fig, ax)

def plot_metric_summary(results_df, case, model_type, metric, save_file, width=8, \
    height=6, save=True, xscale='linear', train_color='cadetblue', \
    test_color='darkseagreen', cv_color='gray', cv_alpha=0.5):
    metric = metric.lstrip('neg_')
    metric_dictionary = {'r2': '$R^2$', 'mean_absolute_error': 'MAE'}

    dataset_size_array = results_df['dataset_size'].to_numpy()
    mean_cv_metric = results_df[f'mean_test_{metric}'].to_numpy()
    std_cv_metric = results_df[f'std_test_{metric}'].to_numpy()
    train_metric = results_df[f'train_{metric}'].to_numpy()
    test_metric = results_df[f'test_{metric}'].to_numpy()
        
    fig, ax = plt.subplots(figsize=(width, height))
    ax.errorbar(dataset_size_array, mean_cv_metric, yerr=std_cv_metric, label='5-fold CV Mean', \
                fmt='.', capsize=3, color=cv_color, zorder=1, alpha=cv_alpha)
    ax.scatter(dataset_size_array, train_metric, label='Train', marker= '.', color=train_color)
    ax.scatter(dataset_size_array, test_metric, label='Test', marker= '.', color=test_color)
    ax.legend(loc='best')
    if metric=='r2':
        ax.plot([0,max(dataset_size_array)], [1,1], color='k', linestyle=(0,(10,5)),zorder=0)
    elif metric=='mean_absolute_error':
        ax.plot([0,max(dataset_size_array)], [0,0], color='k', linestyle=(0,(10,5)),zorder=0)
    
    ax.set_xlabel("Train Dataset Size")
    
    ax.set_ylabel(metric_dictionary[metric])
    ax.set_title(f'{case}: {metric_dictionary[metric]} for {model_type}s and Varying Train Dataset Size')
    ax.set_xscale(xscale)
    
    if save and len(save_file)>0:
        fig.savefig(f"{save_file}.png", bbox_inches="tight")
    elif save and len(save_file)==0:
        fig.savefig(f"{title}.png", bbox_inches="tight")
    else:
        pass
    return (fig, ax)

def plot_metric_best_vs_alt(results_df, alt_model_df, case, model_type, metric, \
    save_file, width=8,height=6, save=True, xscale='linear', \
    train_color='cadetblue', test_color='darkseagreen', cv_color='gray'):
    metric = metric.lstrip('neg_')
    metric_dictionary = {'r2': '$R^2$', 'mean_absolute_error': 'MAE'}

    dataset_size_array = results_df['dataset_size'].to_numpy()
    train_metric = results_df[f'train_{metric}'].to_numpy()
    test_metric = results_df[f'test_{metric}'].to_numpy()
    alt_dataset_size_array = alt_model_df['dataset_size'].to_numpy()
    alt_train_metric = alt_model_df[f'train_{metric}'].to_numpy()
    alt_test_metric = alt_model_df[f'test_{metric}'].to_numpy()
        
    fig, ax = plt.subplots(figsize=(width, height))
    ax.scatter(dataset_size_array, train_metric, label='Train', marker= '.', color=train_color)
    ax.scatter(dataset_size_array, test_metric, label='Test', marker= '.', color=test_color)
    ax.scatter(alt_dataset_size_array, alt_train_metric, label='Alt Train', marker= '*', color=train_color)
    ax.scatter(alt_dataset_size_array, alt_test_metric, label='Alt Test', marker= '*', color=test_color)
    ax.legend(loc='best')
    if metric=='r2':
        ax.plot([0,max(dataset_size_array)], [1,1], color='k', linestyle=(0,(10,5)), zorder=0)
    elif metric=='mean_absolute_error':
        ax.plot([0,max(dataset_size_array)], [0,0], color='k', linestyle=(0,(10,5)), zorder=0)
    
    ax.set_xlabel("Train Dataset Size")
    
    ax.set_ylabel(metric_dictionary[metric])
    ax.set_title(f'{case}: {metric_dictionary[metric]} for {model_type}s and Varying Train Dataset Size')
    ax.set_xscale(xscale)
    
    if save and len(save_file)>0:
        fig.savefig(f"{save_file}.png", bbox_inches="tight")
    elif save and len(save_file)==0:
        fig.savefig(f"{title}.png", bbox_inches="tight")
    else:
        pass
    return (fig, ax)

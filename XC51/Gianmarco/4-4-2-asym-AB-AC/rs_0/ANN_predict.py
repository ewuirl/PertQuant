import operator
import os
import pandas as pd
import sklearn.preprocessing
import sklearn.utils
import numpy as np
from sklearn.ensemble import RandomForestRegressor
import sklearn.metrics as sklm
import sklearn.gaussian_process.kernels
import pickle
import json
import sys
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.kernel_ridge import KernelRidge
from hyperopt import hp, tpe, fmin, Trials
from functools import partial
from molSimplifyAD.retrain.nets import build_ANN, auc_callback, cal_auc, compile_model
from molSimplifyAD.retrain.model_optimization import train_model_hyperopt
import keras
from keras import backend as K
import tensorflow as tf

descriptors = ['CN_IP', 'CN_EA', 'NN_IP', 'NN_EA', 'CN_HOMO', 'CN_LUMO', 'NN_HOMO', 'NN_LUMO', 'lCN_C_mull', 'lCN_N_mull', 'lNN_N1_mull', 'lNN_N2_mull'] # l: ligand only 

def nonnegative_filter(predictions):
    for count, value in enumerate(predictions):
        if value < 0:
            predictions[count] = 0 # set to zero if prediction is negative
    return predictions

def normalize_data(df_train, df_test, fnames, lname, unit_trans=1, debug=False):
    _df_train = df_train.copy().dropna(subset=fnames + lname)
    _df_test = df_test.copy().dropna(subset=fnames + lname)
    X_train, X_test = _df_train[fnames].values, _df_test[fnames].values
    y_train, y_test = _df_train[lname].values * unit_trans, _df_test[lname].values * unit_trans
    if debug:
        print("training data reduced from %d -> %d because of nan." %
              (len(df_train), y_train.shape[0]))
        print("test data reduced from %d -> %d because of nan." %
              (len(df_test), y_test.shape[0]))
    x_scaler = sklearn.preprocessing.StandardScaler()
    x_scaler.fit(X_train)
    X_train = x_scaler.transform(X_train)
    X_test = x_scaler.transform(X_test)
    y_scaler = sklearn.preprocessing.StandardScaler()
    y_scaler.fit(y_train)
    y_train = y_scaler.transform(y_train)
    y_test = y_scaler.transform(y_test)
    for ii, ln in enumerate(lname):
        if debug:
            print(("mean in target: %s, train/test: %f/%f" %
                   (ln, np.mean(y_train[:, ii]), np.mean(y_test[:, ii]))))
    return X_train, X_test, y_train, y_test, x_scaler, y_scaler


def main():

    current_take = '29A' # change as needed

    df = pd.read_csv(f"all_data_take_{current_take}.csv")
    df_train = df[df["type"] == "train"]
    df_val = df[df["type"] == "val"]
    df_train_all = df_train.append(df_val)
    df_test = df[df["type"] == "test"]

    fnames = descriptors
    lname = ['Spectral_Integral'] # change as needed

    X_train, X_test, y_train, y_test, x_scaler, y_scaler = normalize_data(
        df_train_all, df_test, fnames, lname, unit_trans=1, debug=True)    

    # dependencies = {'precision':precision,'recall':recall,'f1':f1}
    model = keras.models.load_model('model_Ir_phosphors.h5', compile=False)

    predictions = y_scaler.inverse_transform(model.predict(X_test))
    predictions = nonnegative_filter(predictions)

    ### Will append the predictions to the dataframe df_test, and save to an Excel file

    df_predictions = pd.DataFrame(predictions, columns = ['Prediction']) # making the predictions into a pandas dataframe, in order to append

    df_predictions.reset_index(drop=True, inplace=True)
    df_test.reset_index(drop=True, inplace=True)
	
    df_predictions = pd.concat([df_predictions, df_test], axis=1) # combining test data with the predictions, column wise

    # Reordering the columns of the dataframe
    column_order = ['ID', 'Prediction', 'Em50_50', 'Lifetime', 'Sigma', 'Spectral_Integral', 'type'] + descriptors
    df_predictions = df_predictions[column_order]

    # Writing it all to a csv

    df_predictions.to_csv(f"predictions_take_{current_take}.csv", index = False)

    # Getting the average % error
    avg_percent_error = 0 
    # will use exact value as denominator
    for index, row in df_predictions.iterrows():
        avg_percent_error += abs(row['Prediction'] - row['Spectral_Integral'])/row['Spectral_Integral'] * 100 # change as needed

    avg_percent_error /= df_predictions.shape[0]

    print(f'avg_percent_error: {avg_percent_error}')

    ### Getting the mean absolute error, and the standard deviation of errors
    errors = []
    for index, row in df_predictions.iterrows():
        errors.append(row['Prediction'] - row['Spectral_Integral']) # change as needed

    import statistics
    abs_errors = [abs(item) for item in errors]
    MAE = statistics.mean(abs_errors)
    stdev = statistics.stdev(errors)

    print(f'The MAE is {MAE:.4f} and the standard deviation of errors is {stdev:.4f}') # :.4f gives four decimal places
     
if __name__ == "__main__":
    main()
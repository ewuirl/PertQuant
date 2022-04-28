from numpy.random import seed
import pandas as pd
import numpy as np
import pickle
import sys, os
from keras import backend as K
import tensorflow as tf
from scipy.spatial import distance_matrix
from sklearn.metrics import pairwise_distances
import json
import sklearn
tf.compat.v1.disable_eager_execution()

seed(1)

# Load in whatever model is of interest here
from keras.models import load_model
model = load_model('../model_Ir_phosphors.h5', compile=False)

# Standard function to normalize data
def normalize_data(df_train, df_test, fnames, lname, unit_trans=1, debug=False):
    _df_train = df_train.copy().dropna(subset=fnames+lname)
    _df_test = df_test.copy().dropna(subset=fnames+lname)
    X_train, X_test = _df_train[fnames].values, _df_test[fnames].values
    y_train, y_test = _df_train[lname].values, _df_test[lname].values
    if debug:
        print("training data reduced from %d -> %d because of nan." % (len(df_train), y_train.shape[0]))
        print("test data reduced from %d -> %d because of nan." % (len(df_test), y_test.shape[0]))
    x_scaler = sklearn.preprocessing.StandardScaler()
    x_scaler.fit(X_train)
    X_train = x_scaler.transform(X_train)
    X_test = x_scaler.transform(X_test)
    y_scaler = sklearn.preprocessing.StandardScaler()
    y_scaler.fit(y_train)
    y_train = y_scaler.transform(y_train)
    y_test = y_scaler.transform(y_test)
    return X_train, X_test, y_train, y_test, x_scaler, y_scaler

descriptors = ['CN_IP', 'CN_EA', 'NN_IP', 'NN_EA', 'CN_HOMO', 'CN_LUMO', 'NN_HOMO', 'NN_LUMO', 'lCN_C_mull', 'lCN_N_mull', 'lNN_N1_mull', 'lNN_N2_mull']

current_take = '29A' # change as needed

# Load in the data, keep relevant features
df = pd.read_csv(f"../all_data_take_{current_take}.csv")
df_train = df[df["type"] == "train"]
df_val = df[df["type"] == "val"]
df_test = df[df["type"] == "test"]
df_train_all = df_train.append(df_val)
features = descriptors
target = ["Spectral_Integral"] # change as needed

# Scale the data
X_train, X_test, y_train, y_test, x_scaler, y_scaler = normalize_data(df_train_all, df_test, features, target, unit_trans=1, debug=False)

# Getting the layer before the last dense layer. This will be useful in next part of code. (can check for yourself using model.summary())
my_layers = model.layers
for i in reversed(range(len(my_layers))):
    if (type(my_layers[i]) is tf.python.keras.layers.core.Dense):
        index_last = i-1
        break

# Define the function for the latent space. This will depend on the model. We want the layer before the last (the one before dense-last).
get_latent = K.function([model.layers[0].input],
                        [model.layers[index_last].output]) 

# Get the latent vectors for the training data first, then the latent vectors for the test data.
training_latent = get_latent([X_train, 0])[0]
design_latent = get_latent([X_test, 0])[0]

# Compute the pairwise distances between the test latent vectors and the train latent vectors to get latent distances.
d1 = pairwise_distances(design_latent,training_latent,n_jobs=30)
df1 = pd.DataFrame(data=d1, index=df_test['ID'].tolist(), columns=df_train_all['ID'].tolist())
df1.to_csv(f'train_test_latent_dists_take_{current_take}.csv')

# # Get train train latent dists.
# d2 = pairwise_distances(training_latent,training_latent,n_jobs=30)
# df2 = pd.DataFrame(data=d2, index=df_train['ID'].tolist(), columns=df_train['ID'].tolist())
# df2.to_csv('train_train_latent_dists.csv')

# Next, will make a dataframe for which has the test complexes as the rows, and the the average latent distance to the 10 nearest train neighbors as the column.
# Will sort the column from greatest to least average distance (measure of uncertainty).

df_content = []

for index, row in df1.iterrows(): # iterate through the test complexes
    row_copy = row.copy()
    latent_distances = []

    # gathering the 10 nearest neighbors
    for i in range(10):
        min_column = row_copy.idxmin() # name of the train complex with the smallest latent distance to this test complex
        latent_distance = row_copy.at[min_column] # latent distance to this train complex
        row_copy = row_copy.drop(labels = min_column) # dropping the minimum distance train complex in order to get the next smallest in the subsequent for loop pass
        latent_distances.append(latent_distance)

    average_distance = np.mean(latent_distances)
    df_content.append([index, average_distance]) # add new row to the Pandas DataFrame with the desired index name (the test complex)

average_dist_df = pd.DataFrame(df_content, columns=['ID', 'avg_lat_dist_10']) # This dataframe will hold the average latent space distances
    # avg_lat_dist_10: average latent space distance across then 10 nearest neighbors in latent space

# Sorting average_dist_df now.
average_dist_df = average_dist_df.sort_values(by=['avg_lat_dist_10'])
average_dist_df.to_csv(f'train_test_avg_lat_dist_take_{current_take}.csv')
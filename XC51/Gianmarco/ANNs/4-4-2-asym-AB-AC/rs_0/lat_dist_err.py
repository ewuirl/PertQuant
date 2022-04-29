# Goal of this script: to make an Excel sheet with the columns ID, avg_lat_dist_10, and absolute error
# To be run after ann_predict_take_{TAKE}.py and UQ_take_{TAKE}.py

import pandas as pd
import numpy as np
from statistics import mean

current_take = '29A' # change as needed

avg_lat_df = pd.read_csv(f'train_test_avg_lat_dist_take_{current_take}.csv')
predictions_df = pd.read_csv(f'../predictions_take_{current_take}.csv')

df_content = []

for index, row in predictions_df.iterrows():
	current_ID = row['ID']
	current_abs_err = np.absolute(row['Prediction'] - row['Spectral_Integral']) # note: change the target property as needed
	
	my_row = avg_lat_df[avg_lat_df['ID'] == current_ID] # desired row of the avg_lat_dist_10 Pandas dataframe. Row corresponding to the current test complex 
	current_avg_lat_dist_10 = my_row['avg_lat_dist_10'].iloc[0]

	df_content.append([current_ID, current_avg_lat_dist_10, current_abs_err])

new_df = pd.DataFrame(df_content, columns=['ID', 'avg_lat_dist_10', 'abs_err'])

new_df.to_csv(f'latent_vs_error_take_{current_take}.csv')

## Next, check what the MAE is when using different avg_lat_dist_10 cutoffs
min_lat_dist = min(new_df['avg_lat_dist_10'])
max_lat_dist = max(new_df['avg_lat_dist_10'])
# print(f'min_lat_dist: {min_lat_dist}')
# print(f'max_lat_dist: {max_lat_dist}')

cutoff_info = [] 

spacings = np.linspace(min_lat_dist, max_lat_dist, 50) # 50 spacers
for cutoff in spacings:
	temp_df = new_df[new_df['avg_lat_dist_10'] <= cutoff]
	MAE = mean(temp_df['abs_err'])
	num_under = len(temp_df['abs_err']) # number of complexes within the UQ cutoff

	cutoff_info.append([cutoff, MAE, num_under])

cutoff_df = pd.DataFrame(cutoff_info, columns=['cutoff', 'MAE', 'number of complexes'])
cutoff_df.to_csv(f'cutoff_statistics_take_{current_take}.csv')

my_mean = mean(new_df["avg_lat_dist_10"])
my_std = np.std(new_df["avg_lat_dist_10"])
print(f'The average avg_lat_dist_10 (i.e. our UQ metric) is {my_mean}')
print(f'The standard deviation of the UQ metric is {my_std}')
print(f'The mean + 1 std is {my_mean + my_std}')
print(f'The mean + 2 std is {my_mean + 2 * my_std}')
print(f'The mean + 3 std is {my_mean + 3 * my_std}')
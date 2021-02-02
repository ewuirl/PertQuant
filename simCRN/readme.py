# Import functions
from multivariate_reg import read_eq_data_file
from multivariate_reg import convert_np2df
from multivariate_reg import get_stats
from multivariate_reg import Z_normalize_data
from multivariate_reg import prep_data
from multivariate_reg import plot_data
from multivariate_reg import plot_data
from multivariate_reg import plot_predict
from multivariate_reg import plot_true_and_pred
from multivariate_reg import plot_error_hist
from multivariate_reg import plot_true_v_error
from multivariate_reg import plot_residuals
from multivariate_reg import subset

# 1. Import the data, eg

Ci_all_array, Am_array, Cmin, Cmax, Ai = read_eq_data_file('filename.txt')

# 2. Prepare the data (converts to pandas data frame and randomly splits it into
# the data into testing and training arrays.)

Am_df, Ci_df, train_Am_df, train_Ci_df, test_Am_df, test_Ci_df = \
prep_data(Am_array, Ci_all_array, 0.8)

# 3. Normalize the data (I normalize input and output because they're both very small)

train_Am_stats = get_stats(train_Am_df)
normed_train_Am_df = Z_normalize_data(train_Am_df, train_Am_stats)
normed_test_Am_df = Z_normalize_data(test_Am_df, train_Am_stats)
train_Ci_stats = get_stats(train_Ci_df)
normed_train_Ci_df = Z_normalize_data(train_Ci_df, train_Ci_stats)
normed_test_Ci_df = Z_normalize_data(test_Ci_df, train_Ci_stats)

# Extra stuff ------------------------------------------------------------------
# Plotting functions
# 1. plot the data as input vs output variables
plot_data(normed_train_Am_df, normed_train_Ci_df, title='title')

# The rest are to be used with predictions and errors, so it'll throw an error
# since there's no model.

# Calculate test predictions and error
test_predictions = model.predict(normed_test_Am_df)
error = test_predictions - normed_test_Ci_df
scaled_mae = error.abs().mean()*train_Ci_stats["std"]

# 2. Plot the predicted vs true values
plot_predict(test_predictions, normed_test_Ci_df, lims=[-2,2])

# 3. Plot the predicted vs input values, and true vs input values
plot_true_and_pred(normed_test_Ci_df, test_predictions, normed_test_Am_df)

# 4. Plot an error histogram
plot_error_hist(error)

# 5. Plot true values vs the error
plot_true_v_error(normed_test_Ci_df, error)

# 6. Plot the residuals (input vs error)
green_list_all = [(0,0),(1,1),(2,2),(3,3)]
plot_residuals(normed_test_Am_df, error, green_list = green_list_all)


# Subset function
# This function takes a random subset of the data. Eg below it takes a subset of
# half the dataset and divides it into a training set and testing set
train_Am_df_1000, train_Ci_df_1000, test_Am_df_1000, test_Ci_df_1000, \
large_test_Am_df_1000, large_test_Ci_df_1000 = subset(Am_df, Ci_df, 0.5, 0.8)

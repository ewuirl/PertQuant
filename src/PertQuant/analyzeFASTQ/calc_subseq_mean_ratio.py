import numpy as np

# See A-B Ratio Discrepancy Subsequence Mean Ratios.ipynb

def calc_mean_subseq_count(count_array, x_array):
    min_len = x_array[0]
    max_len = x_array[-1]
    binned_counts = []
    mean_arr = np.zeros(int(max_len-min_len+1))
    std_dev_arr = np.zeros(int(max_len-min_len+1))
    current_len = min_len
    storage_list = []
    for i in range(len(x_array)):
        if current_len == x_array[i]:
            storage_list.append(count_array[i])
        else:
            binned_counts.append(storage_list)
            storage_list = [count_array[i]]
            current_len += 1
            if i == len(x_array)-1:
                binned_counts.append(storage_list)
            else:
                pass
    for i in range(len(binned_counts)):
        mean_arr[i] = np.mean(binned_counts[i])
        std_dev_arr[i] = np.std(binned_counts[i])
    return(mean_arr, std_dev_arr)

def add_means(mean_arr1, std_dev_arr1, mean_arr2, std_dev_arr2):
    summed_arr = mean_arr1 + mean_arr2
    summed_std_dev_arr = np.sqrt(std_dev_arr1 ** 2.0 + std_dev_arr2 ** 2.0)
    return (summed_arr, summed_std_dev_arr)

def ratio_means(mean_arr1, std_dev_arr1, mean_arr2, std_dev_arr2):
    ratio_arr = mean_arr1/mean_arr2
    ratio_std_dev_arr = ratio_arr * np.sqrt((std_dev_arr1/mean_arr1)** 2.0 + (std_dev_arr2/mean_arr2) ** 2.0)
    return (ratio_arr, ratio_std_dev_arr)

def mean_ratio(ratio_arr, ratio_std_dev_arr):
    ratio_mean = np.mean(ratio_arr)
    ratio_mean_std_dev = np.mean(ratio_std_dev_arr ** 2.0)
    return(ratio_mean, ratio_mean_std_dev)
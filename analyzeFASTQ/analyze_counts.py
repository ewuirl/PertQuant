import numpy as np
import scipy.optimize

def check_min_unique_len(seq_list, seq_len, min_len):
    # Iterate through the sequences and their complements
    for k in range(len(seq_list)):
        # Iterate through the different subsequence lengths
        for i in range(seq_len - min_len + 1):
            # Set the length of the subsequence to examine
            n = i + min_len
            # Iterate through the different subsequences
            for j in range(seq_len + 1 - n):
                # Iterate through the other sequences
                for l in range(len(seq_list)):
                    if l == k:
                        pass
                    else:
                        if seq_list[k][j:j + n] in seq_list[l]:
                            print(f'Found {k}\'s subseq of len {n} in {l}')
                        else:
                            pass

def read_counts(file_name, array):
    with open(file_name, "r") as file:
        index = 0
        while True:
            full_line = file.readline()
            if not full_line:
                break
            line = full_line.rstrip("\n")
            counts = line.split()
            for count in counts:
                array[index] = float(count)
                index += 1

def power_func(x, S, p):
    return S*p ** x

def autolabel(ax, rects, label_height,label_style='sci'):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        if label_style == 'sci':
            ax.annotate('{:1.2e}'.format(height), \
                        xy=(rect.get_x() + rect.get_width() / 2, label_height), \
                        xytext=(0, 3), \
                        textcoords="offset points", \
                        ha='center', va='bottom')
        else:
            ax.annotate('{:.2f}'.format(height), \
                        xy=(rect.get_x() + rect.get_width() / 2, label_height), \
                        xytext=(0, 3), \
                        textcoords="offset points", \
                        ha='center', va='bottom')

# def summed_power_func(x, S, p):
#     global seq_len
#     return S*(seq_len+1-x)*p ** x

def calc_sum(counts_array, sum_array, seq_len, min_len):
    index = 0
    for i in range(seq_len - min_len + 1):
        for j in range(seq_len - min_len - i + 1):
            sum_array[i] += counts_array[index]
            index += 1

def read_time_filt_count_files(folder_path, passfail, time_step, run_length, barcode_ID):
    name_base = f"{folder_path}/fastq_{passfail}/fastq_{passfail}"
    time_step_range = range(time_step, run_length*60+time_step, time_step)
    # Create lists to store count arrays in
    count_list = []
    comp_count_list = []
    for i in range(len(time_step_range)):
        time = time_step_range[i]
        file_name_base = f"_{time}_{barcode_ID}"
        # make arrays to store counts in 
        count_array = np.zeros(231)
        comp_count_array = np.zeros(231)
        # Read in the counts
        read_counts(f"{name_base}{file_name_base}_counts.txt",count_array)
        read_counts(f"{name_base}{file_name_base}_comp_counts.txt",comp_count_array)
        # Add the count arrays into the lists
        count_list.append(count_array)
        comp_count_list.append(comp_count_array)
    return(count_list, comp_count_list)

def sum_list(array_list, seq_len, min_len):
    summed_list = []
    for i in array_list:
        sum_array = np.zeros(seq_len-min_len+1)
        calc_sum(i,sum_array, seq_len, min_len)
        summed_list.append(sum_array)
    return(summed_list)

def time_summed_list(sum_list, seq_len, min_len):
    time_array = sum_list[0]
    time_summed_list = [time_array]
    for i in range(len(sum_list)-1):
        time_array = time_summed_list[i] + sum_list[i+1]
        time_summed_list.append(time_array)
    return(time_summed_list)

def fit_func(func, x_fit, array):
    fit_params, pcov = scipy.optimize.curve_fit(func, x_fit, array)
    fit = func(x_fit, *fit_params)
    perr = np.sqrt(np.diag(pcov))
    return(fit_params, fit, perr)

def fit_array_list(func, x_fit, array_list):
    fit_list = []
    val_arr = np.zeros((len(array_list),4))
    for i in range(len(array_list)):
        fit_params, fit, perr = fit_func(func, x_fit, array_list[i])
        fit_list.append(fit)
        val_arr[i,:] = [fit_params[0], fit_params[1], perr[0], perr[1]]
    return(fit_list, val_arr)

def sum_val_arr(val_arr1, val_arr2, summed=False):
    summed_val_arr = np.zeros((len(val_arr1), 3))
    for i in range(len(val_arr1)):
        summed_val_arr[i,0] = val_arr1[i,0] + val_arr2[i,0]
        if summed:
            summed_val_arr[i,1:] = np.sqrt(val_arr1[i,1:] ** 2.0 + val_arr2[i,1:] ** 2.0)
        else:
            summed_val_arr[i,1:] = np.sqrt(val_arr1[i,2:] ** 2.0 + val_arr2[i,2:] ** 2.0)
    return(summed_val_arr)

def ratio_val_arr(val_arr1, val_arr2, summed=False):
    ratio_arr = np.zeros((len(val_arr1), 2))
    for i in range(len(val_arr1)):
        ratio_arr[i,0] = val_arr1[i,0]/val_arr2[i,0]
        if summed:
            ratio_arr[i,1] = ratio_arr[i,0] * np.sqrt((val_arr1[i,1]/val_arr1[i,0])** 2.0 + \
                                                                 (val_arr2[i,1]/val_arr2[i,0])** 2.0 )
        else:
            ratio_arr[i,1] = ratio_arr[i,0] * np.sqrt((val_arr1[i,2]/val_arr1[i,0])** 2.0 + \
                                                                 (val_arr2[i,2]/val_arr2[i,0])** 2.0 )
    return(ratio_arr)
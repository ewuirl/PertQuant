import numpy as np
import scipy.optimize

def check_min_unique_len(seq_list, seq_len, min_len):
    """
    check_min_unique_len(seq_list, seq_len, min_len)

    This function takes in a list of target sequences and their complements, the
    length of the sequences (assuming the same length), and the minimum 
    subsequence length to check for. If a nonunique subsequence is found, a 
    statement is printed indicating the length of the subsequence, the ID of the
    target, and the index where the subsequence starts in the sequence.

    Arguments:
        seq_list (list): 
        seq_len (int): 
        min_len (int): 

    Returns:
        Nothing. Prints statements if nonunique subsequences are found.
    """
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

def read_binned_counts(count_file_path, array):
    with open(count_file_path, 'r') as file:
        index = 0
        while True:
            line = file.readline()
            if not line:
                break
            line = line.rstrip("\n")
            count_line_list = [int(i) for i in line.split()]
            array[index,:] = np.array(count_line_list)
            index += 1

def power_func(x, S, p):
    return S*p ** x

# def summed_power_func(x, S, p):
#     global seq_len
#     return S*(seq_len+1-x)*p ** x

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

def calc_sum(count_array, sum_array, seq_len, min_len):
    index = 0
    for i in range(seq_len - min_len + 1):
        for j in range(seq_len - min_len - i + 1):
            sum_array[i] += count_array[index]
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
    return(fit_params, fit, pcov, perr)

def sum_fit_func(count_arr, x_fit, seq_len, min_len, func):
    # Create an array to store the summed values in 
    summed_counts_arr = np.zeros(seq_len-min_len+1)
    # Sum the counts
    calc_sum(count_arr, summed_counts_arr, seq_len, min_len)
    # Fit the summed counts
    fit_params, fit, pcov, perr = fit_func(func, x_fit, \
        summed_counts_arr)
    return (summed_counts_arr, fit_params, fit, pcov, perr)

def sum_fit_arr(target_bin_array, x_fit, seq_len, min_len, bin_range, func):
    # Create lists to store data in
    bin_list = []
    summed_counts_list = []
    params_list = []
    fit_list = []
    pcov_list = []
    perr_list = []
    for i in range(len(target_bin_array)):
        count_arr = target_bin_array[i,:]
        # Check to see if counts for the bin are nonzero
        if np.sum(count_arr) > 0:
            # Fit the subsequence counts
            summed_counts_arr, fit_params, fit, pcov, perr = \
            sum_fit_func(count_arr, x_fit, seq_len, min_len, func)
            # Save results to lists
            bin_list.append(bin_range[i])
            summed_counts_list.append(summed_counts_arr)
            params_list.append(fit_params)
            fit_list.append(fit)
            pcov_list.append(pcov)
            perr_list.append(perr)
        # Skip if the subsequence counts are zero
        else:
            pass
    return(bin_list, summed_counts_list, params_list, fit_list, pcov_list, \
        perr_list)

def time_fit_arr(target_bin_array, x_fit, func, is_sum, seq_len, min_len, array_len):
    rows = len(target_bin_array)
    # Create lists to store data in
    if is_sum:
        time_summed_counts_arr = np.zeros((rows,len(x_fit)))        
    else:
        time_summed_counts_arr = np.zeros((rows,array_len))
    params_arr = np.zeros((rows,2))
    fit_list = []
    pcov_list = []
    perr_arr = np.zeros((rows,2))
    count_arr = np.zeros(array_len)

    for i in range(len(target_bin_array)):
        count_arr += target_bin_array[i,:]
        # Check to see if counts for the bin are nonzero
        if is_sum:
            # Fit the subsequence counts
            summed_counts_arr, fit_params, fit, pcov, perr = \
            sum_fit_func(count_arr, x_fit, seq_len, min_len, func)
            # Save summed results to lists
            time_summed_counts_arr[i,:] = summed_counts_arr
        else:
            # Save summed results to lists
            time_summed_counts_arr[i,:] = count_arr
            fit_params, fit, pcov, perr = fit_func(func, x_fit, count_arr)
        # Save fitresults to lists
        params_arr[i,:] = fit_params
        fit_list.append(fit)
        pcov_list.append(pcov)
        perr_arr[i,:] = perr
        
    return(time_summed_counts_arr, params_arr, fit_list, pcov_list, perr_arr)

def sum_binned_target_and_comp(target_params, target_perr, comp_params, comp_perr):
    # Create arrays to store data in 
    tot_S_arr = np.zeros(len(target_params))
    tot_S_err_arr = np.zeros(len(target_params))

    for i in range(len(target_params)):
        tot_S_arr[i] = target_params[i,0] + comp_params[i,0]
        tot_S_err_arr[i] = np.sqrt(target_perr[i,0] ** 2.0 + comp_perr[i,0] ** 2.0)

    return (tot_S_arr, tot_S_err_arr)

def divide_binned_counts(dividend_params, dividend_perr, divisor_params, divisor_perr):
    # Create arrays to store data in 
    ratio_arr = np.zeros(len(dividend_params))
    ratio_err_arr = np.zeros(len(dividend_params))

    for i in range(len(dividend_params)):
        ratio_arr[i] = dividend_params[i]/divisor_params[i]
        ratio_err_arr[i] = ratio_arr[i]*np.sqrt((dividend_perr[i]/dividend_params[i]) ** 2.0 \
            + (divisor_perr[i]/divisor_params[i]) ** 2.0)

    return (ratio_arr, ratio_err_arr)

def fit_array_list(func, x_fit, array_list):
    fit_list = []
    val_arr = np.zeros((len(array_list),4))
    for i in range(len(array_list)):
        fit_params, fit, perr = fit_func(func, x_fit, array_list[i])
        fit_list.append(fit)
        val_arr[i,:] = [fit_params[0], fit_params[1], perr[0], perr[1]]
    return(fit_list, val_arr)

def sum_binned_target_and_comp(target_params, target_perr, comp_params, comp_perr):
    # Create arrays to store data in 
    tot_S_arr = np.zeros(len(target_params))
    tot_S_err_arr = np.zeros(len(target_params))

    for i in range(len(target_params)):
        tot_S_arr[i] = target_params[i,0] + comp_params[i,0]
        tot_S_err_arr[i] = np.sqrt(target_perr[i,0] ** 2.0 + comp_perr[i,0] ** 2.0)

    return (tot_S_arr, tot_S_err_arr)

def divide_binned_counts(dividend_params, dividend_perr, divisor_params, divisor_perr):
    # Create arrays to store data in 
    ratio_arr = np.zeros(len(dividend_params))
    ratio_err_arr = np.zeros(len(dividend_params))

    for i in range(len(dividend_params)):
        ratio_arr[i] = dividend_params[i]/divisor_params[i]
        ratio_err_arr[i] = ratio_arr[i]*np.sqrt((dividend_perr[i]/dividend_params[i]) ** 2.0 \
            + (divisor_perr[i]/divisor_params[i]) ** 2.0)
    return (ratio_arr, ratio_err_arr)

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

# def plot_binned_fits(target, comp, bin_list, x_range, summed_counts_list, params_list, \
#     fit_list, perr_list, bin_arr, bin_step, bin_type, save_folder, sum=True, y_lims=(0,0), is_save=True):
#     for i in range(len(bin_list)):
#         bin_lower_bound = bin_list[i]
#         bin_upper_bound = bin_lower_bound + bin_step
#         fig, ax = plt.subplots()
#         if comp:
#             scatter_color = "k"
#         else:
#             scatter_color = "g"
#         ax.scatter(x_range, summed_counts_list[i], marker='.', color=scatter_color)
#         ax.plot(x_range, fit_list[i], label ="fit: S = {:.0f} $\pm$ {:.0f},\np = {:.2f} $\pm$ {:.1e}".\
#         format(params_list[i][0], perr_list[i][0], params_list[i][1], perr_list[i][1]))
#         ax.set_xticks(np.arange(min(x_fit), max(x_fit)+1, 5.0))
#         ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
#         if y_lims != (0,0):
#             ax.set_ylim(y_lims)
#         else:
#             pass
#         ax.legend(loc='best')
#         fig.suptitle(f'Fit of {target} Subsequence Counts, {bin_type} in [{bin_lower_bound:.2f},{bin_upper_bound:.2f}]')
#         ax.set_xlabel('Length of Subsequence')
#         ax.set_ylabel('Number of Matches')
#         if is_save:
#             if sum:
#                 png_save_name = f"{save_folder}/{target}_sum_fit_{bin_type}-{bin_lower_bound:.2f}".replace(".","-")
#             else:
#                 png_save_name = f"{save_folder}/{target}_fit_{bin_type}-{bin_lower_bound:.2f}".replace(".","-")
#             plt.savefig(png_save_name + ".png")
#         else:
#             pass


# def plot_time_binned_fits(target, comp, bin_list, x_range, time_summed_counts_arr, params_arr, \
#     fit_list, perr_arr, bin_arr, bin_step, bin_type, save_folder, sum=True, y_lims=(0,0), is_save=True):
#     for i in range(len(bin_list)):
#         bin_lower_bound = bin_list[i]
#         bin_upper_bound = bin_lower_bound + bin_step
#         fig, ax = plt.subplots()
#         if comp:
#             scatter_color = "k"
#         else:
#             scatter_color = "g"
#         ax.scatter(x_range, time_summed_counts_arr[i,:], marker='.', color=scatter_color)
#         ax.plot(x_range, fit_list[i], label ="fit: S = {:.0f} $\pm$ {:.0f},\np = {:.2f} $\pm$ {:.1e}".\
#         format(params_arr[i,0], perr_arr[i,0], params_arr[i,1], perr_arr[i,1]))
#         ax.set_xticks(np.arange(min(x_fit), max(x_fit)+1, 5.0))
#         ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
#         if y_lims != (0,0):
#             ax.set_ylim(y_lims)
#         else:
#             pass
#         ax.legend(loc='best')
#         fig.suptitle(f'Fit of {target} Subsequence Counts, {bin_type} in [{bin_lower_bound:.0f},{bin_upper_bound:.0f}]')
#         ax.set_xlabel('Length of Subsequence')
#         ax.set_ylabel('Number of Matches')
#         if is_save:
#             if sum:
#                 png_save_name = f"{save_folder}/{target}_sum_fit_{bin_type}-{bin_upper_bound:.0f}".replace(".","-")
#             else:
#                 png_save_name = f"{save_folder}/{target}_fit_{bin_type}-{bin_upper_bound:.0f}".replace(".","-")
#             plt.savefig(png_save_name + ".png")
#             plt.close()
#         else:
#             pass

def get_subsequence_counts(N, n, min_len, count_arr):
    """
    get_subsequence_counts(N, n, min_len, count_arr)

    This function takes in the total target length N, the target subsequence 
    length n, the minimum subsequence length, and a counts array and creates an 
    array containing all the subsequence counts for targets of length "n", with 
    each row representing the counts for each target subsequence. The counts are
    ordered from smallest subsequence length (min_len) to largest subsequence 
    length (n). For a particular subsequence length, the counts are ordered by 
    the position of start of the subsequence in the target subsequence 
    (eg 0, 1, ...).

    Arguments:
        N (int): 
        n (int): the size of the "target subsequences".
        min_len (int): the minimum subsequence length of the target that was 
            counted.
        count_arr (array): An array of subsequence counts for the target of 
            length N. The counts are ordered from smallest subsequence length 
            (min_len) to largest subsequence length (N). For a particular 
            subsequence length, the counts are ordered by the position of
            start of the subsequence in the target sequence (eg 0, 1, ...).

    Returns:

    """
    # assert (n<N), "The length of the target subsequences must be smaller than the target sequence length."
    assert (n>=min_len), "The length of the target subsequences must be greater than the minimum subsequence length."

    # Figure out the dimensions of the array to store the subsequence counts in
    rows = N-n+1
    cols = int((n-min_len+1)*(n-min_len+2)/2)
    # Make an array to store the counts
    subseq_counts_array = np.zeros((rows,cols))

    # Iterate through all the subsequences of length n
    for i in range(rows):
        # Set the starting point of the indices
        count_index = i
        row_index = -1
        # Iterate through the subsubsequence lengths for each subsequence
        for j in range(n-min_len+1):
            # Bump the index to account for the number of subsubsequences of the
            # original sequence
            if j > 0:
                count_index += N-min_len-j+2
            else:
                pass
            # Iterate through the number of subsubsequences for each subsequence
            for h in range(n-j-min_len+1):
                # Bump the index to account for which subsubsequence is being added
                row_index += 1
                # Add the subsubsequence count to the array
                subseq_counts_array[i, row_index] = count_arr[count_index + h]
    return(subseq_counts_array)

def bin_fit_arr(target_bin_array, x_fit, func, is_sum, seq_len, min_len):
    rows = len(target_bin_array)
    # Create lists to store data in
    if is_sum:
        all_summed_counts_arr = np.zeros((rows,len(x_fit)))        
    else:
        all_summed_counts_arr = np.array([])
    params_arr = np.zeros((rows,2))
    fit_list = []
    pcov_list = []
    perr_arr = np.zeros((rows,2))

    for i in range(len(target_bin_array)):
        count_arr = target_bin_array[i,:]
        # Check to see if counts for the bin are nonzero
        if is_sum:
            # Fit the subsequence counts
            summed_counts_arr, fit_params, fit, pcov, perr = \
            sum_fit_func(count_arr, x_fit, seq_len, min_len, func)
            # Save summed results to lists
            all_summed_counts_arr[i,:] = summed_counts_arr
        else:
            # Save fit results to lists
            fit_params, fit, pcov, perr = fit_func(func, x_fit, count_arr)
        # Save fit results to lists
        params_arr[i,:] = fit_params
        fit_list.append(fit)
        pcov_list.append(pcov)
        perr_arr[i,:] = perr
        
    return(all_summed_counts_arr, params_arr, fit_list, pcov_list, perr_arr)

def make_x_fit_arr(seq_len, min_len):
    # Figure out the maximum nunber of subsequences possible
    max_num_subseqs = seq_len - min_len + 1
    # Calculate the count array length
    array_len = int(max_num_subseqs * (max_num_subseqs+1)/2)
    # Create arrays for fitting
    x_array = np.zeros(array_len)
    # Fill in the x array
    index = 0
    for i in range(seq_len - min_len + 1):
        for j in range(seq_len - min_len - i + 1):
            x_array[index] = i + min_len
            index += 1
    x_fit = np.arange(min_len,seq_len+1)
    return(x_array, x_fit)

def fit_artificial_sequences(target, comp, total_counts_arr, start_len, end_len, N, min_len, \
    fit_save_folder, power_func, summed_power_func, plot=True):
    rows = int(((N-end_len+1)+(N-start_len+1))*(end_len-start_len+1)/2)
    all_params_arr = np.zeros((rows,2))
    all_perr_arr = np.zeros((rows,2))
    all_params_arr_s = np.zeros((rows,2))
    all_perr_arr_s = np.zeros((rows,2))
    index = 0
    for length in np.arange(start_len, end_len+1):
        # Make x range arrays
        subseq_counts = get_subsequence_counts(N, length, min_len, total_counts_arr)
        x_arr, x_fit = make_x_fit_arr(length, min_len)
        # Make the non summed fits
        empty_counts_arr, params_arr, fit_list, pcov_list, perr_arr = \
        bin_fit_arr(subseq_counts, x_arr, power_func, False, length, min_len)
        # Make the summed fits
        summed_counts_arr_s, params_arr_s, fit_list_s, pcov_list_s, perr_arr_s = \
        bin_fit_arr(subseq_counts, x_fit, summed_power_func, True, length, min_len)
        # Save the results
        for i in range(len(params_arr)): 
            all_params_arr[index,:] = params_arr[i,:]
            all_perr_arr[index,:] = perr_arr[i,:]
            all_params_arr_s[index,:] = params_arr_s[i,:]
            all_perr_arr_s[index,:] = perr_arr_s[i,:]
            # Increase the index
            index += 1
        print(f"saved the results of {length}")
        # Save the plots
        if plot:
            plot_binned_fits(length, target, comp, x_arr, subseq_counts, params_arr, \
                fit_list, perr_arr, fit_save_folder, sum=False, y_lims=(0,0), \
                is_save=True)
            plot_binned_fits(length, target, comp, x_fit, summed_counts_arr_s, params_arr_s, \
                fit_list_s, perr_arr_s, fit_save_folder, sum=True, y_lims=(0,0), \
                is_save=True)
        else:
            pass

        print(f"Finished with length {length}.")

    return(all_params_arr, all_perr_arr, all_params_arr_s, all_perr_arr_s)

def sum_by_len_artificial_seq_fits(target_fit_arr0, target_error_arr0, comp_fit_arr1, \
    comp_error_arr1, start_len, end_len, N, min_len):
    # Create arrays for storing the average fit values in.
    sum_param_arr = np.zeros(np.shape(target_fit_arr0))
    sum_err_arr = np.zeros(np.shape(target_error_arr0))

    # Make an index for storing values into the arrays
    save_index = 0
    target_index = 0
    comp_index = -1
    # Iterate through the fit values and average them
    for i in np.arange(start_len, end_len + 1):
        comp_index += N - i + 1
        # Add up the param values + squared values
        for j in range(N - i + 1):
            print(f"target:{target_fit_arr0[target_index]}, comp: {comp_fit_arr1[comp_index]}")
            # sum_param_arr[save_index] += target_fit_arr0[target_index,:] + \
            # comp_fit_arr1[comp_index,:]
            # sum_err_arr[save_index] += np.sqrt(target_error_arr0[target_index,:] ** 2.0 + \
            #     comp_error_arr1[comp_index,:] ** 2.0)
            # Increase the read index
            target_index += 1
            # Decrease the comp index
            comp_index -= 1
            # Increase the save index
            save_index += 1
        comp_index += N - i + 1
    return(sum_param_arr, sum_err_arr)
    
def avg_by_len_artificial_seq_fits(fit_arr, fit_error_arr, start_len, end_len, N, min_len):
    # Create arrays for storing the average fit values in.
    array_len = end_len - start_len + 1
    avg_param_arr = np.zeros(array_len)
    avg_err_arr = np.zeros(array_len)
    # Make an index for storing values into the arrays
    save_index = 0
    read_index = 0
    # Iterate through the fit values and average them
    for i in np.arange(start_len, end_len + 1):
        # Create variables to sum param values + squared errors
        sum_param = 0
        square_err = 0
        # Add up the param values + squared values
        for j in range(N - i + 1):
            sum_param += fit_arr[read_index,0]
            square_err += fit_error_arr[read_index,0] ** 2.0
            # Increase the read index
            read_index += 1

        # Take the average
        avg_param_arr[save_index] = sum_param / (N - i + 1)
        avg_err_arr[save_index] = np.sqrt(square_err) / (N - i + 1)

        # Increase the save index
        save_index += 1

    return(avg_param_arr, avg_err_arr)

def ratio_by_len_artificial_seq_fits(num_fit_arr, num_err_arr, den_fit_arr, \
    den_err_arr, start_len, end_len, N):
    # Make arrays to store values
    array_len = 0
    len_range = np.arange(start_len, end_len + 1)
    for i in len_range:
        array_len += (N - i + 1) ** 2.0
    array_len = int(array_len)
    ratio_arr = np.zeros((array_len, 2))
    ratio_err_arr = np.zeros((array_len, 2))
    avg_ratio_arr = np.zeros((len(len_range), 2))
    avg_ratio_err_arr = np.zeros((len(len_range), 2))

    save_index = 0
    read_index = 0
    # Iterate through the different strand lengths
    for i in len_range:
        # Iterate through the numerator fits for length i
        for j in range(N - i + 1):
            numerator = num_fit_arr[read_index+j,:]
            num_error = num_err_arr[read_index+j,:]

            # Iterate through the denominator fits for length i
            for k in range(N - i + 1):
                denominator = den_fit_arr[read_index+k,:]
                den_error = den_err_arr[read_index+k,:]

                # Calculate the ratio
                ratio = numerator/denominator
                # Calculate the ratio error
                ratio_err = ratio * np.sqrt((num_error/numerator)**2.0 + \
                    (den_error/denominator)**2.0)
                # Save the results
                ratio_arr[save_index] = ratio
                ratio_err_arr[save_index] = ratio_err
                # Increase the save index
                save_index += 1
        # Increase the read index
        read_index += N - i + 1

    # Calculate the averages
    avg_save_index = 0
    avg_read_index = 0
    for i in len_range:
        tot_ratio = 0
        square_ratio_err = 0
        
        for j in range(int((N - i + 1) ** 2.0)):
            tot_ratio += ratio_arr[avg_read_index]
            square_ratio_err += ratio_err_arr[avg_read_index] ** 2.0
            # Increase the read index
            avg_read_index += 1

        avg_ratio_arr[avg_save_index] = tot_ratio / ((N - i + 1) ** 2.0)
        avg_ratio_err_arr[avg_save_index] = np.sqrt(square_ratio_err) / ((N - i + 1) ** 2.0)
        # Increase the save index
        avg_save_index += 1

    return(ratio_arr, ratio_err_arr, avg_ratio_arr, avg_ratio_err_arr)


def get_box_plot_data(ratio_arr, start_len, end_len, N):
    box_plot_list = []
    read_index = 0
    for i in np.arange(start_len, end_len + 1):
        array_index = 0
        len_ratio_arr = np.zeros(int((N-i+1)**2.0))
        for j in range(int((N-i+1)**2.0)):
            len_ratio_arr[array_index] = ratio_arr[read_index,0]
            # Increase the indices
            read_index += 1
            array_index += 1
        # Add the array to the list
        box_plot_list.append(len_ratio_arr)
    return(box_plot_list)




# # Test get_subsequence_counts
# count_array = np.arange(int((8-1)*(7+1)/2))
# sub_array = get_subsequence_counts(8,4,2,count_array)
# print(sub_array)

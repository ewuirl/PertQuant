import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt

def check_min_unique_len(seq_list, seq_len, min_len):
    """
    check_min_unique_len(seq_list, seq_len, min_len)

    This function takes in a list of target sequences and their complements, the
    length of the sequences (assuming the same length), and the minimum 
    subsequence length to check for. If a nonunique subsequence is found, a 
    statement is printed indicating the length of the subsequence, the ID of the
    target, and the index where the subsequence starts in the sequence.

    Arguments:
        seq_list (list): A list of target sequences and their complements. 
        seq_len (int): The length of the target sequences
        min_len (int): The minimum subsequence length to check for

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
    """
    read_counts(file_name, array)

    This function takes in the name of a subsequence count file (txt file) and 
    reads in the subsequence counts into an array of the appropriate size. All
    of the subsequence counts are stored on the same line.

    Arguments:
        file_name (str): A string representing the path to the subsequence count
            file to read in.
        array (arr): An array to store the subsequence counts in.

    Returns:
        Nothing. Fills in the provided array.
    """
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
    """
    read_binned_counts(count_file_path, array)

    This function takes in the name of a binned subsequence count file (txt file)
    and reads in the subsequence counts into an array of the appropriate size.

    Arguments:
        count_file_path (str): A string representing the path to the subsequence
            count file to read in. Each line in the txt file represented the
            counts for a particular bin.
        array (arr): An array to store the subsequence counts in.

    Returns:
        Nothing. Fills in the provided array.
    """
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
    """
    power_func(x, S, p)

    This function is used to fit the expected number of a particular subsequence
    of length x based on the true total number of strands of the entire sequence
    and the probability of calling a base correctly.

    Arguments:
        x (int): The subsequence length
        S (float): The true total number strands of the entire sequence.
        p (float): The probability of calling a base correctly.

    Returns:
        E[# of x] (float): The expected number of instances of a particular
            subsequence of length x.
    """
    return S*p ** x

# def summed_power_func(x, S, p):
#     """
#     summed_power_func(x, S, p)(x, S, p)

#     This function is used to fit the expected number of subsequences of length x
#     based on the true total number of strands of the entire sequence and the 
#     probability of calling a base correctly. 

#     This function is commented out 
#     because it relies on a global argument, which does not work well when 
#     imported.

#     Arguments:
#         x (int): The subsequence length
#         S (float): The true total number strands of the entire sequence.
#         p (float): The probability of calling a base correctly.

#     Global Arguments:    
#     seq_len (int): The length of the sequence.

#     Returns:
#         E[# of x] (float): The expected number of instances of subsequences of
#             length x.
#     """
#     global seq_len
#     return S*(seq_len+1-x)*p ** x

def autolabel(ax, rects, label_height, label_style='sci'):
    """
    autolabel(ax, rects, label_height,label_style='sci')

    This function attaches a text label above each bar in *rects*, displaying 
    its height.

    Arguments:
        ax (matplotlib ax): 
        rects (barplot rect): 
        label_height (int):
        label_style (str): Defaults to 'sci'
    """
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
    """
    calc_sum(count_array, sum_array, seq_len, min_len)

    This function takes in a count array (not binned), sums the subsequence
    counts by subsequence length, and then saves these summed counts into a 
    provided sum array.

    Arguments:
        count_array (arr): An array with subsequence counts.
        sum_array (arr): An array to store the summed subsequence counts in.
        seq_len (int): The target sequence length.
        min_len (int): The minimum subsequence length.

    Returns:
        Nothing. Stores the summed subsequence counts in the provided sum_array.

    """
    index = 0
    for i in range(seq_len - min_len + 1):
        for j in range(seq_len - min_len - i + 1):
            sum_array[i] += count_array[index]
            index += 1

def read_time_filt_count_files(folder_path, passfail, time_step, run_length, barcode_ID):
    """
    read_time_filt_count_files(folder_path, passfail, time_step, run_length, barcode_ID)

    Deprecated.
    This function reads in time filtered count files.

    Arguments:
        folder_path (str): A path to folder containing fastq pass/fail folders.
        passfail (str): "pass" or "fail", depending on whether time filtered
            pass or fail counts are being read.
        time_step (int): The time step in minutes that counts were filtered with.
        run_length (int): The run length in hours.
        barcode_ID (int): The ID number of the target sequence.

    Returns:
        count_list (list): a list of subsequence count arrays binned by time step
        comp_count_list (list): a list of complementary subsequence count arrays
            binned by time step
    """
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
    """
    sum_list(array_list, seq_len, min_len)

    Deprecated.
    This function takes a list of binned subsequence arrays, the sequence length,
    and the minimum subsequence length. It sums the counts of each binned array 
    by subsequence length and returns a list of binned, summed arrays.
    
    Arguments:
        array_list (list): a list of binned subsequence count arrays
        seq_len (int): The target sequence length.
        min_len (int): The minimum subsequence length.

    Returns:
        summed_list (list): A list of binned arrays where counts with the same
            subsequence length are summed together.
    """
    summed_list = []
    for i in array_list:
        sum_array = np.zeros(seq_len-min_len+1)
        calc_sum(i,sum_array, seq_len, min_len)
        summed_list.append(sum_array)
    return(summed_list)

def time_summed_list(sum_list, seq_len, min_len):
    """
    time_summed_list(sum_list, seq_len, min_len)

    Deprecated.
    This function uses a list of time binned, summed arrays and sums them over
    time to create a list of arrays of the total subsequence counts over time.
    
    Arguments:
        sum_list (list): a list of binned, summed subsequence count arrays
        seq_len (int): The target sequence length.
        min_len (int): The minimum subsequence length.
    
    Returns:
        time_summed_list (list): A list of arrays of the counts of different
            subsequence lengths over time.
    """
    time_array = sum_list[0]
    time_summed_list = [time_array]
    for i in range(len(sum_list)-1):
        time_array = time_summed_list[i] + sum_list[i+1]
        time_summed_list.append(time_array)
    return(time_summed_list)

def fit_func(func, x_fit, array):
    """
    fit_func(func, x_fit, array)

    This function takes a function to fit and fits it to the provided data: an 
    array of input values and an array of output values. Fitting is done with
    scipy.optimize.curve_fit

    Arguments:
        func (function): The function to fit.
        x_fit (arr): An array of input data values, eg subsequence lengths.
        array (arr): An array of output data values, eg subsequence counts.

    Returns:
        fit_params (array): An array of the estimated parameter values.
        fit (array): An array of the fitted function (using the fit parameters 
            and x_fit as the input data values).
        pcov (array): The covariance array.
        perr (array): An array of one standard deviation errors for the parameters.
    """
    fit_params, pcov = scipy.optimize.curve_fit(func, x_fit, array)
    fit = func(x_fit, *fit_params)
    perr = np.sqrt(np.diag(pcov))
    return(fit_params, fit, pcov, perr)

def sum_fit_func(count_arr, x_fit, seq_len, min_len, func):
    """
    sum_fit_func(count_arr, x_fit, seq_len, min_len, func)

    This function takes a function to fit and fits it to the provided data: an 
    array of input values and an array of output values summed by subsequence
    length. Fitting is done with scipy.optimize.curve_fit. 

    Arguments:
        count_arr (arr): An array of output data values, eg subsequence counts.
            NOT summed.
        x_fit (arr): An array of input data values, eg subsequence lengths.
        seq_len (int): The target sequence length.
        min_len (int): The minimum subsequence length.
        func (function): The function to fit, aka summed_power_func.

    Returns:
        summed_counts_arr: The subsequence counts summed by subsequence length.
        fit_params (array): An array of the estimated parameter values.
        fit (array): An array of the fitted function (using the fit parameters 
            and x_fit as the input data values).
        pcov (array): The covariance array.
        perr (array): An array of one standard deviation errors for the parameters.
    """
    # Create an array to store the summed values in 
    summed_counts_arr = np.zeros(seq_len-min_len+1)
    # Sum the counts
    calc_sum(count_arr, summed_counts_arr, seq_len, min_len)
    # Fit the summed counts
    fit_params, fit, pcov, perr = fit_func(func, x_fit, \
        summed_counts_arr)
    return (summed_counts_arr, fit_params, fit, pcov, perr)

def time_fit_arr(target_bin_array, x_fit, func, is_sum, seq_len, min_len, array_len):
    """
    This function takes an array of time binned subsequence counts. Each row is a 
    separate bin. It sums the subsequence count data over time and fits this 
    summed data. If is_sum is True, it also sums the data by subsequence length 
    and fits the summed data.

    Arguments:
        target_bin_array (array): An array of time binned subsequence counts. 
            Each row is a separate bin. 
        x_fit (arr): An array of input data values, eg subsequence lengths.
        func (function): The function to fit, eg power_func or summed_power_func.
        is_sum (bool): If True, sums the data by subsequence length and fits the
            data using sum_fit_func. If False, fits the data using fit_func.
        seq_len (int): The target sequence length.
        min_len (int): The minimum subsequence length.
        array_len (int): The number of subsequence counts in each bin.
        func (function): The function to fit, aka summed_power_func.

    Returns:
        time_summed_counts_arr (arr): An array of the subsequence counts summed
            over time. Each row is a different time bin.
        params_arr (arr): A array of the fitted parameters for the 
            binned data. Each row is a separate bin.
        fit_list (list): A list of fitted function output arrays calculated 
            using the fitted parameters
        pcov_list (list): A list of the covariance arrays.
        perr_arr (arr): An array of one standard deviation errors for
            the fitted parameters. Each row is a separate bin.
    """
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

def sum_target_and_comp(target_params, target_perr, comp_params, comp_perr):
    """
    sum_target_and_comp(target_params, target_perr, comp_params, comp_perr)

    This function takes total sequence instance estimates (S) for a target 
    sequence and its complement and sums them together to get an estimate of the
    true total number of instances of that target. It also calculates the error 
    of this estimate using one standard deviation error estimates for the 
    target and complementary sequence values of S.

    Arguments:
        target_params (arr): An array of the fitted parameters for the target
            sequence.
        target_perr (arr): An array of the one standard deviation error estimates
            for the target sequence.
        comp_params (arr): An array of the fitted parameters for the complementary
            sequence.
        comp_perr (arr): An array of the one standard deviation error estimates
            for the complementary sequence.

    Returns:
        tot_S (float): An estimate the total number of instances of the target 
            sequence.
        tot_S_err (float): Estimated error of the total number of instances
            of the target sequence.
    """
    # Add target and complement S together
    tot_S = target_params[0] + comp_params[0]
    # Add the errors in quadrature
    tot_S_err = np.sqrt(target_perr[0] ** 2.0 + comp_perr[0] ** 2.0)

    return (tot_S, tot_S_err)

def divide_counts(dividend_params, dividend_perr, divisor_params, divisor_perr):
    """
    divide_counts(dividend_params, dividend_perr, divisor_params, divisor_perr)

    This function takes in total sequence instance estimates (S) and estimated
    one standard deviation errors for two sequences and computes the ratio of S
    for these sequences, and the ratio error.

    Arguments:
        dividend_params (float): The dividend value of S, the estimated total 
            instances of a sequence.
        dividend_perr (float):The estimated one standard deviation error of the 
            dividend value of S.
        divisor_params (float):The divisor value of S, the estimated total 
            instances of a sequence.
        divisor_perr (float):The estimated one standard deviation error of the 
            divisor value of S.
    Returns:
        ratio (float): The calculated ratio of the provided values of S.
        ratio_err (float): The estimated ratio error.
    """    
    ratio = dividend_params/divisor_params
    ratio_err = ratio*np.sqrt((dividend_perr/dividend_params) ** 2.0 \
        + (divisor_perr/divisor_params) ** 2.0)
    return (ratio, ratio_err)

def sum_binned_target_and_comp(target_params, target_perr, comp_params, comp_perr):
    """
    sum_binned_target_and_comp(target_params, target_perr, comp_params, comp_perr)

    This function takes in binned arrays of total sequence instance estimates (S)
    for a target sequence and its complement and sums them together to get an 
    estimate of the true total number of instances of that target. It also 
    calculates the error of this estimate using one standard deviation error 
    estimates for the target and complementary sequence values of S. Each row
    represents a separate bin.

    Arguments:
        target_params (arr): An array of the fitted parameters for the target
            sequence. Each row is a separate bin.
        target_perr (arr): An array of the one standard deviation error estimates
            for the target sequence. Each row is a separate bin.
        comp_params (arr): An array of the fitted parameters for the complementary
            sequence. Each row is a separate bin.
        comp_perr (arr): An array of the one standard deviation error estimates
            for the complementary sequence. Each row is a separate bin.

    Returns:
        tot_S_arr (arr): An array of the binned estimates of the total number of
            instances of the target sequence.
        tot_S_err_arr (arr): An array of the estimated errors of the total number
            of instances of the target sequence.
    """
    # Create arrays to store data in 
    tot_S_arr = np.zeros(len(target_params))
    tot_S_err_arr = np.zeros(len(target_params))

    for i in range(len(target_params)):
        tot_S_arr[i] = target_params[i,0] + comp_params[i,0]
        tot_S_err_arr[i] = np.sqrt(target_perr[i,0] ** 2.0 + comp_perr[i,0] ** 2.0)

    return (tot_S_arr, tot_S_err_arr)

def divide_binned_counts(dividend_params, dividend_perr, divisor_params, divisor_perr):
    """
    divide_binned_counts(dividend_params, dividend_perr, divisor_params, divisor_perr)

    This function takes in total sequence instance estimates (S) and estimated
    one standard deviation errors for two sequences and computes the ratio of S
    for these sequences, and the ratio error.

    Arguments:
        dividend_params (arr): An array of the fitted parameters for the 
            dividend sequence.
        dividend_perr (arr): An array of the one standard deviation error 
            estimates for the dividend sequence.
        divisor_params (arr):An array of the fitted parameters for the divisor 
            sequence.
        divisor_perr (arr):An array of the one standard deviation error estimates
            for the divisor sequence.
    Returns:
        ratio_arr (arr): An array containing the calculated ratio(s) of the 
            provided values of S.
        ratio_err_arr (arr): An array containing the estimated ratio error(s).
    """
    # Create arrays to store data in 
    ratio_arr = np.zeros(len(dividend_params))
    ratio_err_arr = np.zeros(len(dividend_params))

    for i in range(len(dividend_params)):
        ratio_arr[i] = dividend_params[i]/divisor_params[i]
        ratio_err_arr[i] = ratio_arr[i]*np.sqrt((dividend_perr[i]/dividend_params[i]) ** 2.0 \
            + (divisor_perr[i]/divisor_params[i]) ** 2.0)
    return (ratio_arr, ratio_err_arr)

def plot_time_binned_fits(target, comp, bin_arr, x_range, time_summed_counts_arr, params_arr, \
    fit_list, perr_arr, bin_step, save_folder, is_sum=True, y_lims=(0,0), is_save=True):
    """
    plot_time_binned_fits(target, comp, bin_arr, x_range, time_summed_counts_arr, params_arr, 
    fit_list, perr_arr, bin_step, save_folder, is_sum=True, y_lims=(0,0), is_save=True)

    This function plots the time binned count data together with the fitted function.
    Each time bin gets its own plot, and if is_save is True, plots are saved in 
    the specified save folder as png's. If sum is True, the names of the saved 
    files include "sum" in them. This function will overwrite saved plot files 
    if files with the same name already exist.

    Arguments:
        target (str): Name of the target, eg "Target A" or "Target A Complement". 
            This is used in the title of the plots.
        comp (bool): The markers of the scatter plot of the count data 
            are plotted in black if True, and green if False
        bin_arr (arr): The range of bin values. This should range from 0 to end
            run time. 
        x_range (arr): The input x array, eg x_array (not summed) or x_fit (for
            summed data).
        time_summed_counts_arr (arr): An array of the subsequence counts summed
            over time. Each row is a different time bin.
        params_arr (arr): A array of the fitted parameters for the 
            binned data. Each row is a separate bin.
        fit_list (list): A list of fitted function output arrays calculated 
            using the fitted parameters
        perr_arr (arr): An array of one standard deviation errors for
            the fitted parameters. Each row is a separate bin.
        bin_step (num): The time bin step size.
        save_folder (str): The path to the folder to save the plots in.

    Optional Arguments:
        sum (bool): If True, and is_save is True, edits the plot save file
            names to include "sum".
        y_lims (tuple): A tuple of the y axis limits (bottom, top). Defaults to
            (0,0), which sets the minimum y axis value to 0.
        is_save (bool): If True, saves the plots. Defaults to True.

    Returns:
        Nothing. Plots of the time binned count data and the fitted functions
        are made and optionally saved.

    """
    bin_type = "tstep"
    for i in range(len(bin_arr)):
        bin_lower_bound = bin_arr[i]
        bin_upper_bound = bin_lower_bound + bin_step
        fig, ax = plt.subplots()
        if comp:
            scatter_color = "k"
        else:
            scatter_color = "g"
        ax.scatter(x_range, time_summed_counts_arr[i,:], marker='.', color=scatter_color)
        ax.plot(x_range, fit_list[i], label ="fit: S = {:.0f} $\pm$ {:.0f},\np = {:.2f} $\pm$ {:.1e}".\
        format(params_arr[i,0], perr_arr[i,0], params_arr[i,1], perr_arr[i,1]))
        ax.set_xticks(np.arange(min(x_fit), max(x_fit)+1, 5.0))
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        if y_lims != (0,0):
            ax.set_ylim(y_lims)
        else:
            pass
        ax.legend(loc='best')
        fig.suptitle(f'Fit of {target} Subsequence Counts, {bin_type} in [{bin_lower_bound:.0f},{bin_upper_bound:.0f}]')
        ax.set_xlabel('Length of Subsequence')
        ax.set_ylabel('Number of Matches')
        if is_save:
            if is_sum:
                png_save_name = f"{save_folder}/{target}_sum_fit_{bin_type}-{bin_upper_bound:.0f}".replace(".","-")
            else:
                png_save_name = f"{save_folder}/{target}_fit_{bin_type}-{bin_upper_bound:.0f}".replace(".","-")
            plt.savefig(png_save_name + ".png", bbox_inches="tight")
            plt.close()
        else:
            pass

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
        N (int): The length of the full sequence.
        n (int): The length of the "target subsequences".
        min_len (int): the minimum subsequence length of the target that was 
            counted.
        count_arr (array): An array of subsequence counts for the target of 
            length N. The counts are ordered from smallest subsequence length 
            (min_len) to largest subsequence length (N). For a particular 
            subsequence length, the counts are ordered by the position of
            start of the subsequence in the target sequence (eg 0, 1, ...).

    Returns:
        subseq_counts_array (arr): An array of subsubsequence counts for all the
            different "target subsequences" of length n. Each row represents 
            the counts for a particular "target subsequence". 
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
    """
    bin_fit_arr(target_bin_array, x_fit, func, is_sum, seq_len, min_len)

    This function takes an array of binned subsequence counts and an array of
    input values and fits the provided function to this binned data. If is_sum
    is True, it sums the data by subsequence length and fits the summed data. 

    Arguments:
        target_bin_array (array): An array of binned subsequence counts. 
            Each row is a separate bin. 
        x_fit (arr): An array of input data values, eg subsequence lengths.
        func (function): The function to fit, eg power_func or summed_power_func.
        is_sum (bool): If True, sums the data by subsequence length and fits the
            data using sum_fit_func. If False, fits the data using fit_func.
        seq_len (int): The target sequence length.
        min_len (int): The minimum subsequence length.
        func (function): The function to fit, eg power_func or summed_power_func.

    Returns:
        all_summed_counts_arr (arr): If is_sum is True, this is an array of the 
            subsequence counts summed by subsequence length. Each row is a 
            different time bin. If is_sum is False, this is an empty array.
        params_arr (arr): A array of the fitted parameters for the 
            binned data. Each row is a separate bin.
        fit_list (list): A list of fitted function output arrays calculated 
            using the fitted parameters
        pcov_list (list): A list of the covariance arrays.
        perr_arr (arr): An array of one standard deviation errors for
            the fitted parameters. Each row is a separate bin.
    """
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
    """
    make_x_fit_arr(seq_len, min_len)

    This function makes input arrays to use for fitting with count data summed by
    subsequence length, and data not summed in this way.

    Arguments:
        seq_len (int): The target sequence length.
        min_len (int): The minimum subsequence length.

    Returns:
        x_array (arr): An array of input values to use with unsummed count data.
        x_fit (arr): An array of input values to use with count data summed by
            subsequence length.
    """
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

def plot_binned_artificial_fits(subseq_len, target, comp, x_range, binned_counts_arr, \
    params_arr, fit_list, perr_arr, save_folder, sum=True, y_lims=(0,0), is_save=True):
    """
    plot_binned_artificial_fits(subseq_len, target, comp, x_range, binned_counts_arr, 
    params_arr, fit_list, perr_arr, save_folder, sum=True, y_lims=(0,0), is_save=True)
    
    This function plots the count data of "target subsequences" of the provided
    subsequence length, together with the fitted function. Each subsequence gets
    its own plot, and if is_save is True, plots are saved in the specified save 
    folder as png's. If sum is True, the names of the saved files include "sum" 
    in them. This function will overwrite saved plot files if files with the same
    name already exist.

    Arguments:
        subseq_len (int): The length of the "target subsequences" being plotted.
        target (str): Name of the target, eg "Target A" or "Target A Complement". 
            This is used in the title of the plots.
        comp (bool): The markers of the scatter plot of the count data 
            are plotted in black if True, and green if False
        x_range (arr): The input x array, eg x_array (not summed) or x_fit (for
            summed data).
        binned_counts_arr (arr): An array of the subsequence counts for the 
            "target subsequences". Each row is a different subsequence.
        params_arr (arr): A array of the fitted parameters for the 
            binned data. Each row is a separate bin.
        fit_list (list): A list of fitted function output arrays calculated 
            using the fitted parameters
        perr_arr (arr): An array of one standard deviation errors for
            the fitted parameters. Each row is a separate bin.
        save_folder (str): The path to the folder to save the plots in.

    Optional Arguments:
        sum (bool): If True, and is_save is True, edits the plot save file
            names to include "sum".
        y_lims (tuple): A tuple of the y axis limits (bottom, top). Defaults to
            (0,0), which sets the minimum y axis value to 0.
        is_save (bool): If True, saves the plots. Defaults to True.

    Returns:
        Nothing. Plots of the target subsequence count data and the fitted functions
        are made and optionally saved.
    """
    for i in range(len(binned_counts_arr)):
        fig, ax = plt.subplots()
        if comp:
            scatter_color = "k"
        else:
            scatter_color = "g"
        ax.scatter(x_range, binned_counts_arr[i,:], marker='.', color=scatter_color)
        ax.plot(x_range, fit_list[i], label ="fit: S = {:.0f} $\pm$ {:.0f},\np = {:.2f} $\pm$ {:.1e}".\
        format(params_arr[i,0], perr_arr[i,0], params_arr[i,1], perr_arr[i,1]))
        ax.set_xticks(np.arange(min(x_fit), 25, 5.0))
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        if y_lims != (0,0):
            ax.set_ylim(y_lims)
        else:
            pass
        ax.legend(loc='best')
        fig.suptitle(f'Fit of Target {target} Len {subseq_len} v{i}')
        ax.set_xlabel('Length of Subsequence')
        ax.set_ylabel('Number of Matches')
#         plt.tight_layout()
        if is_save:
            if sum:
                png_save_name = f"{save_folder}/{target}_len_{subseq_len}_sum_fit_{i}".replace(".","-")
            else:
                png_save_name = f"{save_folder}/{target}_len_{subseq_len}_fit_{i}".replace(".","-")
            plt.savefig(png_save_name + ".png", bbox_inches="tight")
            plt.close()
        else:
            pass

def fit_artificial_sequences(target, comp, total_counts_arr, start_len, end_len, N, min_len, \
    fit_save_folder, power_func, summed_power_func, plot=True):
    """
    fit_artificial_sequences(target, comp, total_counts_arr, start_len, end_len, N, min_len, 
    fit_save_folder, power_func, summed_power_func, plot=True)

    This function takes start (smallest) and end lengths for "target subsequences",
    picks out the subsequence count data for these subsequences lengths, fits
    the provided functions power_func and summed_power_func, and optionally
    plots the data and fits.

    Arguments:
        target (str): Name of the target, eg "Target A" or "Target A Complement". 
            This is used in the title of the plots.
        comp (bool): The markers of the scatter plot of the count data 
            are plotted in black if True, and green if False 
        total_counts_arr (arr): An array of the subsequence counts.
        start_len (int): The minimum "target subsequence" length to analyze.
        end_len (int): The maximum "target subsequence" length to analyze.
        N (int): The length of the sequence.
        min_len (int): The minimum subsequence length.
        fit_save_folder (str): The path to the folder to save the plots in.
        power_func (func): The power function to fit.
        summed_power_func (func): The summed power function to fit.

    Optional Arguments:
        plot (bool): If True, saves plots of the "target subsequence" count data
            and the fitted functions. Defaults to True.

    Returns:
        all_params_arr (arr): An array containing the power_func fitted 
            parameters for all of the "target subsequences". Each row represents
            results for a different "target subsequence".
        all_perr_arr (arr): An array containing the estimated error for the 
            power_func fitted parameters for all of the "target subsequences". 
            Each row represents results for a different "target subsequence".
        all_params_arr_s (arr): An array containing the summed_power_func fitted 
            parameters for all of the "target subsequences". Each row represents
            results for a different "target subsequence".
        all_perr_arr_s (arr): An array containing the estimated error for the 
            summed_power_func fitted parameters for all of the "target 
            subsequences". Each row represents results for a different "target 
            subsequence".
    """
    assert (start_len>=min_len), "The length of the target subsequences must be greater than the minimum subsequence length."
    assert (N>=end_len), "The length of the target subsequences must be smaller than the total sequence length."
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
            plot_binned_artificial_fits(length, target, comp, x_arr, subseq_counts, \
                params_arr, fit_list, perr_arr, fit_save_folder, sum=False, \
                y_lims=(0,0), is_save=True)
            plot_binned_artificial_fits(length, target, comp, x_fit, \
                summed_counts_arr_s, params_arr_s, fit_list_s, perr_arr_s, \
                fit_save_folder, sum=True, y_lims=(0,0), is_save=True)
        else:
            pass

        print(f"Finished with length {length}.")

    return(all_params_arr, all_perr_arr, all_params_arr_s, all_perr_arr_s)

def sum_by_len_artificial_seq_fits(target_fit_arr, target_error_arr, comp_fit_arr, \
    comp_error_arr, start_len, end_len, N):
    """
    sum_by_len_artificial_seq_fits(target_fit_arr, target_error_arr, comp_fit_arr, \
    comp_error_arr, start_len, end_len, N)

    This function takes arrays of all the "target subsequence" fitted parameters
    and estimated error, and sums the parameters and errors of the corresponding
    complementary subsequence from another set of provided arrays.

    Arguments:
        target_fit_arr (arr): An array containing fitted parameters for all of 
            the "target subsequences". Each row represents results for a 
            different "target subsequence".
        target_error_arr (arr): An array containing the estimated error for the 
            fitted parameters for all of the "target subsequences". Each row 
            represents results for a different "target subsequence".
        comp_fit_arr (arr): An array containing fitted parameters for the
            complements of all the "target subsequences". Each row represents 
            results for a different "target subsequence".
        comp_error_arr (arr): An array containing the estimated error for the 
            fitted parameters for the complements of all the "target 
            subsequences". Each row represents results for a different "target 
            subsequence".
        start_len (int): The minimum "target subsequence" length analyzed.
        end_len (int): The maximum "target subsequence" length analyzed.
        N (int): The length of the sequence.

    Returns:
        sum_param_arr (arr): An array containing summed target and complement
            fitted parameters for all of the "target subsequences". Each row 
            represents results for a different "target subsequence".
        sum_err_arr (arr): An array containing error estimates of the summed 
            target and complement fitted parameters for all of the "target 
            subsequences". Each row represents results for a different "target 
            subsequence".
    """
    # Create arrays for storing the summed fit values in.
    sum_param_arr = np.zeros(np.shape(target_fit_arr))
    sum_err_arr = np.zeros(np.shape(target_error_arr))

    # Make an index for storing values into the arrays
    save_index = 0
    target_index = 0
    comp_index = -1
    # Iterate through the fit values and average them
    for i in np.arange(start_len, end_len + 1):
        comp_index += N - i + 1
        # Add up the param values + squared values
        for j in range(N - i + 1):
            sum_param_arr[save_index] += target_fit_arr[target_index,:] + \
            comp_fit_arr[comp_index,:]
            sum_err_arr[save_index] += np.sqrt(target_error_arr[target_index,:] ** 2.0 + \
                comp_error_arr[comp_index,:] ** 2.0)
            # Increase the read index
            target_index += 1
            # Decrease the comp index
            comp_index -= 1
            # Increase the save index
            save_index += 1
        comp_index += N - i + 1
    return(sum_param_arr, sum_err_arr)
    
def avg_by_len_artificial_seq_fits(fit_arr, fit_error_arr, start_len, end_len, N):
    """
    avg_by_len_artificial_seq_fits(fit_arr, fit_error_arr, start_len, end_len, N)

    This function takes arrays of all the "target subsequence" fitted parameters
    and estimated error, and averages the parameters and errors by subseqeuence
    length.

    Arguments:
        fit_arr (arr): An array containing fitted parameters for all of 
            the "target subsequences". Each row represents results for a 
            different "target subsequence".
        fit_error_arr (arr): An array containing the estimated error for the 
            fitted parameters for all of the "target subsequences". Each row 
            represents results for a different "target subsequence".
        start_len (int): The minimum "target subsequence" length analyzed.
        end_len (int): The maximum "target subsequence" length analyzed.
        N (int): The length of the sequence.

    Returns:
        avg_param_arr (arr): An array containing fitted parameters for the 
            "target subsequences", averaged by subsequence length. Each row 
            represents results for a different target subsequence length.
        avg_err_arr (arr): An array containing error estimates for the 
            "target subsequences", averaged by subsequence length. Each row 
            represents results for a different target subsequence length.
            subsequence".
    """
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
    """
    ratio_by_len_artificial_seq_fits(num_fit_arr, num_err_arr, den_fit_arr, \
    den_err_arr, start_len, end_len, N)

    This function takes subsequence arrays of "target subsequence" fitted parameters
    and estimated errors for two sequences, and calculates the ratio of the two
    sequences using parameters fitted for subsequences of the same length. It 
    also computes the errors of these ratios. It also computes the average ratio
    for each subsequence length.

    Arguments:
        num_fit_arr (arr): An array containing fitted parameters for all of 
            the "target subsequences" that will be used as the numerators in the
            ratios. Each row represents results for a different "target subsequence".
        num_err_arr (arr): An array containing the estimated error for the 
            fitted parameters for all of the "target subsequences" that will be 
            used as the numerators in the ratios. Each row represents results 
            for a different "target subsequence".
        den_fit_arr (arr): An array containing fitted parameters for all of 
            the "target subsequences" that will be used as the denominators in the
            ratios. Each row represents results for a different "target subsequence".
        den_err_arr (arr): An array containing the estimated error for the 
            fitted parameters for all of the "target subsequences" that will be 
            used as the denominators in the ratios. Each row represents results 
            for a different "target subsequence".
        start_len (int): The minimum "target subsequence" length analyzed.
        end_len (int): The maximum "target subsequence" length analyzed.
        N (int): The length of the sequence.

    Returns:
        ratio_arr (arr): An array containing ratios of the parameters for the 
            "target subsequences". Each row represents results for a different 
            target subsequence.
        ratio_err_arr (arr): An array containing ratio error estimates for the 
            "target subsequences". Each row represents results for a different 
            target subsequence.
        avg_ratio_arr (arr): An array containing ratios of the parameters for the 
            "target subsequences", averaged by subsequence length. Each row 
            represents results for a different target subsequence length.
        avg_ratio_err_arr (arr): An array containing ratio error estimates for 
            the "target subsequences", averaged by subsequence length. Each row 
            represents results for a different target subsequence length.
            subsequence".
    """
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
    """
    get_box_plot_data(ratio_arr, start_len, end_len, N)

    This function takes an array of "target subsequence" parameter ratios, pulls
    out the count ratios and puts them into arrays by subsequence length, and 
    adds these arrays to a list. This list of arrays can then be used to plot 
    box plots of the ratios by subsequence length.
    
    Arguments:
        ratio_arr (arr): An array containing ratios of the parameters for the 
            "target subsequences". Each row represents results for a different 
            target subsequence.
        start_len (int): The minimum "target subsequence" length analyzed.
        end_len (int): The maximum "target subsequence" length analyzed.
        N (int): The length of the sequence.

    Returns:
        box_plot_list (list): A list of arrays of "target subsequence" count 
            ratios, where each array contains ratios for a different subsequence
            length.
    """
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


def read_count_reads(count_read_file_name, time_step, run_length):
    """
    read_count_reads(count_read_file_name, time_step, run)

    This function takes in a file of read counts binned over time, and reads this
    data into an array.

    Arguments:
        count_read_file_name (str): The path to the file containing the read
            count data.
        time_step (int): The time step in minutes that reads were filtered with.
        run_length (int): The run length in hours.

    Returns:
        reads_array (arr): An array containing the total number of read counts
            over time.
    """
    time_range = np.arange(time_step,run_length*60+ time_step,time_step)
    reads_array = np.zeros(len(time_range))
    with open(count_read_file_name, 'r') as count_read_file:
        line = count_read_file.readline().rstrip('\n')
        read_counts = line.split()
    
    reads_array[0] = float(read_counts[0])
    for i in range(len(read_counts)-1):
        reads_array[i+1] = float(read_counts[i+1])+reads_array[i]
    return reads_array

def drop_subseq_counts(count_array, drop_list):
    
    ndims = np.ndim(count_array)
    if ndims == 2:
        orig_rows, orig_cols = np.shape(count_array)
    elif ndims == 1:
        orig_cols = len(count_array)
    else:
        print("count array has more than 2 dimensions")
    new_cols = orig_cols - len(drop_list)
    if ndims == 2:
        new_count_arr = np.zeros((orig_rows, new_cols))
    else:
        new_count_arr = np.zeros(new_cols)
    start = 0
    save = 0
    for drop_col in np.sort(drop_list):
        num_cols = drop_col - start
        if ndims == 2:
            new_count_arr[:,save:save+num_cols] = count_array[:,start:drop_col]
        else:
            new_count_arr[save:save+num_cols] = count_array[start:drop_col]
        start = drop_col + 1
        save = save + num_cols
    if start < orig_cols:
        if ndims == 2:
            new_count_arr[:,save:] = count_array[:,start:]
        else:
            new_count_arr[save:] = count_array[start:]
    return(new_count_arr)

def check_extended_sequence(subseq_len, extension, extended_seq, extended_seq_ID, target_list, sandwich):
    clear = True
    match_dict = {}
    front_check_list = []
    end_check_list = []
    if subseq_len <= len(extension):
        for i in range(subseq_len-1):
            front_check_list.append(extended_seq[len(extension)-subseq_len+i+1:len(extension)+i+1])
            end_check_list.append(extended_seq[-len(extension)-subseq_len+i+1:-len(extension)+i+1])
    else:
        for i in range(len(extension)):
            front_check_list.append(extended_seq[i:subseq_len+i])
            end_check_list.append(extended_seq[-subseq_len-i:len(extended_seq)-i])
    for i in range(len(front_check_list)):
        front_subseq = front_check_list[i]
        end_subseq = end_check_list[i]
        for j in range(len(target_list)):
            if sandwich == "sandwich" or sandwich == "front":
                front_loc = target_list[j].find(front_subseq)
                if front_loc >= 0:
                    print(f"Found subseq {front_subseq} in target {target_list[j]}")
                    match_dict.setdefault(front_subseq, []).append((subseq_len, "f", extended_seq_ID, j, front_loc))
                    clear = False
                else:
                    pass
            else:
                pass
            if sandwich == "sandwich" or sandwich == "end":
                end_loc = target_list[j].find(end_subseq)
                if end_loc >= 0:
                    print(f"Found subseq {end_subseq} in target {target_list[j]}")
                    match_dict.setdefault(end_subseq, []).append((subseq_len, "b", extended_seq_ID, j, end_loc))
                    clear = False
                else:
                    pass
            else:
                pass
    return clear, match_dict

def update_dict(dict1, dict2):
    for key in dict2.keys():
        val = dict2[key]
        dict1.setdefault(key, [])
        for entry in val:
            dict1[key].append(entry)
    return dict1

def add_dicts(dict_list):
    newdict = {}
    for dictionary in dict_list:
        for key in dictionary.keys():
            val = dictionary[key]
            newdict.setdefault(key, [])
            for entry in val:
                newdict[key].append(entry)
    return newdict

def extended_sequence_find_min(subseq_len, extensions, extended_seq, extended_seq_ID, target_list):
    success = 0
    check_len = subseq_len - 1
    match_dict = {}
    front_ext = extensions[0]
    end_ext = extensions[1]
    while success != 2 and check_len < len(target_list[0]):
        check_len += 1
        front_success, front_subseq_match_dict = check_extended_sequence(check_len, front_ext, extended_seq, \
                                                                         extended_seq_ID, target_list, "front")
        update_dict(match_dict, front_subseq_match_dict)
        end_success, end_subseq_match_dict = check_extended_sequence(check_len, end_ext, extended_seq, \
                                                                     extended_seq_ID, target_list, "end")
        update_dict(match_dict, end_subseq_match_dict)
        success = front_success + end_success
    return (check_len, success, match_dict)

def check_extended_sequence_list(subseq_len, extension, extended_seq_list, target_list, sandwich):
    success_arr = np.zeros(len(extended_seq_list))
    check_len = subseq_len - 1
    match_dict = {}
    while sum(success_arr) < 4 and check_len < len(target_list[0]):
        check_len += 1
        for j in range(len(extended_seq_list)):
            success_arr[j], subseq_match_dict = check_extended_sequence(check_len, extension, \
                                                                        extended_seq_list[j], j, target_list, \
                                                                        sandwich)
            match_dict.update(subseq_match_dict)
    return (check_len, success_arr, match_dict)

def dict_to_drop_lists(num_seq, match_dict):
    drop_lists = []
    for i in range(num_seq):
        drop_lists.append([])
    for key in match_dict.keys():
        val_list = match_dict[key]
        for match_tuple in val_list:
            seq_index = match_tuple[3]
            drop_index = match_tuple[4]
            if drop_index not in drop_lists[seq_index]:
                drop_lists[seq_index].append(drop_index)
            else:
                pass
    return drop_lists

def find_outliers(count_arr, sequence, seq_len, min_len):
    outlier_list = []
    arr_index = 0
    means = np.zeros(seq_len-min_len+1)
    std_devs = np.zeros(seq_len-min_len+1)
    for i in range(seq_len - min_len + 1):
        # Set the length of the subsequence to examine
        n = i + min_len
        # Iterate through the different subsequences
        subseq_count_arr = np.zeros(seq_len + 1 - n)
        for j in range(seq_len + 1 - n):
            subseq_count_arr[j] = count_arr[arr_index]
            arr_index += 1
        if n < seq_len:
            mean = np.mean(subseq_count_arr)
            std_dev = np.std(subseq_count_arr)
            means[i] = mean
            std_devs[i] = std_dev
            Z_arr = (subseq_count_arr - mean)/std_dev
            print(f"Subsequence len {n}, mean = {mean:.2f}, std_dev = {std_dev:.2f}")
            for j in range(seq_len + 1 - n):
                if abs(Z_arr[j]) >= 3:
                    print(f"{sequence[j:j+n]},\tpos = {j}:{j+n},\t{subseq_count_arr[j]},\tZ = {Z_arr[j]:.2f}")
                    outlier_list.append((arr_index-(seq_len + 1 - n)+j, Z_arr[j], sequence[j:j+n]))
                else:
                    pass
        else:
            means[-1] = subseq_count_arr[0]
            std_devs[-1] = 0.0
            print(f"Subsequence {n},\tmean = {subseq_count_arr[0]}")
    return(outlier_list, means, std_devs)

def calc_mean_subseq_count(count_array, x_array):
    """
    Drops nonunique subsequence counts
    """
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

# # Test get_subsequence_counts
# count_array = np.arange(int((8-1)*(7+1)/2))
# sub_array = get_subsequence_counts(8,4,2,count_array)
# print(sub_array)

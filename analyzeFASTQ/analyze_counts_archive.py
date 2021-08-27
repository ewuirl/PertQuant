def fit_array_list(func, x_fit, array_list):
    """
    fit_array_list(func, x_fit, array_list)

    Deprecated.
    This function takes a list of count subsequence arrays and fits the provided
    function to the data.

    Arguments:
        func (function): The function to fit.
        x_fit (arr): An array of input data values, eg subsequence lengths.
        array_list (list): A list of arrays containing subsequence count data.
    
    Returns:
        fit_list (list): A list of fitted function output arrays calculated 
            using the fitted parameters
        val_arr (arr): An array containing the fit parameters and their estimated
            one standard deviation errors. Each row consists of the following 
            values: [param_0, param_1, err_0, err_1]
    """
    fit_list = []
    val_arr = np.zeros((len(array_list),4))
    for i in range(len(array_list)):
        fit_params, fit, perr = fit_func(func, x_fit, array_list[i])
        fit_list.append(fit)
        val_arr[i,:] = [fit_params[0], fit_params[1], perr[0], perr[1]]
    return(fit_list, val_arr)

def sum_val_arr(val_arr1, val_arr2, summed=False):
    """
    Deprecated.
    """
    summed_val_arr = np.zeros((len(val_arr1), 3))
    for i in range(len(val_arr1)):
        summed_val_arr[i,0] = val_arr1[i,0] + val_arr2[i,0]
        if summed:
            summed_val_arr[i,1:] = np.sqrt(val_arr1[i,1:] ** 2.0 + val_arr2[i,1:] ** 2.0)
        else:
            summed_val_arr[i,1:] = np.sqrt(val_arr1[i,2:] ** 2.0 + val_arr2[i,2:] ** 2.0)
    return(summed_val_arr)

def ratio_val_arr(val_arr1, val_arr2, summed=False):
    """
    Deprecated.
    """
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

def sum_fit_arr(target_bin_array, x_fit, seq_len, min_len, bin_range, func):
    """
    sum_fit_arr(target_bin_array, x_fit, seq_len, min_len, bin_range, func)
    
    Deprecated.
    This function takes an array of binned subsequence counts. Each row is a 
    separate bin. For bins with nonzero subsequence counts, it sums the data in 
    each bin by subsequence length and fits the summed data.

    Arguments:
        target_bin_array (array):
        x_fit (arr): An array of input data values, eg subsequence lengths.
        seq_len (int): The target sequence length.
        min_len (int): The minimum subsequence length.
        bin_range (array):
        func (function): The function to fit, aka summed_power_func.

    Returns:
        bin_list (list): A list of the bins that had nonzero subsequence counts.
        summed_counts_list (list): A list of binned arrays of subsequence counts 
            summed by subsequence length.
        params_list (list): A list of arrays of the fitted parameters for the 
            binned data.
        fit_list (list): A list of fitted function output arrays calculated 
            using the fitted parameters
        pcov_list (list): A list of the covariance arrays.
        perr_list (list): A list of arrays of one standard deviation errors for
            the fitted parameters.
    """
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

def plot_binned_fits(target, comp, bin_arr, x_range, summed_counts_list, params_list, \
    fit_list, perr_list, bin_step, bin_type, save_folder, sum=True, y_lims=(0,0), is_save=True):
    """
    Deprecated.
    """
    for i in range(len(bin_list)):
        bin_lower_bound = bin_list[i]
        bin_upper_bound = bin_lower_bound + bin_step
        fig, ax = plt.subplots()
        if comp:
            scatter_color = "k"
        else:
            scatter_color = "g"
        ax.scatter(x_range, summed_counts_list[i], marker='.', color=scatter_color)
        ax.plot(x_range, fit_list[i], label ="fit: S = {:.0f} $\pm$ {:.0f},\np = {:.2f} $\pm$ {:.1e}".\
        format(params_list[i][0], perr_list[i][0], params_list[i][1], perr_list[i][1]))
        ax.set_xticks(np.arange(min(x_fit), max(x_fit)+1, 5.0))
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        if y_lims != (0,0):
            ax.set_ylim(y_lims)
        else:
            pass
        ax.legend(loc='best')
        fig.suptitle(f'Fit of {target} Subsequence Counts, {bin_type} in [{bin_lower_bound:.2f},{bin_upper_bound:.2f}]')
        ax.set_xlabel('Length of Subsequence')
        ax.set_ylabel('Number of Matches')
        if is_save:
            if sum:
                png_save_name = f"{save_folder}/{target}_sum_fit_{bin_type}-{bin_lower_bound:.2f}".replace(".","-")
            else:
                png_save_name = f"{save_folder}/{target}_fit_{bin_type}-{bin_lower_bound:.2f}".replace(".","-")
            plt.savefig(png_save_name + ".png", bbox_inches="tight")
        else:
            pass
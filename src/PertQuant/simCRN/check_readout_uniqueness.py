import numpy as np 
import time
import pandas as pd 

def get_pert_assoc_matrix(assoc_matrix, D, B):
    '''get_pert_assoc_matrix(assoc_matrix, D, B)

    This function takes in a symmetric association matrix for a CRN, the number
    of readout species D, the number of bead species B, and returns the 
    perturbation association matrix.

    Arguments:
        assoc_matrix (arr): A symmetric matrix containing all the association
            constants for the CRN. The ith row and column represent the same
            species. The readout species take the first D rows/columns, the 
            bead species take the next B rows/columns, and the target species
            take the next T rows/columns. 
        D (int): The number of readout species in the CRN.
        B (int): The number of bead species in the CRN.

    Returns: 
        pert_assoc_matrix (arr): An array containing the association constants 
            for the reactions between the targets and the readout/bead species 
            in the CRN. Each row represents a different target species. Each 
            column represents a different readout or bead species. The order is 
            set by the order of the species in the association matrix. For the 
            columns, all the readout species precede the bead species.
    '''
    return assoc_matrix[D+B:,:D+B]

def find_connected_graphs(assoc_matrix):
    '''find_connected_graphs(assoc_matrix)

    This function takes in a symmetric association matrix for a CRN and returns 
    a matrix that shows all the connected graphs in the network.

    Arguments:
        assoc_matrix (arr): A symmetric matrix containing all the association
            constants for the CRN. The ith row and column represent the same
            species.

    Returns: 
        connected_graphs (arr): An array containing the connected graphs in the
            CRN. Each row represents a distinct graph in the network. Each 
            column represents a different species. The order is set by the
            order of the species in the association matrix (D->B->T). The value
            corresponding to each species in a group indicates how many 
            reactions it takes part in. 
    '''

    if np.sum(assoc_matrix) > 0:
        # Make a copy of the association matrix
        mutable_matrix = assoc_matrix.copy()
        graphs = []
        rows, cols = np.shape(mutable_matrix)
        # Iterate through the rows
        for i in range(rows):
            # Check if the row has nonzero values
            if np.sum(mutable_matrix[i,:]) > 0:
                # Make a copy of the original association matrix row
                connected_graph = assoc_matrix[i,:].copy()
                connected_graph[i] += 1
                # Set the corresponding row in the mutable matrix to 0
                mutable_matrix[i,:]=np.zeros(cols)
                # Iterate through the columns
                for j in range(cols):
                    # Check if the column is nonzero
                    if connected_graph[j]>0:
                        # Add the rows that also have a nonzero value in this column
                        for l in range(rows):
                            if mutable_matrix[l,j] > 0:
                                connected_graph += mutable_matrix[l,:] != 0
                                # Get rid of the row when it's been added
                                mutable_matrix[l,:] = np.zeros(cols)
                graphs.append(connected_graph)
        return(np.vstack(graphs))
    else:
        return np.array([])

def check_readout_connectivity(assoc_matrix, D, B, T):
    '''check_readout_connectivity(assoc_matrix, D, B, T)

    This function takes in a symmetric association matrix for a CRN, the number
    of readout species D, and the number of bead species B. It returns True if 
    all the readout species are connected to at least one target species.

    Arguments:
        assoc_matrix (arr): A symmetric matrix containing all the association
            constants for the CRN. The ith row and column represent the same
            species.
        D (int): The number of readout species in the CRN.
        B (int): The number of bead species in the CRN.
        T (int): The number of target species in the CRN.

    Returns: 
        connected (bool): True if all the readout species are connected 
            to at least one target species. False if not.
    '''
    # Start with a quick check for direct connections
    direct_connection = assoc_matrix[:D,D+B:]
    summed_direct_connection = np.sum(direct_connection,axis=1)
    # Count how many readout strands are directly connected
    num_connected_D = len(summed_direct_connection[summed_direct_connection != 0])
    sufficient_connectivity = num_connected_D >= T
    if sufficient_connectivity:
        pass
    # Make connected graphs if the readouts aren't all directly linked 
    else:
        connected_readout = np.zeros(D)
        connected_graphs = find_connected_graphs(assoc_matrix)
        if np.shape(connected_graphs)[0] > 0:
            # Identify which graphs have a target in them
            targets = np.sum(connected_graphs[:,D+B:],axis=1)
            for i in range(len(targets)):
                if targets[i]>0:
                    connected_readout += connected_graphs[i,:D]
                else:
                    pass
            num_connected_D = np.sum(connected_readout.astype(bool))
            # Check that enough targets are connected
            sufficient_connectivity = num_connected_D >= T    
        else:
            pass
    return (num_connected_D, sufficient_connectivity)

def check_readout_connectivity_and_detectability(assoc_matrix, D, B, T):
    '''check_readout_connectivity_and_detectability(assoc_matrix, D, B, T)

    This function takes in a symmetric association matrix for a CRN, the number
    of readout species D, and the number of bead species B. It returns True if 
    all the readout species are connected to at least one target species.

    Arguments:
        assoc_matrix (arr): A symmetric matrix containing all the association
            constants for the CRN. The ith row and column represent the same
            species.
        D (int): The number of readout species in the CRN.
        B (int): The number of bead species in the CRN.
        T (int): The number of target species in the CRN.

    Returns: 
        connected (bool): True if all the readout species are connected 
            to at least one target species. False if not.
    '''

    connected_readout = np.zeros((2,D))
    connected_graphs = find_connected_graphs(assoc_matrix)
    connected_result_array = np.zeros(3)
    if np.shape(connected_graphs)[0] > 0:
        # Identify which graphs have a target in them
        targets = np.sum(connected_graphs[:,D+B:],axis=1)
        for i in range(len(targets)):
            if targets[i]>0:
                connected_readout[0,:] += connected_graphs[i,:D]
            else:
                pass
        # Check how many readout strands are connected to targets
        connected_result_array[0] = np.sum(connected_readout[0,:].astype(bool))

        # Identify which graphs have a bead in them
        beads = np.sum(connected_graphs[:,D:D+B],axis=1)
        for i in range(len(beads)):
            if beads[i]>0:
                connected_readout[1,:] += connected_graphs[i,:D]
            else:
                pass
        # Check how many readout strands are connected to beads
        connected_result_array[1] = np.sum(connected_readout[1,:].astype(bool))
        
        # Check how many readout strands are connected to targets and beads
        T_B_connected = np.logical_and(connected_readout[0,:],connected_readout[1,:])
        connected_result_array[2] = np.sum(T_B_connected)
    else:
        pass
    sufficiently_connected_array = connected_result_array >= T
    return (connected_result_array, sufficiently_connected_array)

def generate_random_matrix(min_K, max_K, D,B,T, rng):
    num_species = D+B+T
    # Generate array of random ints
    rand_int_array = rng.integers(min_K, high=max_K, \
        size=(num_species,num_species), endpoint=True)
    # Make the matrix symmetric
    return np.triu(rand_int_array)+np.transpose(np.triu(rand_int_array,k=1))

def make_base_mat_DB_one(D, B, T):
    n_mat = D+B+T
    base_mat = np.zeros((n_mat,n_mat))
    DB_mat = np.identity(D)
    base_mat[:D,D:D+B] = DB_mat
    base_mat[D:D+B,:D] = DB_mat
    return(base_mat)


def gen_stats_rand_DB(rng, D, B, T, base_mat, n_test, min_K, max_K_array, verbose=True):
    # Store  data
    stat_array = np.zeros((len(max_K_array),5),dtype='int')
    # Save the max K
    stat_array[:,0] = max_K_array
    n_mat = D+B+T
    for i, max_K in enumerate(max_K_array):
        if verbose:
            print(f"K = {max_K}, {i+1}/{len(max_K_array)}")
        else:
            pass
        for j in range(n_test):
            # Fix D-B as K = 1
            assoc_mat = base_mat.copy()
            # Generate D-T interactions
            rand_int_array = rng.integers(min_K, high=max_K, size=(D,T), endpoint=True)
            assoc_mat[0:D,D+B:] = rand_int_array
            assoc_mat[D+B:,0:D] = rand_int_array.transpose()
            # get perturbation association matrix
            pert_assoc_mat = get_pert_assoc_matrix(assoc_mat, D, B)
            # get and store the rank
            rank = np.linalg.matrix_rank(pert_assoc_mat)
            can_differentiate = rank==T
            stat_array[i,1] += can_differentiate
            # check the linkage
            connected_D, is_linked = check_readout_connectivity(assoc_mat, D, B, T)
            stat_array[i,2] += connected_D
            stat_array[i,3] += is_linked
            # uniqueness
            stat_array[i,4] += can_differentiate * is_linked
    return stat_array

def gen_stats_sparse_rand_DB(rng, D, B, T, n_test, sparsity_array, min_K, max_K, verbose=True, DB_array=[]):
    # Store  data
    if np.size(DB_array) > 0:
        stat_array = np.zeros((len(sparsity_array),5),dtype='float')
    else:
        stat_array = np.zeros((len(sparsity_array),9),dtype='float')
    # Save the varying variable
    stat_array[:,0] = sparsity_array
    n_mat = D+B+T
    for i, zero_rate in enumerate(sparsity_array):
        if verbose:
            print(f"K = {zero_rate}, {i+1}/{len(sparsity_array)}")
        else:
            pass
        for j in range(n_test):
            # Generate association constants
            rand_int_array = rng.integers(min_K, high=max_K, size=(n_mat,n_mat), endpoint=True)
            # Use zero rate to determine which association constants = 0
            zero_array = rng.choice([0,1],size=(n_mat,n_mat),p=[zero_rate,1-zero_rate])
            assoc_triu_mat = np.multiply(np.triu(rand_int_array),zero_array)
            # Use the D-B association constants, if provided
            if np.size(DB_array) > 0:
                assoc_triu_mat[0:D,D:D+B] = DB_array
            else:
                pass
            # Make the association matrix symmetric
            assoc_mat = assoc_triu_mat+np.triu(assoc_triu_mat,k=1).transpose()
            # get perturbation association matrix
            pert_assoc_mat = get_pert_assoc_matrix(assoc_mat, D, B)
            # get and store the rank
            rank = np.linalg.matrix_rank(pert_assoc_mat)
            can_differentiate = rank==T
            stat_array[i,1] += can_differentiate
            # check the linkage and uniqueness
            if np.size(DB_array) > 0:
                # check the linkage to T strands
                connected_D, is_linked = check_readout_connectivity(assoc_mat, D, B, T)
                stat_array[i,2] += connected_D
                stat_array[i,3] += is_linked
                # uniqueness
                stat_array[i,4] += can_differentiate * is_linked
            else:
                # check the linkage of readout to T and B strands
                connected_result_array, sufficiently_connected_array = \
                check_readout_connectivity_and_detectability(assoc_mat, D, B, T)
                stat_array[i,2:5] += connected_result_array
                stat_array[i,5:8] += sufficiently_connected_array
                # uniqueness
                stat_array[i,8] += can_differentiate * sufficiently_connected_array[2]            
    return stat_array

def vary_sparsity(rng, D, B, T, n_test, sparsity_array, min_K, max_K, \
    verbose=True, DB_array=[]):
    start=time.perf_counter()
    stat_array = gen_stats_sparse_rand_DB(rng, D, B, T, n_test, sparsity_array, \
        min_K, max_K, verbose=verbose, DB_array=DB_array)
    end=time.perf_counter()
    print("Time elapsed: {} second(s)".format(end-start))
    print("Percentage of identifiable arrays")
    print(stat_array[:,1]*100/n_test)
    if np.size(DB_array) > 0:
        print("Percent of linked arrays")
        print(stat_array[:,3]*100/n_test)
        print("Average number of linked readout strands")
        print(stat_array[:,2]/n_test)
        print("Percent of unique arrays")
        print(stat_array[:,4]*100/n_test)
    else:
        print("Percent of T-linked arrays")
        print(stat_array[:,5]*100/n_test)
        print("Average number of T-linked readout strands")
        print(stat_array[:,2]/n_test)
        print("Percent of B-linked arrays")
        print(stat_array[:,6]*100/n_test)
        print("Average number of T-linked readout strands")
        print(stat_array[:,3]/n_test)
        print("Percent of T- and B-linked arrays")
        print(stat_array[:,7]*100/n_test)
        print("Average number of T-linked readout strands")
        print(stat_array[:,4]/n_test)
        print("Percent of unique arrays")
        print(stat_array[:,8]*100/n_test)
    return stat_array

def generate_additional_sparsity(new_sparsity_array, old_sparsity_array, rng, \
    D, B, T, n_test, min_K, max_K, verbose=True, DB_array=[], rtol=1e-5):
    sparsity_list = []
    for zero_rate in new_sparsity_array:
        in_old_sparsity_array = False
        for old_zero_rate in old_sparsity_array:
            if np.isclose(zero_rate,old_zero_rate, rtol=rtol):
                in_old_sparsity_array = True
            else:
                pass
        if in_old_sparsity_array:
            pass
        else:
            sparsity_list.append(zero_rate)
    sparsity_array = np.array(sparsity_list)
    stat_array = vary_sparsity(rng, D, B, T, n_test, sparsity_array, min_K, max_K, \
                             verbose=verbose, DB_array=DB_array)
    return(stat_array)

def save_stat_array_text(array, file_name, column_names, write_mode="w"):
    with open(f"{file_name}.txt", write_mode) as save_file:
        if write_mode=="w":
            # Write column names
            save_file.write(" ".join(column_names))
        else:
            pass
        # Write data
        for i in range(np.shape(array)[0]):
            save_file.write("\n"+" ".join(array[i,:].astype(str)))

def append_additional_data(new_stat_array, old_stat_array, file_name, columns, save=False):
    combined_stat_array = np.concatenate((old_stat_array, new_stat_array))
    stat_df = pd.DataFrame(combined_stat_array,columns=columns)
    # Sort by varying parameter
    stat_df.sort_values(columns[0], inplace=True)
    stat_array = stat_df.to_numpy()
    if save:
        save_stat_array_text(stat_array, file_name, columns)
    else:
        pass
    return(stat_array)

def read_stat_array(file_name):
    with open(file_name, 'r') as read_file:
        lines = read_file.readlines()
        header = lines[0].rstrip("\n").split(" ")
        stat_array = np.zeros((len(lines)-1,len(header)))
        for i in range(len(lines)-1):
            line_list = lines[i+1].rstrip("\n").split(" ")
            stat_array[i,:] = np.array(line_list,dtype='float')
    return(stat_array)


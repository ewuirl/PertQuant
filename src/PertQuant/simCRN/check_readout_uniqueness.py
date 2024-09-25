import numpy as np 

def get_pert_assoc_matrix(assoc_matrix, D, T):
    '''get_pert_assoc_matrix(assoc_matrix, D, T)

    This function takes in a symmetric association matrix for a CRN, the number
    of readout species D, the number of target species T, and returns the 
    perturbation association matrix.

    Arguments:
        assoc_matrix (arr): A symmetric matrix containing all the association
            constants for the CRN. The ith row and column represent the same
            species.
        D (int): The number of readout species in the CRN.
        T (int): The number of target species in the CRN.

    Returns: 
        pert_assoc_matrix (arr): An array containing the association constants 
            for the reactions between the targets and the readout/bead species 
            in the CRN. Each row represents a different target species. Each 
            column represents a different readout or bead species. The order is 
            set by the order of the species in the association matrix. For the 
            columns, all the readout species precede the bead species.
    '''
    # get D-T association constants
    D_assoc = assoc_matrix[D:D+T,:D]
    # get T-B association constants
    B_assoc = assoc_matrix[D:D+T,D+T:]
    return np.concatenate((D_assoc,B_assoc),axis=1)

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
            order of the species in the association matrix.
    '''
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
                            connected_graph += mutable_matrix[l,:]
                            # Get rid of the row when it's been added
                            mutable_matrix[l,:] = np.zeros(cols)
            graphs.append(connected_graph)
    return(np.vstack(graphs))

def check_readout_connectivity(assoc_matrix, D, T):
    '''check_readout_connectivity(assoc_matrix, D, T)

    This function takes in a symmetric association matrix for a CRN, the number
    of readout species D, the number of target species T. It returns True if all
    the readout species are connected to at least one target species.

    Arguments:
        assoc_matrix (arr): A symmetric matrix containing all the association
            constants for the CRN. The ith row and column represent the same
            species.
        D (int): The number of readout species in the CRN.
        T (int): The number of target species in the CRN.

    Returns: 
        connected (bool): True if all the readout species are connected 
            to at least one target species. False if not.
    '''
    # Start with a quick check for direct connections
    direct_connection = assoc_matrix[:D,D:D+T]
    summed_direct_connection = np.sum(direct_connection,axis=1)
    if np.all(summed_direct_connection):
        return True
    # Make connected graphs if the readouts aren't all directly linked 
    else:
        connected_readout = np.zeros(D)
        connected_graphs = find_connected_graphs(assoc_matrix)
        # Identify which graphs have a target in them
        targets = np.sum(connected_graphs[:,D:D+T],axis=1)
        for i in range(len(targets)):
            if targets[i]>0:
                connected_readout += connected_graphs[:D,i]
            else:
                pass
        return np.all(connected_readout)

def generate_matrix(min_K, max_K, D,B,T):
    num_species = D+B+T
    # Generate array of random ints
    rand_int_array = rng.integers(min_K, high=max_K, \
        size=(num_species,num_species), endpoint=True)
    # Make the matrix symmetric
    return np.triu(rand_int_array)+np.transpose(np.triu(rand_int_array,k=1))
import numpy as np
import random as rand
from PertQuant.simCRN.ml_nupack import generate_strand

def gen_detector_strands(orthogonal_strand_file, save_file_name, n_strands, \
    len_addition):
    
    with open(orthogonal_strand_file, 'r') as read_file:
        with open(save_file_name, 'w') as save_file:
            lines = read_file.readlines()

            assert n_strands <= len(lines)

            for i in range(n_strands):
                strand = lines[i].rstrip('\n')
                addition = generate_strand(len_addition)
                save_file.write(strand+addition)
                if i < n_strands-1:
                    save_file.write('\n')


if __name__ == "__main__":
    # Initialize the random number generator
    rand.seed()

    orthogonal_strand_file = '20-base_orthogonal_strands.txt'

    n_strands = 20
    len_addition = 20
    save_file_name = f'{n_strands}_detector_strands.txt'

    gen_detector_strands(orthogonal_strand_file, save_file_name, n_strands, \
    len_addition)
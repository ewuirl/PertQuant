import numpy as np
import random as rand
from PertQuant.simCRN.ml_nupack import gen_complement

def gen_bead_strands(detector_strand_file, save_file_name, n_strands, \
    comp):
    
    with open(detector_strand_file, 'r') as read_file:
        with open(save_file_name, 'w') as save_file:
            lines = read_file.readlines()

            assert n_strands <= len(lines)

            for i in range(n_strands):
                detector_strand = lines[i].rstrip('\n')
                bead_strand = gen_complement(detector_strand, comp=comp, exact=False)
                save_file.write(bead_strand)
                if i < n_strands-1:
                    save_file.write('\n')

if __name__ == "__main__":
    # Initialize the random number generator
    rand.seed()

    detector_strand_file = '20_detector_strands.txt'
    n_strands = 20
    comp = 0.25
    comp_str = f'{comp}'.replace('.','-')
    save_file_name = f'{n_strands}_bead_strands_comp-{comp_str}.txt'

    gen_bead_strands(detector_strand_file, save_file_name, n_strands, \
    comp)
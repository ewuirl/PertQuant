import numpy as np
import random as rand
from PertQuant.simCRN.ml_nupack import generate_strand

if __name__ == "__main__":
    # Initialize the random number generator
    rand.seed()

    n_strands = 20
    length=40
    GC_content = 0.5
    save_file_name = f'{n_strands}_target_strands.txt'

    with open(save_file_name, 'w') as save_file:
        for i in range(n_strands):
            strand = generate_strand(length, GC_content=GC_content)
            save_file.write(strand)
            if i <= n_strands:
                save_file.write('\n')
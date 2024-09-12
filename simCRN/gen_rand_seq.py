import random as rand 
from ml_nupack import generate_strand


length = 20
num_strands = 100
strand_count = 0
strand_list = []

while strand_count < num_strands:
    strand = generate_strand(length, GC_content=rand.random())
    if strand not in strand_list:
        strand_list.append(strand)
        strand_count += 1
    else:
        pass

with open("100-strands.in", 'w') as file:
    for strand in strand_list:
        file.write(f"{strand}\n")
# enough = False
import numpy as np
from PertQuant.simCRN.ml_nupack import gen_complement
from PertQuant.simCRN.ml_nupack import generate_strand
from PertQuant.analyzeFASTQ.countmatches import check_list_min_unique_len
phos_seq = "/5Phos/"
# restriction_enzyme_seq = "TGCA"
sticky_end = "TGCA"
filler_seq = 'GTTCACCAGTTTGTCAAATCACCGCCACGGCGCTGCGCTTGCGTTGTATAACTAATA'

# filler_seq = generate_strand(57)
# print(filler_seq)
# isok = restriction_enzyme_seq in filler_seq
# print(isok)


# # Version 1 (long barcode)
# with open("bc25mer.240k.fasta", 'r') as file:
# 	for i in range(50):
# 		line = file.readline()
# 		if i % 2 == 0:
# 			pass
# 		else:
# 			seq = line.rstrip('\n')
# 			if restriction_enzyme_seq in seq:
# 				pass
# 			else:
# 				seq = seq + filler_seq + restriction_enzyme_seq
# 				print(seq)
# 				# comp_seq = gen_complement(seq)
# 				# print(comp_seq)

# # Version 2 (mega ligation)
# num_strands = 16
# with open("bc25mer.240k.fasta", 'r') as file:
# 	read_count = 0
# 	strand_count = 0
# 	while strand_count < num_strands:
# 		line = file.readline()
# 		read_count += 1
# 		if read_count % 2 == 1:
# 			pass
# 		else:
# 			seq = line.rstrip('\n')
# 			if restriction_enzyme_seq in seq:
# 				pass
# 			else:
# 				strand_count += 1
# 				comp_seq = phos_seq + gen_complement(seq)
# 				seq = phos_seq + restriction_enzyme_seq + seq + restriction_enzyme_seq
# 				print("Strand #{}".format(strand_count))
# 				print(seq)
# 				print(comp_seq)
# 	print("Read count {}".format(read_count))
# 	print("Strand count {}".format(strand_count))

# Version 3 (mega ligation + RE Digest experiments)
# Objectives:
# (1) 2 orthogonal strands that do not contain the SbfI-Hf RE recognition site:
# CCTGCAGG
# (2) terminal bases have A/T then G/C 
# (3) min unique sequence length = 5

restriction_enzyme_seq = "CCTGCAGG"
min_len = 5
num_strands = 50
num_pairs = 20
strand_list = []
comp_list = []
orthogonal_file = "/Users/emilywu/Dropbox (MIT)/Emily <-> Ashwin/PerturbativeQuantification/bc25mer.240k.fasta"
# # Find sequences that satisfy (1) and (2) 
# with open(orthogonal_file, 'r') as file:
# 	read_count = 0
# 	strand_count = 0
# 	while strand_count < num_strands:
# 		line = file.readline()
# 		read_count += 1
# 		if read_count % 2 == 1:
# 			pass
# 		else:
# 			seq = line.rstrip('\n')
# 			if restriction_enzyme_seq in seq:
# 				pass
# 			else:
# 				if seq[-1] == 'T' or seq[-1] == 'A' or \
# 				seq[0] == 'T' or seq[0] == 'A' or \
# 				seq[-2] == 'G' or seq[-2] == 'C' or \
# 				seq[1] == 'G' or seq[1] == 'C':
# 					pass
# 				else:
# 					strand_list.append(seq)
# 					comp = gen_complement(seq)
# 					comp_list.append(comp)
# 					strand_count += 1
# 	print("Read count {}".format(read_count))
# 	print("Strand count {}".format(strand_count))


# # # Figure out which strands have a unique sequence length of at least 5
# # success_list = []
# # pair_count = 0
# # for i in range(len(strand_list)-1):
# # 	for j in range(len(strand_list)-1-i):
# # 		if pair_count < num_pairs:
# # 			strand_i = strand_list[i]
# # 			comp_i = comp_list[i]
# # 			strand_j = strand_list[i+j+1]
# # 			comp_j = comp_list[i+j+1]
# # 			test_list = [strand_i, comp_i, strand_j, comp_j]
# # 			# Check that the ends match
# # 			if strand_i[:2] == strand_j[:2] and strand_i[-2:] == strand_j[-2:]:
# # 				match = 0
# # 			elif strand_i[:2] == comp_j[:2] and strand_i[-2:] == comp_j[-2:]:
# # 				match = 1
# # 			else:
# # 				match = -1
# # 			# If the ends match for the strands or strand + comp check for the min
# # 			# unique length
# # 			if match >= 0:
# # 				unique = check_list_min_unique_len(test_list, min_len)
# # 				if unique:
# # 					success_list.append((match,i,i+j+1))
# # 					pair_count += 1
# # 				else:
# # 					pass 
# # 		else:
# # 			break

# # if len(success_list) > 0:
# # 	for i in range(len(strand_list)):
# # 		print(f"{i}: {strand_list[i]}")
# # 	for success in success_list:
# # 		print(success)
# # 		if success[0] == 0:
# # 			print(strand_list[success[1]])
# # 			print(strand_list[success[2]])
# # 		elif success[0] == 1:
# # 			print(strand_list[success[1]])
# # 			print(comp_list[success[2]])
# # 		else:
# # 			print("I shouldn't have gotten here")
# # else:
# # 	print("No candidates found.")


strand_11 = "CAATTCATCCATATTGCACCGTGAG"
comp_11 = gen_complement(strand_11)
strand_39 = "CAGTGATTTTCGGCGGGCTCTAAAG"
comp_39 = gen_complement(strand_39)

print(phos_seq+sticky_end+strand_11+restriction_enzyme_seq+comp_39)
print(phos_seq+sticky_end+strand_39+restriction_enzyme_seq+comp_11)

# sticky_end = 'TGCA'
# with open("20-20-20-asym-AB-AC.in", 'r') as file:
# 	lines = file.readlines()
# 	for i in range(len(lines)):
# 		line = lines[i]
# 		if sticky_end in line[-4:]:
# 			print("found in {}".format(i))
# 		else:
# 			pass

from PertQuant.simCRN.ml_nupack import gen_complement
import time
import random as rand 
import concurrent.futures
import os
import datetime as dt
import numpy as np

rand.seed()

barcode = 'TATGAGGACGAATCTCCCGCTTATA'
target = 'ATGGTCGAATCAAGGGGAGG'
sticky_end = 'TGCA'
barcode_comp = gen_complement(barcode)
target_comp = gen_complement(target)
num_files = 5
num_lines = 1000

# Open the file
def generate_barcode_test_fastq(file_name):
	with open(file_name, 'w') as file:
		for i in range(num_lines):
			# Decide how many barcodes to add to sequence
			n = rand.randrange(2,21)
			# Decide whether to use sequence or complement
			comp = rand.random()

			# Write sequence ID
			file.write('@sequenceid{}\n'.format(i))

			# Write sequence
			if comp >= 0.3:
				for j in range(n):
					file.write(sticky_end)
					file.write(barcode)
				
			# Write sequence complement
			else:
				for j in range(n):
					file.write(sticky_end)
					file.write(barcode_comp)
			# Write the last sticky end
			file.write(sticky_end + '\n')

			# Write '+'
			file.write('+\n')

			# Write quality score
			if i < num_lines-1:
				file.write('qualityscore{}\n'.format(i))
			else:
				file.write('qualityscore{}'.format(i))

# Test serial generation
# # start timing
# start = time.perf_counter()

# for i in range(num_files):
# 	file_name = 'barcode_test{}.fastq'.format(i)
# 	generate_barcode_test_fastq(file_name)

# # end timing
# end = time.perf_counter()
# print('Serial Time elapsed: {} second(s)'.format(end-start))


# Test parallel process result
def generate_barcode_test_fastq_p(i):
	with open('barcode_test{}.fastq'.format(i), 'w') as file:
		for i in range(num_lines):
			# Decide how many barcodes to add to sequence
			n = rand.randrange(2,21)
			# Decide whether to use sequence or complement
			comp = rand.random()

			# Write sequence ID
			file.write('@sequenceid{}\n'.format(i))

			# Write sequence
			if comp >= 0.3:
				for j in range(n):
					file.write(sticky_end)
					file.write(barcode)
				
			# Write sequence complement
			else:
				for j in range(n):
					file.write(sticky_end)
					file.write(barcode_comp)
			# Write the last sticky end
			file.write(sticky_end + '\n')

			# Write '+'
			file.write('+\n')

			# Write quality score
			if i < num_lines-1:
				file.write('qualityscore{}\n'.format(i))
			else:
				file.write('qualityscore{}'.format(i))

# # Test multiprocessing
# # start timing
# start = time.perf_counter()

# with concurrent.futures.ProcessPoolExecutor() as executor:
# 	executor.map(generate_barcode_test_fastq_p, range(num_files))

# # end timing
# end = time.perf_counter()
# print('Multiprocess Time elapsed: {} second(s)'.format(end-start))

# # Test multithreading
# # start timing
# start = time.perf_counter()

# with concurrent.futures.ThreadPoolExecutor() as executor:
# 	executor.map(generate_test_fastq_p, range(num_files))

# # end timing
# end = time.perf_counter()
# print('Time elapsed: {} second(s)'.format(end-start))


# Open the file
def generate_test_fastq_p(i):
	with open('test{}.fastq'.format(i), 'w') as file:
		for i in range(num_lines):
			# Decide how many barcodes to add to sequence
			n = rand.randrange(2,21)
			# Decide whether to use sequence or complement
			comp = rand.random()

			# Write sequence ID
			file.write('@sequenceid{}\n'.format(i))

			# Write barcode
			if comp < 0.25:
				print('writing barcode')
				for j in range(n):
					file.write(sticky_end)
					file.write(barcode)
			
			# Write barcode complement
			elif 0.25 <= comp < 0.5:
				print('writing barcode_comp')
				for j in range(n):
					file.write(sticky_end)
					file.write(barcode_comp)

			# Write barcode and target complement
			elif 0.5 <= comp < 0.75:
				print('writing barcode, target_comp')
				# Decide to start with barcode or target
				if rand.random() < 0.5:
					for j in range(n):
						file.write(sticky_end)
						if i % 2 == 0:
							file.write(barcode)
						else:
							file.write(target_comp)
				# Start with target
				else:
					for j in range(n):
						file.write(sticky_end)
						if i % 2 == 0:
							file.write(target_comp)
						else:
							file.write(barcode)
			# Write barcode complement and target
			else:
				print('writing barcode_comp, target')
				# Decide to start with barcode or target
				if rand.random() < 0.5:
					for j in range(n):
						file.write(sticky_end)
						if i % 2 == 0:
							file.write(barcode_comp)
						else:
							file.write(target)
				# Start with target
				else:
					for j in range(n):
						file.write(sticky_end)
						if i % 2 == 0:
							file.write(target)
						else:
							file.write(barcode_comp)

			# Write the last sticky end
			file.write(sticky_end + '\n')

			# Write '+'
			file.write('+\n')

			# Write quality score
			if i < num_lines-1:
				file.write('qualityscore{}\n'.format(i))
			else:
				file.write('qualityscore{}'.format(i))

# # Generate test files with multi_processing
# # start timing
# start = time.perf_counter()

# with concurrent.futures.ProcessPoolExecutor() as executor:
# 	executor.map(generate_test_fastq_p, range(num_files))

# # end timing
# end = time.perf_counter()
# print('Multiprocess Time elapsed: {} second(s)'.format(end-start))

# generate_test_fastq_p(0)

# Generate test fastq with multiple targets
def make_Q_score(sequence):
	Q_score = ""
	for base in sequence:
		base_Q_score = rand.randint(33,126)
		Q_score += chr(base_Q_score)
	return Q_score

def datetime_2_str(datetime):
	time_str = f"{datetime.year}-{datetime.month:02d}-{datetime.day:02d}T{datetime.hour:02d}:{datetime.minute:02d}:{datetime.second:02d}Z"
	return time_str

def generate_ratio_test_fastq_p(i):
	with open(f'{save_folder}/test_{i}.fastq', 'w') as file:
		for i in range(num_lines):
			# Decide how many barcodes to add to sequence
			n = rand.randrange(2,21)
			# Decide whether to use sequence or complement
			which_target = rand.random()
			comp = rand.random()

			# Write sequence ID
			file.write(f'{std_sequence_id}\n')

			# Write sequence
			sequence = "TGCAGG"
			for j in range(n):
				# Decide which target
				if which_target <= ratio_threshold:
					target = 0
				else:
					target = 1
				# Decide target or complement
				if comp <= 0.5:
					sequence += target_list[target]
				else:
					sequence += comp_list[target]
				sequence += sticky_end
			# Remove the last 2 bases (simulating RE enzyme cut)
			sequence = sequence[:-2]			
			# Write the sequence
			file.write(f'{sequence}\n')

			# Write '+'
			file.write('+\n')

			# Write quality score
			Q_score = make_Q_score(sequence)
			if i < num_lines-1:
				file.write(f'{Q_score}\n')
			else:
				file.write(Q_score)

def generate_timed_ratio_test_fastq_p(i, num_reads, read_time, read_rate, read_step):
	with open(f'{save_folder}/test_{i}.fastq', 'w') as file:
		for i in range(num_lines):
			# Decide how many barcodes to add to sequence
			n = rand.randrange(2,21)
			# Decide whether to use sequence or complement
			which_target = rand.random()
			comp = rand.random()

			# Write sequence ID
			if num_reads % read_rate == 0 and num_reads != 0:
				read_time += dt.timedelta(minutes=read_step)
			else:
				pass
			time_str = datetime_2_str(read_time)
			file.write(f"@a631bd03-5b9d-44a5-9c01-af5ce3a7c77e runid=6bb769124712ca33f52fe35403bbb34fce14140a read=4 ch=84 start_time={time_str} flow_cell_id=AAA000 protocol_group_id=test_analyzeFASTQ sample_id={folder_name}\n")

			# Write sequence
			sequence = "TGCAGG"
			for j in range(n):
				# Decide which target
				if which_target <= ratio_threshold:
					target = 0
				else:
					target = 1
				# Decide target or complement
				if comp <= 0.5:
					sequence += target_list[target]
				else:
					sequence += comp_list[target]
				sequence += sticky_end
			# Remove the last 2 bases (simulating RE enzyme cut)
			sequence = sequence[:-2]			
			# Write the sequence
			file.write(f'{sequence}\n')

			# Write '+'
			file.write('+\n')

			# Write quality score
			Q_score = make_Q_score(sequence)
			if i < num_lines-1:
				file.write(f'{Q_score}\n')
			else:
				file.write(Q_score)
			# Increase the read count
			num_reads += 1
	return (num_reads, read_time)

ratio = 4 	# change this
folder_name = "1-2ratio" # change this
num_files = 24 # change this
ratio_threshold = 1/(1+ratio)
now = dt.datetime.now()

# Set up the save folder
outer_folder = "/Users/emilywu/OneDrive - Massachusetts Institute of Technology/Minion/test_analyzeFASTQ"
# save_folder = "/Users/emilywu/OneDrive - Massachusetts Institute of Technology/Minion/test_analyzeFASTQ/20211006_1734_MC-110826_0_AAA000_1-1ratio"
save_folder = f"{outer_folder}/{now.year}{now.month:02d}{now.day:02d}_{now.hour:02d}{now.minute:02d}_MC-110826_0_AAA000_{folder_name}"
if not os.path.exists(save_folder):
    os.makedirs(save_folder)
else:
	pass
# Set up generation information
target_list = ["CAATTCATCCATATTGCACCGTGAG", "CAGTGATTTTCGGCGGGCTCTAAAG"]	# barcodes A and B
comp_list = []
for target in target_list:
	comp_list.append(gen_complement(target))
std_sequence_id = f"@a631bd03-5b9d-44a5-9c01-af5ce3a7c77e runid=6bb769124712ca33f52fe35403bbb34fce14140a read=4 ch=84 start_time=2021-09-24T17:37:56Z flow_cell_id=AAA000 protocol_group_id=test_analyzeFASTQ sample_id={folder_name}"
sticky_end = "CCTGCAGG"


# Generate test files with multi_processing
# # start timing
# start = time.perf_counter()
# with concurrent.futures.ProcessPoolExecutor() as executor:
# 	executor.map(generate_ratio_test_fastq_p, range(num_files))
# end timing
# end = time.perf_counter()
# print('Multiprocess Time elapsed: {} second(s)'.format(end-start))

# start timing
start = time.perf_counter()
# Generate test files with different time
read_time = now
read_rate = 100
read_step = 5
num_reads = 0
for i in range(num_files):
	num_reads, read_time = generate_timed_ratio_test_fastq_p(i, num_reads, read_time, read_rate, read_step)
end = time.perf_counter()
print('Time elapsed: {} second(s)'.format(end-start))
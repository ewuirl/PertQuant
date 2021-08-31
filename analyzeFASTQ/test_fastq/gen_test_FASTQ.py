import sys
sys.path.append('../')
from simCRN.ml_nupack import gen_complement
import time
import random as rand 
import concurrent.futures

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

generate_test_fastq_p(0)
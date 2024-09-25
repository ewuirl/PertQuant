import sys
sys.path.append('../')
from simCRN.ml_nupack import gen_complement
import time
import functools
import concurrent.futures

barcode = 'TATGAGGACGAATCTCCCGCTTATA'
sticky_end = 'TGCA'
num_files = 100
file_names = []
for i in range(num_files):
	file_names.append('barcode_test{}.fastq'.format(i))

# This function reads an input file and saves results to an output file
def pick_save_file(save_file):
	"""
	This function returns a function count_barcodes(file_name) which will save 
	its results to a file using the provided save file name.
	Arguments:
		save_file (str): A string representing the name of the file to save 
		the results of count_barcodes to.
			ex: test.txt
	Returns:
		A function count_barcodes(file_name), which will save its results to
		the file name provided by save_file.
	"""	
	def check_barcode_counts(file_name):
		"""
		This function reads in a fastq file, calculates and counts the number of 
		barcodes in each sequence, figures out if these numbers match, and 
		writes the results to a save file. Each line is tab separated as 
		follows:
		@sequenceid    calculated barcodes    counted barcodes    match (True/False)
		Arguments:
			file_name (str): A string representing the name of the fastq file to
			read.
				ex: test.fastq
		"""
		with open(file_name, 'r') as file:
			with open(save_file, 'a') as results:
				i = 0
				# for i in range(10):
				while True:
					# Read the sequence ID
					seq_ID = file.readline().rstrip('\n')
					# If we've reached the end of the file, break the loop
					if not seq_ID:
						break
					# Read in the sequence and strip the new line
					seq = file.readline().rstrip('\n')
					file.readline()
					Q_score = file.readline().rstrip('\n')

					# Figure out how many barcodes are in the sequence
					num_barcodes = int((len(seq)-4)/29)
					# Split up the sequences
					seq_list = seq.split(sticky_end)
					barcode_count = len(seq_list) - 2
					# Write results to the save file
					results.write(seq_ID + '\t')
					results.write('{}\t'.format(num_barcodes))
					results.write('{}\t'.format(barcode_count))
					results.write('{}\n'.format(num_barcodes==barcode_count))
					# Increase the counter
					i += 1
	return check_barcode_counts

# process_file = pick_save_file('barcode_test_serial.txt')
# process_file_mp = pick_save_file('barcode_test_mp.txt')

# # Serial 
# # start timing
# start = time.perf_counter()

# for file_name in file_names:
# 	process_file(file_name)

# end = time.perf_counter()
# print('Serial time elapsed: {} second(s)'.format(end-start))


# We actually just want the sequence and counts
def count_barcodes(file_name):
	"""
	This function reads in a fastq file, calculates and returns the number of 
	barcode instances in it.
	Arguments:
		file_name (str): A string representing the name of the fastq file to
		read.
			ex: test.fastq
	Returns:
		num_barcodes (int): The number of times the barcode showed up in the 
		file
	"""
	with open(file_name, 'r') as file:
		# Save the number of barcodes
		num_barcodes = 0
		# for i in range(10):
		while True:
			# Read the sequence ID
			seq_ID = file.readline().rstrip('\n')
			# If we've reached the end of the file, break the loop
			if not seq_ID:
				break
			# Read in the sequence and strip the new line
			seq = file.readline().rstrip('\n')
			file.readline()
			Q_score = file.readline().rstrip('\n')

			# Figure out how many barcodes are in the sequence
			# Update the count
			num_barcodes += int((len(seq)-4)/29)		
	return num_barcodes


# Serial 
# start timing
start = time.perf_counter()

serial_count = 0
for file_name in file_names:
	serial_count += count_barcodes(file_name)

end = time.perf_counter()
print('Total barcode count: {}'. format(serial_count))
print('Serial time elapsed: {} second(s)'.format(end-start))

# Multiprocessed
# start timing
start = time.perf_counter()

parallel_count = 0
with concurrent.futures.ProcessPoolExecutor() as executor:
	results = executor.map(count_barcodes, file_names)
	print(results)
	for result in results:
		parallel_count += result

end = time.perf_counter()
print('Total barcode count: {}'. format(parallel_count))
print('Time elapsed: {} second(s)'.format(end-start))
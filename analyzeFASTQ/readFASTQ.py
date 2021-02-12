import config
import time
import functools
import concurrent.futures
import numpy as np
import jellyfish as jf

sticky_end = config.sticky_end
barcode = config.barcode
barcodec = config.barcodec
barcode1 = config.barcode1
target = config.target
targetc = config.targetc
target1 = config.target1

barcode_dict = config.barcode_dict
target_dict = config.target_dict

sticky_len = config.sticky_len
target_len = config.target_len
barcode_len = config.barcode_len

file_name = config.file_name

# # Create dictionary of barcodes
# with open('barcodes.txt', 'r') as barcode_file:
# 	# Create empty dictionary
# 	barcode_dict = {}
# 	barcodec_dict = {}
# 	# Create a counter 
# 	i = 0
# 	# Read lines
# 	while True:
# 		# Read line
# 		line = barcode_file.readline()
# 		# Break if the line is empty
# 		if not line:
# 			break
# 		# Remove the newline
# 		barcode = line.rstrip('\n')
# 		# Add the sequence and complement to the dictionary
# 		barcode_dict['barcode{}'.format(i)] = barcode
# 		barcodec_dict['barcode{}c'.format(i)] = gen_complement(seq)
		
# 		# Increase the counter
# 		i += 1

def analyze_file(file_name):
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

def analyze_sequence(seq, seq_id, count_array, error_list, tol):
	head = seq[:4]
	tail = seq[-4:]
	# If the ends match the sticky ends, remove them
	if head == sticky_end and tail == sticky_end:
		print('Stripping sticky ends')
		seq = seq[4:-4]
		print('stripped sequence: {}'.format(seq))
	# Else create a string of the fastq file name, sequence id, error type, and 
	# sequence to add to the error list
	else:
		print('sticky end error')
		error_list.append('{}\t{}\tend\t{}'.format(file_name, seq_id, seq))
		return('end_error')
	# Try splitting the sequence up into 
	print('splitting sequence')
	seq_list = seq.split(sticky_end)
	print('This is seq_list: {}'.format(seq_list))
	check_list = True
	# Set barcode ID flag to -1 (not identified):
	barcode_ID = -1
	# Check to make sure the fragments are all the correct sizes
	for fragment in seq_list:
		if len(fragment) == target_len or len(fragment) == barcode_len:
			pass
		# If a fragment isn't the right size, set the check_list variable to 
		# False
		else:
			check_list = False
	# If check_list is True, identify each fragment
	print('Check_list is {}'.format(check_list))
	if check_list:
		for fragment in seq_list:
			print('Checking sequence: {}'.format(fragment))
			# Check if the fragment is in the target dictionary
			if len(fragment) == target_len:
				print('Checking target')
				target_score = check_dict(fragment, target_dict, tol)
				# If there's a match, increase the count in the array
				if target_score[0] >= tol:
					print('Found target match')
					count_array[target_score[1]] += 1
				# If there's no match, create a string of the fastq file name, 
				# sequence ID, error type, fragment sequence, and best match 
				# score to add to the error list
				else:
					print('Didn\'t find target match')
					error_list.append('{}\t{}\tfrag\t{}\t{}'.format(file_name, \
						seq_id, fragment, target_score[0]))
			# Check if the fragment is in the barcode dictionary if the barcode
			# hasn't been ID'd yet
			elif len(fragment) == barcode_len and barcode_ID == -1:
				print('Finding barcode')
				barcode_score = check_dict(fragment, barcode_dict, tol)
				# If there's a match, set the barcode ID
				if barcode_score[0] >= tol:
					print('Found barcode')
					barcode_ID = barcode_score[1]
				else:
					print('Didn\'t find barcode')
					pass
			else:
				print('Passing on barcode')
				pass
	# If check list is false, use a sliding check
	else:
		print('sliiiiiide')
		# sliding check

	if barcode_ID < 0:
		error_list.append('{}\t{}\tbarcode\t{}'.format(file_name, seq_id, seq))
	else:
		pass	
	return barcode_ID


def check_dict(fragment, seq_dict, tol):
	""" Checks to see if a fragment is in the provided dictionary of sequences.
	If there is no perfect match, it calculates the Levenshtein distances 
	between the fragment and the dictionary, normalizes the distance to a score
	between [0,1] ([worst, perfect]) to see if there is a match with a score 
	high enough to pass the provided tolerance. The normalized score is: 
		1.0 - levenshtein_distance/max_levenshtein distance

	Arguments:
		fragment (str): a string specifying the fragment sequence to analyze
		seq_dict (dict): a dictionary of sequences to compare the fragment 
			against
		tol (float): a float between [0,1] to specify the accepted tolerance.
			1.0 = perfect match, 0.0 = worst match.

	Returns:
		(score, value) (tuple):
			score (float): the normalized levenshtein distance of the best 
				matching key in the dictionary
			value (int): the value of the best matching key in the dictionary
	"""
	# First check if there's a perfect match in the dictionary
	check_frag = seq_dict.get(fragment)
	# If there is, return the value of the corresponding match
	if check_frag != None:
		# print('found a perfect match')
		return (1.0, check_frag)
	# If not, calculate the Levenshtein distance between the fragment and the 
	# keys in the dictionary
	else:
		# print('calculating Levenshtein distances')
		# Max Levenshtein distance is the length of the longer sequence (same)
		max_dist = len(fragment)
		# Store the score in a tuple (lev_score, index, key)
		score = (0.0, -1)
		for key in seq_dict.keys():
			# print('trying key: {}'.format(seq_dict.get(key)))
			# Calculate the Levenshtein distance
			lev_dist = jf.levenshtein_distance(fragment, key)
			# Normalize the Levenshtein distance to [0, 1]
			lev_score = 1.0 - lev_dist / max_dist
			# If the score is greater than the tolerance, return the index of 
			# the corresponding key
			if lev_score >= tol:
				# print('found match within tolerance')
				return (lev_score, seq_dict.get(key))
			# If not, check if the score is better than the last saved score.
			else:
				# Save the score and index if the score is better
				if lev_score > score[0]:
					# print('found better match')
					score = (lev_score, seq_dict.get(key))
				else:
					# print('found worse match')
					pass
		# Getting to this point means a match wasn't found
		return score





def index_analyze_sequence(seq):
	pass



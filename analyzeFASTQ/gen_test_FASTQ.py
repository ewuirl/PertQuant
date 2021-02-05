import sys
sys.path.append('../')
from simCRN.ml_nupack import gen_complement
import time
import random as rand

barcode = 'TATGAGGACGAATCTCCCGCTTATA'
sticky_end = 'TGCA'
barcode_comp = gen_complement(barcode)
print(barcode_comp)
num_lines = 1000

# tim-starting
start = time.perf_counter()

# Open the file
with open('barcode_test.fastq', 'w') as file:
	for i in range(num_lines):
		# Decide how many barcodes to add to sequence
		n = rand.randrange(2,21)
		# Decide whether to use sequence or complement
		comp = rand.random()
		print(comp)

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

end = time.perf_counter()
print('Time elapsed: {} second(s)'.format(end-start))


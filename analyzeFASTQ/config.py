import sys
sys.path.append('../')
from simCRN.ml_nupack import gen_complement
# Define sequences
barcode = 'TATGAGGACGAATCTCCCGCTTATA'
barcodec = gen_complement(barcode)
barcode1 = 'GGTCTTGACAAACGTGTGCTTGTAC'
target = 'ATGGTCGAATCAAGGGGAGG'
targetc = gen_complement(target)
target1 = 'GGGAGCCTCTATCGTGTTGT'
sticky_end = 'TGCA'
file_name = 'test_file'

# Fill barcode dict with barcodes and comp
barcode_dict = {barcode: 0, barcode1: 1}
i = 0
for key in list(barcode_dict.keys()):
	comp = gen_complement(key)
	barcode_dict[comp] = i
	i += 1
target_dict = {target: 0, target1: 1}
i = 0

# Fill target dict with targets and comp
for key in list(target_dict.keys()):
	comp = gen_complement(key)
	target_dict[comp] = i
	i += 1

sticky_barcode_dict = {}
i = 0
for key in list(barcode_dict.keys()):
    slide_key = sticky_end + key + sticky_end
    sticky_barcode_dict[slide_key] = i

sticky_target_dict = {}
i = 0
for key in list(target_dict.keys()):
    slide_key = sticky_end + key + sticky_end
    sticky_target_dict[slide_key] = i

# Set lengths
sticky_len = 4
target_len = 20
barcode_len = 25
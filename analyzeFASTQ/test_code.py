import sys
sys.path.append('../')
from simCRN.ml_nupack import gen_complement
import config
from analyzeFASTQ.readFASTQ import check_match
from analyzeFASTQ.readFASTQ import analyze_sequence
from analyzeFASTQ.readFASTQ import replace_check
import numpy as np
import jellyfish as jf

sticky_end = config.sticky_end
barcode = config.barcode
barcodec = config.barcodec
barcode1 = config.barcode1
barcodec1 = gen_complement(barcode1)
target = config.target
targetc = config.targetc
target1 = config.target1
targetc1 = gen_complement(target1)

barcode_dict = config.barcode_dict
target_dict = config.target_dict

sticky_len = config.sticky_len
target_len = config.target_len
barcode_len = config.barcode_len


# barcode = 'TATGAGGACGAATCTCCCGCTTATA'
# sticky_end = 'TGCA'

# seq = sticky_end + barcode + sticky_end

# print(seq)
# print(seq.replace(sticky_end, 'hohoho'))
# print(seq.replace(sticky_end, ''))

# seq = sticky_end + 'cats' + sticky_end + 'dogs' + sticky_end

# print(seq)
# print(seq.replace(sticky_end, 'hohoho'))
# print(seq)
# print(seq.replace(sticky_end, ''))
# print(seq)
# print(seq[:4])
# print(seq[-4:])

# print(np.zeros((5,2)))


# print(barcode_dict.get(barcode))
# print(barcode_dict.get(target1))

# if barcode_dict.get(barcode) == None:
#     print('barcode not here')

# if barcode_dict.get(target) == None:
#     print('target not here')

barcode_switch = 'AATGAGGACGAATCTCCCGCTTATA'
barcode_ins = 'GTATGAGGACGAATCTCCCGCTTATA'
barcode_del = 'ATGAGGACGAATCTCCCGCTTATA'
barcode_swap = 'ATTGAGGACGAATCTCCCGCTTATA'
barcode1_switch = 
barcode1_ins = 
barcode1_del = 
barcode1_swap = 

# start = time.perf_counter()
# for i in range(100):
#     lev_dist_switch = jf.levenshtein_distance(barcode, barcode_switch)
#     lev_dist_ins = jf.levenshtein_distance(barcode, barcode_ins)
#     lev_dist_del = jf.levenshtein_distance(barcode, barcode_del)
#     lev_dist_swap = jf.levenshtein_distance(barcode, barcode_swap)
# end = time.perf_counter()
# print('Time elapsed{}'.format(end-start))
# print('Switch: {}'.format(lev_dist_switch))
# print('Insert: {}'.format(lev_dist_ins))
# print('Delete: {}'.format(lev_dist_del))
# print('Swap: {}'.format(lev_dist_swap))

# lev_dist_perf = 1-jf.levenshtein_distance(barcode, barcode)/max(len(barcode), len(barcode))
# lev_dist_switch = 1-jf.levenshtein_distance(barcode, barcode_switch)/max(len(barcode), len(barcode_switch))
# lev_dist_ins = 1-jf.levenshtein_distance(barcode, barcode_ins)/max(len(barcode), len(barcode_ins))
# lev_dist_del = 1-jf.levenshtein_distance(barcode, barcode_del)/max(len(barcode), len(barcode_del))
# lev_dist_swap = 1-jf.levenshtein_distance(barcode, barcode_swap)/max(len(barcode), len(barcode_swap))
# print('Perfect: {}'.format(lev_dist_perf))
# print('Switch: {}'.format(lev_dist_switch))
# print('Insert: {}'.format(lev_dist_ins))
# print('Delete: {}'.format(lev_dist_del))
# print('Swap: {}'.format(lev_dist_swap))

# lev_dist_perf = 1-jf.levenshtein_distance(barcode, barcode)/25.0
# lev_dist_switch = 1-jf.levenshtein_distance(barcode, barcode_switch)/25.0
# lev_dist_ins = 1-jf.levenshtein_distance(barcode, barcode_ins)/25.0
# lev_dist_del = 1-jf.levenshtein_distance(barcode, barcode_del)/25.0
# lev_dist_swap = 1-jf.levenshtein_distance(barcode, barcode_swap)/25.0
# print('Perfect: {}'.format(lev_dist_perf))
# print('Switch: {}'.format(lev_dist_switch))
# print('Insert: {}'.format(lev_dist_ins))
# print('Delete: {}'.format(lev_dist_del))
# print('Swap: {}'.format(lev_dist_swap))

# start = time.perf_counter()
# for i in range(100):
#     dam_lev_dist_switch = jf.damerau_levenshtein_distance(barcode, barcode_switch)
#     dam_lev_dist_ins = jf.damerau_levenshtein_distance(barcode, barcode_ins)
#     dam_lev_dist_del = jf.damerau_levenshtein_distance(barcode, barcode_del)
#     dam_lev_dist_swap = jf.damerau_levenshtein_distance(barcode, barcode_swap)
# end = time.perf_counter()
# print('Time elapsed{}'.format(end-start))
# print('Dam Switch: {}'.format(dam_lev_dist_switch))
# print('Dam Insert: {}'.format(dam_lev_dist_ins))
# print('Dam Delete: {}'.format(dam_lev_dist_del))
# print('Dam Swap: {}'.format(dam_lev_dist_swap))

# dam_lev_dist_perf = 1-jf.damerau_levenshtein_distance(barcode, barcode)/max(len(barcode), len(barcode))
# dam_lev_dist_switch = 1-jf.damerau_levenshtein_distance(barcode, barcode_switch)/max(len(barcode), len(barcode_switch))
# dam_lev_dist_ins = 1-jf.damerau_levenshtein_distance(barcode, barcode_ins)/max(len(barcode), len(barcode_ins))
# dam_lev_dist_del = 1-jf.damerau_levenshtein_distance(barcode, barcode_del)/max(len(barcode), len(barcode_del))
# dam_lev_dist_swap = 1-jf.damerau_levenshtein_distance(barcode, barcode_swap)/max(len(barcode), len(barcode_swap))
# print('dam Perfect: {}'.format(dam_lev_dist_perf))
# print('dam Switch: {}'.format(dam_lev_dist_switch))
# print('dam Insert: {}'.format(dam_lev_dist_ins))
# print('dam Delete: {}'.format(dam_lev_dist_del))
# print('dam Swap: {}'.format(dam_lev_dist_swap))

# # Test check_dict
# perf = check_dict(barcode, barcode_dict, 0.9)
# print('perf result: {}'.format(perf))
# switch = check_dict(barcode_switch, barcode_dict, 0.95)
# print('switch result: {}'.format(switch))
# swap = check_dict(barcode_swap, barcode_dict, 0.95)
# print('swap result: {}'.format(swap))

# Test analyze_sequence
# Create test sequences
barbar = sticky_end + barcode + sticky_end + barcode + sticky_end + barcode + sticky_end
barbarc = sticky_end + barcodec + sticky_end + barcodec + sticky_end + barcodec + sticky_end
barbar1 = sticky_end + barcode1 + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end
barbarc1 = sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end
bartar = sticky_end + barcode + sticky_end + targetc + sticky_end + barcode + sticky_end 
bartartar = sticky_end + barcode + sticky_end + targetc + sticky_end + barcode + sticky_end + targetc1 + sticky_end + barcode + sticky_end 
bartarc = sticky_end + barcodec + sticky_end + target + sticky_end + barcodec + sticky_end
bartartarc = sticky_end + barcodec + sticky_end + target + sticky_end + barcodec + sticky_end + target1 + sticky_end + barcodec + sticky_end
no_sticky_end = barcode + sticky_end + barcode + sticky_end + barcode
wrong_size = sticky_end + barcode + 'c' + sticky_end + barcode + sticky_end + barcode + sticky_end
diff_barcodes = sticky_end + barcode + sticky_end + barcode1 + sticky_end + barcode + sticky_end


test_sequence = wrong_size
tol = 0.9
count_array = np.zeros(2)
error_list = []
barbarid = analyze_sequence(test_sequence, '@seqid', count_array, error_list, tol)
print('count')
print(count_array)
print('error list:')
print(error_list)
print(barbarid)


# # Test replace_check
# test_sequence = barbar
# replaceseq = replace_check(test_sequence, barcode_dict, target_dict, 0.9)
# print('test sequence: {}'.format(test_sequence))
# print('replace seq: {}'.format(replaceseq))

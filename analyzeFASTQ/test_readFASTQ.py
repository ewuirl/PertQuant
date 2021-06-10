import sys
sys.path.append('../')
from simCRN.ml_nupack import gen_complement
import config
from analyzeFASTQ.readFASTQ import analyze_barcoded_seq
from analyzeFASTQ.readFASTQ import analyze_seq
from analyzeFASTQ.readFASTQ import make_dictionary
from analyzeFASTQ.readFASTQ import calc_tol
import numpy as np
import jellyfish as jf
import random as rand

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
sticky_barcode_dict = config.sticky_barcode_dict
sticky_target_dict = config.sticky_target_dict

sticky_len = config.sticky_len
target_len = config.target_len
barcode_len = config.barcode_len
sticky_barcode_len = int(2*sticky_len + barcode_len)
sticky_target_len = int(2*sticky_len + target_len)
barcode_lengths = (sticky_len, barcode_len, sticky_barcode_len) 
lengths = (sticky_len, barcode_len, target_len, sticky_barcode_len, \
    sticky_target_len) 

# Calculate the tolerances
target_tol_dist = 2
barcode_tol_dist = 2
target_tol = calc_tol(target_tol_dist, target_len)
barcode_tol = calc_tol(barcode_tol_dist, barcode_len)
sticky_target_tol = calc_tol(target_tol_dist, sticky_target_len)
sticky_barcode_tol = calc_tol(barcode_tol_dist, sticky_barcode_len)
# Group the tolerances into a tuple
barcode_tols = (barcode_tol, sticky_barcode_tol)
tols = (barcode_tol, target_tol, sticky_barcode_tol, sticky_target_tol)


barcode_ins = barcode + 'i'
barcode_del = barcode[1:]
barcode_swap = barcode[1] + barcode[0] + barcode[2:]
barcode_switch = 'j' + barcode[1:]
rand_b = 12
rand_barcode_ins = barcode[0:rand_b] + 'i' + barcode[rand_b:]
rand_barcode_del = barcode[0:rand_b] + barcode[rand_b+1:]
rand_barcode_swap = barcode[0:rand_b] + barcode[rand_b+1] \
+ barcode[rand_b] + barcode[rand_b+2:]
rand_barcode_switch = barcode[0:rand_b] + 'j' + barcode[rand_b+1:]

barcodec_ins = barcodec + 'i'
barcodec_del = barcodec[1:]
barcodec_swap = barcodec[1] + barcodec[0] + barcodec[2:]
barcodec_switch = 'j' + barcodec[1:]
rand_barcodec_ins = barcodec[0:rand_b] + 'i' + barcodec[rand_b:]
rand_barcodec_del = barcodec[0:rand_b] + barcodec[rand_b+1:]
rand_barcodec_swap = barcodec[0:rand_b] + barcodec[rand_b+1] \
+ barcodec[rand_b] + barcodec[rand_b+2:]
rand_barcodec_switch = barcodec[0:rand_b] + 'j' + barcodec[rand_b+1:]

barcode1_ins = barcode1 + 'i'
barcode1_del = barcode1[1:]
barcode1_swap = barcode1[1] + barcode1[0] + barcode1[2:]
barcode1_switch = 'j' + barcode1[1:]
rand_barcode1_ins = barcode1[0:rand_b] + 'i' + barcode1[rand_b:]
rand_barcode1_del = barcode1[0:rand_b] + barcode1[rand_b+1:]
rand_barcode1_swap = barcode1[0:rand_b] + barcode1[rand_b+1] \
+ barcode1[rand_b] + barcode1[rand_b+2:]
rand_barcode1_switch = barcode1[0:rand_b] + 'j' + barcode1[rand_b+1:]

barcodec1_ins = barcodec1 + 'i'
barcodec1_del = barcodec1[1:]
barcodec1_swap = barcodec1[1] + barcodec1[0] + barcodec1[2:]
barcodec1_switch = 'j' + barcodec1[1:]
rand_barcodec1_ins = barcodec1[0:rand_b] + 'i' + barcodec1[rand_b:]
rand_barcodec1_del = barcodec1[0:rand_b] + barcodec1[rand_b+1:]
rand_barcodec1_swap = barcodec1[0:rand_b] + barcodec1[rand_b+1] \
+ barcodec1[rand_b] + barcodec1[rand_b+2:]
rand_barcodec1_switch = barcodec1[0:rand_b] + 'j' + barcodec1[rand_b+1:]

target_ins = target + 'i'
target_del = target[1:]
target_swap = target[1] + target[0] + target[2:]
target_switch = 'j' + target[1:]
rand_t = 12
rand_target_ins = target[0:rand_t] + 'i' + target[rand_t:]
rand_target_del = target[0:rand_t] + target[rand_t+1:]
rand_target_swap = target[0:rand_t] + target[rand_t+1] \
+ target[rand_t] + target[rand_t+2:]
rand_target_switch = target[0:rand_t] + 'j' + target[rand_t+1:]

targetc_ins = targetc + 'i'
targetc_del = targetc[1:]
targetc_swap = targetc[1] + targetc[0] + targetc[2:]
targetc_switch = 'j' + targetc[1:]
rand_targetc_ins = targetc[0:rand_t] + 'i' + targetc[rand_t:]
rand_targetc_del = targetc[0:rand_t] + targetc[rand_t+1:]
rand_targetc_swap = targetc[0:rand_t] + targetc[rand_t+1] \
+ targetc[rand_t] + targetc[rand_t+2:]
rand_targetc_switch = targetc[0:rand_t] + 'j' + targetc[rand_t+1:]

target1_ins = target1 + 'i'
target1_del = target1[1:]
target1_swap = target1[1] + target1[0] + target1[2:]
target1_switch = 'j' + target1[1:]
rand_target1_ins = target1[0:rand_t] + 'i' + target1[rand_t:]
rand_target1_del = target1[0:rand_t] + target1[rand_t+1:]
rand_target1_swap = target1[0:rand_t] + target1[rand_t+1] \
+ target1[rand_t] + target1[rand_t+2:]
rand_target1_switch = target1[0:rand_t] + 'j' + target1[rand_t+1:]

targetc1_ins = targetc1 + 'i'
targetc1_del = targetc1[1:]
targetc1_swap = targetc1[1] + targetc1[0] + targetc1[2:]
targetc1_switch = 'j' + targetc1[1:]
rand_targetc1_ins = targetc1[0:rand_t] + 'i' + targetc1[rand_t:]
rand_targetc1_del = targetc1[0:rand_t] + targetc1[rand_t+1:]
rand_targetc1_swap = targetc1[0:rand_t] + targetc1[rand_t+1] \
+ targetc1[rand_t] + targetc1[rand_t+2:]
rand_targetc1_switch = targetc1[0:rand_t] + 'j' + targetc1[rand_t+1:]

sticky_del = sticky_end[:3]
sticky_ins = sticky_end + 'i'
sticky_swap = sticky_end[1] + sticky_end[0] + sticky_end[2:]
sticky_switch = 'j' + sticky_end[:3]

def test_barcoded_sequence(sequence, test_name, tols, answer):
    true_count, true_error, has_error, true_ID = answer
    count_array = np.zeros(2)
    error_list = []
    barcode_ID, is_error = analyze_barcoded_seq(sequence, '@seqid', 'testfile', \
        sticky_end, lengths, barcode_dict, target_dict, sticky_barcode_dict, \
        sticky_target_dict, count_array, error_list, tols)
    if np.array_equal(count_array, true_count):
        pass
    else:
        passed = 0
    if true_error == len(error_list):
        pass
    else:
        passed = 0
    if true_is_error == is_error:
        pass
    else:
        passed = 0
    if barcode_ID == true_ID:
        pass
    else:
        passed = 0

    return passed

def test_barcoded_sequence_print(sequence, test_name, tols, answer):
    true_count, true_error, has_error, true_ID = answer
    count_array = np.zeros(2)
    error_list = []
    barcode_ID, is_error = analyze_barcoded_seq(sequence, '@seqid', 'testfile', \
        sticky_end, lengths, barcode_dict, target_dict, sticky_barcode_dict, \
        sticky_target_dict, count_array, error_list, tols)
    passed = 1
    if np.array_equal(count_array, true_count):
        pass
    else:
        passed = 0
        print('Test: {} failed counts'.format(test_name))
        print('True counts: {}'.format(true_count))
        print('Calc counts: {}'.format(count_array))
    if true_error == len(error_list):
        pass
    else:
        passed = 0
        print('Test: {} failed # errors'.format(test_name))
        print('True error: {}'.format(true_error))
        print('Calc error: {}'.format(len(error_list)))
        print('{} errors: {}'.format(test_name, error_list))
    if has_error == is_error:
        pass
    else:
        print('Test: {} failed has error'.format(test_name))
        print('True has error: {}'.format(has_error))
        print('Calc has error: {}'.format(is_error))
        passed = 0
    if barcode_ID == true_ID:
        pass
    else:
        passed = 0
        print('Test: {} failed barcode'.format(test_name))
        print('True barcode: {}'.format(true_ID))
        print('Calc barcode: {}'.format(barcode_ID))
    return passed

def test_sequence(sequence, test_name, tols, answer):
    true_count, true_error, has_error, true_ID = answer
    count_array = np.zeros(2)
    error_list = []
    is_error = analyze_seq(sequence, '@seqid', 'testfile', \
        sticky_end, barcode_lengths, barcode_dict, sticky_barcode_dict, count_array, \
        error_list, tols)
    if np.array_equal(count_array, true_count):
        pass
    else:
        passed = 0
    if true_error == len(error_list):
        pass
    else:
        passed = 0
    if true_is_error == is_error:
        pass
    else:
        passed = 0

    return passed

def test_sequence_print(sequence, test_name, tols, answer):
    true_count, true_error, has_error, true_ID = answer
    count_array = np.zeros(2)
    error_list = []
    is_error = analyze_seq(sequence, '@seqid', 'testfile', \
        sticky_end, barcode_lengths, barcode_dict, sticky_barcode_dict, count_array, \
        error_list, tols)
    passed = 1
    if np.array_equal(count_array, true_count):
        pass
    else:
        passed = 0
        print('Test: {} failed counts'.format(test_name))
        print('True counts: {}'.format(true_count))
        print('Calc counts: {}'.format(count_array))
    if true_error == len(error_list):
        pass
    else:
        passed = 0
        print('Test: {} failed # errors'.format(test_name))
        print('True error: {}'.format(true_error))
        print('Calc error: {}'.format(len(error_list)))
        print('{} errors: {}'.format(test_name, error_list))
    if has_error == is_error:
        pass
    else:
        print('Test: {} failed has error'.format(test_name))
        print('True has error: {}'.format(has_error))
        print('Calc has error: {}'.format(is_error))
        passed = 0
    return passed

class testseq:
    def __init__(self, test_name, sequence, true_count, true_error, has_error, \
        true_ID):
        self.test_name = test_name
        self.sequence = sequence
        self.true_count = true_count
        self.true_error = true_error
        self.has_error = has_error
        self.true_ID = true_ID

# Perfect test cases
bar0 = testseq('barcode 0', \
    sticky_end + barcode + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.zeros(2), 0, False, 0)

barc0 = testseq('barcode 0 comp', \
    sticky_end + barcodec + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 0, False, 0)

bar1 = testseq('barcode 1', \
    sticky_end + barcode1 + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.zeros(2), 0, False, 1)

barc1 = testseq('barcode 1 comp', \
    sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.zeros(2), 0, False, 1)

bar0tar0 = testseq('barcode 0 and target 0', \
    sticky_end + barcode + sticky_end + targetc + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0tar1 = testseq('barcode 0 and target 1', \
    sticky_end + barcode + sticky_end + targetc1 + sticky_end + barcode + sticky_end, \
    np.array([0,1]), 0, False, 0)

bar0tar01 = testseq('barcode 0 and targets 0,1', \
    sticky_end + barcode + sticky_end + targetc + sticky_end + barcode + sticky_end + targetc1 + sticky_end + barcode + sticky_end, \
    np.array([1,1]), 0, False, 0)

barc0tar0 = testseq('barcode 0 comp and target 0', \
    sticky_end + barcodec + sticky_end + target + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0tar1 = testseq('barcode 0 comp and target 1', \
    sticky_end + barcodec + sticky_end + target1 + sticky_end + barcodec + sticky_end, \
    np.array([0,1]), 0, False, 0)

barc0tar01 = testseq('barcode 0 comp and targets 0,1', \
    sticky_end + barcodec + sticky_end + target + sticky_end + barcodec + sticky_end + target1 + sticky_end + barcodec + sticky_end, \
    np.array([1,1]), 0, False, 0)

bar1tar0 = testseq('barcode 1 and target 0', \
    sticky_end + barcode1 + sticky_end + targetc + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1tar1 = testseq('barcode 1 and target 1', \
    sticky_end + barcode1 + sticky_end + targetc1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

bar1tar01 = testseq('barcode 1 and targets 0,1', \
    sticky_end + barcode1 + sticky_end + targetc + sticky_end + barcode1 + sticky_end + targetc1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,1]), 0, False, 1)

barc1tar0 = testseq('barcode 1 comp and target 0', \
    sticky_end + barcodec1 + sticky_end + target + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1tar1 = testseq('barcode 1 comp and target 1', \
    sticky_end + barcodec1 + sticky_end + target1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

barc1tar01 = testseq('barcode 1 comp and targets 0,1', \
    sticky_end + barcodec1 + sticky_end + target + sticky_end + barcodec1 + sticky_end + target1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,1]), 0, False, 1)

# Approximate cases
# Approximate cases: Barcode 0 only, 1 barcode error
bar0_del = testseq('barcode 0 deletion', \
    sticky_end + barcode_del + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.zeros(2), 0, False, 0)

bar0_ins = testseq('barcode 0 insertion', \
    sticky_end + barcode_ins + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.zeros(2), 0, False, 0)

bar0_switch = testseq('barcode 0 switch', \
    sticky_end + barcode_switch + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.zeros(2), 0, False, 0)

bar0_swap = testseq('barcode 0 swap', \
    sticky_end + barcode_swap + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.zeros(2), 0, False, 0)

barc0_del = testseq('barcode 0 comp deletion', \
    sticky_end + barcodec_del + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 0, False, 0)

barc0_ins = testseq('barcode 0 comp insertion', \
    sticky_end + barcodec_ins + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 0, False, 0)

barc0_switch = testseq('barcode 0 comp switch', \
    sticky_end + barcodec_switch + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 0, False, 0)

barc0_swap = testseq('barcode 0 comp swap', \
    sticky_end + barcodec_swap + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 0, False, 0)

# Approximate cases: Barcode 1 only, 1 barcode error
bar1_del = testseq('barcode 1 deletion', \
    sticky_end + barcode1_del + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.zeros(2), 0, False, 1)

bar1_ins = testseq('barcode 1 insertion', \
    sticky_end + barcode1_ins + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.zeros(2), 0, False, 1)

bar1_switch = testseq('barcode 1 switch', \
    sticky_end + barcode1_switch + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.zeros(2), 0, False, 1)

bar1_swap = testseq('barcode 1 swap', \
    sticky_end + barcode1_swap + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.zeros(2), 0, False, 1)

barc1_del = testseq('barcode 1 comp deletion', \
    sticky_end + barcodec1_del + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.zeros(2), 0, False, 1)

barc1_ins = testseq('barcode 1 comp insertion', \
    sticky_end + barcodec1_ins + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.zeros(2), 0, False, 1)

barc1_switch = testseq('barcode 1 comp switch', \
    sticky_end + barcodec1_switch + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.zeros(2), 0, False, 1)

barc1_swap = testseq('barcode 1 comp swap', \
    sticky_end + barcodec1_swap + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.zeros(2), 0, False, 1)

# Approximate cases: Barcode 0 only, 1 random barcode error
bar0_del_rand = testseq('barcode 0 deletion rand', \
    sticky_end + rand_barcode_del + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.zeros(2), 0, False, 0)

bar0_ins_rand = testseq('barcode 0 insertion rand', \
    sticky_end + rand_barcode_ins + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.zeros(2), 0, False, 0)

bar0_switch_rand = testseq('barcode 0 switch rand', \
    sticky_end + rand_barcode_switch + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.zeros(2), 0, False, 0)

bar0_swap_rand = testseq('barcode 0 swap rand', \
    sticky_end + rand_barcode_swap + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.zeros(2), 0, False, 0)

barc0_del_rand = testseq('barcode 0 comp deletion rand', \
    sticky_end + rand_barcodec_del + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 0, False, 0)

barc0_ins_rand = testseq('barcode 0 comp insertion rand', \
    sticky_end + rand_barcodec_ins + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 0, False, 0)

barc0_switch_rand = testseq('barcode 0 comp switch rand', \
    sticky_end + rand_barcodec_switch + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 0, False, 0)

barc0_swap_rand = testseq('barcode 0 comp swap rand', \
    sticky_end + rand_barcodec_swap + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 0, False, 0)

# Approximate cases: Barcode 1 only, 1 random barcode error
bar1_del_rand = testseq('barcode 1 deletion rand', \
    sticky_end + rand_barcode1_del + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.zeros(2), 0, False, 1)

bar1_ins_rand = testseq('barcode 1 insertion rand', \
    sticky_end + rand_barcode1_ins + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.zeros(2), 0, False, 1)

bar1_switch_rand = testseq('barcode 1 switch rand', \
    sticky_end + rand_barcode1_switch + sticky_end + barcode1 + sticky_end + barcode1+ sticky_end, \
    np.zeros(2), 0, False, 1)

bar1_swap_rand = testseq('barcode 1 swap rand', \
    sticky_end + rand_barcode1_swap + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.zeros(2), 0, False, 1)

barc1_del_rand = testseq('barcode 1 comp deletion rand', \
    sticky_end + rand_barcodec1_del + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.zeros(2), 0, False, 1)

barc1_ins_rand = testseq('barcode 1 comp insertion rand', \
    sticky_end + rand_barcodec1_ins + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.zeros(2), 0, False, 1)

barc1_switch_rand = testseq('barcode 1 comp switch rand', \
    sticky_end + rand_barcodec1_switch + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.zeros(2), 0, False, 1)

barc1_swap_rand = testseq('barcode 1 comp swap rand', \
    sticky_end + rand_barcodec1_swap + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.zeros(2), 0, False, 1)

# Approximate cases: Barcode 0 + Target errors
bar0_tar0_del = testseq('barcode 0 target 0 deletion', \
    sticky_end + targetc_del + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_tar0_ins = testseq('barcode 0 target 0 insertion', \
    sticky_end + targetc_ins + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_tar0_switch = testseq('barcode 0 target 0 switch', \
    sticky_end + targetc_switch + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_tar0_swap = testseq('barcode 0 target 0 swap', \
    sticky_end + targetc_swap + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_tar0_del = testseq('barcode 0 comp target 0 deletion', \
    sticky_end + target_del + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_tar0_ins = testseq('barcode 0 comp target 0 insertion', \
    sticky_end + target_ins + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_tar0_switch = testseq('barcode 0 comp target 0 switch', \
    sticky_end + target_switch + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_tar0_swap = testseq('barcode 0 comp target 0 swap', \
    sticky_end + target_swap + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_tar1_del = testseq('barcode 0 target 1 deletion', \
    sticky_end + targetc1_del + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([0,1]), 0, False, 0)

bar0_tar1_ins = testseq('barcode 0 target 1 insertion', \
    sticky_end + targetc1_ins + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([0,1]), 0, False, 0)

bar0_tar1_switch = testseq('barcode 0 target 1 switch', \
    sticky_end + targetc1_switch + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([0,1]), 0, False, 0)

bar0_tar1_swap = testseq('barcode 0 target 1 swap', \
    sticky_end + targetc1_swap + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([0,1]), 0, False, 0)

barc0_tar1_del = testseq('barcode 0 comp target 1 deletion', \
    sticky_end + target1_del + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([0,1]), 0, False, 0)

barc0_tar1_ins = testseq('barcode 0 comp target 1 insertion', \
    sticky_end + target1_ins + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([0,1]), 0, False, 0)

barc0_tar1_switch = testseq('barcode 0 comp target 1 switch', \
    sticky_end + target1_switch + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([0,1]), 0, False, 0)

barc0_tar1_swap = testseq('barcode 0 comp target 1 swap', \
    sticky_end + target1_swap + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([0,1]), 0, False, 0)

# Approximate cases: Barcode + Target errors
bar1_tar0_del = testseq('barcode 1 target 0 deletion', \
    sticky_end + targetc_del + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_tar0_ins = testseq('barcode 1 target 0 insertion', \
    sticky_end + targetc_ins + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_tar0_switch = testseq('barcode 1 target 0 switch', \
    sticky_end + targetc_switch + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_tar0_swap = testseq('barcode 1 target 0 swap', \
    sticky_end + targetc_swap + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)
    
barc1_tar0_del = testseq('barcode 1 comp target 0 deletion', \
    sticky_end + target_del + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_tar0_ins = testseq('barcode 1 comp target 0 insertion', \
    sticky_end + target_ins + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_tar0_switch = testseq('barcode 1 comp target 0 switch', \
    sticky_end + target_switch + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_tar0_swap = testseq('barcode 1 comp target 0 swap', \
    sticky_end + target_swap + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_tar1_del = testseq('barcode 1 target 1 deletion', \
    sticky_end + targetc1_del + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

bar1_tar1_ins = testseq('barcode 1 target 1 insertion', \
    sticky_end + targetc1_ins + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

bar1_tar1_switch = testseq('barcode 1 target 1 switch', \
    sticky_end + targetc1_switch + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

bar1_tar1_swap = testseq('barcode 1 target 1 swap', \
    sticky_end + targetc1_swap + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

barc1_tar1_del = testseq('barcode 1 comp target 1 deletion', \
    sticky_end + target1_del + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

barc1_tar1_ins = testseq('barcode 1 comp target 1 insertion', \
    sticky_end + target1_ins + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

barc1_tar1_switch = testseq('barcode 1 comp target 1 switch', \
    sticky_end + target1_switch + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

barc1_tar1_swap = testseq('barcode 1 comp target 1 swap', \
    sticky_end + target1_swap + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

# Approximation cases: Barcode 0 + Target random errors
bar0_tar0_del_rand = testseq('barcode 0 target 0 rand deletion', \
    sticky_end + rand_targetc_del + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_tar0_ins_rand = testseq('barcode 0 target 0 rand insertion', \
    sticky_end + rand_targetc_ins + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_tar0_switch_rand = testseq('barcode 0 target 0 rand switch', \
    sticky_end + rand_targetc_switch + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_tar0_swap_rand = testseq('barcode 0 target 0 rand swap', \
    sticky_end + rand_targetc_swap + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_tar0_del_rand = testseq('barcode comp 0 target 0 rand deletion', \
    sticky_end + rand_target_del + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_tar0_ins_rand = testseq('barcode comp 0 target 0 rand insertion', \
    sticky_end + rand_target_ins + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_tar0_switch_rand = testseq('barcode comp 0 target 0 rand switch', \
    sticky_end + rand_target_switch + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_tar0_swap_rand = testseq('barcode comp 0 target 0 rand swap', \
    sticky_end + rand_target_swap + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_tar1_del_rand = testseq('barcode 0 target 1 rand deletion', \
    sticky_end + rand_targetc1_del + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([0,1]), 0, False, 0)

bar0_tar1_ins_rand = testseq('barcode 0 target 1 rand insertion', \
    sticky_end + rand_targetc1_ins + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([0,1]), 0, False, 0)

bar0_tar1_switch_rand = testseq('barcode 0 target 1 rand switch', \
    sticky_end + rand_targetc1_switch + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([0,1]), 0, False, 0)

bar0_tar1_swap_rand = testseq('barcode 0 target 1 rand swap', \
    sticky_end + rand_targetc1_swap + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([0,1]), 0, False, 0)

barc0_tar1_del_rand = testseq('barcode comp 0 target 1 rand deletion', \
    sticky_end + rand_target1_del + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([0,1]), 0, False, 0)

barc0_tar1_ins_rand = testseq('barcode comp 0 target 1 rand insertion', \
    sticky_end + rand_target1_ins + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([0,1]), 0, False, 0)

barc0_tar1_switch_rand = testseq('barcode comp 0 target 1 rand switch', \
    sticky_end + rand_target1_switch + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([0,1]), 0, False, 0)

barc0_tar1_swap_rand = testseq('barcode comp 0 target 1 rand swap', \
    sticky_end + rand_target1_swap + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([0,1]), 0, False, 0)

# Approximation cases: Barcode 1 + Target random errors
bar1_tar0_del_rand = testseq('barcode 1 target 0 rand deletion', \
    sticky_end + rand_targetc_del + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_tar0_ins_rand = testseq('barcode 1 target 0 rand insertion', \
    sticky_end + rand_targetc_ins + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_tar0_switch_rand = testseq('barcode 1 target 0 rand switch', \
    sticky_end + rand_targetc_switch + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_tar0_swap_rand = testseq('barcode 1 target 0 rand swap', \
    sticky_end + rand_targetc_swap + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_tar0_del_rand = testseq('barcode 1 comp target 0 rand deletion', \
    sticky_end + rand_target_del + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_tar0_ins_rand = testseq('barcode 1 comp target 0 rand insertion', \
    sticky_end + rand_target_ins + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_tar0_switch_rand = testseq('barcode 1 comp target 0 rand switch', \
    sticky_end + rand_target_switch + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_tar0_swap_rand = testseq('barcode 1 comp target 0 rand swap', \
    sticky_end + rand_target_swap + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_tar1_del_rand = testseq('barcode 1 target 1 rand deletion', \
    sticky_end + rand_targetc1_del + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

bar1_tar1_ins_rand = testseq('barcode 1 target 1 rand insertion', \
    sticky_end + rand_targetc1_ins + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

bar1_tar1_switch_rand = testseq('barcode 1 target 1 rand switch', \
    sticky_end + rand_targetc1_switch + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

bar1_tar1_swap_rand = testseq('barcode 1 target 1 rand swap', \
    sticky_end + rand_targetc1_swap + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

barc1_tar1_del_rand = testseq('barcode 1 comp target 1 rand deletion', \
    sticky_end + rand_target1_del + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

barc1_tar1_ins_rand = testseq('barcode 1 comp target 1 rand insertion', \
    sticky_end + rand_target1_ins + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

barc1_tar1_switch_rand = testseq('barcode 1 comp target 1 rand switch', \
    sticky_end + rand_target1_switch + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

barc1_tar1_swap_rand = testseq('barcode 1 comp target 1 rand swap', \
    sticky_end + rand_target1_swap + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,1]), 0, False, 1)

# Approximate cases: Barcode 0 random error + Random target 0 deletion
bar0_del_tar0_del_rand = testseq('barcode 0 deletion target 0 rand deletion', \
    sticky_end + rand_target_del + sticky_end + rand_barcode_del + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_ins_tar0_del_rand = testseq('barcode 0 insertion target 0 rand deletion', \
    sticky_end + rand_target_del + sticky_end + rand_barcode_ins + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_switch_tar0_del_rand = testseq('barcode 0 switch target 0 rand deletion', \
    sticky_end + rand_target_del + sticky_end + rand_barcode_switch + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_swap_tar0_del_rand = testseq('barcode 0 swap target 0 rand deletion', \
    sticky_end + rand_target_del + sticky_end + rand_barcode_swap + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_del_tar0_del_rand = testseq('barcode 0 comp deletion target 0 rand deletion', \
    sticky_end + rand_targetc_del + sticky_end + rand_barcode_del + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_ins_tar0_del_rand = testseq('barcode 0 comp insertion target 0 rand deletion', \
    sticky_end + rand_targetc_del + sticky_end + rand_barcode_ins + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_switch_tar0_del_rand = testseq('barcode 0 comp switch target 0 rand deletion', \
    sticky_end + rand_targetc_del + sticky_end + rand_barcode_switch + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_swap_tar0_del_rand = testseq('barcode 0 comp swap target 0 rand deletion', \
    sticky_end + rand_targetc_del + sticky_end + rand_barcode_swap + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

# Approximate cases: Barcode 0 random error + Random target 0 insertion
bar0_del_tar0_ins_rand = testseq('barcode 0 deletion target 0 rand insertion', \
    sticky_end + rand_target_ins + sticky_end + rand_barcode_del + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_ins_tar0_ins_rand = testseq('barcode 0 insertion target 0 rand insertion', \
    sticky_end + rand_target_ins + sticky_end + rand_barcode_ins + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_switch_tar0_ins_rand = testseq('barcode 0 switch target 0 rand insertion', \
    sticky_end + rand_target_ins + sticky_end + rand_barcode_switch + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_swap_tar0_ins_rand = testseq('barcode 0 swap target 0 rand insertion', \
    sticky_end + rand_target_ins + sticky_end + rand_barcode_swap + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_del_tar0_ins_rand = testseq('barcode 0 comp deletion target 0 rand insertion', \
    sticky_end + rand_targetc_ins + sticky_end + rand_barcode_del + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_ins_tar0_ins_rand = testseq('barcode 0 comp insertion target 0 rand insertion', \
    sticky_end + rand_targetc_ins + sticky_end + rand_barcode_ins + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_switch_tar0_ins_rand = testseq('barcode 0 comp switch target 0 rand insertion', \
    sticky_end + rand_targetc_ins + sticky_end + rand_barcode_switch + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_swap_tar0_ins_rand = testseq('barcode 0 comp swap target 0 rand insertion', \
    sticky_end + rand_targetc_ins + sticky_end + rand_barcode_swap + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

# Approximate cases: Barcode 0 random error + Random target 0 switch
bar0_del_tar0_switch_rand = testseq('barcode 0 deletion target 0 rand switch', \
    sticky_end + rand_target_switch + sticky_end + rand_barcode_del + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_ins_tar0_switch_rand = testseq('barcode 0 insertion target 0 rand switch', \
    sticky_end + rand_target_switch + sticky_end + rand_barcode_ins + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_switch_tar0_switch_rand = testseq('barcode 0 switch target 0 rand switch', \
    sticky_end + rand_target_switch + sticky_end + rand_barcode_switch + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_swap_tar0_switch_rand = testseq('barcode 0 swap target 0 rand switch', \
    sticky_end + rand_target_switch + sticky_end + rand_barcode_swap + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_del_tar0_switch_rand = testseq('barcode 0 comp deletion target 0 rand switch', \
    sticky_end + rand_targetc_switch + sticky_end + rand_barcode_del + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_ins_tar0_switch_rand = testseq('barcode 0 comp insertion target 0 rand switch', \
    sticky_end + rand_targetc_switch + sticky_end + rand_barcode_ins + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_switch_tar0_switch_rand = testseq('barcode 0 comp switch target 0 rand switch', \
    sticky_end + rand_targetc_switch + sticky_end + rand_barcode_switch + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_swap_tar0_switch_rand = testseq('barcode 0 comp swap target 0 rand switch', \
    sticky_end + rand_targetc_switch + sticky_end + rand_barcode_swap + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

# Approximate cases: Barcode 0 random error + Random target 0 swap
bar0_del_tar0_swap_rand = testseq('barcode 0 deletion target 0 rand swap', \
    sticky_end + rand_target_swap + sticky_end + rand_barcode_del + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_ins_tar0_swap_rand = testseq('barcode 0 insertion target 0 rand swap', \
    sticky_end + rand_target_swap + sticky_end + rand_barcode_ins + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_switch_tar0_swap_rand = testseq('barcode 0 switch target 0 rand swap', \
    sticky_end + rand_target_swap + sticky_end + rand_barcode_switch + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

bar0_swap_tar0_swap_rand = testseq('barcode 0 swap target 0 rand swap', \
    sticky_end + rand_target_swap + sticky_end + rand_barcode_swap + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_del_tar0_swap_rand = testseq('barcode 0 comp deletion target 0 rand swap', \
    sticky_end + rand_targetc_swap + sticky_end + rand_barcode_del + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_ins_tar0_swap_rand = testseq('barcode 0 comp insertion target 0 rand swap', \
    sticky_end + rand_targetc_swap + sticky_end + rand_barcode_ins + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_switch_tar0_swap_rand = testseq('barcode 0 comp switch target 0 rand swap', \
    sticky_end + rand_targetc_swap + sticky_end + rand_barcode_switch + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

barc0_swap_tar0_swap_rand = testseq('barcode 0 comp swap target 0 rand swap', \
    sticky_end + rand_targetc_swap + sticky_end + rand_barcode_swap + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 0, False, 0)

# Approximate cases: Barcode 1 random error + Random target 0 deletion
bar1_del_tar0_del_rand = testseq('barcode 1 deletion target 0 rand deletion', \
    sticky_end + rand_target_del + sticky_end + rand_barcode1_del + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_ins_tar0_del_rand = testseq('barcode 1 insertion target 0 rand deletion', \
    sticky_end + rand_target_del + sticky_end + rand_barcode1_ins + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_switch_tar0_del_rand = testseq('barcode 1 switch target 0 rand deletion', \
    sticky_end + rand_target_del + sticky_end + rand_barcode1_switch + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_swap_tar0_del_rand = testseq('barcode 1 swap target 0 rand deletion', \
    sticky_end + rand_target_del + sticky_end + rand_barcode1_swap + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_del_tar0_del_rand = testseq('barcode 1 comp deletion target 0 rand deletion', \
    sticky_end + rand_targetc_del + sticky_end + rand_barcode1_del + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_ins_tar0_del_rand = testseq('barcode 1 comp insertion target 0 rand deletion', \
    sticky_end + rand_targetc_del + sticky_end + rand_barcode1_ins + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_switch_tar0_del_rand = testseq('barcode 1 comp switch target 0 rand deletion', \
    sticky_end + rand_targetc_del + sticky_end + rand_barcode1_switch + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_swap_tar0_del_rand = testseq('barcode 1 comp swap target 0 rand deletion', \
    sticky_end + rand_targetc_del + sticky_end + rand_barcode1_swap + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

# Approximate cases: Barcode 1 random error + Random target 0 insertion
bar1_del_tar0_ins_rand = testseq('barcode 1 deletion target 0 rand insertion', \
    sticky_end + rand_target_ins + sticky_end + rand_barcode1_del + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_ins_tar0_ins_rand = testseq('barcode 1 insertion target 0 rand insertion', \
    sticky_end + rand_target_ins + sticky_end + rand_barcode1_ins + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_switch_tar0_ins_rand = testseq('barcode 1 switch target 0 rand insertion', \
    sticky_end + rand_target_ins + sticky_end + rand_barcode1_switch + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_swap_tar0_ins_rand = testseq('barcode 1 swap target 0 rand insertion', \
    sticky_end + rand_target_ins + sticky_end + rand_barcode1_swap + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_del_tar0_ins_rand = testseq('barcode 1 comp deletion target 0 rand insertion', \
    sticky_end + rand_targetc_ins + sticky_end + rand_barcode1_del + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_ins_tar0_ins_rand = testseq('barcode 1 comp insertion target 0 rand insertion', \
    sticky_end + rand_targetc_ins + sticky_end + rand_barcode1_ins + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_switch_tar0_ins_rand = testseq('barcode 1 comp switch target 0 rand insertion', \
    sticky_end + rand_targetc_ins + sticky_end + rand_barcode1_switch + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_swap_tar0_ins_rand = testseq('barcode 1 comp swap target 0 rand insertion', \
    sticky_end + rand_targetc_ins + sticky_end + rand_barcode1_swap + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

# Approximate cases: Barcode 1 random error + Random target 0 switch
bar1_del_tar0_switch_rand = testseq('barcode 1 deletion target 0 rand switch', \
    sticky_end + rand_target_switch + sticky_end + rand_barcode1_del + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_ins_tar0_switch_rand = testseq('barcode 1 insertion target 0 rand switch', \
    sticky_end + rand_target_switch + sticky_end + rand_barcode1_ins + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_switch_tar0_switch_rand = testseq('barcode 1 switch target 0 rand switch', \
    sticky_end + rand_target_switch + sticky_end + rand_barcode1_switch + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_swap_tar0_switch_rand = testseq('barcode 1 swap target 0 rand switch', \
    sticky_end + rand_target_switch + sticky_end + rand_barcode1_swap + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_del_tar0_switch_rand = testseq('barcode 1 comp deletion target 0 rand switch', \
    sticky_end + rand_targetc_switch + sticky_end + rand_barcode1_del + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_ins_tar0_switch_rand = testseq('barcode 1 comp insertion target 0 rand switch', \
    sticky_end + rand_targetc_switch + sticky_end + rand_barcode1_ins + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_switch_tar0_switch_rand = testseq('barcode 1 comp switch target 0 rand switch', \
    sticky_end + rand_targetc_switch + sticky_end + rand_barcode1_switch + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_swap_tar0_switch_rand = testseq('barcode 1 comp swap target 0 rand switch', \
    sticky_end + rand_targetc_switch + sticky_end + rand_barcode1_swap + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

# Approximate cases: Barcode 1 random error + Random target 0 swap
bar1_del_tar0_swap_rand = testseq('barcode 1 deletion target 0 rand swap', \
    sticky_end + rand_target_swap + sticky_end + rand_barcode1_del + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_ins_tar0_swap_rand = testseq('barcode 1 insertion target 0 rand swap', \
    sticky_end + rand_target_swap + sticky_end + rand_barcode1_ins + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_switch_tar0_swap_rand = testseq('barcode 1 switch target 0 rand swap', \
    sticky_end + rand_target_swap + sticky_end + rand_barcode1_switch + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

bar1_swap_tar0_swap_rand = testseq('barcode 1 swap target 0 rand swap', \
    sticky_end + rand_target_swap + sticky_end + rand_barcode1_swap + sticky_end + barcode1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_del_tar0_swap_rand = testseq('barcode 1 comp deletion target 0 rand swap', \
    sticky_end + rand_targetc_swap + sticky_end + rand_barcode1_del + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_ins_tar0_swap_rand = testseq('barcode 1 comp insertion target 0 rand swap', \
    sticky_end + rand_targetc_swap + sticky_end + rand_barcode1_ins + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_switch_tar0_swap_rand = testseq('barcode 1 comp switch target 0 rand swap', \
    sticky_end + rand_targetc_swap + sticky_end + rand_barcode1_switch + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

barc1_swap_tar0_swap_rand = testseq('barcode 1 comp swap target 0 rand swap', \
    sticky_end + rand_targetc_swap + sticky_end + rand_barcode1_swap + sticky_end + barcodec1 + sticky_end, \
    np.array([1,0]), 0, False, 1)

# Approximation cases: Sticky end errors
bar0_sticky_del = testseq('barcode 0 sticky end deletion', \
    sticky_end + barcodec + sticky_del + barcodec + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 0, False, 0)

bar0_sticky_ins = testseq('barcode 0 sticky end insertion', \
    sticky_end + barcodec + sticky_ins + barcodec + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 0, False, 0)

bar0_sticky_switch = testseq('barcode 0 sticky end switch', \
    sticky_end + barcodec + sticky_switch + barcodec + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 0, False, 0)

bar0_sticky_swap = testseq('barcode 0 sticky end swap', \
    sticky_end + barcodec + sticky_swap + barcodec + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 0, False, 0)

# Error cases
mixed_bar = testseq('mixed barcode', \
    sticky_end + barcode + sticky_end + barcode1 + sticky_end + barcode + sticky_end, \
    np.zeros(2), 1, True, -1)

mixed_barc = testseq('mixed barcode comp', \
    sticky_end + barcodec + sticky_end + barcodec1 + sticky_end + barcodec + sticky_end, \
    np.zeros(2), 1, True, -1)

mixed_bar_tar = testseq('mixed barcode with target 0', \
    sticky_end + barcode + sticky_end + targetc + sticky_end + barcode1 + sticky_end + barcode + sticky_end, \
    np.array([1,0]), 1, True, -1)

mixed_barc_tar = testseq('mixed barcode comp with target 0', \
    sticky_end + barcodec + sticky_end + target + sticky_end + barcodec1 + sticky_end + barcodec + sticky_end, \
    np.array([1,0]), 1, True, -1)

no_bar = testseq('no barcode', \
    sticky_end + target + sticky_end + target + sticky_end + target + sticky_end, \
    np.array([3,0]), 1, True, -2)

unID_subseq = testseq('unidentified subsequence', \
    sticky_end + 'whatisthis' + sticky_end + barcodec + sticky_end + target + sticky_end, \
    np.array([1,0]), 1, True, 0)

unID_endseq = testseq('unidentified end subsequence', \
    sticky_end + barcodec + sticky_end + barcodec + sticky_end + target + sticky_end \
    + 'whatisthis', np.array([1,0]), 1, True, 0)

short = testseq('short', 'whatisthis', np.zeros(2), 2, True, -2)

unidentified_seq = testseq('unidentified sequence', \
    'whatisthiswhatisthiswhatisthiswhatisthiswhatisthis', np.array([0,0]), 2, True, -2)

perfect_seq_list = [bar0, barc0, bar1, barc1, \
bar0tar0, bar0tar1, bar0tar01, \
barc0tar0, barc0tar1, barc0tar01, \
bar1tar0, bar1tar1, bar1tar01, \
barc1tar0, barc1tar1, barc1tar01] 

approx_seq_list = [bar0_del, bar0_ins, bar0_switch, bar0_swap, \
barc0_del, barc0_ins, barc0_switch, barc0_swap, \
bar1_del, bar1_ins, bar1_switch, bar1_swap, \
barc1_del, barc1_ins, barc1_switch, barc1_swap, \
bar0_del_rand, bar0_ins_rand, bar0_switch_rand, bar0_swap_rand, \
barc0_del_rand, barc0_ins_rand, barc0_switch_rand, barc0_swap_rand, \
bar1_del_rand, bar1_ins_rand, bar1_switch_rand, bar1_swap_rand, \
barc1_del_rand, barc1_ins_rand, barc1_switch_rand, barc1_swap_rand, \
bar0_tar0_del, bar0_tar0_ins, bar0_tar0_switch, bar0_tar0_swap, \
barc0_tar0_del, barc0_tar0_ins, barc0_tar0_switch, barc0_tar0_swap, \
bar0_tar1_del, bar0_tar1_ins, bar0_tar1_switch, bar0_tar1_swap, \
barc0_tar1_del, barc0_tar1_ins, barc0_tar1_switch, barc0_tar1_swap, \
bar1_tar0_del, bar1_tar0_ins, bar1_tar0_switch, bar1_tar0_swap, \
barc1_tar0_del, barc1_tar0_ins, barc1_tar0_switch, barc1_tar0_swap, \
bar1_tar1_del, bar1_tar1_ins, bar1_tar1_switch, bar1_tar1_swap, \
barc1_tar1_del, barc1_tar1_ins, barc1_tar1_switch, barc1_tar1_swap, \
bar0_tar0_del_rand, bar0_tar0_ins_rand, bar0_tar0_switch_rand, bar0_tar0_swap_rand, \
barc0_tar0_del_rand, barc0_tar0_ins_rand, barc0_tar0_switch_rand, barc0_tar0_swap_rand, \
bar0_tar1_del_rand, bar0_tar1_ins_rand, bar0_tar1_switch_rand, bar0_tar1_swap_rand, \
barc0_tar1_del_rand, barc0_tar1_ins_rand, barc0_tar1_switch_rand, barc0_tar1_swap_rand, \
bar1_tar0_del_rand, bar1_tar0_ins_rand, bar1_tar0_switch_rand, bar1_tar0_swap_rand, \
barc1_tar0_del_rand, barc1_tar0_ins_rand, barc1_tar0_switch_rand, barc1_tar0_swap_rand, \
bar1_tar1_del_rand, bar1_tar1_ins_rand, bar1_tar1_switch_rand, bar1_tar1_swap_rand, \
barc1_tar1_del_rand, barc1_tar1_ins_rand, barc1_tar1_switch_rand, barc1_tar1_swap_rand, \
bar0_del_tar0_del_rand, bar0_ins_tar0_del_rand, bar0_switch_tar0_del_rand, bar0_swap_tar0_del_rand, \
barc0_del_tar0_del_rand, barc0_ins_tar0_del_rand, barc0_switch_tar0_del_rand, barc0_swap_tar0_del_rand, \
bar0_del_tar0_ins_rand, bar0_ins_tar0_ins_rand, bar0_switch_tar0_ins_rand, bar0_swap_tar0_ins_rand, \
barc0_del_tar0_ins_rand, barc0_ins_tar0_ins_rand, barc0_switch_tar0_ins_rand, barc0_swap_tar0_ins_rand, \
bar0_del_tar0_switch_rand, bar0_ins_tar0_switch_rand, bar0_switch_tar0_switch_rand, bar0_swap_tar0_switch_rand, \
barc0_del_tar0_switch_rand, barc0_ins_tar0_switch_rand, barc0_switch_tar0_switch_rand, barc0_swap_tar0_switch_rand, \
bar0_del_tar0_swap_rand, bar0_ins_tar0_swap_rand, bar0_switch_tar0_swap_rand, bar0_swap_tar0_swap_rand, \
barc0_del_tar0_swap_rand, barc0_ins_tar0_swap_rand, barc0_switch_tar0_swap_rand, barc0_swap_tar0_swap_rand, \
bar1_del_tar0_del_rand, bar1_ins_tar0_del_rand, bar1_switch_tar0_del_rand, bar1_swap_tar0_del_rand, \
barc1_del_tar0_del_rand, barc1_ins_tar0_del_rand, barc1_switch_tar0_del_rand, barc1_swap_tar0_del_rand, \
bar1_del_tar0_ins_rand, bar1_ins_tar0_ins_rand, bar1_switch_tar0_ins_rand, bar1_swap_tar0_ins_rand, \
barc1_del_tar0_ins_rand, barc1_ins_tar0_ins_rand, barc1_switch_tar0_ins_rand, barc1_swap_tar0_ins_rand, \
bar1_del_tar0_switch_rand, bar1_ins_tar0_switch_rand, bar1_switch_tar0_switch_rand, bar1_swap_tar0_switch_rand, \
barc1_del_tar0_switch_rand, barc1_ins_tar0_switch_rand, barc1_switch_tar0_switch_rand, barc1_swap_tar0_switch_rand, \
bar1_del_tar0_swap_rand, bar1_ins_tar0_swap_rand, bar1_switch_tar0_swap_rand, bar1_swap_tar0_swap_rand, \
barc1_del_tar0_swap_rand, barc1_ins_tar0_swap_rand, barc1_switch_tar0_swap_rand, barc1_swap_tar0_swap_rand, \
bar0_sticky_del, bar0_sticky_ins, bar0_sticky_switch, bar0_sticky_swap]

error_seq_list = [mixed_bar, mixed_barc, mixed_bar_tar, mixed_barc_tar, \
no_bar, unID_subseq, unID_endseq, short, unidentified_seq] 

# Perfect test cases
o_bar0 = testseq('barcode 0', \
    sticky_end + barcode + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_barc0 = testseq('barcode 0 comp', \
    sticky_end + barcodec + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_bar1 = testseq('barcode 1', \
    sticky_end + barcode1 + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_barc1 = testseq('barcode 1 comp', \
    sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_bar0bar1 = testseq('barcode 0 barcode 1', \
    sticky_end + barcode + sticky_end + barcode1 + sticky_end + barcode + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_barc0barc1 = testseq('barcode comp 0 barcode 1 comp', \
    sticky_end + barcodec1 + sticky_end + barcodec + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 0)

# Approximate cases
# Approximate cases: Barcode 0 only, 1 barcode error
o_bar0_del = testseq('barcode 0 deletion', \
    sticky_end + barcode_del + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_bar0_ins = testseq('barcode 0 insertion', \
    sticky_end + barcode_ins + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_bar0_switch = testseq('barcode 0 switch', \
    sticky_end + barcode_switch + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_bar0_swap = testseq('barcode 0 swap', \
    sticky_end + barcode_swap + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_barc0_del = testseq('barcode 0 comp deletion', \
    sticky_end + barcodec_del + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_barc0_ins = testseq('barcode 0 comp insertion', \
    sticky_end + barcodec_ins + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_barc0_switch = testseq('barcode 0 comp switch', \
    sticky_end + barcodec_switch + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_barc0_swap = testseq('barcode 0 comp swap', \
    sticky_end + barcodec_swap + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([3,0]), 0, False, 0)

# Approximate cases: Barcode 1 only, 1 barcode error
o_bar1_del = testseq('barcode 1 deletion', \
    sticky_end + barcode1_del + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_bar1_ins = testseq('barcode 1 insertion', \
    sticky_end + barcode1_ins + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_bar1_switch = testseq('barcode 1 switch', \
    sticky_end + barcode1_switch + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_bar1_swap = testseq('barcode 1 swap', \
    sticky_end + barcode1_swap + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_barc1_del = testseq('barcode 1 comp deletion', \
    sticky_end + barcodec1_del + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_barc1_ins = testseq('barcode 1 comp insertion', \
    sticky_end + barcodec1_ins + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_barc1_switch = testseq('barcode 1 comp switch', \
    sticky_end + barcodec1_switch + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_barc1_swap = testseq('barcode 1 comp swap', \
    sticky_end + barcodec1_swap + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

# Approximate cases: Barcode 0 only, 1 random barcode error
o_bar0_del_rand = testseq('barcode 0 deletion rand', \
    sticky_end + rand_barcode_del + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_bar0_ins_rand = testseq('barcode 0 insertion rand', \
    sticky_end + rand_barcode_ins + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_bar0_switch_rand = testseq('barcode 0 switch rand', \
    sticky_end + rand_barcode_switch + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_bar0_swap_rand = testseq('barcode 0 swap rand', \
    sticky_end + rand_barcode_swap + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_barc0_del_rand = testseq('barcode 0 comp deletion rand', \
    sticky_end + rand_barcodec_del + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_barc0_ins_rand = testseq('barcode 0 comp insertion rand', \
    sticky_end + rand_barcodec_ins + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_barc0_switch_rand = testseq('barcode 0 comp switch rand', \
    sticky_end + rand_barcodec_switch + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_barc0_swap_rand = testseq('barcode 0 comp swap rand', \
    sticky_end + rand_barcodec_swap + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([3,0]), 0, False, 0)

# Approximate cases: Barcode 1 only, 1 random barcode error
o_bar1_del_rand = testseq('barcode 1 deletion rand', \
    sticky_end + rand_barcode1_del + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_bar1_ins_rand = testseq('barcode 1 insertion rand', \
    sticky_end + rand_barcode1_ins + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_bar1_switch_rand = testseq('barcode 1 switch rand', \
    sticky_end + rand_barcode1_switch + sticky_end + barcode1 + sticky_end + barcode1+ sticky_end, \
    np.array([0,3]), 0, False, 1)

o_bar1_swap_rand = testseq('barcode 1 swap rand', \
    sticky_end + rand_barcode1_swap + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_barc1_del_rand = testseq('barcode 1 comp deletion rand', \
    sticky_end + rand_barcodec1_del + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_barc1_ins_rand = testseq('barcode 1 comp insertion rand', \
    sticky_end + rand_barcodec1_ins + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_barc1_switch_rand = testseq('barcode 1 comp switch rand', \
    sticky_end + rand_barcodec1_switch + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

o_barc1_swap_rand = testseq('barcode 1 comp swap rand', \
    sticky_end + rand_barcodec1_swap + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([0,3]), 0, False, 1)

# Approximate cases: Barcode 0 + Barcode 1 errors
o_bar0_bar1_del = testseq('barcode 0 barcode 1 deletion', \
    sticky_end + barcode1_del + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_bar0_bar1_ins = testseq('barcode 0 barcode 1 insertion', \
    sticky_end + barcode1_ins + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_bar0_bar1_switch = testseq('barcode 0 barcode 1 switch', \
    sticky_end + barcode1_switch + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_bar0_bar1_swap = testseq('barcode 0 barcode 1 swap', \
    sticky_end + barcode1_swap + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_barc0_barc1_del = testseq('barcode 0 comp barcode 1 deletion', \
    sticky_end + barcodec1_del + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_barc0_barc1_ins = testseq('barcode 0 comp barcode 1 insertion', \
    sticky_end + barcodec1_ins + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_barc0_barc1_switch = testseq('barcode 0 comp barcode 1 switch', \
    sticky_end + barcodec1_switch + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_barc0_barc1_swap = testseq('barcode 0 comp barcode 1 swap', \
    sticky_end + barcodec1_swap + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([2,1]), 0, False, 0)

# Approximate cases: Barcode + barcode errors
o_bar1_bar0_del = testseq('barcode 1 barcode 0 deletion', \
    sticky_end + barcode_del + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_bar0_ins = testseq('barcode 1 barcode 0 insertion', \
    sticky_end + barcode_ins + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_bar0_switch = testseq('barcode 1 barcode 0 switch', \
    sticky_end + barcode_switch + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_bar0_swap = testseq('barcode 1 barcode 0 swap', \
    sticky_end + barcode_swap + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)
    
o_barc1_barc0_del = testseq('barcode 1 comp barcode 0 comp deletion', \
    sticky_end + barcodec_del + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_barc0_ins = testseq('barcode 1 comp barcode 0 comp insertion', \
    sticky_end + barcodec_ins + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_barc0_switch = testseq('barcode 1 comp barcode 0 comp switch', \
    sticky_end + barcodec_switch + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_barc0_swap = testseq('barcode 1 comp barcode 0 comp swap', \
    sticky_end + barcodec_swap + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

# Approximation cases: Barcode 0 + barcode random errors
o_bar0_bar1_del_rand = testseq('barcode 0 barcode 1 rand deletion', \
    sticky_end + rand_barcode1_del + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_bar0_bar1_ins_rand = testseq('barcode 0 barcode 1 rand insertion', \
    sticky_end + rand_barcode1_ins + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_bar0_bar1_switch_rand = testseq('barcode 0 barcode 1 rand switch', \
    sticky_end + rand_barcode1_switch + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_bar0_bar1_swap_rand = testseq('barcode 0 barcode 1 rand swap', \
    sticky_end + rand_barcode1_swap + sticky_end + barcode + sticky_end + barcode + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_barc0_barc1_del_rand = testseq('barcode 0 comp barcode 1 comp rand deletion', \
    sticky_end + rand_barcodec1_del + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_barc0_barc1_ins_rand = testseq('barcode 0 comp barcode 1 comp rand insertion', \
    sticky_end + rand_barcodec1_ins + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_barc0_barc1_switch_rand = testseq('barcode 0 comp barcode 1 comp rand switch', \
    sticky_end + rand_barcodec1_switch + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([2,1]), 0, False, 0)

o_barc0_barc1_swap_rand = testseq('barcode 0 comp barcode 1 comp rand swap', \
    sticky_end + rand_barcodec1_swap + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([2,1]), 0, False, 0)

# Approximation cases: Barcode 1 + barcode random errors
o_bar1_bar0_del_rand = testseq('barcode 1 barcode 0 rand deletion', \
    sticky_end + rand_barcode_del + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_bar0_ins_rand = testseq('barcode 1 barcode 0 rand insertion', \
    sticky_end + rand_barcode_ins + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_bar0_switch_rand = testseq('barcode 1 barcode 0 rand switch', \
    sticky_end + rand_barcode_switch + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_bar0_swap_rand = testseq('barcode 1 barcode 0 rand swap', \
    sticky_end + rand_barcode_swap + sticky_end + barcode1 + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_barc0_del_rand = testseq('barcode 1 comp barcode comp 0 rand deletion', \
    sticky_end + rand_barcodec_del + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_barc0_ins_rand = testseq('barcode 1 comp barcode comp 0 rand insertion', \
    sticky_end + rand_barcodec_ins + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_barc0_switch_rand = testseq('barcode 1 comp barcode comp 0 rand switch', \
    sticky_end + rand_barcodec_switch + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_barc0_swap_rand = testseq('barcode 1 comp barcode comp 0 rand swap', \
    sticky_end + rand_barcodec_swap + sticky_end + barcodec1 + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)


# Approximate cases: Barcode 1 random error + Random barcode 0 deletion
o_bar1_del_bar0_del_rand = testseq('barcode 1 deletion barcode 0 rand deletion', \
    sticky_end + rand_barcode_del + sticky_end + rand_barcode1_del + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_ins_bar0_del_rand = testseq('barcode 1 insertion barcode 0 rand deletion', \
    sticky_end + rand_barcode_del + sticky_end + rand_barcode1_ins + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_switch_bar0_del_rand = testseq('barcode 1 switch barcode 0 rand deletion', \
    sticky_end + rand_barcode_del + sticky_end + rand_barcode1_switch + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_swap_bar0_del_rand = testseq('barcode 1 swap barcode 0 rand deletion', \
    sticky_end + rand_barcode_del + sticky_end + rand_barcode1_swap + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_del_barc0_del_rand = testseq('barcode 1 comp deletion barcode comp 0 rand deletion', \
    sticky_end + rand_barcodec_del + sticky_end + rand_barcodec1_del + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_ins_barc0_del_rand = testseq('barcode 1 comp insertion barcode comp 0 rand deletion', \
    sticky_end + rand_barcodec_del + sticky_end + rand_barcodec1_ins + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_switch_barc0_del_rand = testseq('barcode 1 comp switch barcode comp 0 rand deletion', \
    sticky_end + rand_barcodec_del + sticky_end + rand_barcodec1_switch + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_swap_barc0_del_rand = testseq('barcode 1 comp swap barcode comp 0 rand deletion', \
    sticky_end + rand_barcodec_del + sticky_end + rand_barcodec1_swap + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

# Approximate cases: Barcode 1 random error + Random barcode 0 insertion
o_bar1_del_bar0_ins_rand = testseq('barcode 1 deletion barcode 0 rand insertion', \
    sticky_end + rand_barcode_ins + sticky_end + rand_barcode1_del + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_ins_bar0_ins_rand = testseq('barcode 1 insertion barcode 0 rand insertion', \
    sticky_end + rand_barcode_ins + sticky_end + rand_barcode1_ins + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_switch_bar0_ins_rand = testseq('barcode 1 switch barcode 0 rand insertion', \
    sticky_end + rand_barcode_ins + sticky_end + rand_barcode1_switch + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_swap_bar0_ins_rand = testseq('barcode 1 swap barcode 0 rand insertion', \
    sticky_end + rand_barcode_ins + sticky_end + rand_barcode1_swap + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_del_barc0_ins_rand = testseq('barcode 1 comp deletion barcode comp 0 rand insertion', \
    sticky_end + rand_barcodec_ins + sticky_end + rand_barcodec1_del + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_ins_barc0_ins_rand = testseq('barcode 1 comp insertion barcode comp 0 rand insertion', \
    sticky_end + rand_barcodec_ins + sticky_end + rand_barcodec1_ins + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_switch_barc0_ins_rand = testseq('barcode 1 comp switch barcode comp 0 rand insertion', \
    sticky_end + rand_barcodec_ins + sticky_end + rand_barcodec1_switch + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_swap_barc0_ins_rand = testseq('barcode 1 comp swap barcode comp 0 rand insertion', \
    sticky_end + rand_barcodec_ins + sticky_end + rand_barcodec1_swap + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

# Approximate cases: Barcode 1 random error + Random barcode 0 switch
o_bar1_del_bar0_switch_rand = testseq('barcode 1 deletion barcode 0 rand switch', \
    sticky_end + rand_barcode_switch + sticky_end + rand_barcode1_del + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_ins_bar0_switch_rand = testseq('barcode 1 insertion barcode 0 rand switch', \
    sticky_end + rand_barcode_switch + sticky_end + rand_barcode1_ins + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_switch_bar0_switch_rand = testseq('barcode 1 switch barcode 0 rand switch', \
    sticky_end + rand_barcode_switch + sticky_end + rand_barcode1_switch + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_swap_bar0_switch_rand = testseq('barcode 1 swap barcode 0 rand switch', \
    sticky_end + rand_barcode_switch + sticky_end + rand_barcode1_swap + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_del_barc0_switch_rand = testseq('barcode 1 comp deletion barcode comp 0 rand switch', \
    sticky_end + rand_barcodec_switch + sticky_end + rand_barcodec1_del + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_ins_barc0_switch_rand = testseq('barcode 1 comp insertion barcode comp 0 rand switch', \
    sticky_end + rand_barcodec_switch + sticky_end + rand_barcodec1_ins + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_switch_barc0_switch_rand = testseq('barcode 1 comp switch barcode comp 0 rand switch', \
    sticky_end + rand_barcodec_switch + sticky_end + rand_barcodec1_switch + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_swap_barc0_switch_rand = testseq('barcode 1 comp swap barcode comp 0 rand switch', \
    sticky_end + rand_barcodec_switch + sticky_end + rand_barcodec1_swap + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

# Approximate cases: Barcode 1 random error + Random barcode 0 swap
o_bar1_del_bar0_swap_rand = testseq('barcode 1 deletion barcode 0 rand swap', \
    sticky_end + rand_barcode_swap + sticky_end + rand_barcode1_del + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_ins_bar0_swap_rand = testseq('barcode 1 insertion barcode 0 rand swap', \
    sticky_end + rand_barcode_swap + sticky_end + rand_barcode1_ins + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_switch_bar0_swap_rand = testseq('barcode 1 switch barcode 0 rand swap', \
    sticky_end + rand_barcode_swap + sticky_end + rand_barcode1_switch + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_bar1_swap_bar0_swap_rand = testseq('barcode 1 swap barcode 0 rand swap', \
    sticky_end + rand_barcode_swap + sticky_end + rand_barcode1_swap + sticky_end + barcode1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_del_barc0_swap_rand = testseq('barcode 1 comp deletion barcode comp 0 rand swap', \
    sticky_end + rand_barcodec_swap + sticky_end + rand_barcodec1_del + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_ins_barc0_swap_rand = testseq('barcode 1 comp insertion barcode comp 0 rand swap', \
    sticky_end + rand_barcodec_swap + sticky_end + rand_barcodec1_ins + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_switch_barc0_swap_rand = testseq('barcode 1 comp switch barcode comp 0 rand swap', \
    sticky_end + rand_barcodec_swap + sticky_end + rand_barcodec1_switch + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

o_barc1_swap_barc0_swap_rand = testseq('barcode 1 comp swap barcode comp 0 rand swap', \
    sticky_end + rand_barcodec_swap + sticky_end + rand_barcodec1_swap + sticky_end + barcodec1 + sticky_end, \
    np.array([1,2]), 0, False, 1)

# Approximation cases: Sticky end errors
o_bar0_sticky_del = testseq('barcode 0 sticky end deletion', \
    sticky_end + barcodec + sticky_del + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_bar0_sticky_ins = testseq('barcode 0 sticky end insertion', \
    sticky_end + barcodec + sticky_ins + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_bar0_sticky_switch = testseq('barcode 0 sticky end switch', \
    sticky_end + barcodec + sticky_switch + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([3,0]), 0, False, 0)

o_bar0_sticky_swap = testseq('barcode 0 sticky end swap', \
    sticky_end + barcodec + sticky_swap + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([3,0]), 0, False, 0)

# Error cases
o_unID_subseq = testseq('unidentified subsequence', \
    sticky_end + 'whatisthis' + sticky_end + barcodec + sticky_end + barcodec + sticky_end, \
    np.array([2,0]), 1, True, 0)

o_unID_endseq = testseq('unidentified end subsequence', \
    sticky_end + barcodec + sticky_end + barcodec + sticky_end + barcodec + sticky_end \
    + 'whatisthis', np.array([3,0]), 1, True, 0)

o_short = testseq('short', 'whatisthis', np.zeros(2), 1, True, -2)

o_unidentified_seq = testseq('unidentified sequence', \
    'whatisthiswhatisthiswhatisthiswhatisthiswhatisthis', np.array([0,0]), 1, True, -2)

o_perfect_seq_list = [o_bar0, o_barc0, o_bar1, o_barc1, o_bar0bar1, o_barc0barc1]

o_approx_seq_list = [o_bar0_del, o_bar0_ins, o_bar0_switch, o_bar0_swap, \
o_barc0_del, o_barc0_ins, o_barc0_switch, o_barc0_swap, 
o_bar1_del, o_bar1_ins, o_bar1_switch, o_bar1_swap, \
o_barc1_del, o_barc1_ins, o_barc1_switch, o_barc1_swap, \
o_bar0_del_rand, o_bar0_ins_rand, o_bar0_switch_rand, o_bar0_swap_rand, \
o_barc0_del_rand, o_barc0_ins_rand, o_barc0_switch_rand, o_barc0_swap_rand, \
o_bar1_del_rand, o_bar1_ins_rand, o_bar1_switch_rand, o_bar1_swap_rand, \
o_barc1_del_rand, o_barc1_ins_rand, o_barc1_switch_rand, o_barc1_swap_rand, \
o_bar0_bar1_del, o_bar0_bar1_ins, o_bar0_bar1_switch, o_bar0_bar1_swap, \
o_barc0_barc1_del, o_barc0_barc1_ins, o_barc0_barc1_switch, o_barc0_barc1_swap, \
o_bar1_bar0_del, o_bar1_bar0_ins, o_bar1_bar0_switch, o_bar1_bar0_swap, \
o_barc1_barc0_del, o_barc1_barc0_ins, o_barc1_barc0_switch, o_barc1_barc0_swap, \
o_bar0_bar1_del_rand, o_bar0_bar1_ins_rand, o_bar0_bar1_switch_rand, o_bar0_bar1_swap_rand, \
o_barc0_barc1_del_rand, o_barc0_barc1_ins_rand, o_barc0_barc1_switch_rand, o_barc0_barc1_swap_rand, \
o_bar1_bar0_del_rand, o_bar1_bar0_ins_rand, o_bar1_bar0_switch_rand, o_bar1_bar0_swap_rand, \
o_barc1_barc0_del_rand, o_barc1_barc0_ins_rand, o_barc1_barc0_switch_rand, o_barc1_barc0_swap_rand, \
o_bar1_del_bar0_del_rand, o_bar1_ins_bar0_del_rand, o_bar1_switch_bar0_del_rand, o_bar1_swap_bar0_del_rand, \
o_barc1_del_barc0_del_rand, o_barc1_ins_barc0_del_rand, o_barc1_switch_barc0_del_rand, o_barc1_swap_barc0_del_rand, \
o_bar1_del_bar0_ins_rand, o_bar1_ins_bar0_ins_rand, o_bar1_switch_bar0_ins_rand, o_bar1_swap_bar0_ins_rand, \
o_barc1_del_barc0_ins_rand, o_barc1_ins_barc0_ins_rand, o_barc1_switch_barc0_ins_rand, o_barc1_swap_barc0_ins_rand, \
o_bar1_del_bar0_switch_rand, o_bar1_ins_bar0_switch_rand, o_bar1_switch_bar0_switch_rand, o_bar1_swap_bar0_switch_rand, \
o_barc1_del_barc0_switch_rand, o_barc1_ins_barc0_switch_rand, o_barc1_switch_barc0_switch_rand, o_barc1_swap_barc0_switch_rand, \
o_bar1_del_bar0_swap_rand, o_bar1_ins_bar0_swap_rand, o_bar1_switch_bar0_swap_rand, o_bar1_swap_bar0_swap_rand, \
o_barc1_del_barc0_swap_rand, o_barc1_ins_barc0_swap_rand, o_barc1_switch_barc0_swap_rand, o_barc1_swap_barc0_swap_rand, \
o_bar0_sticky_del, o_bar0_sticky_ins, o_bar0_sticky_switch, o_bar0_sticky_swap] 

o_error_seq_list = [o_unID_subseq, o_unID_endseq, o_short]

# Run tests for analyze_barcoded_seq
def run_tests(test_type, testseq_list, tol, barcoded=True):
    tests_passed = 0
    for testseq in testseq_list:
        test_name = testseq.test_name
        sequence = testseq.sequence
        true_count = testseq.true_count
        true_error = testseq.true_error
        has_error = testseq.has_error
        true_ID = testseq.true_ID
        answer = (true_count, true_error, has_error, true_ID)
        if barcoded:
            tests_passed += test_barcoded_sequence_print(sequence, test_name, tols, answer)
        else:
            tests_passed += test_sequence_print(sequence, test_name, barcode_tols, answer)
    print('Passed {}/{} {} tests'.format(tests_passed, \
        len(testseq_list), test_type))

# run_tests('perfect', perfect_seq_list, 0.9)
# run_tests('approx', approx_seq_list, 0.9)
# run_tests('error', error_seq_list, 0.9)

# run_tests('barcode only perfect', o_perfect_seq_list, 0.9, barcoded=False)
# run_tests('barcode only approx', o_approx_seq_list, 0.9, barcoded=False)
# run_tests('barcode only error', o_error_seq_list, 0.9, barcoded=False)

print(barcode)
print(barcodec)
print(barcode1)
print(barcodec1)
# Run a specifict test for analyze_barcoded_seq
def run_test(testseq, tol, barcoded=True):
    test_name = testseq.test_name
    sequence = testseq.sequence
    true_count = testseq.true_count
    true_error = testseq.true_error
    has_error = testseq.has_error
    true_ID = testseq.true_ID
    answer = (true_count, true_error, has_error, true_ID)
    if barcoded:
        test_barcoded_sequence_print(sequence, test_name, tols, answer)
    else:
        test_sequence_print(sequence, test_name, barcode_tols, answer)

# run_test(unidentified_seq, 0.9)

# # Test make_dictionary
# # No errors, should work fine
# test_dict, test_sticky_dict, test_seq_len, test_sticky_seq_len, test_num_seq = \
# make_dictionary('test_barcodes.txt', sticky_end)
# # Empty last line, should work fine
# test_dict, test_sticky_dict, test_seq_len, test_sticky_seq_len, test_sticky_len = \
# make_dictionary('test_barcodes_empty_end.txt', sticky_end)
# # Empty first line, should throw empty line error
# test_dict, test_sticky_dict, test_seq_len, test_sticky_seq_len, test_sticky_len = \
# make_dictionary('test_barcodes_empty_begin.txt', sticky_end)
# # Empty middle line, should throw empty line error
# test_dict, test_sticky_dict, test_seq_len, test_sticky_seq_len, test_sticky_len = \
# make_dictionary('test_barcodes_empty_line.txt', sticky_end)
# # Difference in line length, should throw value error
# test_dict, test_sticky_dict, test_seq_len, test_sticky_seq_len, test_num_seq = \
# make_dictionary('test_barcodes_wrong_len.txt', sticky_end)
# print('True dict')
# print(barcode_dict)
# print('Test dict')
# print(test_dict)
# print('True sticky dict')
# print(sticky_barcode_dict)
# print('Test sticky dict')
# print(test_sticky_dict)
# print(test_seq_len)
# print(test_sticky_seq_len)
# print(test_num_seq)

# with open('approx.fastq', 'w') as file:
#     for i in range(len(approx_seq_list)):
#         testseq = approx_seq_list[i]
#         file.write('@{}\n'.format(testseq.test_name))
#         file.write('{}\n'.format(testseq.sequence))
#         file.write('+\n')
#         file.write('qualityscore{}\n'.format(i))

# with open('test_perturb.txt', 'w') as file:
#     file.write('1e-6\t0\n')
#     file.write('0\t1.5e-6\n')


# # Check Levenshtein distances
# barcode_del_lv_b = 'ATGAGGACGAATCTCCCGCTTATAd'
# barcode_del_lv_b2 = 'TGAGGACGAATCTCCCGCTTATAdd'
# barcode_del_lv = 'TATGAGGACGATCTCCCGCTTATAd'
# barcode_ins_lv_b = 'iTATGAGGACGAATCTCCCGCTTAT'
# barcode_ins_lv = 'TATGAGGACGAATCiTCCCGCTTAT'
# barcode_swap_lv_b = 'ATTGAGGACGAATCTCCCGCTTATA'
# barcode_swap_lv = 'TATGAGGACGAACTTCCCGCTTATA'

# sticky_barcode_del_lv_b = 'GCATATGAGGACGAATCTCCCGCTTATATGCAd'
# sticky_barcode_del_lv = 'TGCATATGAGGACGATCTCCCGCTTATATGCAd'
# sticky_barcode_ins_lv_b = 'iTGCATATGAGGACGAATCTCCCGCTTATATGC'
# sticky_barcode_ins_lv = 'TGCATATGAGGACGAATCiTCCCGCTTATATGC'
# sticky_barcode_swap_lv_b = 'TGCAATTGAGGACGAATCTCCCGCTTATATGCA'
# sticky_barcode_swap_lv = 'TGCATATGAGGACGAACTTCCCGCTTATATGCA'

# target_del_lv_b = 'TGGTCGAATCAAGGGGAGGd'
# target_del_lv = 'ATGGTCGAATAAGGGGAGGd'
# target_ins_lv_b = 'iATGGTCGAATCAAGGGGAG'
# target_ins_lv = 'ATGGTCGAATCiAAGGGGAG'
# target_swap_lv_b = 'TAGGTCGAATCAAGGGGAGG'
# target_swap_lv = 'ATGGTCGAACTAAGGGGAGG'
# sticky_target_del_lv_b = 'TGCATGGTCGAATCAAGGGGAGGTGCAd'
# sticky_target_del_lv = 'TGCAATGGTCGAATAAGGGGAGGdTGCA'
# sticky_target_ins_lv_b = 'iTGCAATGGTCGAATCAAGGGGAGGTGC'
# sticky_target_ins_lv = 'TGCAATGGTCGAATCiAAGGGGAGGTGC'
# sticky_target_swap_lv_b = 'TGCATAGGTCGAATCAAGGGGAGGTGCA'
# sticky_target_swap_lv = 'TGCAATGGTCGAACTAAGGGGAGGTGCA'

# sticky_barcode = 'TGCATATGAGGACGAATCTCCCGCTTATATGCA'
# sticky_target = 'TGCAATGGTCGAATCAAGGGGAGGTGCA'

# barcode_del_lv_b_score = 1-jf.levenshtein_distance(barcode_del_lv_b, barcode)/barcode_len
# barcode_del_lv_b2_score = 1-jf.levenshtein_distance(barcode_del_lv_b2, barcode)/barcode_len
# barcode_del_lv_score = 1-jf.levenshtein_distance(barcode_del_lv, barcode)/barcode_len
# barcode_ins_lv_b_score = 1-jf.levenshtein_distance(barcode_ins_lv_b, barcode)/barcode_len
# barcode_ins_lv_score = 1-jf.levenshtein_distance(barcode_ins_lv, barcode)/barcode_len
# barcode_swap_lv_b_score = 1-jf.levenshtein_distance(barcode_swap_lv_b, barcode)/barcode_len
# barcode_swap_lv_score = 1-jf.levenshtein_distance(barcode_swap_lv, barcode)/barcode_len
# sticky_barcode_del_lv_b_score = 1-jf.levenshtein_distance(sticky_barcode_del_lv_b, sticky_barcode)/sticky_barcode_len
# sticky_barcode_del_lv_score = 1-jf.levenshtein_distance(sticky_barcode_del_lv, sticky_barcode)/sticky_barcode_len
# sticky_barcode_ins_lv_b_score = 1-jf.levenshtein_distance(sticky_barcode_ins_lv_b, sticky_barcode)/sticky_barcode_len
# sticky_barcode_ins_lv_score = 1-jf.levenshtein_distance(sticky_barcode_ins_lv, sticky_barcode)/sticky_barcode_len
# sticky_barcode_swap_lv_b_score = 1-jf.levenshtein_distance(sticky_barcode_swap_lv_b, sticky_barcode)/sticky_barcode_len
# sticky_barcode_swap_lv_score = 1-jf.levenshtein_distance(sticky_barcode_swap_lv, sticky_barcode)/sticky_barcode_len
# target_del_lv_b_score = 1-jf.levenshtein_distance(target_del_lv_b, target)/target_len
# target_del_lv_score = 1-jf.levenshtein_distance(target_del_lv, target)/target_len
# target_ins_lv_b_score = 1-jf.levenshtein_distance(target_ins_lv_b, target)/target_len
# target_ins_lv_score = 1-jf.levenshtein_distance(target_ins_lv, target)/target_len
# target_swap_lv_b_score = 1-jf.levenshtein_distance(target_swap_lv_b, target)/target_len
# target_swap_lv_score = 1-jf.levenshtein_distance(target_swap_lv, target)/target_len
# sticky_target_del_lv_b_score = 1-jf.levenshtein_distance(sticky_target_del_lv_b, sticky_target)/sticky_target_len
# sticky_target_del_lv_score = 1-jf.levenshtein_distance(sticky_target_del_lv, sticky_target)/sticky_target_len
# sticky_target_ins_lv_b_score = 1-jf.levenshtein_distance(sticky_target_ins_lv_b, sticky_target)/sticky_target_len
# sticky_target_ins_lv_score = 1-jf.levenshtein_distance(sticky_target_ins_lv, sticky_target)/sticky_target_len
# sticky_target_swap_lv_b_score = 1-jf.levenshtein_distance(sticky_target_swap_lv_b, sticky_target)/sticky_target_len
# sticky_target_swap_lv_score = 1-jf.levenshtein_distance(sticky_target_swap_lv, sticky_target)/sticky_target_len

# print('barcode')
# print(f'del b score:{barcode_del_lv_b_score}')
# print(f'del b2 score:{barcode_del_lv_b2_score}')
# print(f'del score:{barcode_del_lv_score}')
# print(f'ins b score:{barcode_ins_lv_b_score}')
# print(f'ins score:{barcode_ins_lv_score}')
# print(f'swap b score:{barcode_swap_lv_b_score}')
# print(f'swap score:{barcode_swap_lv_score}')
# print('sticky_barcode')
# print(f'del b score:{sticky_barcode_del_lv_b_score}')
# print(f'del score:{sticky_barcode_del_lv_score}')
# print(f'ins b score:{sticky_barcode_ins_lv_b_score}')
# print(f'ins score:{sticky_barcode_ins_lv_score}')
# print(f'swap b score:{sticky_barcode_swap_lv_b_score}')
# print(f'swap score:{sticky_barcode_swap_lv_score}')
# print('target')
# print(f'del b score:{target_del_lv_b_score}')
# print(f'del score:{target_del_lv_score}')
# print(f'ins b score:{target_ins_lv_b_score}')
# print(f'ins score:{target_ins_lv_score}')
# print(f'swap b score:{target_swap_lv_b_score}')
# print(f'swap score:{target_swap_lv_score}')
# print('sticky_target')
# print(f'del b score:{sticky_target_del_lv_b_score}')
# print(f'del score:{sticky_target_del_lv_score}')
# print(f'ins b score:{sticky_target_ins_lv_b_score}')
# print(f'ins score:{sticky_target_ins_lv_score}')
# print(f'swap b score:{sticky_target_swap_lv_b_score}')
# print(f'swap score:{sticky_target_swap_lv_score}')

# print(1-1/barcode_len)
# print(1-2/barcode_len)
# print(1-3/barcode_len)
# print(1-4/barcode_len)

# test_tol = 2
# print(f'target tol: {calc_tol(test_tol,target_len)}')
# print(f'barcode tol: {calc_tol(test_tol,barcode_len)}')
# print(f'sticky target tol: {calc_tol(test_tol,sticky_target_len)}')
# print(f'sticky barcode tol: {calc_tol(test_tol,sticky_barcode_len)}')

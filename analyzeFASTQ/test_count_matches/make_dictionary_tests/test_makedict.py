from PertQuant.analyzeFASTQ.countmatches import make_dictionary

one_line_file = "test_barcodes_one_line.txt"
one_line_dict_results = make_dictionary(one_line_file, is_sticky=False)
print(one_line_dict_results[0])

one_line_newline_file = "test_barcodes_one_line_newline.txt"
one_line_newline_dict_results = make_dictionary(one_line_newline_file, is_sticky=False)
print(one_line_newline_dict_results[0])

multi_line_file = "test_barcodes_multi_line.txt"
multi_line_dict_results = make_dictionary(multi_line_file, is_sticky=False)
print(multi_line_dict_results[0])

multi_line_newline_file = "test_barcodes_empty_end.txt"
multi_line_newline_dict_results = make_dictionary(multi_line_newline_file, is_sticky=False)
print(multi_line_newline_dict_results[0])
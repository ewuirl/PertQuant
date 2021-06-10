from ont_fast5_api.fast5_interface import get_fast5_file

def print_all_raw_data(fast5_filepath):
    with get_fast5_file(fast5_filepath, mode="r") as f5:
        reads = f5.get_reads()
        # raw_data = reads[0].get_raw_data()
        # print(read.read_id, raw_data)
        for read in f5.get_reads():
            raw_data = read.get_raw_data()
            print(read.read_id, raw_data)

fast5_filepath = "/Users/emilywu/OneDrive - Massachusetts Institute of Technology/Minion/Barcode0ML/Barcode0/20210301_1939_MC-110826_0_AFO090_cb0e839b/fast5_pass/AFO090_pass_45c91af1_0.fast5" # This can be a single- or multi-read file
print_all_raw_data(fast5_filepath)
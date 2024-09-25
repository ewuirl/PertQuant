from PertQuant.simCRN.ml_nupack import generate_strand
from PertQuant.simCRN.ml_nupack import gen_complement
from PertQuant.analyzeFASTQ.countmatches import check_list_min_unique_len

length = 25
GC_content = 0.5
min_unique_len = 5
sticky_end = "TGCA"
RE_site = "CCTGCAGG"
barcodeA = "CAATTCATCCATATTGCACCGTGAG"
barcodeAc = gen_complement(barcodeA, comp=1)



i = -1
for GC_content in [0.30, 0.50, 0.70]:
    found_seq = False 
    while not found_seq:
        i += 1
        print(f"Testing strand {i}")
        strand = generate_strand(length, GC_content=GC_content)
        strand_comp = gen_complement(strand, comp=1)
        # Make sure the RE site is not in the sequence
        if RE_site in strand:
            pass
        else:
            target_list = [barcodeA, barcodeAc, strand, strand_comp]
            # Check the minimum unique length
            unique, target_target_list, target_barcode_list, barcode_barcode_list \
            = check_list_min_unique_len(target_list, min_unique_len, barcode_list=[], \
            record=False)

            if unique:
                GC_frac = 0
                for base in strand:
                    if base == "G" or base == "C":
                        GC_frac += 1
                    else:
                        pass
                found_seq = True
                print(f"Target GC content: {GC_content}")
                print(f"Actual GC content: {GC_frac/length}")
                print(f"Found strand: {strand}")
                print(f"Complement: {strand_comp}")
                print(f"Reverse Complement: {strand_comp[::-1]}")

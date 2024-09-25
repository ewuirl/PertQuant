# Check if the sticky end is in any of the sequences used in simulation
strand_file = '/Users/emilywu/Dropbox (MIT)/Emily <-> Ashwin/Nupack/20-20-20-asym-AB-AC/20-20-20-asym-AB-AC.in'
sticky_end = 'TGCA'
with open(strand_file, 'r') as file:
    lines = file.readlines()
    for i in range(len(lines)):
        if sticky_end in lines[i]:
            print('Sticky end in line{}'.format(i))
            seq_list = lines[i].split(sticky_end)
            print(seq_list)
        else:
            pass
    print('Done checking strands')

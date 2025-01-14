import numpy as np

def read_array(i, lines):
    array_list = []
    is_array = True
    while is_array and i <len(lines):
        line = lines[i]
        try:
            row = np.array(lines[i].rstrip('\n').split(), dtype='float')
            array_list.append(row)
        except:
            break
        i+=1
    array = np.array(array_list)
    return i, array

def read_eq_data_settings(settings_file, quiet=0):
    settings_dict = {}
    int_list = ['N', 'M', 'L', 'N_runs']
    float_list = ['Cmin', 'Cmax']
    bool_list = ['self']
    with open(settings_file, 'r') as file:
        lines = file.readlines()
        i = 0
        while i < len(lines): 
            line = lines[i]
            # Read in non-array setting
            if '=' in line:
                line_list = line.rstrip('\n').split(" = ")
                key = line_list[0]
                if key in int_list:
                    setting = int(line_list[1])
                elif key in float_list:
                    setting = float(line_list[1])
                elif key in bool_list:
                    setting = bool(line_list[1])
                else:
                    setting = line_list[1]
                settings_dict[key] = setting
                if quiet == 0:
                    print(f'{key} = {setting}')
                i+=1
            # Read in an array
            elif 'array' in line:
                key = line.rstrip('\n')
                i+=1
                i, array = read_array(i, lines)
                settings_dict[key] = array
                if quiet == 0:
                    print(f'{key}')
                    print(f'{array}')
            else:
                i+=1
    
    # Check N, M, L, and N_runs                
    N = settings_dict['N']
    M = settings_dict['M']
    assert N > 0, 'N must be > 0'
    assert M > 0, 'M must be > 0'
    keys = settings_dict.keys()
    if 'L' in keys:
        L = settings_dict['L']
        assert L > 0, 'L must be > 0'
    assert settings_dict['N_runs'] > 0, 'N_runs must be > 0'

    # Check Cmin and Cmax
    assert 0 <= settings_dict['Cmin'] < settings_dict['Cmax'], '0 <= Cmin < Cmax'
    
    # Flatten initial concentration arrays and check sizes
    settings_dict['Ai_array'] = settings_dict['Ai_array'].flatten()
    settings_dict['Bi_array'] = settings_dict['Bi_array'].flatten()
    assert len(settings_dict['Ai_array']) == N, f'Ai_array must be length {N}'
    assert len(settings_dict['Bi_array']) == M, f'Bi_array must be length {M}'

    # Check Karray sizes
    assert settings_dict['KAB_array'].shape == (N,M), f'KAB_array must be shape ({N}, {M})'
    if 'KBC_array' in keys:
        assert settings_dict['KBC_array'].shape == (M,L), f'KBC_array must be shape ({M}, {L})'
    else:
        settings_dict['KBC_array'] = np.zeros((M,L))
    if 'KAC_array' in keys:
        assert settings_dict['KAC_array'].shape == (N,L), f'KAC_array must be shape ({N}, {L})'
    else:
        settings_dict['KAC_array'] = np.zeros((N,L))
    if 'KAA_array' in keys:
        assert settings_dict['KAA_array'].shape == (N,N), f'KAA_array must be shape ({N}, {N})'
    else:
        settings_dict['KAA_array'] = np.zeros((N,N))
    if 'KBB_array' in keys:
        assert settings_dict['KBB_array'].shape == (M,M), f'KBB_array must be shape ({M}, {M})'
    else:
        settings_dict['KBB_array'] = np.zeros((M,M))
    if 'KCC_array' in keys:
        assert settings_dict['KCC_array'].shape == (L,L), f'KCC_array must be shape ({L}, {L})'
    else:
        settings_dict['KCC_array'] = np.zeros((L,L))

    return settings_dict
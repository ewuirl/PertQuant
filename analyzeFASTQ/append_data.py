
with open(file_name, 'a') as file:
    file.write('# N = {:d} \n'.format(N))
    file.write('# M = {:d} \n'.format(M))
    file.write('# L = {:d} \n'.format(L))
    file.write('# N_runs = {:d} \n'.format(N_runs))
    file.write('# Cmin = {} \n'.format(Cmin))
    file.write('# Cmax = {} \n'.format(Cmax))
    file.write('# Ai = {} \n'.format(Ai))

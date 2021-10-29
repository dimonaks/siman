#!/usr/bin/env python3
""" 
Include:
1. runBash(cmd)
2. CalcResults
3. interstitial()
4. out_for_paper()
5. shift_analys(st)
6. write_geo(st)
"""

import subprocess as SP
import shutil as S
from siman.small_functions import makedir
from siman import header


def make_vaspkit_kpoints(type_calc='', type_lattice='', poscar='', list_kpoints=[], kpoints_density=0.02, k_cutoff=0.015, num_points=6, folder_vaspkit=''):

    makedir(folder_vaspkit)

    if type_lattice == 'hex':
        initial_list_kpoints = [((0.000000, 0.000, 0.000, 'G'),(0.33333, 0.33333, 0.000, 'K')),
                                ((0.000000, 0.000, 0.000, 'G'),(0.50000, 0.00000, 0.000, 'M')),
                                ((0.000000, 0.000, 0.000, 'G'),(0.50000, 0.00000, 0.000, 'X')),
                                ((0.000000, 0.000, 0.000, 'G'),(0.00000, 0.50000, 0.000, 'Y'))]
    
    S.copy(poscar, folder_vaspkit+'/POSCAR')
    f = open(folder_vaspkit+'info', 'w')
    f.write(poscar+' = POSCAR')
    f.close()

    if type_calc == 'effective_mass':
        full_list_kpoints = initial_list_kpoints + list_kpoints

        f = open(folder_vaspkit+'/VPKIT.in', 'w')
        f.write('1       # "1" for pre-process (generate KPOINTS), "2" for post-process(calculate m*)\n')
        f.write(str(num_points)+'      #  number of points for quadratic function fitting.\n')
        f.write(str(k_cutoff)+'   # k-cutoff, unit Å-1.\n')
        f.write(str(len(full_list_kpoints))+'       # number of tasks for effective mass calculation\n')
        for i in full_list_kpoints:
            f.write('{0:10.7f} {1:10.7f} {2:10.7f} {3:10.7f} {4:10.7f} {5:10.7f}'.format(i[0][0], i[0][1], i[0][2], 
                                                                                         i[1][0], i[1][1], i[1][2])+3*' '+i[0][3]+'->'+i[1][3]+'\n')
        f.close()

        s = SP.Popen(header.PATH2VASPKIT+' -task 912 -kpr '+str(kpoints_density), cwd=folder_vaspkit, shell=True)
        
        s.wait()

        print('File '+folder_vaspkit+'/KPOINTS was built successfully')

def insert_0weight_kpoints(ibzkpt='', list_kpoints=[], folder_kpoints=''):

    f = open(ibzkpt)
    l = f.readlines()
    f.close()

    l[1] = str(int(l[1])+len(list_kpoints))+'\n'

    l_new = l + list_kpoints

    makedir(folder_kpoints)
    f = open(folder_kpoints+'/KPOINTS', 'w')
    f.writelines(l_new)
    f.close()

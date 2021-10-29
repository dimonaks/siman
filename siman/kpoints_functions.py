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
    """
    The function is used to prepare KPOINTS using the "VASPKIT" package
    ###INPUT:
        
        * type_calc (str)             - quantity, which is calculated 
        * type_lattice (str)          - one of common lattices such as bcc, fcc, hcp etc.
        * poscar (str)                - POSCAR file with structure, for which the KPOINTS file is built
        * list_kpoints (list)         - list of additional kpoints to add into the KPOINTS file with zero weights;
                                        has format [((a1,a2,a3,label_a),(b1,b2,b3,label_b)),...], where "a" and "b" are 
                                        the coordinates of points in the reciprocal space, "label" is the name of the point
                                        These two points determine the direction in the reciprocal space.
        * kpoints_density (float)     - density of k-points with non-zero weights
        * k_cutoff (float)            - distance along the choosen direction in the reciprocal space in A^{-1} unit
        * num_points (int)            - number of k-points along the choosen direction in the reciprocal space
        * folder_vaspkit (str)        - name of folder to make the KPOINTS file

    ###SOURCE

    For more details see:

    https://vaspkit.com/tutorials.html#effective-mass 

    ###RETURN:
        
        None

    """

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
    """
    The function is used to add k-points to the current IBZKPT file 
    and convert it to new KPOINTS file
    ###INPUT:
        
        * ibzkpt (str)             - path to IBZKPT file to insert k-points with zero weight 
        * list_kpoints (list)      - list of k-points with relative coordinates; 
                                     it has the following form ['a1, a2, a3 0\n',...]
        * folder_kpoints (str)     - name of folder to make the resulting KPOINTS file

    ###SOURCE

        None

    ###RETURN:
        
        None

    """
    f = open(ibzkpt)
    l = f.readlines()
    f.close()

    l[1] = str(int(l[1])+len(list_kpoints))+'\n'

    l_new = l + list_kpoints

    makedir(folder_kpoints)
    f = open(folder_kpoints+'/KPOINTS', 'w')
    f.writelines(l_new)
    f.close()

#!/usr/bin/env python3
""" 
Copyright (c) Siman Development Team.
Distributed under the terms of the GNU License.
Author: Kartamyshev A.I. (Darth Feiwante) + Aksyonov D.A
"""
import math
import numpy as np 
import subprocess as SP
import shutil as S
from siman.small_functions import makedir
from siman.geo import calc_recip_vectors
from siman import header




def ngkpt2kspacings(ngkpt, rprimd):
    """
    Calculate kspacing from ngkpt and rprimd (A)

    INPUT:
        - ngkpt (list of 3 int) - k-point grid
        - rprimd (list of lists*3) - three primitive vectors in A


    Author: Aksyonov D.A.
    """
    kspacing = []

    recip = calc_recip_vectors(rprimd)

    for i in 0, 1, 2:
        if ngkpt[i] != 0:
            a = np.linalg.norm( recip[i] ) / ngkpt[i]
            # kspacing.append(red_prec(a))
            kspacing.append(a)
        else:
            printlog('Warning! ngkpt = 0 in siman.geo.calc_kspacings()', imp='y')

    return  kspacing

def kspacing2ngkpt(kspacing, rprimd):
    """
    Calculate number of k-points for the given k-spacing in A^-1; The units align with VASP convention.
    INPUT:
        - kspacing (float) - spacing in reciprocal space in A^-1
        - rprimd (list of lists*3) - three primitive vectors in A

    RETURN:
        - ngkpt (list*3 int)

    TODO:
        remove old small_functions.calc_ngkpt with different interface and use this instead

    Author: Aksyonov D.A.
    """
    to_ang_local = 1
    
    N_from_kspacing = []
    recip = calc_recip_vectors(rprimd)

    for i in 0, 1, 2:
        N_from_kspacing.append( math.ceil( (np.linalg.norm(recip[i]) / to_ang_local) / kspacing) )

    print('Actual k-spacings are {:.2f}, {:.2f}, {:.2f}'.format(*ngkpt2kspacings(N_from_kspacing, rprimd)) )

    return N_from_kspacing



def make_vaspkit_kpoints(type_calc='', type_lattice='', poscar='', list_kpoints=[], kpoints_density=0.02, k_cutoff=0.015, num_points=6, folder_vaspkit=''):
    """
    The function is used to prepare KPOINTS using the "VASPKIT" package
    INPUT:
        
        - type_calc (str)             - quantity, which is calculated 
            - "effective mass" - prepare the KPOINTS file to calculate the effective mass along chosen directions in the reciprocal space
        - type_lattice (str)          - one of common lattices such as bcc, fcc, hcp etc.
            - "hcp" - hexagonal close-packed lattice
        - poscar (str)                - POSCAR file with structure, for which the KPOINTS file is built
        - list_kpoints (list of tuples)         - list of additional kpoints to add into the KPOINTS file with zero weights;
                                                  has format [((a1,a2,a3,label_a),(b1,b2,b3,label_b)),...], where "a" and "b" are 
                                                  the relative coordinates of points (floats) in the reciprocal space, "label" (str) is the name of the point
                                                  These two points determine the direction in the reciprocal space.
        - kpoints_density (float)     - density of k-points with non-zero weights
        - k_cutoff (float)            - distance along the choosen direction in the reciprocal space in A^{-1} unit
        - num_points (int)            - number of k-points along the choosen direction in the reciprocal space
        - folder_vaspkit (str)        - name of folder to make the KPOINTS file
    RETURN:
        None
    SOURCE:
        For more details see:
        https://vaspkit.com/tutorials.html#effective-mass 
    TODO:
        - Add more common lattices such as bcc, fcc etc.
        - Add another types of calculations requiring the specific KPOINTS files

    Author: Kartamyshev A.I.
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
    INPUT:
        
        - ibzkpt (str)             - path to IBZKPT file to insert k-points with zero weight 
        - list_kpoints (list)      - list of k-points with relative coordinates; 
                                     it has the following form ['a1, a2, a3 0\n',...]
        - folder_kpoints (str)     - name of folder to make the resulting KPOINTS file
    RETURN:
        None
    SOURCE
        None
    TODO
        Some improvements

    Author: Kartamyshev A.I.
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

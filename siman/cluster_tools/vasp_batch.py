#!/usr/bin/python
"""
NAME
        vasp_batch.py - create VASP batch script for several INCARs

SYNOPSIS
        make_batch.py path_to_vasp_folder

DESCRIPTION
        Programm creates batch script for sequential run of VASP for several INCARs;

        Provide vasp folder as input argument
        INCARS should be like 1.INCAR, 2.INCAR, 3.INCAR, x.INCAR, ..., n.INCAR; n - is the last one


        The batch script is created in incar.batch folder of vasp folder
        The POSCAR is renamed to 1.POSCAR
        The OUTCARS are renamed to 1.x.OUTCAR
        The OUTCAR obtained with n.INCAR is renamed to 1.OUTCAR

REQUIRE
        numpy, siman

AUTHOR
        Aksyonov Dmitry, Skoltech, Moscow

"""
import sys, glob, os
# print sys.argv
if len(sys.argv) != 2:
    print(__doc__)
    sys.exit()

sys.path.append('/usr/lib64/python2.7/site-packages/numpy')
import numpy as np
sys.path.append('/home/aksenov/Simulation_wrapper/siman2') #path to siman library
sys.path.append('/home/aksenov/tools/siman2') #path to siman library
from set_functions import InputSet
from calc_manage import add_loop
import header





vaspfolder = sys.argv[1]
incars = glob.glob(vaspfolder+'/*INCAR*')
try:
    kpoint = glob.glob(vaspfolder+'/*KPOINTS*')[0]
except:
    kpoint = None

poscar = glob.glob(vaspfolder+'/*POSCAR*')[0]
potcar = glob.glob(vaspfolder+'/*POTCAR*')[0]

# print(incars)

sets = []
for incar in incars:
    ise = os.path.basename(incar).split('.')[0]
    s = InputSet(ise, potcar)
    s.read_incar(incar)
    sets.append(s)
    # s.printme()
def to_int(s):
    try:
        ise = int(s.ise)
    except:
        ise = s.ise
        print('Attention! non-integer prefix of Incar')
    return ise


sets = sorted(sets, key = to_int)
print('Order of INCARS:')
for s in sets:
    print(s.ise)


start_set = sets[0]
start_set.remove_last_wave = False
start_set.kpoints_file = kpoint
start_set.set_sequence = sets[1:]
header.varset['batch'] = start_set

header.final_vasp_clean = False
header.warnings = False
header.copy_to_cluster_flag = False
add_loop('incar', 'batch', 1, input_geo_file = poscar, it_folder = vaspfolder, cluster_home = '' )


print('\n\n\n\nPlease run:   \n   sbatch -p AMG '+vaspfolder+'/incar.batch/incar.batch.run')

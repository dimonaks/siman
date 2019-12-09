#!/bin/python
"""  
Program Polaron hop by Aksyonov Dmitry, Skoltech, Moscow
Multiset is not supported yet
"""

import sys, json, os, glob, copy
print('Python version is', sys.version)
from shutil import copyfile
import random

from os.path import expanduser
home = expanduser("~")
sys.path.append(home+'/tools/') # for numpy libraries
import numpy as np

from siman import header
from siman.monte_functions import metropolis
from siman.header import runBash, printlog
from siman.header import ALKALI_ION_ELEMENTS as AM
from siman.header import TRANSITION_ELEMENTS as TM
from siman.classes import CalculationVasp, Structure
from siman.inout import read_poscar
from siman.functions import invert, update_incar
from siman.analysis import suf_en
from siman.monte import vasp_run
from siman.geo import interpolate
debug2 = 0


def copy_vasp_files(v):
    """
    Tool to treat vasp files after calculation
    v (int) - version
    """
    for file in ['OUTCAR', 'CONTCAR', 'CHGCAR']:
        copyfile(file, str(v)+'.'+file)
        if 'CHGCAR' in file:
            runBash('gzip -f '+str(v)+'.'+file)



if __name__ == "__main__":


    debug = 0

    header.warnings = 'yY'
    # header.warnings = 'neyY'
    header.verbose_log = 1

    printlog('\n\n\nStarting Polaron hop script!')
    
    """0. Read configuration file """
    if os.path.exists('conf.json'):
        with open('conf.json', 'r') as fp:
            params = json.load(fp)
    else:
        printlog('Warning! no configuration file conf.json, exiting')
        sys.exit()
        params = {}

    vasprun_command = params.get('vasp_run') or 'vasp'
    images = params.get('images') or 3 # number of images


    """1. Calculate (relax) initial and final positions """
    copyfile('1.POSCAR', 'POSCAR')
    cl1 = vasp_run(3, 'Start position ')
    copy_vasp_files(v)


    copyfile('2.POSCAR', 'POSCAR')
    cl2 = vasp_run(3, 'End position ')
    copy_vasp_files(v)

    """2. Create intermediate steps"""
    interpolate(cl1.end, cl2.end, images, 3)

    """3. Calculate energies of intermediate steps"""
    update_incar(parameter = 'NSW', value = 0)
    
    for v in range(3, 3+images):
        copyfile(str(v)+'.POSCAR', 'POSCAR')
        vasp_run(3, 'End position ')
        copy_vasp_files(v)


    printlog('PH simulation finished!', imp = 'y')
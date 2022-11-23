#!/bin/python
"""  
Program Polaron hop by Aksyonov Dmitry, Skoltech, Moscow
Multiset is not supported yet
Version naming 
AR - atomic relaxation
SP - single point
1 - polaron and deformation at start position (AR)
2 - polaron and deformation at final position (AR)
21 - polaron and deformation at start position (SP)
21+images - polaron at start position, deformation at final position (SP)
22:22+images - intermediate between start and end (SP) 
42 - polaron and deformation at final position (SP)
42+images - polaron at final position, deformation at start position (SP)
43:43+images -  intermediate deformation between end and start  (SP)


"""

import sys, json, os, glob, copy
print('Python version is', sys.version)
from shutil import copyfile
import random

from os.path import expanduser
home = expanduser("~")
print(home)
sys.path.append(home+'/tools/') # for numpy libraries
sys.path.append(home+'/tools/siman/') # 
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
    for file in ['OUTCAR', 'CONTCAR', 'CHGCAR', 'OSZICAR']:
        copyfile(file, str(v)+'.'+file)
        if 'CHGCAR' in file:
            runBash('gzip -f '+str(v)+'.'+file)


def vasp_step(v, msg, rm = 0):
    printlog('Calculating '+msg+ ' point!\n', imp = 'y')
    copyfile(str(v)+'.POSCAR', 'POSCAR')
    cl = vasp_run(3, msg, vasprun_command = vasprun_command)
    copy_vasp_files(v)
    if rm:
        runBash('rm CHGCAR CHG WAVECAR')
    return cl

if __name__ == "__main__":


    debug = 0

    header.warnings = 'yY'
    # header.warnings = 'neyY'
    header.verbose_log = 1

    printlog('\n\n\nStarting Polaron hop script!\n', imp = 'y')
    
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
    mode   = params.get('mode') or 'inherit' # number of images
    magmom = params.get('magmom') or None


    for i in [1,2]:
        path_i = params.get(f'path{i}')
        for pp in ['charge','output']:
            path_ii = path_i[pp]
            ftype = path_ii.split('.')[-1]              
            runBash(f'cp {home}/{path_ii} {i}.{ftype}')


    cl1 = CalculationVasp(output = '1.OUTCAR')
    cl1.read_results(show = 'fo')

    cl2 = CalculationVasp(output = '2.OUTCAR')
    cl2.read_results(show = 'fo')


    update_incar(parameter = 'NSW', value = 0, run = 1, write = 0)
    runBash('cp 1.CHGCAR  CHGCAR')
    interpolate(cl1.end, cl2.end, images, 21, omit_edges = 0)
    for v in range(21, 21+images):
        vasp_step(v, 'Intermediate position'+str(v), rm = 0 )

    runBash('rm CHGCAR WAVECAR;  cp 2.CHGCAR CHGCAR')

    interpolate(cl2.end, cl1.end, images, 42, omit_edges = 0)
    for v in range(42, 42+images):
        vasp_step(v, 'Intermediate position'+str(v), rm = 0 )




    runBash('rm CHG WAVECAR')
    printlog('PH simulation finished!', imp = 'y')

#!/bin/python
"""  
Program shiftk by Aksyonov Dmitry, Skoltech, Moscow



"""

import sys, json, os, glob, copy
from shutil import copyfile
import random

from os.path import expanduser
home = expanduser("~")
sys.path.append(home+'/tools/') # for siman library

import numpy as np
from siman import header
from siman.header import runBash, printlog
from siman.classes import CalculationVasp, Structure
from siman.inout import read_poscar
from siman.functions import invert
from siman.monte import vasp_run
from siman.kmesh.kshifts import generate_uniform_random_shifts_bcc_sym, equidistant_shifts



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
    print('vasprun_command:', vasprun_command)
    cl = vasp_run(3, msg, vasprun_command = vasprun_command + ' > vasp.log')
    copy_vasp_files(v)
    if rm:
        runBash('rm CHGCAR CHG WAVECAR')
    return cl

def write_energy(shiftk, e):
    with open('ENERGIES', 'a') as f:
        f.write('{:.3f} {:.3f} {:.3f}    {:.9f}\n'.format(shiftk[0], shiftk[1], shiftk[2], e))





if __name__ == "__main__":


    debug = 0

    header.warnings = 'yY'
    # header.warnings = 'neyY'
    header.verbose_log = 1
    printlog('Python version is', sys.version)

    printlog('\n\n\nStarting shiftk script!\n', imp = 'y')
    
    """0. Read configuration file """
    if os.path.exists('shiftk_conf.json'):
        with open('shiftk_conf.json', 'r') as fp:
            params = json.load(fp)
    else:
        printlog('Warning! no configuration file shiftk_conf.json, exiting')
        sys.exit()
        params = {}

    vasprun_command = params.get('vasp_com') or 'vasp_std'
    nk = params.get('nk')# number of k-points
    mk = params.get('mk')# size of grid for shifts
    m = params.get('m')# number of shifts for sobol
    mode   = params.get('shift_type') or 'eqd' # mode type: equidistant, random, provided
    v   = params.get('v') or 1 # version
    chgcar_name   = params.get('chgcar') 
    if chgcar_name is None:
        printlog('Warning! params[\'chgcar\'] is none ', imp = 'y')
    else:
        printlog('However I am looking just for CHGCAR in the current folder as the file was expected to be renamed', imp = 'y')
        chgcar_name = 'CHGCAR'
        if os.path.exists(chgcar_name):
            if not os.path.exists('shiftk_CHGCAR'):
                copyfile(chgcar_name, 'shiftk_CHGCAR') # save it with new name to protect from overwrite
            chgcar_name = 'shiftk_CHGCAR'
        else:
            printlog('Warning! CHGCAR is missing', imp = 'y')


    cl0 = CalculationVasp().deserialize('shiftk.pickle')
    printlog('Selected mode is ', mode, imp = 'y')



    if mode == 'eqd':
        kshifts = equidistant_shifts(mk)
    elif mode == 'sobol':
        kshifts = params.get('kshifts')

    elif mode =='random':
        kshifts = generate_uniform_random_shifts_bcc_sym(
            m=mk,
            min_dist=0.08*3/mk,
            attempt_limit=4000)



    if 0:
        cl0.make_kpoints_file("KPOINTS")
        cl0.set.params['LSORBIT'] = None
        cl0.set.params['ICHARG'] = None
        cl0.make_incar('INCAR')
        cl = vasp_step(v, 'chargcar calc ', 0 )
        cl0.set.params['LSORBIT'] = '.TRUE.'
        cl0.set.params['ICHARG'] = 11
        cl0.make_incar('INCAR')
        chgcar_name = 'shiftk_CHGCAR'
        copyfile('CHGCAR', chgcar_name )

    runBash('rm ENERGIES')

    for shiftk in kshifts:
        
        if chgcar_name:
            copyfile(chgcar_name, 'CHGCAR')
        
        cl0.set.shiftk = shiftk
        print(cl0.set.ngkpt)
        cl0.make_kpoints_file("KPOINTS")
        # sys.exit()
        cl = vasp_step(v, 'k - shift '+str(shiftk), 1 )
        write_energy(shiftk, cl.e0)


    runBash('rm CHG WAVECAR')
    printlog('kshift simulation finished!', imp = 'y')
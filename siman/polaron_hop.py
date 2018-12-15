#!/bin/python
"""  
Program Polaron hop by Aksyonov Dmitry, Skoltech, Moscow

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
from siman.functions import invert
from siman.analysis import suf_en

debug2 = 0


def vasp_run(n, des):
    #allows to run vasp several times, here fireworks can be used
    #n - number of attempts
    #des - description of run
    for i in range(n): # max three attempts
        
        if not debug2:
            out = runBash(vasprun_command)
            printlog(des, 'attempt', i,'out is', out)
        
        cl = CalculationVasp(output = 'OUTCAR')
        out = cl.read_results(show = 'fo')
        printlog('Results are', imp = 'y')
        printlog(out, imp = 'y')

        status = check(cl)
        if status == 0:
            break
        else:
            if os.path.exists('CONTCAR'):
                copyfile('CONTCAR', 'POSCAR')
            else:
                printlog('No CONTCAR was found. No further attempts to run VASP', imp = 'y')
                break

    return cl   



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


    """1. Calculate (relax) initial and final positions """
    copyfile('1.POSCAR', 'POSCAR')

    cl_new = vasp_run(3, 'mcstep '+str(i_mcstep))


    printlog('PH simulation finished!', imp = 'y')
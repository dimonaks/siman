#!/bin/python
"""  
Programm Monte-Carlo
"""

import sys, json, os, glob
print('Python version is', sys.version)
from shutil import copyfile
import random

from os.path import expanduser
home = expanduser("~")
sys.path.append(home+'/tools') # for numpy libraries
import numpy as np

import header
from monte_functions import metropolis
from header import runBash, ALKALI_ION_ELEMENTS, TRANSITION_ELEMENTS, printlog
from classes import CalculationVasp, Structure


debug = 0

def check(cl, exit = 0):
    if hasattr(cl, 'e0'):
        printlog('outcar is ok, continue', imp  = 'y')
        out = 0
    else:
        printlog('outcar is broken ', imp  = 'y')
        out = 1
        if exit:
            printlog('exiting...', imp  = 'y')
            sys.exit()
    return out


header.warnings = 'neyY'


"""0. Read configuration file """
# params = read_monte_params()
with open('monte.json', 'r') as fp:
    params = json.load(fp)
vasp_run = params['vasp_run']
print(vasp_run)

nmcstep = params['mcsteps']

print('Total number of steps is', nmcstep)






"""1. Run initial calculation"""
if debug:
    cl = CalculationVasp()
    cl.read_poscar('1.POSCAR')
    cl.end = cl.init
    lastn = 0
else:
    files = glob.glob('*.pickle') #get list of calculated files
    if files:
        numbers = [int(file.split('.')[0]) for file in files]
        lastn = max(numbers)
        last_file = str(lastn)+'.pickle'
        printlog('Last calculation file is ', last_file)
    else:
        lastn = 0
        last_file = None

    if last_file:
        cl = CalculationVasp().deserialize(last_file)
        printlog('Successfully deserialized')
    else:
        out = runBash(vasp_run)
        print('first run', out)
        cl = CalculationVasp(output = 'OUTCAR')
        cl.read_results()
        check(cl, exit = 1)
        cl.serialize('0')
        with open('ENERGIES', 'w') as f:
            f.write('{:.5f}\n'.format(cl.e0))



st = cl.end



t = params['thickness']
z2 = st.get_surface_pos()[1]


for i_mcstep in range(1+lastn, 1+lastn+nmcstep):

    """3. Exchange two atoms"""
    #choose two atoms 
    alk = st.get_specific_elements(ALKALI_ION_ELEMENTS, z_range = [z2-t, z2])
    tra = st.get_specific_elements(TRANSITION_ELEMENTS, z_range = [z2-t, z2])
    # print(alk, tra)
    st_new_init = st.swap_atoms(random.choice(alk), random.choice(tra))

    if debug:
        st_new_init.write_poscar('POSCAR-'+str(i_mcstep))

    else:
        """4. Write new structure and calculate energy  """
        st_new_init.write_poscar('POSCAR')
    
        out = runBash(vasp_run)
        print('mcstep '+str(i_mcstep), out)

        cl_new = CalculationVasp(output = 'OUTCAR')
        cl_new.read_results()

        if check(cl):
            print('unlucky configuration, trying another ... ')
            with open('ENERGIES', 'a') as f:
                f.write('0\n')
            continue

        """5. Check if to accept new structure  """
        printlog('Energies before and after are ', cl.e0, cl_new.e0, imp = 'y')
        with open('ENERGIES', 'a') as f:
            f.write('{:.5f}\n'.format(cl_new.e0))
        
        if metropolis(cl.e0, cl_new.e0):
            cl = cl_new
            st = cl_new.end
            cl.serialize(str(i_mcstep))
            copyfile('CONTCAR', 'CONTCAR_last')
            copyfile('OUTCAR', 'OUTCAR_last')


print('MC simulation finished!')
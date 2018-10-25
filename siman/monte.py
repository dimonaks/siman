#!/bin/python
"""  
Program Monte-Carlo by Aksyonov Dmitry, Skoltech, Moscow

"""

import sys, json, os, glob
print('Python version is', sys.version)
from shutil import copyfile
import random

from os.path import expanduser
home = expanduser("~")
sys.path.append(home+'/tools/') # for numpy libraries
import numpy as np

from siman import header
from siman.monte_functions import metropolis
from siman.header import runBash, ALKALI_ION_ELEMENTS, TRANSITION_ELEMENTS, printlog
from siman.classes import CalculationVasp, Structure
from siman.inout import read_poscar


def check(cl, exit = 0):
    # return 0 if ok, return 1 if failed
    if hasattr(cl, 'e0'):
        printlog('outcar is ok', imp  = 'y')
        out = 0
    else:
        printlog('outcar is broken ', imp  = 'y')
        out = 1
        if exit:
            printlog('exiting...', imp  = 'y')
            sys.exit()
    return out


def vasp_run(n, des):
    #allows to run vasp several times, here fireworks can be used
    #n - number of attempts
    #des - description of run
    for i in range(n): # max three attempts
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

    """0. Read configuration file """
    # params = read_monte_params()
    if os.path.exists('monte.json'):
        with open('monte.json', 'r') as fp:
            params = json.load(fp)
    else:
        printlog('Warning! no configuration file monte.json, exiting')
        sys.exit()
        params = {}

    vasprun_command = params.get('vasp_run') or 'vasp'
    nmcstep = params.get('mcsteps') or 2 # minimum two steps are done
    thickness = params.get('thickness') or 6 # minimum layer 
    temperature = params.get('temp') or 1




    printlog('\n\n\nStarting Monte-Carlo script!')
    printlog('Command to run vasp', vasprun_command)
    printlog('Total number of steps is', nmcstep)
    printlog('Thickness is ', thickness)
    printlog('Temperature is ', temperature, 'K')






    """1. Run initial calculation"""
    if debug:
        cl = CalculationVasp()
        cl.read_poscar('1.POSCAR')
        cl.end = cl.init
        last_number = 0
    else:
        
        files_yes = glob.glob('*yes.pickle') #get list of calculated yes files
        files_all = glob.glob('*.pickle') #get list of all calculated files
        

        """Find last yes pickle file"""
        if files_yes:
            yes_numbers = [int(file.split('-')[0]) for file in files_yes]
            all_numbers = [int(file.split('-')[0]) for file in files_all]
            last_yes_n = max(yes_numbers)
            last_number = max(all_numbers)
            last_yes_file = str(last_yes_n)+'-yes.pickle'
            printlog('Last calculation file is ', last_yes_file, imp = 'y')
        
        else:
            last_number = 0
            last_yes_file = None

        """Read last pickle file or run vasp """
        
        if last_yes_file:
            cl = CalculationVasp().deserialize(last_yes_file)
            printlog('Successfully deserialized')

            if 0:
                #to correct mistake on first step and add void to pickle
                if 'xvoid' in params and len(params['xvoid']) >0: # read void coordinates
                    cl.end = cl.end.add_atoms(params['xvoid'], 'void')
                    printlog('I found', len(params['xvoid']), 'voids in config file. Added to structure.')

                cl.serialize('0-yes')

        
        else:
            cl = vasp_run(3, 'first run')
            if 'xvoid' in params and len(params['xvoid']) >0: # read void coordinates
                cl.end = cl.end.add_atoms(params['xvoid'], 'void')
                printlog('I found', len(params['xvoid']), 'voids in config file. Added to structure.')

            cl.serialize('0-yes')
            copyfile('OSZICAR', 'OSZICAR-0')
            copyfile('CONTCAR', 'CONTCAR-0')

            with open('ENERGIES', 'w') as f:
                f.write('{:5d}  {:.5f}\n'.format(0, cl.e0))





    """Define rest required parameters"""
    st = cl.end
    z2 = st.get_surface_pos()[1]
    printlog('Position of top surface is {:3.2f}'.format(z2) )

    xcart_voids = st.get_specific_elements([300], fmt = 'x')
    # printlog
    if xcart_voids:
        voidz = [300]
        printlog('Voids were extracted from st, adding them to z group for Monte-Carlo', xcart_voids)
    else:
        voidz = None



    for i_mcstep in range(1+last_number, 1+last_number+nmcstep):
        printlog('---------------------------------', imp = 'y')
        printlog('\n\n\n\nMonte-Carlo step = ', i_mcstep, imp = 'y')


        """3. Exchange two atoms"""
        #choose two atoms 
        z_groups = [ALKALI_ION_ELEMENTS, TRANSITION_ELEMENTS]
        if voidz:
            z_groups.append(voidz)


        printlog('All Z groups are ', z_groups)
        # sys.exit()
        gr1 = random.choice(z_groups)
        z_groups.remove(gr1)
        gr2 = random.choice(z_groups)

        printlog('Chosen Z groups are ', gr1, gr2)

        # print(st.get_elements_z())
        # sys.exit()
        alk = st.get_specific_elements(gr1, z_range = [z2 - thickness, z2])
        tra = st.get_specific_elements(gr2, z_range = [z2 - thickness, z2])
        
        
        if len(alk) == 0 or len(tra) == 0:
            printlog('Attention, alk or tra are too small:', alk, tra, 'trying another')
            continue
        
        printlog('Two groups of atom numbers to swap are', alk, tra)
        # sys.exit()

        at1 = random.choice(alk)
        if at1 in tra:
            tra.remove(at1)

        at2 = random.choice(tra)

        els = st.get_elements()
        st_new_init = st.swap_atoms(at1, at2)
        printlog('I swapped', at1, els[at1],  'and', at2, els[at2], imp = 'y' )

        if voidz:
            xcart_voids = st_new_init.get_specific_elements([300], fmt = 'x')
            printlog('The following voids after changes were extracted from st:', xcart_voids)

        if debug:
            st_new_init.write_poscar('POSCAR-'+str(i_mcstep))
            st_new_init = read_poscar(st_new_init, 'POSCAR-'+str(i_mcstep))
            st_new_init = st_new_init.add_atoms(xcart_voids, 'void')
            # xcart_voids = st_new_init.get_specific_elements([300], fmt = 'x')
            # print('After writing and reading the voids are ', xcart_voids)
            st = st_new_init


        else:
            """4. Write new structure and calculate energy  """

            st_new_init.write_poscar('POSCAR') #here voids are lost
            
            cl_new = vasp_run(3, 'mcstep '+str(i_mcstep))

            if check(cl_new):
                printlog('{:5d} is unlucky configuration, trying another ... '.format(i_mcstep), imp = 'y')
                continue



            """5. Check if to accept new structure  """
            printlog('Energies before and after are {:3.3f} and {:3.3f}, dE = {:3.3f}'.format(cl.e0, cl_new.e0, cl_new.e0 - cl.e0), imp = 'y')
            with open('ENERGIES', 'a') as f:
                f.write('{:5d}  {:.5f}\n'.format(i_mcstep, cl_new.e0))
            
            if metropolis(cl.e0, cl_new.e0, temperature):
                cl = cl_new

                if voidz:
                    #insert voids
                    cl_new.end = cl_new.end.add_atoms(xcart_voids, 'void') # here voids are inserted back

                
                cl.serialize(str(i_mcstep)+'-yes')
                copyfile('CONTCAR', 'CONTCAR_last')
                copyfile('OUTCAR', 'OUTCAR_last')

                st = cl_new.end
            else:
                cl_new.serialize(str(i_mcstep)+'-no')
            copyfile('OSZICAR', 'OSZICAR-'+str(i_mcstep))
            copyfile('CONTCAR', 'CONTCAR-'+str(i_mcstep))

    copyfile('OUTCAR_last', 'OUTCAR')
    printlog('MC simulation finished!', imp = 'y')
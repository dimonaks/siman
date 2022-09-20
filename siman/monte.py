#!/bin/python
"""  
Program Monte-Carlo by Aksyonov Dmitry, Skoltech, Moscow

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


def check_poscar(filename):
    # return 0 if ok, return 1 if failed
    try:
        cl = CalculationVasp()
        cl.read_poscar(filename)
        print(cl.init.natom)
        status = 0 # good
    except:
        status = 1
    return status



def vasp_run(n, des, vasprun_command = None):
    #allows to run vasp several times, here fireworks can be used
    #n - number of attempts
    #des - description of run
    for i in range(n): # max three attempts
        printlog(des, 'attempt', i)
        if not debug2:
            out = runBash(vasprun_command)
            printlog('out is', out)
        
        cl = CalculationVasp(output = 'OUTCAR')
        out = cl.read_results(show = 'fo')
        printlog('Results are', imp = 'y')
        printlog(out, imp = 'y')

        status = check(cl)
        if status == 0:
            break
        else:
            if os.path.exists('CONTCAR'):
                if check_poscar('CONTCAR') == 0:
                    copyfile('CONTCAR', 'POSCAR')
                else:
                    printlog('CONTCAR is broken. No further attempts to run VASP', imp = 'y')
                    break

            else:
                printlog('No CONTCAR was found. No further attempts to run VASP', imp = 'y')
                break

    return cl   



def initial_run(xcart_voids, ):
    """1. Run initial calculation"""
    if debug:
        cl = CalculationVasp()
        cl.read_poscar('1.POSCAR')
        cl.end = cl.init
        if xcart_voids: # read void coordinates
            cl.end = cl.end.add_atoms(xcart_voids, 'void')
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
            xcart_voids = cl.end.get_specific_elements([300], fmt = 'x') # extract voids form the last calculation
        
        else:
            cl = vasp_run(3, 'first run', vasprun_command = vasprun_command)
            
            if xcart_voids: # read void coordinates
                cl.end = cl.end.add_atoms(xcart_voids, 'void')
                printlog('I found', len(xcart_voids), 'voids in config file. Added to structure.')

            if check(cl) == 0:
                cl.serialize('0-yes')
                copyfile('OUTCAR', 'OUTCAR-0')
                copyfile('OSZICAR', 'OSZICAR-0')
                copyfile('CONTCAR', 'CONTCAR-0')
                copyfile('OUTCAR', 'OUTCAR_last')
                copyfile('CONTCAR', 'CONTCAR_last')
                with open('ENERGIES', 'w') as f:
                    f.write('{:5d}  {:.5f}\n'.format(0, cl.e0))
            else:
                printlog('Calculation is broken, no data was saved, exiting ...', imp = 'y')
                sys.exit()



    if debug2:
        sys.exit()

    return cl, last_number

def get_zr_range(st, thickness, zr):

    red_thick = thickness/np.linalg.norm(st.rprimd[2])
    # z_range = [z2 - thickness, z2]
    zr_range = [zr - red_thick, zr]
    printlog('zr_range is ', zr_range)
    return zr_range


def exchange_atoms(st, xcart_voids, zr = None, condition = None, params = None):
    """
    Swap two atoms 

    INPUT:
        - st (Structure) - input structure
        - xcart_voids (list) - xcart of voids
        - zr (float) - reduced position of the surface
        - condition (str) - possible additional conditions:
            - 'no_surface_TM' - do not make swaps which reduce oxygen coordination of transition metals
            - 'max_avdist_increase' - maximum allowed increase of TM-O distance after swapping; 
                    (for example larger than 0.5 A allows to exclude swaps to surface) 

        - params (dict) - params from json


    COMMENTS:
        voidz - list with z voids; actually either None or [300]
        z_groups = params.get('z_groups') # z_groups for exchange, list of lists; e.g. [[300], [42]] if None than all groups are used

    """

    if xcart_voids:
        voidz = [300]
        printlog('Voids were extracted from st, adding them to z group for Monte-Carlo', xcart_voids)
    else:
        voidz = None

    if params.get('z_groups'):
        z_groups = params['z_groups']
    else:
        z_groups = [AM, TM]
        if voidz:
            z_groups.append(voidz)


    printlog('All Z groups are ', z_groups)
    # sys.exit()

    if params.get('thickness'):
        zr_range = get_zr_range(st, thickness, zr)
    else:
        zr_range = None

    for i in range(100): # try 100 attempts until the condition is satisfied, otherwise terminate
        z_groups_cp = copy.deepcopy(z_groups)
        # print('z_groups_cp', z_groups_cp)
        gr1 = random.choice(z_groups_cp)
        z_groups_cp.remove(gr1)
        gr2 = random.choice(z_groups_cp)

        printlog('Chosen Z groups are ', gr1, gr2)

        # print(st.get_elements_z())
        # sys.exit()
        # nn1 = st.get_specific_elements(gr1, z_range = z_range) # atom numbers
        # nn2 = st.get_specific_elements(gr2, z_range = z_range)

        nn1 = st.get_specific_elements(gr1, zr_range = zr_range) # atom numbers
        nn2 = st.get_specific_elements(gr2, zr_range = zr_range)
        
        
        if len(nn1) == 0 or len(nn2) == 0:
            printlog('Attention, nn1 or nn2 are too small:', nn1, nn2, 'trying another')
            # print(st.get_elements())
            print(gr1, gr2, zr_range)
            print([st.xred[i] for i in st.get_specific_elements([300]) ])
            # sys.exit()
            continue
        
        printlog('Two groups of atom numbers to swap are', nn1, nn2)
        # sys.exit()

        at1 = random.choice(nn1)
        if at1 in nn1:
            nn1.remove(at1)

        at2 = random.choice(nn2)

        els = st.get_elements()
        st_new_init = st.swap_atoms(at1, at2)
        printlog('I swapped', at1+1, els[at1],  'and', at2+1, els[at2], imp = 'y' )
    

        #condition check-up
        if condition == 'no_surface_TM':
            elsz = st_new_init.get_elements_z()
            z1 = elsz[at1]
            z2 = elsz[at2]
            if (z1 in TM and z2 in TM) or (z1 not in TM and z2 not in TM):
                break # nothing to do
            # elif z1 in TM or z2 in TM:
            if z1 in TM:
                atTM = at1
            else:
                atTM = at2
            printlog('I found that one swapping atom is transition metal', atTM, els[atTM], 'checking coordination')

            # nO1 = st.nn(atTM,          6, only = [8], from_one = 0)['el'].count('O')
            # nO2 = st_new_init.nn(atTM, 6, only = [8], from_one = 0)['el'].count('O')

            av1 = st.nn(atTM,          6, only = [8], from_one = 0, silent = 1)['av(A-O,F)']
            av2 = st_new_init.nn(atTM, 6, only = [8], from_one = 0, silent = 1)['av(A-O,F)']
            
            # printlog('The oxygen-TM average', av1, av2, imp = 'y')

            if av2 > av1+0.5:
                printlog('Surface TM detected, the TM-O average distances before and after are {:.2f} {:.2f} A. Trying another swap.'.format(av1, av2), imp = 'y')
            else:
                printlog('TM-O av. dist before and after are {:.2f} {:.2f} A. Good, accepted'.format(av1, av2), imp = 'y')
                break


            # if nO1 == nO2:
            #     printlog('The oxygen coordination of TM after swap is the same, accepting', nO1, nO2)
            #     break
            # else:
            #     printlog('Warning! The oxygen coordination of TM was reduced, trying another step:', nO1, nO2)



    else:
        printlog('exchange_atoms(): The given condition on atom swapping cant be satisfied! exitting', imp = 'y' )
        sys.exit()



    return st_new_init


def exchange_with_external(st, zr, thickness, external = None):
    """
    external (dict) - {'Ni':['Li']} - atoms from external reservior which can replace existing elements, in this example Ni can replace Li
    """
    
    printlog('Starting exchange with external')

    zr_range = get_zr_range(st, thickness, zr)

    # get  external 
    # print(external.keys())
    el_ext = random.choice(list(external.keys()))

    # get element to change
    el_int = random.choice(external[el_ext])


    nn1 = st.get_specific_elements([invert(el_int)], zr_range = zr_range)

    nn2 = st.get_specific_elements([invert(el_ext)], zr_range = zr_range) # just to know good TM-O distance

    j = random.choice(nn2)

    av2 = st.nn(j, 6, only = [8], from_one = 0, silent = 1)['av(A-O,F)'] # in slab

    if len(nn1) == 0:
        printlog('All atoms were replaced, exiting ', imp = 'y')
        sys.exit()


    for k in range(10):
        i = random.choice(nn1)
        av1 = st.nn(i,          6, only = [8], from_one = 0, silent = 1)['av(A-O,F)']
        if av1 < av2+0.4:
            printlog('The chosen position is compared to the existing in slab for {:s}: {:.2f} {:.2f} A'.format(el_ext, av1, av2), imp = 'y')
            break
        else:
            printlog('Trying ', i)
    else:
        printlog('No more good options, trying all', imp = 'y')

    st_new = st.replace_atoms([i], el_ext)
    printlog('Atom', el_int, i, 'was replaced by ', el_ext, 'from external reservoir')

    return st_new 

if __name__ == "__main__":


    debug = 0

    header.warnings = 'yY'
    # header.warnings = 'neyY'
    header.verbose_log = 1

    printlog('\n\n\nStarting Monte-Carlo script!')
    

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
    thickness = params.get('thickness') # minimum thickness of layer below the surface where exchanges are performed. If None than full cell is considered active
    temperature = params.get('temp') or 1
    xcart_voids = params.get('xvoid')



    if params.get('external'):
        if not params.get('chem_pot'):
            printlog('Error! no chem_pot parameter detected in external mode')
        cl_bulk = CalculationVasp().deserialize_json('bulk.json')
        printlog('Mode with external chemical potential. Reading of bulk structure is OK', cl_bulk.e0)




    printlog('Command to run vasp', vasprun_command)
    printlog('Total number of steps is', nmcstep)
    printlog('Thickness is ', thickness)
    printlog('Temperature is ', temperature, 'K')


    cl, last_number = initial_run(xcart_voids, )
    st = cl.end


    """Determine surface position"""
    z2 = st.get_surface_pos()[1]
    zr2 = st.get_surface_pos(reduced = True)[1]
    printlog('Position of top surface is {:3.2f} {:3.2f}'.format(z2, zr2) )

    # printlog




    """Start Monte-Carlo"""

    for i_mcstep in range(1+last_number, 1+last_number+nmcstep):
        printlog('---------------------------------', imp = 'y')
        printlog('\n\n\n\nMonte-Carlo step = ', i_mcstep, imp = 'y')




        """3. The section where changes are done """

        if params.get('external'):
            st_new_init = exchange_with_external(st, zr2, thickness, external = params.get('external'))
        else:
            st_new_init = exchange_atoms(st, xcart_voids, zr = zr2, condition = 'no_surface_TM', params = params)








        if xcart_voids:
            xcart_voids = st_new_init.get_specific_elements([300], fmt = 'x')
            printlog('The following voids after changes were extracted from st:', xcart_voids)

        if debug:
            st_new_init.write_poscar('POSCAR-'+str(i_mcstep))
            st_new_init = read_poscar(st_new_init, 'POSCAR-'+str(i_mcstep))
            # print(xcart_voids)
            # sys.exit()
            if xcart_voids:
                st_new_init = st_new_init.add_atoms(xcart_voids, 'void')
            # xcart_voids = st_new_init.get_specific_elements([300], fmt = 'x')
            # print('After writing and reading the voids are ', xcart_voids)
            st = st_new_init
            cl = CalculationVasp()
            cl_new = CalculationVasp()
            cl.end = st
            cl_new.end = st_new_init
            cl.e0 = random.random()
            cl_new.e0 = random.random()
        else:
            """4. Write new structure and calculate energy  """

            st_new_init.write_poscar('POSCAR') #here voids are lost
            
            cl_new = vasp_run(3, 'mcstep '+str(i_mcstep), vasprun_command = vasprun_command)

            if check(cl_new):
                printlog('{:5d} is unlucky configuration, trying another ... '.format(i_mcstep), imp = 'y')
                continue



        """5. Check if to accept new structure  """
        if params.get('external'):
            # print('')
            gamma, E1         = suf_en(cl,     cl_bulk, chem_pot = params.get('chem_pot'), return_diff_energy =1)
            gamma_new, E2     = suf_en(cl_new, cl_bulk, chem_pot = params.get('chem_pot'),  return_diff_energy =1)

        else:
            E1 = cl.e0
            E2 = cl_new.e0



        printlog('Energies before and after are {:3.3f} and {:3.3f}, dE = {:3.3f}'.format(E1, E2, E2 - E1), imp = 'y')
        with open('ENERGIES', 'a') as f:
            f.write('{:5d}  {:.5f}\n'.format(i_mcstep, cl_new.e0))


        
        if metropolis(E1, E2, temperature):
            cl = cl_new

            if xcart_voids:
                #insert voids
                cl_new.end = cl_new.end.add_atoms(xcart_voids, 'void') # here voids are inserted back

            if not debug:
                cl.serialize(str(i_mcstep)+'-yes')
                copyfile('CONTCAR', 'CONTCAR_last')
                copyfile('OUTCAR', 'OUTCAR_last')
                if os.path.exists('LOCPOT'):
                    copyfile('LOCPOT', 'LOCPOT_last')

            st = cl_new.end

            printlog('The step was accepted', imp = 'y')

        else:
            printlog('The step was rejected', imp = 'y')

            if not debug:
                cl_new.serialize(str(i_mcstep)+'-no')

        if not debug:
            copyfile('OSZICAR', 'OSZICAR-'+str(i_mcstep))
            copyfile('CONTCAR', 'CONTCAR-'+str(i_mcstep))
            copyfile('OUTCAR', 'OUTCAR-'+str(i_mcstep))
            if os.path.exists('LOCPOT'):

                copyfile('LOCPOT', 'LOCPOT-'+str(i_mcstep))

    if not debug:
        copyfile('OUTCAR_last', 'OUTCAR')
        copyfile('CONTCAR_last', 'CONTCAR')
        if os.path.exists('LOCPOT_last'):
        
            copyfile('LOCPOT_last', 'LOCPOT')
    
    printlog('MC simulation finished!', imp = 'y')
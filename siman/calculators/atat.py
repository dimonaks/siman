# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
"""
This module contains interface tools to work with ATAT software

Contributors: 
Aksyonov Dmitry
"""
import shutil, re, sys
from siman.small_functions import calc_ngkpt, list2string
from siman.header import printlog
from siman import header

def setup_atat_step1(varset, setlist, input_st, params):
    """
    Make necessary adjusments for ATAT 

    INPUT:
        - varset (dict) - dict of sets from internal database
        - setlist (list) - list of sets for current calculation
        - input_st (Structure) - current structure
        - params (dict) - dict for add()

    """

    ks = varset[setlist[0]].vasp_params['KSPACING']
    input_st = input_st.reorder_element_groups(order = 'alphabet') # required for correct work of ATAT
    # input_st.printme()
    # sys.exit()

    N = calc_ngkpt(input_st.get_recip(), ks)
    # print(N)
    KPPRA = N[0]*N[1]*N[2]*input_st.natom # here probably symmetry should taken into account? not really good working

    if ks == 0.3:
        KPPRA = 1200
    else:
        KPPRA = 1200
        printlog('Warning! KPPRA = 1200 for all kspacings! please contact developer to improve or modify this code in calc_manage.py')

    printlog('Atat mode, KPPRA is', KPPRA, imp = 'Y')
    # sys.exit()
    exclude_new = []
    exclude = params['atat'].get('exclude_atoms_n')

    subatom = params['atat'].get('subatom') or None
    # print(exclude)
    # sys.exit()
    if exclude:
        for i in exclude:
            exclude_new.append(input_st.old_numbers.index(i)) # required since the structure was reordered
        # exclude = exclude_new
        params['atat']['exclude_atoms_n'] = exclude_new

    params['update_set_dic']={'add_nbands':None, 'mul_nbands':None, 'USEPOT':'PAWPBE', 'KPPRA':KPPRA, 
    'MAGATOM':list2string(input_st.magmom), 
    'MAGMOM':None,
    'SUBATOM':subatom,
    'DOSTATIC':''}


    rm_files = 'rm maps_is_running pollmach_is_running; \n'

    maps_keys = params['atat'].get('maps_keys') or '-d' # default
    # print(maps_keys)
    # sys.exit()

    if len(params['atat']['active_atoms']) == 1:
        maps = 'maps '
    else:
        maps = 'mmaps '

    header.atat_run_command = rm_files+maps+maps_keys+'&\npollmach runstruct_vasp mpirun\n'

    return



def write_lat_in(st, params, dirpath):
    file = dirpath +  'lat.in'
    to_ang = 1
    rprimd = st.rprimd
    with open(file, 'w') as f:
        for i in 0, 1, 2:
            f.write('{:10.6f} {:10.6f} {:10.6f}\n'.format(rprimd[i][0]*to_ang,rprimd[i][1]*to_ang,rprimd[i][2]*to_ang) )
            # f.write("\n")
        f.write(' 1 0 0\n 0 1 0\n 0 0 1\n')

        active_atoms = params['atat']['active_atoms'] #dict
        
        exclude = params['atat'].get('exclude_atoms_n') or []
        # print(exclude)
        # sys.exit()
        subs = []

        active_numbers = []
        subs_dict    = {}
        active_elements = []
        for el in set(st.get_elements()):
            for elan in active_atoms:
                
                if el not in elan:
                    continue


                # print(elan)
                active_elements.append(el)
                elann = re.split('(\d+)',elan)
                printlog('Exctracting symmetry position of Na ', elann)
                if len(elann) > 2:
                    ela = elann[0]
                    isym   = int(elann[1])
                else:
                    ela = elann[0]
                    isym = None
                
                if isym:
                    natoms = st.determine_symmetry_positions(el)
                    a_numbers = natoms[isym-1]
                else:
                    a_numbers = st.get_numbers(el)
        
                for i in a_numbers:
                    subs_dict[i] = active_atoms[elan]
                
                active_numbers.extend(a_numbers)
        # print(active_numbers)
        # sys.exit()


        # st.printme()
        # sys.exit()


        for i, el in enumerate(st.get_elements()):
            if i not in exclude and i in active_numbers:
                # print(el, i, st.xred[i])
                subs.append(subs_dict[i])
            else:
                subs.append(None)

        printlog('Substitutions are', subs,)
        # print('subs=',subs)

        # print(st.magmom)
        # sys.exit()
        if len(st.magmom) == 0 or None in st.magmom:
            write_magmom = False
            magmom = st.natom*['']
        else:
            magmom = st.magmom
            write_magmom = True

        # print('magmom', magmom, st.magmom,st.natom)
        # sys.exit()

        for x, el, sub, m in zip(st.xred, st.get_elements(), subs, magmom):
            # if el == 'O':
            #     m = 0
            # print(m)
            if m and abs(m) < 0.1:
                m = 0
            if write_magmom:
                f.write('{:10.6f} {:10.6f} {:10.6f} {:s}'.format(*x, el))
            else:
                f.write('{:10.6f} {:10.6f} {:10.6f} {:s}'.format(*x, el))
            
            if sub:
                f.write(','+sub)
            f.write("\n")

    #highlight active atoms
    st_h = st.replace_atoms(active_numbers, 'Pu')
    st_h.write_poscar()
    printlog('ATAT files sucesfully created. Check active atoms using file above.', imp = 'y')



    return file



def write_atat_crange(params, dirpath):
    """
    Write file with crange

        - params (dict) - dict from add() with 'atat' keyword

    """
    filename = dirpath+'/crange.in' 
    with open(filename, 'w') as f:
        for constr in params['atat']['constraints']:
            f.write(constr+'\n')

    return filename




def setup_atat_step2(cl, params, list_to_copy):



    lat_in = write_lat_in(cl.init, params, cl.dir)

    crange_in = write_atat_crange(params, cl.dir)

    list_to_copy.append(lat_in)
    list_to_copy.append(crange_in)

    wrap = cl.dir+'/vasp.wrap'
    shutil.copyfile(cl.dir+'/INCAR', wrap)

    with open(wrap, "r") as fwr:
        wrap_cont = fwr.readlines()

    with open(wrap, "w") as fwr:
        fwr.write("[INCAR]\n"+"".join(wrap_cont[1:])) # remove first line with SYSTEM tag - not working in ATAT for some reason

    list_to_copy.append(wrap)

    return
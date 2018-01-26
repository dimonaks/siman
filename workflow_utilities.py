from __future__ import division, unicode_literals, absolute_import 
import sys, copy, re, os

#external
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import header
from header import print_and_log as printlog
from header import calc, db
from picture_functions import fit_and_plot
from small_functions import merge_dics
from calc_manage import add_loop, name_mod_supercell, res_loop, inherit_icalc, push_figure_to_archive
from neb import add_neb
from classes import Calculation
from analysis import calc_redox,  matrix_diff
from geo import create_deintercalated_structure
from geo import remove_one_atom
from geo import create_replaced_structure, create_antisite_defect3, determine_symmetry_positions, image_distance
from inout import write_occmatrix







def prepare(it_new, opt_vol, it_folder, ise, cl, st_type, option):

    if not it_folder:
        if option == None:
            option = ''
        it_folder = header.struct_des[cl.id[0]].sfolder+'/'+option
    if st_type == 'end':
        st = cl.end
    elif st_type == 'init':
        st = cl.init


    if not ise:
        if opt_vol:
            ise = '4uis'
        else:
            ise = '1u'

    if opt_vol:
        id_new = (it_new+'.su', ise, 100)
    else:
        id_new = (it_new, ise, 1)

    return id_new, st, it_folder



def make_defect(cl, el, st_type = 'end', option = 'vac', pos = None, ise = None, opt_vol = 0, 
    suf = '', it_folder = None,
    el_rep = '', pos_rep = 1, pos_rep2 = None, polaron_pos = None, occ_matrix = None,
    up = 0, fit = 0,  outcar = None, only_read = 0,
    compat1 = False, add_loop_arg = {}):
    """
    Function allow to create point defects and run them
	previous name: make_vacancy()


    cl - starting Calculation 
    st_type - starting structure of cl: 'init' or 'end' 
    el - element to be removed or replaced
    
    option -
        'vac'  - make vacancy
        'rep'  - replace one atom with Ti, was used for V-Ti project
        'pair' - make vacancy -Ti complex for V-Ti project 

    pos - unique position of el if non-eqivalent atoms exist - for vac
    ise - new set
    opt_vol (bool) - optimize volume

    suf (str) - mannually added suffix
    it_folder - mannually provided it_folder

    up (bool) - [ 0, 1 ] update current calculation
    fit = 0,  outcar = None, only_read = 0 - flow control as usual

    polaron_pos - choose polaron position
    occ_matrix - list of lists see format in classes


    compat1 - compatability with previous calculations, which were used for Na2FePO4F project


    TODO: rename to ?_point_defects()
    """





    if pos == None:
        pos = ''
    
    if polaron_pos == None:
        pol_suf = '' 
    else:
        pol_suf = '.p'+str(polaron_pos) # polaron suffix



    ssuf = el+str(pos)+el_rep+pol_suf+suf
    if 'su' in cl.id[0] and not 'su.' in cl.id[0]:
        it_new = cl.id[0].replace('su', option) + ssuf
    else:
        it_new = cl.id[0] + option + ssuf

    if compat1: # no element in name
        it_new = cl.id[0].replace('su', 'vac')+str(pos)

    id_new, st, it_folder = prepare(it_new, opt_vol, it_folder, ise, cl, st_type, option)

    if not only_read and (up or id_new not in calc):
        # it_new

        if 'vac' in option:
            st_del1, i_del = remove_one_atom(st, el, pos)
            st_vis  = st.replace_atoms([i_del], 'U')
            st_vis.name=it_new+'_visual'
            st_vis.write_xyz()

            #possible polaron positions
            tr = st.get_transition_elements(fmt = 'z') 
            i_tr = st.get_transition_elements(fmt = 'n') 
            # dist = []
            max_d = 0
            for i in i_tr:
                d1, d2 = image_distance(st.xcart[i], st.xcart[i_del], st.rprimd)
                # print(i+1, d1, d2)
                if d1 < d2:
                    if d1 > max_d:
                        max_d = d1
                        i_max_d = i

            print('The longest distance to transition metal in current supercell is ', max_d, 'A for atom', i_max_d+1, st.get_elements()[i_max_d])


            numb = st.nn(i_del, from_one = 0, n = len(tr)+5, only = list(set(tr)))['numbers'][1:]
            
            printlog('Choose polaron position starting from 1 using *polaron_pos*', imp = 'y')
            if polaron_pos:
                i_pol = numb[polaron_pos-1]
                printlog('atom', i_pol+1, st.get_elements()[i_pol], 'is chosen', imp ='y')
                # print(numb)
                # sys.exit()
                #take_occupation matrices from cl
                print('substitution occupation matrix of atom', i_pol+1)
                occ_matrices = copy.deepcopy(cl.occ_matrices)
                occ_matrices[i_pol] = occ_matrix
                # print(pd.DataFrame(cl.occ_matrices[i_pol]))
                occfile = write_occmatrix(occ_matrices, cl.dir+'/occ/')
                # print(occfile)

                # sys.exit()



        elif 'rep' in option:
            st_del1 = st.replace_atoms([pos_rep], el)
            print('Atom', str(pos_rep),st.get_elements()[pos_rep],' replaced with', el, )
            print(st_del1.get_elements()[pos_rep])
            st_del1.name+=it_new
        
        elif 'pair' in option:
            st_del1 = st.replace_atoms([1], el_rep)
            if 'pair2' in option:
                st_del1 = st_del1.replace_atoms([pos_rep2], el_rep)

            st_del1 = remove_one_atom(st_del1, el, pos)
            print('Atom 1 replaced with', el,'and atom removed' )
            st_del1.name+=it_new


        st_del1.write_xyz()
        



        if opt_vol:
            it = add_loop(it_new, ise, 1, calc_method = 'uniform_scale',
             scale_region = (-4, 4), inherit_option = 'inherit_xred', input_st = st_del1, it_folder = it_folder, 
             params = {'occmatrix':occfile}, **add_loop_arg)            
        else:
            it = add_loop(it_new, ise, 1, input_st = st_del1, it_folder = it_folder, 
            params = {'occmatrix':occfile}, **add_loop_arg)             





    else:
        if opt_vol and fit:
            res_loop(it_new+'.su', ise, list(range(1,8))+[100], analys_type = 'fit_a', show = 'fitfo', up = '2', choose_outcar = outcar)
        else:
            res_loop(*id_new, up = '2', choose_outcar = outcar, show = 'fo')
            # calc[it_new+'.su', ise, 100].end.jmol()

        cl_v = calc[id_new]
        if '4' not in cl.state:
            cl.res()
        calc_redox(cl_v, cl)
        # print(cl_v.end.vol, cl.end.vol)
        dE = None


        if option == 'vac':
            print('Evac = {:3.2f} eV'.format(cl_v.e0 - cl.e0/cl.end.natom*cl_v.end.natom))    
        
        elif option == 'rep':
            diffE, diffV = matrix_diff(cl_v, cl)

            print('Esol = {:3.2f} eV'.format(diffE))    

        elif 'pair' in option:
            ''
            cl_bulk = cl
            cl_pair = cl_v
            it = cl.id[0]
            if 'V54' in it:
                it = it.replace('.su', '.')
            id_vac = (it + 'vacV', cl.id[1], 1)
            id_sol = (it + 'repTi', cl.id[1], 1)
            cl_vac  = calc[id_vac]
            cl_sol  = calc[id_sol]

            # print(id_vac, id_sol)

            # if  option == 'pair':
            dE = e_bind(cl_bulk, cl_vac, cl_sol, cl_pair)
            print('Ecomplex = {:3.2f} eV'.format(dE))    

            # elif '2' in option:


        return {'dE':dE, 'N':cl_v.end.natom, 'Name':cl.id}



def process_modified(cl, mod_dic = None, scale_region = (-4,4), opt_vol = 1, fit = 0,  st_type = 'end', name = None, el_new = None, run = 0, ise = None, it_folder = None, mode = None, add_loop_arg = None):
    """
    inherited from create_charges - functionality is extended
    The utility allows to (contrlolled by mode parameter):
    1) create charged cells by removing specific atoms provided in del_dic
    2) replace specific atoms 

    add_loop
    res_loop


    mode - 
        delete
        remove
        None

    mod_dic - dic of configurations with atom numbers starting from 1

    """
    # if not del_dic:
    if add_loop_arg == None:
        add_loop_arg = {}

    if mod_dic == None:
        mod_dic = {1:1}

    for key in mod_dic:
        mod_pos = mod_dic[key]
        if mode:
            mod_pos = [p-1 for p in mod_pos]

        if mode:
            suf = '.'+mode[0]+str(key)
        else:
            suf = ''
        
        it_new = cl.id[0] + suf
        
        id_new, stA, it_folder = prepare(it_new, opt_vol, it_folder, ise, cl, st_type, mode)

        if run: 
            if mode == 'delete':
                st = stA.remove_atoms(mod_pos)
            elif mode == 'replace':
                st = stA.replace_atoms(atoms_to_replace = mod_pos, el_new = el_new)
            else:
                st = stA

            st.name+=suf
            st.write_xyz()
            # sys.exit()
            if opt_vol:
                add_loop(it_new, id_new[1], 1, calc_method = 'uniform_scale',
                 scale_region = scale_region, inherit_option = 'inherit_xred', input_st = st, it_folder = it_folder, **add_loop_arg)            
            else:
                add_loop(it_new, id_new[1], 1, input_st = st, it_folder = it_folder, **add_loop_arg)             




        else:
            if 'check_job' in add_loop_arg:
                cj = add_loop_arg['check_job']
            else:
                cj = None
            if opt_vol and fit:
                # res_loop(id_new[0], id_new[1], list(range(1,8))+[100], analys_type = 'fit_a', show = 'fitfo', up = '2', choose_outcar = None)
                res_loop(id_new[0], id_new[1], list(range(1,8))+[100], analys_type = 'fit_a', show = 'fitfo', up = '2', choose_outcar = None, check_job = cj)
            else:
                res_loop(*id_new, up = '1', choose_outcar = None, show = 'fo')
        
    return

from __future__ import division, unicode_literals, absolute_import 
import sys, copy, re, os

#external
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import header
from header import print_and_log as printlog
from header import calc
from picture_functions import fit_and_plot
from small_functions import merge_dics
from calc_manage import add_loop, name_mod_supercell, res_loop, inherit_icalc, push_figure_to_archive, add_to_database
from neb import add_neb
from classes import Calculation
from analysis import calc_redox,  matrix_diff
from geo import create_deintercalated_structure
from geo import remove_one_atom
from geo import create_replaced_structure, create_antisite_defect3, determine_symmetry_positions

def make_defect(cl, el, st_type = 'end', option = 'vac', pos = None, ise = None, opt_vol = 0, 
    suf = '', it_folder = None,
    el_rep = '', pos_rep = 1, pos_rep2 = None, 
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

    pos - unique position of el if non-eqivalent atoms exist
    ise - new set
    opt_vol (bool) - optimize volume

    suf (str) - mannually added suffix
    it_folder - mannually provided it_folder

    up (bool) - [ 0, 1 ] update current calculation
    fit = 0,  outcar = None, only_read = 0 - flow control as usual


    compat1 - compatability with previous calculations, which were used for Na2FePO4F project


    TODO: rename to ?_point_defects()
    """





    if pos == None:
        pos = ''
    
    if 'su' in cl.id[0] and not 'su.' in cl.id[0]:
        it_new = cl.id[0].replace('su', option)+el+str(pos)+el_rep+suf
    else:
        it_new = cl.id[0]+option+el+str(pos)+el_rep+suf


    if compat1: # no element in name
        it_new = cl.id[0].replace('su', 'vac')+str(pos)

    if not it_folder:
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

    if not only_read and (up or id_new not in calc):
        

        if 'vac' in option:
            st_del1 = remove_one_atom(st, el, pos)
        
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
            add_loop(it_new, ise, 1, calc_method = 'uniform_scale',
             scale_region = (-4, 4), inherit_option = 'inherit_xred', input_st = st_del1, it_folder = it_folder, **add_loop_arg)            
        else:
            add_loop(it_new, ise, 1, input_st = st_del1, it_folder = it_folder, **add_loop_arg)             



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

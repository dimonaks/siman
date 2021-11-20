from __future__ import division, unicode_literals, absolute_import 
import sys, copy, re, os

#external
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from siman import header
from siman.header import printlog, calc, db, struct_des
from siman.picture_functions import fit_and_plot
from siman.small_functions import merge_dics
from siman.calc_manage import add_loop, name_mod_supercell, res_loop, inherit_icalc, push_figure_to_archive
from siman.neb import add_neb
from siman.classes import Calculation
from siman.analysis import calc_redox,  matrix_diff
from siman.geo import create_replaced_structure, create_antisite_defect3, determine_symmetry_positions, image_distance, remove_one_atom, replic, local_surrounding, create_deintercalated_structure
from siman.inout import write_occmatrix, write_xyz
from siman.functions import invert







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
    up = 0, fit = 0,  outcar = None, only_read = 0, Eref = 0,
    compat1 = False, add_loop_arg = {}):
    """
    Function allow to create point defects and run them
	previous name: make_vacancy()


    cl - starting Calculation 
    st_type - starting structure of cl: 'init' or 'end' 
    el - element to be removed or replaced
    
    option -
        'vac'  - make vacancy
        'rep'  - replace one atom with 'el_rep', 
        'pair' - make vacancy -Ti complex for V-Ti project 

    pos - unique position of el if non-eqivalent atoms exist - for vac
    pos_rep - number of position to replace from 0

    ise - new set
    opt_vol (bool) - optimize volume

    suf (str) - mannually added suffix
    it_folder - mannually provided it_folder

    up (bool) - [ 0, 1 ] update current calculation
    fit = 0,  outcar = None, only_read = 0 - flow control as usual

    polaron_pos - choose polaron position
    occ_matrix - list of lists see format in classes


    compat1 - compatability with previous calculations, which were used for Na2FePO4F project


    Eref - reference energy for solution energy

    TODO: rename to ?_point_defects()
    """

    from siman.project_funcs import e_bind




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
    occfile = None
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
            i_max_d = None
            for i in i_tr:
                d1, d2 = image_distance(st.xcart[i], st.xcart[i_del], st.rprimd)
                # print(i+1, d1, d2)
                if d1 < d2:
                    if d1 > max_d:
                        max_d = d1
                        i_max_d = i
            if i_max_d is not None:
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
            else:
                occfile = None


        elif 'rep' in option:
            st_del1 = st.replace_atoms([pos_rep], el_rep)
            print('Atom', str(pos_rep),st.get_elements()[pos_rep],' replaced with', el_rep, )
            print(st_del1.get_elements()[pos_rep])
            st_del1.name=it_new
        
        elif 'pair' in option:
            st_del1 = st.replace_atoms([1], el_rep)
            if 'pair2' in option:
                st_del1 = st_del1.replace_atoms([pos_rep2], el_rep)

            st_del1 = remove_one_atom(st_del1, el, pos)
            print('Atom 1 replaced with', el,'and atom removed' )
            st_del1.name=it_new


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
        if not hasattr(cl_v, 'e0'):
            printlog('Warning', cl_v.id, 'is bad')
            return
        calc_redox(cl_v, cl)
        # print(cl_v.end.vol, cl.end.vol)
        dE = None


        if option == 'vac':
            cl_v.res()
            cl.res()
            print('Evac = {:3.2f} eV'.format(cl_v.e0 - cl.e0/cl.end.natom*cl_v.end.natom))    
        
        elif option == 'rep':
            diffE, diffV = matrix_diff(cl_v, cl)

            print('Esol = {:3.2f} eV'.format(diffE - Eref))    

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





def create_segregation_cases(it, ise, verlist, dist_gb, gbpos = None, ise_new = None, option = None, 
    precip_folder = None, use_init = False, precision = None):
    """
    Written for Ti-Fe project.
    Allows to create segregation by substituting atoms;
    dist_gb - distance from gb inside which the atoms are included
    

    option = 'precip'- adding additional impurities to the already existing at gb. Please use 'precip_folder'
    use_init - allows to use initial structure.


    !Warning PBC are not used in determination of seg positions 

    """


    # hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
    # try:
    #     if hstring != header.history[-1]: header.history.append( hstring  )
    # except:
    #     header.history.append( hstring  )
    def write_local(cl, it_new, it_new_path, el_sub, main_path ): 
        cl.version = v
        # it_new_path   
        
        cl.name = it_new
        cl.des = 'Obtained from end state of '+str((it, ise, v))+' by substitution of one atom near gb with '+el_sub+' impurity '
        path_new_geo = it_new_path+"/"+it_new+"/"+it_new+'.imp.'+el_sub+'.'+str(cl.version)+'.'+'geo'
        cl.init.name = it_new+".init."+str(cl.version)
        xyzpath = it_new_path+"/"+it_new
        cl.path["input_geo"] = path_new_geo
        print(path_new_geo)
        cl.write_geometry("init", cl.des, override = 1)
        write_xyz(cl.init, xyzpath)
        return it_new_path

    if 0:
        res_loop(it, ise, verlist, up = 0)

    cl = header.calc[(it, ise, verlist[0])]

    znucl_sub = 3 #atom to be added


    """1. Create list of atoms near gb to substitue"""
    cl.gbpos = gbpos
    # print cl.gbpos
    seg_pos_list = [] # numbers of segregation positions
    
    if use_init:
        st = cl.init
    else:
        st = cl.end


    # print(st.xcart)
    if len(st.xcart) == 0:
        print('Warning!, xcart is empty', cl.id, )
    for i, x in enumerate(st.xcart):
        z_cur = st.znucl[st.typat[i]-1]
        print ('z_cur', z_cur)
        if z_cur == znucl_sub:
            printlog('Skipping znucl_sub atom\n')
            continue
        if abs(x[0] - gbpos)< dist_gb:
            print('adding possible seg position')
            seg_pos_list.append(i)


    """2. Substitue"""
    el_sub    = invert(znucl_sub)
    base_name = it
    main_path = header.struct_des[it].sfolder
    based_on = it+'.'+ise
    des_list = []
    add_list = []

    cl_list = []
    sumr_list = []

    i = 0
    # print(seg_pos_list)
    for j, replace_atom in enumerate(seg_pos_list): #

        v = verlist[0] # the first version from list is used
        cl = calc[(it, ise, v)]
        cl.gbpos = gbpos

        new = copy.deepcopy(cl)
        if use_init:
            new.end = new.init
        else:
            new.init = new.end #replace init structure by the end structure
        
        if 1: #atom substitution !make function; see TODO
            
            if znucl_sub not in new.init.znucl:
                new.init.znucl.append(znucl_sub)
                new.init.ntypat+=1
                new.init.typat[replace_atom] = new.init.ntypat
            else:
                ind = new.init.znucl.index(znucl_sub)
                new.init.typat[replace_atom] = ind + 1
            new.init.nznucl = []
            for typ in range(1, new.init.ntypat+1):
                new.init.nznucl.append(new.init.typat.count(typ) )
            printlog("Impurity with Z="+str(znucl_sub)+" has been substituted in "+new.name+"\n\n")
            
            it_new        = base_name+el_sub+'is'+str(i+1) #interface substitution
            
            if option == 'precip':
                it_new_path = precip_folder

            else:
                it_new_path = main_path+'/'+base_name+'_segreg'
            

        #Check if configuration is unique
        add = 1
        # for cl in cl_list:
        st = new.end
        st_replic = replic(st, (2,2,2))
        st_replic = replic(st_replic, (2,2,2), -1) #replic in negative direction also
        sumr = local_surrounding(st.xcart[replace_atom], st_replic, n_neighbours =  6) # sum of distances to surrounding atoms
        print ("sumr", sumr)
        for ad_sumr in sumr_list:
            if abs(ad_sumr-sumr) < precision :
                add = 0
                printlog("The void is non-equivalent; skipping\n")


        if add:
            i+=1
            sumr_list.append(sumr)
            # cl_list.append(new)
            write_local(new, it_new, it_new_path, el_sub, main_path) 
            for v in verlist[1:]: #write files; versions are scaled
                cl = calc[(it, ise, v)]
                rprimd_scaled = cl.end.rprimd 
                new_scaled    = copy.deepcopy(new)
                new_scaled.init.rprimd = copy.deepcopy(rprimd_scaled)
                new_scaled.init.xred2xcart()
                write_local(new_scaled, it_new, it_new_path, el_sub, main_path) 
            #create names
            des_list.append(  "struct_des['{0:s}'] = des('{1:s}', 'segregation configurations; made from {2:s}'   )".format(it_new, it_new_path, based_on)      ) 
            add_list.append(  "add_loop('"+it_new+"','"+ise_new+"',"+"range(1,6)"+", up = 'up1', it_folder = '"+it_new_path+"')"  )




    for d in des_list:
        print (d)


    for d in add_list:
        print (d)



    return








def optimize_wrapper(cl, ise, add = 0, show_fit = 1, params = None):
    #wrapper for optimization function


    up_res = 'up1'
    readfiles = 1
    check_job = 1
    id_res = (cl.id[0]+'.su', ise, 100)

    if add:
        add_loop(*cl.id, ise_new = ise,  up = 'up2', calc_method = 'uniform_scale',  
            scale_region = (-4, 4), input_st = cl.end, show = '', run = 0,  inherit_option = 'inherit_xred', 
            it_folder = cl.sfolder+'/scaled/', params = params) 
    else:

        if show_fit:
            res_loop(*id_res[0:2],list(range(0+1,0+8))+[100], up = up_res, readfiles= readfiles, analys_type = 'fit_a', show = 'fitfo')
        else:
            res_loop(*id_res, up = up_res, show = 'fo', check_job = check_job)




def run_wrapper(sts, ise = None, add = 0, cl = None, suf = 'w',  it_folder = None, cls = None, ngkpt = None, acc = None, ise1= None, acc2 = None, ise2 = None, params = None):
    """
    Add Several  structures

    params - pass to add_loop


    if add == 0:
        read results

    RETURN
    cl with lowest energy

    """
    
    if params is None:
        params = {}

    folder = suf.replace('.', '')
    # print(folder)
    # sys.exit()
    if ise1 is None:
        ise1 = ise

    if cls is None:
        cls = [cl]*len(sts)
    energies = []
    for i, st, cl_i in zip(range(len(sts)),sts,cls) : 

        itn = cl_i.id[0]+suf+ str(i)
        # del header.struct_des[itn]
        # continue
        if add:
            # add_loop(itn, ise, 1, show = 'fo', up = 'up2', input_st = st,  ngkpt = ngkpt, it_folder = cl.sfolder+'/'+folder+'/', **params ) #
            add_loop(itn, ise, 1, show = 'fo', up = 'up2', input_st = st,  ngkpt = ngkpt, it_folder = cl.sfolder, **params ) #
        
        else:
            ''
            
            if acc:
                if acc2:
                    db[itn+'.ifc', ise1, 1].run(ise2, show = 'fo', iopt = 'full_chg', add  = 0, up = 'up1', ngkpt = ngkpt)
                    cln = db[itn+'.ifc.ifc', ise2, 1]

                else:
                    # db[itn, ise, 1].run(ise1, show = 'fo', iopt = 'full_chg', add  = 0, up = 'up1', ngkpt = ngkpt)
                    cln = db[itn+'.ifc', ise1, 1]
                    cln.res(choose_outcar=0, show = 'fo')
                    suf_acc = '.ifc'
            else:
                suf_acc = ''

                res_loop(itn, ise, 1, up = 'up1')
                cln = db[itn, ise, 1]
            
            if hasattr(cln, 'e0'):
                energies.append(cln.e0)

    for i, e in enumerate(energies):
        print(i, e)

    if not add and len(energies)>0:
        i_min = energies.index(min(energies))
        print('Minimum energy is for ', i_min, energies[i_min])
        itn = cl.id[0]+suf+ str(i_min)+suf_acc
        db[itn, ise, 1].res()

        return db[itn, ise, 1]
    else:
        return None
from __future__ import division, unicode_literals, absolute_import, print_function

import copy, sys, os
import numpy as np
from operator import itemgetter
import shutil
from tabulate import tabulate

from siman import header
from siman.header import print_and_log, printlog, runBash
from siman.calc_manage import add_loop, res_loop, add_des, inherit_ngkpt
from siman.functions import  return_atoms_to_cell, push_to_server, invert
from siman.inout import write_xyz
from siman.small_functions import is_list_like, makedir, list2string
from siman.classes import CalculationVasp, cd
from siman.geo import xcart2xred, xred2xcart, local_surrounding, replic, determine_symmetry_positions
from siman.impurity import find_pores, determine_voids, determine_unique_voids










def determine_unique_final(st_pores, sums, avds, x_m):

    final_table = []
    insert_positions = []


    """Please determine unique positions with similar distances taking into acc PBC!!!
    below is incorrect
    """
    # crude_prec = 1
    # sums_crude = np.unique(sums.round(crude_prec))
    # print_and_log('The unique voids based on the sums:', 
    #     '\nwith 0.01 A prec:',np.unique(sums.round(2)),
    #     '\nwith 0.1  A prec:',sums_crude,
    #     imp ='y')
    # print_and_log('Based on crude criteria only', len(sums_crude),'types of void are relevant', imp = 'y') 

    # xcart_unique = []
    # avds_unique  = []
    # sums_unique  = []

    # for i, s in enumerate(sums_crude):
    #     index_of_first =  np.where(sums.round(crude_prec)==s)[0][0]
    #     xcart_unique.append(st_pores.xcart[index_of_first])
    #     avds_unique.append(avds[index_of_first] )
    #     sums_unique.append(sums[index_of_first] )

    # st_pores_unique = copy.deepcopy(st_pores)

    # st_pores_unique.xcart = xcart_unique
    # st_pores_unique.xcart2xred()

    sur = local_surrounding(x_m, st_pores, n_neighbours = len(st_pores.xcart), control = 'atoms', periodic  = True)

    # print('neb.determine_unique_final(): sur',  sur)

    print_and_log(
    'I can suggest you '+str (len(sur[0]) )+' end positions.', imp = 'y' )
    

    for i, (x, d, ind) in enumerate( zip(sur[0], sur[3], sur[2])):
        # if i == 0:
        #     continue
        final_table.append([i, np.array(x).round(2), round(d, 2), avds[ind], sums[ind] ]  )

    print_and_log( tabulate(final_table, headers = ['void #', 'Cart.', 'Dist', 'Dev.', 'Sum'], tablefmt='psql'), imp = 'Y' )

    return sur






def add_neb(starting_calc = None, st = None, st_end = None,
    it_new = None, ise_new = None, i_atom_to_move = None, 
    up = 'up2',
    search_type = 'vacancy_creation',
    init_neb_geo_fld = None,
    images  = None, r_impurity = None, 
    calc_method = ['neb'], 
    inherit_option  = None, mag_config = None, i_void_start = None, i_void_final = None, 
    atom_to_insert = None,
    atom_to_move   = None,
    rep_moving_atom = None,
    end_pos_types_z = None,
    replicate = None,
    it_new_folder = None,
    it_folder = None,
    inherit_magmom = False,
    x_start = None, xr_start = None,
    x_final = None, xr_final = None,
    upload_vts = False,
    center_on_moving = True,
    run = False, add_loop_dic = None, old_behaviour = None, params = None
     ):


    """
    Prepare needed files for NEB
    Provides several regimes controlled by *search_type* flag:
        - existing_voids - search for voids around atom and use them as a final position 
        - vacancy_creation - search for neighbors of the same type and make a vacancy as a start position
        - interstitial_insertion - search for two neighboring voids; use them as start and final positions
                                    by inserting atom *atom_to_insert*
        - None - just use st and st2 as initial and final
        - external - folder with subfolders 00, 01, ... in VASP format is provided using *init_neb_geo_fld*

    INPUT:
        - starting_calc (Calculation) - Calculation object with structure
        - st (Structure) - structure, can be used instead of Calculation
            - it_new (str) - name for calculation
        - st_end (Structure) - final structure

        - i_atom_to_move (int) - number of atom for moving starting from 0;
        - *mag_config* (int ) - choose magnetic configuration - allows to obtain different localizations of electron
        - *replicate* (tuple 3*int) - replicate cell along rprimd
        - i_void_start,  i_void_final (int) - position numbers of voids (or atoms) from the suggested lists
        - atom_to_insert  (str) - element name of atom to insert
        - atom_to_move (str) - element name of atom to move
        - it_new_folder or it_folder  (str) - section folder
        - inherit_option (str) - passed only to add_loop
        - inherit_magmom (bool) - if True than magmom from starting_calc is used, else from set

        - end_pos_types_z (list of int) - list of Z - type of atoms, which could be considered as final positions in vacancy creation mode

        - calc_method (list)
            - 'neb'
            - 'only_neb' - run only footer

        - x_start, x_final (array) - explicit xcart coordinates of moving atom for starting and final positions, combined with atom_to_insert
        - xr_start, xr_final (array) - explicit xred
        - rep_moving_atom (str)- replace moving atom by needed atom - can be useful than completly different atom is needed. 

        - upload_vts (bool) - if True upload Vasp.pm and nebmake.pl to server
        - run (bool)  - run on server
        - init_neb_geo_fld - folder with initial geometry files in VASP format, automatically switch on 'external' mode 

        - old_behaviour (str) - choose naming behavior before some date in the past for compatibility with your projects
            '020917'
            '261018' - after this moment new namig convention applied if end_pos_types_z is used

        - add_loop_dic - standart parameters of add()
        - params (dic) - provide additional parameters to add() # should be removed

    RETURN:
        None

    DEPENDS:

    TODO
    1. Take care of manually provided i_atom_to_move in case of replicate flag using init_numbers 
    2. For search_type == None x_m and x_del should be determined for magnetic searching and for saving their coordinates
    to struct_des; now their just (0,0,0) 


    """
    naming_conventions209 = True # set False to reproduce old behavior before 2.09.2017
    if old_behaviour == '020917':
        naming_conventions209 = False #


    # print('atom_to_insert', atom_to_insert)
    # sys.exit()

    calc = header.calc
    struct_des = header.struct_des
    varset = header.varset
    
    if not add_loop_dic:
        add_loop_dic = {}

    if not end_pos_types_z:
        end_pos_types_z = []
        end_pos_types_z = sorted(end_pos_types_z)


    if not hasattr(calc_method, '__iter__'):
        calc_method = [calc_method]


    if init_neb_geo_fld:

        printlog('I use folder with input files for neb in VASP format:', init_neb_geo_fld, imp = 'y')
        if st or starting_calc:
            printlog('Warning! st and starting_calc are ignored')
        neb_folders = [ f.name for f in os.scandir(init_neb_geo_fld) if f.is_dir() ]
        images = len(neb_folders) - 2
        search_type = 'external'

    elif starting_calc and st:
        printlog('Warning! both *starting_calc* and *st* are provided. I use *starting_calc*')
        st = copy.deepcopy(starting_calc.end)

    elif starting_calc:
        st = copy.deepcopy(starting_calc.end)
        printlog('I use *starting_calc*')


    elif st:
        ''
        printlog('I use *st*')

    else:
        printlog('Error! no input structure. Use either *starting_calc* or *st*')

    corenum = add_loop_dic.get('corenum')
    # print(corenum)
    # sys.exit()

    if corenum == None:
        if images == 3:
            corenum = 15
        elif images == 5:
            corenum = 15
        elif images == 7:
            corenum = 14
        else:
            printlog('add_neb(): Error! number of images', images,'is unknown to me; please provide corenum!')

    # print(corenum)
    # sys.exit()


    # print(atom_to_insert)
    # sys.exit()


    if corenum:
        # header.corenum = corenum
        ''
    else:
        corenum = header.CORENUM

    if corenum % images > 0:
        print_and_log('Error! Number of cores should be dividable by number of IMAGES', images, corenum)


    if not ise_new:
        ise_new = starting_calc.id[1]
        printlog('I use', ise_new, 'as ise_new', imp = 'y')


    name_suffix = ''
    st_pores = []

    name_suffix+='n'+str(images)


    """Replicate cell """
    if replicate:
        print_and_log('You have chosen to replicate the structure by', replicate)

        st = replic(st, mul = replicate)
        name_suffix += str(replicate[0])+str(replicate[1])+str(replicate[2])



    printlog('Search type is ', search_type, imp = 'y')
    
    if search_type == 'external':
        from siman.calc_manage import smart_structure_read
        st1 = smart_structure_read(init_neb_geo_fld+'/'+neb_folders[0]+'/POSCAR')
        st2 = smart_structure_read(init_neb_geo_fld+'/'+neb_folders[-1]+'/POSCAR')

        x_m = (0,0,0) # please improve, temporary
        x_del = (0,0,0)
        it = None
        obtained_from = 'from external folder '+init_neb_geo_fld




    elif search_type == None:
        
        if st_end == None:
            printlog('Error! You have provided search_type == None, st_end should be provided!')

        st1 = st
        st2 = st_end

        x_m = (0,0,0)
        x_del = (0,0,0)




    else:

        """1. Choose  atom (or insert) for moving """

        if is_list_like(xr_start):
            x_start = xred2xcart([xr_start], st.rprimd)[0]
            # print('atom_to_insert', atom_to_insert)
            # sys.exit()

            st1, i_m = st.add_atoms([x_start], atom_to_insert, return_ins = 1)
            x_m = x_start
            # i_m = st1.find_atom_num_by_xcart(x_start)
            # print(st1.get_elements()[i_m])
            # sys.exit()


            if i_atom_to_move:
                nn = str(i_atom_to_move+1)
            else:
                nn = str(i_void_start)


            name_suffix+=atom_to_insert+nn
            write_xyz(st1, file_name = st.name+'_manually_start')
            printlog('Start position is created manually by adding xr_start', xr_start, x_start)
            type_atom_to_move = atom_to_insert
            el_num_suffix = ''


        else:

            atoms_to_move = []
            atoms_to_move_types = []
            
            # print('d', i_atom_to_move)
            # sys.exit()

            if i_atom_to_move:
                typ = st.get_elements()[i_atom_to_move]
                printlog('add_neb(): atom', typ, 'will be moved', imp = 'y')
                atoms_to_move.append([i_atom_to_move, typ, st.xcart[i_atom_to_move]])
                atoms_to_move_types.append(typ)
    
                if naming_conventions209:
                    name_suffix+=typ+str(i_atom_to_move+1)



            else:
                #try to find automatically among alkali - special case for batteries
                for i, typ, x in zip(range(st.natom), st.get_elements(), st.xcart): 
                    if typ in ['Li', 'Na', 'K', 'Rb', 'Mg']:
                        atoms_to_move.append([i, typ, x])
                        if typ not in atoms_to_move_types:
                            atoms_to_move_types.append(typ)

            if  atoms_to_move:
                # print(atom_to_move)
                # sys.exit()
                if not atom_to_move:
                    atom_to_move = atoms_to_move_types[0] # taking first found element
                    if len(atoms_to_move_types) > 1:
                        printlog('Error! More than one type of atoms available for moving detected', atoms_to_move_types,
                            'please specify needed atom with *atom_to_move*')

                type_atom_to_move = atom_to_move #atoms_to_move[0][1]

                # printlog('atom ', type_atom_to_move, 'will be moved', imp ='y')

                if i_atom_to_move:
                    printlog('add_neb(): *i_atom_to_move* = ', i_atom_to_move, 'is used', imp ='y')
                    numbers = [[i_atom_to_move]]
                    i_void_start = 1
                else:
                    printlog('add_neb(): determine_symmetry_positions ...', imp ='y')

                    numbers = determine_symmetry_positions(st, atom_to_move)

                # print(numbers)
                # sys.exit()
                if len(numbers)>0:
                    printlog('Please choose position using *i_void_start* :', [i+1 for i in range(len(numbers))],imp = 'y' )
                    printlog('*i_void_start* = ', i_void_start)
                    i_m = numbers[i_void_start-1][0]
                    printlog('Position',i_void_start,'chosen, atom:', i_m+1, type_atom_to_move, imp = 'y' )
                
                else:
                    i_m = numbers[0][0]

                
                x_m = st.xcart[i_m]

                el_num_suffix =  type_atom_to_move +str(i_m+1)
                atom_to_insert = atom_to_move

                st1 = st
            # elif atom_to_replace:
            #     num = st.get_specific_elements(atom_to_replace)

            #     if len(n)>0:
            #         printlog('Please choose position using *i_void_start* :', [i+1 for i in range(len(num))],imp = 'y' )
            #         printlog('*i_void_start* = ', i_void_start)
            #         i_m = num[i_void_start-1]
            #         printlog('Position',i_void_start,'chosen, atom to replace:', i_m+1, atom_to_replace, imp = 'y' )
            #         sys.exit()


            else:

                print_and_log('No atoms to move found, you probably gave me deintercalated structure', important = 'y')
                
                st_pores, sums, avds = determine_voids(st, r_impurity, step_dec = 0.1, fine = 2)
              

                insert_positions = determine_unique_voids(st_pores, sums, avds)

                print_and_log('Please use *i_void_start* to choose the void for atom insertion from the Table above:', 
                    end = '\n', imp = 'Y')


                if i_void_start == None:
                    sys.exit()
                if atom_to_insert == None:
                    printlog('Error! atom_to_insert = None')

                st = st.add_atoms([insert_positions[i_void_start],], atom_to_insert)

                name_suffix+='i'+str(i_void_start)

                i_m = st.natom-1
                x_m = st.xcart[i_m]


                search_type = 'existing_voids'
                type_atom_to_move = atom_to_insert
                el_num_suffix = ''

                st1 = st



        """2. Choose final position"""

        if is_list_like(xr_final):
            x_final = xred2xcart([xr_final], st.rprimd)[0]
            

            #old
            #check if i_atom_to_move should be removed
            # st2 = st1.del_atom(i_m)
            # st2 = st2.add_atoms([x_final], atom_to_insert)

            #new
            st2 = st1.mov_atoms(i_m, x_final)


            
            # st1.printme()
            # st2.printme()
            # sys.exit()


            x_del = x_final 
            search_type = 'manual_insertion'
            name_suffix+= 'v'+str(i_void_final) 
            write_xyz(st2, file_name = st.name+'_manually_final')
            printlog('Final position is created manually by adding xr_final', xr_final, x_del)



        elif search_type == 'existing_voids':
            #Search for voids around choosen atoms

            if not st_pores: 
                st_pores, sums, avds = determine_voids(st, r_impurity)
            
            sur = determine_unique_final(st_pores, sums, avds, x_m)

            print_and_log('Please choose *i_void_final* from the Table above:', end = '\n', imp = 'Y')


            if i_void_final == None:
                sys.exit()

            x_final = sur[0][i_void_final] #
            
            printlog('You chose:', np.array(x_final).round(2), end = '\n', imp = 'Y')


            x_del = x_final #please compare with vacancy creation mode

            write_xyz(st.add_atoms([ x_final], 'H'), replications = (2,2,2), file_name = st.name+'_possible_positions2_replicated')
            
            print_and_log('Choosing the closest position as end', important = 'n')

            st1 = st

            st2 = st.mov_atoms(i_m, x_final)
            
            name_suffix += el_num_suffix+'e'+str(i_void_final)+atom_to_insert

            st1 = return_atoms_to_cell(st1)
            st2 = return_atoms_to_cell(st2)

            write_xyz(st1, file_name = st1.name+name_suffix +'_start')

            write_xyz(st2, file_name = st2.name+name_suffix +'_final')


        elif search_type == 'vacancy_creation':
            #Create vacancy by removing some neibouring atom of the same type 
            
            print_and_log('You have chosen vacancy_creation mode of add_neb tool', imp= 'Y')

            print_and_log( 'Type of atom to move = ', type_atom_to_move, imp = 'y')
            # print 'List of left atoms = ', np.array(st.leave_only(type_atom_to_move).xcart)

            final_pos_z = end_pos_types_z or [invert(type_atom_to_move)] # by default only moving atom is considered
            end_pos_types_el = [invert(z) for z in end_pos_types_z]


            sur = local_surrounding(x_m, st, n_neighbours = 14, control = 'atoms', 
                only_elements = final_pos_z,
                periodic  = True) #exclude the atom itself

            # print(x_m)
            # print(sur)

            # st.nn()
            end_pos_n = sur[2][1:]
            print_and_log(
            'I can suggest you '+str (len(end_pos_n) )+' end positions. The distances to them are : ',np.round(sur[3][1:], 2), ' A\n ',
            'They are ', [invert(z) for z in final_pos_z], 'atoms, use *i_void_final* to choose required: 1, 2, 3 ..', imp = 'y')

            i_sym_final_l = []
            for j in end_pos_n:
                for i, l in enumerate(numbers):
                    if j in l: 
                        i_sym_final_l.append(i+1)
            printlog('Their symmetry positions are ', i_sym_final_l, imp = 'y')

            # sys.exit()


            if not i_void_final:
                printlog('Changing i_void_final: None -> 1', imp = 'y')
                i_void_final = 1 #since zero is itself
            chosen_dist = sur[3][i_void_final]
            print_and_log('Choosing position ', i_void_final, 'with distance', round(chosen_dist, 2), 'A', imp = 'y')

            # print(end_pos_n)
            i_sym_final = 0
            n_final = sur[2][i_void_final]
            for i, l in enumerate(numbers):
                if n_final in l:
                    i_sym_final = i+1
            printlog('It is symmetrically non-equiv position #', i_sym_final, imp = 'y')

            # sys.exit()

            header.temp_chosen_dist = chosen_dist


            if old_behaviour == '261018':
                name_suffix += el_num_suffix+'v'+str(i_void_final)
            else:

                name_suffix += el_num_suffix+'v'+str(i_void_final)+list2string(end_pos_types_el, joiner = '')

                # print(name_suffix)
                # sys.exit()


            x_del = sur[0][i_void_final]
            printlog('xcart of atom to delete', x_del)
            i_del = st.find_atom_num_by_xcart(x_del)
            # print(x_del)
            # print(st.xcart)
            # for x in st.xcart:
            #     if x[0] > 10:
            #         print(x)


            print_and_log( 'number of atom to delete = ', i_del, imp = 'y')
            if i_del == None:
                printlog('add_neb(): Error! I could find atom to delete!')

            # print st.magmom
            # print st1.magmom


            # try:
            if is_list_like(xr_start):
                st2 = st1.mov_atoms(i_m, x_del) # i_m and sur[0][neb_config] should coincide
                # i_del = st1.find_atom_num_by_xcart(x_del)
                
                st1 = st1.del_atom(i_del) 

            else:
                print_and_log('Making vacancy at end position for starting configuration', imp = 'y')
                st1 = st.del_atom(i_del) 


                print_and_log('Making vacancy at start position for final configuration', important = 'n')
                st2 = st.mov_atoms(i_m, x_del) # i_m and sur[0][neb_config] should coincide
            # except:
                # st2 = st

            st2 = st2.del_atom(i_del) # these two steps provide the same order



    """Checking correctness of path"""
    #if start and final positions are used, collisions with existing atoms are possible 
    if is_list_like(xr_start) and is_list_like(xr_final):
        printlog('Checking correctness')
        st1, _, _ = st1.remove_close_lying()

        stt = st1.add_atoms([x_final,], 'Pu')
        stt, x, _ = stt.remove_close_lying(rm_both = True) # now the final position is empty for sure; however the order can be spoiled
        # print(st._removed)
        if stt._removed:
            st1 = stt # only if overlapping was found we assign new structure


        st2, _, _ = st2.remove_close_lying(rm_first = stt._removed)
        stt = st2.add_atoms([x_start,], 'Pu')
        stt, x, _ = stt.remove_close_lying(rm_both = True) # now the start position is empty for sure
        if stt._removed:
            st2 = stt

        print(st2.get_elements())
        # sys.exit()


    elif is_list_like(xr_final) and not is_list_like(xr_start) or is_list_like(xr_start) and not is_list_like(xr_final):
        printlog('Attention! only start or final position is provided, please check that everything is ok with start and final states!!!')





    """ Determining magnetic moments  """
    vp = varset[ise_new].vasp_params



    if 'ISPIN' in vp and vp['ISPIN'] == 2:
        print_and_log('Magnetic calculation detected. Preparing spin modifications ...', imp = 'y')
        cl_test = CalculationVasp(varset[ise_new])
        cl_test.init = st1
        # print 'asdfsdfasdfsadfsadf', st1.magmom
        if inherit_magmom and hasattr(st, 'magmom') and st.magmom and any(st.magmom):
            print_and_log('inherit_magmom=True: You have chosen MAGMOM from provided structure', imp = 'y')
            name_suffix+='mp' #Magmom from Previous
        else:
            cl_test.init.magmom = None
            print_and_log('inherit_magmom=False or no magmom in input structure : MAGMOM will be determined  from set', imp = 'y')
            name_suffix+='ms' #Magmom from Set


        cl_test.actualize_set() #find magmom for current structure

        st1.magmom = copy.deepcopy(cl_test.init.magmom)
        st2.magmom = copy.deepcopy(cl_test.init.magmom)

        # sys.exit()
        # print_and_log('The magnetic moments from set:')
        # print cl_test.init.magmom
        if search_type not in [None, 'external']:  # for None not implemented; x_m should be determined first for this
            #checking for closest atoms now only for Fe, Mn, Ni, Co
            sur   = local_surrounding(x_m, st1, n_neighbours = 3, control = 'atoms',
            periodic  = True, only_elements = header.TRANSITION_ELEMENTS)

            dist = np.array(sur[3]).round(2)
            numb = np.array(sur[2])
            a = zip(numb, dist )
            # a=  np.array(a)
            # print a[1]
            # a = np.apply_along_axis(np.unique, 1, a)
            # print a
            def unique_by_key(elements, key=None):
                if key is None:
                    # no key: the whole element must be unique
                    key = lambda e: e
                return list ( {key(el): el for el in elements}.values() )

            # print a
            mag_atoms_dists = unique_by_key(a, key=itemgetter(1))
            # print (mag_atoms_dists)
            # a = unique_by_key(a, key=itemgetter(1))
            print_and_log( 'I change spin for the following atoms:\ni atom     dist\n', np.round(mag_atoms_dists, 2) , imp = 'y' )
            # print 'I have found closest Fe atoms'
            muls = [(1.2, 0.6), (0.6, 1.2)]
            mag_moments_variants = []
            for mm in muls:
                mags = copy.deepcopy(cl_test.init.magmom)
                # print mags
                for a, m in zip(mag_atoms_dists, mm):
                    # print t[1]
                    mags[a[0]] = mags[a[0]]*m
                mag_moments_variants.append(mags)

            print_and_log( 'The list of possible mag_moments:', imp = 'y' )
            for i, mag in enumerate(mag_moments_variants):
                print_and_log( i, mag)

            print_and_log( 'Please use *mag_config* arg to choose desired config' , imp = 'y' )


        if mag_config != None:

            st1.magmom = copy.deepcopy(mag_moments_variants[mag_config])
            st2.magmom = copy.deepcopy(mag_moments_variants[mag_config])

            name_suffix+='m'+str(mag_config)

            print_and_log('You have chosen mag configuration #',mag_config,imp = 'y')

    else:
        print_and_log('Non-magnetic calculation continue ...')















    """3. Add to struct_des, create geo files, check set, add_loop """

    if starting_calc:
        it = starting_calc.id[0]
        it_new = it+'v'+str(starting_calc.id[2])+'.'+name_suffix

        if not it_new_folder:
            it_new_folder = struct_des[it].sfolder+'/neb/'
        obtained_from = str(starting_calc.id) 

        if not ise_new:
            print_and_log('I will run add_loop() using the same set', important = 'Y')
            ise_new = cl.id[1]

    elif st:
        if not it_new:
            printlog('Error! please provide *it_new* - name for your calculation', important = 'Y')


        it = None
        it_new+='.'+name_suffix
        obtained_from = st.name
        
        if not ise_new:
            printlog('Error! please provide *ise_new*', important = 'Y')

        if not it_new_folder and not it_folder:
            
            printlog('Error! please provide *it_new_folder* - folder for your new calculation', important = 'Y')
        if it_folder:
            it_new_folder = it_folder

    if rep_moving_atom:
        it_new += 'r'+rep_moving_atom

    if it_new not in struct_des:
        add_des(struct_des, it_new, it_new_folder, 'Automatically created and added from '+obtained_from  )




    print_and_log('Creating geo files for starting and final configurations (versions 1 and 2) ', important = 'y')

    # if starting_calc:
    #     cl = copy.deepcopy(starting_calc)
    # else:

    cl = CalculationVasp()

    #write start position
    if search_type  not in [None, 'external']:
        struct_des[it_new].x_m_ion_start = x_m
        struct_des[it_new].xr_m_ion_start = xcart2xred([x_m], st1.rprimd)[0]

        # st1, _, _ = st1.remove_close_lying()
        # st2, _, _ = st2.remove_close_lying()
        print('Trying to find x_m', x_m)
        i1 = st1.find_atom_num_by_xcart(x_m, prec = 0.45, )

        # sys.exit()
        print('Trying to find x_del', x_del)

        i2 = st2.find_atom_num_by_xcart(x_del, prec = 0.45, )

        if rep_moving_atom: #replace the moving atom by required
            st1 = st1.replace_atoms([i1], rep_moving_atom)
            st2 = st2.replace_atoms([i2], rep_moving_atom)
        else:
            #allows to make correct order for nebmake.pl
            st1 = st1.replace_atoms([i1], type_atom_to_move)
            st2 = st2.replace_atoms([i2], type_atom_to_move)

        i1 = st1.find_atom_num_by_xcart(x_m, prec = 0.45) # the positions were changed # check if this is correct
        i2 = st2.find_atom_num_by_xcart(x_del, prec = 0.45)


    cl.end = st1
    ver_new = 1
    cl.version = ver_new
    cl.path["input_geo"] = header.geo_folder + struct_des[it_new].sfolder + '/' + \
        it_new+"/"+it_new+'.auto_created_starting_position_for_neb_'+str(search_type)+'.'+str(ver_new)+'.'+'geo'
    
    cl.write_siman_geo(geotype = 'end', description = 'Starting conf. for neb from '+obtained_from, override = True)


    #write final position

    struct_des[it_new].x_m_ion_final = x_del
    struct_des[it_new].xr_m_ion_final = xcart2xred([x_del], st2.rprimd)[0]



    cl.end = st2
    ver_new = 2
    cl.version = ver_new
    cl.path["input_geo"] = header.geo_folder + struct_des[it_new].sfolder + '/' + \
        it_new+"/"+it_new+'.auto_created_final_position_for_neb_'+str(search_type)+'.'+str(ver_new)+'.'+'geo'
    
    cl.write_siman_geo(geotype = 'end', description = 'Final conf. for neb from '+obtained_from, override = True)

    if not rep_moving_atom and search_type  not in [None, 'external']:
        st1s = st1.replace_atoms([i1], 'Pu')
        st2s = st2.replace_atoms([i2], 'Pu')
    else:
        st1s = copy.deepcopy(st1)
        st2s = copy.deepcopy(st2)
    
    if center_on_moving and search_type  not in [None, 'external']:

        vec = st1.center_on(i1)
        st1s = st1s.shift_atoms(vec)
        st2s = st2s.shift_atoms(vec)
        write_xyz(st1s, file_name = it_new+'_start')
        write_xyz(st2s, file_name = it_new+'_end')


    st1s.write_poscar('xyz/POSCAR1')
    st2s.write_poscar('xyz/POSCAR2')
    # print(a)
    # runBash('cd xyz; mkdir '+it_new+'_all;'+"""for i in {00..04}; do cp $i/POSCAR """+ it_new+'_all/POSCAR$i; done; rm -r 00 01 02 03 04')
    
    with cd('xyz'):
        if init_neb_geo_fld:
            printlog('Check neb files at', init_neb_geo_fld, imp = 'y')
        else:
           a = runBash(header.PATH2NEBMAKE+' POSCAR1 POSCAR2 3')
           print(a)
           dst = it_new+'_all'
           makedir(dst+'/any')
           for f in ['00', '01', '02', '03', '04']:
              shutil.move(f+'/POSCAR', dst+'/POSCAR'+f)
              shutil.rmtree(f)





    #prepare calculations
    # sys.exit()








    #Check if nebmake is avail
    # if int(runBash('ssh '+cluster_address+' test -e '+project_path_cluster+'/tools/vts/nebmake.pl; echo $?') ):

    #     ''
    #     print_and_log('Please upload vtsttools to ',cluster_address, project_path_cluster+'/tools/vts/')
    #     raise RuntimeError

    #     copy_to_server(path_to_wrapper+'/vtstscripts/nebmake.pl', to = project_path_cluster+'/tools/',  addr = cluster_address)
    # if  int(runBash('ssh '+cluster_address+' test -e '+project_path_cluster+'/tools/Vasp.pm; echo $?') ):
    #     copy_to_server(path_to_wrapper+'/vtstscripts/Vasp.pm', to = project_path_cluster+'/tools/',  addr = cluster_address)





    inherit_ngkpt(it_new, it, varset[ise_new])

    if run:
        add_loop_dic['run'] = run

    add_loop_dic['corenum'] = corenum

    if 'params' not in add_loop_dic:
        add_loop_dic['params'] = {}

    if init_neb_geo_fld:
        add_loop_dic['params']['init_neb_geo_fld'] = init_neb_geo_fld
    # print(add_loop_dic)
    add_loop(it_new, ise_new, verlist = [1,2], up = up, calc_method = calc_method, savefile = 'oc', 
        inherit_option = inherit_option, n_neb_images = images, 
        **add_loop_dic  )
    

    if upload_vts:
        siman_dir = os.path.dirname(__file__)
        # print(upload_vts)
        push_to_server([siman_dir+'/cluster_tools/nebmake.pl', siman_dir+'/cluster_tools/Vasp.pm'], to = header.cluster_home+'/tools/vts',  addr = header.cluster_address)
    
    else:
        print_and_log('Please be sure that vtsttools are at',header.cluster_address, header.cluster_home+'/tools/vts/', imp = 'Y')

    printlog('add_neb finished')
    return it_new 

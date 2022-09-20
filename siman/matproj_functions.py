import os, sys

from siman import header

from siman.header import db
from siman.geo import  create_surface2, stoichiometry_criteria
from siman.calc_manage import  get_structure_from_matproj, smart_structure_read
from siman.calc_manage   import (clean_history_file, prepare_run,  manually_remove_from_struct_des, update_des, inherit_icalc, add_loop, res_loop, complete_run, add_des )

try:
    # pmg config --add VASP_PSP_DIR $VASP_PSP_DIR MAPI_KEY $MAPI_KEY
    from pymatgen.ext.matproj import MPRester
    from pymatgen.io.vasp.inputs import Poscar
    from pymatgen.io.cif import CifParser
    pymatgen_flag = True 
except:
    print('matproj_functions.py: pymatgen is not available')
    pymatgen_flag = False 

from siman.classes import CalculationVasp


import csv
from siman.analysis import suf_en


# Author: A.Boev


###############################################3
### read and write pmg csv file
#############################################
def read_matproj_info(path):
    """
    """

    data = []
    with open(path+'.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile)

        for d in reader:
            data.append(d)

    return data

# def read_matproj_info_1(path):
#     """
#     Addition of MP_Compound() objects to db
#     """
#     from siman.classes import MP_Compound
#     # data = []
#     compound_list = []
#     with open(path+'.csv', newline='') as csvfile:
#         reader = csv.DictReader(csvfile)

#         for d in reader:
#             mp_name = d['pretty_formula']+'.MP'
#             mp = MP_Compound()
#             mp.name = mp_name
#             mp.pretty_formula = d['pretty_formula']
#             mp.material_id = d['material_id']
#             mp.elements = d['elements']
#             mp.sg_symbol = d['spacegroup.symbol']
#             mp.sg_crystal_str = d['spacegroup.crystal_system']
#             mp.band_gap = d['band_gap']
#             mp.price_per_gramm = d['price_per_gramm']
#             mp.total_magnetization = d['total_magnetization']
            
#             db[mp_name] = mp
#             # print(db[mp_name].pretty_formula)
#             compound_list.append(mp_name)

#     return compound_list


def read_matproj_list(path, new_list = 0, new_object = 0):
    """
    Addition of MP_Compound() objects to db
    """
    from siman.classes import MP_Compound
    # data = []
    compound_list = []
    with open(path+'.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile)

        for d in reader:
            mp_name = d['pretty_formula']+'.MP'
            compound_list.append(mp_name)
            if new_list:
                if new_object:
                    mp = MP_Compound()
                else:
                    mp = db[mp_name]
                mp.name = mp_name
                mp.pretty_formula = d['pretty_formula']
                mp.material_id = d['material_id']
                mp.elements = d['elements']
                mp.sg_symbol = d['spacegroup.symbol']
                mp.sg_crystal_str = d['spacegroup.crystal_system']
                mp.band_gap = d['band_gap']
                try:
                    mp.price_per_gramm = d['price_per_gramm']
                except KeyError:
                    mp.price_per_gramm = None
                mp.total_magnetization = d['total_magnetization']
                mp.formation_energy_per_atom = d['formation_energy_per_atom']
                mp.e_above_hull = d['e_above_hull']
                try:
                    mp.icsd_ids = d['icsd_ids']
                except KeyError:
                    mp.icsd_ids = None
                db[mp_name] = mp



    return compound_list



def write_matproj_info(data, path, properties):
    """
    """
    with open(path+'.csv', 'w', newline='') as csvfile:
            fieldnames = properties
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for st_name in sorted(data.keys()):
                d_i = data[st_name]
                # print(d_i)
                for ss in d_i.keys():
                    if type(d_i[ss]) == float:
                        d_i[ss] = format(d_i[ss],'.2f')

                writer.writerow(d_i)



def write_MP_compound(compound_list, path, properties):
    """
    """
    if not properties:
        properties = ['material_id', 'pretty_formula', 'sg_symbol', 'sg_crystal_str', 'formation_energy_per_atom',  
                    'band_gap', 'total_magnetization', 'e_above_hull', 'price_per_gramm', 'bulk_status_scale', 'e_cohesive',  'e_cohesive_MP', 'ec_es',  'suf_en', 'icsd_ids']


    with open(path+'.csv', 'w', newline='') as csvfile:
            fieldnames = properties
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            string = {}
            for st_name in compound_list:
                # print(dir(db[st_name]))
                for i in properties:
                    try:
                        string[i] = eval('db[st_name].'+i)
                    except AttributeError:
                        None
                # print(string)


                # d_i = data[st_name]
                # # print(d_i)
                # for ss in d_i.keys():
                #     if type(d_i[ss]) == float:
                #         d_i[ss] = format(d_i[ss],'.2f')

                writer.writerow(string)



#############################################
#############################################


def get_matproj_info(criteria, properties, filename = None, price = 0, element_price = None, only_stable = 0, exclude_list = [], all_phases = 0):
    """
    Retrieve a list of structures with information from materials project database (MatProj) according to a set of constraints

    INPUT:
        criteria (dict)   - dictionary with conditions for query method to choose required structures
                            see examples for detailed description of this dict
        properties (list) - list of retrieved properties (MatProj keys) for each structure; e.g. m.get_data('mp-12957') 
        filename (str)    - name of file with retrieved structures and info
        price (bool)      - use prices as one of the constraints
        element_price (dict) - dictionary with prices of elements {element: price per kg}
        only_stable (bool)   - retrievely only experimentally stable compounds from MatProj database (highlighted with green in web interface)
        exclude_list (list)  - list of compound pretty formulas, which won't be included into the final list of structures
        all_phases (bool)    - let to write down all the polymorphs for each composition (in test)
    

    SIDE:
        writes file *filename*.cvs

    RETURN:
        data (dict) - dictionaly with retrieved structures


    """

    import pymatgen
   
    data = {}

    if 'e_above_hull' not in properties:
        properties.append('e_above_hull')


    with MPRester(header.pmgkey) as m:
            properties.append('elements')

            results = m.query(criteria= criteria, properties=properties)

            for string_to_write in results:


                if string_to_write['pretty_formula'] not in exclude_list:

                    if price:

                        if 'price_per_gramm' not in properties:
                            properties.append('price_per_gramm')
                        cost = 0
                        for p in string_to_write['elements']:
                            if p in element_price.keys():
                                # print(element_price[p])
                                cost+= element_price[p]/1000
                            else:
                                print('Warning! Element',p, 'price was not found')
                                cost = 'No data'
                                break


                        string_to_write.update([('price_per_gramm', cost)])



                    st_name = string_to_write['pretty_formula']

                    if not only_stable:

                        if st_name not in data.keys():         
                            data.update([(st_name, string_to_write)])
                        else:
                            if all_phases:
                                n = list(data.keys()).count(st_name)
                                st_name = st_name + '.' + str(n)
                                data.update([(st_name, string_to_write)])
                            else:

                                if string_to_write['e_above_hull'] < data[st_name]['e_above_hull']:
                                    data.update([(st_name, string_to_write)])
                    else:
                        try:
                            entries = m.get_entries_in_chemsys(string_to_write['elements'])
                            print(st_name)

                            for i in m.get_stability(entries):

                                if i['entry_id'] in i['decomposes_to'][0]['material_id']:
                                    data.update([(st_name, string_to_write)])
                                    # print(i)
                        except pymatgen.ext.matproj.MPRestError:
                            properties.append('warnings')
                            string_to_writeupdate([('warnings', 'Stability MPRestError')]) 
                            data.update([(st_name, string_to_write)])

    # write data in csv file
    if filename:
        write_matproj_info(data, filename, properties)
        

    return data



def get_matproj_chem_pot(atom_list =None):
    """
    """

    at_mol = ['H', 'N', 'O', 'F','Cl', 'Br', 'I', 'At']
    atom_list_n = []
    if atom_list:
        for i in atom_list:
            if i in at_mol:
                i+='2'
            atom_list_n.append(i)




    data_out = {}

    if atom_list:
        criteria_string = {'elements': {'$in': atom_list},'nelements':1,'formation_energy_per_atom':{'$lt': 1e-8}}
    else:
        criteria_string = {'nelements':1,'formation_energy_per_atom':{'$lt': 1e-8}}

    properties = ['pretty_formula', 'energy_per_atom']



    data_mp = get_matproj_info(criteria_string, properties, name = 0, price = 0)  #list of pure elements

    for i in data_mp.keys():

        data_out[i] = data_mp[i]['energy_per_atom']

    if len(atom_list_n) == 1:
        return data_out[atom_list_n[0]]
    
    else:
        return data_out







###################################################################################################
######################################################################################################3


def get_matproj_st(mat_in_list, folder = 'geo/'):
    """
    check downloaded POSCAR files in geo/ folder
    if not POSCAR of some structure - download it from Mat Proj

    mat_in_list - data dict for any structure from MP,  result of get_fata('mp-...')
    """
    
    name = mat_in_list['material_id']+'.POSCAR'
    if name not in os.listdir(folder):
        os.chdir(folder)
        st = get_structure_from_matproj(mat_proj_id = mat_in_list['material_id'], it_folder = folder)
        os.chdir('..')
    else:
        st = smart_structure_read(folder+name)
        # print('ok')
    return st




###################################################################################################
###################################################################################################
"""
For list calculations I use some nomenclature:


bulk cl has a name: compound.crystal_symmetry.set.version  - for example Ag2O.cubic.8.1
suf cl has a name: compound.crystal_symmetry.spacegroup_symbol.facet.facet_version.set.version  - for example Ag2O.cubic.Pn-3m.110.0.9s.1


magnetic sets have the same name with usually set with addition 'm':
'8' and '8m' static for bulk
'9' and '9m' at relax
'9s' and '9sm' at relax for suf
0mboxn for boxed atoms

"""



def calc_bulk_list(data_list, spacegroup = '', ise_nomag = '8', ise_mag = '8m', status = None, up =None, corenum = 1):
    """
    function can add or res for set of bulk calculations
    
    data_list = read_pmg_info() of some csv file

    ise  - set of calculation. Usually use '8' for nonmag calc and '8m' for mag calc
    status = add or res
    """

    spacegroup = '.'+spacegroup


    for i in data_list:
        st = get_matproj_st(i)


        if float(i['total_magnetization']): 
            mag_flag = 1
            ise = ise_mag
        else: 
            mag_flag = 0
            ise = ise_nomag




        if status == 'add': 
            add_loop(i['pretty_formula']+spacegroup, ise, 1, it_folder = 'bulk', input_st = st, corenum = corenum)
        if status == 'res': 

            try:

                res_loop(i['pretty_formula']+spacegroup, ise, 1, it_folder = 'bulk', up = up)


            except KeyError:
                print(i['total_magnetization'])
                # add_loop(i['pretty_formula']+spacegroup, ise, 1, it_folder = 'bulk', input_st = st, corenum = corenum)












def calc_suf_list_sg(data_list, sg, suf_list, flag = 0):

    """
    to run a set of slabs for every str from matproj data_list
    sg - 'spacegroup.symbol'
    """


    # print('start')
    for i in data_list:

        
        if i['spacegroup.symbol'] == sg:
            print(i['spacegroup.symbol'])
            st_name = i['pretty_formula']+'.'+i['spacegroup.crystal_system']
            print(st_name)

            if float(i['total_magnetization']): 
                mag_flag = 1
                ise = '9sm'
                st_bulk = db[st_name, '8m', 1].end
            else: 
                mag_flag = 0
                ise = '9s'
                st_bulk = db[st_name, '8', 1].end
            

            try:
                for surface in suf_list:
                    slabs = create_surface2(st_bulk, miller_index = surface,  min_slab_size = 10, min_vacuum_size = 10, surface_i = 0,
                      symmetrize = True)
                    
                    s_i = 0
                    if len(slabs):
                        for sl_i in slabs:
                            st = st_bulk
                            sl = st.update_from_pymatgen(sl_i)

                            if stoichiometry_criteria(sl, st_bulk):

                                if flag == 'add':
                                    add_loop(st_name+'.cubic.'+i['spacegroup.symbol']+'.'+ str(surface[0])+str(surface[1])+str(surface[2]) +'.'+str(s_i), ise, 1 , 
                                        input_st = sl,  it_folder = 'slab/'+st_name+'.cubic.'+i['spacegroup.symbol'], up = 'up2')
                                elif flag == 'res':
                                    res_loop(st_name+'.cubic.'+i['spacegroup.symbol']+'.'+ str(surface[0])+str(surface[1])+str(surface[2]) +'.'+str(s_i), ise, 1,  up = 'up2')
                                else: 
                                    print(st_name + '\n', surface, s_i+1, ' from ', len(slabs), ' slabs\n', sl.natom, ' atoms'  )
                                s_i+=1
                                if s_i == 3: break

                            else:
                                print('Non-stoichiometric slab')

                    else:
                        print('\nWarning!  Zero slabs constructed!\n')

            except AttributeError:
                print('\nWarning!  Bulk calculation of {} has some problems!!\n'.format(st_name)) 


#######################################################################
### Energy calculations
###########################################################################

def coh_en_list(it_bulk, ise_bulk = '8', ise_box = '0mboxn'):
        """
        Is waiting
        """

        cl_bulk = db[it_bulk,ise_bulk,1]#.jmol()
        e_bulk = cl_bulk.energy_sigma0
        n_at_sum = cl_bulk.end.natom
        el_list = cl_bulk.end.get_elements()

        e_at_sum = 0
        for el in el_list:
            it = '.'.join([el,'box'])
            cl_at = db[it, ise_box, 1]
            e_at = cl_at.energy_sigma0

            e_at_sum+=e_at

        e_coh = (e_at_sum-e_bulk)/n_at_sum
        print('{}  \t\tE_coh = {} eV'.format(it_bulk.split('.')[0], round(e_coh,1)))

        return e_coh



def suf_en_list(it_suf, it_bulk, ise_suf, ise_bulk, hkl_list):
    """
    function calculates surface enegries using suf_en()
    cl_name - pretty_formula
    sg - symmetry group
    mag- magnetic or not
    #

    return suf_en_pack = [surface_index, suf_en, max_f, max_f_bulk ]
    """

    suf_en_pack = []


    cl_bulk = db[it_bulk, ise_bulk, 1]



    for hkl in hkl_list:

        for i in range(0,8):
        
            it = '{}.{}.{}'.format(it_suf, hkl, str(i))
            
            if it in header.struct_des:

                # try:

                    from siman.geo import create_surface2
                    # slab = create_surface2(cl_bulk.end, miller_index = [int(hkl[0]),int(hkl[1]),int(hkl[2])],  min_slab_size = 10, min_vacuum_size = 10, 
                    #     surface_i = i, return_one = 1,
                    #               symmetrize = True)

                    try:
                        stoich_cr = stoichiometry_criteria(db[it,ise_suf,1].end, cl_bulk.end)
                    except (AttributeError, KeyError):
                        # db[it,ise_suf,1].res()
                        continue
                    if stoich_cr:


                    # try:
                        # db[it, ise_suf, 1].res()
                            gamma = suf_en(db[it, ise_suf, 1], cl_bulk, silent = 1)
                            max_f = str(max(db[it, ise_suf, 1].maxforce_list[-1]))
                            max_f_bulk = str(max(cl_bulk.maxforce_list[-1]))


                    # except KeyError:
                    else:
                        gamma = None
                        max_f = "Non_st"
                        max_f_bulk = str(max(cl_bulk.maxforce_list[-1]))


                    surface ='{}.{}'.format(hkl, str(i))

                    if not gamma:
                        gamma = 0.
                        # surface+='.Non-st'
                        # max_f = "None"
                        # max_f_bulk = "None"
                    print(surface, 'E_suf = ',round(gamma,1), ' eV')
                    suf_en_pack.append([surface, round(gamma,2), max_f, max_f_bulk ])

                # except AttributeError: 
                #     continue

    # print(suf_en_pack)
    if not suf_en_pack:
        suf_en_pack.append(['None',0,"None", 'None'])

    return suf_en_pack






def energy_list(data_list, sg_symbol = None, suf_en = 0, coh_en = 0, ise_bulk = '8', ise_box = '0mboxn', ise_suf = '9s', hkl_list = None):

    """
    Is waiting
    """
    if not hkl_list:
        hkl_list = ['110','111']



    out = {'coh_en': None, 'suf_en': None,}
    coh_en_f = []
    suf_en_f = []

    for i in data_list:

        mat_i = i['pretty_formula']
        sg_i = i['spacegroup.crystal_system']
        sg_sym_i = i['spacegroup.symbol']
        ise_b=ise_bulk
        ise_s=ise_suf

        if not sg_symbol:
            stat = 1
        elif sg_sym_i == sg_symbol:
            stat = 1
        else:
            stat = 0


        if stat:            
            print('\n\n','.'.join([mat_i,sg_i,sg_sym_i] ),'\n')

            if float(i['total_magnetization']) > 1e-8:
                ise_b+='m'
                ise_s+='m'


            name_bulk = '.'.join([mat_i,sg_i])
            name_suf  = '.'.join([mat_i,sg_i,sg_sym_i])

            if coh_en:
                e_c = coh_en_list(it_bulk = name_bulk, ise_bulk = ise_bulk, ise_box = ise_box)
                coh_en_f.append(e_c)

            if suf_en:
                e_s = suf_en_list(it_suf = name_suf, it_bulk = name_bulk, ise_suf = ise_suf, ise_bulk = ise_bulk, hkl_list = hkl_list)
                suf_en_f.append(e_s)

    out['coh_en'] = coh_en_f
    out['suf_en'] = suf_en_f
    return out





#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################

def calc_suf_mat(self, ise = '9sm', bulk_cl_name = None,  it_folder = 'bulk/', status = 'create', suf_list = None, 
            min_slab_size = 12, only_stoich = 1, symmetrize = True, corenum = 4, cluster = 'cee', conv = 1, 
            create = 1, reset_old = 0, suffix = None): 
    '''
    '''
    from siman.header import db
    from siman.geo import create_surface2, stoichiometry_criteria,stoichiometry_criteria2
    from siman.calc_manage   import add_loop, res_loop


    if not bulk_cl_name:
        bulk_cl_name = ['it','ise', '1']
        
    if not suf_list:
        suf_list = [[1,1,0], [1,1,1]]


    st_bulk = db[self.bulk_cl_scale].end

    if conv:
        st_bulk = st_bulk.get_conventional_cell()

    st_name = '.'.join([self.pretty_formula, self.sg_crystal_str, self.sg_symbol])
    st_name = st_name.replace('/','_')
    # st_bulk.printme()
    try:
        # print(self.suf)
        self.suf.keys()
    except AttributeError:
        self.suf = {}
        self.suf_status = {}

    try:
        self.suf_cl.keys()
    except AttributeError:
        self.suf_cl = {}

    try:
        self.suf_en.keys()
    except AttributeError:
        self.suf_en = {}

    if reset_old:
        self.suf = {}
        self.suf_status = {}
        self.suf_cl = {}
        self.suf_en = {}




    ##################### create slab
    if status == 'create':


        try:
            for surface in suf_list:
                print(surface)

                if self.suf == {} or ''.join([str(surface[0]),str(surface[1]),str(surface[2])]) not in self.suf.keys():

                    slabs = create_surface2(st_bulk, miller_index = surface,  min_slab_size = min_slab_size, min_vacuum_size = 15, surface_i = 0,
                      symmetrize = symmetrize, lll_reduce = 1, primitive = 1)
                    for sl_i in slabs:
                                st = st_bulk
                                sl = st.update_from_pymatgen(sl_i)
                                stoichiometry_criteria2(st,sl)
                    s_i = 0
                    if create:
                        if len(slabs):
                            for sl_i in slabs:
                                st = st_bulk
                                sl = st.update_from_pymatgen(sl_i)
                                stoichiometry_criteria2(st,sl)
                                suf_name = str(surface[0])+str(surface[1])+str(surface[2]) +'.'+str(s_i)

                                if self.suf_status == 'Zero slabs constructed':
                                    self.suf_status = {}

                                self.suf[suf_name] = sl
                                self.suf_status[suf_name] = 'created'
                                print(st_name + '\n', surface, s_i+1, ' from ', len(slabs), ' slabs\n', sl.natom, ' atoms'  )

                                s_i+=1

                        else:
                            self.suf_status = 'Zero slabs constructed'
                            print('\nWarning!  Zero slabs constructed!\n')

        except AttributeError:
            self.suf_status = 'Bulk calculation has some problems'
            print('\nWarning!  Bulk calculation of {} has some problems!!\n'.format(st_name)) 



    ##################### add slab calc

    if status == 'add':
        # for surface in suf_list:
            for suf_name in self.suf_status.keys():
                if self.suf_status[suf_name] in ['created']:
                    sl = self.suf[suf_name]
                    
                    if only_stoich:

                        if stoichiometry_criteria2(sl, st_bulk):
                                # print(st_name+'.sl.'+ suf_name, ise, 1, suf_name)

                                add_loop(st_name+'.sl.'+ suf_name+suffix, ise, 1 , 
                                    input_st = sl,  it_folder = 'slabs_new/'+st_name, up = 'up2', cluster = cluster, corenum = corenum)
                                self.suf_status[suf_name] = 'added'
                                self.suf_cl[suf_name] = '.'.join([st_name+'.sl.'+ suf_name, ise, '1'])

                        else:
                            print('Non-stoichiometric slab')

                    else:
                        None
                else:
                    print('\nThis slab hasn\'t been created yet!\n')
                    print(self.suf_status[suf_name])

    



    ##################### res slab calc

    if status == 'res':
        # print('1111')
        print(self.suf.keys(), self.suf_status)
        for suf_name in self.suf.keys():
            s_hkl = suf_name.split('.')[0]
            # print(s_hkl)
            key = []
            delta = 0
            for i in range(0,3):
                i+=delta
                if s_hkl[i] != '-':
                    key.append(int(s_hkl[i]))
                else:
                    key.append(-int(s_hkl[i+1]))
                    delta = 1
            # hkl = [int(s_hkl[0]), int(s_hkl[1]), int(s_hkl[2])]
            hkl = key
            if self.suf_status[suf_name] in ['added', 'calculated'] and hkl in suf_list:

                try:
                    try:
                        res_loop(st_name+'.sl.'+ suf_name+suffix, ise, 1 , up = 'up2')
                    except ValueError:
                        self.suf_status[suf_name] = 'Error'
                        self.suf_en[suf_name] = 'Error'
                        break
                except KeyError:
                    break

                self.suf_status[suf_name] = 'calculated'

                from siman.analysis import suf_en 

                try:
                    print(self.suf_cl[suf_name])
                except KeyError:
                    self.suf_cl[suf_name] = '.'.join([st_name+'.sl.'+ suf_name, ise, '1'])

                try:
                    e = '-'
                    if self.suf_cl[suf_name] in db.keys():
                        e = suf_en(db[self.suf_cl[suf_name]],db[self.bulk_cl_scale])
                    # self.suf_en = {}
                    self.suf_en[suf_name] = e
                except AttributeError:
                    self.suf_en[suf_name] = 'Error'

        if self.suf_status == 'Zero slabs constructed':
            self.suf_en = 'Zero slabs constructed'

    return


def calc_suf_stoich_mat(self, ise = '9sm', bulk_cl_name = None,  it_folder = 'bulk/', status = 'create', suf_list = None, 
        min_slab_size = 12, only_stoich = 1, corenum = 4, cluster = 'cee', conv = 1, create = 1, reset_old_polar = 0, 
        reset_old_nonpolar = 0, suffix = None, polar = 1, nonpolar = 1): 
    from siman.header import db
    from siman.geo import create_surface2, stoichiometry_criteria,stoichiometry_criteria2, symmetry_criteria, calc_k_point_mesh
    from siman.calc_manage   import add_loop, res_loop


    if not bulk_cl_name:
        bulk_cl_name = ['it','ise', '1']
        
    if not suf_list:
        suf_list = [[1,1,0], [1,1,1]]


    st_bulk = db[self.bulk_cl_scale].end.get_conventional_cell()

    # if conv:
    #     st_bulk = st_bulk.get_conventional_cell()

    st_name = '.'.join([self.pretty_formula, self.sg_crystal_str, self.sg_symbol])
    st_name = st_name.replace('/','_')
    # st_bulk.printme()



    try:
        # print(self.suf)
        self.suf_pol.keys()
    except AttributeError:
        self.suf_pol = {}
        self.suf_pol_status = {}
    try:
        # print(self.suf)
        self.suf_nonpol.keys()
    except AttributeError:
        self.suf_nonpol = {}
        self.suf_nonpol_status = {}

    try:
        self.suf_pol_cl.keys()
    except AttributeError:
        self.suf_pol_cl = {}
    try:
        self.suf_nonpol_cl.keys()
    except AttributeError:
        self.suf_nonpol_cl = {}

    try:
        self.suf_en_pol.keys()
    except AttributeError:
        self.suf_en_pol = {}
    try:
        self.suf_en_nonpol.keys()
    except AttributeError:
        self.suf_en_nonpol = {}

    if reset_old_polar:
        self.suf_pol = {}
        self.suf_pol_status = {}
        self.suf_pol_cl = {}
        self.suf_en_pol = {}

    if reset_old_nonpolar:
        self.suf_nonpol = {}
        self.suf_nonpol_status = {}
        self.suf_nonpol_cl = {}
        self.suf_en_nonpol = {}




    ##################### create slab
    if status == 'create':


        try:
            for surface in suf_list:
                print(surface)

                if self.suf_pol == {} or ''.join([str(surface[0]),str(surface[1]),str(surface[2])]) not in self.suf_pol.keys():

                    ###############polar
                    if polar:
                        print('Polar_suf')
                        slabs_polar = create_surface2(st_bulk, miller_index = surface,  min_slab_size = min_slab_size, min_vacuum_size = 15, surface_i = 0,
                          symmetrize = 0, lll_reduce = 1, primitive = 1)
                        for sl_i in slabs_polar:
                                    st = st_bulk
                                    sl = st.update_from_pymatgen(sl_i)
                                    print(' Stoichoimetry:',stoichiometry_criteria2(st,sl), ' Eqvuivalent:', sl_i.have_equivalent_surfaces())
                                    
                    ############nonpolar
                    if nonpolar:
                        print('Nonpolar_suf')
                        slabs_nonpolar = create_surface2(st_bulk, miller_index = surface,  min_slab_size = min_slab_size, min_vacuum_size = 15, surface_i = 0,
                          symmetrize = 1, lll_reduce = 1, primitive = 1)
                        for sl_i in slabs_nonpolar:
                                    st = st_bulk
                                    sl = st.update_from_pymatgen(sl_i)
                                    
                                    print(' Stoichoimetry:',stoichiometry_criteria2(st,sl), ' Eqvuivalent:', sl_i.have_equivalent_surfaces())
                    
                    if create:

                        ###############polar
                        if polar:
                            if len(slabs_polar):
                                s_i = 0
                                n = 0
                                for sl_i in slabs_polar:
                                    st = st_bulk
                                    sl = st.update_from_pymatgen(sl_i)
                                    if stoichiometry_criteria2(st,sl) and sl_i.have_equivalent_surfaces() == False:
                                        # if sl.natom < 31:
                                            suf_name = str(surface[0])+str(surface[1])+str(surface[2]) +'.'+str(s_i)

                                            # if self.suf_pol_status == 'Zero slabs constructed':
                                            #     self.suf_pol_status = {}

                                            self.suf_pol[suf_name] = sl
                                            self.suf_pol_status[suf_name] = 'created'
                                            n+=1

                                    s_i+=1
                                print(st_name + ' Polar\n', surface, n, ' from ', len(slabs_polar), ' slabs\n', sl.natom, ' atoms\n'  )

                            else:
                                suf_name = str(surface[0])+str(surface[1])+str(surface[2]) +'.'+str(0)
                                self.suf_pol_status[suf_name] = 'Zero slabs constructed'
                                print('\nWarning!  Zero slabs constructed!\n')

                        ###############nonpolar
                        if nonpolar: 
                            if len(slabs_nonpolar):
                                s_i = 0
                                n=0
                                for sl_i in slabs_nonpolar:
                                    st = st_bulk
                                    sl = st.update_from_pymatgen(sl_i)
                                    if stoichiometry_criteria2(st,sl) and sl_i.have_equivalent_surfaces() == True:
                                    

                                        # if sl.natom < 31:
                                            suf_name = str(surface[0])+str(surface[1])+str(surface[2]) +'.'+str(s_i)

                                            # if self.suf_nonpol_status == 'Zero slabs constructed':
                                            #     self.suf_nonpol_status = {}

                                            self.suf_nonpol[suf_name] = sl
                                            self.suf_nonpol_status[suf_name] = 'created'

                                            n+=1
                                    s_i+=1
                                print(st_name + ' Nonolar\n', surface, n, ' from ', len(slabs_nonpolar), ' slabs\n', sl.natom, ' atoms\n'  )

                            else:
                                suf_name = str(surface[0])+str(surface[1])+str(surface[2]) +'.'+str(0)

                                self.suf_nonpol_status[suf_name] = 'Zero slabs constructed'
                                print('\nWarning!  Zero slabs constructed!\n')

        except AttributeError:
            self.suf_pol_status = 'Bulk calculation has some problems'
            self.suf_nonpol_status = 'Bulk calculation has some problems'
            print('\nWarning!  Bulk calculation of {} has some problems!!\n'.format(st_name)) 



    ##################### add slab calc

    if status == 'add':
        for suf_name in self.suf_pol_status.keys():
            if self.suf_pol_status[suf_name] != 'Zero slabs constructed':
                suffix = '.polar'

                if self.suf_pol_status[suf_name] in ['created']:
                    sl = self.suf_pol[suf_name]
                    
                    if only_stoich:

                        if stoichiometry_criteria2(sl, st_bulk):
                                print(sl.rprimd)
                                ngkpt=calc_k_point_mesh(sl.rprimd, kspacing = 0.2)
                                ngkpt=ngkpt[:2]+(1,)
                                add_loop(st_name+'.sl.'+ suf_name+suffix, ise, 1 , ngkpt = ngkpt,
                                    input_st = sl,  it_folder = 'slabs_new/'+st_name, up = 'up2', cluster = cluster, corenum = corenum)
                                self.suf_pol_status[suf_name] = 'added'
                                self.suf_pol_cl[suf_name] = '.'.join([st_name+'.sl.'+ suf_name+suffix, ise, '1'])

                        else:
                            print('Non-stoichiometric slab')

                    else:
                        None
                else:
                    print('\nThis slab hasn\'t been created yet!\n')
                    print(self.suf_pol_status[suf_name])

        
            print(self.suf_nonpol_status)
        for suf_name in self.suf_nonpol_status.keys():
            if self.suf_nonpol_status[suf_name] != 'Zero slabs constructed':
                suffix = '.nonpolar'

                if self.suf_nonpol_status[suf_name] in ['created']:
                    sl = self.suf_nonpol[suf_name]
                    
                    if only_stoich:

                        if stoichiometry_criteria2(sl, st_bulk):
                                # print(st_name+'.sl.'+ suf_name, ise, 1, suf_name)
                                ngkpt=calc_k_point_mesh(sl.rprimd, kspacing = 0.2)
                                ngkpt=ngkpt[:2]+(1,)
                                add_loop(st_name+'.sl.'+ suf_name+suffix, ise, 1 , ngkpt = ngkpt,
                                    input_st = sl,  it_folder = 'slabs_new/'+st_name, up = 'up2', cluster = cluster, corenum = corenum)
                                self.suf_nonpol_status[suf_name] = 'added'
                                self.suf_nonpol_cl[suf_name] = '.'.join([st_name+'.sl.'+ suf_name+suffix, ise, '1'])

                        else:
                            print('Non-stoichiometric slab')

                    else:
                        None
                else:
                    print('\nThis slab hasn\'t been created yet!\n')
                    print(self.suf_nonpol_status[suf_name])


    ##################### res slab calc

    if status == 'res':
    #     # print('1111')
        print(self.suf_pol.keys(), self.suf_pol_status)
        suffix = '.polar'
        for suf_name in self.suf_pol.keys():
            s_hkl = suf_name.split('.')[0]
    #         # print(s_hkl)
            key = []
            delta = 0
            for i in range(0,3):
                i+=delta
                if s_hkl[i] != '-':
                    key.append(int(s_hkl[i]))
                else:
                    key.append(-int(s_hkl[i+1]))
                    delta = 1
            # hkl = [int(s_hkl[0]), int(s_hkl[1]), int(s_hkl[2])]
            hkl = key
            if self.suf_pol_status[suf_name] in ['added', 'calculated'] and hkl in suf_list:

                try:
                    try:

                        res_loop(st_name+'.sl.'+ suf_name+suffix, ise, 1 , up = 'up2')
                    except ValueError:
                        self.suf_pol_status[suf_name] = 'Error'
                        self.suf_en_pol[suf_name] = 'Error'
                        break
                except KeyError:
                    break

                self.suf_pol_status[suf_name] = 'calculated'

                from siman.analysis import suf_en 

                try:
                    print(self.suf_pol_cl[suf_name])
                except KeyError:
                    self.suf_pol_cl[suf_name] = '.'.join([st_name+'.sl.'+ suf_name, ise, '1'])

                try:
                    e = '-'
                    if db[self.suf_pol_cl[suf_name]]:# in db.keys():
                        e = suf_en(db[self.suf_pol_cl[suf_name]],db[self.bulk_cl_scale])
                    # self.suf_en = {}
                    try:
                        self.suf_en_pol[suf_name] = round(float(e),3)

                    except ValueError:
                        self.suf_en_pol[suf_name] = e
                except AttributeError:
                    self.suf_en_pol[suf_name] = 'Error'

        if self.suf_pol_status == 'Zero slabs constructed':
            self.suf_en_pol = 'Zero slabs constructed'
                    

        print(self.suf_nonpol.keys(), self.suf_nonpol_status)
        suffix = '.nonpolar'
        for suf_name in self.suf_nonpol.keys():
            s_hkl = suf_name.split('.')[0]
    #         # print(s_hkl)
            key = []
            delta = 0
            for i in range(0,3):
                i+=delta
                if s_hkl[i] != '-':
                    key.append(int(s_hkl[i]))
                else:
                    key.append(-int(s_hkl[i+1]))
                    delta = 1
            # hkl = [int(s_hkl[0]), int(s_hkl[1]), int(s_hkl[2])]
            hkl = key
            if self.suf_nonpol_status[suf_name] in ['added', 'calculated'] and hkl in suf_list:

                try:
                    try:

                        res_loop(st_name+'.sl.'+ suf_name+suffix, ise, 1 , up = 'up2')
                    except ValueError:
                        self.suf_nonpol_status[suf_name] = 'Error'
                        self.suf_en_nonpol[suf_name] = 'Error'
                        break
                except KeyError:
                    break

                self.suf_nonpol_status[suf_name] = 'calculated'

                from siman.analysis import suf_en 

                try:
                    print(self.suf_nonpol_cl[suf_name])
                except KeyError:
                    self.suf_nonpol_cl[suf_name] = '.'.join([st_name+'.sl.'+ suf_name, ise, '1'])

                try:
                    e = '-'
                    if db[self.suf_nonpol_cl[suf_name]]:# in db.keys():
                        e = suf_en(db[self.suf_nonpol_cl[suf_name]],db[self.bulk_cl_scale])
                    # self.suf_en = {}
                    try:
                        self.suf_en_nonpol[suf_name] = round(float(e),3)

                    except ValueError:
                        self.suf_en_nonpol[suf_name] = e


                except AttributeError:
                    self.suf_en_nonpol[suf_name] = 'Error'

        if self.suf_nonpol_status == 'Zero slabs constructed':
            self.suf_en_nonpol = 'Zero slabs constructed'

def add_relax_mat(self, ise, status = 0):
    from siman.header import db
    from siman.analysis import suf_en
    from siman.geo import  calc_k_point_mesh

    if status == 'add':
        try:
            self.suf_pol_cl_new.keys()
        except AttributeError:
            self.suf_pol_cl_new = {}
        try:
            self.suf_nonpol_cl_new.keys()
        except AttributeError:
            self.suf_nonpol_cl_new = {}

        for suf in self.suf_pol_cl.keys():
            # for cl_name in self.suf_pol_cl[suf]:
                cl_name = self.suf_pol_cl[suf]
                try:
                    db[cl_name].e0
                    ngkpt=calc_k_point_mesh(db[cl_name].end.rprimd, kspacing = 0.2)
                    ngkpt=ngkpt[:2]+(1,)
                    db[cl_name].run(ise, ngkpt = ngkpt)
                    name = cl_name.split('.')[:-2]
                    name.extend(['ifn',ise, '1'])
                    name = '.'.join(name)
                    # print(name)
                    self.suf_pol_cl_new[suf] = name
                except AttributeError:
                    None

        for suf in self.suf_nonpol_cl.keys():
            # for cl_name in self.suf_nonpol_cl[suf]:
                cl_name = self.suf_nonpol_cl[suf]
                try:
                    db[cl_name].e0
                    ngkpt=calc_k_point_mesh(db[cl_name].end.rprimd, kspacing = 0.2)
                    ngkpt=ngkpt[:2]+(1,)
                    db[cl_name].run(ise, ngkpt = ngkpt)
                    name = cl_name.split('.')[:-2]
                    name.extend(['ifn',ise, '1'])
                    name = '.'.join(name)
                    # print(name)
                    self.suf_nonpol_cl_new[suf] = name
                except AttributeError:
                    None

    if status == 'res':
        bulk_new_name = self.bulk_cl_scale.split('.')[:-2]
        bulk_new_name.extend(['ifn','8newm', '100'])
        bulk_new_name = '.'.join(bulk_new_name)
        self.bulk_cl_scale_new = bulk_new_name

        for suf in self.suf_pol_cl_new.keys():
            try:
                db[self.suf_pol_cl_new[suf]].res()
                suf_en(db[self.suf_pol_cl_new[suf]], db[self.bulk_cl_scale_new])
            except AttributeError:
                None
        for suf in self.suf_nonpol_cl_new.keys():
            try:
                db[self.suf_nonpol_cl_new[suf]].res()
                suf_en(db[self.suf_nonpol_cl_new[suf]], db[self.bulk_cl_scale_new])
            except AttributeError:
                None
def move_suf_en_mat(self, only_show = 1):
    from siman.header import db
    from siman.geo import stoichiometry_criteria2, symmetry_criteria
    print('\n\n',self.pretty_formula)
    # print(self.suf_en)
    st_bulk = db[self.bulk_cl_scale].end
    if only_show == 0:
        self.suf_nonpol_cl = {}
        self.suf_en_nonpol = {}
        self.suf_pol_cl = {}
        self.suf_en_pol = {}

    if self.suf_en != 'Zero slabs constructed':
        for suf_name in self.suf_en.keys():
            try:
                if self.suf_en[suf_name] != 'Error':
                    print('\n\n\n\n',suf_name, round(self.suf_en[suf_name], 2), 'J/m2')
                else:
                    print('\n\n\n\n',suf_name, self.suf_en[suf_name])

                sl = db[self.suf_cl[suf_name]].init
                if stoichiometry_criteria2(sl, st_bulk, silent = 0):
                    if symmetry_criteria(sl):
                        if only_show == 0:
                            self.suf_nonpol_cl[suf_name] = self.suf_cl[suf_name]
                            self.suf_en_nonpol[suf_name] = self.suf_en[suf_name]
                    else:
                        # print('Non-symmetric')
                        if only_show == 0:
                            self.suf_pol_cl[suf_name] = self.suf_cl[suf_name]
                            self.suf_en_pol[suf_name] = self.suf_en[suf_name]

            except KeyError or AttributeError:
                None



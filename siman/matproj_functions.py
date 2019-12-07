import os, sys

from siman.header import db
from siman.geo import  create_surface2
from siman.calc_manage import  get_structure_from_matproj, smart_structure_read
from siman.calc_manage   import (clean_history_file, prepare_run,  manually_remove_from_struct_des, update_des, inherit_icalc, add_loop, res_loop, complete_run, add_des )

try:
    # pmg config --add VASP_PSP_DIR $VASP_PSP_DIR MAPI_KEY $MAPI_KEY
    from pymatgen.ext.matproj import MPRester
    from pymatgen.io.vasp.inputs import Poscar
    from pymatgen.io.cif import CifParser
    pymatgen_flag = True 
except:
    print('pymatgen is not available')
    pymatgen_flag = False 

from siman.classes import CalculationVasp


import csv




def get_pmg_info2(criteria, properties, name = 'table', price = 0, element_price = None):

    #get data from  materials project server


    # criteria - string with condition for query method to choice some structures from MatProj
    # properties - list of MatProj structure keys in  m.get_data('mp-12957') to write
    # price - logical 
    # element_price - dict {element: price per kg}

   
    data = {}


    with MPRester(header.pmgkey) as m:
            properties.append('elements')

            results = m.query(criteria= criteria, properties=properties)

            for string_to_write in results:

                if price:

                    if 'price_per_gramm' not in properties:
                        properties.append('price_per_gramm')
                    cost = 0
                    for p in element_price.keys():
                        if p in string_to_write['elements']:
                            # print(element_price[p])
                            cost+= element_price[p]/1000


                    string_to_write.update([('price_per_gramm', cost)])



                st_name = string_to_write['pretty_formula']


                if st_name not in data.keys():         
                    data.update([(st_name, string_to_write)])
                else:
                    if string_to_write['e_above_hull'] < data[st_name]['e_above_hull']:
                        data.update([(st_name, string_to_write)])
                        

    # write data in csv file
    with open(name+'.csv', 'w', newline='') as csvfile:
        fieldnames = properties
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for st_name in sorted(data.keys()):
            # print(data[st_name])
            writer.writerow(data[st_name])

    return data


###############################################3
### read and write pmg csv file
#############################################
def read_pmg_info(path):

    data = []
    with open(path+'.csv', newline='') as csvfile:
        reader = csv.DictReader(csvfile)

        for d in reader:
            data.append(d)

    return data

def write_pmg_info(data, path):

    with open(path+'.csv', 'w', newline='') as csvfile:
        fieldnames = data[0].keys()
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i in data:
            # print(data[st_name])
            writer.writerow(i)

###################################################################################################
######################################################################################################3


def get_mat_st(mat_in_list, folder = 'geo/'):
    #check downloaded POSCAR files in geo/ folder
    #if not POSCAR of some structure - download it from Mat Proj

    # mat_in_list - data dict for any structure from MP,  result of get_fata('mp-...')

    
    name = mat_in_list['material_id']+'.POSCAR'
    if name not in os.listdir(folder):
        os.chdir(folder)
        st = get_structure_from_matproj(mat_proj_id = mat_in_list['material_id'], it_folder = folder)
        os.chdir('..')
    else:
        st = smart_structure_read(folder+name)
        # print('ok')
    return st




def calc_bulk_list(data_list, spacegroup = '', ise = '8', ise_mag = '8m', status = None, corenum = 1):

    # function can add or res for set of bulk calculations
    
    # data_list = read_pmg_info() of some csv file

    # ise  - set of calculation. Usually use '8' for nonmag calc and '8m' for mag calc
    #status = add or res

    spacegroup = '.'+spacegroup
    for i in data_list:
        st = get_mat_st(i)


        if float(i['total_magnetization']) > 1e-3: 
            mag_flag = 1
            ise = ise_mag
        else: 
            mag_flag = 0
            ise = ise


        # if mag_flag:
        if status == 'add': 
            add_loop(i['pretty_formula']+spacegroup, ise, 1, it_folder = 'bulk', input_st = st, corenum = corenum)
        if status == 'res': 
            res_loop(i['pretty_formula']+spacegroup, ise, 1, it_folder = 'bulk')





def stoichiometry_criteria(st1,st2):

    natom1 = st1.get_natom()
    natom2 = st2.get_natom()

    tra1 = st1.get_transition_elements()
    tra2 = st2.get_transition_elements()
    ntra1 = len(tra1)
    if ntra1 == 0: 
        ntra1 = natom1
    ntra2 = len(tra2)
    if ntra2 == 0: 
        ntra2 = natom2
    rat1 = natom1/ntra1
    rat2 = natom2/ntra2
    mul = ntra1/ntra2

    if rat1 == rat2:
        return 1
    else:
        return 0





def calc_suf_list_sg(data_list, sg, suf_list, flag = 0):
    # print('start')
    for i in data_list:

        
        if i['spacegroup.symbol'] == sg:
            print(i['spacegroup.symbol'])
            st_name = i['pretty_formula']
            print(st_name)

            if float(i['total_magnetization']) > 1e-3: 
                mag_flag = 1
                ise = '9sm'
                st_bulk = db[st_name, '8m', 1].end
            else: 
                mag_flag = 0
                ise = '9s'
                st_bulk = db[st_name, '8', 1].end
            

            for surface in suf_list:
                slabs = create_surface2(st_bulk, miller_index = surface,  min_slab_size = 10, min_vacuum_size = 10, surface_i = 0,
                  symmetrize = True)
                
                s_i = 0
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


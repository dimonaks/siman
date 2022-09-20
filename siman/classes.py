# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
from __future__ import division, unicode_literals, absolute_import, print_function
import itertools, os, copy, math, glob, re, shutil, sys, pickle, gzip, shutil, random
import re, io, json
import pprint

from textwrap import wrap

import numpy as np

#additional packages
try:
    from tabulate import tabulate
except:
    print('tabulate is not avail')
    tabulate = None
try:
    import pandas as pd
except:
    print('pandas is not avail')


from siman import header

try:
    import pymatgen
    header.pymatgen_flag = True
except:
    print('classes.py: pymatgen is not available')
    header.pymatgen_flag = False

if header.pymatgen_flag:
    from pymatgen.io.cif import CifWriter
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.core.surface import Slab
    from pymatgen.core.composition import Composition


# import matplotlib.pyplot as plt

#siman packages

from siman.header import printlog, print_and_log, runBash, plt

from siman import set_functions
# from siman.small_functions import return_xred, makedir, angle, is_string_like, cat_files, grep_file, red_prec, list2string, is_list_like, b2s, calc_ngkpt, setting_sshpass
from siman.small_functions import return_xred, makedir, angle, is_string_like, cat_files, grep_file, red_prec, list2string, is_list_like, b2s, calc_ngkpt, setting_sshpass
from siman.functions import (read_vectors, read_list, words, read_string,
     element_name_inv, invert, calculate_voronoi, 
    get_from_server, push_to_server, run_on_server, smoother, file_exists_on_server, check_output)
from siman.inout import write_xyz, write_lammps, read_xyz, read_poscar, write_geometry_aims, read_aims_out, read_vasp_out
from siman.geo import (image_distance, replic, calc_recip_vectors, calc_kspacings, xred2xcart, xcart2xred, 
local_surrounding, local_surrounding2, determine_symmetry_positions, remove_closest, remove_vacuum, make_neutral, 
rms_between_structures, rms_between_structures2)
from siman.set_functions import InputSet, aims_keys


"""
Classes used in siman
TODO:
1. Please  combine calculate_nbands(), calc_kspacings(), magmom filling in make incar  with actualize_set() 
2. split make_incar_and_copy_all() into make_incar() and copy_calc_files_to_cluster()
3. split .read_results() into download_output_files() and .parse_outcar() and .analyze_output()
4. outcar name should be returned by write_sge_script in all cases and used, now only in u_ramping
5. write_sge() - в режиме inherit_option = continue - предыдущие outcar резервируются только один
раз с префиксом prev, повторный запуск перезатрет их, поэтому нужно писать спец код
типа
if test -f prev1.outcar
    cp name prev2+name
чтобы она находила prev3 с максимальным числом, и к этому числу прибавляла единицу для нового файла

NEW:
Calculation Structure():
    *self.magmom* (list) - magnetic moments for each ion in structure; has higher preference than self.set.magnetic_moments which
    include only moments for atom types


"""


'aliases required for old databases to work correctly:'
from siman.core.structure import Structure
from siman.core.calculation import Calculation
from siman.calculators.vasp import CalculationVasp

# class Calculation()


class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


class empty_struct(): # here it is needed for back-compatability for reading old databases
    def __init__(self):
        pass



class Description():
    """
    Objects of this class include just folder and description of specific calculation.
    Mostly was needed for manual addition of new calculations

    self.ngkpt_dict_for_kspacings (dict of lists) - the key is kspacing; the dict
    contains k-meshes 
    for all calculations
    based on this geometry structure.
    can be useful for fine tuning of k-mesh for specific kspacing.

    """
    def __init__(self, sectionfolder = "forgot_folder", description = "forgot_description"):
        self.des = description
        self.sfolder = sectionfolder
        self.ngkpt_dict_for_kspacings = {} #the key is kspacing







class MP_Compound():
    """This class includes information about chemical compounds from MatProj and next operations (bulk calc, slab construction etc.)

    db key is 'pretty_formula.MP': ('AgC.MP')
    """
    def __init__(self):
        self.pretty_formula = ""
        self.material_id = "material_id"
        self.elements = []
        self.sg_symbol =''
        self.sg_crystal_str = ''
        self.band_gap = None
        self.e_above_hull = None
        self.icsd_ids = None
        self.total_magnetization = None
        self.price_per_gramm = None


        self.bulk_cl = None
        self.bulk_status = 'Unknown'


    def copy(self):
        return copy.deepcopy(self)






    def calc_bulk(self, ise, bulk_cl_name = ['it','ise', '1'], it_folder = 'bulk/', status = 'add'):
        from siman.header import db
        from siman.calc_manage   import add_loop, res_loop

        name = '.'.join(bulk_cl_name)

        it = '.'.join([self.pretty_formula, self.sg_crystal_str])
        st = self.get_st()
        self.bulk_cl = name


        if status == 'add':
            # if bulk_cl_name[0] not in header.struct_des:

                add_loop(it,ise,1, input_st = st, it_folder = it_folder, override = 1)
                self.bulk_status = 'run'

        if status == 'res':
            res_loop(it,ise,1)
            try: 
                if max(db[name].maxforce_list[-1]) > 50:
                    self.bulk_status = 'big max_f'
                else:
                    self.bulk_status = 'calculated'

            except AttributeError:
                    self.bulk_status = 'Error!'
                    print(name, '\tUnfinished calculation!!!\n\n')

        if status == 'add_scale':
            # if bulk_cl_name[0] not in header.struct_des:

                add_loop(it,ise,1, input_st = st, calc_method = 'uniform_scale',  scale_region = (-9, 5), n_scale_images = 10, it_folder = it_folder)
                self.bulk_status_scale = 'run_scale'
        
        if status == 'res_scale':
            # if bulk_cl_name[0] not in header.struct_des:

                name_scale = '.'.join([it,'su',ise,'100'])
                self.bulk_cl_scale = name_scale

                try:
                    # res_loop(it+'.su',ise,list(range(1,11))+[100], up = 'up2', show = 'fit', analys_type = 'fit_a')
                    res_loop(it+'.su',ise,[100], up = 'up2')
                except ValueError:
                    self.bulk_status_scale = 'Error!'
                    print('\n\nValueError!!!\n\n')
                    return
                try: 
                    if max(db[name_scale].maxforce_list[-1]) > 50:
                        self.bulk_status_scale = 'big max_f'
                    else:
                        self.bulk_status_scale = 'calculated'

                except AttributeError:
                    self.bulk_status_scale = 'Error!'
                    print(name_scale, '\tUnfinished calculation!!!\n\n')


    def calc_suf(self, **argv):
        from siman.matproj_functions import calc_suf_mat
        calc_suf_mat(self, **argv)

    def calc_suf_stoich(self, **argv):
        from siman.matproj_functions import calc_suf_stoich_mat
        calc_suf_stoich_mat(self, **argv)

    def add_relax(self, **argv):
        from siman.matproj_functions import add_relax_mat
        add_relax_mat(self, **argv)
    def move_suf_en(self, **argv):
        from siman.matproj_functions import move_suf_en_mat
        move_suf_en_mat(self, **argv)


    def get_st(self, folder = 'geo/'):
        """
        check downloaded POSCAR files in geo/ folder
        if not POSCAR of some structure - download it from Mat Proj

        mat_in_list - data dict for any structure from MP,  result of get_fata('mp-...')
        """

        from siman.calc_manage import  get_structure_from_matproj, smart_structure_read
        
        
        name = self.material_id+'.POSCAR'
        # st = get_structure_from_matproj(mat_proj_id = self.material_id, it_folder = folder)

        if name not in os.listdir(folder):
            os.chdir(folder)
            st = get_structure_from_matproj(mat_proj_id = self.material_id, it_folder = folder)
            os.chdir('..')
        else:
            st = smart_structure_read(folder+name)
            # print('ok')
        return st

    def e_cohesive_calc(self, e_box):
        from siman.header import db
        '''
        return cohesive energy

        e_box - dict{element: energy_of_element_in_box}
        '''
        
        # print(self.bulk_cl_scale)
        try:
            cl_bulk = db[self.bulk_cl_scale]
            e_bulk = cl_bulk.energy_sigma0
            n_at_sum = cl_bulk.end.natom
            el_list = cl_bulk.end.get_elements()

            e_at_sum = 0
            for el in el_list:
                e_at = e_box[el]
                e_at_sum+=e_at

            e_coh = (e_at_sum-e_bulk)/n_at_sum
            print('{}  \t\tE_coh = {} eV'.format(self.pretty_formula, round(e_coh,1)))
            self.e_cohesive = round(e_coh,2)
            
        except AttributeError:
            self.e_cohesive = None


    def e_cohesive_from_MP(self):


        try:
            from pymatgen.ext.matproj import MPRester
            from pymatgen.io.vasp.inputs import Poscar
            from pymatgen.io.cif import CifParser
            pymatgen_flag = True 
        except:
            print('pymatgen is not available')
            pymatgen_flag = False 

        with MPRester(header.pmgkey) as m:
            material_id = self.material_id
            ec = round(m.get_cohesive_energy(material_id, per_atom = 1),2)
            self.e_cohesive_MP = ec
            print(self.pretty_formula, ec)


    def calc_ec_es(self, ev = 0):
        from siman.header import db
        '''
        
        '''

        ec_es = []

        try:
            print(self.e_cohesive)
            for i in self.suf_en.keys():
                # print(self.suf_en)
                if self.suf_en[i] !='Error':
                    # print(self.suf_en[i])
                    if ev:
                        ec_es.append(round(float(self.e_cohesive)/float(self.suf_en[i]/ header.eV_A_to_J_m), 2))
                    else:
                        ec_es.append(round(float(self.e_cohesive)/float(self.suf_en[i]), 2))
                else:
                    ec_es.append('None')

            self.ec_es = ec_es
        except AttributeError:
            self.ec_es = None
        

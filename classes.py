# -*- coding: utf-8 -*- 
#Copyright Aksyonov D.A
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

# import pymatgen
# sys.exit()
from siman import header

try:
    import pymatgen
    from pymatgen.io.cif import CifWriter
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.core.surface import Slab
    from pymatgen.core.composition import Composition
    header.pymatgen_flag = True
except:
    print('pymatgen is not avail')
    header.pymatgen_flag = False

# import matplotlib.pyplot as plt

#siman packages

from siman.header import printlog, print_and_log, runBash, plt

from siman import set_functions
from siman.small_functions import return_xred, makedir, angle, is_string_like, cat_files, grep_file, red_prec, list2string, is_list_like, b2s, calc_ngkpt, setting_sshpass
from siman.functions import (read_vectors, read_list, words, read_string,
     element_name_inv, invert, calculate_voronoi, update_incar, 
    get_from_server, push_to_server, run_on_server, smoother, file_exists_on_server, check_output)
from siman.inout import write_xyz, write_lammps, read_xyz, read_poscar, write_geometry_aims, read_aims_out, read_vasp_out
from siman.geo import (image_distance, replic, calc_recip_vectors, calc_kspacings, xred2xcart, xcart2xred, 
local_surrounding, local_surrounding2, determine_symmetry_positions, )
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






class Structure():
    """This class includes only structure related information such as primitive vectors, coordinates, forces and so on"""
    def __init__(self):
        self.name = ""
        self.des = ''
        self.hex_a = None
        self.hex_c = None
        self.gbpos = None
        self.rprimd = [np.zeros((3)) for i in [0,1,2] ]
        self.xcart = []
        self.xred = []
        self.magmom = []
        self.select = [] # flags for selective dynamics
        self.vel    = [] # velocities, could be empty
        self.predictor = None # string , predictor-corrector information to continiue MD run
        # self.pos = write_poscar()

    def copy(self):
        return copy.deepcopy(self)

    def new(self):
        return Structure()

    def selective_all(self):
        """
        allow all atoms to relax
        """
        st = copy.deepcopy(self)
        # if hasattr(self, 'select'):
        st.select = []
        # print(st.natom)
        for i in range(st.natom):
            st.select.append([True,True,True])
            # print('9')
        # print(st.select)
        return st 

    def check_selective(self):
        """
        check if some atoms should be frozen
        """
        selective_dyn = False
        if self.select is not None:
            for sel in self.select:
                # print (not all(sel))
                # print(sel)
                if not all(sel):
                    selective_dyn = True
        return selective_dyn

    def selective_byCompare(self, st2, tol = 0.0001, freeze = 'present' ):
        """
        set selective dynamics falgs in the calling structure (self) by freezing the atoms that are present (or missing) in the supplied structure (st2)
        tol - tolerance for finding the same atoms in the two structures
        freeze = 'present' or 'missing'

        TODO:
        read st2 from file
        optional save flag-changed atoms to cif (to check the result)

        """
        st1 = copy.deepcopy(self)
        if freeze == 'missing':
            flag_change = [True, True, True]
            flag_default = [False, False, False]
        else:
            flag_change = [False, False, False]
            flag_default = [True, True, True]
            if freeze != 'present':
                print('Warning! incorrect \'freeze\' argument, \'presnt\' is used')

        tol = tol ** 2  # ^2 tol instead sqrt(dist)

        totNatom1 = len(st1.typat)
        totNatom2 = len(st1.typat)
        if totNatom2 > totNatom1:
            printlog('Warning! struct to compare with has more atoms than original struct!')

        st1.select = [flag_default] * totNatom1

        natom1 = {} # dictionary of atom types (by znucl); each entry contains list of atom nubers of certain type
        for i, atype in enumerate(st1.typat):
            if st1.znucl[atype-1] in natom1:
                natom1[st1.znucl[atype-1]].append(i)
            else:
                natom1[st1.znucl[atype-1]] = [i]

        natom2 = {}  # dictionary of atom types (by znucl); each entry contains list of atom nubers of certain type
        for i, atype in enumerate(st2.typat):
            if st2.znucl[atype-1] in natom2:
                natom2[st2.znucl[atype-1]].append(i)
            else:
                natom2[st2.znucl[atype-1]] = [i]

        for ztype in natom2:
            for i2 in natom2[ztype]:
                for k, i1 in enumerate(natom1[ztype]):
                    dist = np.sum(np.square(st1.xred[i1] - st2.xred[i2])) # square distance between compared atoms
                    if dist < tol:
                        st1.select[i1] = flag_change # change SD flags for if the atoms are identical
                        del (natom1[ztype][k])
        return st1


    def fix_atoms(self, numbers = None):
        #fix all atom with numbers in numbers list
        
        #TODO allow choice of axis of fix 
        st = copy.deepcopy(self)
        
        if len(st.select) != st.natom:
            printlog('Warning! number of select atom is not equal to total number of atoms. First I make all moveable')
            st = st.selective_all()

        # print(st.select)
        for i in numbers:
            st.select[i] = [False, False, False]

        return st

    def fix_layers(self, xred_range = None, xcart_range = None, highlight = False):
        """
        fix atoms in layers normal to R3

        xred_range (list) [from, to]
        highlight - replace with Pu to check

        """
        st = copy.deepcopy(self)

        r3 = np.linalg.norm(st.rprimd[2])
        xcr = xcart_range
        if xcr:
            xred_range = [xcr[0]/r3, xcr[1]/r3]
            printlog('xcart_range converted to xred', xred_range)



        if not hasattr(st, 'select') or len(st.select)==0:
            st = st.selective_all()
        fixed = []
        
        # print(st.select)

        for i, xr in enumerate(st.xred):
            if xred_range[0]  < xr[2] < xred_range[1]:
                # print(xred_range[0],xr[2],xred_range[1])
                st.select[i] = [False, False, False] # fix
                fixed.append(i)
        
        st.name+='_fix'
        if highlight:
            st_h = st.replace_atoms(fixed, 'Pu')
            st_h.write_poscar()

        return st





    def get_layers_pos(self, xred_range):
        #return layer positions along vector 3 in xred_range
        st = self
        zred_req = []

        for i, xr in enumerate(st.xred):
            z =xr[2]
            if xred_range[0]  < z < xred_range[1]:
                if len(zred_req) == 0:
                    zred_req.append(z)

                m = min(np.abs(np.array(zred_req) - z))
                # print()
                if m > 0.05/np.linalg.norm(st.rprimd[2]): #tolerance 0.1 A
                    zred_req.append(xr[2])
        z_unique = sorted(zred_req)
        return z_unique


    def get_slice(self, zr_range):
        #return element numbers from the zr_range
        #slab should start from bottom
        st = self
        z2 = st.get_surface_pos()[1]
        nn = []
        for i, xr in enumerate(st.xred):
            if zr_range[0] < xr[2] <= zr_range[1]:
                nn.append(i)

        return nn



    def xcart2xred(self,):
        self.xred = xcart2xred(self.xcart, self.rprimd)
        self.natom = len(self.xred)
    def update_xred(self,):
        self.xred = xcart2xred(self.xcart, self.rprimd)
        self.natom = len(self.xred)


    def xred2xcart(self,):
        self.xcart = xred2xcart(self.xred, self.rprimd)

    def update_xcart(self,):
        self.xcart = xred2xcart(self.xred, self.rprimd)

    def exchange_axes(self, i1_r, i2_r):
        """
        
        """
        st = copy.deepcopy(self)
        r = copy.deepcopy(st.rprimd)

        st.rprimd[i1_r] = r[i2_r]
        st.rprimd[i2_r] = r[i1_r]

        st.update_xred()

        return st




    def get_volume(self):
        self.vol = np.dot( self.rprimd[0], np.cross(self.rprimd[1], self.rprimd[2])  ); #volume
        return self.vol

    def get_recip(self):
        """Calculate reciprocal vectors"""
        self.recip = calc_recip_vectors(self.rprimd)
        return self.recip

    def get_nznucl(self):
        """list of numbers of atoms of each type, order is as in typat and znucl
            updated directly
        """
        self.nznucl = []
        for typ in range(1,self.ntypat+1):
            self.nznucl.append(  self.typat.count(typ) )
        return self.nznucl

    def get_elements(self):
        #return list of elements names
        # print(self.typat)
        return [element_name_inv(self.znucl[t-1]) for t in self.typat]
    
    def get_el_name(self, i):
        #return name of element
        return self.get_elements()[i]

    def get_elements_z(self):
        #return list of elements names
        return [self.znucl[t-1] for t in self.typat]
    def get_el_z(self, i):
        #return name of element
        return self.get_elements_z()[i]
    def get_natom(self):
        #get number of real atoms, excluding voids
        # print([z for z in self.get_elements_z() if z != 300])
        return len([z for z in self.get_elements_z() if z != 300])

    def get_elements_zval(self):
        #return list with number of valence electrons for each atom
        zvals = []
        for z in self.get_elements_z():
            i = self.znucl.index(z)
            zv = self.zval[i]
            zvals.append(zv)
        return zvals




    def el_diff(self, st2, mul = 1, silent = 0):
        """
        Determine difference in number of atoms between two structures self and st2
        mul (int) - allows to compare supercells; self.natom = mul * st2.natom
        
        RETURN:
        dict[key], where key is element and the value is difference in number of atoms of this element
        """
        st1 = self

        els1 = st1.get_elements()
        els2 = st2.get_elements()
        uniqe_elements = list(set(els1+els2))
        el_dif = {} # difference of elements between slab and normalized by transition metals bulk phase
        for el in uniqe_elements:
            dif = els1.count(el) - mul * els2.count(el)
            if not float(dif).is_integer():
                printlog('Error! difference of atom numbers is not integer for element ', el, 'something is wrong')
            if abs(dif) > 0:
                el_dif[el] = int(dif) / mul

        if not silent:
            print('The following elements are off-stoicheometry in the slab', el_dif, 'please provide corresponding chemical potentials')
        return el_dif






    def get_total_number_electrons(self):
        zvals = self.get_elements_zval()
        return int(sum(zvals))


    def determine_symmetry_positions(self, element):
        from siman.geo import determine_symmetry_positions

        return determine_symmetry_positions(self, element)


    def get_maglist(self):
        """
        return bool list of which  elements are magnetic (here all transition metals are searched!)
        and dictionary with numbers in lists for each transition metal
        
        RETURN:
            ifmaglist (list of bool) - magnetic or not
            mag_numbers (dict of int) - for each element using z, their numbers
        """

        ifmaglist = []
        zlist = self.get_elements_z()
        mag_numbers = {}
        for i, z in enumerate(zlist): #
            if z in header.TRANSITION_ELEMENTS:
                if z not in mag_numbers:
                    mag_numbers[z] = []
                ifmaglist.append(True)
                mag_numbers[z].append(i)
            else:
                ifmaglist.append(False)
        ifmaglist = np.array(ifmaglist)

        return ifmaglist, mag_numbers

    def show_mag(self, i):
        # show magmom of i atoms, starting from 1
        i-=1
        m = self.magmom[i]
        el = self.get_elements()[i]
        print('Mag for {:s} = {:.2f}'.format(el, m))

        return 

    def group_magmom(self, tol = 0.3):
        """
        Group magmom according to values 
        tol - tolerance according to which the magmom are grouped

        """
        st = self
        lengths = list(np.around(st.magmom, 2))
        unique = []
        unique.append(lengths[0])
        groups_size = {}
        groups_nums = {}
        for l in lengths[1:]:
            if min(np.abs(unique-l)) > tol:
                unique.append(l)
        # print('lengths', lengths)
        # print('unique bonds are', unique)
        for u in unique:
            groups_size[u] = 0
            groups_nums[u] = []
            for i, l in enumerate(lengths):
                if abs(l-u) < tol:
                    groups_size[u] += 1
                    groups_nums[u].append(i)
        return groups_size, groups_nums


    def get_mag_tran(self, to_ox = None, silent = 0):
        #show formatted mag moments of transition metals 
        #to_ox - convert to oxidation state, substract from to_ox
        # if to_ox is negative, then m-to_ox
        l, mag_numbers = self.get_maglist()

        keys = list(mag_numbers.keys())#[0]
        print('The following TM are found:', keys)

        mag = list(np.array(self.magmom)[l])
        magnetic = None
        magnetic_all = []
        for key in keys:
        
            # print(mag)
            magnetic = mag[:len(mag_numbers[key])]
            mag = mag[len(mag_numbers[key]):]
            # print(magnetic)

            s = ' '.join(['{:5.2f} ']*len(magnetic))
            
            s0 = ' '.join(['{:5d} ']*len(magnetic))

            # print(*mag_numbers[key])
            if not silent:
                print('\n Znucl:  ', key)
                # print(' '+s0.format(*mag_numbers[key]))
                # print(magnetic)
                print(' '+s.format(*magnetic))
            magnetic_all += magnetic
            if to_ox:
                if to_ox > 0:
                    ox = [to_ox-abs(m) for m in magnetic]
                else:
                    ox = [abs(m)+to_ox for m in magnetic]
                s2 = ' '.join(['{:5.1f}+']*len(magnetic))
                if not silent:
                    print(s2.format(*ox))
                    print('Average {:5.1f}+'.format(sum(ox)/len(ox)))


        return magnetic_all







    def set_magnetic_config(self, element, moments):
        #set magnetic configuration based on symmetry non-equivalent positions
        # element (str) - elements for which moments should be set 
        # moments (list) - list of moments for non-equivalent positions - the same order as in determine
        #                  the length should be the same as number of unique postions

        st = copy.deepcopy(self)

        if not hasattr(st, 'magmom') or None in st.magmom:
            magmom = [0.6]*st.natom
        else:
            magmom = st.magmom

        pos = st.determine_symmetry_positions(element)

        for j, p in enumerate(pos):
            for i in p:
                magmom[i] = moments[j]

        # print(magmom)
        st.magmom = magmom

        return st

       


    def convert2pymatgen(self, oxidation = None, slab = False, chg_type = 'ox'):
        """
        oxidation (dict) - {'Ti':'Ti3+'}
        if self.charges exist then it is used to update oxidation states of atoms
        chg_type - 'ox' - oxidation states calculated from self.charges 
                   'norm' - normalized self.charges
                    'tot' - just self.charges
                    'pot' - charges from potentials zval, number of valence electrons
                    'pm' - guess oxidation states using pymatgen, s

        slab - if True return slab object is returned - limited functional is implemented
        """
        from siman.analysis import calc_oxidation_states
        from siman.analysis import set_oxidation_states_guess

        st = self

        site_properties = {}

        if hasattr(self, 'magmom') and any(self.magmom):
            site_properties["magmom"] = self.magmom







        if oxidation is None:
            elements = self.get_elements()
        else:
            elements = [oxidation[el] for el in self.get_elements()]

        if slab:
            # print(dir(pymatgen.core))
            pm = Slab(self.rprimd, elements, self.xred, 
                miller_index = [0,0,1], oriented_unit_cell = None, shift = None, scale_factor = None,reorient_lattice = False,
                site_properties = site_properties)
        else:
            # print(elements)
            # print(len(self.xred))
            # print(site_properties)

            pm = pymatgen.Structure(self.rprimd, elements, self.xred, site_properties = site_properties)

        oxi = None
        

        if chg_type == 'pot':
            
            printlog('Using zval as charges', imp = '')
            oxi = [z*-1 for z in self.get_elements_zval()           ]
            print(oxi)
            pm.add_oxidation_state_by_site(oxi)
        elif chg_type == 'pm':
            # oxi = set_oxidation_states(st)
            pm.add_oxidation_state_by_guess()
            oxi = None
        
        else:
            if hasattr(self, 'charges') and any(self.charges):

                # chg_type = 'ox' # 'norm', 'tot'

                # print(chg_type)
                if chg_type == 'norm': #normalize charges
                    t = sum(self.charges)/len(self.charges)
                    oxi = [c-t for c in self.charges ]
                    # print(t)
                
                elif chg_type == 'ox':
                    # print(self.charges)
                    oxi = calc_oxidation_states(st = self)
                    # print(chg)
                    # sys.exit()

                elif chg_type == 'tot':
                    oxi = self.charges
            if hasattr(self, 'oxi_state') and self.oxi_state and any(self.oxi_state):
                #simply use predefined oxi_state
                
                oxi = self.oxi_state


                if 0: #check total charge
                    # st = st.copy()
                    chg = copy.copy(chg)
                    tot = sum(chg)
                    print('Total charge is ', tot, 'I subtract it uniformly from each atom')
                    d = tot/self.natom
                    chg = [c - d for c in chg]




        if oxi and not oxidation:
            # print(oxi)
            try:
                pm.add_oxidation_state_by_site(oxi)
            except:
                printlog('Warning! oxidation states were not added')




        return pm


    # pm = convert2pymatgen



    def get_pm_composition(self):
        ''
        pm = self.convert2pymatgen()
        cm = Composition(pm.formula)
        return cm


    def get_reduced_formula(self):
        ''
        pm = self.convert2pymatgen()
        cm = Composition(pm.formula)
        # cm = Composition(self.get_elements())
        return cm.reduced_formula

    def get_name(self):
        from siman.small_functions import latex_chem
        return latex_chem(self.get_reduced_formula())

    def get_formula(self):
        ''
        # pm = self.convert2pymatgen()
        # cm = Composition(pm.formula)
        # cm = Composition(self.get_elements())
        return self.convert2pymatgen().formula

    def get_reduced_composition(self):
        ''
        pm = self.convert2pymatgen()
        cm = Composition(pm.formula)
        # cm = Composition(self.get_elements())
        return cm.reduced_composition

    def get_reduced_formula_and_factor(self):
        pm = self.convert2pymatgen()
        cm = Composition(pm.formula)
        # cm = Composition(self.get_elements())
        return cm.get_reduced_formula_and_factor()

    def get_fractional_composition(self):
        pm = self.convert2pymatgen()

        cm = Composition(pm.formula)
        # cm = Composition(self.get_elements())
        return cm.fractional_composition





    def update_types(self, elements):
        # update typat, ntypat, znucl, nznucl from elements - list of elements names
        st = copy.deepcopy(self)
        st.ntypat = len(set(elements))

        st.typat = []
        st.znucl = []
        unique = []
        curtyp = 0
        types = {}
        for el in elements:
            if el not in unique:
                curtyp += 1
                types[el] = curtyp
                if is_string_like(el):
                    z = invert(el)
                else:
                    z = el

                st.znucl.append(z)
                # nznucl.append(0)
                unique.append(el)
            st.typat.append(types[el])
        
        st.get_nznucl()

        # print(st.ntypat, st.typat, st.nznucl, st.znucl)

        return st

    def update_from_pymatgen(self, stpm):
        """
        stpm - pymatgen structure
        update the current structure from pymatgen structure
        stil experimental!!!!!

        TODO:


        """
        st = copy.deepcopy(self)
        if hasattr(stpm, '_lattice'):
            st.rprimd = [np.array(vec) for vec in stpm._lattice._matrix]
        else:
            st.rprimd = [[10,0,0],[0,10,0],[0,0,10]] #temporary workaround
        # for site in stpm._sites:
            # print(dir(site))
        if hasattr(stpm._sites[0], '_frac_coords'):
            st.xred   = [np.array(site._frac_coords) for site in stpm._sites]
            st.update_xcart()
        else:
            st.xcart   = [np.array(site.coords) for site in stpm._sites]
            st.xred = [None]
        # print(elements)

        s = stpm._sites[0]

        if 'magmom' in s.properties:
            for i in stpm._sites:
                if 'magmom' not in i.properties:
                    i.properties['magmom'] = 0
            st.magmom = [site.properties['magmom'] for site in stpm._sites]
        else:
            st.magmom = [None]

        # print( dir(s._lattice) )
        # print( dir(s) )
        # print( s.properties )
        # print( s._properties )
        # print( dir(s.specie) )
        # print( s.specie.name )
        # sys.exit()
        elements = [s.specie.name for s in stpm._sites]
        # print(s.specie.oxi_state)
        if hasattr(s.specie, 'oxi_state'):
            oxi_state = [s.specie.oxi_state for s in stpm._sites]
            st.oxi_state = oxi_state
        
        # if hasattr(s.specie, 'oxi_state'):
        #     charges = [s.specie.oxi_state for s in stpm._sites]
        #     st.charges = charges


        # print(st.charges)
        # else:
        #     charges = [None]
        # print(elements)
        # print(charges)
        # sys.exit()

        st = st.update_types(elements)

        # print(len(st.typat))
        st.natom = len(st.typat)
        # sys.exit()

        if st.natom != len(st.xcart):
            printlog('Error! number of atoms was changed, please improve this method')

        # print(st.xcart)


        st.name+='_from_pmg'
        return st


    def rotate(self, axis, angle):
        """
        Rotate relative to cartesian coordinates
        axis - list of 3 elements, [0,0,1] in  cartesian coordinates
        angle in degrees
        """
        from pymatgen.transformations.standard_transformations import RotationTransformation
        
        st = copy.deepcopy(self)
        rot = RotationTransformation(axis, angle)
        stpm = st.convert2pymatgen()
        stpmr1 = rot.apply_transformation(stpm)
        st_r1 = st.update_from_pymatgen(stpmr1)
        st_r1 = st_r1.return_atoms_to_cell()
        return st_r1

    def align_with_axes(self):
        """
        Align with coordinate axes
        First vector with x
        Third vector with z

        """
        from siman.small_functions import normal, angle
        r_orig = self.rprimd
        st = self.copy()
        r = st.rprimd

        st = st.rotate(normal(r[0], [1,0,0]), angle(r[0], [1,0,0])) #combine with x and y
        r = st.rprimd

        if angle(r[2], [0,0,1]) > 1:
            st = st.rotate(normal(r[2], [0,0,1]), angle(r[2], [0,0,1])) #combine with x and y
            r = st.rprimd

        with np.printoptions(precision=2, suppress=True):
            print('rprimd before:\n',np.array(r_orig))
            print('rprimd after :\n',np.array(r))
        return st

    def invert_axis(self, axis):
        #invert one vector
        st = copy.deepcopy(self)

        st.rprimd[axis] *= -1
        st.update_xred()
        st = st.return_atoms_to_cell()
        return st
    def invert_xred(self, axis):
        #invert xred coordinates along one vector axis
        st = copy.deepcopy(self)

        for i in range(st.natom):
            st.xred[i][axis] = 1 - st.xred[i][axis]
        st.update_xcart()
        # st = st.return_atoms_to_cell()
        return st

    def mirror(self, axis):
        #mirror along vector
        j = axis
        st = copy.deepcopy(self)
        for i, x in enumerate(st.xred):
            st.xred[i][j] = -x[j]
        st.update_xcart()
        st = st.return_atoms_to_cell()
        return st

    def sizes(self):
        #return sizes along x, y, and z

        # for x in self.xcart:
        xyz = list(map(list, zip(*self.xcart)))

        dx = max(xyz[0]) - min(xyz[0]) 
        dy = max(xyz[1]) - min(xyz[1]) 
        dz = max(xyz[2]) - min(xyz[2]) 
        return dx,dy, dz
    def rprimd_len(self):
        #return vector lengths
        r = self.rprimd
        n = np.linalg.norm
        return n(r[0]), n(r[1]), n(r[2])






    def pvec(self):
        """
        print primitive vectors in formatted way
        """
        with np.printoptions(precision=3, suppress=True):
            print( np.array(self.rprimd) )

        return


    def add_z(self, z):
        # method appends additional height to the cell
        # negative value of z appends to remove vacuum layer
        st = copy.deepcopy(self)

        st.rprimd[2][2] += z
        for i in st.xcart:
            i[2] += z
        st.update_xred()
        st = st.return_atoms_to_cell()
        return st


    def get_oxi_states(self, typ = 'charges'):
        """
        Create and return list of oxidation states from charges and valences
        self.charges should exist as full charges (e.g. from Bader analysis)
        
        INPUT:
            typ (str) 
                'charges' - from charges and zval
                'guess'   - from guess
        
        RETURN
            oxi (list) - list of oxidation states for each atom
        """
        st = self
        if typ == 'charges':
            printlog('Using zval as reference', imp = '')
            st = self
            oxi = []
            for j, z_val, el in zip(range(st.natom), st.get_elements_zval(), st.get_elements()):
                oxi.append( z_val - st.charges[j] )
            # st.oxi_state = oxi
        
        elif typ == 'guess':
            pm = st.convert2pymatgen()
            pm.add_oxidation_state_by_guess()
            st = st.update_from_pymatgen(pm)
            oxi = st.oxi_state


        return oxi


    def generate_charge_orders(self, el, states = None, x = 0.5):
        """
        Method generates charge order for provided ion, oxidation states, and ratio (not realized yet)
        
        INPUT:
            el (str) - element with charge order
            states (tuple) - two possible charge states  (e.g. +2 and +4)
            x (float) - concentration of ions with state[0] charge state 
        RETURN:
            oxi_states (list of lists) - list of oxi_state lists
        """


        def order(ls, i):
            """
            Find recursivly all possible orderings for the given x
            ls - initial list of atoms 
            i - index in ls  

            """
            # print(i)
            for s in 1,-1:
                
                ls[i] = s
                
                if i < len(ls)-1:
                
                    order(ls, i+1)
                
                else:
                    # print (ls.count(-1)/tot - x)
                    if abs(ls.count(-1)/tot - x ) < 0.001:
                        orderings.append(copy.deepcopy(ls) )  
            return

        st = self

        iels = st.get_specific_elements([invert(el)])
        oxi_state = st.oxi_state

        oxi_states = []
        orderings = []
        tot = len(iels)
        ls = [0]*tot
        # print(ls)
        order(ls, 0)
        
        print('Total number of charge orderings for x=',x,'is',len(orderings))

        for order in orderings:
            atoms_with_minor = [i for i, s in enumerate(order) if s < 0]
            # atoms_with_major = [i for i, s in enumerate(order) if s > 0]
            # print(atoms_with_minor)
            for iloc in range(tot):
                i = iels[iloc]
                if iloc in atoms_with_minor:
                    oxi_state[i] = states[0]
                else:
                    oxi_state[i] = states[1]
            oxi_states.append(copy.copy(oxi_state))
            # print(oxi_state[0:8])
        return oxi_states

    def get_conventional_cell(self):
        """
        return conventional cell 
        """
        st_mp = self.convert2pymatgen()

        # st_test = self.update_from_pymatgen(st_mp)
        # st_test.printme()

        sf = SpacegroupAnalyzer(st_mp, ) #symprec = 0.1

        sc = sf.get_conventional_standard_structure() # magmom are set to None

        # print(sc)
        st = self.update_from_pymatgen(sc)

        # print(st.rprimd)
        # print(len(st.xcart))
        # print(st.ntypat)

        return st




    def get_primitive_cell(self, international_monoclinic = True):
        """
        return primitive cell 
        """
        st_mp = self.convert2pymatgen()

        sf = SpacegroupAnalyzer(st_mp, ) #symprec = 0.1

        sc = sf.get_primitive_standard_structure(international_monoclinic=international_monoclinic) # magmom are set to None

        st = self.update_from_pymatgen(sc)

        return st

    def get_refined_structure(self, ):
        """
        Get the refined structure based on detected symmetry. 
        The refined structure is a conventional cell setting with atoms moved to the expected symmetry positions.
        """
        st_mp = self.convert2pymatgen()

        sf = SpacegroupAnalyzer(st_mp, ) #symprec = 0.1

        sc = sf.get_refined_structure() # magmom are set to None
        # sc = sf.get_symmetrized_structure() # magmom are set to None


        st = self.update_from_pymatgen(sc)

        return st


    def get_symmetry_operations(self, symprec = 0.1):
        """
        Return symmetry operations as a list of SymmOp objects. 
        By default returns fractional coord symmops. But cartesian can be returned too.
        """
        st_mp = self.convert2pymatgen()

        sf = SpacegroupAnalyzer(st_mp, symprec = symprec)

        sym_op = sf.get_symmetry_operations(cartesian=False)

        return sym_op




    def get_surface_pos(self, reduced = False):
        """
        Allows to return positions of bottom and top surfaces (edge atoms) in cartesian
        assumed normal to R3
        small number is added or subtracted to/from edge atom to overcome further numericall errors
            
            reduced - reduced coordinates

        returns list of two z coordinates [bottom, top]
        
        """
        st = self

        z1 = 100
        z2 = -100
        z = []
        # print(st.xcart)
        if reduced:
            l = st.xred
        else:
            l = st.xcart
        for x in l:
            if z1 > x[2]:
                z1 = x[2]
            if z2 < x[2]:
                z2 = x[2]

        # z1-=0.01
        # z2+=0.01
        printlog('Surfaces are ', z1, z2)

        z.append(z1)
        z.append(z2)

        return z


    def get_surface_atoms(self, element = None, surface = 0, surface_width = 0.5 ):
        """
        return numbers of surface atoms
            elememt (str) - which element is interesting? can be CoO to take two elements; if None than all elements are takes
            surface_width - which atoms to consider as surface 
            surface (int) - 0 or 1 - one of two surfaces; surface 1 has smaller z coordinate
            TODO - 
            now only along third vector, make more general
            PBC may work incorrect
        """

        st = self
        surface_atoms = [[],[]]
        
        z = st.get_surface_pos(reduced  = True)


        els = st.get_elements()

        suf_width_red = surface_width/np.linalg.norm(st.rprimd[2])


        for i, x in enumerate(st.xred):
            el = els[i]

            if element is None or el in element:
                # print(x[2])
                if z[0] <= x[2] < z[0]+suf_width_red:
                    surface_atoms[0].append(i)

                if z[1] - suf_width_red  < x[2] <= z[1]:
                    surface_atoms[1].append(i)


        return surface_atoms[surface]


    def cover_surface(self, el, d = 1):
        """
        Covers both surfaces of the slab with element el, 
        each surface atoms on top positions
        
        the third vector should be parallel to z axis

        el - element name 
        d - distance in A

        """


        suf0_at = self.get_surface_atoms(surface = 0)
        suf1_at = self.get_surface_atoms(surface = 1)
        # print(suf_at)
        st = self.copy()

        coords = []
        for i in suf0_at:
            coords.append(st.xcart[i]-np.array([0,0,d]))
        for i in suf1_at:
            coords.append(st.xcart[i]+np.array([0,0,d]))

        st = self.add_atoms(coords, element = el)

        return st 


    def  if_surface_atom(self, i):
        """
        Based on reduced coordinates, 
            i - number of atom from zero 
        """
        st = self
        suf = st.get_surface_pos(reduced = 1)
        sufd = min(abs(st.xred[i][2]-suf[0]),abs(st.xred[i][2]-suf[1]))
        if sufd < 0.5/np.linalg.norm(st.rprimd[2]):
            printlog('TM is on surface, sufd', sufd , imp = 'y')
            return True


    def get_surface_area(self):
        """
        currently should be normal to rprim[0] and rprim[1]
        """
        st = self
        # print()
        return np.linalg.norm( np.cross(st.rprimd[0] , st.rprimd[1]) )


    def printme(self):
        print(self.convert2pymatgen())
        return 


    def get_space_group_info(self, symprec = None):
        
        default = 0.01
        if not symprec:
            symprec = default
        # print(symprec)
        if hasattr(self, 'spg') and symprec == default:
            spg = self.spg
        else:
            # print('get')
            p = self.convert2pymatgen()
            spg = p.get_space_group_info(symprec, angle_tolerance=5.0)
        # p = self.convert2pymatgen()

        # print(p.get_symmetry_operations(symprec))


        return spg

    def sg(self,symprec = None, silent = 0):
        try:
            s = self.get_space_group_info(symprec)
        except:
            s = 'error'
        if not silent:
            print(s)
        return s

    def get_angles(self):
        R = self.rprimd
        alpha = angle(R[1], R[2])
        beta  = angle(R[0], R[2])
        gamma = angle(R[0], R[1])
        return alpha, beta, gamma

    # def __str__(self):
    #     # print(self.convert2pymatgen())
    #     return 


    def get_element_xred(self, element):
        """
        Get xred of *element* first occurance
        """
        i = self.get_elements().index(element)
        return self.xred[i]
    
    def get_element_xcart(self, element):
        """
        Get xred of *element* first occurance
        """
        i = self.get_elements().index(element)
        
        return self.xcart[i]

    def get_transition_elements(self, fmt = 'names'):
        """Returns list of transition elements (chemical names or z) in the structure
        fmt - 
            'names'
            'z'
            'n' - numbers of atoms
        """
        el = self.get_elements()
        tra = []
        ns = []
        for i, e in enumerate(el):
            for t in header.TRANSITION_ELEMENTS:
                if e == invert(t):
                    tra.append(e)
                    ns.append(i)
        
        if fmt == 'z':
            tra = [invert(t) for t in tra]
        elif fmt == 'n':
            tra = ns
        return tra


    def get_specific_elements(self, required_elements = None, fmt = 'n', z_range = None, zr_range = None):
        """Returns list of specific elements (chemical names. z, or numbers from 0) in the structure
        required_elements - list of elements z of interest
        z_range - (2 index tuple) range of z coordinates in A: only atoms from z1 to z2 are taken
        fmt - format of output
            'names'
            'z'
            'n' - numbers of atoms
            'x' - xcart
            'xr' - xred
        


        """
        el = self.get_elements()
        tra = []
        ns = []
        r3 = np.linalg.norm(self.rprimd[2])

        if zr_range is None and z_range is not None:
            zr_range = [z_range[0]/r3, z_range[1]/r3]

        if zr_range:
            def additional_condition(xr):
                return zr_range[0] < xr <= zr_range[1]
        else:
            def additional_condition(xr):
                return True

        xcart = []
        xred = []
        for i, e, xc, xr in zip( range(self.natom), el, self.xcart, self.xred ):
            Z = invert(e)
            if Z in required_elements and additional_condition(xr[2]):
                tra.append(e)
                ns.append(i)
                xcart.append(xc)
                xred.append(xr)
        
        if fmt == 'z':
            tra = [invert(t) for t in tra]
        elif fmt == 'n':
            tra = ns
        elif fmt == 'x':
            tra = xcart
        elif fmt == 'xr':
            tra = xred

        return tra


    def get_dipole(self, ox_states = None, chg_type = 'ox'):
        """ Return dipole moment in e*A calculated by pymatgen
        ox_states (dict) - oxidation states of elements
        chg_type (str) - type of charges  if provided in self.charges; see description of self.convert2pymatgen()
        
        If you need to convert e*A to debye (D), use 1 D = 0.20819434 eA = 3.33564×10−30 Cm; 1 eA = 4.8 D
        """
        slab = self.convert2pymatgen(slab = 1, oxidation = ox_states, chg_type = chg_type)
        return slab.dipole


    def add_atoms(self, atoms_xcart = None, element = 'Pu', return_ins = False, selective = None, atoms_xred = None, mag = None):
        """
        appends at the end if element is new. Other case insertered according to VASP conventions
        Updates ntypat, typat, znucl, nznucl, xred, magmom and natom
        atoms_xcart (list of ndarray)
        atoms_xred (list of coordinate lists) - if provided both, this has higher priority 

        selective (list of lists) - selective dynamics

        mag magnetic moment of added atoms, if None, than 0.6 is used
            magmom is appended with 0.6, 
            please improve me! by using the corresponding list of magmoms


        if return_ins:
            Returns Structure(), int - place of insertion of first atom
        else:
            Structure()
       


        """

        printlog('self.add_atoms(): adding atom ', element, imp = 'n')

        st = copy.deepcopy(self)

        # print(st.select)
        if selective: 
            if not hasattr(st, 'select') or st.select is None or len(st.select) == 0:
                st = st.selective_all()
        # print(st.select)
        # sys.exit()

        if atoms_xred is not None:
            atoms_xcart = xred2xcart(atoms_xred, st.rprimd)


        natom_to_add = len(atoms_xcart)
        if natom_to_add == 0:
            printlog('self.add_atoms(): Nothing to add, continiue')
            return st


        st.natom+=natom_to_add

        # print(element)
        # sys.exit()
        if type(element) is not str:
            printlog('Error! element is', element, 'but should be element name')

        el_z_to_add = element_name_inv(element)

        if hasattr(st, 'magmom') and any(st.magmom) or mag:
            magmom_flag = True
        else:
            magmom_flag = False

        if mag is None:
            mag = 0.6

        if el_z_to_add not in st.znucl:
            
            st.znucl.append( el_z_to_add )
            
            st.nznucl.append(natom_to_add)            

            st.ntypat+=1

            # print (st.typat)
            if st.typat:
                typ = max(st.typat)+1
            else:
                typ = 1
        
            st.xcart.extend(atoms_xcart)
            
            st.typat.extend( [typ]*natom_to_add )

            if selective is not None:
                st.select.extend(selective)
            elif hasattr(st, 'select') and st.select and len(st.select) > 0:
                # printlog('adding default selective', imp = 'y')

                st.select.extend( [[True,True,True] for i in range(natom_to_add)] )
            else:
                ''

            if magmom_flag:
                # print(mag, natom_to_add)

                # print('add_atoms: magmom = ',st.magmom)
                st.magmom.extend( [mag]*natom_to_add  )
                # print('add_atoms: magmom = ',st.magmom)

            j_ins = self.natom # first one

        else:
            i = st.znucl.index(el_z_to_add)
            
            st.nznucl[i]+=natom_to_add
            
            typ = i+1

            # print(el_z_to_add, typ)
            for j, t in enumerate(st.typat):
                if t == typ:
                    j_ins = j+1 # place to insert




            st.xcart[j_ins:j_ins] = atoms_xcart
            
            st.typat[j_ins:j_ins] = [typ]*natom_to_add

            # print(st.select)
            # sys.exit()
            if selective is not None:
                printlog('adding selective', imp = '')

                st.select[j_ins:j_ins] = selective
            elif hasattr(st, 'select') and  st.select and len(st.select) > 0:
                printlog('adding default selective', imp = '')
                st.select[j_ins:j_ins] = [[True,True,True] for i in range(natom_to_add)]
            else:
                ''

            if magmom_flag:
                st.magmom[j_ins:j_ins] =  [mag]*natom_to_add


        st.xcart2xred()
        


        # print(st.select)



        if return_ins:
            return st, j_ins
        else:
            return st


    def add_atom(self, xr = None, element = 'Pu', xc = None, selective = None):
        """
        
        allows to add one atom using reduced coordinates or cartesian
        xr - reduced
        xc - cartesian
        element - element name
        """

        if xr is not None:
            ''
            xc = xred2xcart([xr], self.rprimd)[0]
        
        elif xc is not None:
            ''

        else:
            ''
            printlog('Error! Provide reduced *xr* or cartesian *xc* coordinates!')


        if selective is not None:
            selective = [selective]

        st = self.add_atoms([xc], element = element, selective = selective)
        return st 







    def reorder_for_vasp(self, inplace = False):
        """
        
        Group and order atoms by atom types; consistent with VASP
        return st
        """
        ''
        if inplace:
            st = self
        else:
            st = copy.deepcopy(self)
        nt = range(st.ntypat)
        zxred  = [[] for i in nt]
        zxcart = [[] for i in nt]
        ztypat = [[] for i in nt]
        zmagmom= [[] for i in nt]
        ziat =   [[] for i in nt]
        i = 0
        # print(st.ntypat)
        for t, xr, xc in zip(st.typat, st.xred, st.xcart):
            # print ("t ", t, xr)
            zxred[ t-1].append(xr)
            zxcart[t-1].append(xc)
            ztypat[t-1].append(t)
            ziat[t-1].append(i)
            i+=1

        st.nznucl = [len(typat) for typat in ztypat]

        st.xcart = [item for sublist in zxcart for item in sublist]
        st.xred  = [item for sublist in zxred for item in sublist]
        st.typat = [item for sublist in ztypat for item in sublist]
        original_numbers  = [item for sublist in ziat for item in sublist]
        
        st.perm = [original_numbers.index(i) for i in range(st.natom)] # show the initial order of atoms; starting from 0

        if hasattr(st, 'magmom') and any(st.magmom):
            for t, m in zip(st.typat, st.magmom):
                zmagmom[t-1].append(m)
            st.magmom = [item for sublist in zmagmom for item in sublist]

        else:
            st.magmom = [None]

        # print(st.get_elements())

        # print(st.perm)

        

        return st


    def reorder_element_groups(self, order = None, inplace = False):
        """
        
        Group and order atoms by atom types; consistent with VASP
        order (list) -required order e.g. ['O', 'Li'] or str 'alphabet'

        return st
        """

        st = copy.deepcopy(self)

        # for z in st.znucl:
        #     print(z)
        
        typat = []
        xcart = []
        magmom = []
        znucl = []
        # st.write_poscar()


        els = st.get_elements()

        if 'alphabet' in order:
            order =  list(sorted(set(els)))
        # print(unique_sorted)
        # sys.exit()

        old_numbers = []
        t = 1
        for el in order:
            if el not in els:
                printlog('Error! Check *order* list')
            
            znucl.append( invert(el) )

            for i in range(st.natom):
                # old_numbers
                el_i = els[i]
                if el_i not in order:
                    printlog('Error! Check *order* list')

                if el_i == el:
                    # print(el)
                    typat.append(t)
                    xcart.append(st.xcart[i])
                    old_numbers.append(i)
                    if None not in st.magmom:
                        magmom.append(st.magmom[i])
            t+=1

        if len(magmom) == 0:
            magmom = [None]

        st.old_numbers = old_numbers
        st.xcart = xcart
        st.magmom = magmom
        st.typat = typat
        st.znucl = znucl
        st.update_xred()
        st.name+='_r'
        # st.write_poscar()

        return st

    def reorder(self, new_order):
        """
        Reorder according to new_order list, len(new_order) should be equal to natom

        """
        st = self.copy()
        els = self.get_elements()
        st = st.remove_atoms(atoms_to_remove = els,  clear_magmom=0)
        # print(st.natom)

        magmom_flag = False
        if len(self.magmom) == self.natom:
            magmom_flag = True
        # print(st.magmom)

        for i in new_order:
            x = self.xcart[i]
            el = els[i]

            if magmom_flag:
                m = self.magmom[i]
            else:
                m = None
            # print(m)
            st = st.add_atoms([x], el, mag = m)

        return st



    def permutate_to_ref(self, st_ref):
        """
        Permutate atom numbers of self according to st
        Structures should have the same amount of atoms and be quite similar
        the self can have one extra type of atoms

        #TODO: now extra atom types could be only in self structure

        """

        st = self.copy()
        els = st.get_elements()
        els_ref = st_ref.get_elements()
        # print(els, els_ref)
        extra = list(set(els)-set(els_ref))[0] # only one is implemented currently
        new_order = []
        for el1, x1 in zip(els_ref, st_ref.xcart):
            i,s,d = st.find_closest_atom(x1)
            # print(el1, st.get_elements()[i])
            new_order.append(i)
        
        for i in st.get_specific_elements([1]):
            # print(i, els[i])
            new_order.append(i)

        # sys.exit()


        if len(new_order) != st.natom:
            printlog('Error! something is wrong with number of atoms')

        # print('ref ', st_ref.get_elements())

        # print('init', st.get_elements())
        st = st.reorder(new_order)
        # print('reor', st.get_elements())

        return st




    def del_atom(self, iat):
        """
        Now can delete only one atom with number iat (int), starting from 0. 
        Takes care of magmom, ntypat, typat, znucl, nznucl, xred and natom
        Returns Structure()
        """


        # print_and_log('Warning! Method del_atoms() was not carefully tested ')
        st = copy.deepcopy(self)
        # print(st.nznucl)

        i = iat

        typ = st.typat[i]

        # printlog('del_atom(): I remove atom ',  st.get_elements()[i], imp = 'n')
        del st.typat[i]
        del st.xred[i]
        del st.xcart[i]

        # print ('Magmom deleted?')
        # print(st.magmom)
        if hasattr(st, 'magmom') and any(st.magmom):
            del st.magmom[i]
            # print ('Yes!')
        else:
            ''
            # print ('No!')



        st.natom-=1

        if typ in st.typat:
            st.nznucl[typ-1]-=1
        else:
            del st.nznucl[typ-1]
            del st.znucl[typ-1]
            st.ntypat-=1

            # for i, n in enumerate(st.nznucl):
            #     typ = i+1
            #     st.typat = [typ, ]

            for i, t in enumerate(st.typat):
                if t > typ:
                    st.typat[i]-=1


        # print(st.nznucl)

        return st







    def mov_atoms(self, iat = None, to_x = None, to_xr = None, relative = False):
        """
        Move one atom to xcart position *to_x*
        relative (bool) - if shift is relative

        """
        st = copy.deepcopy(self)
        
        if to_xr is not None:
            if relative:
                st.xred[iat] += to_xr
            else:
                st.xred[iat] = to_xr
            
            st.xred2xcart()



        else:
            if relative:
                st.xcart[iat] += to_x
            else:
                st.xcart[iat] = to_x
            
            st.xcart2xred()

        return st

    def swap_atoms(self, iat1, iat2):
        
        st = copy.deepcopy(self)
        els = st.get_elements()
        # printlog('You choose', els[iat1], 'and', els[iat2])

        x1 = st.xcart[iat1]
        st.xcart[iat1] = st.xcart[iat2]
        st.xcart[iat2] = x1
        st.xcart2xred()
        
        return st


    def leave_only(self, atom_type = None):
        #Remove all atoms except *atom_type*(str, mendeleev element name)
        
        print_and_log('Starting leave_only()', imp = 'n')

        st = copy.deepcopy(self)
        
        print_and_log('    N of atoms before = ',st.natom, imp = 'n')


        z = element_name_inv(atom_type)

        new_xred = []
        new_magmom = []
        
        if hasattr(st, 'magmom') and any(st.magmom):
            for t, xr, m in zip(st.typat, st.xred, st.magmom):
                if st.znucl[t-1] == z:
                    new_xred.append(xr)
                    new_magmom.append(m)
        else:
            for t, xr in zip(st.typat, st.xred):
                if st.znucl[t-1] == z:
                    new_xred.append(xr)


        st.magmom = new_magmom

        st.xred = new_xred

        st.natom = len(new_xred)

        st.ntypat = 1

        st.typat = [1]*st.natom

        st.znucl = [z,]

        st.nznucl = [st.natom,]

        st.xcart = xred2xcart(st.xred, st.rprimd)

        # print st.xred

        print_and_log('    N of atoms after  = ',st.natom, imp = 'n')


        return st


    def get_numbers(self, element):
        "return numbers of specific element "
        return [i for i, el in enumerate(self.get_elements()) if el == element]


    def remove_atoms(self, atoms_to_remove, from_one = 0, clear_magmom  = 1 ):
        """
        remove atoms either of types provided in *atoms_to_remove* or having numbers provided in *atoms_to_remove*, starting from 0
        st (Structure)
        atoms_to_remove (list) - list of element names or numbers
        from_one (int)- if 1 the numbers of atoms in provided list are starting from one
        clear_magmom - by default magmom is cleared
        """
        st = copy.deepcopy(self)
        # print(st.nznucl)

        numbers = list(range(st.natom))


        atom_exsist = True

        while atom_exsist:

            for i, (n, el) in enumerate(  zip(numbers, st.get_elements()) ):
                # print(i)

                if el in atoms_to_remove or n+from_one in atoms_to_remove:
                    # print(n)
                    # atoms_to_remove.remove(i)
                    st = st.del_atom(i)
                    del numbers[i]

                    break
            else:
                atom_exsist = False
        printlog('remove_atoms(): Atoms', atoms_to_remove, 'were removed')
        if clear_magmom:
            st.magmom = [None]

        # print(st.nznucl)

        # print(st.get_elements())
        return st


    def del_layers(self, xred_range = None, xcart_range = None):
        """
        remove atoms normal to R3 in given range

        xred_range (list) [from, to]
        highlight - replace with Pu to check

        """
        st = copy.deepcopy(self)
        # print(st.nznucl)
        r3 = np.linalg.norm(st.rprimd[2])
        # print('r3',r3 )
        xcr = xcart_range
        if xcr:
            xred_range = [xcr[0]/r3, xcr[1]/r3]
            printlog('xcart_range converted to xred', xred_range)
            
            # print('xcart_range converted to xred',xcart_range, xred_range)

        dels = []
        for i, xr in enumerate(st.xred):
            if xred_range[0]  < xr[2] <= xred_range[1]:
                # print(xred_range[0], xr[2], xred_range[1])
                dels.append(i)
        # print(dels)
        st = st.remove_atoms(dels)
        st.name+='_del'
        # print(st.nznucl)
        
        return st


    def remove_at_in_zrange(self, z_range, del_range = 0):
        # z_range - at 2 remove
        # del_range: if true  - remove atoms in z_range, if false - remove atoms out of given z_range

        st = copy.deepcopy(self)
        at2remove = []
        
        for i in range(0, st.natom):

            if (z_range[1]<st.xcart[i][2] or st.xcart[i][2]<z_range[0]) and not del_range:
                at2remove.append(i)

            if z_range[0]<=st.xcart[i][2]<=z_range[1]  and del_range:
                at2remove.append(i)

        st_new = st.remove_atoms(at2remove)
        return st_new  



    def replace_atoms(self, atoms_to_replace, el_new, mag_new = None, silent = 1, mode = 1):
        """
        atoms_to_replace - list of atom numbers starting from 0
        el_new - new element periodic table short name
        mag_new - new magnetic moment
        mode 
            1 - old behaviour, numbering is not conserved
            2 - numbering is conserved if el_new already exists in self

        TODO:
        Now if el_new already exists in structure, numbering is conserved,
        otherwise numbering is not conserved.
        Make numbering conservation in case when new element is added
        Both modes can be useful, as the first case is compat with VASP


        """
        st = copy.deepcopy(self)

        numbers = list(range(st.natom))
        z_new = invert(el_new)
        atom_exsist = True
        if silent:
            warn = 'n'
        else:
            warn = 'Y'


        while atom_exsist:


            for i, (n, el) in enumerate(  zip(numbers, st.get_elements()) ):
                # print(i)

                if n in atoms_to_replace:
                    xcart = st.xcart[i]

                    if mode == 2:
                        if st.get_elements().count(el) == 1:
                            printlog('Error! The functions replace_atoms() in mode == 2 works incorrectly if one atom of type ', el)
                        if el_new not in st.get_elements():
                            st.znucl.append(z_new)
                            # print(st.nznucl)
                            st.ntypat+=1
                        
                        it = st.znucl.index(z_new)+1
                        
                        st.typat[n] = it
                        # print('mag_new',mag_new)
                        # sys.exit()
                        st.magmom[n] = mag_new
                        # print(it, z_new, st.typat)
                        # print(st.get_elements())
                        # sys.exit()
                        del numbers[i]

                        printlog('replace_atoms(): atom', i, el, 'replaced with', st.get_elements()[n], '; Atom number stayed the same' , imp = warn)

                    else:
                        # atom number is changed, since new typat is added
                        el_rep = st.get_elements()[i]
                        st = st.del_atom(i)
                        del numbers[i]
                        # print(mag_new)
                        st = st.add_atoms([xcart], element = el_new, mag = mag_new)
                        printlog('replace_atoms(): atom', i, el_rep, 'replaced with', el_new, '; Atom number was changed', imp = warn)
                        # print('replace_atoms(): mag, magmom', mag_new, st.magmom)
                    break
            else:
                atom_exsist = False
        st.get_nznucl()
        # printlog('remove_atoms(): Atoms', atoms_to_remove, 'were removed')

        # print(st.get_elements())
        return st


    def replace_atoms2(self, el_old, el_new, concentration):
        """
        Replace atoms using random  

        el_old - element to replace

        el_new - new element

        concentration - part of atoms el_old to replace by el_new. Number from 0 to 1.
        """
        import random

        st = copy.deepcopy(self)

        numbers = list(range(st.natom))

        nums = []


        for i, (n, el) in enumerate(  zip(numbers, st.get_elements()) ):
            if el == el_old: nums.append(n)
            # print(nums)
        n_replace = int(len(nums)*concentration)
        c = float(n_replace/len(nums))
        if c != concentration: print('\n\nAttention! This concentraiton is impossible. Real concentration is - ', c) 
        print('\nI have found {} from {} random atoms of {} to replace by {} \n'.format(n_replace, len(nums), el_old, el_new ))
        random.shuffle(nums)

        atoms2replace = nums[0:n_replace]
        # print(num2replace)
        st = st.replace_atoms(atoms2replace, el_new)

        return st



    def remove_part(self, element, new_conc):
        """
        element to remove
        new_conc <1 - new concentration of element atoms (part of unity)
        """
        st = copy.deepcopy(self)


        numb = self.get_numbers(element)

        nat_el = int(np.ceil((len(numb)*new_conc)))
        printlog('New number of ', element, 'atoms is ', nat_el, imp = 'y')

        del_num = numb[nat_el:len(numb)]



        return st.remove_atoms(del_num)




    def add_vacuum(self, vector_i, thickness):
        """
        Allows to add or remove vacuum along one of the rprimd vectors
        vector_i (int) - index of vector along which vacuum should be added 0, 1, 2
        thickness (float) - thickness of added (positive) or removed (negative) vacuum 


        TODO:
        add capability to add vacuum normal to surface in case of non-orthogonal cells

        """
        st = copy.deepcopy(self)
        v = st.rprimd[vector_i]
        v_l = np.linalg.norm(v)
        new_len = v_l+thickness

        st.rprimd[vector_i]*=new_len/v_l

        st.update_xred()
        st.name+='_vac'
        # st.write_xyz()
        return st


    # def sum_of_coord(self):
    #     sumx = 0
    #     for x in self.xcart:
    #         sumx+=x
    #     sumx/=len(self.xcart)
    #     return sumx

    def return_atoms_to_cell(self, shift = 0):
        #shift - shift from the end of vectors between 0 and 1 - allows to collect atoms close to origin

        st = copy.deepcopy(self)
        bob = 0-shift; upb = 1-shift;
        n = 0 
        # print st.xred
        for xr in st.xred:
            for j in 0,1,2:
                if xr[j]  < bob:  
                    xr[j] = xr[j] - int(xr[j]) + 1 #allows to account that xr can be more than 2
                if xr[j]  >= upb:  
                    # print(xr[j], int(xr[j]))
                    xr[j] = xr[j] - int(xr[j]) 
        # n+=1
        # zmin = 100
        # for xr in st.xred:
        #     if xr[2]<zmin: zmin = xr[2]
        # if zmin < 0:
        #     for xr in st.xred:
        #         xr[2] = xr[2]-zmin



        st.xcart = xred2xcart(st.xred, st.rprimd)

        # print_and_log(str(n)+" atoms were returned to cell.\n")
        #print st.xred
        return st





    def combine(self, st_list, only_numbers = None):
        """
        Combine several structures into one
        using reduced coordinates
        """
        st_b = self.copy()

        if only_numbers is None:
            only_numbers = []

        for i, st in enumerate(st_list):
            # print(i)

            for j, xr, el in zip(list(range(st.natom)), st.xred, st.get_elements() ):
                if j in only_numbers:
                    st_b = st_b.add_atom(xr, el)


        return st_b




    def combine_atoms(self, d = 0.1):
        """
        Combine close-lying atoms into one
        d (float) - all atoms with d less then d are combined into one in Angstrom
        """
        st = self
        # copy()
        remove_list = []
        for i, x1 in enumerate(st.xcart):
            for ii, x2 in enumerate(st.xcart[i+1:]):
                j = ii+i+1
                dx = st.distance(x1=x1, x2=x2)
                xlist = []
                if dx < d:
                    xlist.append(x2)
                    remove_list.append(i)
                    # print(dx)
            if xlist:
                x_new = sum(xlist)/len(xlist)
                # print(x_new)
                st.xcart[i] = x_new
                print('Atom i= {:n}, {:s} replaced with average position'.format(i, st.get_elements()[i]) )

        st = st.remove_atoms(remove_list)



        return st






    def perturb(self, d=0.1):
        """
        d is distance
        """
        st = self.copy()
        pm = st.convert2pymatgen()
        pm.perturb(d)
        st = st.update_from_pymatgen(pm)
        return st

    def find_atom_num_by_xcart(self, x_tar, prec = 1e-6, search_by_xred = 1 ):
        """take into account periodic conditions

        search_by_xred - use xred to search with PBC
        prec - difference in atomic positions in A; transformed to reduced difference using the longest vector
        
        TODO:
        make normal function that treats periodic boundary conditions normally!!!
            done please test

        """

        [xr_tar] = xcart2xred([x_tar], self.rprimd)
        printlog('find_atom_num_by_xcart(): xr_tar = ', xr_tar)
        #PBC!!!
        if 1:
            for i in [0,1,2]:
                if xr_tar[i] < 0:
                    xr_tar[i]+= 1
                    
                if xr_tar[i] >= 1:
                    xr_tar[i]-= 1
        
        printlog('find_atom_num_by_xcart(): xr_tar after periodic = ', xr_tar)

        # print(xr_tar)
        # print(self.rprimd)
        [x_tar] = xred2xcart([xr_tar], self.rprimd)
        
        printlog('find_atom_num_by_xcart(): x_tar after periodic = ', x_tar)
        # sys.exit()
        # print(x_tar)
        self = self.return_atoms_to_cell() # please solve the problem, as neb not always works correctly!

        # for i, x in enumerate(self.xcart): # xcart
        #     if np.linalg.norm(x-x_tar) < prec:
        #         printlog('Atom', i+1, 'corresponds to', x_tar)
        #         return i
        
        vmax = max(self.vlength)
        prec_xr = prec/vmax
        print('Reduced precision is', prec_xr, '')

        i_matched = None
        d_min = 100

        for i, x in enumerate(self.xred): # xred
            # d = np.linalg.norm(x-xr_tar)
            # print(d)
            # if d >= 0.5:
                # d = abs(d-1)
            dv = []
            for j in 0,1,2:
                di = abs(x[j]-xr_tar[j])
                if di >= 0.5:
                    di=di-1
                dv.append(di)
            dvl = np.linalg.norm(dv)

            # print(dvl)
            if dvl < prec_xr:
                print('Difference is ', dvl*vmax, 'A', 'for atom', i)
                if dvl < d_min:
                    i_matched = i
                    d_min = dvl
                    d_xc = np.linalg.norm(self.xcart[i]-x_tar)

                    printlog('Atom', i, self.xcart[i], 'corresponds to the requested atom with difference', d_xc, 'A',  imp = 'Y')



        if i_matched is None:
            printlog('Attention, atom ', x_tar, 'was not found' )



        return i_matched



    def shift_atoms(self, vector_red = None, vector_cart = None, return2cell = 1):
        """
        Shift all atoms according to *vector_red* or *vector_cart*
        Use *return2cell* if atoms coordinates should be inside cell
        """
        st = copy.deepcopy(self)
        if vector_cart is not None:
            vec_cart = np.array(vector_cart)
            for xc in st.xcart:
                xc+=vec_cart
            st.update_xred()
            
        elif vector_red is not None:
            vec = np.array(vector_red)
            for xr in st.xred:
                xr+=vec
            st.xred2xcart()

        
        if return2cell:
            st = st.return_atoms_to_cell()
        return st


    def shake_atoms(self, amplitude = 0.1, el_list = None, ):
        """
        Randomly shake atoms around the lattice 
        amplitude (float) - maximum shift in A it is multiplied by random vector
        el_list (list of int) - shake only el atoms, None - shake all atoms
        """
        st = self.copy()

        ru = random.uniform
        for i, el in enumerate(st.get_elements()):
            # print(el, el_list)
            if el_list is None or el in el_list:
                rand_vec = amplitude*np.array([ru(-1, 1),ru(-1, 1),ru(-1, 1)])
                # print(rand_vec) 

                st.xcart[i]+=rand_vec
        st.update_xred()

        return st

    def replic(self, *args, **kwargs):

        return replic(self, *args, **kwargs)

    def image_distance(self, *args, **kwargs):

        return image_distance(*args, **kwargs)


    def distance(self, i1=None, i2=None, x1=None, x2 = None, coord_type = 'xcart'):
        """
        Shortest distance between two atoms acounting PBC, from 0
        i1 and i2 override x1 and x2

        coord_type - only when x1 and x2 are provided

        """
        # print(self.xcart)
        if i1:
            x1 = self.xcart[i1]
        if i2:
            x2 = self.xcart[i2]
        return image_distance(x1, x2, self.rprimd, coord_type = coord_type)[0]

    def remove_close_lying(self, rm_both = 0, rm_first = 0, tol = 0.4):
        """
        rm_both (bool) - if True than remove both atoms, if False than only the second atom is removed
        rm_first (bool) - if True than the first atom of two overlapping is removed, otherwise the second atom is removed
        tol (float) - atoms separated by distance less than *tol* A are removed

        PBC is realized through image_distance

        SIDE:
            write _removed field to returned st

        TODO:
            Works incorrectly for more than one overlap !!!


        """
        st = copy.deepcopy(self)    
        removed = False
        x1_del = []
        x2_del = []
        numbers = range(st.natom)
        count = 0
        for i, x1 in enumerate(st.xcart):
            for j, x2 in list(zip(numbers, st.xcart))[i+1:]:
                # print(i,j)
                # if all(x1 == x2):
                #     continue
                # if all(x1 == x1_del) or (x2 == x2_del):
                #     continue
                if self.image_distance(x1, x2, r = self.rprimd)[0] < tol:
                    count+=1
                    if count > 1:
                        raise RuntimeError # for detecting multiple overlaps please make this function more universal - removing not by numbers, but by coordinates or more intelligent- see remove_atoms()
                    x1_del = x1
                    x2_del = x2
                    if rm_both:
                        printlog('remove_close_lying(): Atoms', i,j, 'of types ', st.get_elements()[i], st.get_elements()[j], 'are removed')

                        st = st.remove_atoms([i, j])
                        removed = True
                    else:
                        printlog('remove_close_lying(): Atom', j, 'of type ', st.get_elements()[j], 'is removed')
                        
                        if rm_first:
                            st = st.remove_atoms([i]) # the existing atom is removed
                        else:
                            st = st.remove_atoms([j]) # the added atom is removed
                        
                        removed = True
        st._removed = removed

        return st, x1_del, x2_del


    def remove_close_lying2(self, tol = 0.4):
        """
        Very fast removal, linear with respect to number of atoms
        Support multiple overlaps
        Makes a mesh with cubes and determine wich atoms appeared in the mesh 

        Leaves only the first entry, the rest are removed

        
        TODO:
        A problem may occur if two-close lying atoms goes into neibouring bins
            This may be solved by making three additional meshes, 
            where all atoms are shifted by (tol/2, 0, 0), (0, tol/2, 0), (0, 0, tol/2)

        Problem with PBC, not taken into account
            can be solved by making mesh for atoms with (tol, tol, tol) shift 

        Then all overlaps detected on any mesh are treated 

        Replace all overlapping atoms by their center of gravity.


        """
        st = copy.deepcopy(self)    
        
        st = st.return_atoms_to_cell()
        vl = np.array(st.rprimd_len() )
        t = tol/vl
        # print(meshes)
        # sys.exit()
        # shifts_xc = [[0,0,0], [tol, tol, tol], [tol/2, 0, 0], [0, tol/2, 0], [0, 0, tol/2]]
        # shifts = [[0,0,0], [0.5,0.5,0.5], [t[0]/2, 0, 0], [0, t[1]/2, 0], [0, 0, t[2]/2]]
        # shifts = [[0,0,0], [1.5*t[0], 0, 0], [0, 1.5*t[1], 0], [0, 0, 1.5*t[2]]]
        shifts = []
        space = np.linspace(1.9,2.9,3)
        if 1:
            for s1 in space:
                for s2 in space:
                    for s3 in space:
                        s = np.array([s1,s2,s3])
                        # print(s*t) 
                        shifts.append(s*t)
        
        else:
            for i in 0,1,2:
                v = [0,0,0]
                for s in space:
                    v[i] = s
                    # print(v*t)    
                    shifts.append(v*t)

        # sys.exit()


        # print(t)
        print('Number of shifts is ', len(shifts) )


        meshes = [{} for i in shifts]

        
        NML = (vl/tol).astype(int) # mesh sizes

        # print('Mesh sizes are', NML, )#np.array([1/31, 1/27, 1/36])*NML)
        for i, xr in enumerate(st.xred):
            for m, s in zip(meshes, shifts):
                # print(xr)
                xrs = return_xred(xr+s)
                # print(xrs,'\n')
                p = (xrs*NML).astype(int)
                pos = str(p[0])+' '+str(p[1])+' '+str(p[2])
                # print()
                # print(pos)

                if pos in m:
                    m[pos].append(i)
                else:
                    m[pos] = [i]
            # print('\n')
        rem_lists = []
        for i, m in enumerate(meshes):
            rem_list_m = []
            for key in m:
                if len(m[key]) > 1:
                    # x = sum([st.xcartxr for ]  )
                    rem_list_m.extend(m[key][1:])
                    # print(x)
            # print('For mesh', i, 'the list of atoms to remove is ',rem_list_m )
            rem_lists.append(rem_list_m)

        # print( list(set(rem_lists[0]).symmetric_difference(set(rem_lists[1]))) )
        # print( list(set(rem_lists[2]).symmetric_difference(set(rem_lists[3]))) )
        rem_flat_list = [item for sublist in rem_lists for item in sublist]
        # print(rem_flat_list)
        # print(list(set(rem_flat_list)))
        nat_before = st.natom
        st = st.remove_atoms(rem_flat_list)
        print(nat_before- st.natom, 'atoms were removed')

        return st


    def remove_closest(self, el, nn = 6, n = 2, x = 0.0):
        """
        Remove closest lying atoms of type el  
        st (Structure) - input structure 
        el (int array) - list of elements to remove
        nn (int) - number of closest atoms 
        n (int array) - number of removing atoms 
        x (float array) - relative number of removing atoms 
        author - A. Burov 

        """
        st = copy.deepcopy(self)
        
        atoms = st.get_specific_elements(required_elements = el, fmt = 'n', z_range = None, zr_range = None)
        if (x != 0):
            itr = x
        elif (x_rel != 0):
            itr = int(len(atoms)*x_rel)
        else:
            itr = 0
        atoms_removed = [] 
        for i in range(n):
            dist_min = 1e3 
            idx_min = -1  
            for atom_idx in atoms:
                dist = st.nn(atom_idx, n, from_one = 0, only=[3], silent = 1)['dist'][1:]
                dist_cur = sum(dist)/len(dist)
                if (dist_cur < dist_min):
                    dist_min, idx_min = dist_cur, atom_idx
            st = st.remove_atoms([idx_min], from_one = 0, clear_magmom  = 1)
            atoms_removed.append(idx_min)
            print("Atoms were removed: {}".format(i+1))
        print("Atoms with indicies {} were removed".format(atoms_removed))
        print("The final reduced formula is {}".format(st.get_reduced_formula()))
        return st



    def find_closest_atom(self, xc = None, xr = None):
        """
        Find closest atom in structure to xc (cartesian) or xr (reduced) coordinate

        RETURN:
        i shifts, and dist
        """
        if xc is not None:
            x = np.asarray(xc)
            coord_type = 'xcart'
            coords = self.xcart
        if xr is not None:
            x = np.asarray(xr)
            coord_type = 'xred'
            coords = self.xred

        # abs_shifts = [np.linalg.norm(x-x1) for x1 in self.xcart]
        
        abs_shifts = [self.distance(x1 = x, x2 = x1, coord_type = coord_type) for x1 in coords]
        # print(sorted(abs_shifts))
        i = np.argmin(abs_shifts)
        return i, abs_shifts[i], self.distance(x1 = x, x2 = coords[i], coord_type = coord_type)

    def find_closest_neighbor(self,i_at):
        #find closest atom in structure to i_at

        x = self.xcart[i_at]
        abs_shifts = []
        for x1 in self.xcart:
            if list(x1) != list(x):
                abs_shifts.append(np.linalg.norm(x-x1))
        i = np.argmin(abs_shifts)
        return i, abs_shifts[i], x - self.xcart[i]

    def nn(self, i, n = 6, ndict = None, only = None, silent = 0, 
        from_one = True, more_info = 0, oxi_state = 0, print_average = 0):
        """
        show neigbours

        INPUT:
        i - number of central atom, from 1 or 0 (from_one = True or False)
        n - number of neigbours to return
        ndict (dic) - number of specific neigbour atoms to take into account e.g ndict = {8:3} - 3 oxygen atoms will be considered
        only - list of interesting z neighbours

        more_info - return more output - takes time

        from_one - if True, strart first atom from 1, otherwise from 0

        oxi_state (bool) - if 1 then showing oxidation state as well

        print_average (bool) - print more

        RETURN
            dict with the following keys:
            'av(A-O,F)'
            'numbers'
            'dist'
            'xcart'
            'st' - surrounding


        Important:
            'numbers' from 0 in the new version!!!!!


        """

        if from_one:
            i -= 1
            mod = 1
        else:
            mod = 0 # for table 
        st = self
        zn = st.znucl
        x = st.xcart[i]
        out_or = local_surrounding(x, st, n, 'atoms', True, only_elements = only)
        # out =  (xcart_local, typat_local, numbers, dlist )
        # print(out_or)
        out = list(out_or)
        # out[0] = list(itertools.chain.from_iterable(out[0]))
        out[1] = [invert(zn[o-1]) for o in out[1]]
        numbers = copy.copy(out[2]) 
        out[2] = [o+mod for o in out[2]]

        out_tab = [range(0, len(out[2])), out[2], out[1], out[3]]

        tab = np.asarray(out_tab).T.tolist()

 
        # df = pd.DataFrame(tab)
        # print(df)
        if  silent:
            imp = ''
        else:
            imp = 'Y'
        printlog('Neighbors around atom', i+mod, st.get_elements()[i],':', imp = imp)
        # if not silent:
        
        headers = ['nn', 'No.', 'El', 'Dist, A']
        if oxi_state:
            headers.append('Oxi state')
            i = 0 
            oxi = st.get_oxi_states()
            for t in tab:
                i_at = numbers[i]
                t.append(oxi[i_at])
                i+=1


        if tabulate:
            printlog( tabulate(tab[1:], headers = headers, tablefmt='psql', floatfmt=".2f"), imp = imp )
        else:
            printlog(tab[1:], imp = imp )

        info = {}
        info['numbers'] = out_or[2]
        info['dist'] = out_or[3]
        info['xcart'] = out_or[0]


        el = st.get_elements()
        info['el'] = [el[i] for i in out_or[2]]
        info['av(A-O,F)'] = local_surrounding(x, st, n, 'av', True, only_elements = [8,9], round_flag = 0)

        if more_info:
            info['avsq(A-O,F)'] = local_surrounding2(x, st, n, 'avsq', True, only_elements = [8,9])
            info['avharm(A-O,F)'] = local_surrounding2(x, st, n, 'avharm', True, only_elements = [8,9])
            info['avdev(A-O,F)'], _   = local_surrounding2(x, st, n, 'av_dev', True, only_elements = [8, 9])
            info['sum(A-O,F)'] = local_surrounding2(x, st, n, 'sum', True, only_elements = [8,9])

        t = set(out_or[2])
        s = set(range(st.natom)) 
        d = s.difference(t) 
        # d = d.remove(i)
        # print(t)
        # print(i)
        # print(d)
        st_left = st.remove_atoms(d)
        st_left.name+='_loc'
        # sys.exit()
        st_left.dlist = out_or[3] # distances to neighbours
        st_left.ellist = info['el'] # types of neighbours
        info['st'] = st_left

        if ndict:
            info['av(A-O)']   = local_surrounding(x, st, ndict[8], 'av', True, only_elements = [8])
            info['avdev(A-O)'], _   = local_surrounding(x, st, ndict[8], 'av_dev', True, only_elements = [8])
            info['min(A-O)'], _ ,info['max(A-O)']    = local_surrounding(x, st, ndict[8], 'mavm', True, only_elements = [8])
            atoms = local_surrounding(x, st, ndict[8], 'atoms', True, only_elements = [8])
            info['Onumbers'] = atoms[2][1:] # exclude first, because itself!
            # print(info['Onumbers'])

        if print_average:
            print('av(A-O,F)', info['av(A-O,F)'])

        return info

    def check_JT(self, criteria = 0.03):
        #Check Yan-Teller effect
        #check average TM-O distance in the cell and find the outstanding bonds 
        #return TM-O dist list
        #criteria - value in % of average bond length when bond is outstanding
        #

        tra = self.get_transition_elements()
        
        if len(tra): 
            print('Starting...\n\n I ve obtained  %i TM atoms \n\n\n'%len(tra))
        else:
            print('Starting...\n\n I ve obtained  no TM atoms \n\n\n')
            return 
        


        el = self.get_elements()

        aver_list = []
        dist_list = []


        # print(self.nn(1, silent = 1)['dist'][1:],self.nn(1, silent = 1)['numbers'][1:])

        for i in range(0, len(el)):
            d = []
            if el[i] in tra:
                dist = self.nn(i+1, silent = 1)['dist'][1:]
                numbers = self.nn(i+1, silent = 1)['numbers'][1:]
                # print(numbers)
                n = self.nn(i+1, silent = 1)['numbers'][0]
                for k in range(0,len(dist)):
                    if el[numbers[k]] == 'O':
                        d.append(dist[k])
                        dist_list.append([round(dist[k],4),n,numbers[k]])
                aver_list.append(round(np.mean(d),2))
        # print(dist_list)
        # print(aver_list)
        aver_distance = round(np.mean(aver_list),2)
        print('Average TM-O bond length is %s A \n'%aver_distance)
       
        k = 0

        min_dist = []
        max_dist = []


        for i in dist_list:
            if (el[i[1]] == 'O' or el[i[2]] == 'O'):
                if i[0] > aver_distance*(1+criteria): 
                    max_dist.append(i[0])    
                    # print('Outstanding bond length %.4s between %s (%s) and %s (%s) \n'%(i[0],i[1], el[i[1]],i[2], el[i[2]]))
                    k = 1

                if i[0] < aver_distance*(1-criteria):
                    min_dist.append(i[0])    
                    # print('Outstanding bond length %.4s between %s (%s) and %s (%s) \n'%(i[0],i[1], el[i[1]],i[2], el[i[2]]))
                    k = 1

        if k:
            maxd = round(np.mean(max_dist),2)
            mind = round(np.mean(min_dist),2)

            if maxd and mind: 
                print('Jahn-Teller effect is found\n Average min TM-O length is %s \n Average max TM-O length is %s \n'%(mind, maxd) )


        if not k: print('Ok! None outstanding bonds found\n')

        return dist_list


    def find_unique_topologies(self, el1, el2, nn = 6, tol = 0.5, told = 0.005, tolmag = 0.4, write_loc = 0):

        """
        Looks for unique topologies
        Currently only octahedral and pentahedral are realized

        el1, el2 (str) - elements that forms topology
        nn (int) - number of neighbours for topology analysis
        tol (float) - tolerance for unique centers defined by deviation, mA
        told (float) - tolerance for distances applied for grouping bonds, A
        tolmag (float) - tolerance for magnetic moments works with tol
        write_loc (int) - write local topology


        """

        def group_bonds(lengths, tol):
            #
            lengths = list(np.around(lengths, 2))
            unique = []
            unique.append(lengths[0])
            groups = {}
            for l in lengths[1:]:
                if min(np.abs(unique-l)) > tol:
                    unique.append(l)
            # print('lengths', lengths)
            # print('unique bonds are', unique)
            for u in unique:
                groups[u] = 0
                for l in lengths:
                    if abs(l-u) < tol:
                        groups[u] += 1
            return groups


        st = self
        z1 = invert(el1)
        z2 = invert(el2)
        n1 = self.get_specific_elements([z1])


        unique_centers = [] # numbers of unique topology centers 
        unique_deviations = []
        unique_magmoms = []
        av_dev5 = 0
        for i in n1:
            x = st.xcart[i]

            av_dev, _   = local_surrounding2(x, st, nn, 'av_dev', True, only_elements = [z2], round_flag = 0 )
            # if av_dev > 100:
                #probably surface atom, deviation is too much
            mag = st.magmom[i]
            print('Deviation for atom {:d} is {:.1f}'.format(i, av_dev) )
            if len(unique_centers) == 0:
                unique_centers.append(i)
                unique_deviations.append(av_dev)
                unique_magmoms.append(mag)
                continue
            # print(unique_centers)
            # print(av_dev, min(np.abs(np.array(unique_deviations-av_dev))))
            # print(np.array(unique_magmoms-mag)) 
            if min(np.abs(np.array(unique_deviations-av_dev))) < tol and min(np.abs(np.array(unique_magmoms)-mag)) < tolmag:
                continue
            else:
                unique_centers.append(i)
                unique_magmoms.append(mag)
                unique_deviations.append(av_dev)
        
        # pretty = pprint.PrettyPrinter(width=30)

        print('Unique centers are ', unique_centers,'. number, deviation6, deviation5, magmom, and topology of polyhedra and  for each:')
        for i, d in zip(unique_centers, unique_deviations):
            dic = st.nn(i, only = [z2], from_one = 0, silent = 1)
            lengths = dic['dist'][1:]
            av = dic['av(A-O,F)']
            if d > 100:
                x = st.xcart[i]
                av_dev5, _   = local_surrounding2(x, st, 5, 'av_dev', True, only_elements = [z2], round_flag = 0 )
                st.name+=str(i)
                if write_loc:
                    st.write_xyz(show_around=i+1, analysis = 'imp_surrounding', only_elements = [z2])
            # print(lengths)
            groups = group_bonds(lengths, told)
            print( '{:2d} | {:4.1f} | {:4.1f} | {:4.1f} :'.format(i, d, av_dev5, st.magmom[i]))
            print(groups, 'av={:.2f} \n'.format(av))

        return

    def center(self, reduced = 0):
        #return cartesian or reduced center of the cell
        if reduced:
            center = np.sum(self.xred, 0)/self.natom
        else:
            center = np.sum(self.xcart, 0)/self.natom

        return center

    def center_on(self, i):
        #calc vector which alows to make particular atom in the center 
        if i and i < len(self.xred):
            x_r = self.xred[i]
            center = np.sum(self.xred, 0)/self.natom
            # print(center)
            # sys.exit()
            # print(x_r)
            dv = center - x_r
            # print(dv)
            # print(dv+x_r)
        else:
            dv = None
        return dv







    def localize_polaron(self, i, d, nn = 6):
        """
        Localize small polaron at transition metal by adjusting TM-O distances
        i - number of transition atom, from 0
        d - shift in angstrom; positive increase TM-O, negative reduce TM-O distance
        nn - number of neigbours

        """
        # nn

        st = copy.deepcopy(self)
        TM = st.get_el_z(i)
        TM_name = st.get_el_name(i)
        if TM not in header.TRANSITION_ELEMENTS:
            printlog('Warning! provided element ', TM_name, 'is not a transition metal, I hope you know what you are doing. ')

        silent = 1
        if 'n' in header.warnings or 'e' in header.warnings:
            silent = 0
        # silent = 0


        dic = st.nn(i, nn, from_one = 0, silent = silent, only = [8,9])
        printlog('Average TM-O distance before localization is {:.2f}'.format(dic['av(A-O,F)']), imp = '')

        #updated xcart
        xc = st.xcart[i]
        for j, x in zip(dic['numbers'][1:], dic['xcart'][1:]):
            x1 = st.xcart[j]
            
            # print(xc, x)
            v = xc-x
            # print(v)
            vn = np.linalg.norm(v)
            mul = d/vn
            # print(vn, mul)
            dv = v * mul
            # print(st.xcart[j])
            st.xcart[j] = st.xcart[j] -  dv 
            # print(st.xcart[j])

        st.update_xred()

        dic = st.nn(i, nn, from_one = 0, silent = silent, only = [8,9])
        printlog('Average TM-O distance after localization is {:.2f}'.format(dic['av(A-O,F)']), imp = '')

        st.name+='pol'+str(i+1)

        return st


    def localize_polaron_dist(self, i_center, d, nn = 6, axis = None, direction = None, mode = 'axis_expand'):
        """
        
        localization small polaron at transition metal by adjusting TM-O distances
        with distortions
        i_center - number of transition atom, from 0
        d - shift in angstrom; positive increase TM-O, negative reduce TM-O
            or shift along *direction* 
        nn - number of neighbors
        axis - axis of octahedra, 0, 1, 2,
        direction - vector to shift central atom in reduced coordinates

            Axes of octahedra are determined relative to cartesian coordinates
        mode 
            'axis_expand'
            'shift_center'

        TODO
        Make it more general to include any ligands; now only O and F are supported
        Make for other topologies, now tested only for octa

        """
        st = copy.deepcopy(self)
        TM = st.get_el_z(i_center)
        TM_name = st.get_el_name(i_center)
        if TM not in header.TRANSITION_ELEMENTS:
            printlog('Warning! provided element ', TM_name, 'is not a transition metal, I hope you know what you are doing. ')

        silent = 1
        if 'n' in header.warnings or 'e' in header.warnings:
            silent = 0
        # silent = 0

        np.set_printoptions(formatter={'float': '{: 6.2f}'.format})

        dic = st.nn(i_center, nn, from_one = 0, silent = silent)
        av = dic['av(A-O,F)']
        printlog('Average TM-O distance before localization is {:.2f}'.format(av), imp = '')

        xc = st.xcart[i_center]
        ligand_xcart = [x-xc for x in  dic['xcart'][1:]]
        # print(np.array(ligand_xcart))
        ligand_order     = copy.copy(list(dic['numbers'][1:]))



        # find octahedron axes
        pairs = []# first, second, and third pair are along first, second, and third octahedron axes
        order = []
        # if id(x1) in map(id, checked):
        i=0
        for i1, x1 in zip(ligand_order, ligand_xcart):
            for i2, x2 in zip(ligand_order[i+1:], ligand_xcart[i+1:]):
                # print( x1+x2 )
                ssum = np.linalg.norm(x1+x2 ) # for ideal octa should be zero for axis
                if ssum < av/2: #should work even for highly distorted octahedra
                    pairs.append(x1)
                    pairs.append(x2)
                    order.append(i1)
                    order.append(i2)
                    # print(av, ssum)
            i+=1
        if len(pairs) < len(ligand_xcart):
            #only two axes detected; i.e. pyramid; the third axis is determined relative to the center
            for i1, x1 in zip(ligand_order, ligand_xcart):
                # if x1 in pairs:
                if id(x1) in map(id, pairs):
                    continue
                else:
                    pairs.append(x1)
                    pairs.append(xc-xc)
                    order.append(i1)
                    order.append(dic['numbers'][0])

        # print(np.array(pairs))

        # check the order of pairs; to have first vector in positive xy quater
        # and third vector #determine pair along z

        pairs_new = [0]* len(pairs)
        order_new = [0]* len(order)
        # print(xc)
        # print(np.array(dic['xcart'][1:]))
        # print(np.array(pairs))
        # print(order)
        iz = np.argmax(np.abs(np.array(pairs).dot([0,0,1]) ))//2
        pairs_new[4] = copy.copy(pairs[iz*2])
        pairs_new[5] = copy.copy(pairs[iz*2+1])
        order_new[4] = order[iz*2]
        order_new[5] = order[iz*2+1]
        del pairs[iz*2+1] # 
        del pairs[iz*2]
        del order[iz*2+1] # 
        del order[iz*2]        
        ix = np.argmax(np.array(pairs).dot([1,0,0]) )//2
        iy = np.argmax(np.array(pairs).dot([0,1,0]) )//2
        
        q1 = sum(np.sign(pairs[0][0:2]))-sum(np.sign(pairs[1][0:2])) # for positive quater should be 4, for negative zero
        q2 = sum(np.sign(pairs[2][0:2]))-sum(np.sign(pairs[3][0:2]))
        # print(q1, q2)

        if q1 > q2:
            pairs_new[0:4] = order
            order_new[0:4] = order
        else:
            #swap axis 
            pairs_new[0:2] = pairs[2:4]
            pairs_new[2:4] = pairs[0:2]
            order_new[0:2] = order[2:4]
            order_new[2:4] = order[0:2]

        # print(order)
        # print(order_new)


        if mode == 'shift_center':
            #shift along lattice vectors
            #a12 
            ''
            v = np.dot( direction, st.rprimd) # cart
            vn = np.linalg.norm(v)
            dv = v/vn*d
            # print(dv)
            st.xcart[i_center] = st.xcart[i_center] + dv

        if mode == 'axis_expand':
            #expand or shring alond one of the axes by d
            # axis = 0
            # print(order_new[axis*2:axis*2+2])
            for i in order_new[axis*2:axis*2+2]:
                x = st.xcart[i]
                print(x)
                v = xc-x
                vn = np.linalg.norm(v)
                if vn < 0.1:
                    mul = 0 # central atom will not move! exatly what we need
                else:
                    mul = d/vn
                dv = v * mul
                st.xcart[i] = st.xcart[i] -  dv 
                print(st.xcart[i] )

        # sys.exit()

        st.update_xred()

        dic = st.nn(i, 6, from_one = 0, silent = silent)
        printlog('Average TM-O distance after localization is {:.2f}'.format(dic['av(A-O,F)']), imp = '')

        st.name+='pol'+str(i+1)

        return st


    def make_polarons(self, atoms, pol_type = 'hole', mag = None, silent = 1):
        """
        create polarons
        """
        st = self.copy()
        for i in atoms:
            st = st.localize_polaron(i, d=-0.1, nn = 6)
            st.magmom[i] = mag


        return st

    def ewald(self, ox_st = None, site = None):
        # ox_st 
        #   # 1 - oxidation states from guess
            # 2 - from potential
            # None - from charge
        # site if provided (from 0), than site energy is printed
        from pymatgen.analysis.ewald import EwaldSummation
        # from siman.analysis import set_oxidation_states
        st = self
        if ox_st == 1:
            # st = set_oxidation_states(st)
            # st.printme()
            stpm = st.convert2pymatgen(chg_type = 'pm')
            # print('The following oxi states were set', st.oxi_state)
        if ox_st == 2:
            stpm = st.convert2pymatgen(chg_type = 'pot')

        else:
            stpm = st.convert2pymatgen()
        
        ew = EwaldSummation(stpm)
        if site is not None:
            site_e = 2*ew.get_site_energy(site)
            print('Energy for site ', st.get_elements()[site], site_e)

            return ew.total_energy,  site_e
        else:
            return ew.total_energy


    def write_espresso(self, filename = None, shift = None):
        st = copy.deepcopy(self)
        st = st.remove_atoms(['void']) # remove voids
        if shift:
            st = st.shift_atoms(shift)
        if not filename:
            filename = ('xyz/espresso_'+st.name).replace('.', '_')

        makedir(filename)

        printlog('Starting writing Quantum Espresso', filename)

        with io.open(filename,'w', newline = '') as f:
            f.write('ATOMIC_POSITIONS\n')
            for el, x in zip(st.get_elements(), st.xred):
                f.write(" {:2s}   {:12.10f}  {:12.10f}  {:12.10f} \n".format(el, x[0], x[1], x[2]) )




        return


    def write_poscar(self, filename = None, coord_type = 'dir', vasp5 = True, charges = False, energy = None, selective_dynamics = False, shift = None):
        """
        write 

        charges (bool) - write charges, self.charges should be available
        energy - write total energy

        selective dynamics - 
            if at least one F is found than automatically switched on
            !works only for coord_type = 'dir' and charges = False
            None - not written

        shift - shift atoms
        
        NOTE
        #void element type is not written to POSCAR

        TODO
            selective_dynamics for coord_type = 'cart'
            Velocity and predictor are not reordered; Can be used only for continiue MD
        """




        st = copy.deepcopy(self)
        st = st.remove_atoms(['void']) # remove voids
        if shift:
            st = st.shift_atoms(shift)


        to_ang = 1
        rprimd = st.rprimd
        xred   = st.xred
        xcart  = st.xcart
        typat = st.typat  
        znucl = st.znucl
        els   = st.get_elements()


        # print(st.convert2pymatgen())

        # print()
        try:
            select = st.select
        except:
            st = st.selective_all()
            select = st.select

        # print(select)

        if selective_dynamics is False:
            selective_dynamics = st.check_selective()

        # print(selective_dynamics)


        if not filename:
            filename = os.getcwd()+('/xyz/POSCAR_'+st.name).replace('.', '_')

        makedir(filename)

        printlog('Starting writing POSCAR', filename, 'Vasp5:', vasp5)

        # print 
        """1. Generate correct nznucl and zxred and zxcart"""
        zxred  = [[] for i in znucl]
        zxcart = [[] for i in znucl]
        zchar = [[] for i in znucl]
        zelem = [[] for i in znucl]
        zselect = [[] for i in znucl]
        zmagmom = [[] for i in znucl] # not used for the moment
        ziatom = [[] for i in znucl]
        iatom = range(st.natom)
        # nznucl = []
        if len(typat) != len(xred) or len(xred) != len(xcart):
            raise RuntimeError
        

        # print(xred)
        # print(typat)

        for t, xr, xc, el, i in zip(typat, xred, xcart, els, iatom ):
            # print ("t ", t-1, xr)
            zxred[ t-1].append(xr)
            zxcart[t-1].append(xc)
            zelem[t-1].append(el)
            ziatom[t-1].append(i)
        
        if selective_dynamics:
            for t, s in zip(typat, select):
                zselect[t-1].append(s)

        # print(zselect)
        # print(zxred)


        # print(charges)
        if charges:
            charg = self.charges

            for t, ch in zip(typat, charg):
                zchar[t-1].append(ch)

        if st.magmom:
            '' # if needed put it here

        poscar_atom_order = []
        for iatom in ziatom:
            for i in iatom:
                poscar_atom_order.append(i)


        self.poscar_atom_order = poscar_atom_order

        nznucl = [len(xred) for xred in zxred]
        # print(nznucl)
        # sys.exit()

        elnames = [element_name_inv(z) for z in znucl]



        with io.open(filename,'w', newline = '') as f:
            """Writes structure (POSCAR) in VASP format """
            if energy:
                energy_string = 'e0='+str(energy)+' ; '
            else:
                energy_string = ''

            f.write('i2a=['+list2string(elnames).replace(' ', ',') + '] ; '+energy_string)
            
            if hasattr(self, 'tmap'):
                f.write('tmap=[{:s}] ; '.format(list2string(st.tmap).replace(' ', ',') ))

            # print(self.name)
            f.write(self.name)


            f.write("\n{:18.15f}\n".format(1.0))
            
            for i in 0, 1, 2:
                f.write('{:10.6f} {:10.6f} {:10.6f}'.format(rprimd[i][0]*to_ang,rprimd[i][1]*to_ang,rprimd[i][2]*to_ang) )
                f.write("\n")

            if vasp5:
                for el in elnames:
                    f.write(el+' ')
                f.write('\n')


            for n in nznucl:    
                f.write(str(n)+' ')
                # print(str(n)+' ')
            f.write('\n')

            if selective_dynamics:
                f.write("Selective dynamics\n")



            if "car" in coord_type:
                print_and_log("Warning! Cartesian regime of coordination may be obsolete and incorrect !!!", imp = 'Y')
                f.write("Cartesian\n")
                for xcart in zxcart:
                    for x in xcart:
                        f.write(str(x[0]*to_ang)+" "+str(x[1]*to_ang)+" "+str(x[2]*to_ang))
                        f.write("\n")

                

            elif "dir" in coord_type:
                f.write("Direct\n")
                
                if charges:
                    for xred, elem, char in zip(zxred, zelem, zchar, ):
                        for x, el, ch in zip(xred, elem, char):
                            f.write("  {:12.10f}  {:12.10f}  {:12.10f}  {:2s}  {:6.3f}\n".format(x[0], x[1], x[2], el, ch) )
                elif selective_dynamics:
                    for xred, select in zip(zxred, zselect):
                        for x, s in zip(xred, select):
                            # print(x,s)
                            f.write("  {:19.16f}  {:19.16f}  {:19.16f}  {:s} {:s} {:s}\n".format(x[0], x[1], x[2], b2s(s[0]), b2s(s[1]), b2s(s[2])) )
                else:
                    for xred  in zxred :
                        for x in xred :
                            f.write("  {:19.16f}  {:19.16f}  {:19.16f}\n".format(x[0], x[1], x[2]))

            


            elif 'None' in coord_type:
                pass

            else:
                print_and_log("Error! The type of coordinates should be 'car' or 'dir' ")
                raise NameError



            # print('write_poscar(): predictor:\n', st.predictor)
            if hasattr(st, 'vel') and len(st.vel)>0:
                printlog("Writing velocity to POSCAR ", imp = 'y')
                # f.write("Cartesian\n")
                f.write("\n")
                for v in st.vel:
                    f.write( '  {:18.16f}  {:18.16f}  {:18.16f}\n'.format(v[0]*to_ang, v[1]*to_ang, v[2]*to_ang) )

            if hasattr(st, 'predictor') and st.predictor:
                printlog("Writing predictor POSCAR ", imp = 'y')
                f.write("\n")
                f.write(st.predictor)


        

        f.close()
        # if os.getcwd() not in filename:
        #     print('rep', str(os.getcwd()), filename)
        #     path = os.getcwd()+'/'+filename
        # else:
        path = filename
        print_and_log("POSCAR was written to", path, imp = 'y')
        return path



    def write_cif(self, filename = None, mcif = False, symprec = 0.1, write_prim = 0):
        """
        Find primitive cell and write it in cif format
        
        filename (str) - name of produced cif file
        mcif (bool) - if True, than write mcif file with magnetic moments included, primitive cell is not supported
        symprec (float) - symmetry precision, symprec = None allows to write the structure as is
        write_prim (bool) - convert structure to primitive 
        

        """
        
        if mcif:
            m = 'm'
        else:
            m = ''

        if filename == None:
            filename = os.getcwd()+'/cif/'+self.name


        makedir(filename)

        st_mp = self.convert2pymatgen()

        # print(st_mp)

        try:
            sg_before =  st_mp.get_space_group_info() 


            sf = SpacegroupAnalyzer(st_mp, symprec = symprec)

            st_mp_prim = sf.find_primitive()

            sg_after = st_mp_prim.get_space_group_info()

        except:
            sg_before = [None]
            sg_after = [None]
            st_mp_prim = None
            printlog('Warning! could not analyze space group')

        if sg_before[0] != sg_after[0]:
            printlog('Attention! the space group was changed after primitive cell searching', sg_before, sg_after)
            printlog('I will save supercell in cif Pay attention that CifWriter can symmetrize and change vectors. Also use *write_prim*')
            # st_mp_prim = st_mp
            # symprec = 0.001
            # symprec = None

        if mcif:
            cif = CifWriter(st_mp, symprec = symprec, write_magmoms=mcif)
        else:
            if st_mp_prim:
                cif_prim = CifWriter(st_mp_prim, symprec = symprec, )
            
            cif = CifWriter(st_mp, symprec = symprec, )

        
        cif_name =  filename+'.'+m+'cif'
        cif_prim_name =  filename+'_prim.'+m+'cif'
        

        cif.write_file( cif_name  )
        if write_prim:
            cif_prim.write_file( cif_prim_name  )
        
        printlog('Writing cif', cif_name, imp = 'y')

        return cif_name

    def write_xyz(self, *args, **kwargs):
        #see description for write_xyz()
        return write_xyz(self, *args, **kwargs)



    def write_lammps(self, *args, **kwargs):
        return write_lammps(self, *args, **kwargs)


    def read_xyz(self, *args, **kwargs):
        
        # print(self.perm)
        return read_xyz(self, *args, **kwargs)




    def jmol(self, shift = None, r = 0, show_voids = False, rep = None, program = 'jmol'):
        """open structure in Jmol or vesta
        
        INPUT:
        shift (list) - shift vector  in reduced coordinates
        r (int ) - parameter
            0 - open POSCAR
            1 - open OUTCAR to see optimization steps
            2 - open mcif to see magnetic moments
            3 - xyz
        show_voids (bool) - replace voids (z = 300) with Po to visualize them
        rep  (list 3*int) - replicate along vectors
        program - 
            'jmol'
            'vesta'
        
        """
        st = copy.deepcopy(self)

        if rep:
            st = st.replic(rep)

        if shift:
            st = st.shift_atoms(shift)

        if show_voids:
            atom_numbers = st.get_specific_elements([300])
            st = st.replace_atoms(atom_numbers, 'Po')

        # filename, _ = st.write_xyz()
        if r == 1:
            filename = st.outfile 
        elif r == 2:
            filename = st.write_cif(mcif = 1)
        elif r == 3:
            filename = st.write_xyz()[0]
        else:
            filename = st.write_poscar(vasp5 = 1)
        
        # print(r, filename)
        # sys.exit()
        if 'jmol' in program :
            runBash(header.PATH2JMOL+' '+filename, detached = True)
        elif 'vesta' in program:
            runBash(header.PATH2VESTA+' '+filename, detached = True)

        return

    def vesta(self, *args, **kwargs):
        kwargs['program'] = 'vesta'
        self.jmol(*args, **kwargs)



    @property
    def vlen(self):
        #return vector lengths
        r = self.rprimd
        n = np.linalg.norm
        return n(r[0]), n(r[1]), n(r[2])


class Calculation(object):
    """Main class of siman. Objects of this class contain all information about first-principles calculation
        List of important fields:
            - init (Structure)
            - end  (Structure)
            - occ_matrices (dict) - occupation matrices, number of atom (starting from 0) is used as key


    """
    def __init__(self, inset = None, iid = None, output = None):
        #super(CalculationAbinit, self).__init__()
        self.name = "noname"
        if inset:
            self.set = copy.deepcopy(inset)
        else:
            self.set = InputSet()
        
        # if self.set.set_sequence:



        self.init = Structure()
        self.end = Structure()
        self.children = [] # inherited calculations 
        self.state = "0.Initialized"
        self.path = {
        "input":None,
        "input_geo":None,
        "potential":None,
        "output":output}
        self.calc_method = None #
        self.prev = [] # list of previous calculations
        if iid:
            self.id = iid
            self.name = str(iid[0])+'.'+str(iid[1])+'.'+str(iid[2])
        else:
            self.id = (output,'0', 1)
        header.db[self.id] = self
        self.cluster_address = ''
        self.project_path_cluster = ''
    def get_path(self,):
        path = os.path.dirname(os.getcwd()+'/'+self.path['output'])
        print( path)
        return path

    def read_geometry(self, filename = None):
        """Reads geometrical data from filename file in abinit format"""
        if self.path["input_geo"] == None:
            self.path["input_geo"] = filename



        with open(filename,"r") as file:
            #For large files can be time consuming
            memfile = file.read()
            gen_words = memfile.split()

            self.des = '';  
            for line in memfile.splitlines():
                if 'des' in line: 
                    # print line; 
                    self.des = line.split('des ')[1]+';'
                
                self.build = empty_struct()                
                if 'BEGIN BUILD INFORMATION' in line:
                    print_and_log("File contain build information! Start to read", imp = 'n')
                    # self.build = Structure()
                    # # self.build.rprimd = None
                    # # self.build.xred = None
                    # # self.build.xcart = None
                    # # self.build.des = None
                    # # self.build.name = None

                    self.build.calctype = read_list("calctype", 1, str, gen_words)[0]
                    self.build.a_c_conv = read_list("a_c_conv", 4,float, gen_words)
                    self.build.build_natom = read_list("build_natom", 1, int, gen_words)[0]
                    self.build.build_acell = read_vectors("build_acell", 1, gen_words)
                    self.build.build_rprim = read_vectors("build_rprim", 3, gen_words)
                    self.build.build_xred = read_vectors("build_xred", self.build.build_natom, gen_words)
                    self.build.build_ntypat = read_list("build_ntypat", 1, int, gen_words)[0]
                    self.build.build_typat = read_list("build_typat", self.build.build_natom, int, gen_words)
                    self.build.build_znucl = read_list("build_znucl", self.build.build_ntypat, float, gen_words)
                    self.build.hkl1 = read_vectors("hkl1", 1, gen_words)
                    self.build.uvw1 = read_vectors("uvw1", 1, gen_words)
                    self.build.uvw2 = read_vectors("uvw2", 1, gen_words)
                    self.build.uvw3 = read_vectors("uvw3", 1, gen_words)
                    self.build.mul  = read_vectors("mul", 1, gen_words)
                    self.build.nadded = read_list("nadded", 1, int, gen_words)[0] #total number of added atoms after building structure
                    self.build.listadded = read_list("listadded", self.build.nadded, int, gen_words) #list of added atoms corresponding to xred 

                    print_and_log("Build information has been read")




            self.init = Structure()

            #sys.exit()
            self.useable = 0
            #Read total number of atoms
           # nznucl, since We will have more impurity atoms of different types
           #  command="""grep -w -m 1 "natom " """+filename
#             s1=runBash(command)
            self.natom = read_list("natom", 1, int, gen_words)[0]
            # print command
            # print s1
#             self.natom_str = s1
#             if s1=='':
#                 self.natom = 0
#                 print_and_log( """Warning! In filename """+filename+""" not found natom! set to zero.
#                 It is very likely that other parameters was not 
#                 found too, Calculation completely unusable!!!""")
#                 raise RuntimeError
#             else:
#                 self.natom=int(s1.split()[1]) 

            self.acell = read_list("acell", 3, float, gen_words)
            self.rprim = read_vectors("rprim", 3, gen_words)

            self.rprimd = copy.deepcopy( self.rprim )
            for i in 0,1,2:
                self.rprimd[i] = self.rprim[i] * self.acell[i]         #Calculate rprimd
            
            self.vol = np.dot( self.rprimd[0], np.cross(self.rprimd[1], self.rprimd[2])  ); #volume
                      
            self.recip = calc_recip_vectors(self.rprimd) #Determine reciprocal vectors


            self.ntypat = read_list("ntypat", 1, int, gen_words)[0]
            self.typat = read_list("typat", self.natom, int, gen_words)
            if 0 in self.typat:
                print_and_log('Error; 0 in typat is not allowed')
                raise RuntimeError

            self.nznucl = []

            for typ in range(1,self.ntypat+1):
                self.nznucl.append(  self.typat.count(typ) )


            self.znucl = read_list("znucl", self.ntypat, float, gen_words)
            self.xcart = read_vectors("xcart", self.natom, gen_words)
            
            self.xred = read_vectors("xred", self.natom, gen_words)
            #print self.xred
            # print(self.xcart)
            if self.xred is [None]:
                print_and_log("Convert xcart to xred")
                self.xred = xcart2xred(self.xcart, self.rprimd)
            
            if self.xcart is [None]:
                print_and_log("Convert xred to xcart")
                self.xcart = xred2xcart(self.xred, self.rprimd)

            self.hex_a = read_list("hex_a", 1, float, gen_words)[0]
            self.hex_c = read_list("hex_c", 1, float, gen_words)[0]
            self.len_units = read_list("len_units", 1, str, gen_words)[0]

            self.version = read_list("version", 1, int, gen_words)[0]

            self.gbpos = read_list("gbpos", 1, float, gen_words)[0]



            self.init.hex_a = self.hex_a
            self.init.hex_c = self.hex_c
            self.init.gbpos = self.gbpos
            self.init.name = self.name+'.init'
            self.init.xcart = self.xcart 
            self.init.xred = self.xred
            self.init.rprimd = self.rprimd
            self.init.recip = self.recip
            self.init.vol = self.vol
            self.init.znucl = self.znucl
            self.init.nznucl = self.nznucl 
            self.init.typat = self.typat
            self.init.ntypat = self.ntypat 
            self.init.natom = self.natom 



            vel = read_vectors("vel", self.natom, gen_words)
            if vel[0] is not None: 
                self.init.vel = vel

            #read magnetic states; name of vasp variable
            curset = self.set
#             if hasattr(curset, 'magnetic_moments') and curset.magnetic_moments and ('ISPIN' in curset.vasp_params.keys()) and curset.vasp_params['ISPIN'] == 2:
            self.init.magmom = read_list("magmom", self.natom, float, gen_words)
            # if magmom[0] is not None:
            #     self.init.magmom = magmom


            select = read_vectors("select", self.natom, gen_words, type_func = lambda a : int(a), lists = True )
            if None not in select: 
                self.init.select = select

            predictor_length = read_list("pred_length", 1, int, gen_words)[0]
            # print('pred_length', predictor_length)
            # sys.exit()
            if predictor_length:
                predictor = read_string('predictor', predictor_length, memfile)
                # print('predictor', predictor)
                self.init.predictor = predictor




            self.state = "1.Geometry has been read"



        #file.close();

        print_and_log( "If no warnings, geometry has been succesfully read from file "+filename+" \n")

        return





    def write_geometry(self, geotype = "init", description = "", override = False, atomic_units = 0):
        """Writes geometrical data in custom siman format bases on abinit format to self.path["input_geo"]"""
        geo_dic = {}
        geofile = self.path["input_geo"]
        geo_exists = os.path.exists(geofile)
        # print (os.path.exists(geofile))
        
        if atomic_units:
            en = 1/header.to_eV
            le = 1/header.to_ang
        else:
            en = 1
            le = 1            



        if geo_exists:
            if override:
                print_and_log("Warning! File "+geofile+" was replaced"); 
            else: 
                print_and_log("Error! File "+geofile+" exists. To replace it set parameter override"); 
                return False
                #raise RuntimeError
        # print "geofile name, classes:",  geofile
        # print "folder :",  os.path.dirname(geofile)
        if not os.path.exists(os.path.dirname(geofile)):
            os.makedirs(os.path.dirname(geofile))

        if geotype == "init": #write initial structure
            st = self.init
        elif geotype == "end":
            st = self.end
            # if not hasattr(st, 'natom'):  st.natom = self.init.natom
            # if not hasattr(st, 'ntypat'): st.ntypat = self.init.ntypat
            # if not hasattr(st, 'typat'): st.typat = self.init.typat
            # if not hasattr(st, 'znucl'): st.znucl = self.init.znucl 
        else: print_and_log("Error! Unknown geotype \n"); raise RuntimeError                                  

        if st.natom != len(st.xred) != len(st.xcart) != len(st.typat) or len(st.znucl) != max(st.typat): 
            print_and_log("Error! write_geometry: check your arrays.", imp = 'Y')
            raise RuntimeError

        # print (st.magmom)
        # sys.exit()
        with open(self.path["input_geo"],"w", newline = '') as f:
            f.write("des "+description+"\n")
            f.write("len_units "+self.len_units+"\n")
            
            if hasattr(st, 'hex_a'):
                f.write("hex_a "+str(st.hex_a)+"\n")
            if hasattr(st, 'hex_c'):
                f.write("hex_c "+str(st.hex_c)+"\n")

            # try: self.gbpos
            # except AttributeError:
            #     self.gbpos = None
            if hasattr(st, 'gbpos'):
                f.write("gbpos "+str(st.gbpos)+"\n")
            
            if hasattr(self, 'version'):

                f.write("version "+str(self.version)+"\n")
            
            try: 
                st.magmom
            except AttributeError:
                st.magmom = [None]
            # print st.magmom 
            # sys.exit()
            if len(st.typat) != len(st.xred) or len(st.xred) != len(st.xcart):
                printlog('Error! Check sizes of your atom lists')


            if len(st.magmom) > 0 and not None in st.magmom:
                mag_str = ' '.join(np.array(st.magmom).astype(str))
                f.write("magmom "+'\n'.join( wrap(mag_str) ) +"\n")
                if len(st.typat) != len(st.magmom):
                    printlog('Error! Check size of your magmom list')


            f.write("acell 1 1 1\n")

            f.write("natom  " +str(st.natom) +"\n")
            
            f.write("ntypat " +str(st.ntypat) +"\n")
            
            f.write("znucl  ")
            for z in st.znucl:
                f.write(str(z)+" ")
            
            f.write("\ntypat  ")
            i = 0
            for t in st.typat:
                f.write("%i "%(t)  ); i+=1;
                if i >= 20: f.write("\n"); i = 0;

            f.write("\nrprim  ")
            for v in st.rprimd:
                f.write("%.12f %.12f %.12f \n"%(v[0]*le, v[1]*le, v[2]*le)  )

            f.write("xred  ")
            #print st.xred
            if len(st.xred) != st.natom: print_and_log("Warning! write_geometry(): xred is empty or overfull\n");raise RuntimeError 
            for v in st.xred:
                f.write("%.12f %.12f %.12f \n"%(v[0], v[1], v[2])  )

            f.write("xcart  ")
            if len(st.xcart) != st.natom: 
                print_and_log("Warning! write_geometry(): xcart is empty or overfull, I make it from xred\n");#raise RuntimeError
                st.xcart = xred2xcart(st.xred, st.rprimd) 
            for v in st.xcart:
                f.write("%.12f %.12f %.12f \n"%(v[0]*le, v[1]*le, v[2]*le)  )

            if hasattr(st, 'select') and len(st.select) > 0 and not None in st.select:
                f.write("\nselect  ")
                for v in st.select:
                    # print(type(v[0]))
                    # sys.exit()
                    for i in 0,1,2:
                        if v[i] == 'T' or v[i] is True:
                            v[i] = 1
                        else:
                            v[i] = 0
                    f.write("{:d} {:d} {:d}\n".format(v[0], v[1], v[2])  )
                    # print(v)
            if hasattr(st, 'vel') and len(st.vel) > 0:
                f.write("\nvel ")
                for v in st.vel:
                    f.write("{:18.16f} {:18.16f} {:18.16f}\n".format(v[0], v[1], v[2])  )

            if hasattr(st, 'predictor') and st.predictor:
                
                f.write("\npred_length "+str(len(st.predictor)))
                f.write("\npredictor ")
                f.write(st.predictor)

            #Write build information
            try:
                self.build
            except AttributeError:
                pass
            else:
                f.write("\n\n\n#BEGIN BUILD INFORMATION!!!\n")
                #print self.build.__dict__
                for name in self.build.__dict__:
                    val = getattr(self.build, name)
                    if val == None or val == [None]: continue
                    if hasattr( val, '__iter__'):
                        temp = " ".join( map(str,val))
                        temp = temp.replace('[','')
                        temp = temp.replace(']','')
                        #print temp
                        f.write("%s %s\n" %(name,temp ))  # iterable
                    else:
                        f.write("%s %s\n" %(name, val)  )                    # not iterable
                f.write("\n#END BUILD INFORMATION!!!\n")


        return True


    def write_siman_geo(self, *args, **kwargs):
        """
        Please rename write_geometry() to write_abinit_geo() everywhere and transfer the code here
        """
        return self.write_geometry(*args, **kwargs)

    def serialize(self, filename):
        """
        save as pickle object, return path
        """
        file = filename+'.pickle'
        makedir(file)
        with open(file, 'wb') as f:
            pickle.dump(self, f, 2)
        return file

    def deserialize(self, filename, encoding = ''):
        
        # import chardet  
        # with open(filename, 'rb') as f:
        #     result = chardet.detect(f.read(10000))  
        # print(result)
        # sys.exit()
        with open(filename, 'rb') as f:
            if encoding:
                self = pickle.load(f, encoding = encoding)
            else:
                self = pickle.load(f, )
        # printlog('Calculation object succesfully read from ', filename)
        return self


    def serialize_json(self, filename):
        """
        save in json object - works
        the problem is how to decode correctly
        """
        cl = copy.deepcopy(self)
        for st in cl.init, cl.end:
            st.xcart = [list(xc) for xc in st.xcart]
            st.xred = [list(xc) for xc in st.xred]
            st.rprimd = [list(xc) for xc in st.rprimd]
            st.recip = [list(xc) for xc in st.recip]
        for mat in cl.occ_matrices:
            cl.occ_matrices[mat] = [list(line) for line in cl.occ_matrices[mat]]

        cl.ldauu = list(cl.ldauu)

        file = filename+'.json'
        makedir(file)
        # print(cl.__dict__)
        # print(cl.e0)
        with open(file, 'w') as f:
            json.dump(cl, f, default=lambda o: o.__dict__, 
            sort_keys=True, indent=4)
        


        # print(cl.__dict__)

        return file

    def deserialize_json(self, filename):
        """
        limited support, should be generalized
        """
        with open(filename, 'r') as fp:
            d = json.load(fp,) # works incorrect

        cl = CalculationVasp()
        
        sup = {}
        ats = ['set', 'init', 'end']
        for attr in ats: 
            sup[attr] = d[attr]
            del d[attr]
        cl.__dict__.update(d)
        # for attr in sup: 
            # print(sup[attr])
            # setattr(cl, attr+'.__dict__', sup[attr])
        cl.set.__dict__ = sup['set']
        cl.init.__dict__ = sup['init']
        cl.end.__dict__ = sup['end']
        # print(cl.end.rprimd)
        return cl
    def get_kpoints_density(self):
        """
        Number of k-points per atom
        """
        print(self.NKPTS*self.end.natom) #KPPRA - k-points per reciprocal atom? 




    def copy(self, id = None):
        # make entry with new id
        clcopy = copy.deepcopy(self)
        if id is not None:
            header.db[id] = clcopy
            header.db[id].id = id
        return clcopy

    def jmol(self, *args, **kwargs):
        self.end.jmol(*args, **kwargs)
    def poscar(self):
        self.end.write_poscar()

    def me(self):
        self.end.printme()
    def gmt(self, *args, **kwargs):
        return self.end.get_mag_tran(*args, **kwargs)


    def isggau(self):
        #check if calculation is gga+u
        #TODO - please finish
        ''


    def mag_diff(self, cl2, dm_skip = 0.5, el = 'FeNiCoVMnO', more = 0):

        """
        rms difference of magmom, skippting large deviations, due to defects
        dm_skip (float) - skip differences larger than this 
        el (str) - only for this element, could be several 
        more (bool) - show more info about mag diff at each pos
        the order of elements should be the same!!!

        """
        m1 = self.end.magmom
        m2 = cl2.end.magmom
        el1 = self.end.get_elements()
        el2 = cl2.end.get_elements()
        
        ms = 0
        tot=0
        maxdm = 0
        for i in range(len(m1)):
            if el1[i] != el2[i]:
                print('Warinig! el1 is not equal el2 for i=', i, el1[i], el2[i])
            
            dm = abs(m1[i]-m2[i])
            if dm > dm_skip:
                print('For i=', i, el1[i], 'dm= {:0.3f} muB'.format(dm), ', which is larger than dm_skip =', dm_skip, '; probably defect, skipping')
                continue
            if el and el1[i] not in el:
                continue

            if more:
                print('For ', i, el1[i], el2[i], 'dm= {:0.3f} muB'.format(dm))

            if dm > maxdm:
                maxdm = dm
            ms+= (dm)**2
            tot+=1
        rms = (ms/tot)**0.5

        if el:
            print('For '+el+' atoms RMS difference is {:0.3f} muB; dE is {:0.3f} eV'.format(rms, self.e0-cl2.e0))
            print('For '+el+' atoms max difference is {:0.3f} muB; dE is {:0.3f} eV'.format(maxdm, self.e0-cl2.e0))

        else:
            print('For rest atoms RMS difference is {:0.3f} muB; dE is {:0.3f} eV'.format(rms, self.e0-cl2.e0))
        mag1 = sum(self.end.magmom)
        mag2 = sum(cl2.end.magmom)
        mag_abs1 = sum([abs(m) for m in self.end.magmom])
        mag_abs2 = sum([abs(m) for m in cl2.end.magmom])

        el = 'NiCoO'
        w= 1
        st1 = self.end
        st2 = cl2.end
        suf_at1 = self.end.get_surface_atoms(el, surface = 0, surface_width=w)+self.end.get_surface_atoms(el, surface = 1, surface_width=w)
        suf_at2 = cl2.end.get_surface_atoms(el, surface = 0, surface_width=w)+cl2.end.get_surface_atoms(el, surface = 1, surface_width=w)
        # print(suf_at2)
        mag_sufabs1 = sum([abs(st1.magmom[i]) for i in suf_at1 ])
        mag_sufabs2 = sum([abs(st2.magmom[i]) for i in suf_at2  ])

        # self.end.jmol(r=2)

        print('Absolute magnetizations {:0.1f} muB, {:0.1f} muB'.format(mag_abs1, mag_abs2) )
        print('Absolute suf magnetizat {:0.1f} muB, {:0.1f} muB'.format(mag_sufabs1, mag_sufabs2) )
        print('Diff of summed magnetizations = {:0.1f} muB, total = {:0.1f} muB, absolute = {:0.1f} muB and abs suf = {:0.1f} muB'.format(mag1-mag2, self.mag_sum[-1][0]-cl2.mag_sum[-1][0], mag_abs1 - mag_abs2, mag_sufabs1 - mag_sufabs2) )


        return rms




    def occ_diff(self, cl2, li_at1 = None, li_at2 = None):
        """
        difference bettween occupation matricies for atoms li_at1 from self and li_at2 from cl2

        li_at1  -  list of atom numbers from self
        li_at2  -  list of atom numbers from cl2
        both lists should have the same length and the differences are taken between items at the same positions in lists


        otherwise

        self and cl2 should be commensurate, ideally having completly the same order of atoms
        first five for spin up
        then five for spin down



        """
        if li_at1 or li_at2:
            TM1 = li_at1
            TM2 = li_at2

        else:
            TM1 = self.end.get_transition_elements(fmt = 'n')
            TM2 = self.end.get_transition_elements(fmt = 'n')
        # print(TM)
        max_diff = 0.01
        nodiff  = True
        for i_at1, i_at2 in zip(TM1, TM2):
            occ1 = self.occ_matrices.get(i_at1)
            occ2 = cl2.occ_matrices.get(i_at2)

            if not occ1:
                print('Warning! no', i_at1, 'in self, skipping')
                continue
            if not occ2:
                print('Warning! no', i_at2, 'in cl2, skipping')
                continue

            occ1 = np.array(occ1)
            occ2 = np.array(occ2)

            # print(occ1-occ2)
            docc = occ1-occ2
            l05 = len(docc)//2

            # print(occ1[0:l05])
            det1 = np.linalg.det(docc[0:l05])
            det2 = np.linalg.det(docc[l05:])
            # m1 = np.matrix.max(np.matrix(docc))
            m1 = max(docc.min(), docc.max(), key=abs)
            # print(det1, det2, m1)
            df = pd.DataFrame(docc).round(5)

            if abs(m1) > max_diff:
                nodiff = False
                printlog('max diff larger than ', max_diff, 'was detected', imp = 'y')
                printlog('For atom ', i_at1, 'and atom', i_at2,  'max diff is ', '{:0.2f}'.format(m1), imp = 'y')
                printlog(tabulate(df, headers = ['dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2'], floatfmt=".2f", tablefmt='psql'),end = '\n', imp = 'Y'  )
        if nodiff:
            printlog('No diffs larger than', max_diff, '; Last matrix:', imp = 'y')
            printlog(tabulate(df, headers = ['dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2'], floatfmt=".2f", tablefmt='psql'),end = '\n', imp = 'Y'  )
        else:
            printlog('No more diffs', imp = 'y')


        return




    def dos(self, isym = None, el = None, i_at = None, iatoms = None,  *args, **kwargs):
        """
        Plot dos either for self or for children with dos
        isym (int) - choose symmetry position to plot DOS,
        otherwise use 
        i_at - number of atom from 0
        iatoms - list of atom numbers (from 0) to make one plot with several dos
        el - element for isym, otherwise first TM is used
        orbitals

        """
        from siman.header import db
        from siman.dos_functions import plot_dos
        # print(self.children)


        pm = kwargs
        x_nbins = pm.get('x_nbins')
        ylim = pm.get('ylim') or (-6,7)
        xlim = pm.get('xlim') or (-8,6)
        fontsize = pm.get('fontsize') or 13
        ver_lines = pm.get('ver_lines')
        corner_letter = pm.get('corner_letter')
        orbitals = pm.get('orbitals')
        efermi_origin = pm.get('efermi_origin')
        nsmooth = pm.get('nsmooth') or 0.0001
        linewidth = pm.get('linewidth') or 0.8
        efermi_shift = pm.get('efermi_shift') or 0
        labels = pm.get('labels')
        image_name = pm.get('image_name')
        fig_format = pm.get('fig_format') or 'pdf'

        if efermi_origin is None:
            efermi_origin = 1


        if corner_letter is None:
            corner_letter = 1
        # print(corner_letter)
        # sys.exit()

        ifdos = False
        if hasattr(self, 'children'):
            for id in self.children:
                # print(s[1])
                if 'dos' in id[1]:
                    printlog('Child with DOS set is found', id, imp = 'y')
                    id_dos = id
                    ifdos = True

                    break
            else:
                ifdos = False 
        if not ifdos:
            printlog('No children were found, using self', self.id, imp = 'y')
            id = self.id

        cl = db[id]

        cl.res()

        if orbitals is None:
            orbitals = ['d', 'p6']

        st = cl.end
        if isym is not None:

            if el:
                n = st.get_specific_elements(required_elements = [invert(el)], fmt = 'n')
            else:
                n = st.get_transition_elements(fmt = 'n')

            iTM = n[0]
            el = st.get_elements()[iTM]
            pos = determine_symmetry_positions(st, el)
            iTM = pos[isym][0]
            print('Choosing ', isym, 'atom ',iTM)
        else:
            iTM = i_at

        if not iatoms:
            #just one plot
            plot_dos(cl,  iatom = iTM+1,  
            dostype = 'partial', orbitals = orbitals, 
            labels = labels, 
            nsmooth = nsmooth, 
            image_name = image_name, 
            # invert_spins = invert_spins,
            efermi_origin = efermi_origin,
            efermi_shift = efermi_shift,
            show = 0,  plot_param = {
            'figsize': (6,3), 
            'linewidth':linewidth, 
            'fontsize':fontsize,
            'ylim':ylim, 'ver':1, 'fill':1,
            # 'ylim':(-1,1), 
            'ver_lines':ver_lines,
            'xlim':xlim, 
            'x_nbins':x_nbins,
            # 'xlim':(-0.5,0.1), 
            'dashes':(5,1), 'fig_format':fig_format, 'fontsize':fontsize})

        if iatoms:


            if fontsize:
                # header.mpl.rcParams.update({'font.size': fontsize+4})
                # fontsize = 2
                SMALL_SIZE = fontsize
                MEDIUM_SIZE = fontsize
                BIGGER_SIZE = fontsize

                header.mpl.rc('font', size=SMALL_SIZE)          # controls default text sizes
                header.mpl.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
                header.mpl.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
                header.mpl.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
                header.mpl.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
                header.mpl.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
                header.mpl.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


            color_dicts = [None, {'s':'k', 'p':'r', 'p6':'#FF0018', 'd':'g'}]

            total = len(iatoms)*1
            # letters = ['(a)', '(b)', '(c)', '(d)']*10
            letters = [str(i) for i in iatoms]
            font = 8
            fig, axs = plt.subplots(total,1,figsize=(6,total*3))    
            fig.text(0.03, 0.5, 'PDOS (states/atom/eV)', size = font*1.8, ha='center', va='center', rotation='vertical')
        
            i = 0
            first = 0
            last = 0
            hide_xlabels = 1
            xlabel = None
            ylabel = None
            # i_last = 1

            for iat in iatoms:
            # for cl, iat in zip([RbVsd, KVsd, Vsd], [13, 61, 53]):

                ax = axs[i]

                if corner_letter:
                    letter = letters[i]
                else:
                    letter = None
                
                # print(letter)
                # sys.exit()
                if i == 0:
                    first = True
                    last = False
                if i == total-1:
                    last = True
                    hide_xlabels = 0
                    xlabel = "Energy (eV)"
                plot_dos(cl,  iatom = iat+1,  efermi_origin = efermi_origin,
                dostype = 'partial', orbitals = orbitals, 
                labels = ['', ''], 
                nsmooth = 1, 
                color_dict = color_dicts[i%2],
                image_name = image_name, 

                # invert_spins = invert_spins,
                show_gravity = (1, 'p6', (-10, 10)), 
                show = 0,  plot_param = {
                # 'figsize': (6,3), 
                'linewidth':linewidth, 
                'fontsize':fontsize, 'legend_fontsize':font+3,
                'first':first, 'last':last, 'ax':ax, 'pad':1, 'hide_xlabels':hide_xlabels,
                'xlabel':xlabel, 'ylabel':ylabel,
                'corner_letter':letter,
                'ylim':ylim, 'ver':1, 'fill':1,
                # 'ylim':(-1,1), 
                'ver_lines':ver_lines,
                'xlim':xlim, 
                'x_nbins':x_nbins,
                # 'xlim':(-0.5,0.1), 
                'dashes':(5,1), 'fig_format':fig_format})

                i+=1



        return 


    def plot_locpot(self, filename = None):
        'plot LOCPOT'
        from siman.chg.chg_func import chg_at_z_direct
        from siman.picture_functions import fit_and_plot

        z_coord1, elst1 =  chg_at_z_direct(self, filetype = 'LOCPOT', plot = 0)

        if filename:
            show = False
            filename='figs/'+filename
        else:
            show = True

        fit_and_plot(pot=(z_coord1, elst1, '-b', ),
            xlabel = 'Z coordinate, $\AA$', 
            ylabel = 'Potential, eV', legend = None, fontsize = 12,
            show = show, hor_lines = [{'y':elst1[0]}],
            filename = filename
            )









    def add_new_name(self, idd):
        """
        
        just adding new key in database for that calculation
        warning cl.id is updated; old name in db remains
        the calculation folder remains the same
        
        idd - key
        """
        if idd in header.db:
            printlog('Error! ',idd,'already used in database! Choose another name')
            
        header.db[idd] = self
        header.struct_des[idd[0]] = header.struct_des[self.id[0]]

    def check_kpoints(self, ngkpt = None):
        """
        The method updates init.ngkpt and ngkpt_dict_for_kspacings !!! as well provides possible options for it
        TODO probably the method should transfered to Structure?
        Attention: the order should be the same as in make_kpoints_file
        """
        struct_des = header.struct_des
        # to_ang_local = header.to_ang
        to_ang_local = 1
        
        # try:
        #     if "Ang" in self.len_units:
        #         to_ang_local = 1
        #         #print "units angs"
        # except AttributeError:
        #     print_and_log("Warning! no len_units for "+self.name+" calculation, I use Bohr \n") 
        
        N_from_kspacing = []

        it = self.id[0]


        if not hasattr(struct_des[it], 'ngkpt_dict_for_kspacings'): #compatibiliy issues
            struct_des[it].ngkpt_dict_for_kspacings = {}

        ngkpt_dict = struct_des[it].ngkpt_dict_for_kspacings

        # if self.set.kpoints_file  == False:#self.set.vasp_params['KSPACING']:
        #     N = N_from_kspacing
        kspacing = self.set.vasp_params['KSPACING']
        # print(kspacing)
        # sys.exit()
        # print (struct_des)
        if ngkpt:
            N = ngkpt

        elif kspacing in ngkpt_dict:
            N = ngkpt_dict[kspacing]
            printlog('check_kpoints(): k-points will be used from *ngkpt_dict* of',it, N)
        
        elif self.set.ngkpt:
            N = self.set.ngkpt
            printlog('check_kpoints(): k-points will be used from set.ngkpt of',self.set.ise)


        elif kspacing:
            # print(self.init.rprimd)
            self.init.recip = self.init.get_recip()
            N_from_kspacing = calc_ngkpt(self.init.recip, kspacing)

            N = N_from_kspacing
            printlog('check_kpoints(): k-points are determined from kspacing',kspacing)

        elif self.set.kpoints_file:
            print_and_log("K-points file was provided", self.set.kpoints_file)
            N = None

        else:
            # print(self.dir)
            N = None
            if self.set.periodic:
                print_and_log("Error! check_kpoints(): no information about k-points for periodic calculation\n")



        self.init.ngkpt = N

        if kspacing != None and kspacing not in ngkpt_dict:
            ngkpt_dict[kspacing] = N
            printlog('check_kpoints(): I added ',N,'as a k-grid for',kspacing,'in struct_des of', it)


        print_and_log("check_kpoints(): Kpoint   mesh is: ", N, imp = 'Y')


        if not hasattr(struct_des[it], 'ngkpt_dict_for_kspacings') or  kspacing not in struct_des[it].ngkpt_dict_for_kspacings:
            print_and_log('Several other options instead of automatically determined ngkpt = ',N,np.array(self.calc_kspacings(N) ).round(2), ':', end = '\n', imp = 'y')
            print_and_log('ngkpt              |    actual kspacings       ', end = '\n', imp = 'y' )
            
            if N:
                for ngkpt in itertools.product([N[0]-1, N[0], N[0]+1], [N[1]-1, N[1], N[1]+1], [N[2]-1, N[2], N[2]+1]):
                    print_and_log(ngkpt, np.array(self.calc_kspacings(ngkpt) ).round(2), end = '\n', imp = 'y' )

            # user_ngkpt = input('Provide ngkpt:')
            # print(user_ngkpt)
            # sys.exit()

        else:
            print_and_log("check_kpoints(): The actual k-spacings are ", np.array(self.calc_kspacings(N) ).round(2), imp = 'Y')
        return N


    def calc_kspacings(self, ngkpt = None, sttype = 'init'):
        """Calculates reciprocal vectors and kspacing from ngkpt"""
        # to_ang_local = header.to_ang
        # try:
        #     if "Ang" in self.len_units:
        #         to_ang_local = 1
        #         #print "units angs"
        # except AttributeError:
        #     print_and_log("Warning! no len_units for "+self.name+" calculation, I use Bohr \n")
        

        if sttype == 'init':
            st = self.init
        if sttype == 'end':
            st = self.end 


        self.kspacing = []
        st.kspacings = []

        if not ngkpt:
            ngkpt = self.set.ngkpt

        k = [0,0,0]

        if ngkpt:
            k = calc_kspacings(ngkpt, st.rprimd)
            self.kspacing = copy.deepcopy(k)
            st.kspacing   = copy.deepcopy(k)

        return  k

    def actualize_set(self, curset = None, params = None):
        """
        Makes additional processing of set parameters, which also depends on calculation
    
        adding parameters for atat

        """


        #check if some parameters should be filled according to number of species
        #make element list
        el_list = [element_name_inv(el) for el in self.init.znucl]
        if not curset:
            curset = self.set
        vp = curset.vasp_params

        # print(['LDAU'])
        # print(vp)
        # print(vp['LDAU'])

        if 'LDAUL' in vp and vp['LDAUL'] is not None: 
            # print(vp['LDAU'])
            # if 
            for key in ['LDAUL', 'LDAUU', 'LDAUJ']:
                # print( vp[key])
                try:
                    if set(vp[key].keys()).isdisjoint(set(el_list)): #no common elements at all
                        print_and_log('\n\n\nAttention! The '+str(key)+' doesnt not contain values for your elements! Setting to zero\n\n\n')
                        # raise RuntimeError

                    new = []
                    for el in el_list:
                        
                        if el in vp[key]:
                            val = vp[key][el]
                            
                            if 'S' in el_list:  # use another value in the format of Fe/S
                                kk = el+'/S' 
                                if kk in vp[key]:
                                    val = vp[key][kk]
                    





                        else:
                            if key == 'LDAUL':
                                val = -1
                            else:
                                val =  0

                        new.append(val)
                    
                    vp[key] = new
                
                except AttributeError:
                    printlog('Error! LDAU* were not processed')
                    pass

        """Process magnetic moments"""
        if self.calc_method and 'afm_ordering' in self.calc_method:
            self.init.magmom = [None]



        # print(hasattr(self.init, 'magmom') and hasattr(self.init.magmom, '__iter__') and not None in self.init.magmom)
        # print(self.init.magmom)
        # print(None in self.init.magmom)
        # sys.exit()
        if hasattr(self.init, 'magmom') and hasattr(self.init.magmom, '__iter__') and not None in self.init.magmom and bool(self.init.magmom):

            print_and_log('actualize_set(): Magnetic moments are determined from self.init.magmom:',self.init.magmom, imp = 'y')

        elif hasattr(curset, 'magnetic_moments') and curset.magnetic_moments:
            print_and_log('actualize_set(): Magnetic moments are determined using siman key "magnetic_moments" and corresponding dict in set', end = '\n')
            print_and_log('curset.magnetic_moments = ', curset.magnetic_moments)
            
            mag_mom_other = 0.6 # magnetic moment for all other elements
            magmom = []
            for iat in range(self.init.natom):
                typ = self.init.typat[iat]
                el  = el_list[typ-1]
                if el in curset.magnetic_moments:
                    magmom.append(curset.magnetic_moments[el])
                else:
                    magmom.append(mag_mom_other)
            

            #convert magmom to vasp ordering

            zmagmom = [[] for x in range(0,self.init.ntypat)]

            # print zmagmom

            for t, m in zip(self.init.typat, magmom):
                # print "t, m = ", t, m
                zmagmom[t-1].append(m)
                # print t-1, zmagmom[3]

            # print 'sdfsdf', zmagmom[3], 
            poscar_ordered_magmom = [m for mag in zmagmom for m in mag  ]
            # sys.exit()
               
            vp['MAGMOM'] = poscar_ordered_magmom

            #check possible antiferromagnetic configurations:
            spec_mom_is = []
            for i, m in enumerate(magmom):
                if m != mag_mom_other: #detected some specific moment
                    spec_mom_is.append(i)

            if len(spec_mom_is) % 2 == 0 and len(spec_mom_is) > 0:
                print_and_log('Number of elements is even! trying to find all antiferromagnetic orderings:', imp = 'y')
                ns = len(spec_mom_is); 
                number_of_ord = int(math.factorial(ns) / math.factorial(0.5 * ns)**2)
                
                if number_of_ord > 10000:
                    printlog('Attention! Too much orderings (1000), skipping ...')
                else:
                    nords = 71
                    use_each = number_of_ord // nords  # spin() should be improved to find the AFM state based on the number of configuration 
                    if use_each == 0:
                        use_each = 1

                    if number_of_ord > nords:
                        print_and_log('Attention! Number of orderings is', number_of_ord, ' more than', nords, ' - I will check only each first ', imp = 'y')
                # else:

                    ls = [0]*len(spec_mom_is)
                    # print ls
                    orderings = []
                    



                    def spin(ls, i):
                        """
                        Find recursivly all possible orderings
                        ls - initial list of mag moments
                        i - index in ls  

                        """
                        # nonlocal i_current
                        if len(orderings) < nords:

                            for s in 1,-1:
                                
                                ls[i] = s
                                
                                if i < len(ls)-1:
                                
                                    spin(ls, i+1)
                                
                                else:
                                    if sum(ls) == 0:
                                        i_current['a']+=1  
                                        # print (i_current)

                                        if 1: #  i_current % use_each == 0:  # every use_each will be calculated; two slow even for sampling!
                                            orderings.append(copy.deepcopy(ls) )  
                                            # print (i_current)
                        return

                    i_current = {'a':0}
                    spin(ls, 0)

                    mag_orderings = []
                    mag_orderings.append(magmom)
                    printlog('Only '+str(nords)+' orderings equally sampled along the whole space are checked !')




                    for j, order in enumerate(orderings):
                        
                        # if j >nords: # old behaviour - just first ten orderings were checked
                        #     break

                        new_magmom = copy.deepcopy(magmom)
                        for i, s in zip(spec_mom_is, order):
                            # print i
                            new_magmom[i] = s * magmom[i]
                        

                        printlog(j, new_magmom,)
                        
                        mag_orderings.append(new_magmom)

                    # print orderings
                    print_and_log('Total number of orderings is ', len(orderings),imp = 'y')
                    
                    if self.calc_method and 'afm_ordering' in self.calc_method:
                        self.magnetic_orderings = mag_orderings
                  
            self.init.magmom = magmom # the order is the same as for other lists in init

        
        elif 'MAGMOM' in vp and vp['MAGMOM']: #just add * to magmom tag if it is provided without it
            print_and_log('Magnetic moments from vasp_params["MAGMOM"] are used\n')
            
            # if "*" not in vp['MAGMOM']:
            #     vp['MAGMOM'] = str(natom) +"*"+ vp['MAGMOM']
        


        # print (self.init.magmom, 'asdfaaaaaaaaaaaa')
        
        # sys.exit()

        # number of electrons

        if vp.get('MAGATOM') is not None: # for ATAT
            # print (vp['MAGATOM'])
            del vp['MAGMOM']
            # self.init.magmom = [None]
            # sys.exit()

        if self.calculator == 'aims':
            if None not in self.init.magmom:
                ''
                # vp['default_initial_moment'] = 0.6 # per atom - not good, since for different elements you need different moments




        return

    def write_sge_script(self, input_geofile = "header", version = 1, option = None, 
        prevcalcver = None, savefile = None, schedule_system = None,
        output_files_names = None,
        mode = None,
        batch_script_filename = None):
        """Without arguments writes header, else adds sequence of calculatios
            option - the same as inherit_option, 'inherit_xred' - control inheritance, or 'master' - run serial on master 
            prevcalcver - ver of previous calc; for first none
            savefile - 'cdawx', where c-charge, d-dos, a- AECCAR, w-wavefile, x-xml
            schedule_system - type of job scheduling system:'PBS', 'SGE', 'SLURM', 
                'none' - just run without any system
            mode - 
                body
                footer
        """


        # print 'Starting write_sge()', input_geofile
        varset = header.varset
        
        f = open(batch_script_filename,'a', newline = '') #

        # print(savefile)
        # sys.exit()


        def prepare_input(prevcalcver = None, option = None, input_geofile = None, name_mod_prev = '', write = True, 
            curver = None, copy_poscar_flag = True):
            """1. Input files preparation

                curver - current version
            """  


            if write:
                # if not 'only_neb' in self.calc_method:
                precont = str(prevcalcver)+name_mod_prev+'.CONTCAR' #previous contcar
                if option == 'inherit_xred' and prevcalcver:
                    if copy_poscar_flag:

                        f.write('grep -A '+str(self.init.natom)+ ' "Direct" '+precont+' >> '+input_geofile+ ' \n')

                if copy_poscar_flag: # only for first set 
                    if option == 'continue': #test for the case of sequence set - OK
                        ''
                        precont = str(curver)+name_mod_prev+'.CONTCAR ' #previous contcar
                        preout  = str(curver)+name_mod_prev+'.OUTCAR ' #previous outcar
                        f.write("cp "+precont+" POSCAR  # inherit_option = continue\n")
                        f.write("cp "+preout+'prev.'+preout+" # inherit_option = continue\n")
                        f.write('mv CHGCAR prev.CHGCAR   # inherit_option = continue\n')
                    
                    else:
                        f.write("cp "+input_geofile+" POSCAR\n")
        



            return





        def run_command(option, name, parrallel_run_command, condition = False, write = True):
            """2. write commands for running vasp. condition = true allows override additional conditions""" 

            if write:
                # if not condition:
                #     condition = (not 'only_neb' in self.calc_method)

                # if condition:


                if option == 'master':
                    f.write("vasp >"+name+".log\n")

                elif 'monte' in self.calc_method:
                    f.write("python "+header.cluster_home+'/'+ header.cluster_tools+'/siman/monte.py > monte.log\n')

                elif 'polaron' in self.calc_method:
                    f.write("python "+header.cluster_home+'/'+ header.cluster_tools+'/siman/polaron.py > polaron.log\n')

                elif 'atat' in  self.calc_method:
                    f.write('maps -d&\npollmach runstruct_vasp mpirun\n')


                else:
                    f.write(parrallel_run_command +" >"+name+".log\n")
                



                f.write("sleep 20\n")
            return


        def mv_files_according_versions(savefile, v, name_mod = '', write = True, rm_chg_wav = 'cw',
            ):    
            """3. Out files saving block
                
                rm_chg_wav - if True than CHGCAR and WAVECAR are removed

                savefile (str) - key, which determines what files should be saved
                    'o' - OUTCAR
                    'i' - INCAR
                    'v' - CHG
                    'c' - CHGCAR
                    'p' - PARCHG
                    'l' - LOCPOT
                    'd' - DOSCAR
                    'a' - AECCAR0, AECCAR2
                    'x' - vasprun.xml
                    't' - XDATCAR
                    'z' - OSZICAR
                    'w' - WAVECAR

            """   
            printlog('The value of savefile is', savefile)
            
            if 'polaron' in self.calc_method:
                write = 0 # not needed, since files are automatically saved by python script

            pre = v + name_mod
            contcar = pre+'.CONTCAR'

            if write:

                if "o" in savefile:

                    f.write("mv OUTCAR "          + v + name_mod +  ".OUTCAR\n")
                    f.write("mv CONTCAR "         + contcar+'\n')

                if "i" in savefile:
                    f.write("cp INCAR "           + v + name_mod +  ".INCAR\n")
                
                if "v" in savefile: # v means visualization chgcar
                    chg  = pre + '.CHG'
                    f.write("mv CHG "+chg+"\n")
                    f.write("gzip -f "+chg+"\n")
                
                if 'c' in savefile: # 
                    fln = 'CHGCAR'
                    chgcar  = pre +'.'+fln
                    f.write('cp '+fln+' '+chgcar+'\n') #use cp, cause it may be needed for other calcs in run
                    f.write('gzip -f '+chgcar+'\n')                

                if 'p' in savefile: # 
                    fln = 'PARCHG'
                    parchg  = pre +'.'+fln
                    f.write('cp '+fln+' '+parchg+'\n') #use cp, cause it may be needed for other calcs in run
                    f.write('gzip -f '+parchg+'\n') 

                if 'l' in savefile: # 
                    fln = 'LOCPOT'
                    locpot  = pre +'.'+fln
                    f.write('cp '+fln+' '+locpot+'\n') #use cp, cause it may be needed for other calcs in run
                    f.write('gzip -f '+locpot+'\n') 


                # else:
                #     f.write("rm CHG \n") #file can be used only for visualization


                if "d" in savefile:
                    fln = 'DOSCAR'
                    doscar  = pre +'.'+fln
                    f.write('mv '+fln+' '+doscar+'\n')
                    f.write('gzip -f '+doscar+'\n')                
                


                if "a" in savefile:
                    f.write("mv AECCAR0 "     + v + name_mod + ".AECCAR0\n")
                    f.write("mv AECCAR2 "     + v + name_mod + ".AECCAR2\n")
                
                if 'x' in savefile:
                    f.write("mv vasprun.xml " + v + name_mod + ".vasprun.xml\n")

                if 't' in savefile:
                    f.write("mv XDATCAR " + v + name_mod + ".XDATCAR\n")

                if 'z' in savefile:
                    f.write("mv OSZICAR " + v + name_mod + ".OSZICAR\n")
               
               
               
                if 'w' in savefile:
                    fln = 'WAVECAR'
                    wavecar  = pre +'.'+fln
                    # f.write("mv WAVECAR "     + v + name_mod + ".WAVECAR\n")
                    f.write('mv '+fln+' '+wavecar+'\n') #
                    f.write('gzip -f '+wavecar+'\n')  
                    rm_chg_wav = rm_chg_wav.replace('w','')
                # else:
                #     f.write("rm WAVECAR\n")


                if 'c' in rm_chg_wav:
                    f.write("rm CHGCAR   # rm_chg_wav flag\n") #file is important for continuation
                if 'w' in rm_chg_wav:
                    ''
                    f.write("rm WAVECAR  # rm_chg_wav flag\n") #
                if 'v' in rm_chg_wav: #chgcar for visualization
                    ''
                    f.write("rm CHG   # rm_chg_wav flag\n") #


            return contcar


        def analysis_script(write = True):
            #now only for u-ramping
            if write:
                f.write("touch ENERGIES\n")                

                for outcar in self.associated_outcars:
                    f.write("grep 'energy  without entropy' "+outcar+" | awk '{print $7}' >> ENERGIES\n")


        def name_mod_U_last():
            name_mod_last = 'U'+str(
                        update_incar(parameter = 'LDAUU', 
                            u_ramp_step = self.set.u_ramping_nstep-1, write = False, f = f, st = self )).replace('.','') #used to det last U

            return name_mod_last


        if schedule_system == 'SGE':
            # parrallel_run_command = "mpirun -x PATH vasp" # MPIE
            parrallel_run_command = header.vasp_command
        elif schedule_system in ['PBS', 'PBS_bsu', 'none']:
            # parrallel_run_command = "mpiexec --prefix /home/aleksenov_d/mpi/openmpi-1.6.3/installed vasp" bsu cluster
            # parrallel_run_command = "mpirun  vasp_std" #skoltech cluster
            parrallel_run_command = header.vasp_command #skoltech cluster
        
        elif schedule_system == 'SLURM':
            # parrallel_run_command = "prun /opt/vasp/bin/vasp5.4.1MPI"
            parrallel_run_command = header.vasp_command
        else:
            raise RuntimeError


        run_name = batch_script_filename     
        job_name = self.id[0]+"."+self.id[1]
        neb_flag = ('neb' in self.calc_method or 'only_neb' in self.calc_method)

        if hasattr(self.set, 'set_sequence') and self.set.set_sequence and any(self.set.set_sequence):
            sets = [self.set]+[se for se in self.set.set_sequence]
        else:
            sets = [self.set]




 
        def write_body(v = None, savefile = None, set_mod = '', copy_poscar_flag = True,
            final_analysis_flag = True, penult_set_name = None, curset = None):
            """
            set_mod (str) - additional modification of names needed for *set_sequence* regime, should be '.setname'
            """
            if 'only_neb' in self.calc_method:
                write = False
                write_poscar = False
            else:
                write = True
                write_poscar = True                

            #neb
            if 'neb' in self.calc_method: 
                if write: 
                    f.write("#NEB run, start and final configurations, then IMAGES:\n") 
                update_incar(parameter = 'IMAGES', value = 0, write = write, f = f, st = self) # start and final runs

            
            if 0: #experimental preliminary non-magnetic run
                ''
                #     if self.set.vasp_params['ISPIN'] == 2:
                #         print_and_log('Magnetic calculation detected; For better convergence',
                #          'I add first non-magnetic run', imp = 'Y')
                #         write = True
                #         name_mod_last = '.'+'NM'
                #         name_mod = '.NM'

                #         if write: 
                #             f.write("#Preliminary non-magnetic run:\n")  
                #         prepare_input(prevcalcver = prevcalcver, option = option,
                #          input_geofile = input_geofile, name_mod_prev = name_mod_last, write = write, curver = version)

                #         update_incar(parameter = 'ISPIN', value = 1, write = write) #
                        
                #         run_command(option = option, name = self.name+name_mod, parrallel_run_command = parrallel_run_command, write = write)

                #         if write:
                #             f.write("cp CONTCAR POSCAR  #prepare for basic run\n")
                #             write_poscar = False  

                #         mv_files_according_versions('co', v, name_mod = name_mod, write = write, rm_chg_wav = '')

                #         update_incar(parameter = 'ISPIN', value = 2, write = write) #




            if 'u_ramping' in self.calc_method:

                if write: 
                    f.write("#U-ramping run:\n")  

                # name_mod_last = '.'+name_mod_U_last()
                name_mod_last = '.U00' #since u_ramp starts from u = 00, it is more correct to continue from 00
                if penult_set_name:
                    name_mod_last += '.'+penult_set_name #however, for multiset run, the structure for u=00 exists only
                                                      #for penult set or maybe even for the first only set
                    # print (name_mod_last, penult_set_name)
                    # sys.exit()
                    
                # print 'prevcalver', prevcalcver

                if write and copy_poscar_flag: 
                    f.write("rm CHGCAR    #u-ramp init from scratch\n")                

                prepare_input(prevcalcver = prevcalcver, option = option,
                 input_geofile = input_geofile, name_mod_prev = name_mod_last, write = write_poscar, curver = version,
                 copy_poscar_flag = copy_poscar_flag)

                if copy_poscar_flag:
                    usteps = range(self.set.u_ramping_nstep)
                else:
                    usteps = [self.set.u_ramping_nstep-1]  # now it the case of sequence_set for contin sets only the last U is used

                u_last = 100
                for i_u in usteps:

                    u = update_incar(parameter = 'LDAUU', u_ramp_step = i_u, write = write, f = f , st = self)
                    if u == u_last:
                        continue
                    name_mod   = '.U'+str(u).replace('.', '')+set_mod
                   
                    run_command(option = option, name = self.name+name_mod, parrallel_run_command = parrallel_run_command, write = write)

                    if write: 
                        if copy_poscar_flag:

                            f.write("cp CONTCAR POSCAR   #u-ramp preparation\n")                

            # print(savefile)
                    contcar_file = mv_files_according_versions(savefile, v, name_mod = name_mod, write = write, rm_chg_wav = '')
                


                    self.associated_outcars.append( v + name_mod +  ".OUTCAR"  )
                    u_last = u
                
                if final_analysis_flag:
                    rm_chg_wav = 'w' #The wavcar is removed for the sake of harddrive space
                
                else:
                    rm_chg_wav = ''

                if curset.save_last_wave:
                    save_last = 'cw'
                else:
                    save_last = 'c'


                mv_files_according_versions(savefile = save_last, v=v, name_mod = name_mod, rm_chg_wav = rm_chg_wav) #save more files for last U
                

                analysis_script(write = write)
                # print self.associated





            elif 'afm_ordering' in self.calc_method:

                #Comment - inherit_xred option is not available here
                f.write("rm CHGCAR\n")                
                if not savefile: 
                    savefile = 'o'

                for i, magmom in enumerate(self.magnetic_orderings):

                    name_mod   = '.AFM'+str(i)+set_mod

                    update_incar(parameter = 'MAGMOM', value = magmom, write = write, f = f, st = self)

                    prepare_input(prevcalcver = prevcalcver, option = option, input_geofile = input_geofile,
                        copy_poscar_flag = copy_poscar_flag)
                    
                    run_command(option = option, name = self.name+name_mod, parrallel_run_command = parrallel_run_command)

                    contcar_file = mv_files_according_versions(savefile, v, name_mod = name_mod)
                
                    self.associated_outcars.append( v + name_mod +  ".OUTCAR"  )

                analysis_script()
            

            else: #simple run
                
                if not savefile: 
                    savefile = 'vco'

                if write: 
                        f.write("#Basic run:\n")  

                name_mod   = set_mod
                name_mod_last = ''

                prepare_input(prevcalcver = prevcalcver, option = option, name_mod_prev = name_mod_last,
                    input_geofile = input_geofile, write = write_poscar, curver = version,
                    copy_poscar_flag = copy_poscar_flag)

                run_command(option = option, name = self.name+name_mod, parrallel_run_command = parrallel_run_command, write = write)

                if final_analysis_flag:
                    rm_chg_wav = 'w' #The wavcar is removed for the sake of harddrive space
                
                else:
                    rm_chg_wav = ''


                contcar_file = mv_files_according_versions(savefile, v, write = write, name_mod = name_mod, rm_chg_wav = rm_chg_wav)

                self.associated_outcars.append( v + name_mod +  ".OUTCAR"  )

            return contcar_file
        


        def write_footer(set_mod = '', run_tool_flag = True, savefile = None, final_analysis_flag = True):
            """footer"""
            
            def u_ramp_prepare():
                if 'u_ramping' in self.calc_method:
                    u = update_incar(parameter = 'LDAUU', u_ramp_step = self.set.u_ramping_nstep-1, write = False, f = f, st = self)
                    name_mod   = '.U'+str(u).replace('.', '')
                    # name_mod_last = name_mod_U_last()+'.'
                    name_mod_last = '.'+'U00' #since u_ramp starts from u = 00, it is more correct to continue from 00
                
                else:
                    name_mod_last = ''
                    name_mod   = ''                

                return name_mod, name_mod_last

            def u_ramp_loop(ver_prefix = '', subfolders = None, run_name_prefix = None, set_mod = ''):

                if not subfolders:
                    subfolders = ['.']

                

                if run_tool_flag:
                    usteps = range(self.set.u_ramping_nstep)
                else:
                    usteps = [self.set.u_ramping_nstep-1]  # now it the case of sequence_set for contin sets only the last U is used


                u_last = 100

                for i_u in usteps:


                    u = update_incar(parameter = 'LDAUU', u_ramp_step = i_u, write = 1,  f = f, st = self)
                    if u == u_last:
                        continue

                    name_mod   = ver_prefix+'U'+str(u).replace('.', '')+set_mod

                    
                    run_command(option = option, name = run_name_prefix+'.'+name_mod, 
                        parrallel_run_command = parrallel_run_command, write = True)
                    
                    u_last = u


                    for n_st in subfolders:

                        f.write('cp '+n_st+'/CONTCAR '+n_st+'/POSCAR'+'               #u_ramp_loop()\n' )
                        f.write('cp '+n_st+'/OUTCAR  '+n_st+'/'+name_mod+'.OUTCAR'+'  #u_ramp_loop()\n' )
                        contcar = name_mod+'.CONTCAR'
                        f.write('cp '+n_st+'/CONTCAR  '+n_st+'/'+contcar+'            #u_ramp_loop()\n' )

                        # self.associated_outcars.append( v + name_mod +  ".OUTCAR"  )


                return contcar

            subfolders = None
            contcar_file = None
            
            if neb_flag:
                print_and_log('Writing scripts for NEB method', important = 'n')
                nim = self.set.vasp_params['IMAGES']
                nim_str = str(nim)

                subfolders = []
                for n in range(1, nim+1):
                    if n < 10:
                        n_st = '0'+str(n)
                    elif n < 100:
                        n_st = str(n)
                    subfolders.append(n_st)


                name_mod, name_mod_last = u_ramp_prepare()


                start = '1'+name_mod+'.OUTCAR '
                final = '2'+name_mod+'.OUTCAR '



                f.write("\n\n#Starting NEB script \n")

                if option and 'continue' in option:
                    prevout = name_mod_last+'OUTCAR '

                    for n_st in subfolders:
                        f.write('cp '+n_st+'/'+prevout+n_st+'/'+'prev.'+prevout+'  # inherit_option = continue\n' )
                        f.write('cp '+n_st+'/CONTCAR '+n_st+'/POSCAR  # inherit_option = continue\n')
                        f.write('mv '+n_st+'/CHGCAR '+n_st+'/prev.CHGCAR   # inherit_option = continue\n')

                else:
                    
                    if run_tool_flag:
                        f.write('export PATH=$PATH:'+header.cluster_home+'/tools/vts/\n') #header.project_path_cluster

                        f.write('nebmake.pl '+ start.replace('OUT','CONT') + final.replace('OUT','CONT') + nim_str +' \n')


                if nim+1 < 10: 
                    nim_plus_one_str = '0'+str(nim+1)

                if run_tool_flag:
                    f.write('cp '+start +  '00/OUTCAR\n')
                    f.write('cp '+final +  nim_plus_one_str + '/OUTCAR\n' )


                update_incar(parameter = 'IMAGES', value = nim, write  =1, f  = f , st = self)


                if 'u_ramping' in self.calc_method:


                    contcar_file = u_ramp_loop(subfolders = subfolders, 
                        run_name_prefix = self.name+'.n_'+nim_str, 
                        set_mod = set_mod)
              

                else:

                    run_command(option = option, name = self.name+set_mod+'.n_'+nim_str+name_mod, 
                    parrallel_run_command = parrallel_run_command, write = True)
                    # print(set_mod)
                    # sys.exit()
                    if '.' in set_mod and set_mod[0] == '.':
                        set_mod_loc = set_mod[1:]
                    else:
                        set_mod_loc = set_mod

                    name_mod   = set_mod_loc
                    if name_mod:
                        contcar = name_mod+'.CONTCAR'
                        outcar  = name_mod+'.OUTCAR'
                        for n_st in subfolders:
                            f.write('cp '+n_st+'/OUTCAR  '+n_st+'/'+outcar  +'  #sequence set: save file\n' )
                            f.write('cp '+n_st+'/CONTCAR  '+n_st+'/'+contcar+'  #sequence set: save file\n' )
                    else:
                        contcar = 'CONTCAR'

                    contcar_file = contcar






                if final_analysis_flag:
                    # f.write('export PATH=$PATH:'+header.cluster_home+'/tools/gnuplot/bin/ \n')
                    # f.write(header.cluster_home+'/tools/vts/nebresults.pl  \n')
                    f.write('find . -name WAVECAR -delete\n')
                    f.write('find . -name PROCAR -delete\n')
                # for n in range



            # print (calc[id].calc_method )
            # sys.exit()
            if 'uniform_scale' in self.calc_method or 'c_scale' in self.calc_method:
                # print (input_geofile)
                name_mod = set_mod
                
                if run_tool_flag:
                    f.write("\n\n#Starting fitting tool \n")
                    outputs = [ os.path.basename(out) for out in output_files_names ]
                    # f.write('export PYTHONPATH=$PYTHONPATH:'+CLUSTER_PYTHONPATH+'\n')
                    # f.write('/home/aksenov/tools/fit_tool.py '+list2string(outputs)+'\n' )
                    f.write('python '+header.cluster_home+'/tools/fit_tool.py '+list2string(outputs)+'\n' )
                    

                    f.write('cp 100.POSCAR POSCAR \n')
                
                if 'u_ramping' in self.calc_method:
                    

                    contcar_file = u_ramp_loop(ver_prefix = '100.', run_name_prefix = self.id[0]+'.fitted', set_mod = set_mod)

                else:
                    if final_analysis_flag:
                        rm_chg_wav = 'w' #The wavcar is removed for the sake of harddrive space
                    
                    else:
                        rm_chg_wav = ''




                    run_command(option = option, name = self.id[0]+'.'+self.id[1]+'.100'+name_mod+'.fitted', 
                        parrallel_run_command = parrallel_run_command, write = True)

                    # print(final_analysis_flag)
                    # sys.exit()

                    contcar_file = mv_files_according_versions(savefile, '100', name_mod = name_mod, write = True, rm_chg_wav = rm_chg_wav)

                # sys.exit()


            #clean at the end
            if final_analysis_flag: 
                if header.final_vasp_clean:
                    f.write('rm LOCPOT CHGCAR CHG PROCAR DOSCAR OSZICAR PCDAT REPORT XDATCAR vasprun.xml\n')
                f.write('rm RUNNING\n')



            return contcar_file, subfolders




        nsets = len(sets)
        footer_flag = not set(self.calc_method).isdisjoint(['uniform_scale', 'neb', 'only_neb' ])




        if mode == "body": #control part of script
            self.associated_outcars = []

        penult_set_name = None

        for k, curset in enumerate(sets):
            
            if nsets > 1: #the incar name is modified during creation only if more than 1 set is detected
                if mode == 'body' or footer_flag:
                    f.write('\n#sequence set: '+curset.ise+' \n')
                    f.write('cp '+curset.ise+'.INCAR  INCAR\n')
                    if hasattr(curset, 'savefile') and len(curset.savefile) > 0:
                        savefile = curset.savefile 


                penult_set_name = sets[-2].ise
            

            if k < nsets-1:
                set_mod = '.'+curset.ise
                final_analysis_flag = False
            else: #last set
                set_mod = '' # the last step do not use modifications of names 
                final_analysis_flag = True #for footer



            if k == 0: # additional control of prepare_input routine and footer
                copy_poscar_flag = True # the flag is also used to detect first set
                run_tool_flag = True
            else:
                copy_poscar_flag = False
                run_tool_flag  = False

            if mode == "body":
                
                contcar_file = write_body( v = str(version), savefile = savefile, 
                    set_mod = set_mod, copy_poscar_flag = copy_poscar_flag, 
                    final_analysis_flag = final_analysis_flag, penult_set_name = penult_set_name, curset = curset)

                

            elif mode == 'footer':
                if copy_poscar_flag: 
                    f.write('\n#Footer section: \n')


                # print(savefile)
                # sys.exit()
                contcar_file, subfolders = write_footer(set_mod = set_mod, run_tool_flag = run_tool_flag, savefile = savefile,
                 final_analysis_flag = final_analysis_flag)

            
            if k < nsets-1 and contcar_file:
                if 'o' in savefile:
                    if neb_flag and mode == 'footer':
                        for n_st in subfolders:
                            f.write('cp '+n_st+'/'+contcar_file+' '+n_st+'/POSCAR  # sequence_set: preparation of input geo for next neb set\n' )
                    else:
                        f.write('cp '+contcar_file+' POSCAR  #sequence_set: preparation of input geo for next set\n')










        if hasattr(self, 'associated_outcars') and  self.associated_outcars:
            out = self.associated_outcars[-1]
        else:
            out = None
        # print 'write_sge() out=', out
        f.close()
        
        return  out#return OUTCAR name
    



    def make_run(self, schedule_system, run_name):
        """Generate run file

        INPUT:
            schedule_system - 
        """

        with open('run','a', newline = '') as f:

            if schedule_system == 'SGE':
                #'qsub -pe 'mpi*' NCORES -l CLUSTER_TAG script.parallel.sh' for mpi-jobs which should run on CLUSTER_TAG (cmmd or cmdft)
                #IMPORTANT: NCORES must be a multiple of 8(24) on cmdft(cmmd). 
                # f.write("qsub -pe 'mpi*' "+str(header.corenum)+" "+header.queue+" "+run_name+"\n") #str(self.set.np) #-l cmmd; on MPIE
                
                f.write("qsub "+" "+run_name+"\n") 
            
                # f.write('sleep 5\n')
                # runBash('chmod +x run')
            
            elif schedule_system in ['PBS', 'PBS_bsu']:
                if header.PATH2PROJECT == '':
                    header.PATH2PROJECT = '.'

                f.write("cd "+header.PATH2PROJECT+'/'+self.dir+"\n")
                f.write("qsub "+run_name.split('/')[-1]+"\n") 
                f.write("cd -\n")
                f.write('sleep 1\n')                        
            elif schedule_system in ['none']:
                if header.PATH2PROJECT == '':
                    header.PATH2PROJECT = '.'

                f.write("cd "+header.PATH2PROJECT+'/'+self.dir+"\n")
                f.write('./'+run_name.split('/')[-1]+"\n") 
                f.write("cd -\n")
                # f.write('sleep 1\n')     

            
            elif schedule_system == 'SLURM':
                f.write("squeue\n") 
                f.write("sbatch " + run_name+"\n") 
                # f.write("sbatch -p AMG " + run_name+"\n") 
            else:
                printlog('Error! Unknown schedule_system', schedule_system)
                



        printlog("\nRun file created\n")     
        return

    def calculate_nbands(self, curset, path_to_potcar = None, params = None):
        """Should be run after add_potcar()
            updates set, including number of electrons
        """
        #1 add additional information to set
        if not curset:
            curset = self.set
        vp = curset.vasp_params
        st = copy.deepcopy(self.init)
        st = st.remove_atoms(['void']) # remove voids

        if path_to_potcar:
            # path_to_potcar = self.dir+'/POTCAR'
            self.init.zval = []
            # print path_to_potcar
            for line in open(path_to_potcar,'r'):
                if "ZVAL" in line:
                    # print line
                    self.init.zval.append(float(line.split()[5]))
            
            try: 
                curset.add_nbands
            except AttributeError: 
                curset.add_nbands = None

            if curset.add_nbands != None:
                tve =0
                for i in range(st.ntypat):
                    # print self.init.zval
                    tve += self.init.zval[i] * st.nznucl[i] #number of electrons 
                    # print(self.init.zval[i], self.init.nznucl[i])
                nbands_min = math.ceil(tve / 2.)
                self.nbands = int ( round ( nbands_min * curset.add_nbands ) )
                # print(self.nbands)
                

                vp['NBANDS'] = self.nbands
                printlog('I found that at least', nbands_min, ' bands are required. I will use', self.nbands, 'bands; add_nbands = ', curset.add_nbands)





            if 'LSORBIT' in vp and vp['LSORBIT']:
                # print (vp)
                printlog('SOC calculation detected; increasing number of bands by two', imp = 'Y')
                vp['NBANDS']*=2




            if params and 'charge' in params and params['charge']:
                vp['NELECT'] = int(tve - params['charge'])


        else:
            printlog('Attention! No path_to_potcar! skipping NBANDS calculation')

        return

    def show_force(self,):
        force_prefix = ' tot '

        printlog("\n\nMax. F."+force_prefix+" (meV/A) = \n{:};".format(np.array([m[1] for m in self.maxforce_list ])[:]  ), imp = 'Y'  )

    def check_job_state(self):
        #check if job in queue or Running

        cl = self
        if header.check_job == 1:
            job_in_queue = ''
            if hasattr(cl,'schedule_system'):


                check_string =  cl.id[0]+'.'+cl.id[1]
                if 'SLURM' in cl.schedule_system:


                    job_in_queue = check_string in run_on_server("squeue -o '%o' ", cl.cluster['address'])
                    printlog(cl.id[0]+'.'+cl.id[1], 'is in queue or running?', job_in_queue)

                elif 'PBS' in cl.schedule_system:
                    job_in_queue = check_string in run_on_server("qstat -x ", cl.cluster['address'])

                elif 'SGE' in cl.schedule_system:
                    job_in_queue = check_string in run_on_server("qstat -xml ", cl.cluster['address'])
                
                elif 'none' in cl.schedule_system:
                    job_in_queue = ''
                    
                else:
                    print_and_log('Attention! unknown SCHEDULE_SYSTEM='+cl.schedule_system+'; Please teach me here! ', imp = 'y')
                    job_in_queue = ''


            if file_exists_on_server(os.path.join(cl.dir, 'RUNNING'), addr = cl.cluster['address']) and job_in_queue: 
                
                cl.state = '3. Running'

            elif job_in_queue:
                
                cl.state = '3. In queue'
       
            else:
                ''
                if '3' in cl.state:
                    cl.state = '2. Unknown'

        else:
            cl.state = '2. Unknown'




        return cl.state 



    def get_file(self, filetype = '', nametype = '', up = 'up1', root = 0):
        """
        allow to get any file of type filetype 
        cl - (Calculation) 
        filetype (str) - 'CHG', 'CHGCAR', etc just the name of file in calculation folder
        nametype (str) - 'asoutcar' - update filetype to OUTCAR format
        up (str) - control flag 
            '1' - do not update
            '2' - update

        root - root calculation folder location of file

        Comment
            initially used for chg files - rename!
        """

        setting_sshpass(self)
        # print(filetype)

        if nametype == 'asoutcar':
            path_to_file = self.path['output'].replace('OUTCAR',filetype)
        else:
            if root:
                path_to_file = self.dir +'/'+ filetype
            else:
                path_to_file = os.path.dirname(self.path['output']) +'/'+ filetype
        if 'CHGCAR' in filetype:
            self.path['chgcar'] = path_to_file
            self.path['charge'] = path_to_file
        elif 'xml' in filetype:
            self.path['xml'] = path_to_file


        # print(self.cluster_address)
        # print(self.project_path_cluster+'/')
        # sys.exit()
        if hasattr(self, 'cluster'):
            address = self.cluster['address']
        if header.override_cluster_address:
            if hasattr(self, 'cluster') and self.cluster.get('name'):
                clust = header.CLUSTERS[self.cluster['name']]
            else:
                printlog('Youve chosen to override cluster_address, but name of cluster is None, trying default', imp = 'Y')
                clust = header.CLUSTERS[header.DEFAULT_CLUSTER]

            self.project_path_cluster = clust['homepath']
            address = clust['address']


        path2file_cluster = self.project_path_cluster+'/'+path_to_file

        # print(self.project_path_cluster)
        # sys.exit()

        if os.path.exists(path_to_file) and '2' not in up: 
            out = None
        else:
            # printlog('File', path_to_file, 'was not found. Trying to update from server')
            out = get_from_server(path2file_cluster, os.path.dirname(path_to_file), addr = address)


        if out:
            printlog('File', path2file_cluster, 'was not found, trying archive:',header.PATH2ARCHIVE, imp = 'Y')
            # printlog('Charge file', path_to_file, 'was not found')
            try:
                pp = self.project_path_cluster.replace(self.cluster_home, '') #project path without home
            except:
                pp = ''
            # print(pp)
            path_to_file_scratch = header.PATH2ARCHIVE+'/'+pp+'/'+path_to_file

            out = get_from_server(path_to_file_scratch, os.path.dirname(path_to_file), addr = self.cluster['address'])
            
            if out:
                printlog('File', path_to_file_scratch, 'was not found', imp = 'Y')
                path_to_file = None
           
        printlog('File', path_to_file, ' was download', imp = 'e')
        
        return path_to_file

    def run_on_server(self, command, addr = None):
        setting_sshpass(self)
        if addr is None:
            addr = self.cluster['address']
        out = run_on_server(command, addr)
        return out

    def update_name(self): 
        self.name = str(self.id[0])+'.'+str(self.id[1])+'.'+str(self.id[2])
        return self.name




    @property
    def sfolder(self):
        self._x = header.struct_des[self.id[0]].sfolder
        return self._x

    def e0_fu(self, n_fu = None):
        # please improve
        #n_fu - number of atoms in formual unit
        if n_fu:
            n1 = self.end.natom/n_fu
            print('e0_fu: Normalization by provided n_fu', n_fu)

        else:
            self.end.get_nznucl()
            n1 = self.end.nznucl[0]
            print('e0_fu: Normalization by element z=', self.end.znucl[0])

        e0_fu = self.e0/n1
        print('e0_fu: e0_fu=',e0_fu)
        
        return e0_fu

    @property
    def e0_at(self,):
        return self.e0/self.end.natom




class CalculationAbinit(Calculation):
    """docstring for CalculationAbinit"""
    pass


class CalculationAims(Calculation):
    """object for Aims code """
    def __init__(self, inset = None, iid = None, output = None):
        super(CalculationAims, self).__init__(inset, iid, output)
        self.len_units = 'Angstrom'
        self.calculator = 'aims'

    def write_structure(self, name_of_output_file, type_of_coordinates = 'dir', option = None, prevcalcver = None, path = None, state = 'init'):

        if path == None: 
            path = self.dir
        
        if state == 'init':
            st  = self.init
        elif state == 'end':
            st  = self.end
        else: 
            raise RuntimeError 
        
        filename = os.path.join(path, 'geometry.in')

        makedir(filename)

        write_geometry_aims(st, filename, coord_type = type_of_coordinates, periodic = self.set.periodic)


    def add_potcar(self):

        d = self.dir

        incar = d+'control.in'

        with open(self.set.path_to_potcar, 'r') as f:
            fil = f.read()

        with open(incar, 'w') as f:
            f.write(fil)

        self.path['potcar'] = self.set.path_to_potcar

    def make_incar(self):
        d = self.dir
        
        incar = d+'control.in'
        with open(incar, 'r') as f:
            fil = f.read()
        vp = self.set.params
        
        N = self.check_kpoints()
        # print(N)
        # self.exit()
        if N:
            vp['k_grid'] = list2string(N)

        with open(incar, 'w') as f:
            f.write(vp['universal'])
            f.write('\n')
            for key in vp:
                if key in aims_keys:
                    # print(key, self.set.params[key])
                    if vp[key] is not None:
                        f.write(key+' '+str(vp[key])+'\n')
            f.write(fil)
        
        return [incar]

    def make_kpoints_file(self):
        printlog( "Attention! ngkpt for kpoints file are created from kspacing\n")
        N = self.check_kpoints()
        self.set.ngkpt = N
        return ['']


    def copy_to_cluster(self, list_to_copy, update):
        d = self.dir
        list_to_copy.extend( glob.glob(   os.path.join(d, '*geometry*')  ) )
        
        if "up" in update: #Copy to server
            printlog('Files to copy:', list_to_copy)

            push_to_server(list_to_copy,  self.project_path_cluster +'/'+ self.dir, self.cluster_address)

    def download(self, load):

        path_to_outcar  = self.path["output"]

        # print(path_to_outcar)
        # sys.exit()

        self.get_file(os.path.basename(path_to_outcar), up = load)

        return path_to_outcar

    def read_results(self, load = '', out_type = '', voronoi = None, show = '', choose_outcar = None, alkali_ion_number = None):

        """
        Aims

        choose_outcar - for now is dummy
        alkali_ion_number - for now is dummy
        voronoi - dummy
        """
        cl = self
        filename = cl.download(load) # wrapper for downloading output files


        cl.state = check_output(filename, 'Have a nice day', load)
        
        if "4" in cl.state:

            outstr = read_aims_out(cl, out_type, show)
            
            printlog(outstr)

        else:
            
            printlog('Status of calculation is', cl.state, 'continiue', imp = 'y')
            outstr = cl.state
        

        return outstr




class CalculationVasp(Calculation):
    """Methods for calculations made using VASP DFT code"""
    def __init__(self, inset = None, iid = None, output = None):
        super(CalculationVasp, self).__init__(inset, iid, output)
        self.len_units = 'Angstrom'
        self.calculator = 'vasp'



    def read_poscar(self, filename, version = None):
        """
        Read POSCAR file using st.read_poscar 
        

        """


        if self.path["input_geo"] == None:
            self.path["input_geo"] = filename
        self.path["poscar"] = filename
        

        self.hex_a = None
        self.hex_c = None
        self.gbpos = None


        #Determine version
        if version:
            self.version = version
        else:
            print_and_log('Trying to find version at the end of filename POSCAR-v ...')
            try:
                ver = int(filename.split('-')[-1])
                print_and_log('OK\n')

            except:
                print_and_log('\nTrying to find version at the begenning of filename v.POSCAR...')

                try:
                    ver = int(os.path.basename(filename).split('.')[0] )
                    print_and_log('OK\n')
               
                except:
                    print_and_log('No version, using 1\n')
                    ver = 1

            self.version = ver

        self.init = Structure()
        self.init = read_poscar(self.init, filename)
        self.des = self.init.des

        return
    



    def write_structure(self, name_of_output_file, type_of_coordinates = 'dir', option = None, prevcalcver = None, path = None, state = 'init'):
        """Generates POSCAR file
           type_of_coordinates - 'direct' (xred) or 'cartesian' (xcart)
           option -inheritance option
           prevcalcver - ver of first calc in calc list; for first None
           state - 'init' or 'end' 
        """
        #units
        # try:
        #     if "ang" in self.len_units or "Ang" in self.len_units: 
        #         global to_ang; to_ang = 1.0; print_and_log("Conversion multiplier to_ang is "+str(to_ang) )
        # except AttributeError:
        #     pass

        if option == 'inherit_xred' and 'car' in type_of_coordinates: 
            raise RuntimeError 

        if option == 'inherit_xred' and prevcalcver: 
            type_of_coordinates = 'None' # do not write xred or xcart if they will be transfered on cluster
        
        if path == None: 
            path = self.dir
        
        if state == 'init':
            st  = self.init
        elif state == 'end':
            st  = self.end
        else: 
            raise RuntimeError 
        
        filename = os.path.join(path, name_of_output_file)

        makedir(filename)

        st.write_poscar(filename, coord_type = type_of_coordinates)



        return



    def add_potcar(self):
        """Should be run for the first calculation only"""
        #Add POTCAR

        path_to_potcar = os.path.join(self.dir, 'POTCAR')
        potcar_files   = []

        if hasattr(self.set, 'path2pot' ) and self.set.path2pot:
            path2pot = self.set.path2pot
        else:
            path2pot = header.PATH2POTENTIALS
        printlog('Potentials from ', path2pot, 'are taken')

        if self.set.potdir:
            # print (self.set.potdir)
            for z in self.init.znucl:
                if z == 300:
                    continue # skip voids
                potcar_files.append(os.path.join(path2pot, self.set.potdir[ int(z) ], 'POTCAR') )

            printlog("POTCAR files:", potcar_files)        
            # print(path_to_potcar)            
            cat_files(potcar_files, path_to_potcar)

        


        elif self.set.path_to_potcar:
            printlog('Attention! set.path_to_potcar is used !', self.set.path_to_potcar)
            shutil.copyfile(self.set.path_to_potcar, path_to_potcar)
            printlog('POTCAR was copied to', path_to_potcar)
            path_to_potcar = self.set.path_to_potcar


        else:
            printlog('Error! set.potdir and set.path_to_potcar are empty; no POTCAR was not created!')
            path_to_potcar = None
        
        self.path['potcar'] = path_to_potcar

        return path_to_potcar












    def make_incar(self):
        """Makes Incar file for current calculation and copy all
        TO DO: there is no need to send all POSCAR files; It is enothg to send only one. However for rsync its not that crucial
        """
        #print "Begin make---------------------------------------------"
        
        
        #Generate incar
        varset = header.varset
        d = self.dir
        natom = self.init.natom
        poscar_atom_order = self.init.poscar_atom_order # order of atoms in POSCAR, can be different from init!!!! used for magmom
        incar_list = []

        setseq = [self.set]
        
        if hasattr(self.set, 'set_sequence') and self.set.set_sequence:
            for s in self.set.set_sequence:
                setseq.append(s)


        nsets = len(setseq)
        for i, curset in enumerate(setseq):

            if nsets == 1:
                name_mod = ''
            else:
                name_mod = curset.ise+'.'

            incar_filename = d+name_mod+'INCAR'
            vp = curset.vasp_params
            

            with open(incar_filename,'w', newline = '') as f:

                f.write( 'SYSTEM = ')
                if hasattr(self.init, 'perm'):
                    f.write( 'perm=[{:s}] ; '.format( list2string([i+1 for i in self.init.perm]).replace(' ', ',') )) #write permuations
                f.write( '{:s}\n'.format(self.des) )


                for key in sorted(vp):
    
                    if key == 'SYSTEM':
                        ''
                    elif key == 'MAGMOM' and hasattr(self.init, 'magmom') and self.init.magmom and any(self.init.magmom): #
                        mag = self.init.magmom
                        magmom_aligned_with_poscar = [mag[i] for i in poscar_atom_order ]
                        f.write('MAGMOM = '+list2string(magmom_aligned_with_poscar)+"\n") #magmom from geo file has higher preference
                   
                    elif vp[key] == None:
                        ''

                    elif key == 'KSPACING' and self.set.kpoints_file: #attention! for k-points only the base set is used!!
                        '' # since VASP has higher priority of KSPACING param, it should not be written 

                    elif is_list_like(vp[key]):
                        lis = vp[key]
                        f.write(key + " = " + ' '.join(['{:}']*len(lis)).format(*lis) + "\n")
                   
                    else:
                        f.write(key+" = "+str( vp[key] ) +"\n")
               

                f.write("\n")




            print_and_log(incar_filename, "was generated\n")
        
            incar_list.append(incar_filename)
        

        return incar_list




    
    def make_kpoints_file(self):

        struct_des = header.struct_des
        #Generate KPOINTS
        kspacing = self.set.vasp_params['KSPACING']

        filename = os.path.join(self.dir, "KPOINTS")

        it = self.id[0]



        if hasattr(self.set, 'k_band_structure') and self.set.k_band_structure:
            k = self.set.k_band_structure
            printlog('Writing k-points file for band structure calculation.', imp = 'y')
            
            with open(filename, 'w', newline = '') as f:
                f.write('k-points along high symmetry lines\n')
                f.write('{:} ! intersections\n'.format(k[0]))
                f.write('Line-mode\n')
                f.write('rec\n') # now only reciprocal are supported
                ps= k[1]
                for pn in k[2:]:
                    # pn  = next(k)
                    f.write('{:6.3f} {:6.3f} {:6.3f} ! {:s}\n'.format(ps[1], ps[2], ps[3], ps[0]) ) 
                    f.write('{:6.3f} {:6.3f} {:6.3f} ! {:s}\n\n'.format(pn[1], pn[2], pn[3], pn[0]) ) 
                    ps = pn





        elif self.set.kpoints_file:
            if self.set.kpoints_file == True:

                print_and_log( "You said to generate KPOINTS file.\n")
                self.calc_kspacings()
                #Generate kpoints file

                #
                if hasattr(struct_des[it], 'ngkpt_dict_for_kspacings') and kspacing in struct_des[it].ngkpt_dict_for_kspacings:
                    N =    struct_des[it].ngkpt_dict_for_kspacings[kspacing]
                    print_and_log( 'Attention! ngkpt = ',N, 
                        ' is adopted from struct_des which you provided for it ',it, ' and kspacing = ',kspacing)
                    nk1 = N[0]; nk2 = N[1]; nk3 = N[2]
                    self.set.ngkpt = N

                elif self.set.ngkpt:
                    nk1 = self.set.ngkpt[0]; nk2 = self.set.ngkpt[1]; nk3 = self.set.ngkpt[2];
                    print_and_log( "Attention! ngkpt was used for kpoints file\n")

                
                elif kspacing:
                    print_and_log( "Attention! ngkpt for kpoints file are created from kspacing; ngkpt is empty\n")
                    N = self.check_kpoints()
                    self.set.ngkpt = N
                    nk1 = N[0]; nk2 = N[1]; nk3 = N[2]
                
                else:
                    print_and_log( "Error! could not find information about k-points\n")



                with open(filename,'w', newline = '') as f:

                    f.write("Automatic Mesh\n") #Comment
                    f.write("0 \n")#Number of points; 0-Auto
                    if 'KGAMMA' in self.set.vasp_params and self.set.vasp_params["KGAMMA"] in (1,'.TRUE.', 'True', '1'): 
                        f.write("Gamma\n")
                    else: 
                        f.write("Monkhorst Pack\n")
                    f.write('%i %i %i \n'%(nk1, nk2, nk3) )
                    f.write("0 0 0\n") # optional shift

                print_and_log( "KPOINTS was generated\n")
            
            else:
                # print()
                shutil.copyfile(self.set.kpoints_file, filename)
                print_and_log( "KPOINTS was copied from"+self.set.kpoints_file+"\n")



        else:
            print_and_log( "This set is without KPOINTS file.\n")
            filename = ''



        return [filename]

    def copy_to_cluster(self, list_to_copy, update):
        d = self.dir
        list_to_copy.append( os.path.join(d, 'POTCAR')  )
        list_to_copy.extend( glob.glob(   os.path.join(d, '*POSCAR*')  ) )
        # list_to_copy.extend( glob.glob(   os.path.join(d, '*.run*'  )  ) )

        if 'OCCEXT' in self.set.vasp_params and self.set.vasp_params['OCCEXT'] == 1:
            list_to_copy.append(  os.path.join(d, 'OCCMATRIX')   )

        
        if "up" in update: #Copy to server
            printlog('Files to copy:', list_to_copy)

            # command = ' mkdir -p {:}'.format( os.path.join(self.project_path_cluster, self.dir)  )

            # run_on_server(command, self.cluster_address)

            push_to_server(list_to_copy,  self.project_path_cluster +'/'+ self.dir, self.cluster_address)


        return








    def plot_energy_force(self, force_type = 'max'):
        # print(self.maxforce)
        if 'max' in force_type:
            force = [m[1] for m in self.maxforce_list ]
            lab = 'Max.'
        elif 'av' in force_type:
            # print(self.average_list)
            force = [m for m in self.average_list ]
            lab = 'Av.'


        # print(maxf)
        plt.plot(force, 1000*(np.array(self.list_e_sigma0)-self.energy_sigma0) , '-o')
        # plt.xlabel('MD step')
        # plt.ylabel('Energy per cell (eV')
        plt.xlabel(lab+' force on atom (meV/$\AA$)')
        plt.ylabel('Energy per cell relative to min (meV)')

        plt.show()
        return

    def plot_energy_step(self,):
        # print(self.maxforce)
        # maxf = [m[1] for m in self.maxforce_list ]
        # print(maxf)
        steps = range(len(self.list_e_sigma0))
        plt.plot(steps, 1000*(np.array(self.list_e_sigma0)-self.energy_sigma0) , '-o')
        # plt.xlabel('MD step')
        # plt.ylabel('Energy per cell (eV')
        plt.xlabel('Step')
        plt.ylabel('Energy per cell relative to min (meV)')

        plt.show()
        return

    def plot_energy_conv(self,):
        # print(self.maxforce)
        # maxf = [m[1] for m in self.maxforce_list ]
        # print(maxf)
        en = self.list_e_conv[10:]
        steps = range(len(en))
        plt.plot(steps, (np.array(en)-self.energy_sigma0) , '-o')
        # plt.xlabel('MD step')
        # plt.ylabel('Energy per cell (eV')
        plt.xlabel('SCF Step')
        plt.ylabel('Energy per cell relative to min (eV)')

        plt.show()
        return




    def read_results(self, load = '', out_type = '', voronoi = False, show = '', choose_outcar = None, alkali_ion_number = None, only_load = False):

        """
        Download and Read VASP OUTCAR file

        ###INPUT:
            - load (str) - 'x' - download xml, o - download outcar and contcar, un - read unfinished
            - show (str) - print additional information
                alkali_ion_number - show mag around this ion
            - choose_outcar - see description in res_loop(), from 1
            - out_type - controls the return string
                see in code, add here
                also controls reading of OUTCAR
                'xcarts' read xcart every relaxation step and write into self.end.list_xcart

            - only_load (bool) - if true - only load the files (used for database)

        ###RETURN:
            ?

        ###DEPENDS:
        TODO:
        please split into outcar parser, downloader, and checker-formatter

        """

        # print (choose_outcar, hasattr(self, 'associated_outcars'), self.associated_outcars)
        join = os.path.join
        dirname = os.path.dirname

        if header.show:
            show +=header.show

        if not hasattr(self, 'dir'):
            self.dir = os.path.dirname(self.path['output'])

        # print(self.associated_outcars)

        # print(choose_outcar)
        # sys.exit()
        if choose_outcar and hasattr(self, 'associated_outcars') and self.associated_outcars and len(self.associated_outcars) >= choose_outcar and len(self.associated_outcars) > 1:
            # print ('associated outcars = ',self.associated_outcars)
            printlog('read_results(): choose_outcar', choose_outcar)

            path_to_outcar = join( dirname(self.path["output"]), self.associated_outcars[choose_outcar-1] )

            printlog(self.associated_outcars)
        else:
            path_to_outcar  = self.path["output"]
        
        # print(path_to_outcar)
        if 'OUTCAR' in path_to_outcar:
            path_to_contcar = path_to_outcar.replace('OUTCAR', "CONTCAR")
            path_to_poscar = path_to_outcar.replace('OUTCAR', "POSCAR")
            path_to_xml     = path_to_outcar.replace('OUTCAR', "vasprun.xml")
        else:
            path_to_contcar = ''
            path_to_xml     = ''


        if self.calc_method  and not set(self.calc_method  ).isdisjoint(  ['u_ramping', 'afm_ordering']):
            '' 
            # print self.associated_outcars
            # lor = self.associated_outcars[-1]
            # path_to_outcar  = self.dir + lor
            # path_to_contcar = self.dir + lor.replace('OUTCAR', 'CONTCAR')
            # path_to_xml     = self.dir + lor.replace('OUTCAR', 'vasprun.xml')         
            # print 'sdf', path_to_outcar, path_to_contcar, path_to_xml

            # energies_str = runBash("ssh "+self.cluster_address+" cat "+self.dir+"ENERGIES")
            
            # print "ssh "+self.cluster_address+" cat "+self.dir+"ENERGIES"
            # print (  energies_str )
            # if not 'cat' in energies_str:
            #     self.associated_energies = [float(e) for e in energies_str.split()]
            
            # self.u_ramping_u_values = np.arange(*self.u_ramping_list)
            # print 'associated_energies:', self.associated_energies
        print_and_log('read_results() path to outcar', path_to_outcar)
        # sys.exit()



        if not os.path.exists(path_to_outcar):
            load = load+'o'




        """Copy from server """

        printlog('The load flag is ', load)

        if 'o' in load and hasattr(self, 'cluster_address'):

            #reduce size of downloadable file by removing occupations: vasp 4 and 5
            command_reduce = """ssh {0:s} nbands=\`grep \\"NBANDS=\\" \{1:s} \| awk \\'{{print \$NF - 1}}\\'\`\; sed -i -e \\"/band No./,+\${{nbands}}d\\" \{1:s} """.format(
                self.cluster['address'], join(self.project_path_cluster, path_to_outcar) )


            # runBash(command_reduce)


            if 'un2' in load:
                out_name  = os.path.basename(path_to_outcar)
                cont_name = os.path.basename(path_to_contcar)
                path_to_outcar = path_to_outcar.replace(out_name, 'OUTCAR')
                path_to_contcar = path_to_contcar.replace(cont_name, 'CONTCAR')
                # self.path['output'] = path_to_outcar

            files = [ self.project_path_cluster+'/'+path_to_outcar, self.project_path_cluster+'/'+path_to_contcar ]
            # print(load)
            # print(files)
            # get_from_server(files = files, to = os.path.dirname(path_to_outcar),  addr = self.cluster_address)
            for file in files:
                self.get_file(os.path.basename(file), up = load)


        if 'x' in load:

            # get_from_server(files = join(self.project_path_cluster, path_to_xml), to = os.path.dirname(path_to_outcar),  
            #     addr = self.cluster_address)
            
            self.get_file(os.path.basename(path_to_xml), up = load)





        if os.path.exists(path_to_outcar):
            outcar_exist   = True
        else:
            outcar_exist   = False
            path_to_zip = path_to_outcar+'.gz'
            if os.path.exists(path_to_zip):
                with gzip.open(path_to_zip, 'rb') as f_in:      # unzip OUTCAR
                    with open(path_to_outcar, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                if os.path.exists(path_to_outcar):
                    outcar_exist = True


        """Start reading """

        self.state = check_output(path_to_outcar, 'General timing', load)
        
        outst = self.state


        if "4" in self.state:
            
            outst  = read_vasp_out(self, load = load, out_type = out_type, show = show, voronoi = voronoi,
                path_to_outcar = path_to_outcar, path_to_contcar = path_to_contcar)


        else:
            cl = self
            try:
                os.rename(cl.path['output'], cl.path['output']+"_unfinished") 
                printlog('read_results():',cl.id, 'is unfinished, continue:', cl.dir, cl.cluster['address'], imp = 'y')
                cl.state = '5. Unfinished'
            except:
                printlog('read_results():',cl.id, 'probably was not submitted:', cl.dir, imp = 'y')


        return outst




    def determine_filenames(self, nametype = 'asoutcar'):
        """
        try to determine correct filenames
        """
        if nametype == 'asoutcar':
            for filetype in 'CHGCAR', 'AECCAR0', 'AECCAR2':
                self.path[filetype.lower()] = self.path['output'].replace('OUTCAR',filetype)
                printlog('determine_filenames()',self.path[filetype.lower()])





    def get_chg_file(self, *args, **kwargs):
        """just wrapper to get chgcar files """
        if 'CHGCAR' in kwargs:
            del kwargs['CHGCAR']
        # print(self.path['charge'])
        # return self.get_file(filetype = str(self.id[2])+'.CHGCAR', **kwargs)
        return self.get_file(filetype = 'CHGCAR', **kwargs)




    def bader(self):
        chgcar = self.path['chgcar']
        acf = self.dir+'/ACF.dat'
        printlog('Bader should be installed', imp = 'Y')
        if not os.path.isfile(acf):
            cwd = os.getcwd()
            os.chdir(self.dir)
            print(runBash('bader ' + os.path.basename(chgcar) ) )
            os.chdir(cwd)

        else:
            charges = []
            with open(acf, 'r') as f:
                for line in f:
                    try:
                        charges.append(round(float(line.split()[4]), 3))
                    except:
                        pass
            print(dict(zip(charges, self.end.get_elements())))



    def get_bader_ACF(self, p = 0):
        #Make bader on server
        #assumes that bader is installed
        from siman.small_functions import bash_chk_file_cmd
        self.res()
        v = str(self.version)
        path = self.project_path_cluster+'/'+self.dir
        ppc = self.project_path_cluster+'/'
        self.determine_filenames()

        # print()
        CHG_scratch_gz  = header.PATH2ARCHIVE+'/'+self.dir+'/'+v+".CHGCAR.gz"

        CHG     = ppc + self.path['chgcar']
        AECCAR0 = ppc + self.path['aeccar0']
        AECCAR2 = ppc + self.path['aeccar2']
        CHGCAR_sum = path+v+".CHGCAR_sum"
        baderlog =  path+v+".bader.log"
        ACF      = path+v+'.ACF.dat'
        self.path['acf'] = ACF



        run_chgsum = "cd "+path+"; ~/tools/vts/chgsum.pl "+AECCAR0+" "+AECCAR2+"; "+\
        "mv CHGCAR_sum "+CHGCAR_sum+";" # on cluster

        command_chg_gunzip = 'gunzip '+CHG+'.gz ' # on cluster
        
        restore_CHG = "rsync "+CHG_scratch_gz+' '+path+' ; gunzip '+CHG+'.gz ' # on cluster
        restore_AEC = "rsync "+header.PATH2ARCHIVE+'/'+self.path['aeccar0']+' '+ header.PATH2ARCHIVE+'/'+self.path['aeccar2']+' '+path # on cluster


        mv = v+".bader.log; mv ACF.dat "+v+".ACF.dat; mv AVF.dat "+v+".AVF.dat; mv BCF.dat "+v+".BCF.dat; cat "+v+".bader.log"
        
        bader_on_sum = "cd "+path+"; ~/tools/bader "+CHG+" -ref "+CHGCAR_sum+" > "+mv
        bader_on_chg = "cd "+path+"; ~/tools/bader "+CHG+" > "+mv #simple
        command_cat_ACF    = " cat "+path+v+".ACF.dat"
        
        no_CHG_sum = bash_chk_file_cmd(CHGCAR_sum) #true if file not exists
        no_CHG     = bash_chk_file_cmd(CHG)
        no_ACF     = bash_chk_file_cmd(ACF)
        no_AECCAR0 = bash_chk_file_cmd(AECCAR0)
        no_AECCAR2 = bash_chk_file_cmd(AECCAR2)
        
        def remote(cmd):
            return run_on_server(cmd, self.cluster['address'])



        #Calculate CHGCAR_sum
        if remote(no_CHG_sum): 
            printlog(  CHGCAR_sum, "does not exist. trying to calculate it ...", imp = 'Y')
            
            if remote(no_AECCAR0) or remote(no_AECCAR2):
                printlog(  AECCAR0, "does not exist, trying to take it from archive ...", imp = 'Y')
                printlog(remote(restore_AEC)+'\n', imp = 'y')

            if not remote(no_AECCAR0):
                printlog( remote(run_chgsum)+'\n', imp = 'Y' ) 
                printlog( remote(" rm "+path+v+".ACF.dat")+'\n', imp = 'Y' ) 
        # sys.exit()

        #Check chgcar
        if remote(no_CHG): #true if file not exists
            printlog( 'Warning! File ', CHG, "does not exist. Checking .gz .. ", imp = 'Y')

            printlog( remote(command_chg_gunzip)+'\n', imp = 'y' ) 

        if remote(no_CHG): #true if file not exists
            printlog( 'Warning! File ', CHG, "does not exist. Trying to restore it from archive .. ", imp = 'Y')

            printlog( remote(restore_CHG)+'\n', imp = 'y' ) 



        def run_bader(command, ):        
            
            if remote(no_ACF): #true if file not exists
                printlog(  ACF, " does not exist. trying to calculate Bader ... ", imp = 'Y')
                printlog( remote(command)+'\n', imp = 'y' ) 
                
            ACF_text = remote(command_cat_ACF)

            return ACF_text

        ACF_text = run_bader(bader_on_sum)

        if 'No such file or directory' in ACF_text:
            printlog('Warning! Probably you have problems with',CHGCAR_sum)
            printlog('Trying to calculate charges from CHGCAR ...', imp = 'Y')
            
            ACF_text = run_bader(bader_on_chg)

        print('ACF_text = ', ACF_text)


        if ACF_text:
            ACF_l = ACF_text.splitlines()
        charges = []
        

        for line in ACF_l:
            try:
                charges.append(round(float(line.split()[4]), 3))
            except:
                pass
        
        self.charges = charges

        # print(charges[[1,2]])
        if p:
            print(list(zip(charges, self.end.get_elements())))


        return charges
        


    def get_occ_mat(self, i):
        """
        return occupation matrix for atom i (from zero)


        TODO: probably it is better to move occ_matrix to structure class and make this method their
        """
        st = self.end
        # i_tran = st.get_transition_elements('n')
        # print(st.get_elements())
        # print(i_tran[21])
        # i_mag = i_tran.index(i)
        # print(i, i_mag)
        # print( self.occ_matrices )

        return self.occ_matrices.get(i)

    def set_occ_mat(self, i, m):
        """
        set occupation matrix m for atom i (from zero)


        TODO: probably it is better to move occ_matrix to structure class and make this method their
        """
        cl = self.copy()
        st = cl.end
        # i_tran = st.get_transition_elements('n')
        # print(st.get_elements())
        # print(i_tran[21])
        # i_mag = i_tran.index(i)
        # print(i_mag)
        cl.occ_matrices[i] = m

        return cl


    def write_occmatrix(self):

        #write occmatrix file to calculation folder
        from siman.inout import write_occmatrix
        # print(self.get_path())
        # sys.exit()

    

        return write_occmatrix(self.occ_matrices, self.get_path())


    def vasp_dipole_center(self):
        """
        Determine dipole center in reduced coordinates 
        based on vasp definition, https://www.vasp.at/wiki/index.php/DIPOL
        as only number of position of minimum charge density is given
        
        TODO
        works only in the case of dipole along z axis, as min_pos is read only for one direction
        
        """

        if hasattr(self, 'ngxf'):
            dc = return_xred([0,0, self.dipole_min_pos[-1]/self.ngxf[2]+0.5])[2]
        else:
            dc = None

        return dc




    def bader_coseg():

        "used in coseg project Ti- C,O" 
        ACF = self.get_bader_ACF()

        ACF = ACF.splitlines()[2:] #list of lines with charges for each atom

        # print ACF[0]
        if self.end.znucl[1] == 8:
            imp_valence_chg = 6 #!!should be taken from potential
        elif self.end.znucl[1] == 6:
            imp_valence_chg = 4 #!!should be taken from potential
        if self.end.znucl[0] == 22:
            mat_valence_chg = 12 #!!should be taken from potential



        local_atoms = local_surrounding(self.xcart[-1], self.end, 6, control = 'atoms')
        numbers = local_atoms[2] # first atom is impurity
        # print numbers
        imp_partial_chg = imp_valence_chg - float(ACF[numbers[0]].split()[4])

        mat_partial_chg = [mat_valence_chg - float(ACF[i].split()[4]) for i in numbers[1:] ]
        print_and_log( "Partial charge of impurity ", imp_partial_chg, imp = 'Y' )
        print_and_log( "Partial charges of neibouring Ti atoms", " ".join("{:.2f}".format(m) for m in mat_partial_chg), imp = 'Y' )
        print_and_log( "Partial charge of matrix", sum(mat_partial_chg), imp = 'Y' )
        
        print_and_log( "Sum of mat and imp charges:", sum(mat_partial_chg)+imp_partial_chg, imp = 'Y' )

        return path_to_file






    def res(self, **argv):
        from siman.calc_manage import res_loop
        # print(argv)
        # sys.exit()
        res_loop(*self.id, **argv)

    def run(self, ise, iopt = 'full_nomag', up = 'up1', vers = None, i_child = -1, add = 0, it_suffix_del = True, *args, **kwargs):
        """
        Wrapper for add_loop (in development).
        On a first run create new calculation. On a second run will try to read results.
        All children are saved in self.children list.
        By default uses self.end structure
        To overwrite existing calculation 
        use combination of parameters: add = 1, up = 'up2'.
        Allows to use all arguments available for add_loop()


        INPUT:
            ise (str) - name of new set available in header.varset

            iopt (str) - inherit_option
                'full_nomag' - full without magmom
                'full' - full with magmom
                'full_chg' - full with magmom and including chg file
            
            up (str) - update key transferred to add_loop and res_loop;
                'up1' - create new calculation if not exist
                'up2' - recreate new calculation overwriting old; for reading results redownload output files

            vers (list of int) - list of version for which the inheritance is done

            i_child (int) - choose number of child in self.children to run res_loop(); can be relevant if more than one
                calculation exists for the same set
            
            add (bool) - 
                1 - overwrite existing children

            it_suffix_del (bool) - needed to be false to use it_suffix with run. Provides compatibility with old behaviour; should be improved


        RETURN:
            cl (Calculation) - new created calculation 


        TODO:
        1. if ise is not provided continue in the same folder under the same name,
        however, it is not always what is needed, therefore use inherit_xred = continue
        """

        add_flag  = add
        if add:
            up = 'up2'

        from siman.calc_manage import add_loop

        if not iopt:
            iopt = 'full'

        if iopt == 'full_nomag':
            suffix = '.ifn'
        if iopt == 'full':
            suffix = '.if'
        if iopt == 'full_chg':
            suffix = '.ifc'




        if 'it_suffix' in kwargs:
            it_suffix = '.'+kwargs['it_suffix']
        else:
            it_suffix = ''
        
        if it_suffix_del:
            if kwargs.get('it_suffix'):
                del kwargs['it_suffix']


        # if self.id[1] != ise:
        if 1:
            if not hasattr(self, 'children'):
                self.children = []

            if not add and len(self.children)>0:
                print('Children were found in self.children:', len(self.children), ' childs, by default reading last, choose with *i_child* ')
                
                idd = None
                for i in self.children:
                    # print(i, ise, i[1], i[1] == ise)
                    # print(i[0], self.id[0]+it_suffix)
                    if self.id[0]+suffix+it_suffix == i[0] and i[1] == ise:
                        # print(i)
                        idd = i
                        # add = True
                        # break

                if idd is None:
                    add = True
                    # idd  = self.children[i_child]
                # print(idd)
                # sys.exit()

                if idd:
                    # print('setaset')
                    cl_son = header.calc[idd]
                    try:
                        res_params = kwargs['params'].get('res_params') or {}
                    except:
                        res_params = {}

                    cl_son.res(up = up, **res_params, **kwargs)
                    child = idd
                    add = 0
                else:
                    child = None
            

            vp = header.varset[ise].vasp_params
            ICHARG_or = 'impossible_value'

            # print(add, len(self.children) )
            # sys.exit()

            if add or len(self.children) == 0:
                

                if iopt  == 'full_chg':
                    if 'ICHARG' in vp and vp['ICHARG'] != 1:
                        printlog('Warning! Inheritance of CHGCAR and ICHARG == 0; I change locally ICHARG to 1')
                        ICHARG_or = vp['ICHARG']
                        vp['ICHARG'] = 1
                    


                if not vers:
                    vers = [self.id[2]]

                idd = self.id
                it_new = add_loop(idd[0],idd[1], vers, ise_new = ise, up = up, inherit_option = iopt, override = 1, *args, **kwargs)
                # it_new = add_loop(*self.id, ise_new = ise, up = up, inherit_option = iopt, override = 1)
                child = (it_new, ise, self.id[2])

                if child not in self.children:
                    self.children.append(child)


                if ICHARG_or != 'impossible_value':
                    vp['ICHARG'] = ICHARG_or  #return original value
 

        return header.calc[child]


    def full(self, ise = None, up = 0, fit = 1, suf = '', add_loop_dic  = None, up_res = 'up1'):
        """
        Wrapper for full optimization
        ise (str) - optimization set; if None then choosen from dict
        up (int) - 0 read results if exist, 1 - update
        fit (int) - 1 or 0
        suf - additional suffix
        """
        from siman.project_funcs import optimize
        if ise is None:
            if 'u' in self.id[1]:
                ise = '4uis'
        st = self.end.copy()
        it = self.id[0]
        child = (it+suf+'.su', ise, 100)
        #st.printme()
        if not hasattr(self, 'children'):
            self.children = []
        if not up and child in self.children:
            optimize(st, it+suf, ise = ise, fit = fit, add_loop_dic = add_loop_dic,up_res = up_res) # read results
        else:
            #run
            optimize(st, it+suf, ise = ise, add = 1, add_loop_dic = add_loop_dic, up_res = up_res)
            self.children.append(child)

        return






    def read_pdos_using_phonopy(self, mode = 'pdos', poscar = '', plot = 1, up = 'up1'):
        """
        mode - 
            pdos
            band
            free - thermal properties, converted to eV!!!
        """

        if plot == 1:
            p = ' -p '
        else:
            p = ''

        from siman.calc_manage import create_phonopy_conf_file, read_phonopy_data
        from siman.header import PATH2PHONOPY as phonopy_command
        
        self.get_file('vasprun.xml', nametype = 'asoutcar', up = up)
        create_phonopy_conf_file(self.end, mp = [10, 10, 10], dim = [1, 1, 1], path = self.dir)
        # create_phonopy_conf_file(self.end, mp = [36, 36, 36], path = self.dir) #almost no difference was found for Na2X
        create_phonopy_conf_file(self.end, path = self.dir, filetype = 'band', dim = [1, 1, 1]) #create band file



        cwd = os.getcwd()


        os.chdir(self.dir)
        print(self.dir)
        out = runBash(phonopy_command+' --fc '+os.path.basename(self.path['xml']))

        printlog('phonopy out: ', out)



        if 'poscar' not in self.path:
            self.path['poscar'] = self.path['output'].replace('OUTCAR','POSCAR')

        if not poscar:
            poscar = os.path.basename(self.path['poscar'])

        if mode == 'pdos':
            # print('phonopy -c '+os.path.basename(self.path['poscar'])+p+'  mesh.conf --readfc ')
            # runBash('phonopy -c '+os.path.basename(self.path['poscar'])+p+' mesh.conf --readfc ')
            print(phonopy_command+' -c '+poscar+p+'  mesh.conf --readfc ')
            print(runBash(phonopy_command+' -c '+poscar+p+' mesh.conf --readfc '))

        from siman.calc_manage import read_phonopy_dat_file

        self.pdos = read_phonopy_dat_file('total_dos.dat')


        #phonons
        

        if mode == 'band':
            print(phonopy_command+' -c '+os.path.basename(self.path['poscar'])+' -p band.conf --readfc ')
            runBash(phonopy_command+' -c '+os.path.basename(self.path['poscar'])+' -p band.conf --readfc ')

        if mode == 'free':
            print(phonopy_command+' -c '+os.path.basename(self.path['poscar'])+' -t -p mesh.conf --readfc ')

            runBash(phonopy_command+' -c '+os.path.basename(self.path['poscar'])+' -t' +p+' mesh.conf --readfc ')


            Trange, func = read_phonopy_data('thermal_properties.yaml', convert = 1)

            self.F = func # free energy function in eV, still for the whole supercell!
            # print(self.id, self.F)
            Trange, func = read_phonopy_data('thermal_properties.yaml', key = 'entropy', convert = 1)
            self.entropy = func/1000 # entropy function in eV/K, still for the whole supercell!
            Trange, func = read_phonopy_data('thermal_properties.yaml', key = 'energy', convert = 1)
            self.Uvib = func # internal energy (phonon, including zero-point) function in eV, still for the whole supercell!



        os.chdir(cwd)

        return



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
        


def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__
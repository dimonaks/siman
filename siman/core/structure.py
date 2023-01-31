# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
import os
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
    header.pymatgen_flag = True
except:
    print('structure.py: pymatgen is not available')
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

    def get_elements_by_el_name(self, el_name):
        #return list of at numbers of el_name
        elements = self.get_elements()
        el_list = [n for n,x in enumerate(elements) if x==el_name]
        # el_list = elements.index(el_name)
        return el_list


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


    def determine_symmetry_positions(self, element, silent = 0):
        from siman.geo import determine_symmetry_positions

        return determine_symmetry_positions(self, element, silent)

    def get_symmetry_positions(self, element):
        #just another name
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

    def show_mag(self, i, from_one = None):
        # show magmom of i atoms, starting from 0
        if from_one is None:
            from_one = 0
        if from_one:
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

            pm = pymatgen.core.Structure(self.rprimd, elements, self.xred, site_properties = site_properties)

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


    def get_oxi_states(self, typ = 'charges', silent = 1,):
        """
        Create and return list of oxidation states from charges and valences
        self.charges should exist as full charges (e.g. from Bader analysis)
        
        INPUT:
            - typ (str) 
                - 'charges' - from charges and zval
                - 'guess'   - from guess
            - silent (bool) - no output

        RETURN
            - oxi (list) - list of oxidation states for each atom

        SIDE
            - fill in self.oxi_state property
        """
        st = self
        if typ == 'charges':
            printlog('Using zval as reference', imp = '')
            st = self
            oxi = []
            for j, z_val, el in zip(range(st.natom), st.get_elements_zval(), st.get_elements()):
                # if el == 'Ni':
                    # print(z_val, st.charges[j])
                oxi.append( z_val - st.charges[j] )
            st.oxi_state = oxi
        
        elif typ == 'guess':
            pm = st.convert2pymatgen()
            pm.add_oxidation_state_by_guess()
            st = st.update_from_pymatgen(pm)
            oxi = st.oxi_state


        return oxi

    def print_oxi(self, el, silent = 0, from_one = 0):
        """
        Print oxidation states for particular element 

        INPUT:
            - el (str) - element
            - from_one (bool) - show atom numbers starting from one
        """
        iel = self.get_specific_elements([invert(el)])

        fo = int(from_one)

        if not silent:
            # print(' '.join( ['{:}{:}={:.1f}'.format(el, i, self.oxi_state[i]) for i in iel] ) )
            print(f'Element is {el}')
            print(f'i_at    Oxidation state')
            for i in iel:
                ox = self.oxi_state[i]
                print(f'{i+fo:4}   {ox:+5.1f} ')


        return [self.oxi_state[i] for i in iel]

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


    def make_neutral(self, *args, **kwargs):
        return make_neutral(self, *args, **kwargs)


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


    def get_el(self, el):
        """
        Get atomic numbers for el
        """
        z = invert(el)
        num = self.get_specific_elements([z,])


        return num


    def get_dipole(self, ox_states = None, chg_type = 'ox'):
        """ Return dipole moment in e*A calculated by pymatgen
        ox_states (dict) - oxidation states of elements
        chg_type (str) - type of charges  if provided in self.charges; see description of self.convert2pymatgen()
        
        If you need to convert e*A to debye (D), use 1 D = 0.20819434 eA = 3.33564×10−30 Cm; 1 eA = 4.8 D
        """
        slab = self.convert2pymatgen(slab = 1, oxidation = ox_states, chg_type = chg_type)
        return slab.dipole


    def add_atoms(self, atoms_xcart = None, element = 'Pu', return_ins = False, selective = None, 
        atoms_xred = None, mag = None, add2end = False):
        """
        Appends atoms of type *element* at the end if the element is new. Otherwise insert according to VASP conventions if *add2end* is False.
        
        Updates ntypat, typat, znucl, nznucl, xred, magmom and natom
        
        INPUT: 

            - atoms_xcart (list of ndarray)
            
            - atoms_xred (list of coordinate lists) - if provided both, this has higher priority 

            - selective (list of lists) - selective dynamics

            - mag magnetic moment of added atoms, if None, than 0.6 is used
                magmom is appended with 0.6, 
                please improve me! by using the corresponding list of magmoms
            
            - add2end (bool) - override default behavior by appending all atoms to the end


        RETURN:

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

        if hasattr(st, 'magmom') and any(st.magmom) or (mag and st.natom == 0) :
            magmom_flag = True
        else:
            magmom_flag = False

        # print(magmom_flag, st.magmom)
        # sys.exit()

        if mag is None:
            mag = 0.6

        if el_z_to_add not in st.znucl or add2end:
            
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


    def add_atom(self, xr = None, element = 'Pu', xc = None, selective = None, add2end = False):
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

        st = self.add_atoms([xc], element = element, selective = selective, add2end = add2end)
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
        # print(st.magmom, i)
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


        if mode in [1,2]:
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
                            # print(st.magmom)
                            if hasattr(st, 'magmom') and any(st.magmom):

                                st.magmom[n] = mag_new
                            # sys.exit()
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
        if mode == 3:
            xcart_replace = [st.xcart[i] for i in atoms_to_replace]
            st = st.remove_atoms(atoms_to_replace)
            st = st.add_atoms(xcart_replace, element = el_new, mag = mag_new)


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





    def combine(self, st_list, only_numbers = None, add2end = False):
        """
        Combine several structures into one
        using reduced coordinates

        INPUT:
            - st_list (list) - list of structures to combine
            - only_numbers (list) - list of atom numbers to be added; if None all atoms are added
            - add2end (bool) - see self.add_atoms()

        """
        st_b = self.copy()

        # if only_numbers is None:
            # only_numbers = []

        for i, st in enumerate(st_list):
            # print(i)

            for j, xr, el in zip(list(range(st.natom)), st.xred, st.get_elements() ):

                if only_numbers is None or j in only_numbers:
                    # print(j)
                    st_b = st_b.add_atom(xr, el, add2end = add2end)


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
        # print('vlength', self.vlength)
        vmax = max(self.vlen)
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
        Shortest distance between two atoms acounting PBC, numbers from 0
        i1 and i2 override x1 and x2

        coord_type - only when x1 and x2 are provided

        """
        # print(self.xcart)
        if i1 is not None:
            x1 = self.xcart[i1]
        if i2 is not None:
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


    def remove_closest(self, *args, **kwargs):
        return remove_closest(self, *args, **kwargs)

    def remove_vacuum(self, *args, **kwargs):
        return remove_vacuum(self, *args, **kwargs)

    def move_edge(self, *args, **kwargs):
        return move_edge(self, *args, **kwargs)

    def find_slab_width(self, *args, **kwargs):
        return find_slab_width(self, *args, **kwargs)

    def find_closest_atom(self, xc = None, xr = None):
        """
        Find closest atom in structure to xc (cartesian) or xr (reduced) coordinate

        TODO: Should be revised!

        RETURN:
            - i (int)
            - shifts ()
            - dist
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


    def rms(self, st, el = None ):
        """
        Get rms difference 

        INPUT:
            - st (Structure) - structure to be compared with
            - el (str) - see rms_between_structures2()
        RETURN:


        """
        rms_between_structures2(self, st, el)

        return




    def nn(self, i, n = 6, ndict = None, only = None, silent = 0, 
        from_one = None, more_info = 0, oxi_state = 0, print_average = 0):
        """
        show neigbours

        INPUT:
        i - number of central atom, from 1 or 0 (from_one = True or False); if header.from_one is False or True than the behaviour is overriden
        n - number of neigbours to return
        ndict (dic) - number of specific neigbour atoms to take into account e.g ndict = {8:3} - 3 oxygen atoms will be considered
        only - list of interesting z neighbours

        more_info - return more output - takes time

        from_one - if True, strart first atom from 1, otherwise from 0

        oxi_state (bool) - if 1 then showing oxidation state as well, turned on automatically if self.oxi_state list is non zero

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
        
        if from_one is None:
            from_one == True # default due to compatability with previous code

            if header.FROM_ONE is not None: #overwrites default behavior
                from_one = header.FROM_ONE

        if header.FROM_ONE != from_one:
            printlog('Warning! provided *from_one* and header.FROM_ONE are different. I am using from header')



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

        if hasattr(st, 'oxi_state') and len(st.oxi_state) > 0:
            oxi_state = 1

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

        print

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
            print('av(A-O,F) = {:5.2f} A'.format( info['av(A-O,F)']))

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
            # print('center', center)
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

    def localize_ox_polaron(self, i, d_a, d_tm, nn = 6, mag = None):
        """
        Localize small polaron at oxygen atom by adjusting distances
        i - number of oxygen, from 0
        d_a - shift in angstrom; positive increase Li-O bonds, negative reduce Li-O bonds 
        d_tm - shift in angstrom; positive increase TM-O bonds, negative reduce TM-O bonds 
        nn - number of neigbours

        """
        # nn

        st = copy.deepcopy(self)
        i_ox = st.get_el_z(i)
        i_name = st.get_el_name(i)
        if i_name != 'O':
            printlog('Warning! provided element ', i_name, 'is not an oxygen.. Choose another mode. ')
            return
        else:

            silent = 1
            if 'n' in header.warnings or 'e' in header.warnings:
                silent = 0
            # silent = 0


            dic = st.nn(i, nn, from_one = 0, silent = silent)
            printlog('Average distances around O before localization is {:.2f}'.format(dic['av(A-O,F)']), imp = '')

            #updated xcart
            xc = st.xcart[i]
            for j, x in zip(dic['numbers'][1:], dic['xcart'][1:]):
                TM = st.get_el_z(j)
                if TM in header.TRANSITION_ELEMENTS:
                    d = d_tm
                else:
                    d = d_a 
                x1 = st.xcart[j]
                v = xc-x
                vn = np.linalg.norm(v)
                mul = d/vn
                dv = v * mul
                st.xcart[j] = st.xcart[j] -  dv 

            st.update_xred()

            dic = st.nn(i, nn, from_one = 0, silent = silent)
            printlog('Average distances around O after localization is {:.2f}'.format(dic['av(A-O,F)']), imp = '')

            st.name+='pol'+str(i+1)
            if mag:
                st.magmom[i] = mag

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

    def get_coordination(self, el, silent = 1):
        """
        Get coordination for atom *el* for all non-equivalent types

        INPUT:

            - el (str) - element name


        RETURN:

            - coordination_list (list of lists of str)

        AUTHOR:

            Aksyonov DA
        """        

        st = self
        numbers = st.determine_symmetry_positions(el, silent = 1)
        i = 1
        els = st.get_elements()
        coords = []
        for nn in numbers:
            elsc = [els[i] for i in st.nn(nn[0], from_one = 0, silent =1)['numbers'][1:]]
            if not silent:
                printlog(f'Coordination of {el}{i} is {elsc}', imp = 'y')
            i+=1
            coords.append(elsc)

        return coords


    def add_types_for_el(self, el):
        """
        Make several types for the same element. Needed to set different U and magnetic moments
        Currently new types are created based on symmetry. New types are added at the end of typat list

        INPUT:

            - el (str) - element name for which new types should be created

        RETURN:

            - st (Structure)
            - coords (list of lists) coordination of detected types of element *el*

        TODO:
            other methods may not work correctly with such structures. Check their 


        AUTHOR:

            Aksyonov DA
        """

        st = self.copy()

        numbers = st.determine_symmetry_positions(el, silent = 1)
        ntn = len(numbers) # new number of types for el
        z = invert(el)
        ntc = st.znucl.count(z) #current number of types for el
        # print(st.ntypat, st.typat)# st.typat, st.znucl)
        # print(numbers)
        nt_add = ntn-ntc
        # coords = None
        coords = self.get_coordination(el, silent = 0)
        if nt_add > 0:
            printlog('Current number of types=', ntc, 'is smaller than the number of non-equivalent positions=', ntn, '; Additional', nt_add, 'types will be added', imp = 'y')
            
            # st.ntypat += nt_add
            nt_last = st.ntypat 


            for nn in numbers[1:]: # skip first
                st.ntypat += 1
                st.znucl.append(z)
                for i in nn:
                    st.typat[i] = st.ntypat
            # print
            # print(st.ntypat, st.typat)# st.typat, st.znucl)
            st.get_nznucl()

            # st.printme()
            # st.write_poscar()

        else:
            printlog('Current number of types=', ntc, 'is equal or larger than number of non-equivalent positions=', ntn, 'nothing is done', imp = 'y')

        return st, coords



    def get_unique_type_els(self, coordination = False):
        """
        Get list of unique type elements 

        INPUT:

            - coordination (bool) - if true than the coordination of the element is given in format 'el/el_coord', where
            el_cood is the closest-lying coordinating element. It is done only if more than one type for one element is found

        RETURN:

            - els (list of str)

        AUTHOR:

            Aksyonov DA

        """

        els_type = [invert(z) for z in self.znucl] # 

        els = list(set(els_type))

        els_typen = {}

        if coordination:
            els_typec = els_type.copy()
            for i, el in enumerate(els_type):
                if els_type.count(el) > 1:
                    if el not in els_typen:
                        els_typen[el] = 0
                    coords = self.get_coordination(el, silent = 1)
                    coord_el = coords[els_typen[el]][0] #currently only first element is used
                    els_typec[i] = els_type[i]+f'/{coord_el}'
                    els_typen[el] +=1



        return els_typec



    def ewald(self, ox_st = None, site = None):
        # ox_st 
        #   # 1 - oxidation states from guess
            # 2 - from potential
            # None - from charge
        # site if provided (from 0), than site energy is printed
        from pymatgen.analysis.ewald import EwaldSummation
        # from siman.analysis import set_oxidation_states

        st = copy.deepcopy(self)

        if ox_st == 1:
            # st = set_oxidation_states(st)
            #st.printme()
            stpm = st.convert2pymatgen(chg_type = 'pm')
            # print('The following oxi states were set', st.oxi_state)
        elif ox_st == 2:
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
            4 - open mcif to see magnetic moments only on oxygen, other magmoms are set as zero
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
        elif r == 4:
            for i in range(0, len(st.xcart)):
                if st.get_el_name(i) != 'O':
                    st.magmom[i] = 0
                else:
                    st.magmom[i] *= 5
            filename = st.write_cif(mcif = 1)
        else:
            filename = st.write_poscar(vasp5 = 1)
        
        # print(r, filename)
        # sys.exit()
        if 'jmol' in program :
            runBash(header.PATH2JMOL+' -j \"background white\" '+filename, detached = True)
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


    def run_vasp(self):

        """
        Convinient wrapper for add_loop()

        TODO: to be finished; It is more reasonable to move this method to class; Then a new calculation is created which requires a structure
        and set; But probably it can be here as well. It will be very simple to use. 

        please make it in such a way that a res command is suggested at the end, required to read the results
        """
        return

# -*- coding: utf-8 -*- 
#Copyright Aksyonov D.A
from __future__ import division, unicode_literals, absolute_import, print_function
import itertools, os, copy, math, glob, re, shutil, sys, pickle, gzip, shutil
import re, io, json
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
from siman.small_functions import makedir, angle, is_string_like, cat_files, grep_file, red_prec, list2string, is_list_like, b2s, calc_ngkpt, setting_sshpass
from siman.functions import (read_vectors, read_list, words,
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

    def get_mag_tran(self, to_ox = None):
        #show formatted mag moments of transition metals 
        #to_ox - convert to oxidation state, substract from to_ox
        # if to_ox is negative, then m-to_ox
        l, mag_numbers = self.get_maglist()

        keys = list(mag_numbers.keys())#[0]
        print('In the present cell are next TM types:', keys)

        mag = list(np.array(self.magmom)[l])

        for key in keys:
        
            # print(mag)
            magnetic = mag[:len(mag_numbers[key])]
            mag = mag[len(mag_numbers[key]):]
            # print(magnetic)

            s = ' '.join(['{:5.2f} ']*len(magnetic))
            
            s0 = ' '.join(['{:5d} ']*len(magnetic))

            # print(*mag_numbers[key])

            print('\n Znucl:  ', key)
            # print(' '+s0.format(*mag_numbers[key]))
            print(' '+s.format(*magnetic))
            if to_ox:
                if to_ox > 0:
                    ox = [to_ox-abs(m) for m in magnetic]
                else:
                    ox = [abs(m)+to_ox for m in magnetic]
                s2 = ' '.join(['{:5.1f}+']*len(magnetic))
                print(s2.format(*ox))
                print('Average {:5.1f}+'.format(sum(ox)/len(ox)))


        return s.format(*magnetic) 

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
            if hasattr(self, 'oxi_state') and any(self.oxi_state):
                #simply use predefined oxi_state
                
                oxi = self.oxi_state


                if 0: #check total charge
                    # st = st.copy()
                    chg = copy.copy(chg)
                    tot = sum(chg)
                    print('Total charge is ', tot, 'I subtract it uniformly from each atom')
                    d = tot/self.natom
                    chg = [c - d for c in chg]




        if oxi:
            pm.add_oxidation_state_by_site(oxi)




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
        st.rprimd = [np.array(vec) for vec in stpm._lattice._matrix]
        # for site in stpm._sites:
            # print(dir(site))
        st.xred   = [np.array(site._frac_coords) for site in stpm._sites]
        

        # print(elements)
        st.update_xcart()

        s = stpm._sites[0]

        if 'magmom' in s.properties:
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

        st.natom = len(st.typat)
        # sys.exit()

        if st.natom != len(st.xred):
            printlog('Error! number of atoms was changed, please improve this method')

        st.name+='_from_pmg'
        return st


    def rotate(self, axis, angle):
        #axis - list of 3 elements, [0,0,1]
        #angle in degrees
        
        from pymatgen.transformations.standard_transformations import RotationTransformation
        
        st = copy.deepcopy(self)
        rot = RotationTransformation(axis, angle)
        stpm = st.convert2pymatgen()
        stpmr1 = rot.apply_transformation(stpm)
        st_r1 = st.update_from_pymatgen(stpmr1)
        st_r1 = st_r1.return_atoms_to_cell()
        return st_r1


    def invert_axis(self, axis):
        #invert one vector
        st = copy.deepcopy(self)

        st.rprimd[axis] *= -1
        st.update_xred()
        st = st.return_atoms_to_cell()
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




    def get_primitive_cell(self):
        """
        return primitive cell 
        """
        st_mp = self.convert2pymatgen()

        sf = SpacegroupAnalyzer(st_mp, ) #symprec = 0.1

        sc = sf.get_primitive_standard_structure() # magmom are set to None

        st = self.update_from_pymatgen(sc)


        return st

    def get_surface_pos(self, reduced = False):
        #allows to return positions of top and bottom surfaces (edge atoms) in cartesian
        #assumed normal to R3
        #small number is added or subtracted to/from edge atom to overcome further numericall errors
        # reduced - reduced coordinations

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


    def get_surface_atoms(self, element, surface = 0, surface_width = 0.5 ):
        #return numbers of surface atoms
        #elememt - which element is interesting?
        #surface_width - which atoms to consider as surface 
        #surface (int) - 0 or 1 - one of two surfaces; surface 1 has lowest z coordinate
        st = self
        surface_atoms = [[],[]]
        
        z = st.get_surface_pos()


        els = st.get_elements()


        for i, x in enumerate(st.xcart):
            el = els[i]
            if el == element:
                # print(x[2])
                if z[0] <= x[2] < z[0]+surface_width:
                    surface_atoms[0].append(i)

                if z[1] - surface_width  < x[2] <= z[1]:
                    surface_atoms[1].append(i)


        return surface_atoms[surface]


    def get_surface_area(self):
        #currently should be normal to rprim[0] and rprim[1]
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
        for i, e, xc, xr in zip( range(self.natom), el, self.xcart, self.xred ):
            Z = invert(e)
            if Z in required_elements and additional_condition(xr[2]):
                tra.append(e)
                ns.append(i)
                xcart.append(xc)
        
        if fmt == 'z':
            tra = [invert(t) for t in tra]
        elif fmt == 'n':
            tra = ns
        elif fmt == 'x':
            tra = xcart

        return tra



    def get_dipole(self, ox_states = None, chg_type = 'ox'):
        #return dipole moment, e*A
        #
        # if you need to convert it to debye (D), use 1 D = 0.20819434 eÅ = 3.33564×10−30; 1 eA = 4.8 D

        slab = self.convert2pymatgen(slab = 1, oxidation = ox_states, chg_type = chg_type)
        return slab.dipole




    def add_atoms(self, atoms_xcart, element = 'Pu', return_ins = False, selective = None, atoms_xred = None):
        """
        appends at the end if element is new. Other case insertered according to VASP conventions
        Updates ntypat, typat, znucl, nznucl, xred, magmom and natom
        atoms_xcart (list of ndarray)
        atoms_xred (list of coordinate lists) - if provided both, this has higher priority (not implemented yet, use xcart!!!)

        selective (list of lists) - selective dynamics

        magmom is appended with 0.6, please improve me! by using other values for magnetic elements 


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

        natom_to_add = len(atoms_xcart)
        if natom_to_add == 0:
            printlog('self.add_atoms(): Nothing to add, continiue')
            return st


        st.natom+=natom_to_add

        el_z_to_add = element_name_inv(element)

        if hasattr(st, 'magmom') and any(st.magmom):
            magmom_flag = True
        else:
            magmom_flag = False


        if el_z_to_add not in st.znucl:
            
            st.znucl.append( el_z_to_add )
            
            st.nznucl.append(natom_to_add)            

            st.ntypat+=1

            typ = max(st.typat)+1
        
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
                st.magmom.extend( [0.6]*natom_to_add  )

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
                st.magmom[j_ins:j_ins] =  [0.6]*natom_to_add


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






    def replace_atoms(self, atoms_to_replace, el_new):
        """
        atoms_to_replace - list of atom numbers starting from 0
        el_new - new element periodic table short name
        """
        st = copy.deepcopy(self)

        numbers = list(range(st.natom))


        atom_exsist = True

        while atom_exsist:


            for i, (n, el) in enumerate(  zip(numbers, st.get_elements()) ):
                # print(i)

                if n in atoms_to_replace:
                    xcart = st.xcart[i]
                    print('replace_atoms(): atom', i, st.get_elements()[i], 'replaced with', el_new)
                    st = st.del_atom(i)

                    st = st.add_atoms([xcart], element = el_new)
                    # print(st.natom)
                    # print(st.get_elements())

                    # print(st.natom)

                    # print(st.get_elements())
                    del numbers[i]

                    break
            else:
                atom_exsist = False
        # printlog('remove_atoms(): Atoms', atoms_to_remove, 'were removed')

        # print(st.get_elements())
        return st


    def replace_atoms2(self, el_old, el_new, concentration):
        """
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




    def add_vacuum(self, vec, thick):
        """
        To be improved
        vector - along which vector, 0, 1, 2
        thick - thickness of vector 


        TODO:
        make thick to be thickness along axis and not vector

        """
        st = copy.deepcopy(self)
        v = st.rprimd[vec]
        v_l = np.linalg.norm(v)
        new_len = v_l+thick

        st.rprimd[vec]*=new_len/v_l

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



    def find_atom_num_by_xcart(self, x_tar, prec = 1e-6):
        """take into account periodic conditions

        TODO:
        make normal function that treats periodic boundary conditions normally!!!
        """

        [xr_tar] = xcart2xred([x_tar], self.rprimd)
        printlog('find_atom_num_by_xcart(): xr_tar = ', xr_tar)
        #PBC!!!
        for i in [0,1,2]:
            if xr_tar[i] < 0:
                xr_tar[i]+= 1
                
            if xr_tar[i] >= 1:
                xr_tar[i]-= 1
        printlog('find_atom_num_by_xcart(): xr_tar after periodic = ', xr_tar)

        # print(xr_tar)
        [x_tar] = xred2xcart([xr_tar], self.rprimd)
        
        printlog('find_atom_num_by_xcart(): x_tar after periodic = ', x_tar)

        # print(x_tar)
        self = self.return_atoms_to_cell()

        for i, x in enumerate(self.xcart):
            if np.linalg.norm(x-x_tar) < prec:
            # if all(x == x_tar):
                printlog('Atom', i+1, 'corresponds to', x_tar)

                return i
        else:
            printlog('Attention, atom ', x_tar, 'was not found' )
        # print self.xcart.index(x_tar)
        # return self.xcart.index(x_tar)


        # print type(atoms_xcart)


    # def 



    def shift_atoms(self, vector_red = None, vector_cart = None):
        """
        Shift all atoms according to *vector_red*
        """
        st = copy.deepcopy(self)
        if vector_cart is not None:
            vec_cart = np.array(vector_cart)
            for xc in st.xcart:
                xc+=vec_cart
            st.update_xred()
            
        else:
            vec = np.array(vector_red)
            for xr in st.xred:
                xr+=vec
            st.xred2xcart()

        

        st = st.return_atoms_to_cell()
        return st



    def replic(self, *args, **kwargs):

        return replic(self, *args, **kwargs)

    def image_distance(self, *args, **kwargs):

        return image_distance(*args, **kwargs)


    def distance(self, i1, i2):
        """
        Shortest distance between two atoms acounting PBC, from 0
        """
        # print(self.xcart)
        x1 = self.xcart[i1]
        x2 = self.xcart[i2]
        return image_distance(x1, x2, self.rprimd)[0]

    def remove_close_lying(self, rm_both = 0, rm_first = 0):
        """
        rm_both (bool) - if True than remove both atoms, if False than only the second atom is removed
        rm_first (bool) - if True than the first atom of two overlapping is removed, otherwise the second atom is removed

        PBC is realized through image_distance

        SIDE:
            write _removed field to returned st

        TODO:
            Works incorrectly for more than one overlap !!!


        """
        st = copy.deepcopy(self)    
        tol = 0.4
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


    def find_closest_atom(self,x):
        #find closest atom in structure to x cartesian coordinate
        #return i and dist
        # for ixs in self.xcart:
        x = np.asarray(x)
        abs_shifts = [np.linalg.norm(x-x1) for x1 in self.xcart]
        # print(sorted(abs_shifts))
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
        info['av(A-O,F)'] = local_surrounding2(x, st, n, 'av', True, only_elements = [8,9], round_flag = 0)
        
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

    def tm_o_distance(self, criteria = 0.03):
        #return average TM-O distance in the cell and list of outstanding bonds 
        #criteria - value in % of average bond length when bond is outstanding

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
                print('Yan-Teller effect is found\n Average min TM-O length is %s \n Average max TM-O length is %s \n'%(mind, maxd) )


        if not k: print('Ok! None outstanding bonds found\n')

        return



    def center(self):
        #return cartesian center of the cell
        return np.sum(self.xcart, 0)/self.natom

    def center_on(self, i):
        #calc vector which alows to make particular atom in the center 
        x_r = self.xred[i]
        center = np.sum(self.xred, 0)/self.natom
        # print(center)
        # sys.exit()
        # print(x_r)
        dv = center - x_r
        # print(dv)
        # print(dv+x_r)
        return dv







    def localize_polaron(self, i, d):
        """
        Localize small polaron at transition metal by adjusting TM-O distances
        i - number of transition atom, from 0
        d - shift in angstrom; positive increase TM-O, negative reduce TM-O
        """
        st = copy.deepcopy(self)
        TM = st.get_el_z(i)
        TM_name = st.get_el_name(i)
        if TM not in header.TRANSITION_ELEMENTS:
            printlog('Warning! provided element ', TM_name, 'is not a transition metal, I hope you know what you are doing. ')

        silent = 1
        if 'n' in header.warnings or 'e' in header.warnings:
            silent = 0
        # silent = 0


        dic = st.nn(i, 6, from_one = 0, silent = silent)
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

        dic = st.nn(i, 6, from_one = 0, silent = silent)
        printlog('Average TM-O distance after localization is {:.2f}'.format(dic['av(A-O,F)']), imp = '')

        st.name+='pol'+str(i+1)

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
            filename = ('xyz/POSCAR_'+st.name).replace('.', '_')

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
                print_and_log("Warning! may be obsolete!!! and incorrect", imp = 'Y')
                f.write("Cartesian\n")
                for xcart in zxcart:
                    for x in xcart:
                        f.write(str(x[0]*to_ang)+" "+str(x[1]*to_ang)+" "+str(x[2]*to_ang))
                        f.write("\n")

                
                if hasattr(self.init, 'vel'):
                    print_and_log("I write to POSCAR velocity as well")
                    f.write("Cartesian\n")
                    for v in self.init.vel:
                        f.write( '  {:18.16f}  {:18.16f}  {:18.16f}\n'.format(v[0]*to_ang, v[1]*to_ang, v[2]*to_ang) )

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
        

        f.close()
        path = os.getcwd()+'/'+filename
        print_and_log("POSCAR was written to", path, imp = 'y')
        return path



    def write_cif(self, filename = None, mcif = False, symprec = 0.1, write_prim = 0):
        """
        Find primitive cell and write it in cif format
        

        mcif (bool) - if True, than write mcif file with magnetic moments included, primitive cell is not supported


        

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


    def jmol(self, shift = None, r = 0, show_voids = False):
        """open structure in Jmol
        
        INPUT:
        shift (list) - shift vector  in reduced coordinates
        r (int ) - parameter
            0 - open POSCAR
            1 - open OUTCAR to see optimization steps
            2 - open mcif to see magnetic moments
            3 - xyz
        show_voids (bool) - replace voids (z = 300) with Po to visualize them
        
        """
        st = copy.deepcopy(self)
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
        runBash(header.PATH2JMOL+' '+filename, detached = True)
        return

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
            self.id = ('0','0',0)
    
    def get_path(self,):
        print( os.path.dirname(os.getcwd()+'/'+self.path['output']))


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
            
            if self.xred == [None]:
                print_and_log("Convert xcart to xred")
                self.xred = xcart2xred(self.xcart, self.rprimd)
            
            if self.xcart == [None]:
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
            # print vel
            if vel[0] != None: 
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
                    f.write("{:d} {:d} {:d}\n".format(v[0], v[1], v[2])  )

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
        self.end.get_mag_tran(*args, **kwargs)


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
            # print(dir(self.init))
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

        if 'LDAU' in vp and vp['LDAU']: 
            # print(vp['LDAU'])

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
            if 'uniform_scale' in self.calc_method:
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
                    f.write('rm PROCAR DOSCAR OSZICAR PCDAT REPORT XDATCAR vasprun.xml\n')
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




            if params and 'charge' in params:
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
            'up1' - do not update
            'up2' - update

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
        address = self.cluster['address']
        if header.override_cluster_address:
            clust = header.CLUSTERS[header.DEFAULT_CLUSTER]
            self.project_path_cluster = clust['homepath']
            address = clust['address']


        path2file_cluster = self.project_path_cluster+'/'+path_to_file

        # print(self.project_path_cluster)
        # sys.exit()

        if os.path.exists(path_to_file) and 'up2' not in up: 
            out = None
        else:
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

    def run(self, ise, iopt = 'full_nomag', up = 'up1', vers = None, i_child = -1, add = 0, *args, **kwargs):
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
                'full_nomag'
                'full'
                'full_chg' - including chg file
            
            up (str) - update key transferred to add_loop and res_loop;
                'up1' - create new calculation if not exist
                'up2' - recreate new calculation overwriting old; for reading results redownload output files

            vers (list of int) - list of version for which the inheritance is done

            i_child (int) - choose number of child in self.children to run res_loop(); can be relevant if more than one
                calculation exists for the same set
            
            add (bool) - 
                1 - overwrite existing children


        RETURN:
            cl (Calculation) - new created calculation 


        TODO:
        1. if ise is not provided continue in the same folder under the same name,
        however, it is not always what is needed, therefore use inherit_xred = continue
        """

        add_flag  = add
        from siman.calc_manage import add_loop

        if not iopt:
            iopt = 'full'



        # if self.id[1] != ise:
        if 1:
            if not hasattr(self, 'children'):
                self.children = []

            if not add and len(self.children)>0:
                print('Children were found:', self.children, 'by defauld reading last, choose with *i_child* ')
                
                idd = None
                for i in self.children:
                    print(i, ise, i[1], i[1] == ise)
                    if i[1] == ise:
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
                else:
                    child = None
            

            vp = header.varset[ise].vasp_params
            ICHARG_or = 'impossible_value'
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

# -*- coding: utf-8 -*- 
#Copyright Aksyonov D.A
from __future__ import division, unicode_literals, absolute_import, print_function
import itertools, os, copy, math, glob, re, shutil, sys, pickle, gzip, shutil
import re, io


#additional packages
try:
    from tabulate import tabulate
except:
    print('tabulate is not avail')
try:
    import pandas as pd
except:
    print('pandas is not avail')

# import pymatgen
# sys.exit()

try:
    import pymatgen
    from pymatgen.io.cif import CifWriter
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.core.surface import Slab
    from pymatgen.core.composition import Composition
    pymatgen_flag = True
except:
    print('pymatgen is not avail')
    pymatgen_flag = False

import numpy as np
# import matplotlib.pyplot as plt


#siman packages
from siman import header

from siman.header import printlog, print_and_log, runBash, plt

from siman.small_functions import makedir, angle, is_string_like, cat_files, grep_file, red_prec, list2string, is_list_like
from siman.functions import (read_vectors, read_list, words,
     element_name_inv, invert, calculate_voronoi,
    get_from_server, push_to_server, run_on_server, smoother, file_exists_on_server)
from siman.inout import write_xyz, write_lammps, read_xyz, read_poscar
from siman.geo import image_distance, replic, calc_recip_vectors, calc_kspacings, xred2xcart, xcart2xred, local_surrounding, determine_symmetry_positions



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




class empty_struct():
    def __init__(self):
        pass

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


    def get_slice(self, thickness):
        #return element numbers from the top part of slab along z
        #slab should start from bottom
        st = self
        z2 = st.get_surface_pos()[1]
        nn = []
        for i, x in enumerate(st.xcart):
            if z2-thickness < x[2] < z2:
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
        """list of numbers of atoms of each type, order is not important
            updated directly
        """
        self.nznucl = []
        for typ in range(1,self.ntypat+1):
            self.nznucl.append(  self.typat.count(typ) )
        return self.nznucl

    def get_elements(self):
        #return list of elements names
        return [element_name_inv(self.znucl[t-1]) for t in self.typat]

    def get_elements_z(self):
        #return list of elements names
        return [self.znucl[t-1] for t in self.typat]

    def get_elements_zval(self):
        #return list with number of valence electrons for each element
        zvals = []
        for z in self.get_elements_z():
            i = self.znucl.index(z)
            zv = self.zval[i]
            zvals.append(zv)
        return zvals

    def determine_symmetry_positions(self, element):
        from siman.geo import determine_symmetry_positions

        return determine_symmetry_positions(self, element)


    def get_maglist(self):
        #return bool list of which  elements are magnetic (here all transition metals are searched!)
        #and dictionary with numbers of each transition metal
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

    def get_mag_tran(self):
        #show mag moments of transition metals 
        l = self.get_maglist()[0]
        # print(l)
        return np.array(self.magmom)[l]

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


        slab - if True return slab object is returned - limited functional is implemented
        """
        from siman.analysis import calc_oxidation_states

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


        if chg_type == 'pot':
            
            printlog('Using zval as charges', imp = '')
            chg = [z*-1 for z in self.get_elements_zval()           ]
            # print(chg)
            pm.add_oxidation_state_by_site(chg)
        else:
            if hasattr(self, 'charges') and any(self.charges):

                chg_type = 'ox' # 'norm', 'tot'

                # print(chg_type)
                if chg_type == 'norm': #normalize charges
                    t = sum(self.charges)/len(self.charges)
                    chg = [c-t for c in self.charges ]
                    # print(t)
                
                elif chg_type == 'ox':
                    chg = calc_oxidation_states(st = self)
                elif chg_type == 'tot':
                    chg = self.charges
                
                if 0: #check total charge
                    # st = st.copy()
                    chg = copy.copy(chg)
                    tot = sum(chg)
                    print('Total charge is ', tot, 'I subtract it uniformly from each atom')
                    d = tot/self.natom
                    chg = [c - d for c in chg]





                pm.add_oxidation_state_by_site(chg)




        return pm







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
        only rprimd, xred and xcart are updated now!!!!!

        TODO:


        """
        st = copy.deepcopy(self)
        st.rprimd = [np.array(vec) for vec in stpm._lattice._matrix]
        st.xred   = [np.array(site._fcoords) for site in stpm._sites]
        

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
            charges = [s.specie.oxi_state for s in stpm._sites]
            st.charges = charges
        
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
        return st_r1


    def invert_axis(self, axis):
        st = copy.deepcopy(self)

        st.rprimd[axis] *= -1
        st.update_xred()
        st = st.return_atoms_to_cell()
        return st











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

    def get_surface_pos(self, ):
        #allows to return positions of top and bottom surfaces (edge atoms) in cartesian
        #assumed normal to R3
        #small number is added or subtracted to/from edge atom to overcome further numericall errors
        st = self

        z1 = 100
        z2 = -100
        z = []
        # print(st.xcart)
        for x in st.xcart:
            if z1 > x[2]:
                z1 = x[2]
            if z2 < x[2]:
                z2 = x[2]

        z1-=0.01
        z2+=0.01
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

        if hasattr(self, 'spg') and symprec == default:
            spg = self.spg
        else:
            p = self.convert2pymatgen()
            spg = p.get_space_group_info(symprec)
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

    def get_specific_elements(self, required_elements = None, fmt = 'n', z_range = None):
        """Returns list of specific elements (chemical names. z, or numbers from 0) in the structure
        required_elements - list of elements z of interest
        z_range - (2 index tuple) range of z coordinates: only atoms from z1 to z2 are taken
        fmt - format of output
            'names'
            'z'
            'n' - numbers of atoms
            'x' - xcart
        


        """
        el = self.get_elements()
        tra = []
        ns = []

        if z_range:
            def additional_condition(x):
                return z_range[0] < x < z_range[1]
        else:
            def additional_condition(x):
                return True

        xcart = []
        for i, e, xc in zip( range(self.natom), el, self.xcart ):
            Z = invert(e)
            if Z in required_elements and additional_condition(xc[2]):
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




    def add_atoms(self, atoms_xcart, element = 'Pu', return_ins = False, selective = None):
        """
        appends at the end if element is new. Other case insertered according to VASP conventions
        Updates ntypat, typat, znucl, nznucl, xred, magmom and natom
        atoms_xcart (list of ndarray)
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
        wrapper 
        allows to add one atom using reduced coordinates or cartesian
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
        order (list) -required order e.g. ['O', 'Li']

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
        t = 1
        for el in order:
            if el not in els:
                printlog('Error! Check *order* list')
            
            znucl.append( invert(el) )

            for i in range(st.natom):
                el_i = els[i]
                if el_i not in order:
                    printlog('Error! Check *order* list')

                if el_i == el:
                    # print(el)
                    typat.append(t)
                    xcart.append(st.xcart[i])
                    magmom.append(st.magmom[i])
            t+=1

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

        printlog('del_atom(): I remove atom ',  st.get_elements()[i], imp = 'n')
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







    def mov_atoms(self, iat = None, to_x = None, relative = False):
        """
        Move one atom to xcart position *to_x*
        relative (bool) - if shift is relative

        """
        st = copy.deepcopy(self)
        
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


    def remove_atoms(self, atoms_to_remove, from_one = 0):
        """
        remove atoms either of types provided in *atoms_to_remove* or having numbers provided in *atoms_to_remove*, starting from 0
        st (Structure)
        atoms_to_remove (list) - list of element names or numbers
        from_one (int)- if 1 the numbers of atoms in provided list are starting from one
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
        xcr = xcart_range
        if xcr:
            xred_range = [xcr[0]/r3, xcr[1]/r3]
            printlog('xcart_range converted to xred', xred_range)

        dels = []
        for i, xr in enumerate(st.xred):
            if xred_range[0]  < xr[2] < xred_range[1]:
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
        el_new - new element
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



    def shift_atoms(self, vector_red):
        """
        Shift all atoms according to *vector_red*
        """
        st = copy.deepcopy(self)
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

    def nn(self, i, n = 6, ndict = None, only = None, silent = 0, from_one = True, more_info = 0):
        """
        show neigbours
        i - number of central atom, from 1 or 0 (from_one = True or False)
        ndict (dic) - number of specific neigbour atoms
        only - list of interesting z neighbours

        more_info - return more output - takes time

        out
            'numbers' from 0 in the new version!!!!!
        """
        import itertools
        from siman.functions import invert
        if from_one:
            i -= 1

        zn = self.znucl
        x = self.xcart[i]
        out_or = local_surrounding(x, self, n, 'atoms', True, only_elements = only)
        # out =  (xcart_local, typat_local, numbers, dlist )

        out = list(out_or)
        # out[0] = list(itertools.chain.from_iterable(out[0]))
        out[1] = [invert(zn[o-1]) for o in out[1]]
        out[2] = [o+1 for o in out[2]]

        out_tab = [range(0, len(out[2])), out[2], out[1], out[3]]

        tab = np.asarray(out_tab).T.tolist()

 
        # df = pd.DataFrame(tab)
        # print(df)
        if not silent:
            print('Neighbors around atom', i+1, self.get_elements()[i],':')
            print( tabulate(tab[1:], headers = ['nn', 'No.', 'El', 'Dist, A'], tablefmt='psql', floatfmt=".2f") )


        info = {}
        info['numbers'] = out_or[2]


        el = self.get_elements()
        info['el'] = [el[i] for i in out_or[2]]
        info['av(A-O,F)'] = local_surrounding(x, self, n, 'av', True, only_elements = [8,9], round_flag = 0)
        
        if more_info:
            info['avsq(A-O,F)'] = local_surrounding(x, self, n, 'avsq', True, only_elements = [8,9])
            info['avdev(A-O,F)'], _   = local_surrounding(x, self, n, 'av_dev', True, only_elements = [8, 9])
            info['sum(A-O,F)'] = local_surrounding(x, self, n, 'sum', True, only_elements = [8,9])

        t = set(out_or[2])
        s = set(range(self.natom)) 
        d = s.difference(t) 
        # d = d.remove(i)
        # print(t)
        # print(i)
        # print(d)
        st_left = self.remove_atoms(d)
        st_left.name+='_loc'
        # sys.exit()
        st_left.dlist = out_or[3] # distances to neighbours
        st_left.ellist = info['el'] # types of neighbours
        info['st'] = st_left

        if ndict:
            info['av(A-O)']   = local_surrounding(x, self, ndict[8], 'av', True, only_elements = [8])
            info['avdev(A-O)'], _   = local_surrounding(x, self, ndict[8], 'av_dev', True, only_elements = [8])
            info['min(A-O)'], _ ,info['max(A-O)']    = local_surrounding(x, self, ndict[8], 'mavm', True, only_elements = [8])
            atoms = local_surrounding(x, self, ndict[8], 'atoms', True, only_elements = [8])
            info['Onumbers'] = atoms[2][1:] # exclude first, because itself!
            # print(info['Onumbers'])

        # print(info)

        return info


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


    def write_poscar(self, filename = None, coord_type = 'dir', vasp5 = True, charges = False, energy = None, selective_dynamics = False):
        """
        write 

        charges (bool) - write charges, self.charges should be available
        energy - write total energy

        selective dynamics - 
            if at least one F is found than automatically switched on
            !works only for coord_type = 'dir' and charges = False

        NOTE
        #void element type is not written to POSCAR

        TODO
            selective_dynamics for coord_type = 'cart'
        """

        def b2s(b):
            #bool to vasp str
            if b:
                s = 'T'
            else:
                s = 'F'

            return s


        st = copy.deepcopy(self)
        st = st.remove_atoms(['void']) # remove voids
 
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


    def jmol(self, shift = None, r = 0):
        # self.write_poscar('CONTCAR', vasp5 = 1)
        #if r == 1 then open outcar
        st = self
        if shift:
            st = st.shift_atoms(shift)


        # filename, _ = st.write_xyz()
        if r == 1:
            filename = st.outfile 
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
        self.set = copy.deepcopy(inset)
        
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
#             #Programm nznucl, since We will have more impurity atoms of different types
#             command="""grep -w -m 1 "natom " """+filename
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
            if hasattr(curset, 'magnetic_moments') and curset.magnetic_moments:
                self.init.magmom = read_list("magmom", self.natom, float, gen_words)
            # self.init.mag_moments 







            self.state = "1.Geometry has been read"



        #file.close();

        print_and_log( "If no warnings, geometry has been succesfully read from file "+filename+" \n")

        return





    def write_geometry(self, geotype = "init", description = "", override = False):
        """Writes geometrical data in custom siman format bases on abinit format to self.path["input_geo"]"""
        geo_dic = {}
        geofile = self.path["input_geo"]
        geo_exists = os.path.exists(geofile)
        # print (os.path.exists(geofile))
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
            

            f.write("version "+str(self.version)+"\n")
            
            try: 
                st.magmom
            except AttributeError:
                st.magmom = [None]
            # print st.magmom 
            # sys.exit()
            if not None in st.magmom:
                f.write("magmom "+' '.join(np.array(st.magmom).astype(str)) +"\n")



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
                f.write("%.12f %.12f %.12f \n"%(v[0], v[1], v[2])  )

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
                f.write("%.12f %.12f %.12f \n"%(v[0], v[1], v[2])  )


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


    def get_kpoints_density(self):
        """
        Number of k-points per atom
        """
        print(self.NKPTS*self.end.natom) #KPPRA - k-points per reciprocal atom? 




    def copy(self):
        return copy.deepcopy(self)

    def jmol(self, *args, **kwargs):
        self.end.jmol(*args, **kwargs)
    def poscar(self):
        self.end.write_poscar()
    def me(self):
        self.end.printme()

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


    @property
    def sfolder(self):
        self._x = header.struct_des[self.id[0]].sfolder
        return self._x




class CalculationAbinit(Calculation):
    """docstring for CalculationAbinit"""
    pass








class CalculationVasp(Calculation):
    """Methods for calculations made using VASP DFT code"""
    def __init__(self, inset = None, iid = None, output = None):
        super(CalculationVasp, self).__init__(inset, iid, output)
        self.len_units = 'Angstrom'



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
            self.init.recip = self.init.get_recip()
            for i in 0, 1, 2:
                N_from_kspacing.append( math.ceil( (np.linalg.norm(self.init.recip[i]) / to_ang_local) / kspacing) )

            N = N_from_kspacing
            printlog('check_kpoints(): k-points are determined from kspacing',kspacing)

        elif self.set.kpoints_file:
            print_and_log("K-points file was provided", self.set.kpoints_file)
            N = None

        else:
            # print(self.dir)
            print_and_log("Error! check_kpoints(): no information about k-points\n")



        self.init.ngkpt = N

        if kspacing != None and kspacing not in ngkpt_dict:
            ngkpt_dict[kspacing] = N
            printlog('check_kpoints(): I added ',N,'as a k-grid for',kspacing,'in struct_des of', it)


        print_and_log("check_kpoints(): Kpoint   mesh is: ", N, imp = 'Y')


        if not hasattr(struct_des[it], 'ngkpt_dict_for_kspacings') or  kspacing not in struct_des[it].ngkpt_dict_for_kspacings:
            print_and_log('Several other options instead of automatically determined ngkpt = ',N,np.array(self.calc_kspacings(N) ).round(2), ':', end = '\n', imp = 'y')
            print_and_log('ngkpt              |    actual kspacings       ', end = '\n', imp = 'y' )
            

            for ngkpt in itertools.product([N[0]-1, N[0], N[0]+1], [N[1]-1, N[1], N[1]+1], [N[2]-1, N[2], N[2]+1]):
                print_and_log(ngkpt, np.array(self.calc_kspacings(ngkpt) ).round(2), end = '\n', imp = 'y' )

            # user_ngkpt = input('Provide ngkpt:')
            # print(user_ngkpt)
            # sys.exit()

        else:
            print_and_log("check_kpoints(): The actual k-spacings are ", np.array(self.calc_kspacings(N) ).round(2), imp = 'Y')
        return N_from_kspacing


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

        if option == 'inherit_xred' and 'car' in type_of_coordinates: raise RuntimeError 

        if option == 'inherit_xred' and prevcalcver: type_of_coordinates = 'None' # do not write xred or xcart if they will be transfered on cluster
        
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


    def calculate_nbands(self, curset, path_to_potcar = None):
        """Should be run after add_potcar()"""
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
        else:
            printlog('Attention! No path_to_potcar! skipping NBANDS calculation')

        return


    def actualize_set(self, curset = None):
        """
        Makes additional processing of set parameters, which also depends on calculation
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

        return

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



    def write_sge_script(self, input_geofile = "header", version = 1, option = None, 
        prevcalcver = None, savefile = None, schedule_system = None,
        output_files_names = None,
        mode = None,
        batch_script_filename = None):
        """Without arguments writes header, else adds sequence of calculatios
            option - the same as inherit_option, 'inherit_xred' - control inheritance, or 'master' - run serial on master 
            prevcalcver - ver of previous calc; for first none
            savefile - 'cdawx', where c-charge, d-dos, a- AECCAR, w-wavefile, x-xml
            schedule_system - type of job scheduling system:'PBS', 'SGE', 'SLURM'
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



        def update_incar(parameter = None, value = None, u_ramp_step = None, write = True):    
            """Modifications of INCAR. Take attention that *parameter* will be changed to new *value*
            if it only already exist in INCAR.  *u_ramp_step*-current step to determine u,
            *write*-sometimes just the return value is needed. 
            Returns U value corresponding to *u_ramp_step*.
            """

            u_step = None
            if parameter == 'LDAUU':
                #Update only non-zero elements of LDAUU with value

                set_LDAUU_list = self.set.vasp_params['LDAUU']
                new_LDAUU_list = copy.deepcopy(set_LDAUU_list)
                
                # print set_LDAUU_list
                u_step = 0.0
                for i, u in enumerate(set_LDAUU_list):
                    if u == 0:
                        continue
                    u_step = np.linspace(0, u, self.set.u_ramping_nstep)[u_ramp_step]
                    u_step = np.round(u_step, 1)
                    # new_LDAUU_list[i] = value
                    new_LDAUU_list[i] = u_step


                new_LDAUU = 'LDAUU = '+' '.join(['{:}']*len(new_LDAUU_list)).format(*new_LDAUU_list)
                
                command = "sed -i.bak '/LDAUU/c\\" + new_LDAUU + "' INCAR\n"
                #print('u_step',u_step)
                #sys.exit()

            elif parameter == 'MAGMOM':

                new_incar_string = parameter + ' = ' + ' '.join(['{:}']*len(value)).format(*value)
                command = "sed -i.bak '/"+parameter+"/c\\" + new_incar_string + "' INCAR\n"

            elif parameter in ['IMAGES', 'ISPIN']:

                new_incar_string = parameter + ' = ' + str(value)
                command = "sed -i.bak '/"+parameter+"/c\\" + new_incar_string + "' INCAR\n"




            if write:
                f.write(command)

            return  u_step #for last element

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
            if write:
                pre = v + name_mod

                contcar = pre+'.CONTCAR'

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
                            u_ramp_step = self.set.u_ramping_nstep-1, write = False)).replace('.','') #used to det last U

            return name_mod_last


        if schedule_system == 'SGE':
            # parrallel_run_command = "mpirun -x PATH vasp" # MPIE
            parrallel_run_command = header.vasp_command
        elif schedule_system in ['PBS', 'PBS_bsu']:
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
                update_incar(parameter = 'IMAGES', value = 0, write = write) # start and final runs

            
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

                    u = update_incar(parameter = 'LDAUU', u_ramp_step = i_u, write = write)
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

                    update_incar(parameter = 'MAGMOM', value = magmom)

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
                    u = update_incar(parameter = 'LDAUU', u_ramp_step = self.set.u_ramping_nstep-1, write = False)
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


                    u = update_incar(parameter = 'LDAUU', u_ramp_step = i_u)
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


                update_incar(parameter = 'IMAGES', value = nim)


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
                    f.write('export PATH=$PATH:'+header.cluster_home+'/tools/gnuplot/bin/ \n')
                    f.write(header.cluster_home+'/tools/vts/nebresults.pl  \n')
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
                set_mod = '' # the last step do no use modifications of names 
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
                f.write('sleep 5\n')                        
            
            elif schedule_system == 'SLURM':
                f.write("squeue\n") 
                f.write("sbatch -p AMG " + run_name+"\n") 
            else:
                printlog('Error! Unknown schedule_system', schedule_system)
                



        printlog("\nRun file created\n")     
        return



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













    def read_results(self, load = '', out_type = '', voronoi = False, show = '', choose_outcar = None, alkali_ion_number = None, only_load = False):

        """
        Download and Read VASP OUTCAR file

        ###INPUT:
            - load (str) - 'x' - download xml, o - download outcar and contcar, un - read unfinished
            - show (str) - print additional information
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


        if choose_outcar and hasattr(self, 'associated_outcars') and self.associated_outcars and len(self.associated_outcars) >= choose_outcar:
            # print ('associated outcars = ',self.associated_outcars)
            printlog('read_results(): choose_outcar', choose_outcar)

            path_to_outcar = join( dirname(self.path["output"]), self.associated_outcars[choose_outcar-1] )

            printlog(self.associated_outcars)
        else:
            path_to_outcar  = self.path["output"]
        

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
                self.cluster_address, join(self.project_path_cluster, path_to_outcar) )
            # runBash(command_reduce)


            if 'un2' in load:
                out_name  = os.path.basename(path_to_outcar)
                cont_name = os.path.basename(path_to_contcar)
                path_to_outcar = path_to_outcar.replace(out_name, 'OUTCAR')
                path_to_contcar = path_to_contcar.replace(cont_name, 'CONTCAR')


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
        if outcar_exist:

            out = grep_file('General timing', path_to_outcar, reverse = True)

            printlog('The grep result of',path_to_outcar, 'is:', out)
            # sys.exit()
            if 'Gen' in out or 'un' in load:
                self.state = '4. Finished'
            else:
                self.state = '5. Broken outcar'



        else:
            self.state = '5. no OUTCAR'
        
        outst = self.state


        if "4" in self.state:
            

            """Try to read xred from CONCAR and calculate xcart"""

            printlog('Path to CONTCAR', path_to_contcar)
            if os.path.exists(path_to_contcar):
                contcar_exist   = True
            else:
                contcar_exist   = False


            printlog('The status of CONTCAR file is', contcar_exist)
            # self.end.update_xred()
        
            if contcar_exist:
                # try:
                self.end = read_poscar(self.end, path_to_contcar, new = False) # read from CONTCAR
                # except:
                contcar_read = True
            else:
                printlog('Attention!, No CONTCAR:', path_to_contcar, '. I use data from outcar')
                contcar_read = False






            read = 1
            if read:
                if 0: #please use this only for linux or create cross-platform way
                    nw = runBash('sed -n "/NPAR = approx SQRT( number of cores)/=" '+path_to_outcar) #remove warinig
                    tmp = path_to_outcar+".tmp"
                    if nw:
                        nw = int(nw)
                        runBash("sed '"+str(nw-11)+","+str(nw+8)+"d' "+path_to_outcar+">"+tmp+";mv "+tmp+" "+path_to_outcar)


                with open(path_to_outcar, 'r') as outcar:
                    
                    printlog("Start reading from "+ path_to_outcar, imp = 'n')
                    outcarlines = outcar.readlines()




                re_lengths = re.compile("length of vectors")
                re_eltime = re.compile("Elapsed time")
                re_nkpts = re.compile("NKPTS")
                iterat = 0
                i_line = 0
                mdstep_prev = 0
                dipol = None
                self.mdstep = 1
                warnings = 0#""
                self.time = 0
                nscflist = []; mdstep_old = 1; niter_old = 0
                maxforce = []; average = [];  gstress =[]
                # mforce = []
                self.list_e_sigma0 = []
                self.list_e_without_entr = []
                # try:
                #     self.end = copy.deepcopy(self.init) # below needed end values will be updated
                # except:
                if not contcar_read:
                    self.end = Structure()

                # if not hasattr(self.end, "natom"): 
                #     self.end.natom = self.natom
                #Structure() #create structure object with end values after calculation
                #self.end.typat = self.typat
                #self.end.znucl = self.znucl
                self.end.name = self.name+'.end'
                self.end.list_xcart = []
                self.energy = empty_struct()

                de_each_md = 0 # to control convergence each md step
                de_each_md_list = []


                nsgroup = None
                magnitudes = []
                self.mag_sum = [] #toatal mag summed by atoms, +augmentation

                tot_mag_by_atoms = [] #magnetic moments by atoms on each step
                tot_chg_by_atoms = []
                tot_mag_by_mag_atoms = []


                ldauu = None
                e_sig0 = 0 #energy sigma 0 every scf iteration
                occ_matrices = {} # the number of atom is the key

                #which kind of forces to use
                if ' CHAIN + TOTAL  (eV/Angst)\n' in outcarlines:
                    force_keyword = 'CHAIN + TOTAL  (eV/Angst)'
                    ff  = (0, 1, 2)
                    force_prefix = ' chain+tot '

                else:
                    force_keyword = 'TOTAL-FORCE'
                    ff  = (3, 4, 5)
                    force_prefix = ' tot '



                # try:
                #     spin_polarized = self.set.spin_polarized # again it will be better to determine this from outcar 
                # except:
                #     spin_polarized = None



                self.potcar_lines = []
                self.stress = None
                self.intstress = None
                spin_polarized = None
                for line in outcarlines:

                    #Check bands

                    # if 'band No.' in line:
                    #     kpoint = float(outcarlines[i_line-1].split()[1])
                    #     lastocc = float(outcarlines[i_line+self.nbands].split()[2])
                    #     lastbandno = outcarlines[i_line+self.nbands].split()[0]
                    #     if lastocc > 0:
                    #         print "Warning!!! at kpoint ", kpoint, " last band No. ",lastbandno, " is not empty ", lastocc

                    if 'TITEL' in line:
                        self.potcar_lines.append( line.split()[2:] )

                    if 'LEXCH  =' in line:
                        # print(line)
                        self.xc_pot = line.split()[2].strip() #xc from potential 

                    if 'GGA     =' in line:
                        # print(line)
                        self.xc_inc = line.split()[2].strip() #xc from incar

                    if 'ions per type =' in line:
                        if not contcar_read:
                            self.end.nznucl = [int(n) for n in line.split()[4:]]
                            self.end.ntypat = len(self.end.nznucl)

                            self.end.natom  = sum(self.end.nznucl)

                            #correction of bug; Take into account that VASP changes typat by sorting impurities of the same type.
                            self.end.typat = []
                            for i, nz in enumerate(self.end.nznucl):
                                for j in range(nz):
                                    self.end.typat.append(i+1)
                            #correction of bug



                            # print(self.potcar_lines)
                            elements = [t[1].split('_')[0] for t in self.potcar_lines]
                            # printlog('I read ',elements, 'from outcar')
                            self.end.znucl = [element_name_inv(el) for el in elements]
                            # print (self.end.znucl)

                        ifmaglist, _ = self.end.get_maglist()


                    if 'ISPIN' in line:
                        if line.split()[2] == '2':
                            spin_polarized = True
                            self.spin_polarized = spin_polarized
                        else:
                            spin_polarized = False
                            self.spin_polarized = False


                    if "TOO FEW BANDS" in line:
                        print_and_log("Warning! TOO FEW BANDS!!!\n\n\nWarning! TOO FEW BANDS!!!\n")



                    #Check W(q)
                    if 'operators is LMAX' in line:
                        lmax = int(line.split()[7])
                        # print 'lmax', lmax
                    if "W(low)/X(q)" in line:
                        kk = 1; 
                        low = []; 
                        high = [];
                        
                        while kk < 100:
                            if 'Optimization' in outcarlines[i_line + kk] or len(outcarlines[i_line + kk].split() ) != 7: 
                                break
                            if 'PSMAXN' in outcarlines[i_line + kk]:
                                # print(line)
                                printlog('Warning! PSMAXN for non-local potential too small')
                                break
                            # print( 'line', outcarlines[i_line + kk])

                            low.append(  float(outcarlines[i_line + kk].split()[4]) )
                            high.append( float(outcarlines[i_line + kk].split()[5]) )
                            kk+=1


                        if any(v > 1e-3 for v in low+high):
                            print_and_log("W(q)/X(q) are too high, check output!\n")
                            print_and_log('Low + high = ', low+high, imp = 'Y' )
                            print_and_log([v > 1e-3 for v in low+high], imp = 'Y' )
                    
                    if "direct lattice vectors" in line:
                        if not contcar_read:
                            for v in 0,1,2:
                                line = outcarlines[i_line+1+v]
                                line = line.replace('-', ' -')
                                # print(line)
                                self.end.rprimd[v] = np.asarray( [float(ri) for ri in line.split()[0:3]   ] )


                        #print self.end.rprimd
                        #print self.rprimd
                    if "POSITION" in line:
                        # if not contcar_exist or out_type == 'xcarts':
                        if not contcar_read or out_type == 'xcarts':
                            local_xcart = []
                            for i in range(self.end.natom):
                                #print outcarlines[i_line+1+i].split()[0:3] 
                                xcart = np.asarray ( 
                                            [   float(x) for x in outcarlines[i_line+2+i].split()[0:3]   ] 
                                        )
                                
                                local_xcart.append( xcart )

                            self.end.xcart = local_xcart

                    
                            if out_type == 'xcarts':
                                self.end.list_xcart.append(local_xcart) #xcart at each step only for dimer

                                #the change of typat is accounted below

                    if "number of electron " in line:
                        # print line
                        # print line.split()[-1]
                        try:
                            self.magn1 = float(line.split()[-1])
                        except:
                            self.magn1 = 0

                    if "augmentation part " in line:
                        try:
                            self.magn2 = float(line.split()[-1])
                        except:
                            self.magn2 = 0


                    if force_keyword in line:
                        # Calculate forces here...
                        forces = []
                        magnitudes = []

                        # print(self.end.select)
                        for j in range(self.end.natom):
                            parts = outcarlines[i_line+j+2].split()
                            # print "parts", parts
                            # sys.exit()
                            if hasattr(self.end, 'select') and self.end.select:
                                # print(float(parts[ff[0]]), self.end.select[j][0])
                                b = []
                                # print (self.end.select)
                                for kkk in 0,1,2:
                                    cur = self.end.select[j][kkk]
                                    # print(cur)
                                    
                                    if cur == False:# or 'F' in cur:
                                        b.append(0)
                                    elif cur == True:# or 'T' in cur:
                                        b.append(1)
                                    else:
                                        b.append(cur)
                                # print(b)
                                x = float(parts[ff[0]]) * b[0]
                                y = float(parts[ff[1]]) * b[1]
                                z = float(parts[ff[2]]) * b[2]
                            else:
                                x = float(parts[ff[0]])
                                y = float(parts[ff[1]])
                                z = float(parts[ff[2]])
                            
                            
                            forces.append([x,y,z])
                            magnitudes.append(math.sqrt(x*x + y*y + z*z))
                        # print('new step:')
                        # for f, s in zip(forces, self.end.select):
                        #     print('{:5.2f} {:5.2f} {:5.2f} {}'.format(*f, s))
                        # sys.exit()
                        average.append( red_prec( sum(magnitudes)/self.end.natom * 1000 ) )
                        imax = np.asarray(magnitudes).argmax()
                        maxforce.append( [imax, round(magnitudes[imax] * 1000)]  )
                        # mforce.append( round(magnitudes[imax] * 1000))
                        
                   

                    #Check total drift
                    if "total drift:" in line:
                        #print line
                        tdrift = [float(d) for d in line.split()[2:5]]
                        #if any(d > 0.001 and d > max(magnitudes) for d in tdrift):
                            #print_and_log("Total drift is too high = "+str(tdrift)+", check output!\n")
                            #pass


                    if "g(Stress)" in line:
                        #print line
                        gstress.append( round( float(line.split()[4])*1000 *100, 3 )  )
                    #if "Total" in line:
                        #gstress.append( red_prec(float(line.split()[4])*1000 *100, 1000 )  )
                    if "volume of cell" in line:
                        try:                     
                            self.end.vol = float(line.split()[4])
                        except ValueError: 
                            print_and_log("Warning! Cant read volume in calc "+self.name+"\n")
                        #print self.vol      

                    if "generate k-points for:" in line: 
                        self.ngkpt = tuple(  [int(n) for n in line.split()[3:]]  )
                        #print self.set.ngkpt

                      # Kohn-Sham hamiltonian: http://en.wikipedia.org/wiki/Kohn%E2%80%93Sham_equations
                      #kinetic energy
                      #+ the external potential + the exchange-correlation energy +
                      #+ Hartree (or Coulomb) energy
                    # print line
                    
                    if  "alpha Z        PSCENC" in line:
                        # print line
                        self.energy.alpha = float(line.split()[-1]) # the electrostatic interaction of the ions in a compensating electron gas.

                    if  "Ewald energy   TEWEN" in line:
                        self.energy.ewald = float(line.split()[-1]) # the electrostatic interaction of the ions in a compensating electron gas.
                        # print self.energy.ewald
                    if  "-1/2 Hartree   DENC" in line or "-Hartree energ DENC" in line:
                        self.energy.hartree = float(line.split()[-1]) #Coulomb electron-electron energy
                        # print self.energy.hartree
                    if  "-V(xc)+E(xc)   XCENC" in line:
                        self.energy.xc = float(line.split()[-1]) # Kohn-Sham exchange-correlation energy
                    if  "PAW double counting" in line:
                        self.energy.pawdc1 = float(line.split()[-2]) #
                        self.energy.pawdc2 = float(line.split()[-1]) #
                    if  "eigenvalues    EBANDS" in line:
                        self.energy.bands = float(line.split()[-1]) # - Kohn Sham eigenvalues - include kinetic energy , but not exactly
                    if  "atomic energy  EATOM" in line:
                        self.energy.atomic = float(line.split()[-1]) #energy of atoms in the box



                    if "energy  without entropy=" in line:
                        #self.energy = float(line.split()[4])
                        self.e_without_entr = float(line.split()[3]) #
                        self.energy_sigma0 = float(line.split()[6]) #energy(sigma->0)
                        self.e0 = self.energy_sigma0
                        self.list_e_sigma0.append(  self.energy_sigma0  )
                        self.list_e_without_entr.append(  self.e_without_entr  )

                        de_each_md_list.append(de_each_md)


                    if "energy without entropy =" in line:
                        e_sig0_prev = e_sig0
                        try:
                            e_sig0 = float(line.split()[7])
                        except:
                            e_sig0 = 0
                        de_each_md = e_sig0_prev - e_sig0

                    if "free  energy   TOTEN  =" in line:
                        #self.energy = float(line.split()[4])
                        self.energy_free = float(line.split()[4]) #F free energy
                    




                    if re_lengths.search(line):
                        self.vlength = [red_prec( float(l),1000 ) for l in outcarlines[i_line + 1].split()[0:3]]
                        #print self.vlength
                    if "in kB" in line:
                        # print(line)
                        line = line.replace('-', ' -')
                        # print(line)
                        if '*' in line:
                            self.stress = [0,0,0] # problem with stresses
                            printlog('Warning! Some problem with *in kB* line of OUTCAR')
                        else:
                            self.stress = [float(i)*100 for i in line.split()[2:]]  # in MPa 
                    if "Total  " in line:
                        # print(line)
                        line = line.replace('-', ' -')
                        try:
                            self.intstress = [int(float(i)*1000) for i in line.split()[1:]] #stress in internal units; can be regarded as forces
                        except:
                            self.intstress = [0,0,0]
                            printlog('Warning! Some problem with *Total * line of OUTCAR')

                    if "external pressure =" in line: 
                        #print iterat
                        self.extpress = float(line.split()[3]) * 100 # in MPa 
                        if self.mdstep == 1 : self.extpress_init = self.extpress

                    if "E-fermi :" in line: 
                        # print line
                        self.efermi = float(line.split()[2]) # in eV


                    if "Elapsed time" in line:
                        self.time = float(line.split()[3])
                    if re_nkpts.search(line):
                        self.NKPTS = int(line.split()[3])
                    if "WARNING" in line:
                        warnings += 1#line


                    if "Subroutine DYNSYM returns" in line and not nsgroup:
                        nsgroup = line.split()[4]#number of space group operations
                    # if nsgroup == None:
                    if "Subroutine GETGRP returns:" in line and not nsgroup:
                        nsgroup = line.split()[4]    


                    if "Iteration" in line:
                        self.mdstep = int(line.split('(')[0].split()[2].strip())
                        iterat +=1
                        # print self.mdstep
                        # print line
                        if mdstep_old != self.mdstep:
                            nscflist.append( niter ) # add to list number of scf iterations during mdstep_old
                        niter = int(line.split(')')[0].split('(')[-1].strip()) #number of scf iterations
                        mdstep_old = self.mdstep


                    if 'number of electron ' in line:
                        # print (line)
                        try:
                            self.mag_sum.append( [float(line.split()[5]), 0])
                        except:
                            pass

                    if 'augmentation part' in line:
                        # print (line)
                        try:
                            self.mag_sum[-1][1]= float(line.split()[4])
                        except:
                            pass

                    if 'total charge ' in line:
                        chg = []
                        for j in range(self.end.natom):
                            chg.append( float(outcarlines[i_line+j+4].split()[4]) )
                        
                        tot_chg_by_atoms.append(np.array(chg))#[ifmaglist])                    


                    if 'magnetization (x)' in line:
                        # print(line)
                        mags = []
                        for j in range(self.end.natom):
                            mags.append( float(outcarlines[i_line+j+4].split()[4]) )
                        
                        tot_mag_by_atoms.append(np.array(mags))#[ifmaglist])
                        # print(ifmaglist)
                        tot_mag_by_mag_atoms.append(np.array(mags)[ifmaglist])
                        # print tot_mag_by_atoms
                        # magnetic_elements
                        # ifmaglist
                        # self.tot_mag_by_atoms = tot_mag_by_atoms



                    if 'LDAUU' in line:
                        ldauu = line


                    if 'onsite density matrix' in line:
                        i_at = int( outcarlines[i_line-2].split()[2]  ) #starting from one
                        l_at = int( outcarlines[i_line-2].split()[8]  )
                        # print (spin_polarized)
                        spin1 = []
                        spin2 = []
                        nm = 2*l_at+1
                        for i in range(nm):
                            line = outcarlines[i_line+4+i]
                            try:
                                spin1.append( np.array(line.split()).astype(float) )
                            except:
                                print_and_log('Warning! Somthing wrong with occ matrix:', line)
                        if spin_polarized:
                            for i in range(nm):
                                # try:
                                line = outcarlines[i_line+7+nm+i]
                                # print(line)
                                line = line.replace('-', ' -')
                                spin2.append( np.array(line.split()).astype(float) )
                                # except:
                                #     printlog('Attention! Could not read spin2, probably no spaces')
                                #     spin2.append(0)        

                        occ_matrices[i_at-1] = spin1+spin2
                        # print (np.array(spin1) )


                    if 'freq' in show:

                        if 'Eigenvectors and eigenvalues of the dynamical matrix' in line:
                            freq = []

                            i = 0
                            while 'ELASTIC MODULI CONTR FROM IONIC RELAXATION' not in line:
                                i+=1
                                line = outcarlines[i_line+i]
                                if 'f  =' in line:
                                    freq.append(float(line.split()[3]) ) #THz
                                    # print(line)


                    if 'TOTAL ELASTIC MODULI' in line:
                        eltensor = []
                        for i in range(9):
                            line = outcarlines[i_line+i]
                            print(line.strip())
                            if i > 2:
                                eltensor.append([float(c)/10 for c in line.split()[1:]])

                        eltensor = np.asarray(eltensor)
                        # print(eltensor)
                        w, v = np.linalg.eig(eltensor)
                        printlog('Eigenvalues are:', w, imp = 'y')
                                # eltensor

                    if 'average eigenvalue GAMMA=' in line:
                        # print(line)
                        gamma = float(line.split()[-1])
                        if gamma > 1 and 'conv' in show:
                            printlog('average eigenvalue GAMMA >1', gamma, imp = 'y')
                        # sys.exit()



                    # if 'DIPCOR: dipole corrections for dipol' in line:
                    if self.mdstep > mdstep_prev:
                        # print(self.mdstep, dipol)
                        mdstep_prev = self.mdstep

                    if 'dipolmoment' in line:
                        dipol = line.split()[1:4]
                        self.dipol = [float(d) for d in dipol]
                        # print(line)

                        # for i in range(1,4):
                        #     line = outcarlines[i_line+i]
                        #     print(line)



                    # if 'irreducible k-points:': in line:
                    #     self.nkpt = int(line.split()[1])




                    i_line += 1
                # sys.exit()
                #Check total drift
                






            try:
                toldfe = self.set.toldfe  # eV
            except:
                toldfe = 0




            max_magnitude = max(magnitudes)
            max_tdrift    = max(tdrift)
            self.maxforce = maxforce[-1][1]
            # if max_magnitude < self.set.toldff/10: max_magnitude = self.set.toldff
            # print 'magn', magnitudes
            # print 'totdr', tdrift
            # print 'max_magnitude', max_magnitude
            try: 
                
                if max_magnitude < self.set.tolmxf: 
                    max_magnitude = self.set.tolmxf
            except:
                ''

            #if any(d > 0.001 and d > max_magnitude for d in tdrift):
            if max_tdrift > 0.001 and max_tdrift > max_magnitude:
                
                printlog( "Total drift is too high! At the end one component is {:2.1f} of the maximum force, check output!\n".format(max_tdrift)  )
                pass
            #else: maxdrift = 
            # print magn
            if tot_mag_by_atoms:
                self.end.magmom = tot_mag_by_atoms[-1].tolist()

            """update xred"""
            self.end.update_xred()









            #print "init pressure = ",self.extpress_init,"; final pressure =",self.extpress
            #print self.end.xred
            #self.vol = np.dot( self.rprimd[0], np.cross(self.rprimd[1], self.rprimd[2])  ); #volume
            nscflist.append( niter ) # add to list number of scf iterations during mdstep_old
            #print "Stress:", self.stress
            v = self.vlength
            self.end.vlength = self.vlength

            s = self.stress
            yznormal = np.cross(self.init.rprimd[1], self.init.rprimd[2])
            #print yznormal
            #print np.cross( yznormal, np.array([1,0,0]) )
            if not hasattr(self.init, 'gbpos'):
                self.init.gbpos = None#for compatability

            self.gbpos = self.init.gbpos #for compatability
            if self.gbpos:
                if any( np.cross( yznormal, np.array([1,0,0]) ) ) != 0: 
                    print_and_log("Warning! The normal to yz is not parallel to x. Take care of gb area\n")
            self.end.yzarea = np.linalg.norm( yznormal )  #It is assumed, that boundary is perpendicular to x


            """Calculate voronoi volume"""
            # print hasattr(self, 'vorovol')
            voro = ''
            if voronoi:# and not hasattr(self, 'vorovol'):#out_type == 'e_seg':
                voro = calculate_voronoi(self)
                calculate_voronoi(self, state = 'init')


            #deal with ldauu
            u_hubbard = 0
            if ldauu: 
                ldauu = np.array(ldauu.split()[7:]).astype(float)
                # print (ldauu)
                #find first non-zero
                self.ldauu = ldauu
                u_hubbard = ( next((u for u in ldauu if u), 0) )
                # print ( np.unique(ldauu)  )
            else:
                self.ldauu = [0]

            #Check if energy is converged relative to relaxation
            e_diff_md = self.energy_sigma0
            if len(self.list_e_sigma0) > 2:
                e_diff_md = (self.list_e_sigma0[-1] - self.list_e_sigma0[-2])*1000 #meV

            e_diff = (e_sig0_prev - e_sig0)*1000 #meV

            if abs(e_diff) > toldfe*1000:
                toldfe_warning = '!'
                print_and_log("Attention!, SCF was not converged to desirable prec", 
                    round(e_diff,3), '>', toldfe*1000, 'meV', imp = 'y')
            else:
                toldfe_warning = ''

            if 'conv' in show:
                for i, de in enumerate(de_each_md_list ):
                    if de/toldfe > 1.01:
                        printlog('Attention! bad SCF convergence {:6.1g} eV for MD step {:}; toldfe = {:6.0g} eV'.format(de, i+1, toldfe))


            #  Construct beatifull table
            #self.a1 = float(v[0])/2 ; self.a2 = float(v[1])/2/math.sqrt(0.75); self.c = float(v[2])  # o1b
            
            try:
                self.a = self.hex_a ; self.c = self.hex_c  # c1b
            except AttributeError:
                self.a  = 0; self.c = 0 #calculations with full relaxation
            if self.a == None or self.a == [None]:
                self.a  = 0; self.c = 0

            j = (35,12,7,7,8,9,14,5,5,20,5,20,8,12,20,8,5,8,8,25,8,4,3)

            d = "|"
            name = ("%s.%s.%s" % (self.id[0],self.id[1], str(self.id[2]) )).ljust(j[0])
            etot = ("%.4f" % ( self.energy_sigma0 )).center(j[1])
            etot1 = ("%.4f" % ( self.energy_sigma0/self.end.natom )).center(j[1])
            # print self.a
            a = ("%.4f" %      ( self.a )      ).center(j[2])
            c = ("%.4f" %      ( self.c )      ).center(j[3])
            time = ("%.3f" % (self.time/3600.)    ).center(j[4])
            itertm = ("%.1f" % (self.time/1./iterat)    ).center(j[5])
            Nmd = ("%1i,%2i,%3i" % (self.mdstep, iterat/self.mdstep, iterat)    ).center(j[6])
            self.iterat = iterat
            War = ("%i" % (warnings)    ).center(j[7])
            #nbands = ("%i" % (self.set.vasp_params["NBANDS"])    ).center(j[8])
            #added = ("%.0f" % ( (self.set.add_nbands - 1) * 100 )    ).center(j[15])
            try:
                kmesh = ("%s" % (str(self.ngkpt) )    ).center(j[8])
                ks = self.calc_kspacings()
                kspacing = ("[%.2f,%.2f,%.2f]" % ( ks[0], ks[1], ks[2] )    ).center(j[9])
                ks1 = ("[%.2f]" % ( ks[0] )    ).center(j[9])
            except:
                kmesh = ''
                ks    = ''
                kspacing = ''
                ks1     = ''

            nkpt = ("%i" % ( self.NKPTS)     ).center(j[10])
            if self.stress:
                istrs = ("[%5i,%5i,%5i] " % ( self.intstress[0],self.intstress[1],self.intstress[2]  )    ).center(j[11])
                strs = ("%.0f,%.0f,%.0f " % ( self.stress[0],self.stress[1],self.stress[2]  )    ).center(j[11])   
                eprs = ("%.0f" % (self.extpress)).center(j[12])
            
            else:
                istrs = ''
                strs  = ''
                eprs =  ''
            try:
                tsm = ("%.0f" % (self.set.tsmear*1000)).center(j[13])
            except:
                tsm = ''

            entrr = ("%.3f" % (   (self.energy_free - self.energy_sigma0)/self.end.natom * 1000    )   ).center(j[14]) #entropy due to the use of smearing

            try:
                npar = ("%i" % (self.set.vasp_params["NPAR"])).center(j[16])
                lpl = ("%s" % (self.set.vasp_params["LPLANE"])).center(j[17])
                ecut = ("%s" % (self.set.vasp_params["ENCUT"]) ).center(j[18]) 
            except:
                npar = ''
                lpl  = ''
                ecut = ''

            # lens = ("%.2f;%.2f;%.2f" % (v[0],v[1],v[2] ) ).center(j[19])
            lens = "{:4.2f};{:4.2f};{:4.2f}".format(v[0],v[1],v[2] ) 
            r1 = ("%.2f" % ( v[0] ) ).center(j[19])            
            vol = ("%.1f" % ( self.end.vol ) ).center(j[20])
            nat = ("%i" % ( self.end.natom ) ).center(j[21])
            try:
                totd = ("%.0f" % (   max_tdrift/max_magnitude * 100      ) ).center(j[22])
            except:
                totd = ''

            nsg = ("%s" % (     nsgroup     ) ).center(j[22])
            Uhu   = " {:3.1f} ".format(u_hubbard)
            ed    = ' {:3.0f}'.format( e_diff)
            edg   = ' {:3.1f} '.format( e_diff_md)
            spg   = ' {:4s} '.format( self.end.sg(silent = 1)[0])
            """Warning! concentrations are calculated correctly only for cells with one impurity atom"""
            #gbcon = ("%.3f" % (     1./self.end.yzarea      ) ).center(j[23]) # surface concentation at GB A-2
            #bcon = ("%.1f" % (     1./self.natom * 100      ) ).center(j[24]) # volume atomic concentration, %
            #outstring_nbands = name+d+Etot+d+a+d+c+d+time+d+itertm+d+Nmd+d+War+d+nbands+d+added+"\\\\"
            #outstring_npar = name+d+Etot+d+a+d+c+d+time+d+itertm+d+Nmd+d+War+d+npar+d+lpl+"\\\\"

            #outstring_stress = name+d+Etot+d+a+d+c+d+time+d+itertm+d+Nmd+d+War+d+istrs+d+eprs

            #outstring_kp_ec = name+d+Etot+d+a+d+c+d+time+d+itertm+d+Nmd+d+War+d+strs+d+eprs+d+kmesh+d+ecut+"\\\\"
            outst_ecut= etot+d+a+d+c+                                            d+time+d+itertm+d+Nmd+d+War+d+ecut+"\\\\"
            outst_kp  = etot+d+a+d+c+                                            d+time+d+itertm+d+Nmd+d+War+d+kmesh+d+kspacing+d+nkpt+"\\\\"            
            outst_ts  = etot+d+a+d+c+                                            d+time+d+itertm+d+Nmd+d+War+d+kmesh+d+tsm+d+entrr+"\\\\"

            outst_all = voro+etot+d+a+d+c+d+lens+d+vol+d+kspacing+d+strs+d+eprs+d+nat+d+time+d+Nmd+d+War+d+totd+d+nsg+"\\\\"
            outst_seg = voro+etot+d+        lens+d+vol+d+ks1     +d+strs+d+eprs+d+nat+d+time+d+Nmd+d+War+d+totd+d+nsg+"\\\\" #for segregation
            outst_coseg=voro+etot+d+                                strs+d+eprs+d+nat+d+time+d+Nmd+d+War+d+totd+d+nsg+"\\\\" #for co-segregation; 
            outst_gbe = voro+etot+               d+vol+d+kspacing+d+strs+d+eprs+d+nat+d+time+d+Nmd+d+War+d+nsg+"\\\\" # For comparing gb energies and volume
            outst_imp = voro+etot+d+a+d+c+d+lens+d+vol+d+kspacing+d+       eprs+d+nat+d+time+d+Nmd+d+War+d+totd+d+nsg+"\\\\" # For comparing impurity energies
            
            outst_cathode = d.join([spg,etot, etot1, lens, vol,nkpt, strs, nat, time, Nmd, War, nsg, Uhu, ed, edg ])
            # print self.end.xred[-1]
            #print outstring_kp_ec
            # print show
            # print 'force' in show


            if 'conv' in show:
                # print('asdf', de_each_md_list)
                # show achived convergence every step with respect to toldfe, should be less than 1
                # np.set_printoptions(linewidth=150, formatter={'float':lambda x: "%3.0f->" % x}) #precision=1,
                np.set_printoptions(precision=0, linewidth=150, )
                printlog('Conv each step, de/toldfe (toldfe = {:.0g} eV) =  \n{:};'.format(toldfe, np.array([de/toldfe for de in de_each_md_list ])), imp = 'Y')
            



            if 'fo' in show:
                # print "Maxforce by md steps (meV/A) = %s;"%(str(maxforce)  )
                print_and_log("\n\nMax. F."+force_prefix+" (meV/A) = \n{:};".format(np.array([m[1] for m in maxforce ])[:]  ), imp = 'Y'  )
                # print "\nAve. F. (meV/A) = \n%s;"%(  np.array(average)  )
                # import inspect
                # print inspect.getargspec(plt.plot).args
                # print plt.plot.__doc__
                if 'p' in show[0]:
                    plt.plot(maxforce, )
                    plt.xlabel('MD step')
                    plt.ylabel('Max. force on atom (meV/$\AA$)')
                    plt.show()
            
            if 'sur' in show:
                self.sumAO = {}
                self.devAO = {}
                for el in 'Li', 'Na', 'Fe', 'O':
                    if el in self.end.get_elements():
                        pos  = determine_symmetry_positions(self.end, el)
                        # print(pos)
                        # sys.exit()
                        # xc = self.end.xcart[pos[0][0]]
                        for ps in pos:
                            print('position', ps[0])
                            xc = self.end.xcart[ps[0]]

                            if el == 'O':
                                neib = 6
                            else:
                                neib = 6
                            sumAO = local_surrounding(xc, self.end, neib, periodic = True, only_elements = [8, 9], control = 'av')#[0]
                            self.devAO[el+'-O'] = local_surrounding(xc, self.end, neib, periodic = True, only_elements = [8, 9], control = 'av_dev')[0]
                            print('d_av '+el+'-O:',sumAO )
                            print('dev_av '+el+'-O:',self.devAO[el+'-O'] )
                            AO = local_surrounding(xc, self.end, neib, periodic = True, only_elements = [8, 9], control = 'mavm')
                            print('d_min, d_avex, d_max: {:4.2f}, {:4.2f}, {:4.2f}'.format(*AO))




                            self.sumAO[el+'-O'] = sumAO
                            if self.id[2] in [1,12]:
                                self.end.write_xyz(show_around_x = xc, nnumber = neib, filename = self.end.name+'_'+el+'-OF'+str(neib), analysis = 'imp_surrounding', only_elements = [8, 9])

                            if el in ['Li', 'Na']:
                                neib = 2
                                sumAO = local_surrounding(xc, self.end, neib, periodic = True, only_elements = [26,], control = 'av')#[0]
                                if self.id[2] in [1,12]:
                                    self.end.write_xyz(show_around_x = xc, nnumber = neib, filename = self.end.name+'_'+el+'-Fe'+str(neib), analysis = 'imp_surrounding', only_elements = [26])
                                
                                print(el+'-Fe',sumAO )
                                self.sumAO[el+'-Fe'] = sumAO


            if 'en' in show:
                    maxf = [m[1] for m in maxforce ]
                    # print(maxf)
                    plt.plot(maxf, 1000*(np.array(self.list_e_sigma0)-self.energy_sigma0) , '-o')
                    # plt.xlabel('MD step')
                    # plt.ylabel('Energy per cell (eV')
                    plt.xlabel('Max. force on atom (meV/$\AA$)')
                    plt.ylabel('Energy per cell relative to min (meV)')

                    plt.show()

            if 'smag' in show:
                # printlog('{:s}'.format([round(m) for m in self.mag_sum]), imp = 'Y' )
                printlog(np.array(self.mag_sum).round(2), imp = 'Y' )

            if 'mag' in show or 'occ' in show:
                from siman.analysis import around_alkali
                numb, dist, chosen_ion = around_alkali(self.end, 4, alkali_ion_number)
                
                #probably not used anymore
                # dist_dic = {}
                # self.dist_numb = zip(dist, numb)
                # for d, n in self.dist_numb:
                #     dist_dic[n] = d 
                #probably not used anymore


            if 'mag' in show and tot_mag_by_atoms:
                print ('\n\n\n')
                # print_and_log
                # print 'Final mag moments for atoms:'
                # print np.arange(self.end.natom)[ifmaglist]+1
                # print np.array(tot_mag_by_atoms)

                # print (tot_mag_by_atoms)
                # if tot_mag_by_atoms:
                # print ('first step ', tot_mag_by_atoms[0][numb].round(3) )
                # print ('first step all ', tot_mag_by_atoms[0][ifmaglist].round(3) )
                # for mag in tot_mag_by_atoms:
                #     print ('  -', mag[numb].round(3) )

                # print ('last  step ', tot_mag_by_atoms[-1][numb].round(3), tot_chg_by_atoms[-1][numb].round(3) )
                mmm = tot_mag_by_atoms[-1][numb].round(3)

                print ('atom:mag  = ', ', '.join('{}:{:4.2f}'.format(iat, m) for iat, m  in zip(  numb+1, mmm   )) )
                if 'a' in show:
                    ''
                    # print ('last  step all', tot_mag_by_atoms[-1][ifmaglist].round(3) )

                    # sys.exit()
                if chosen_ion:
                    printlog ('Dist from 1st found alkali ion ',element_name_inv( chosen_ion[1]),
                        ' to sur. transition met atoms: (Use *alkali_ion_number* to choose ion manually)')
                    print ('atom:dist = ', 
                    ', '.join('{}:{:.2f}'.format(iat, d) for iat, d  in zip(  numb+1, dist   )  ) )

                # plt.plot(np.array(sur[3]).round(2), tot_mag_by_atoms[-1][numb]) mag vs dist for last step
                
                # print ('Moments on all mag atoms:\n', tot_mag_by_atoms[-1][ifmaglist].round(3))
                if 'p' in show:
                    plt.plot(np.array(tot_mag_by_mag_atoms)) # magnetization vs md step
                    plt.show()
                    plt.clf()

            if 'chg' in show:
                self.tot_chg_by_atoms = tot_chg_by_atoms[-1] #save for last step only
                # print(list(zip(self.end.get_elements(), self.tot_chg_by_atoms)))
                els  = self.end.get_elements()
                try:
                    only_el = show.split('.')[1:]
                except:
                    only_el = None
                print('\nMulliken charges are:')
                for el, ch in zip(els, self.tot_chg_by_atoms):
                    if only_el == None or (only_el and el in only_el):
                        print('{:s} {:4.2f};'.format(el, ch), end = ' ')
                print()
            if 'occ' in show:
                ''
                # print (matrices)
                # print (df)
                if chosen_ion:
                    print_and_log('Distances (A) from alkali ion #',chosen_ion[0]+1,' to transition atoms:', 
                        ',  '.join([ '({:}<->{:}): {:.2f}'.format(chosen_ion[0]+1, iat, d) for d, iat in zip(  dist, numb+1  )  ]), imp = 'Y'  )
                
                show_occ_for_atoms = [int(n) for n in re.findall(r'\d+', show)]
                # print (show_occ_for_atom)
                # sys.exit()
                if show_occ_for_atoms:
                    iat = show_occ_for_atoms[0]-1
                    # dist_toi = dist_dic[iat]
                    i_mag_at = iat
                else:
                    i = 0
                    # dist_toi = dist[i]
                    i_mag_at = numb[i]
                # print (st.znucl[st.typat[i_mag_at]-1] )
                l05 = len(occ_matrices[i_mag_at])//2

                df = pd.DataFrame(occ_matrices[i_mag_at]).round(5)

                print_and_log( 'Occ. matrix for atom ', i_mag_at+1, end = '\n', imp = 'Y'  )
                    # ':  ; dist to alk ion is ',  dist_toi, 'A', end = '\n' )
                print_and_log('Spin 1:',end = '\n', imp = 'Y'  )
                print_and_log(tabulate(df[0:l05], headers = ['dxy', 'dyz', 'dz2', 'dxz', 'dx2-y2'], floatfmt=".1f", tablefmt='psql'),end = '\n', imp = 'Y'  )
                # print(' & '.join(['d_{xy}', 'd_{yz}', 'd_{z^2}', 'd_{xz}', 'd_{x^2-y^2}']))
                # print_and_log(tabulate(occ_matrices[i_mag_at][l05:], headers = ['d_{xy}', 'd_{yz}', 'd_{z^2}', 'd_{xz}', 'd_{x^2-y^2}'], floatfmt=".1f", tablefmt='latex'),end = '\n' )
                # print(tabulate(a, tablefmt="latex", floatfmt=".2f"))
                print_and_log('Spin 2:',end = '\n', imp = 'Y'  )
                print_and_log(tabulate(df[l05:], floatfmt=".1f", tablefmt='psql'), imp = 'Y'  )
            self.occ_matrices = occ_matrices
            


            if 'freq' in show:
                dos = [1]*len(freq)
                # from scipy.ndimage.filters import gaussian_filter
                from scipy.signal import butter, lfilter, freqz
                # blurred = gaussian_filter(freq, sigma=7)
                fmin = min(freq)
                fmax = max(freq)
                fw   = fmax-fmin

                finefreq = np.linspace(fmin, fmax, 1000)
                dos = [0]*1000

                # for i in range(1000):
                # print(fw)
                for f in freq:
                    # print(f)
                    i = int( np.round( (f-fmin)/ fw * 999 ,0) )
                    dos[i] = 1
                    # print(i, finefreq[i], f)
                

                def butter_lowpass(cutoff, fs, order=5):
                    nyq = 0.5 * fs
                    normal_cutoff = cutoff / nyq
                    b, a = butter(order, normal_cutoff, btype='low', analog=False)
                    return b, a

                def butter_lowpass_filter(data, cutoff, fs, order=5):
                    b, a = butter_lowpass(cutoff, fs, order=order)
                    y = lfilter(b, a, data)
                    return y

                order = 6
                fs = 30.0       # sample rate, Hz
                cutoff = 3.667  # desired cutoff frequency of the filter, Hz

                y = butter_lowpass_filter(finefreq, cutoff, fs, order)

                plt.plot(finefreq, smoother(smoother(dos,50), 50), '-') 
                plt.savefig('figs/'+str(self.id)+'.eps')
                # plt.show()
                plt.clf()

            # sys.exit()


            printlog("Reading of results completed\n\n", imp = 'n')
            self.end.outfile = path_to_outcar
            

            if pymatgen_flag:
                ''
                # self.end.write_cif(os.path.join(self.dir,self.name))
            

            # print(out_type)
            # sys.exit()
            if   out_type == 'gbe'  : outst = outst_gbe
            elif out_type == 'e_imp': outst = outst_imp
            elif out_type == 'e_seg': outst = outst_seg            
            elif out_type == 'coseg': outst = outst_coseg            
            elif 'ecut' in out_type : outst = outst_ecut
            elif 'kp' in out_type   : outst = outst_kp
            elif 'ts' in out_type   : outst = outst_ts
            
            elif not header.siman_run:
                outst_simple = '|'.join([etot, lens, strs, Nmd])
                # print("Bi2Se3.static.1               |  -20.1543  |    10.27;10.27;10.27    | -680,-680,-657 |   1,13, 13   |    ")
                if header.show_head:
                    printlog("name                          |  energy(eV)|    Vector lenghts (A)   | Stresses (MPa)     | N MD, N SCF   ", end = '\n', imp = 'Y')
                    header.show_head = False
                
                outst = outst_simple
            else: 
                printlog('Output type: outst_cathode')
                outst = outst_cathode
            #else: print_and_log("Uknown type of outstring\n")


            #save cif file


        else:
            # if not hasattr(cl,'energy_sigma0'):
            cl = self
            try:
                os.rename(cl.path['output'], cl.path['output']+"_unfinished") 
                printlog('read_results():',cl.id, 'is unfinished, continue:', cl.dir, cl.cluster_address, imp = 'y')
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



    def get_file(self, filetype = '', nametype = '', up = 'up1'):
        #allow to get any file of type filetype 
        #cl - object of CalculationVasp class
        #filetype (str) - 'CHG', 'CHGCAR', etc just the name of file in calculation folder
        #nametype (str) - 'asoutcar' - update filetype to OUTCAR format
        #up = up1 - donot update
              #up2 - update

        #Comment
            #initially used for chg files - rename!
        if nametype == 'asoutcar':
            path_to_file = self.path['output'].replace('OUTCAR',filetype)
        else:
            path_to_file = os.path.dirname(self.path['output']) +'/'+ filetype

        if 'CHGCAR' in filetype:
            self.path['chgcar'] = path_to_file
        elif 'xml' in filetype:
            self.path['xml'] = path_to_file


        # print(path_to_file)
        # print(self.cluster_address)
        # print(self.project_path_cluster+'/')
        # sys.exit()
        path2file_cluster = self.project_path_cluster+'/'+path_to_file
        if os.path.exists(path_to_file) and 'up2' not in up: 
            out = None
        else:
            out = get_from_server(path2file_cluster, os.path.dirname(path_to_file), addr = self.cluster_address)


        if out:
            printlog('File', path2file_cluster, 'was not found, trying archive:',header.PATH2ARCHIVE, imp = 'Y')
            # printlog('Charge file', path_to_file, 'was not found')
            try:
                pp = self.project_path_cluster.replace(self.cluster_home, '') #project path without home
            except:
                pp = ''
            # print(pp)
            path_to_file_scratch = header.PATH2ARCHIVE+'/'+pp+'/'+path_to_file

            out = get_from_server(path_to_file_scratch, os.path.dirname(path_to_file), addr = self.cluster_address)
            
            if out:
                printlog('File', path_to_file_scratch, 'was not found', imp = 'Y')
                path_to_file = None
           
        printlog('File', path_to_file, ' was download', imp = 'e')
        
        return path_to_file

    def get_chg_file(self, *args, **kwargs):
        """just wrapper to get chgcar files """
        del kwargs['CHGCAR']
        return self.get_file(self, filetype = 'CHGCAR', **kwargs)




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
            return run_on_server(cmd, self.cluster_address)



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


    def check_job_state(self):
        #check if job in queue or Running

        cl = self
        if header.check_job == 1:
            job_in_queue = ''
            if hasattr(cl,'schedule_system'):

                check_string =  cl.id[0]+'.'+cl.id[1]
                if 'SLURM' in cl.schedule_system:

                    job_in_queue = check_string in run_on_server("squeue -o '%o' ", cl.cluster_address)
                    printlog(cl.id[0]+'.'+cl.id[1], 'is in queue or running?', job_in_queue)

                elif 'PBS' in cl.schedule_system:
                    job_in_queue = check_string in run_on_server("qstat -x ", cl.cluster_address)

                else:
                    print_and_log('Attention! unknown SCHEDULE_SYSTEM='+'; Please teach me here! ', imp = 'y')
                    job_in_queue = ''


            if file_exists_on_server(os.path.join(cl.dir, 'RUNNING'), addr = cl.cluster_address) and job_in_queue: 
                
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




    def res(self, **argv):
        from siman.calc_manage import res_loop
        res_loop(*self.id, **argv)

    def run(self, ise, iopt = 'full_nomag', up = 'up1', vers = None, i_child = -1, add = 0, *args, **kwargs):
        """
        Wrapper for add_loop (in development)
        By default inherit self.end
        ise - new ise

        iopt - inherit_option
            'full_nomag'
            'full'
            'full_chg' - including chg file
        vers - list of version for which the inheritance is done

        i_child - choose number of child to run res_loop()
        add - if 1 than add new calculation irrelevant to children
        TODO:
        1. if ise is not provided continue in the same folder under the same name,
        however, it is not always what is needed, therefore use inherit_xred = continue
        """

        from siman.calc_manage import add_loop

        if not iopt:
            iopt = 'full'



        # if self.id[1] != ise:
        if 1:
            if not hasattr(self, 'children'):
                self.children = []

            if not add and len(self.children)>0:
                print('Children were found:', self.children, 'by defauld reading last, choose with *i_child* ')
                
                for i in self.children:
                    # print(i, ise, i[1], i[1] == ise)
                    if i[1] == ise:
                        idd = i
                        # add = True
                        # break
                    else:
                        idd = None

                if idd is None:
                    add = True
                    # idd  = self.children[i_child]
                # print(idd)

                if idd:
                    cl_son = header.calc[idd]
                    
                    cl_son.res(up = up, **kwargs)
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


    def read_pdos_using_phonopy(self, mode = 'pdos', plot = 1, up = 'up1'):
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

        self.get_file('vasprun.xml', nametype = 'asoutcar', up = up)
        create_phonopy_conf_file(self.end, mp = [10, 10, 10], path = self.dir)
        # create_phonopy_conf_file(self.end, mp = [36, 36, 36], path = self.dir) #almost no difference was found for Na2X
        create_phonopy_conf_file(self.end, path = self.dir, filetype = 'band') #create band file



        cwd = os.getcwd()


        os.chdir(self.dir)
        print(self.dir)
        out = runBash('phonopy --fc '+os.path.basename(self.path['xml']))

        printlog('phonopy out: ', out)



        if 'poscar' not in self.path:
            self.path['poscar'] = self.path['output'].replace('OUTCAR','POSCAR')

        if mode == 'pdos':
            print('phonopy -c '+os.path.basename(self.path['poscar'])+p+'  mesh.conf --readfc ')
            runBash('phonopy -c '+os.path.basename(self.path['poscar'])+p+' mesh.conf --readfc ')

        from siman.calc_manage import read_phonopy_dat_file

        self.pdos = read_phonopy_dat_file('total_dos.dat')


        #phonons
        

        if mode == 'band':
            print('phonopy -c '+os.path.basename(self.path['poscar'])+' -p band.conf --readfc ')
            runBash('phonopy -c '+os.path.basename(self.path['poscar'])+' -p band.conf --readfc ')

        if mode == 'free':
            print('phonopy -c '+os.path.basename(self.path['poscar'])+' -t -p mesh.conf --readfc ')

            runBash('phonopy -c '+os.path.basename(self.path['poscar'])+' -t' +p+' mesh.conf --readfc ')


            Trange, func = read_phonopy_data('thermal_properties.yaml', convert = 1)

            self.F = func # free energy function in eV, still for the whole supercell!
            # print(self.id, self.F)


        os.chdir(cwd)

        return

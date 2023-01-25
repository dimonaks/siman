# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
import itertools, os, copy, math, glob, re, shutil, sys, pickle, gzip, shutil, random
import re, io, json
import pprint

from textwrap import wrap

import numpy as np

#additional packages
from siman import header

try:
    from tabulate import tabulate
except:
    print('tabulate is not available')
    tabulate = None
try:
    import pandas as pd
except:
    print('pandas is not available')

try:
    import pymatgen
    header.pymatgen_flag = True
except:
    print('calculation.py: pymatgen is not available')
    header.pymatgen_flag = False

if header.pymatgen_flag:
    from pymatgen.io.cif import CifWriter
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.core.surface import Slab
    from pymatgen.core.composition import Composition

from siman.header import printlog

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


class Calculation(object):
    """Main class of siman. Objects of this class contain all information about first-principles calculation
        List of important fields:
            - init (Structure) - the initial structure used for calculation
            - end  (Structure) - the structure obtained after optimization
            - occ_matrices (dict) - occupation matrices, number of atom (starting from 0) is used as key


    """
    def __init__(self, inset = None, iid = None, output = None):
        #super(CalculationAbinit, self).__init__()
        self.name = "noname"
        if inset:
            self.set = copy.deepcopy(inset)
        else:
            self.set = InputSet(calculator = self.calculator)
        
        # if self.set.set_sequence:


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
        """Reads geometrical data from filename file in abinit format
            should be moved to Structure class
        """
        from siman.classes import empty_struct
        from siman.core.structure import Structure

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
        """Writes geometrical data in custom siman format bases on abinit format 
        to self.path["input_geo"]


        TODO: move to Structure class
        """
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




    def dos(self, isym = None, el = None, i_at = None, iatoms = None, multi = None,  *args, **kwargs):
        """
        Plot DOS either for self or for children with dos

        INPUT:
            - isym (int) - choose symmetry position to plot DOS, otherwise use i_at
            - i_at - number of atom from 0
            - iatoms - list of atom numbers (from 0) to make one plot with several dos
            - el - element for isym, otherwise first TM is used orbitals
            - multi (dict) - special dict to make multiplots, 
                    default is {'first':1, 'last':1, 'ax':None, 'hide_xlabels':False, 'pad':None}; see fit_and_plot()
            - nneighbors

            fit_and_plot arguments can be used

        RETURN:

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
        dostype = pm.get('dostype') or 'partial'
        nneighbors = pm.get('nneighbors') or 6

        if efermi_origin is None:
            efermi_origin = 1


        xlabel = '$E-E_F$, eV'
        if multi is None:
            multi = {'first':1, 'last':1, 'ax':None, 'hide_xlabels':False, 'pad':None}
            ylabel = 'Total DOS, states/eV'
        else:
            ylabel = None
            if multi['hide_xlabels']:
                xlabel = None

        

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




        if not iatoms:
            #just one plot
            plot_dos(cl,  iatom = iTM+1,  
            dostype = dostype, orbitals = orbitals, 
            labels = labels, 
            nsmooth = nsmooth, 
            image_name = image_name, 
            # invert_spins = invert_spins,
            efermi_origin = efermi_origin,
            efermi_shift = efermi_shift,
            neighbors = nneighbors,
            show = 0,  plot_param = {
            'figsize': (6,3), 
            'first':multi['first'], 'last':multi['last'], 'ax':multi['ax'], 'pad':multi['pad'], 'hide_xlabels':multi['hide_xlabels'],
            'xlabel':xlabel, 'ylabel':ylabel,
            'corner_letter':corner_letter,
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

            if corner_letter is None:
                corner_letter = 1

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
                


                # print()
                plot_dos(cl,  iatom = iat+1,  efermi_origin = efermi_origin,
                dostype = 'partial', orbitals = orbitals, 
                neighbors = nneighbors,
                labels = ['', ''], 
                nsmooth = 1, 
                color_dict = color_dicts[i%2],
                image_name = image_name, 

                # invert_spins = invert_spins,
                # show_gravity = (1, 'p6', (-10, 10)), 
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














    def plot_locpot(self, filename = None):
        """
        Plot local electrostatic potential
        """
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
        Try to make it general for different codes
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
        if hasattr(self.set, 'k_band_structure'):
            band_structure = self.set.k_band_structure
        else:
            band_structure = None
        kspacing = self.set.vasp_params['KSPACING']
        # self.set.printme()
        # print(self)
        # print(kspacing)
        # sys.exit()
        # print (self.set.kpoints_file)
        # sys.exit()
        # printlog("Warning! If you use *update_set_dic* in add(), check_kpoints() works for set before update. Please fix.") #done
        if ngkpt:
            N = ngkpt

        elif is_string_like(self.set.kpoints_file):
            print_and_log("External K-points file was provided", self.set.kpoints_file)
            N = None

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
        elif band_structure:
            printlog('check_kpoints(): the following path is used for band structure ',band_structure)
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

                elif 'simple' in cl.schedule_system:
                    print_and_log('For SCHEDULE_SYSTEM='+cl.schedule_system+' please manually run on server! ', imp = 'y')                    


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

    def push_file(self, file, subfolder=''):
        """
        Upload file to server to calculation folder into *subfolder*
        """

        push_to_server([file],  self.project_path_cluster +'/'+ self.dir + '/'+subfolder, self.cluster_address)

        return  

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
                # printlog(out, imp  = 'Y')

        if out:
            if header.PATH2ARCHIVE:
                # printlog('Charge file', path_to_file, 'was not found')
                try:
                    pp = self.project_path_cluster.replace(self.cluster_home, '') #project path without home
                except:
                    pp = ''
                # print(pp)
                path_to_file_scratch = header.PATH2ARCHIVE+'/'+pp+'/'+path_to_file

                out = get_from_server(path_to_file_scratch, os.path.dirname(path_to_file), addr = self.cluster['address'])
                
                if out:
                    printlog('File', path_to_file, ' was downloaded from archive', imp = 'e')

                else:
                    printlog('File', path_to_file_scratch, 'was not found, trying local archive path', imp = 'Y')


            else:
                printlog('project_conf.PATH2ARCHIVE is empty, trying local archive path', imp = 'Y')


        if out: 
            if header.PATH2ARCHIVE_LOCAL:

                path_to_file_AL= header.PATH2ARCHIVE_LOCAL+'/'+path_to_file
                # print(path_to_file_AL)
                if os.path.exists(path_to_file_AL):
                    shutil.copyfile(path_to_file_AL, path_to_file)

                    printlog('Succefully found in local archive and copied to the current project folder', imp = 'Y')
                else:
                    printlog('File', path_to_file_AL,'not found in local archive', imp = 'Y')

            else:
                path_to_file = None
                printlog('project_path.PATH2ARCHIVE_LOCAL is empty', imp = 'Y')



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
            it_suffix_del = False
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
                    if kwargs.get('it_suffix'):
                        del kwargs['it_suffix']
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

    def show_force(self,):
        force_prefix = ' tot '

        printlog("\n\nMax. F."+force_prefix+" (meV/A) = \n{:};".format(np.array([m[1] for m in self.maxforce_list ])[:]  ), imp = 'Y'  )




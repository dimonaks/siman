# -*- coding: utf-8 -*- 

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
from header import *
import header
from functions import (read_vectors, read_list, words, local_surrounding, 
    xred2xcart, xcart2xred, element_name_inv, calculate_voronoi,
    get_from_server, push_to_server, list2string)



class Description():
    """
    Objects of this class include just folder and description of specific calculation.
    Mostly was needed for manual addition of new calculations

    """
    def __init__(self, sectionfolder = "forgot_folder", description = "forgot_description"):
        self.des = description
        self.sfolder = sectionfolder




class empty_struct():
    def __init__(self):
        pass

class Structure():
    """This class includes only structure related information such as primitive vectors, coordinates, forces and so on"""
    def __init__(self):
        self.name = ""
        self.des = ''
        self.rprimd = [np.zeros((3)) for i in 0,1,2]
        self.xcart = []
        self.xred = []

    def xcart2xred(self,):
        self.xred = xcart2xred(self.xcart, self.rprimd)

    def xred2xcart(self,):
        self.xcart = xred2xcart(self.xred, self.rprimd)


    def add_atoms(self, atoms_xcart, element = 'Pu'):
        """
        appends at the end. Takes care of ntypat, typat, znucl, nznucl, xred and natom
        Returns Structure()
        """

        print_and_log('Warning! Method add_atoms() was not carefully tested ')
        print_and_log('Adding atom ', element)

        st = copy.deepcopy(self)

        # print type(atoms_xcart)

        st.xcart.extend(atoms_xcart)
        st.xcart2xred()
        natom_to_add = len(atoms_xcart)
        # print natom_to_add

        st.natom+=natom_to_add
        el_z_to_add = element_name_inv(element)
        print 'el_z_to_add', el_z_to_add

        if el_z_to_add not in st.znucl:
            
            st.znucl.append( el_z_to_add )
            
            st.nznucl.append(natom_to_add)            

            st.ntypat+=1
            typ = max(st.typat)+1
        
        else:
            i = st.znucl.index(el_z_to_add)
            print i
            st.nznucl[i]+=natom_to_add
            typ = i+1
            print typ
            print st.znucl
            print st.nznucl

        # sys.exit()



        st.typat.extend( [typ]*natom_to_add )



        return st

    def del_atoms(self, iat):
        """
        Now can delete only one atom with number iat (int), starting from 0. Takes care of ntypat, typat, znucl, nznucl, xred and natom
        Returns Structure()
        """


        print_and_log('Warning! Method del_atoms() was not carefully tested ')
        st = copy.deepcopy(self)


        i = iat

        typ = st.typat[i]

        del st.typat[i]
        del st.xred[i]
        del st.xcart[i]

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



        return st


    def mov_atoms(self, iat = None, to_x = None):

        st = copy.deepcopy(self)
        st.xcart[iat] = to_x
        st.xcart2xred()

        return st


    def leave_only(self, atom_type = None):
        #Remove all atoms except *atom_type*(str, mendeleev element name)
        print_and_log('Starting leave_only()', imp = 'n')

        st = copy.deepcopy(self)
        
        print_and_log('    N of atoms before = ',st.natom, imp = 'n')


        z = element_name_inv(atom_type)

        new_xred = []
        for t, xr in zip(st.typat, st.xred):
            if st.znucl[t-1] == z:
                new_xred.append(xr)

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


    def find_atom_num_by_xcart(self, x_tar, prec = 1e-10):
        for i, x in enumerate(self.xcart):
            if np.linalg.norm(x-x_tar) < prec:
            # if all(x == x_tar):
                return i
        # print self.xcart.index(x_tar)
        # return self.xcart.index(x_tar)


        # print type(atoms_xcart)


class Calculation():
    """Main class of siman. Objects of this class contain all information about first-principles calculation"""
    def __init__(self, inset = None):
        #super(CalculationAbinit, self).__init__()
        self.name = "noname"
        self.set = copy.deepcopy(inset)
        self.end = Structure()
        self.state = "1.Initialized"
        self.path = {
        "input":None,
        "input_geo":None,
        "potential":None,
        "output":None}
        self.calc_method = None #



    def read_geometry(self,filename = None):
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
                    print "\nFile contain build information! Start to read"
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

                    print "Build information has been read\n"




            self.init = Structure()

            #sys.exit()
            self.useable = 0
            #Read total number of atoms
            #Programm nznucl, since We will have more impurity atoms of different types
            command="""grep -w -m 1 "natom " """+filename
            s1=runBash(command)
            # print command
            # print s1
            self.natom_str = s1
            if s1=='':
                self.natom = 0
                print_and_log( """Warning! In filename """+filename+""" not found natom! set to zero.
                It is very likly that other parameters was not 
                found too, Calculation completly not useable!!!""")
                raise RuntimeError
            else:
                self.natom=int(s1.split()[1]) 

            self.acell = read_list("acell", 3, float, gen_words)
            self.rprim = read_vectors("rprim", 3, gen_words)

            self.rprimd = copy.deepcopy( self.rprim )
            for i in 0,1,2:
                self.rprimd[i] = self.rprim[i] * self.acell[i]         #Calculate rprimd
            #Determine reciprocal vectors
            self.recip = []
            self.vol = np.dot( self.rprimd[0], np.cross(self.rprimd[1], self.rprimd[2])  ); #volume
            #print vol
            self.recip.append(   np.cross( self.rprimd[1], self.rprimd[2] )   )
            self.recip.append(   np.cross( self.rprimd[2], self.rprimd[0] )   )
            self.recip.append(   np.cross( self.rprimd[0], self.rprimd[1] )   )
            for i in 0,1,2:
                self.recip[i] =  self.recip[i] * 2 * math.pi / self.vol;
            #print self.recip

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
                print "Convert xcart to xred"
                self.xred = xcart2xred(self.xcart, self.rprimd)
            
            if self.xcart == [None]:
                print "Convert xred to xcart"
                self.xcart = xred2xcart(self.xred, self.rprimd)

            self.hex_a = read_list("hex_a", 1, float, gen_words)[0]
            self.hex_c = read_list("hex_c", 1, float, gen_words)[0]
            self.len_units = read_list("len_units", 1, str, gen_words)[0]

            self.version = read_list("version", 1, int, gen_words)[0]

            self.gbpos = read_list("gbpos", 1, float, gen_words)[0]




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
            if vel[0] != None: self.init.vel = vel

            #read magnetic states; name of vasp variable
            self.init.magmom = read_list("magmom", self.natom, float, gen_words)
            # self.init.mag_moments 







            self.state = "2.Geometry has been read"



        #file.close();

        print_and_log( "If no warnings, geometry has been succesfully read from file "+filename+" \n")

        return





    def write_geometry(self, geotype = "init", description = "", override = False):
        """Writes geometrical data in custom siman format bases on abinit format to self.path["input_geo"]"""
        geo_dic = {}
        geofile = self.path["input_geo"]
        geo_exists = os.path.exists(geofile)

        if geo_exists:
            if override:
                print_and_log("Warning! File "+geofile+" was replaced\n"); 
            else: 
                print_and_log("Error! File "+geofile+" exists. To replace it set parameter override\n"); 
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
            if not hasattr(st, 'natom'):  st.natom = self.init.natom
            if not hasattr(st, 'ntypat'): st.ntypat = self.init.ntypat
            if not hasattr(st, 'typat'): st.typat = self.init.typat
            if not hasattr(st, 'znucl'): st.znucl = self.init.znucl 
        else: print_and_log("Error! Unknown geotype \n"); raise RuntimeError                                  

        if st.natom != len(st.xred) != len(st.xcart) != len(st.typat) or len(st.znucl) != max(st.typat): 
            print "Error! write_geometry: check your arrays.\n\n" 



        with open(self.path["input_geo"],"w") as f:
            f.write("des "+description+"\n")
            f.write("len_units "+self.len_units+"\n")
            f.write("hex_a "+str(self.hex_a)+"\n")
            f.write("hex_c "+str(self.hex_c)+"\n")
            try: self.gbpos
            except AttributeError:
                self.gbpos = None
            f.write("gbpos "+str(self.gbpos)+"\n")
            f.write("version "+str(self.version)+"\n")
            
            try: st.magmom
            except AttributeError:
                st.magmom = None
            print st.magmom 
            # sys.exit()
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



class CalculationAbinit(Calculation):
    """docstring for CalculationAbinit"""
    pass









class CalculationVasp(Calculation):
    """Methods for calculations made using VASP DFT code"""



    def read_poscar(self, filename):
        """Read POSCAR file
        """
        if self.path["input_geo"] == None:
            self.path["input_geo"] = filename
        
        self.len_units = 'Angstrom'
        self.hex_a = None
        self.hex_c = None
        self.gbpos = None


        #Determine version
        self.version = None
        try:
            print_and_log('Trying to find version at the end of filename POSCAR-v ...')
            ver = int(filename.split('-')[-1])
            # print filename.split('-')[-1]
            # print 'ver', ver
            print_and_log('OK\n')

        except:
            try:
                print_and_log('\nTrying to find version at the begenning of filename v.POSCAR...')
                ver = int(os.path.basename(filename).split('.')[0] )
                print_and_log('OK\n')
            
            except:
                raise RuntimeError    

        # sys.exit()

        self.version = ver
        # except:
        #     pass
        # print self.version, 'version'



        self.init = Structure()
        st = self.init
        with open(filename,'r') as f:
            name = f.readline().strip()
            # print self.name, "self.name"

            self.des = name
            # st.name = self.name

            mul = float( f.readline() )
            # print 'mul', mul


            st.rprimd = []
            for i in 0, 1, 2:
                vec = f.readline().split()
                st.rprimd.append( np.asarray([float(vec[0])*mul, float(vec[1])*mul, float(vec[2])*mul]) )

            st.nznucl = []
            for z in f.readline().split():
                st.nznucl.append( int(z)  )

            type_of_coordinates = f.readline()

            st.xred = []


            elements_list = []

            if "Car" in type_of_coordinates:
                # print "Warning! may be obsolete!!! and incorrect"
                # f.write("Cartesian\n")
                # for xcart in zxcart:
                #     for x in xcart:
                #         f.write(str(x[0]*to_ang)+" "+str(x[1]*to_ang)+" "+str(x[2]*to_ang))
                #         f.write("\n")
                raise RuntimeError
                
            elif "dir" in type_of_coordinates or 'Dir' in type_of_coordinates:
                for nz in st.nznucl:



                    for i in range(nz):
                        vec = f.readline().split()
                        st.xred.append( np.asarray([float(vec[0]), float(vec[1]), float(vec[2])]) )




                        if len(vec) == 4: # elements may be added by pymatgen
                            # print_and_log("Probably elements names are added at the end of coordinates, trying to read them")
                            if vec[3] not in elements_list:
                                elements_list.append(vec[3])
                        




                st.xcart = xred2xcart(st.xred, st.rprimd)

            elif 'None' in type_of_coordinates:
                pass

            else:
                print_and_log("Error! The type of coordinates should be 'car' or 'dir' ")
                raise NameError



            if 'Species order:' in name:
                print_and_log('I detect that input file was generated by cif2cell\n')
                name = name.split(':')[-1]


            if not elements_list:
                elements_list = name.split('!')[0].strip().split()
                print_and_log('I take elements from the first line, The line is '+str(name.split('!'))+' you could use ! to add comment after name+\n')
            


            else:
                print_and_log("Elements names have been taken from the end of coordinates, pymatgen file?\n")



            st.znucl = []
            for elname in elements_list:
                st.znucl.append( element_name_inv(elname) )
            # print_and_log('znucl is ')



            st.natom = len(st.xred)

            st.ntypat = len(st.znucl)

            st.typat = []
            for i, nz in enumerate(st.nznucl):
                for j in range(nz):
                    st.typat.append(i+1)

            #Determine reciprocal vectors
            st.recip = []
            st.vol = np.dot( st.rprimd[0], np.cross(st.rprimd[1], st.rprimd[2])  ); #volume
            #print vol
            st.recip.append(   np.cross( st.rprimd[1], st.rprimd[2] )   )
            st.recip.append(   np.cross( st.rprimd[2], st.rprimd[0] )   )
            st.recip.append(   np.cross( st.rprimd[0], st.rprimd[1] )   )
            for i in 0,1,2:
                st.recip[i] =  st.recip[i] * 2 * math.pi / st.vol;

                # if hasattr(self.init, 'vel'):
                #     print "I write to POSCAR velocity as well"
                #     f.write("Cartesian\n")
                #     for v in self.init.vel:
                #         f.write( '%.12f %.12f %.12f\n'%(v[0]*to_ang, v[1]*to_ang, v[2]*to_ang) )

        print_and_log('The following Z were read = '+ str(self.init.znucl)+'\n')


        print_and_log( "POSCAR was read\n")
        return
    

    def check_kpoints(self):
        to_ang_local = to_ang
        try:
            if "Ang" in self.len_units:
                to_ang_local = 1
                #print "units angs"
        except AttributeError:
            print_and_log("Warning! no len_units for "+self.name+" calculation, I use Bohr \n") 
        N_from_kspacing = []
        # print self.set.vasp_params['KSPACING']
        # print 
        for i in 0, 1, 2:
            N_from_kspacing.append( math.ceil( (np.linalg.norm(self.init.recip[i]) / to_ang_local) / self.set.vasp_params['KSPACING']) )
        #print "Vector length is:", (np.linalg.norm(self.rprimd[0]), "Bohr"

        if self.set.kpoints_file  == False:#self.set.vasp_params['KSPACING']:
            N = N_from_kspacing

        if self.set.ngkpt:
            N = self.set.ngkpt
        else:
            N = N_from_kspacing


        print_and_log("Kpoint   mesh is: "+str(N)+'\n' )
        print_and_log("The actual k-spacings are "+str(self.calc_kspacings(N) )+'\n' )
        return N_from_kspacing


    def write_structure(self, name_of_output_file, type_of_coordinates = 'dir', option = None, prevcalcver = None, path = None, state = 'init'):
        """Generates POSCAR file
           type_of_coordinates - 'direct' (xred) or 'cartesian' (xcart)
           option -inheritance option
           prevcalcver - ver of first calc in calc list; for first None
           state - 'init' or 'end' 
        """
        #units
        try:
            if "ang" in self.len_units or "Ang" in self.len_units: 
                global to_ang; to_ang = 1.0; print_and_log("Conversion multiplier to_ang is "+str(to_ang) )
        except AttributeError:
            pass

        if option == 'inherit_xred' and 'car' in type_of_coordinates: raise RuntimeError 

        if option == 'inherit_xred' and prevcalcver: type_of_coordinates = 'None' # do not write xred or xcart if they will be transfered on cluster
        
        if path == None: path = self.dir
        
        if state == 'init':
            st  = self.init
        elif state == 'end':
            st  = self.end
        else: 
            raise RuntimeError 
        
        rprimd = st.rprimd
        xred   = st.xred
        xcart  = st.xcart
        typat = st.typat  
        znucl = st.znucl

        
        # print 
        """1. Generate correct nznucl and zxred and zxcart"""
        zxred  = [[] for i in znucl]
        zxcart = [[] for i in znucl]
        # nznucl = []
        if len(typat) != len(xred) or len(xred) != len(xcart):
            raise RuntimeError
        for t, xr, xc in zip(typat, xred, xcart):
            # print "t ", t, xr
            zxred[ t-1].append(xr)
            zxcart[t-1].append(xc)
        
        nznucl = [len(xred) for xred in zxred]

        # print znucl, typat
        # print zxred
        # print nznucl
        if not os.path.exists(path):
            log.write( runBash("mkdir -p "+path) )

        with open(path+'/'+name_of_output_file,'w') as f:
            """Writes structure (POSCAR) in VASP format """
            f.write(self.name)
            
            f.write("\n1\n")
            
            for i in 0, 1, 2:
                f.write('  %18.16f %18.16f %18.16f'%(rprimd[i][0]*to_ang,rprimd[i][1]*to_ang,rprimd[i][2]*to_ang))
                f.write("\n")

            for n in nznucl:    
                f.write(str(n)+' ')

            f.write('\n')

            if "car" in type_of_coordinates:
                print "Warning! may be obsolete!!! and incorrect"
                f.write("Cartesian\n")
                for xcart in zxcart:
                    for x in xcart:
                        f.write(str(x[0]*to_ang)+" "+str(x[1]*to_ang)+" "+str(x[2]*to_ang))
                        f.write("\n")

                
                if hasattr(self.init, 'vel'):
                    print "I write to POSCAR velocity as well"
                    f.write("Cartesian\n")
                    for v in self.init.vel:
                        f.write( '  {:18.16f}  {:18.16f}  {:18.16f}\n'.format(v[0]*to_ang, v[1]*to_ang, v[2]*to_ang) )

            elif "dir" in type_of_coordinates:
                f.write("Direct\n")
                for xred in zxred:
                    for x in xred:
                        f.write("  {:18.16f}  {:18.16f}  {:18.16f}\n".format(x[0], x[1], x[2]) )
            elif 'None' in type_of_coordinates:
                pass

            else:
                print_and_log("Error! The type of coordinates should be 'car' or 'dir' ")
                raise NameError
        f.close()
        print_and_log( "POSCAR was generated\n")
        return



    def add_potcar(self):
        """Should be run for the first calculation only"""
        #Add POTCAR
        #try:
        #    shutil.copy2(self.set.potdir[0]+'/POTCAR',self.dir)
        #except IOError:
        #    runBash("gunzip "+self.set.potdir[0]+"/POTCAR.Z")

        path_to_potcar = self.dir+'/POTCAR'
        potcar_files   = ""

        for z in self.init.znucl:
            potcar_files += path_to_potentials+'/'+self.set.potdir[ int(z) ]+"/POTCAR "

        runBash("cat "+potcar_files+" >"+path_to_potcar)

        print_and_log( "POTCAR files: "+potcar_files+"\n")        
        return


    def calculate_nbands(self):
        """Should be run after add_potcar()"""
        #1 add additional information to set

        path_to_potcar = self.dir+'/POTCAR'
        self.init.zval = []
        # print path_to_potcar
        for line in open(path_to_potcar,'r'):
            if "ZVAL" in line:
                # print line
                self.init.zval.append(float(line.split()[5]))
        
        try: self.set.add_nbands
        except AttributeError: self.set.add_nbands = None

        if self.set.add_nbands != None:
            tve =0
            for i in range(self.init.ntypat):
                # print self.init.zval
                tve += self.init.zval[i] * self.init.nznucl[i]
            self.nbands = int ( round ( math.ceil(tve / 2.) * self.set.add_nbands ) )
            self.set.vasp_params['NBANDS'] = self.nbands
        return


    def actualize_set(self):
        """
        Makes additional processing of set parameters, which also depends on calculation
        """


        #check if some parameters should be filled according to number of species
        #make element list
        el_list = [element_name_inv(el) for el in self.init.znucl]
        vp = self.set.vasp_params


        if 'LDAU' in vp and vp['LDAU']:        
            for key in ['LDAUL', 'LDAUU', 'LDAUJ']:
                if set(vp[key].keys()).isdisjoint(set(el_list)): #no common elements at all
                    print_and_log('\n\n\nAttention! The '+str(key)+' doesnt not contain values for your elements! Setting to zero\n\n\n')
                    # raise RuntimeError

                new = []
                for el in el_list:
                    if el in vp[key]:
                        val = vp[key][el]
                    else:
                        if key == 'LDAUL':
                            val = -1
                        else:
                            val =  0

                    new.append(val)
                vp[key] = new


        """Process magnetic moments"""
        # print dir(self.set)
        # print self.set.vasp_params

        # print self.set.vasp_params['MAGMOM']
        if hasattr(self.init, 'magmom') and self.init.magmom:
            print_and_log('Magnetic moments are determined from self.init.magmom:',self.init.magmom, imp = 'y')

        elif hasattr(self.set, 'magnetic_moments') and self.set.magnetic_moments:
            print_and_log('Magnetic moments are determined using siman key "magnetic_moments" and corresponding dict in set\n')
            mag_mom_other = 0.6 # magnetic moment for all other elements
            magmom = []
            for iat in range(self.init.natom):
                typ = self.init.typat[iat]
                el  = el_list[typ-1]
                if el in self.set.magnetic_moments:
                    magmom.append(self.set.magnetic_moments[el])
                else:
                    magmom.append(mag_mom_other)
            

            #convert magmom to vasp ordering

            zmagmom = [[] for x in xrange(0,self.init.ntypat)]

            # print zmagmom

            for t, m in zip(self.init.typat, magmom):
                # print "t, m = ", t, m
                zmagmom[t-1].append(m)
                # print t-1, zmagmom[3]

            # print 'sdfsdf', zmagmom[3], 
            magmom = [m for mag in zmagmom for m in mag  ]
            # sys.exit()



            # print el_list
            # print self.init.typat
            # print magmom
                
            self.set.vasp_params['MAGMOM'] = magmom

            #check possible antiferromagnetic configurations:
            spec_mom_is = []
            for i, m in enumerate(magmom):
                if m != mag_mom_other: #detected some specific moment
                    spec_mom_is.append(i)

            if len(spec_mom_is) % 2 == 0 and len(spec_mom_is) > 0:
                print_and_log('Number of elements is even! trying to find all antiferromagnetic orderings:', imp = 'y')
                ns = len(spec_mom_is); 
                number_of_ord = math.factorial(ns) / math.factorial(0.5 * ns)**2
                if number_of_ord > 100:
                    print_and_log('Attention! Number of orderings is more than 100 - no sense to check them; May be to check several random?', imp = 'y')
                else:

                    ls = [0]*len(spec_mom_is)
                    # print ls
                    orderings = []
                    
                    def spin(ls, i):
                        """
                        Find recursivly all possible orderings
                        """
                        for s in 1,-1:
                            ls[i] = s
                            if i < len(ls)-1:
                                spin(ls, i+1)
                            else:
                                if sum(ls) == 0:
                                    orderings.append(copy.deepcopy(ls) )            
                        return

                    spin(ls, 0)

                    mag_orderings = []
                    mag_orderings.append(magmom)
                    for j, order in enumerate(orderings):
                        # print order

                        new_magmom = copy.deepcopy(magmom)
                        for i, s in zip(spec_mom_is, order):
                            # print i
                            new_magmom[i] = s * magmom[i]
                        print j, new_magmom
                        mag_orderings.append(new_magmom)

                    # print orderings
                    print_and_log('Total number of orderings is ', len(orderings),imp = 'y')
                    
                    if self.calc_method and 'afm_ordering' in self.calc_method:
                        self.magnetic_orderings = mag_orderings

                    


            self.init.magmom = magmom

        
        elif 'MAGMOM' in self.set.vasp_params and self.set.vasp_params['MAGMOM']: #just add * to magmom tag if it is provided without it
            print_and_log('Magnetic moments from vasp_params["MAGMOM"] are used\n')
            
            if "*" not in self.set.vasp_params['MAGMOM']:
                self.set.vasp_params['MAGMOM'] = str(natom) +"*"+ self.set.vasp_params['MAGMOM']
        




        return

    def make_incar_and_copy_all(self, update):
        """Makes Incar file for current calculation and copy all
        TO DO: there is no need to send all POSCAR files; It is enothg to send only one. However for rsync its not that crucial
        """
        #print "Begin make---------------------------------------------"
        
        
        #Generate incar
        vp = self.set.vasp_params
        natom = self.init.natom
        #please make consistent
        # print vp

        with open(self.dir+"INCAR",'w') as f:
            f.write( 'SYSTEM = %s\n\n'%(self.des) )
            for key in sorted(self.set.vasp_params):

                if key == 'MAGMOM' and hasattr(self.init, 'magmom') and self.init.magmom: #
                    f.write('MAGMOM = '+list2string(self.init.magmom)+"\n") #magmom from geo file has higher preference
                    # sys.exit()
                    continue
                
                if self.set.vasp_params[key] == None:
                    continue


                if type(self.set.vasp_params[key]) == list:
                    lis = self.set.vasp_params[key]
                    f.write(key + " = " + ' '.join(['{:}']*len(lis)).format(*lis) + "\n")
               
                else:
                    f.write(key+" = "+str( self.set.vasp_params[key] ) +"\n")
           

            f.write("\n")

        print_and_log( "INCAR was generated\n")

        #Generate KPOINTS
        d = self.dir
        if self.set.kpoints_file:
            if self.set.kpoints_file == True:
                print_and_log( "You said to generate KPOINTS file.\n")
                self.calc_kspacings()
                #Generate kpoints file

                #
                if self.set.ngkpt:
                    nk1 = self.set.ngkpt[0]; nk2 = self.set.ngkpt[1]; nk3 = self.set.ngkpt[2];
                    print_and_log( "Attention! ngkpt was used for kpoints file\n")

                else:
                    print_and_log( "Attention! ngkpt for kpoints file are created from kspacing; ngkpt is empty\n")
                    N = self.check_kpoints()
                    self.set.ngkpt = N
                    nk1 = N[0]; nk2 = N[1]; nk3 = N[2]
                
                with open(self.dir+"KPOINTS",'w') as f:

                    f.write("Automatic Mesh\n") #Comment
                    f.write("0 \n")#Number of points; 0-Auto
                    if self.set.vasp_params["KGAMMA"] == ".TRUE.": 
                        f.write("Gamma\n")
                    else: 
                        f.write("Monkhorst Pack\n")
                    f.write('%i %i %i \n'%(nk1, nk2, nk3) )
                    f.write("0 0 0\n") # optional shift

                print_and_log( "KPOINTS was generated\n")
            else:
                shutil.copyfile(self.set.kpoints_file, self.dir+"KPOINTS")
                print_and_log( "KPOINTS was copied from"+self.set.kpoints_file+"\n")



            list_to_copy = [d+"INCAR",d+"POTCAR",d+"KPOINTS"]

        else:
            print_and_log( "This set is without KPOINTS file.\n")

            list_to_copy = [d+"INCAR",d+"POTCAR"]

            #N = []
            #for i in 0, 1, 2:
                #N.append( math.ceil( (np.linalg.norm(self.recip[i]) / to_ang) / self.set.vasp_params['KSPACING']) )
            #print "Vector length is:", (np.linalg.norm(self.rprimd[0]), "Bohr"
            #print_and_log("Kpoint   mesh is: "+str(N) )
            #print_and_log("The actual k-spacings is "+str(self.calc_kspacings(N) ) )
        #Copy section


        list_to_copy.extend( glob.glob(d+'/*POSCAR*') )
        list_to_copy.extend( glob.glob(d+'/*.run*') )
        string_of_paths = ""

        for nf in list_to_copy:
            string_of_paths += nf+" " #Compose all files for copy in one string

        
        if "up" in update: #Copy to server
            print_and_log("Files to copy: "+string_of_paths+"\n")

            # temp_dir = os.path.dirname(project_path_cluster+self.dir)
            runBash('ssh '+ cluster_address+' "mkdir -p '+project_path_cluster+self.dir+'"') #create directory
            log.write( runBash("rsync -zave ssh "+string_of_paths+" "+cluster_address+":"+project_path_cluster+self.dir)+"\n" )
        #print "End make---------------------------------------------\n\n"



    def write_sge_script(self, input_geofile = "header", version = 1, option = None, prevcalcver = None, savefile = None, schedule_system = None):
        """Without arguments writes header, else adds sequence of calculatios
            option - 'inherit_xred' - control inheritance, or 'master' - run serial on master 
            prevcalcver - ver of previous calc; for first none
            savefile - 'cdawx', where c-charge, d-dos, a- AECCAR, w-wavefile, x-xml
            schedule_system - type of job scheduling system:'PBS', 'SGE', 'SLURM'

        """


        # print 'Starting write_sge()', input_geofile

        def prepare_input(prevcalcver = None, option = None, input_geofile = None, name_mod_prev = '', write = True, curver = None):
            """1. Input files preparation

                curver - current version
            """  


            if write:
                # if not 'only_neb' in self.calc_method:
                precont = str(prevcalcver)+name_mod_prev+'.CONTCAR' #previous contcar
                if option == 'inherit_xred' and prevcalcver:
                    f.write('grep -A '+str(self.init.natom)+ ' "Direct" '+precont+' >> '+input_geofile+ ' \n')

                if option == 'continue':
                    ''
                    precont = str(curver)+name_mod_prev+'.CONTCAR ' #previous contcar
                    preout  = str(curver)+name_mod_prev+'.OUTCAR ' #previous outcar
                    f.write("cp "+precont+" POSCAR  # inherit_option = continue\n")
                    f.write("cp "+preout+'prev.'+preout+" # inherit_option = continue\n")

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

                for i, u in enumerate(set_LDAUU_list):
                    if u == 0: continue
                    u_step = np.linspace(0, u, self.set.u_ramping_nstep)[u_ramp_step]
                    u_step = np.round(u_step, 1)
                    # new_LDAUU_list[i] = value
                    new_LDAUU_list[i] = u_step


                new_LDAUU = 'LDAUU = '+' '.join(['{:}']*len(new_LDAUU_list)).format(*new_LDAUU_list)
                
                command = "sed -i.bak '/LDAUU/c\\" + new_LDAUU + "' INCAR\n"


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
                else:
                    f.write(parrallel_run_command +" >"+name+".log\n")
                
                f.write("sleep 20\n")
            return


        def mv_files_according_versions(savefile, v, name_mod = '', write = True, rm_chg_wav = True):    
            """3. Out files saving block
                
                rm_chg_wav - if True than CHGCAR and WAVECAR are removed

            """   
            
            if write:
                pre = v + name_mod

                if "o" in savefile:

                    f.write("mv OUTCAR "          + v + name_mod +  ".OUTCAR\n")
                    f.write("mv CONTCAR "         + v + name_mod +  ".CONTCAR\n")
                    f.write("cp INCAR "           + v + name_mod +  ".INCAR\n")
                
                if "c" in savefile:
                    chg  = pre + '.CHG'
                    f.write("mv CHG "+chg+"\n")
                    f.write("gzip "+chg+"\n")
                else:
                    f.write("rm CHG \n") #file can be used only for visualization






                if "d" in savefile:
                    f.write("mv DOSCAR "      + v + name_mod + ".DOSCAR\n")
                
                if "a" in savefile:
                    f.write("mv AECCAR0 "     + v + name_mod + ".AECCAR0\n")
                    f.write("mv AECCAR2 "     + v + name_mod + ".AECCAR2\n")
                
                if 'x' in savefile:
                    f.write("mv vasprun.xml " + v + name_mod + ".vasprun.xml\n")
               
                if 'w' in savefile:
                    f.write("mv WAVECAR "     + v + name_mod + ".WAVECAR\n")
                # else:
                #     f.write("rm WAVECAR\n")


                if rm_chg_wav:
                    f.write("rm CHGCAR \n") #file is important for continuation
                    f.write("rm WAVECAR\n")


            return


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
            parrallel_run_command = "mpirun -x PATH vasp"
        elif schedule_system == 'PBS':
            parrallel_run_command = "mpiexec --prefix /home/aleksenov_d/mpi/openmpi-1.6.3/installed vasp"
        
        elif schedule_system == 'SLURM':
            parrallel_run_command = "prun /opt/vasp/bin/vasp5.4.1MPI"
        else:
            raise RuntimeError


        run_name = self.dir+self.id[0]+"."+self.id[1]+'.run'        
        job_name = self.id[0]+"."+self.id[1]
        neb_flag = ('neb' in self.calc_method or 'only_neb' in self.calc_method)

        if input_geofile == "header":
            with open(run_name,'w') as f:
                if schedule_system == 'SGE':
                    f.write("#!/bin/tcsh   \n")
                    f.write("#$ -M aksenov@mpie.de\n")
                    f.write("#$ -m be\n")
                    f.write("#$ -S /bin/tcsh\n")
                    f.write("#$ -cwd \n")
                    f.write("#$ -R y \n")
                    f.write("#$ -o "+self.dir+" -j y\n\n")

                    f.write("cd "+self.dir+"\n")
                    f.write("module load sge\n")
                    f.write("module load vasp/parallel/5.2.12\n\n")


                # import random
                # foo = ['01', '02', '03', '04', '05', '06', '07', '08', '10', '12', '13', '15', '16', '17', '18', '19', '20']
                # print(random.choice(foo)

                if schedule_system == 'PBS':
                    f.write("#!/bin/bash   \n")
                    f.write("#PBS -N "+job_name+"\n")
                    f.write("#PBS -l walltime=99999999:00:00 \n")
                    f.write("#PBS -l nodes=1:ppn="+str(header.corenum)+"\n")
                    f.write("#PBS -r n\n")
                    f.write("#PBS -j eo\n")
                    f.write("#PBS -m bea\n")
                    f.write("#PBS -M dimonaks@gmail.com\n")
                    f.write("cd $PBS_O_WORKDIR\n")
                    f.write("PATH=/share/apps/vasp/bin:/home/aleksenov_d/mpi/openmpi-1.6.3/installed/bin:/usr/bin:$PATH \n")
                    f.write("LD_LIBRARY_PATH=/home/aleksenov_d/lib64:$LD_LIBRARY_PATH \n")

                    # f.write("cd "+self.dir+"\n")

                    # f.write("module load sge\n")
                    # f.write("module load vasp/parallel/5.2.12\n\n")

                if schedule_system == 'SLURM':
                    f.write("#!/bin/bash   \n")
                    f.write("#SBATCH -J "+job_name+"\n")
                    f.write("#SBATCH -t 250:00:00 \n")
                    f.write("#SBATCH -N 1\n")
                    f.write("#SBATCH -n "+str(header.corenum)+"\n")
                    f.write("#SBATCH -o "+cluster_home+self.dir+"sbatch.out\n")
                    f.write("#SBATCH -e "+cluster_home+self.dir+"sbatch.err\n")
                    f.write("#SBATCH --mem-per-cpu=7675\n")
                    f.write("#SBATCH --mail-user=d.aksenov@skoltech.ru\n")
                    f.write("#SBATCH --mail-type=END\n")
                    f.write("cd ~/"+self.dir+"\n")
                    f.write("export OMP_NUM_THREADS=1\n")

                    f.write("module add prun/1.0\n")
                    f.write("module add intel/16.0.2.181\n")
                    f.write("module add impi/5.1.3.181\n")




                # f.write("rm WAVECAR\n")                
                # f.write("rm WAVECAR\n")                

            # f.close()




        elif input_geofile != "footer": #control part of script

            
            v = str(version)


            self.associated_outcars = []




            with open(run_name,'a') as f: #append information about run

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

                #experimental preliminary non-magnetic run
                if 0:
                    if self.set.vasp_params['ISPIN'] == 2:
                        print_and_log('Magnetic calculation detected; For better convergence',
                         'I add first non-magnetic run', imp = 'Y')
                        write = True
                        name_mod_last = '.'+'NM'
                        name_mod = '.NM'

                        if write: 
                            f.write("#Preliminary non-magnetic run:\n")  
                        prepare_input(prevcalcver = prevcalcver, option = option,
                         input_geofile = input_geofile, name_mod_prev = name_mod_last, write = write, curver = version)

                        update_incar(parameter = 'ISPIN', value = 1, write = write) #
                        
                        run_command(option = option, name = self.name+name_mod, parrallel_run_command = parrallel_run_command, write = write)

                        if write:
                            f.write("cp CONTCAR POSCAR  #prepare for basic run\n")
                            write_poscar = False  

                        mv_files_according_versions('co', v, name_mod = name_mod, write = write, rm_chg_wav = False)

                        update_incar(parameter = 'ISPIN', value = 2, write = write) #




                if 'u_ramping' in self.calc_method:

                    if write: 
                        f.write("#U-ramping run:\n")  



                    # name_mod_last = '.'+name_mod_U_last()
                    name_mod_last = '.'+'U00' #since u_ramp starts from u = 00, it is more correct to continue from 00

                    # print name_mod_last
                    # print 'prevcalver', prevcalcver

                    if write: 
                        f.write("rm CHGCAR\n")                

                    prepare_input(prevcalcver = prevcalcver, option = option,
                     input_geofile = input_geofile, name_mod_prev = name_mod_last, write = write_poscar, curver = version)

                    for i_u in range(self.set.u_ramping_nstep):


                        u = update_incar(parameter = 'LDAUU', u_ramp_step = i_u, write = write)

                        name_mod   = '.U'+str(u).replace('.', '')




                        
                        run_command(option = option, name = self.name+name_mod, parrallel_run_command = parrallel_run_command, write = write)


                        if write: 
                            f.write("cp CONTCAR POSCAR\n")                


                        mv_files_according_versions('o', v, name_mod = name_mod, write = write, rm_chg_wav = False)
                    
                        self.associated_outcars.append( v + name_mod +  ".OUTCAR"  )

                        print 'write_sge(): as_outcars=', self.associated_outcars


                    mv_files_according_versions(savefile = 'c', v=v, name_mod = name_mod) #save more files for last U
                    analysis_script(write = write)
                    # print self.associated


                elif 'afm_ordering' in self.calc_method:

                    #Comment - inherit_xred option is not available here
                    f.write("rm CHGCAR\n")                
                    if not savefile: savefile = 'o'

                    for i, magmom in enumerate(self.magnetic_orderings):

                        name_mod   = '.AFM'+str(i)

                        update_incar(parameter = 'MAGMOM', value = magmom)

                        prepare_input(prevcalcver = prevcalcver, option = option, input_geofile = input_geofile)
                        
                        run_command(option = option, name = self.name+name_mod, parrallel_run_command = parrallel_run_command)

                        mv_files_according_versions(savefile, v, name_mod = name_mod)
                    
                        self.associated_outcars.append( v + name_mod +  ".OUTCAR"  )

                    analysis_script()
                



                else:
                    
                    if not savefile: savefile = 'cdox'

                    if write: 
                            f.write("#Basic run:\n")  

                    prepare_input(prevcalcver = prevcalcver, option = option, 
                        input_geofile = input_geofile, write = write_poscar, curver =version)

                    run_command(option = option, name = self.name, parrallel_run_command = parrallel_run_command, write = write)

                    mv_files_according_versions(savefile, v, write = write)




                # f.write("\n")


        else: #footer
            with open(run_name,'a') as f: #append information about run

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


                    if 'u_ramping' in self.calc_method:
                        u = update_incar(parameter = 'LDAUU', u_ramp_step = self.set.u_ramping_nstep-1, write = False)
                        name_mod   = '.U'+str(u).replace('.', '')
                        # name_mod_last = name_mod_U_last()+'.'
                        name_mod_last = '.'+'U00' #since u_ramp starts from u = 00, it is more correct to continue from 00
                    
                    else:
                        name_mod_last = ''
                        name_mod   = ''

                    start = '1'+name_mod+'.OUTCAR '
                    final = '2'+name_mod+'.OUTCAR '

                    f.write("\n\n#Starting NEB script \n")

                    if option and 'continue' in option:
                        prevout = name_mod_last+'OUTCAR '

                        for n_st in subfolders:
                            f.write('cp '+n_st+'/'+prevout+n_st+'/'+'prev.'+prevout+'  # inherit_option = continue\n' )



                    f.write('~/tools/vts/nebmake.pl '+ start.replace('OUT','CONT') + final.replace('OUT','CONT') + nim_str +' \n')
                    
                    f.write('cp '+start+ '00/OUTCAR\n')
                    
                    if nim+1 < 10: 
                        nim_plus_one_str = '0'+str(nim+1)

                    f.write('cp '+final + nim_plus_one_str + '/OUTCAR\n' )

                    update_incar(parameter = 'IMAGES', value = nim)










                    if 'u_ramping' in self.calc_method:



                        
                        for i_u in range(self.set.u_ramping_nstep):


                            u = update_incar(parameter = 'LDAUU', u_ramp_step = i_u)

                            name_mod   = 'U'+str(u).replace('.', '')

                            
                            run_command(option = option, name = self.name+'.images'+nim_str+'.'+name_mod, 
                                parrallel_run_command = parrallel_run_command, write = True)


                            for n_st in subfolders:

                                f.write('cp '+n_st+'/CONTCAR '+n_st+'/POSCAR'+'\n' )
                                f.write('cp '+n_st+'/OUTCAR  '+n_st+'/'+name_mod+'.OUTCAR'+'\n' )

                        
                            # self.associated_outcars.append( v + name_mod +  ".OUTCAR"  )
                    




                    else:

                        run_command(option = option, name = (self.name+'.images'+nim_str), 
                        parrallel_run_command = parrallel_run_command, write = True)










                    f.write('export PATH=$PATH:/home/aksenov/tools/gnuplot/bin/ \n')
                    f.write('~/tools/vts/nebresults.pl  \n')
                    f.write('find . -name WAVECAR -delete\n')
                    f.write('find . -name PROCAR -delete\n')
                    # for n in range


                #clean 
                f.write('rm PROCAR DOSCAR OSZICAR PCDAT REPORT XDATCAR vasprun.xml\n')


            runBash('chmod +x '+run_name)


        if hasattr(self, 'associated_outcars') and  self.associated_outcars:
            out = self.associated_outcars[-1]
        else:
            out = None
        print 'write_sge() out=', out
        
        return  out#return OUTCAR name
    



    def make_run(self, schedule_system = None):
        """Generate run file

        INPUT:
            schedule_system - 
        """
        run_name = self.dir+self.id[0]+"."+self.id[1]+'.run'

        with open('run','a') as f:
            if schedule_system == 'SGE':
                #'qsub -pe 'mpi*' NCORES -l CLUSTER_TAG script.parallel.sh' for mpi-jobs which should run on CLUSTER_TAG (cmmd or cmdft)
                #IMPORTANT: NCORES must be a multiple of 8(24) on cmdft(cmmd). 
                f.write("qsub -pe 'mpi*' "+str(header.corenum)+" "+header.queue+" "+run_name+"\n") #str(self.set.np) #-l cmmd
                f.write('sleep 5\n')
                # runBash('chmod +x run')
            if schedule_system == 'PBS':
                f.write("cd "+self.dir+"\n")
                f.write("qsub "+self.id[0]+"."+self.id[1]+'.run'+"\n") 
                f.write("cd -\n")
                f.write('sleep 5\n')                        
            
            if schedule_system == 'SLURM':
                f.write("sbatch -p AMG " + run_name+"\n") 

                



        self.state = "3. Ready for start"
        log.write("\nRun file created\n")     
        return



    def calc_kspacings(self,ngkpt = []):
        """Calculates reciprocal vectors and kspacing from ngkpt"""
        to_ang_local = to_ang
        try:
            if "Ang" in self.len_units:
                to_ang_local = 1
                #print "units angs"
        except AttributeError:
            print_and_log("Warning! no len_units for "+self.name+" calculation, I use Bohr \n")
        #Determine reciprocal vectors
        if 1: #Already calculated during reading of structure
            self.recip = []
            vol = np.dot( self.init.rprimd[0], np.cross(self.init.rprimd[1], self.init.rprimd[2])  ); #volume
            #print vol
            self.recip.append(   np.cross( self.init.rprimd[1], self.init.rprimd[2] )   )
            self.recip.append(   np.cross( self.init.rprimd[2], self.init.rprimd[0] )   )
            self.recip.append(   np.cross( self.init.rprimd[0], self.init.rprimd[1] )   )
            for i in 0,1,2:
                self.recip[i] =  self.recip[i] * 2 * math.pi / vol;
        #print self.recip
        self.kspacing = []
        if not ngkpt:
            ngkpt = self.set.ngkpt
        k = [0,0,0]
        if ngkpt:
            # print 'ngkpt are ', ngkpt
            # print 'recip are', self.recip
            for i in 0, 1, 2:
                a = np.linalg.norm( self.recip[i] ) / ngkpt[i] / to_ang_local
                self.kspacing.append(red_prec(a))
            k = self.kspacing
        #print "Spacings for %s is [%.2f, %.2f, %.2f]"%(self.set.ngkpt,k[0], k[1], k[2])
        return  k













    def read_results(self, load = '', out_type = '', voronoi = False, show = '', choose_outcar = None):

        """
        Download and Read VASP OUTCAR file

        ###INPUT:
            - load (str) - 'x' - download xml, o - download outcar and contcar
            - show (str) - print additional information
            - choose_outcar - see description in res_loop()


        ###RETURN:
            ?

        ###DEPENDS:


        """

        if choose_outcar:
            path_to_outcar = os.path.dirname(self.path["output"])+ '/'+ self.associated_outcars[choose_outcar-1]

        else:
            path_to_outcar  = self.path["output"]
        
        print 'classes: path to outcar', path_to_outcar
        path_to_contcar = path_to_outcar.replace('OUTCAR', "CONTCAR")
        path_to_xml     = path_to_outcar.replace('OUTCAR', "vasprun.xml")


        if self.calc_method  and not set(self.calc_method  ).isdisjoint(  ['u_ramping', 'afm_ordering']):
            # print self.associated_outcars
            # lor = self.associated_outcars[-1]
            # path_to_outcar  = self.dir + lor
            # path_to_contcar = self.dir + lor.replace('OUTCAR', 'CONTCAR')
            # path_to_xml     = self.dir + lor.replace('OUTCAR', 'vasprun.xml')         
            # print 'sdf', path_to_outcar, path_to_contcar, path_to_xml
            energies_str = runBash("ssh "+self.cluster_address+" cat "+self.dir+"ENERGIES")
            # print "ssh "+self.cluster_address+" cat "+self.dir+"ENERGIES"
            self.associated_energies = [float(e) for e in energies_str.split()]
            # self.u_ramping_u_values = np.arange(*self.u_ramping_list)
            # print 'associated_energies:', self.associated_energies
        
        outcar_exist   = False

        contcar_exist   = False





        if not os.path.exists(path_to_outcar):
            load = 'o'


 
        self.natom = self.init.natom

        """Copy from server """
        if 'o' in load:

            #reduce size of downloadable file by removing occupations: vasp 4 and 5
            command_reduce = """ssh {0:s} nbands=\`grep \\"NBANDS=\\" \{1:s} \| awk \\'{{print \$NF - 1}}\\'\`\; sed -i -e \\"/band No./,+\${{nbands}}d\\" \{1:s} """.format(
                self.cluster_address, self.project_path_cluster+path_to_outcar )
            runBash(command_reduce)


            files = [self.project_path_cluster+path_to_outcar, self.project_path_cluster+path_to_contcar]
            get_from_server(files = files, to = os.path.dirname(path_to_outcar)+'/',  addr = self.cluster_address)



        if 'x' in load:

            get_from_server(files = self.project_path_cluster+path_to_xml, to = os.path.dirname(path_to_outcar)+'/',  
                addr = self.cluster_address)






        if os.path.exists(path_to_contcar):
            contcar_exist   = True

        if os.path.exists(path_to_outcar):
            outcar_exist   = True

        """Start reading """
        if outcar_exist:
            s = runBash("grep 'General timing' "+path_to_outcar)
        else:
            s = 'no OUTCAR'
        if "g" in s:
            self.state = "4. Calculation completed."
        else: 
            self.state+=s
            outst = self.state


        if "4" in self.state:

            #war = runBash('grep "NPAR = approx SQRT( number of cores)" '+path_to_outcar) #remove warinig
            nw = runBash('sed -n "/NPAR = approx SQRT( number of cores)/=" '+path_to_outcar) #remove warinig
                    # print runBash("sed '/from  /, /to  /d' "+path_to_outcar)
                    #print "sed -e '15,34d' "+path_to_outcar
            #print 'grep "NPAR = approx SQRT( number of cores)" '+path_to_outcar
            #print nw
            tmp = path_to_outcar+".tmp"
            if nw:
                nw = int(nw)
                #remove warning
                runBash("sed '"+str(nw-11)+","+str(nw+8)+"d' "+path_to_outcar+">"+tmp+";mv "+tmp+" "+path_to_outcar)
                pass

            with open(path_to_outcar, 'r') as outcar:
                log.write("Start reading from "+ path_to_outcar+" \n")
                outcarlines = outcar.readlines()

            re_lengths = re.compile("length of vectors")
            re_eltime = re.compile("Elapsed time")
            re_nkpts = re.compile("NKPTS")

            i_line = 0
            iterat = 0
            self.mdstep = 1
            warnings = 0#""
            self.time = 0
            nscflist = [] ; mdstep_old = 1; niter_old = 0
            maxforce = []; average = [];  gstress =[]
            # mforce = []
            self.list_e_sigma0 = []
            self.list_e_without_entr = []
            try:
                self.end = copy.deepcopy(self.init) # below needed end values will be updated
            except:
                self.end = Structure()

            if not hasattr(self.end, "natom"): self.end.natom = self.natom
            #Structure() #create structure object with end values after calculation
            #self.end.typat = self.typat
            #self.end.znucl = self.znucl
            self.end.name = self.name+'.end'
            self.end.list_xcart = []
            self.energy = empty_struct()

            nsgroup = 1
            magnitudes = []
            tot_mag_by_atoms = [] #magnetic moments by atoms on each step
            tot_mag_by_mag_atoms = []
            #which atoms to use
            magnetic_elements = [26, 27, 28]
            #Where magnetic elements?
            zlist = [int(self.end.znucl[t-1]) for t in self.end.typat]
            # print zlist
            # i_mag_start = None
            # i_mag_end   = None
            ifmaglist = [] #np.array() #[False]*len(zlist)
            for i, z in enumerate(zlist): #
                if z in magnetic_elements:
                    ifmaglist.append(True)
                else:
                    ifmaglist.append(False)

            ifmaglist = np.array(ifmaglist)
            # print ifmaglist

            # sys.exit()



            # print self.set.vasp_params

            for line in outcarlines:

                #Check bands

                # if 'band No.' in line:
                #     kpoint = float(outcarlines[i_line-1].split()[1])
                #     lastocc = float(outcarlines[i_line+self.nbands].split()[2])
                #     lastbandno = outcarlines[i_line+self.nbands].split()[0]
                #     if lastocc > 0:
                #         print "Warning!!! at kpoint ", kpoint, " last band No. ",lastbandno, " is not empty ", lastocc




                if "TOO FEW BANDS" in line:
                    print_and_log("Warning! TOO FEW BANDS!!!\n\n\nWarning! TOO FEW BANDS!!!\n")



                #Check W(q)
                if 'operators is LMAX' in line:
                    lmax = int(line.split()[7])
                    # print 'lmax', lmax
                if "W(low)/X(q)" in line:
                    kk = 1; low = []; high = [];
                    while kk < 100:
                        if 'Optimization' in outcarlines[i_line + kk] or len(outcarlines[i_line + kk].split() ) != 7: break
                        # print 'line', outcarlines[i_line + kk]

                        low.append(  float(outcarlines[i_line + kk].split()[4]) )
                        high.append( float(outcarlines[i_line + kk].split()[5]) )
                        kk+=1


                    if any(v > 1e-3 for v in low+high):
                        print_and_log("W(q)/X(q) are too high, check output!\n")
                        print 'Low + high = ', low+high
                        print [v > 1e-3 for v in low+high]
                if "direct lattice vectors" in line:
                    for v in 0,1,2:
                        self.end.rprimd[v] = np.asarray( [float(ri) for ri in outcarlines[i_line+1+v].split()[0:3]   ] )
                    #print self.end.rprimd
                    #print self.rprimd
                if "POSITION" in line:
                    if not contcar_exist or out_type == 'dimer':
                        self.end.xcart = [] #clean xcart before filling
                        for i in range(self.init.natom):
                            #print outcarlines[i_line+1+i].split()[0:3] 
                            xcart = np.asarray ( 
                                        [   float(x) for x in outcarlines[i_line+2+i].split()[0:3]   ] 
                                    )
                            
                            self.end.xcart.append( xcart )
                            #self.end.xred.append ( xcart2xred( xcart, self.end.rprimd) )
                
                        if out_type == 'dimer':
                            self.end.list_xcart.append(self.end.xcart) #xcart at each step only for dimer

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

                if "TOTAL-FORCE" in line:
                    # Calculate forces here...
                    forces = []
                    magnitudes = []

                    for j in range(0,self.init.natom):
                        parts = outcarlines[i_line+j+2].split()
                        # print "parts", parts
                        x = float(parts[3])
                        y = float(parts[4])
                        z = float(parts[5])
                        forces.append([x,y,z])
                        magnitudes.append(math.sqrt(x*x + y*y + z*z))
                    average.append( red_prec( sum(magnitudes)/self.init.natom * 1000 ) )
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
                    try:                     self.end.vol = float(line.split()[4])
                    except ValueError: print_and_log("Warning! Cant read volume in calc "+self.name+"\n")
                    #print self.vol      

                if "generate k-points for:" in line: 
                    self.set.ngkpt = tuple(  [int(n) for n in line.split()[3:]]  )
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
                if  "-1/2 Hartree   DENC" in line:
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
                    self.list_e_sigma0.append(  self.energy_sigma0  )
                    self.list_e_without_entr.append(  self.e_without_entr  )




                if "free  energy   TOTEN  =" in line:
                    #self.energy = float(line.split()[4])
                    self.energy_free = float(line.split()[4]) #F free energy
                if re_lengths.search(line):
                    self.vlength = [red_prec( float(l),1000 ) for l in outcarlines[i_line + 1].split()[0:3]]
                    #print self.vlength
                if "in kB" in line:
                    self.stress = [float(i)*100 for i in line.split()[2:]]  # in MPa 
                if "Total  " in line:
                    self.intstress = [int(float(i)*1000) for i in line.split()[1:]] #stress in internal units; can be regarded as forces
                
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
                    self.NKPTS = float(line.split()[3])
                if "WARNING" in line:
                    warnings += 1#line


                if "Subroutine DYNSYM returns" in line:
                    nsgroup = line.split()[4]#number of space group operations



                if "Iteration" in line:
                    self.mdstep = int(line.split('(')[0].split()[2].strip())
                    iterat +=1
                    if mdstep_old != self.mdstep:
                        #print "Stress:", self.stress 
                        nscflist.append( niter ) # add to list number of scf iterations during mdstep_old
                    niter = int(line.split(')')[0].split('(')[-1].strip()) #number of scf iterations
                    mdstep_old = self.mdstep



                if 'magnetization (x)' in line:
                    mags = []
                    for j in range(self.init.natom):
                        mags.append( float(outcarlines[i_line+j+4].split()[4]) )
                    
                    tot_mag_by_atoms.append(np.array(mags))#[ifmaglist])
                    tot_mag_by_mag_atoms.append(np.array(mags)[ifmaglist])
                    # print tot_mag_by_atoms
                    # magnetic_elements
                    # ifmaglist





                i_line += 1
            #Check total drift
            max_magnitude = max(magnitudes)
            max_tdrift    = max(tdrift)
            self.maxforce = maxforce[-1][1]
            # if max_magnitude < self.set.toldff/10: max_magnitude = self.set.toldff
            # print 'magn', magnitudes
            # print 'totdr', tdrift
            # print 'max_magnitude', max_magnitude
            # try 
            if max_magnitude < self.set.tolmxf: max_magnitude = self.set.tolmxf
            #if any(d > 0.001 and d > max_magnitude for d in tdrift):
            if max_tdrift > 0.001 and max_tdrift > max_magnitude:
                
                #print_and_log( ("Total drift is too high! At the end one component is %0.f %% of the maximum force, check output!\n") %(maxdrift)  )
                pass
            #else: maxdrift = 
            # print magn
            """Try to read xred from CONCAR and calculate xcart"""

            #correction of bug; Take into account that VASP changes typat by sorting impurities of the same type.
            self.end.typat = []
            for i, nz in enumerate(self.init.nznucl):
                for j in range(nz):
                    self.end.typat.append(i+1)
            #correction of bug




            #print contcar_exist
            if contcar_exist:
                with open(path_to_contcar, 'r') as contcar:
                    
                    for line in contcar:
                        
                        if "Direct" in line:
                            self.end.xred = []
                            for i in range(self.natom):
                                xr = np.asarray ( [float(x) for x in contcar.next().split()] )
                                self.end.xred.append( xr )






                self.end.xcart = xred2xcart( self.end.xred, self.end.rprimd)
            else: 
                self.end.xred = xcart2xred( self.end.xcart, self.end.rprimd)

            #print "init pressure = ",self.extpress_init,"; final pressure =",self.extpress
            #print self.end.xred
            #self.vol = np.dot( self.rprimd[0], np.cross(self.rprimd[1], self.rprimd[2])  ); #volume
            nscflist.append( niter ) # add to list number of scf iterations during mdstep_old
            #print "Stress:", self.stress
            v = self.vlength
            s = self.stress
            yznormal = np.cross(self.init.rprimd[1], self.init.rprimd[2])
            #print yznormal
            #print np.cross( yznormal, np.array([1,0,0]) )
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


            #  Construct beatifull table
            #self.a1 = float(v[0])/2 ; self.a2 = float(v[1])/2/math.sqrt(0.75); self.c = float(v[2])  # o1b
            
            try:
                self.a = self.hex_a ; self.c = self.hex_c  # c1b
            except AttributeError:
                self.a  = 0; self.c = 0 #calculations with full relaxation
            if self.a == None or self.a == [None]:
                self.a  = 0; self.c = 0

            j = (15,12,7,7,8,9,9,5,5,20,5,20,8,12,20,8,5,8,8,25,8,4,3)

            d = "&"
            name = ("%s.%s.%s" % (self.id[0],self.id[1], str(self.id[2]) )).ljust(j[0])
            etot = ("%.4f" % ( self.energy_sigma0 )).center(j[1])
            # print self.a
            a = ("%.4f" %      ( self.a )      ).center(j[2])
            c = ("%.4f" %      ( self.c )      ).center(j[3])
            time = ("%.2f" % (self.time/3600.)    ).center(j[4])
            itertm = ("%.1f" % (self.time/1./iterat)    ).center(j[5])
            Nmd = ("%1i,%2i,%3i" % (self.mdstep, iterat/self.mdstep, iterat)    ).center(j[6])
            self.iterat = iterat
            War = ("%i" % (warnings)    ).center(j[7])
            #nbands = ("%i" % (self.set.vasp_params["NBANDS"])    ).center(j[8])
            #added = ("%.0f" % ( (self.set.add_nbands - 1) * 100 )    ).center(j[15])
            kmesh = ("%s" % (str(self.set.ngkpt) )    ).center(j[8])
            ks = self.calc_kspacings()
            kspacing = ("[%.2f,%.2f,%.2f]" % ( ks[0], ks[1], ks[2] )    ).center(j[9])
            ks1 = ("[%.2f]" % ( ks[0] )    ).center(j[9])
            nkpt = ("%i" % ( self.NKPTS)     ).center(j[10])
            istrs = ("[%5i,%5i,%5i] " % ( self.intstress[0],self.intstress[1],self.intstress[2]  )    ).center(j[11])
            strs = ("%.0f,%.0f,%.0f " % ( self.stress[0],self.stress[1],self.stress[2]  )    ).center(j[11])   
            eprs = ("%.0f" % (self.extpress)).center(j[12])
            tsm = ("%.0f" % (self.set.tsmear*1000)).center(j[13])
            entrr = ("%.3f" % (   (self.energy_free - self.energy_sigma0)/self.init.natom * 1000    )   ).center(j[14]) #entropy due to the use of smearing
            npar = ("%i" % (self.set.vasp_params["NPAR"])).center(j[16])
            lpl = ("%s" % (self.set.vasp_params["LPLANE"])).center(j[17])
            ecut = ("%s" % (self.set.vasp_params["ENCUT"]) ).center(j[18]) 

            lens = ("%.2f;%.2f;%.2f" % (v[0],v[1],v[2] ) ).center(j[19])
            r1 = ("%.2f" % ( v[0] ) ).center(j[19])            
            vol = ("%.1f" % ( self.end.vol ) ).center(j[20])
            nat = ("%i" % ( self.natom ) ).center(j[21])
            totd = ("%.0f" % (   max_tdrift/max_magnitude * 100      ) ).center(j[22])
            nsg = ("%s" % (     nsgroup     ) ).center(j[22])

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
                        
            # print self.end.xred[-1]
            #print outstring_kp_ec
            # print show
            # print 'force' in show
            # if not hasattr(show, '__iter__')
            #     show = [show]

            if 'fo' in show:
                # print "Maxforce by md steps (meV/A) = %s;"%(str(maxforce)  )
                print "\nMax. F. (meV/A) = \n%s;"%(np.array([m[1] for m in maxforce ])[-50:]  )
                # print "\nAve. F. (meV/A) = \n%s;"%(  np.array(average)  )
                # import inspect
                # print inspect.getargspec(plt.plot).args
                # print plt.plot.__doc__
                if 'p' in show[0]:
                    plt.plot(maxforce, )
                    plt.xlabel('MD step')
                    plt.ylabel('Max. force on atom (meV/$\AA$)')
                    plt.show()
            
            if 'en' in show:
                    maxf = [m[1] for m in maxforce ]
                    plt.plot(maxf, 1000*(np.array(self.list_e_sigma0)-self.energy_sigma0) , )
                    # plt.xlabel('MD step')
                    # plt.ylabel('Energy per cell (eV')
                    plt.xlabel('Max. force on atom (meV/$\AA$)')
                    plt.ylabel('Energy per cell relative to min (meV)')

                    plt.show()

            if 'mag' in show:
                # print_and_log
                # print 'Final mag moments for atoms:'
                # print np.arange(self.end.natom)[ifmaglist]+1
                # print np.array(tot_mag_by_atoms)
                print 'Dist from 1st atom to Fe atoms:, please make me more general'
                sur   = local_surrounding(self.end.xcart[0], self.end, n_neighbours = 4, control = 'atoms', 
                periodic  = True, only_elements = [26,])

                dist = np.array(sur[3]).round(2)
                numb = np.array(sur[2])
                for mag in tot_mag_by_atoms:
                    print mag[numb]

                print np.array(sur[3]).round(2), np.array(sur[2])+1

                self.tot_mag_by_atoms = tot_mag_by_atoms
                plt.plot(np.array(tot_mag_by_mag_atoms))
                plt.show()




            log.write("Reading of results completed\n\n")
            
            if   out_type == 'gbe'  : outst = outst_gbe
            elif out_type == 'e_imp': outst = outst_imp
            elif out_type == 'e_seg': outst = outst_seg            
            elif out_type == 'coseg': outst = outst_coseg            
            elif 'ecut' in out_type : outst = outst_ecut
            elif 'kp' in out_type   : outst = outst_kp
            elif 'ts' in out_type   : outst = outst_ts
            else: outst = outst_all
            #else: print_and_log("Uknown type of outstring\n")





        else:
            # print_and_log("Still no OUTCAR for mystery reason for", self.id)
            print_and_log('OUTCAR not finished for', self.id)
            # raise RuntimeError

        return outst







    def get_chg_file(self, filetype = 'CHGCAR'):
        #cl - object of CalculationVasp class
        path_to_chg = self.dir+str(self.version)+"."+filetype
        if not os.path.exists(path_to_chg): 
            print 'Charge file is downloading'
            log.write( runBash("rsync -zave ssh "+self.cluster_address+":"+self.project_path_cluster+path_to_chg+" "+self.dir)+'\n' ) #CHG
            print path_to_chg, 'was downloaded'
            
        return path_to_chg

    def get_file(self, filename):
        #cl - object of CalculationVasp class
        # filename - standart Vasp file name
        path_to_file = self.dir+str(self.version)+"."+filename
        if not os.path.exists(path_to_file): 
            log.write( runBash("rsync -zave ssh "+self.cluster_address+":"+self.project_path_cluster+path_to_file+" "+self.dir)+'\n' ) #CHG
            print path_to_file, 'was downloaded'
            
        return path_to_file


    def bader_analysis(self):
        #Make bader on server
        #assumes that bader is intalled
        v = str(self.version)
        path = project_path_cluster+self.dir
        CHG     = path+v+".CHG"
        AECCAR0 = path+v+".AECCAR0"
        AECCAR2 = path+v+".AECCAR2"
        CHGCAR_sum = path+v+".CHGCAR_sum"
        baderlog =  path+v+".bader.log"
        command1 = "cd "+path+"; ~/utils/chgsum.pl "+AECCAR0+" "+AECCAR2+"; "+\
        "mv CHGCAR_sum "+CHGCAR_sum+";"
        command2 = \
        "cd "+path+"; ~/utils/bader "+CHG+" -ref "+CHGCAR_sum+" > "+\
        v+".bader.log; mv ACF.dat "+v+".ACF.dat; mv AVF.dat "+v+".AVF.dat; mv BCF.dat "+v+".BCF.dat;"
        # print "ssh "+cluster_address+" '"+command1+"'"
        
        if runBash("ssh "+self.cluster_address+" '[ -e "+   CHGCAR_sum       +""" ] || echo "NO"     ;' """): #true if file not exists
            print  CHGCAR_sum, "not exist. try to calculate it "
            log.write( runBash("ssh "+self.cluster_address+" '"+command1+"'")+'\n' ) 
        
        if runBash("ssh "+self.cluster_address+" '[ -e "+   baderlog       +""" ] || echo "NO"     ;' """): #true if file not exists
            print  baderlog, "not exist. try to calculate Bader "
            log.write( runBash("ssh "+self.cluster_address+" '"+command2+"'")+'\n' ) 
        ACF = runBash("ssh "+self.cluster_address+" 'cat "+path+v+".ACF.dat"  +"'" )
        # print ACF
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
        print "Partial charge of impurity ", imp_partial_chg
        print "Partial charges of neibouring Ti atoms", " ".join("{:.2f}".format(m) for m in mat_partial_chg)
        print "Partial charge of matrix", sum(mat_partial_chg)
        
        print "Sum of mat and imp charges:", sum(mat_partial_chg)+imp_partial_chg

        return path_to_chg

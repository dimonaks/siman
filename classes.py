"""
The most important header. Should be always first
"""
from header import *
import header
from functions import (local_surrounding, image_distance, xred2xcart, xcart2xred, 
write_xyz, element_name_inv, write_lammps, calculate_voronoi, replic,
log_history)

import optparse
import re
import glob


import sys
import colorsys
import shelve



# from ase.utils.eos import EquationOfState
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = "stix"

# print header.calc















def words(fileobj):
    """Generator of words. However does not allow to use methods of list for returned"""
    for line in fileobj:
        for word in line.split():
            yield word

def read_vectors(token, number_of_vectors, list_of_words):
    """Returns the list of numpy vectors for the last match"""

    number_of_matches = list_of_words.count( token )
    if number_of_matches == 0: 
        #print_and_log("Warning token '"+token+"' was not found! return empty\n")
        return [None]

    if number_of_matches > 1:
        print_and_log("Warning token '"+token+"' was found more than one times\n")
        raise RuntimeError


    index = list_of_words.index(token, number_of_matches - 1 )     #Return the index of the last match
    #print list_of_words[index]
    list_of_vectors = []
    vector = np.zeros((3))
    for i in range(number_of_vectors):
        vector[0] = float(list_of_words[index + 1])
        vector[1] = float(list_of_words[index + 2])
        vector[2] = float(list_of_words[index + 3])
        index+=3
        list_of_vectors.append(vector.copy())
    return list_of_vectors


def read_list(token, number_of_elements, ttype, list_of_words):
    """Input is token to find, number of elements to read, type of elements and list of words, 
    where to search
    Returns the list of elements for the last match"""
    

    number_of_matches = list_of_words.count( token )


    
    #if number_of_elements == 0:        raise RuntimeError
    
    if number_of_matches > 1:
        print_and_log("Warning token '"+token+"' was found more than one times\n")
        raise RuntimeError

    if number_of_matches == 0 or number_of_elements == 0: 
        #print_and_log("Warning token '"+token+"' was not found or asked number of elements is zero! set to [None]\n")
        #if ttype == str:
        #    return ['']*number_of_elements
        #else:
        #    return [0]*number_of_elements
        return [None]

    try:
        index = list_of_words.index(token, number_of_matches - 1 )     #Return the index of the last match

    except ValueError: 
        print_and_log("Warning!, token "+token+" was not found. I return [None]!\n")
        return [None]
    
    index+=1 #the position of token value
    list_of_elements = []
    
    #define function dependig on type:

    if   ttype == int  : 
        def convert(a): return int(a)
    elif ttype == float: 
        def convert(a): return float(a)
    elif ttype == str  : 
        def convert(a): return str(a)
    
    #print list_of_words[index], type(list_of_words[index])
    if list_of_words[index] == "None"  : 
        def convert(a): return [None]
    
    #Make convertion
    for i in range(number_of_elements):
        
        list_of_elements.append(    convert(  list_of_words[index]   )     )
        index+=1


    return list_of_elements

def read_database(scratch = False):
    #2. Read database of calculations
    #global history; 
    #global struct_des;
    databasefile = 'calc.s'
    if scratch == True: databasefile =   '/scratch/aksenov/calc.s'
    
    log.write("\nLaunch at "+str( datetime.datetime.today() )+'\n')
    d = shelve.open(databasefile, protocol=1)
    try:
        calc = d['calc']; 
        conv = d['conv']; 
        varset = d['varset']; 
        header.history = d['history']
        #struct_des = d['struct_des']
        #print struct_des 
       #
    except KeyError: 
        try: calc = d['calc'] #dictionary of calculations
        except KeyError:
            log.write( "There is no database of calculations. I create new"); calc = {}
        try: conv = d['conv'] #dictionary of convergence lists
        except KeyError:
            log.write( "There is no dictionary of convergence lists. I create new"); conv = {}   
        try: varset = d['varset'] 
        except KeyError:
            log.write( "There is no dictionary of inputsets. I create new");  varset = {} 
        try: header.history = d['history'] 
        except KeyError:
            log.write( "There is still no history in database. The list is in header module ");  #history = [] 
        #try: struct_des = d['struct_des'] 
        #except KeyError:
            #log.write( "There is still no struct_des in database. The dict is global "); # struct_des = {} 

    d.close()
    #print history

    return calc,conv,varset,sys.getsizeof(d)
def write_database(calc,conv,varset,size_on_start):
    #size_on_finish = sys.getsizeof(dbase)
    #if size_on_finish != size_on_start:
    runBash("cp calc.s calc_copy.s") #create copy before writing
    d = shelve.open('calc.s', protocol=1) #Write database of calculations
    #print struct_des 
    d['calc'] = calc
    d['conv'] = conv
    d['varset'] = varset
    d['history'] = header.history
    #d['struct_des'] = struct_des 

    d.close()
    log.write("\nEnd of work at "+str(datetime.datetime.now())+'\n')
    log.close()
    with  open('history','w') as his:
        #print history
        for i in header.history:
            #print i
            his.write(i+"\n")
    print("\nDatabase was succesfully updated\n")
    return


class Structure():
    """This class includes only structure related information such as primitive vectors, coordinates, forces and so on"""
    def __init__(self):
        self.name = ""
        self.des = ''
        self.rprimd = [np.zeros((3)) for i in 0,1,2]
        self.xcart = []
        self.xred = []


class Calculation():
    """docstring for Calculation"""
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
                    print line; self.des = line.split('des ')[1]+';'
                
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



            self.init = Structure()
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







            self.state = "2.Geometry has been read"



        #file.close();

        print_and_log( "If no warnings, geometry has been succesfully read from file "+filename+" \n")

        return


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
        # try:
        print filename.split('-')[-1], 'in read_poscar'
        self.version = int(filename.split('-')[-1] )
        # except:
        #     pass
        print self.version, 'version'



        self.init = Structure()
        st = self.init
        with open(filename,'r') as f:
            name = f.readline().strip()
            # print self.name, "self.name"

            self.des = name
            # st.name = self.name

            mul = float( f.readline() )
            print 'mul', mul


            st.rprimd = []
            for i in 0, 1, 2:
                vec = f.readline().split()
                st.rprimd.append( np.asarray([float(vec[0])*mul, float(vec[1])*mul, float(vec[2])*mul]) )

            st.nznucl = []
            for z in f.readline().split():
                st.nznucl.append( int(z)  )

            type_of_coordinates = f.readline()

            st.xred = []


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
                st.xcart = xred2xcart(st.xred, st.rprimd)

            elif 'None' in type_of_coordinates:
                pass

            else:
                print_and_log("Error! The type of coordinates should be 'car' or 'dir' ")
                raise NameError


            print '!Name should contain names of elements; Name is ', name.split('!'), 'you could use ! to add comment after name'

            st.znucl = []
            for elname in name.split('!')[0].strip().split():
                st.znucl.append( element_name_inv(elname) )

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

        print self.init.znucl, 'self.init'


        print_and_log( "POSCAR was read\n")
        return



    def write_geometry(self, geotype = "init", description = "", override = False):
        """Writes geometrical data to self.path["input_geo"]"""
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







class CalculationAbinit(Calculation):
    """docstring for CalculationAbinit"""
    pass









class CalculationVasp(Calculation):
    """Methods for calculations made using VASP DFT code"""

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
        print 
        for i in 0, 1, 2:
            N_from_kspacing.append( math.ceil( (np.linalg.norm(self.init.recip[i]) / to_ang_local) / self.set.vasp_params['KSPACING']) )
        #print "Vector length is:", (np.linalg.norm(self.rprimd[0]), "Bohr"

        if self.set.kpoints_file  == False:#self.set.vasp_params['KSPACING']:
            N = N_from_kspacing

        if self.set.ngkpt:
            N = self.set.ngkpt
        else:
            N = N_from_kspacing


        print_and_log("\nKpoint   mesh is: "+str(N) )
        print_and_log("The actual k-spacings are "+str(self.calc_kspacings(N) ) )
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
                global to_ang; to_ang = 1.0; print "Conversion multiplier to_ang is",to_ang 
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

        
        print 
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
            potcar_files += self.set.potdir[ int(z) ]+"/POTCAR "

        runBash("cat "+potcar_files+" >"+path_to_potcar)

        print_and_log( "POTCAR files: "+potcar_files+"\n")        
        return


    def calculate_nbands(self):
        """Should be run after add_potcar()"""
        #1 add additional information to set
        path_to_potcar = self.dir+'/POTCAR'
        self.init.zval = []
        for line in open(path_to_potcar,'r'):
            if "ZVAL" in line:
                self.init.zval.append(float(line.split()[5]))
        try: self.set.add_nbands
        except AttributeError: self.set.add_nbands = None

        if self.set.add_nbands != None:
            tve =0
            for i in range(self.init.ntypat):
                tve += self.init.zval[i] * self.init.nznucl[i]
            self.nbands = int ( round ( math.ceil(tve / 2.) * self.set.add_nbands ) )
            self.set.vasp_params['NBANDS'] = self.nbands
        return


    def make_incar_and_copy_all(self, update):
        """Makes Incar file for current calculation and copy all
        TO DO: there is no need to send all POSCAR files; It is enothg to send only one. However for rsync its not that crucial
        """
        #print "Begin make---------------------------------------------"
        
        
        #Generate incar


        with open(self.dir+"INCAR",'w') as f:
            f.write( 'SYSTEM = %s\n\n'%(self.des) )
            #f.write( 'Other parameters for this Run:\n' )
            for key in sorted(self.set.vasp_params):
                if self.set.vasp_params[key] == None:
                    continue
                #print type(self.set.vasp_params[key])
                if type(self.set.vasp_params[key]) == list:
                    f.write(key+" = "+str(self.set.vasp_params[key][0])+" "+str(self.set.vasp_params[key][1])+"\n")
                else:
                    f.write(key+" = "+str( self.set.vasp_params[key] ) +"\n")
            try:
                if self.set.vasp_params['ISPIN'] == 2 and self.init.magmom[0]: #write magnetic moments of atoms
                    f.write("MAGMOM = ")
                    for m in self.init.magmom:    
                        f.write(str(m)+' ')
                    f.write("\n")
            except KeyError:
                pass

            f.write("\n")
            #f.write( 'Electronic Relaxation 1:\n' )

            #f.write("\n")
            #f.write( 'Ionic Relaxation:\n' )
       
            #f.write("\n")
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
            print_and_log("Files to copy: "+string_of_paths)
            log.write( runBash("rsync -zave ssh "+string_of_paths+" "+cluster_address+":"+project_path_cluster+self.dir)+"\n" )
        #print "End make---------------------------------------------\n\n"



    def write_sge_script(self, input_geofile = "header", version = 1, option = None, prevcalcver = None, savefile = "all"):
        """Without arguments writes header, else adds sequence of calculatios
            option - 'inherit_xred' - control inheritance, or 'master' - run serial on master 
            prevcalcver - ver of previous calc; for first none
            savefile - all, allw - +wavecar

        """

        type = 'SGE'
        type = 'PBS'
        
        if type == 'SGE':
            parrallel_run_command = "mpirun -x PATH vasp"
        if type == 'PBS':
            parrallel_run_command = "mpiexec vasp"



        run_name = self.dir+self.id[0]+"."+self.id[1]+'.run'        
        if input_geofile == "header":
            with open(run_name,'w') as f:
                if type == 'SGE':
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
                    f.write("rm WAVECAR\n")

                if type == 'PBS':
                    f.write("#!/bin/bash   \n")
                    f.write("#PBS -N "+run_name+"\n")
                    f.write("#PBS -l walltime=99999999:00:00 \n")
                    f.write("#PBS -l nodes=1:ppn="+str(header.corenum)+"\n")
                    f.write("#PBS -r n\n")
                    f.write("#PBS -j eo\n")
                    f.write("#PBS -m bea\n")
                    f.write("#PBS -M dimonaks@gmail.com\n")
                    f.write("cd $PBS_O_WORKDIR\n")
                    f.write("PATH=/share/apps/vasp/bin/vasp:/usr/bin/mpiexec:$PATH \n")

                    # f.write("cd "+self.dir+"\n")

                    # f.write("module load sge\n")
                    # f.write("module load vasp/parallel/5.2.12\n\n")
                    f.write("rm WAVECAR\n")                




            # f.close()




        else:
            v = str(version)
            precont = str(prevcalcver)+'.CONTCAR'
            #POSC = str(version)+'.POSCAR'
            with open(run_name,'a') as f: #append information about run
                
                if option == 'inherit_xred' and prevcalcver:
                    f.write('grep -A '+str(self.init.natom)+ ' "Direct" '+precont+' >> '+input_geofile+ ' \n')

                f.write("cp "+input_geofile+" POSCAR\n")
                if option == 'master':
                    f.write("vasp >"+self.name+".log\n")
                else:
                    f.write(parrallel_run_command +" >"+self.name+".log\n")
                
                f.write("sleep 20\n")
                # f.write("mv " +"POSCAR "  + "CONTCAR\n") #test
                
                f.write("mv OUTCAR "  + v+".OUTCAR\n")
                f.write("mv CONTCAR " + v+".CONTCAR\n")
                
                if "all" in savefile:
                    f.write("mv CHG "     + v+".CHG\n")
                    f.write("mv CHGCAR "  + v+".CHGCAR\n")
                    f.write("mv DOSCAR "  + v+".DOSCAR\n")
                    # f.write("mv PROCAR "  + v+".PROCAR\n")
                    f.write("mv AECCAR0 " + v+".AECCAR0\n")
                    f.write("mv AECCAR2 " + v+".AECCAR2\n")
                    if 'w' in savefile:
                        f.write("mv WAVECAR " + v+".WAVECAR\n")
                    if 'x' in savefile:
                        f.write("mv vasprun.xml " + v+".vasprun.xml\n")
                else:

                    # f.write("rm CHG\n")
                    # f.write("rm CHGCAR\n")
                    pass
                f.write("rm WAVECAR\n")
                f.write("\n")

            f.close()
            runBash('chmod +x '+run_name)

            return
    def make_run(self):
        #Generate run file


        type = 'SGE'
        type = 'PBS'

        run_name = self.dir+self.id[0]+"."+self.id[1]+'.run'

        with open('run','a') as f:
            if type == 'SGE':
                #'qsub -pe 'mpi*' NCORES -l CLUSTER_TAG script.parallel.sh' for mpi-jobs which should run on CLUSTER_TAG (cmmd or cmdft)
                #IMPORTANT: NCORES must be a multiple of 8(24) on cmdft(cmmd). 
                f.write("qsub -pe 'mpi*' "+str(header.corenum)+" "+header.queue+" "+run_name+"\n") #str(self.set.np) #-l cmmd
                f.write('sleep 5\n')
                # runBash('chmod +x run')
            if type == 'PBS':
                f.write("cd "+self.dir+"\n")
                f.write("qsub "+self.id[0]+"."+self.id[1]+'.run'+"\n") #str(self.set.np) #-l cmmd
                f.write("cd -\n")
                f.write('sleep 5\n')                        

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

    def read_results(self, load = 0, out_type = '', voronoi = False, show = [] ):

        #Start to read OUTCAR
        # print self.version, 'version'
        path_to_outcar  = self.path["output"]
        path_to_contcar = self.dir+str(self.version)+".CONTCAR"
        path_to_xml = self.dir+str(self.version)+".vasprun.xml"

        # print self.version, type(self.version)
        # print path_to_outcar, path_to_contcar

        # print os.path.exists(path_to_contcar)
        contcar_exist   = False
 
        self.natom = self.init.natom

        """Copy from server """
        if load == 1:
            # print "rsync -ave ssh "+cluster_address+":"+project_path_cluster+path_to_outcar+" "+self.dir
            log.write( runBash("rsync -zave ssh "+cluster_address+":"+project_path_cluster+path_to_outcar+" "+self.dir)+'\n' ) #OUTCAR
            log.write( runBash("rsync -zave ssh "+cluster_address+":"+project_path_cluster+path_to_contcar+" "+self.dir)+'\n' ) #CONTCAR
            log.write( runBash("rsync -zave ssh "+cluster_address+":"+project_path_cluster+path_to_xml+" "+self.dir)+'\n' ) #CONTCAR
            
            if not os.path.exists(path_to_outcar): 
                print_and_log("\nNo OUTCAR file even on server; Continue... \n")
                """Calculate voronoi volume"""
                # print hasattr(self, 'vorovol')
                if voronoi:# and not hasattr(self, 'vorovol'):#out_type == 'e_seg':
                    calculate_voronoi(self)
                return
                # raise RuntimeError

        if not os.path.exists(path_to_outcar) and self.id[2] == 1:
            path_to_outcar = self.dir+"OUTCAR" # For compability only for version 1
            if load == 1:
                log.write( runBash("rsync -zave ssh "+cluster_address+":"+project_path_cluster+path_to_outcar+" "+self.dir)+'\n' )
            log.write("Warning! I have used OUTCAR instead of 1.OUTCAR\n")


        if not os.path.exists(path_to_outcar):
            print_and_log("\nNo OUTCAR file... \n")
            raise RuntimeError

        if os.path.exists(path_to_contcar):
            contcar_exist   = True





        """Start reading """
        s = runBash("grep 'General timing' "+path_to_outcar)
        if "g" in s:
            #print s
            self.state = "4. Calculation completed."
        else: self.state = "3. Partly completed!."
        # print self.state
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
                #print "sed '"+str(nw-11)+","+str(nw+8)+"d' "+path_to_outcar+">"+tmp+";mv "+tmp+" "+path_to_outcar
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
            # print self.set.vasp_params

            for line in outcarlines:

                #Check bands

                # if 'band No.' in line:
                #     kpoint = float(outcarlines[i_line-1].split()[1])
                #     lastocc = float(outcarlines[i_line+self.nbands].split()[2])
                #     lastbandno = outcarlines[i_line+self.nbands].split()[0]
                #     if lastocc > 0:
                #         print "Warning!!! at kpoint ", kpoint, " last band No. ",lastbandno, " is not empty ", lastocc





                #Check W(q)
                if 'operators is LMAX' in line:
                    lmax = int(line.split()[7])
                if "W(low)/X(q)" in line:
                    low = [ float(outcarlines[i].split()[4]) for i in range(i_line+1, i_line+lmax+1)]
                    high = [ float(outcarlines[i].split()[5]) for i in range(i_line+1, i_line+lmax+1)]
                    if any(v > 1e-3 for v in low+high):
                        print_and_log("W(q)/X(q) are too high, check output!\n")

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




                if "TOTAL-FORCE" in line:
                    # Calculate forces here...
                    forces = []
                    magnitudes = []
                    for j in range(0,self.init.natom):
                        parts = outcarlines[i_line+j+2].split()
                        x = float(parts[3])
                        y = float(parts[4])
                        z = float(parts[5])
                        forces.append([x,y,z])
                        magnitudes.append(math.sqrt(x*x + y*y + z*z))
                    average.append( red_prec( sum(magnitudes)/self.init.natom * 1000 ) )
                    imax = np.asarray(magnitudes).argmax()
                    maxforce.append( [imax, round(magnitudes[imax] * 1000)]  )

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




  # alpha Z        PSCENC =       490.96968543
  # Ewald energy   TEWEN  =     -9607.95993345
  # -1/2 Hartree   DENC   =     -3050.68217759
  # -exchange  EXHF       =         0.00000000
  # -V(xc)+E(xc)   XCENC  =       516.30683760
  # PAW double counting   =     10207.08402417   -10005.81713472
  # entropy T*S    EENTRO =         0.00000000
  # eigenvalues    EBANDS =     -1760.14896244
  # atomic energy  EATOM  =     13136.98513889
                  # Kohn-Sham hamiltonian: http://en.wikipedia.org/wiki/Kohn%E2%80%93Sham_equations
                  #kinetic energy
                  #+ the external potential + the exchange-correlation energy +
                  #+ Hartree (or Coulomb) energy
                if  "alpha Z        PSCENC" in line:
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
                    self.e_without_entr = float(line.split()[3]) #energy(sigma->0)
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
                    niter = int(line.split(')')[0].split('(')[-1].strip()) #number of scf iteration
                    mdstep_old = self.mdstep





                i_line += 1
            #Check total drift
            max_magnitude = max(magnitudes)
            max_tdrift    = max(tdrift)
            if max_magnitude < self.set.toldff/10: max_magnitude = self.set.toldff
            #if any(d > 0.001 and d > max_magnitude for d in tdrift):
            if max_tdrift > 0.001 and max_tdrift > max_magnitude:
                
                #print_and_log( ("Total drift is too high! At the end one component is %0.f %% of the maximum force, check output!\n") %(maxdrift)  )
                pass
            #else: maxdrift = 

            """Try to read xred from CONCAR and calculate xcart"""
            #print contcar_exist
            if contcar_exist:
                with open(path_to_contcar, 'r') as contcar:
                    
                    for line in contcar:
                        
                        if "Direct" in line:
                            self.end.xred = []
                            for i in range(self.natom):
                                xr = np.asarray ( [float(x) for x in contcar.next().split()] )
                                self.end.xred.append( xr )

                #print self.end.xred
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
            if 'force' in show:
                print "Maxforce by md steps (meV/A) = %s;"%(str(maxforce)  )
                print "Avforce by md steps = %s;"%(str(average)  )

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




            return outst

        else:
            print_and_log("Still no OUTCAR for mystery reason\n\n")
            # raise RuntimeError

    def get_chg_file(self, filetype = 'CHGCAR'):
        #cl - object of CalculationVasp class
        path_to_chg = self.dir+str(self.version)+"."+filetype
        if not os.path.exists(path_to_chg): 
            print 'Charge file is downloading'
            log.write( runBash("rsync -zave ssh "+cluster_address+":"+project_path_cluster+path_to_chg+" "+self.dir)+'\n' ) #CHG
            print path_to_chg, 'was downloaded'
            
        return path_to_chg

    def get_file(self, filename):
        #cl - object of CalculationVasp class
        # filename - standart Vasp file name
        path_to_file = self.dir+str(self.version)+"."+filename
        if not os.path.exists(path_to_file): 
            log.write( runBash("rsync -zave ssh "+cluster_address+":"+project_path_cluster+path_to_file+" "+self.dir)+'\n' ) #CHG
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
        
        if runBash("ssh "+cluster_address+" '[ -e "+   CHGCAR_sum       +""" ] || echo "NO"     ;' """): #true if file not exists
            print  CHGCAR_sum, "not exist. try to calculate it "
            log.write( runBash("ssh "+cluster_address+" '"+command1+"'")+'\n' ) 
        
        if runBash("ssh "+cluster_address+" '[ -e "+   baderlog       +""" ] || echo "NO"     ;' """): #true if file not exists
            print  baderlog, "not exist. try to calculate Bader "
            log.write( runBash("ssh "+cluster_address+" '"+command2+"'")+'\n' ) 
        ACF = runBash("ssh "+cluster_address+" 'cat "+path+v+".ACF.dat"  +"'" )
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
        print numbers
        imp_partial_chg = imp_valence_chg - float(ACF[numbers[0]].split()[4])

        mat_partial_chg = [mat_valence_chg - float(ACF[i].split()[4]) for i in numbers[1:] ]
        print "Partial charge of impurity ", imp_partial_chg
        print "Partial charges of neibouring Ti atoms", " ".join("{:.2f}".format(m) for m in mat_partial_chg)
        print "Partial charge of matrix", sum(mat_partial_chg)
        
        print "Sum of mat and imp charges:", sum(mat_partial_chg)+imp_partial_chg

        # print command1
        # print command2
        return


        # if not os.path.exists(path_to_chg): 
        #     log.write( runBash("rsync -ave ssh "+cluster_address+":"+project_path_cluster+path_to_chg+" "+self.dir)+'\n' ) #CHG
        #     print path_to_chg, 'was downloaded'
            
        return path_to_chg
"""Functions related to CalculationVasp class"""
def get_chg_file(cl):
    #cl - object of CalculationVasp class
    pass




class InputSet():
    """docstring for InputSet"""
    def __init__(self,ise):
        #super(InputSet, self).__init__()
        self.ise = ise
        self.potdir = {}
        self.units = "abinit"
        self.vasp_params = {}
        self.mul_enaug = 1
        self.history = "Here is my uneasy history( :\n"     
        # self.kpoints_file = False
        # self.use_ngkpt = False


    def update(self):
        c1 = 1; c2 = 1
        if self.units == "abinit":
            c1 = to_eV
            c2 = Ha_Bohr_to_eV_A
        #Update Vasp parameters
        if self.units == "vasp":
            c1 = 1
            c2 = 1
        if self.ecut == None:
            self.vasp_params['ENCUT'] = None
            self.vasp_params['ENAUG'] = None
        else:           
            self.vasp_params['ENCUT'] = self.ecut * c1* self.dilatmx * self.dilatmx
            self.vasp_params['ENAUG'] = self.mul_enaug * self.vasp_params['ENCUT']
        self.vasp_params['SIGMA'] = self.tsmear * c1
        self.vasp_params['EDIFF'] = self.toldfe * c1
        self.vasp_params['NELM'] = self.nstep
        self.vasp_params['NSW'] = self.ntime
        self.vasp_params['EDIFFG'] = -self.tolmxf * c2



    def set_LREAL(self,arg):
        if type(arg) is not str:
            sys.exit("\nset_LREAL error\n")
        old = self.vasp_params['LREAL']    
        self.vasp_params['LREAL'] = arg
        if old == arg:
            print "Warning! You did not change one of your parameters in new set"
            return
        self.history += "LREAL was changed from "+str(old)+" to "+str(arg) + "\n"
        self.update()

    def set_ngkpt(self,arg):
        if type(arg) is not tuple:
            sys.exit("\nset_ngkpt type error\n")
        old = copy.copy(self.ngkpt)     
        self.ngkpt = copy.copy(arg)
        self.kpoints_file = True
        self.vasp_params['KSPACING'] = None
        if old == arg:
            print "Warning! You did not change one of your parameters in new set"
            return
        self.history += "ngkpt was changed from "+str(old)+" to "+str(arg) + " and KPOINTS file was swithed on\n"
        self.update()

    def set_tsmear(self,arg):
        if type(arg) is not float:
            sys.exit("\nset_tsmear error\n")
        old = self.tsmear   
        self.tsmear = arg
        if old == arg:
            print "Warning! You did not change one of your parameters in new set"
            return
        self.history += "tsmear was changed from "+str(old)+" to "+str(arg) + "\n"
        self.update()

    def set_KGAMMA(self,arg):
        if type(arg) is not str:
            sys.exit("\nset_KGAMMA error\n")
        old = self.vasp_params['KGAMMA']    
        self.vasp_params['KGAMMA'] = arg
        if old == arg:
            print "Warning! You did not change one of your parameters in new set"
            return
        self.history += "KGAMMA was changed from "+str(old)+" to "+str(arg) + "\n"
        self.update()

    def add_conv_kpoint(self,arg):
        if type(arg) is not str:
            sys.exit("\nadd_conv_kpoint error\n")
        if arg in self.conv_kpoint:
            print "Warning! You already have this name in list"
            return    
        self.conv_kpoint.append(arg)
        self.history += "Name "+arg+" was added to self.conv_kpoint\n"
        self.update()

    def add_conv_tsmear(self,arg):
        if type(arg) is not str:
            sys.exit("\nadd_conv_tsmear type error\n")
        try:
            self.conv_tsmear[0]
        except AttributeError:
            log.write( "Error! Set "+self.ise+" does not have conv_tsmear, I create new\n")
            self.conv_tsmear = []
        if arg in self.conv_tsmear:
            print "Warning! You already have this name in list"
            return    
        self.conv_tsmear.append(arg)
        self.history += "Name "+arg+" was added to self.conv_tsmear\n"
        self.update()

    def add_conv(self,arg,type_of_conv):
        if type(arg) is not str:
            raise TypeError
        if type_of_conv not in ["kpoint_conv","tsmear_conv","ecut_conv","nband_conv","npar_conv"]:
            raise TypeError
        try:
            self.conv[type_of_conv][0]
        except AttributeError:
            log.write( "Warning! Set "+self.ise+" does not have conv, I create new\n")
            self.conv = {}
        except KeyError:
            log.write( "Warning! Set "+self.ise+" does not have list for this key in conv, I add new\n")
            self.conv[type_of_conv] = []
        except IndexError:
            pass
        if arg in self.conv[type_of_conv]:
            print_and_log( "Warning! You already have name %s in list of conv %s. Nothing done.\n" % \
            (str(arg), str(self.conv[type_of_conv])    ) )
            return    
        self.conv[type_of_conv].append(arg)
        self.history += "Name "+arg+" was added to self.conv["+type_of_conv+"]\n"
        log.write( "Name "+arg+" was added to self.conv["+type_of_conv+"] of set "+self.ise+" \n")
        self.update()


    def set_potential(self,znucl, arg):
        
        if type(arg) is not str:
            sys.exit("\nset_potential error\n")
        
        if znucl in self.potdir:
            if arg == self.potdir[znucl]:
                print_and_log( "Warning! You already have the same potential for "+str(znucl)+" element\n" )
                return    
        # print type(self.potdir)
        self.potdir[znucl] = arg
        self.history += "Potential for "+str(znucl)+" was changed to "+arg+"\n"
        self.update()
        return

    def set_compare_with(self,arg):
        if type(arg) is not str:
            raise TypeError ("\nset_compare_with error\n")
        self.compare_with += arg+" "

    def set_PREC(self,arg):
        if arg not in ["Normal", "Accurate"]:
            raise TypeError
        old = self.vasp_params['PREC']    
        self.vasp_params['PREC'] = arg
        if old == arg:
            print_and_log("Warning! You did not change one of your parameters in "+self.ise+" set\n")
            return
        self.history += "PREC was changed from "+str(old)+" to "+str(arg) + "\n"

    def set_ecut(self,arg):
        old = self.ecut    
        if arg == 'default':
            self.ecut = None
        elif type(arg) not in [float, int]:
            raise TypeError
        else:
            self.ecut = arg
        if old == arg:
            print_and_log("Warning! You did not change one of your parameters in "+self.ise+" set\n")
            return
        self.history += "ecut was changed from "+str(old)+" to "+str(arg) + "\n"
        log.write("ecut was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        self.update()

    def set_dilatmx(self,arg):
        if type(arg) not in [float, int]:
            raise TypeError
        old = self.dilatmx    
        self.dilatmx = arg
        if old == arg:
            print_and_log("Warning! You did not change one of your parameters in "+self.ise+" set\n")
            return
        self.history += "dilatmx was changed from "+str(old)+" to "+str(arg) + "\n"
        log.write("dilatmx was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")

    def set_add_nbands(self,arg):
        name = "add_nbands"  
        if type(arg) not in [float, ]:
            raise TypeError
        try: self.add_nbands
        except AttributeError: self.add_nbands = 1.
        old = self.add_nbands    
        self.add_nbands = arg
        if old == arg:
            print_and_log("Warning! You did not change  "+name+"  in "+self.ise+" set\n")
            return
        self.history += " "+name+"  was changed from "+str(old)+" to "+str(arg) + "\n"
        log.write(" "+name+"  was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        

    def set_nmdsteps(self,arg):
        name = "number of MD steps"
        if type(arg) not in [int, ]:
            raise TypeError
        try: self.ntime
        except AttributeError: self.ntime = 1.
        old = self.ntime    
        self.ntime = arg
        if old == arg:
            print_and_log("Warning! You did not change  "+name+"  in "+self.ise+" set\n")
            return
        self.history += " "+name+"  was changed from "+str(old)+" to "+str(arg) + "\n"
        log.write(" "+name+"  was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        self.update()


    def set_relaxation_type(self,type_of_relaxation):
        name = "Type of relaxation ISIF"
        if type(type_of_relaxation) not in [str, ]:
            raise TypeError
        old = self.vasp_params["ISIF"]
        if "ions" == type_of_relaxation:
            if int(self.ise[0]) != 9:
                print_and_log("Warning! The name of set is uncostintent with relaxation type\n")
                raise TypeError     
            self.vasp_params["ISIF"] = 2
            # self.set_nmdsteps(200)
        elif type_of_relaxation == "full":
            if int(self.ise[0]) != 2:
                print_and_log("Warning! The name of set is uncostintent with relaxation type\n")
                raise TypeError              
            self.vasp_params["ISIF"] = 3
        else:
            print_and_log("Error! Uncorrect type of relaxation\n")
            raise TypeError
        arg = self.vasp_params["ISIF"]
        if old == arg:
            print_and_log("Warning! You did not change  "+name+"  in "+self.ise+" set\n")
            return
        self.history += " "+name+"  was changed from "+str(old)+" to "+str(arg) + "\n"
        log.write(name+"  was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        self.update()                
        #print self.history

    def set_ALGO(self,arg):
        name = "speed of calculation ALGO"
        if arg not in ["Normal","Fast" ]:
            raise TypeError
        old = self.vasp_params["ALGO"]
        self.vasp_params["ALGO"] = arg
        if old == arg:
            print_and_log("Warning! You did not change  "+name+"  in "+self.ise+" set\n")
            return
        self.history += " "+name+"  was changed from "+str(old)+" to "+str(arg) + "\n"
        log.write(name+"  was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        self.update() 


    def set_IBRION(self,arg):
        name = "algorithm of relaxation  IBRION"
        if arg not in ["CG","damped","RMM", "MD" ]:
            raise TypeError
        old = self.vasp_params["IBRION"]
        if "CG" == arg:
            self.vasp_params["IBRION"] = 2
        elif "damped" == arg:
            self.vasp_params["IBRION"] = 3
        elif "RMM" == arg:
            self.vasp_params["IBRION"] = 1
        elif "MD" == arg:
            self.vasp_params["IBRION"] = 0
        else:
            print_and_log("Error! Uncorrect algorithm of relaxation\n")
            raise TypeError
        arg = self.vasp_params["IBRION"]
        if old == arg:
            print_and_log("Warning! You did not change  "+name+"  in "+self.ise+" set\n")
            return
        self.history += " "+name+"  was changed from "+str(old)+" to "+str(arg) + "\n"
        log.write(name+"  was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        self.update()


    def set_LPLANE(self,arg):
        name = "control of parallezation LPLANE"
        if arg not in [".TRUE.",".FALSE." ]:
            raise TypeError
        old = self.vasp_params["LPLANE"]
        self.vasp_params["LPLANE"] = arg
        if old == arg:
            print_and_log("Warning! You did not change  "+name+"  in "+self.ise+" set\n")
            return
        self.history += " "+name+"  was changed from "+str(old)+" to "+str(arg) + "\n"
        log.write(name+"  was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        self.update() 

    def set_NPAR(self,arg):
        name = "control of parallezation NPAR"
        if type(arg) not in [int, ]:
            raise TypeError
        old = self.vasp_params["NPAR"]
        self.vasp_params["NPAR"] = arg
        if old == arg:
            print_and_log("Warning! You did not change  "+name+"  in "+self.ise+" set\n")
            return
        self.history += " "+name+"  was changed from "+str(old)+" to "+str(arg) + "\n"
        log.write(name+"  was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        self.update()


    def set_NELMIN(self,arg):
        name = "minimum number of scf steps NELMIN"
        if type(arg) not in [int, ]:
            raise TypeError
        old = self.vasp_params["NELMIN"]
        self.vasp_params["NELMIN"] = arg
        if old == arg:
            print_and_log("Warning! You did not change  "+name+"  in "+self.ise+" set\n")
            return
        self.history += " "+name+"  was changed from "+str(old)+" to "+str(arg) + "\n"
        log.write(name+"  was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        self.update()

    def set_tolmxf(self,arg):
        name = "relaxation until the forces will be less than tolmxf"
        if type(arg) not in [float, ]:
            raise TypeError
        old = self.tolmxf
        self.tolmxf = arg
        if old == arg:
            print_and_log("Warning! You did not change  "+name+"  in "+self.ise+" set\n")
            return
        self.history += " "+name+"  was changed from "+str(old)+" to "+str(arg) + "\n"
        log.write(name+"  was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        self.update()

    def set_toldfe(self,arg):
        name = "energy tolerance toldfe"
        if type(arg) not in [float, ]:
            raise TypeError
        old = self.toldfe
        self.toldfe = arg
        if old == arg:
            print_and_log("Warning! You did not change  "+name+"  in "+self.ise+" set\n")
            return
        self.history += " "+name+"  was changed from "+str(old)+" to "+str(arg) + "\n"
        log.write(name+"  was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        self.update()

    def set_vaspp(self,token,arg,des="see manual"):

        if token in ("ISMEAR",):
            if type(arg) not in [int, None, ]:
                raise TypeError
        if token in ("KSPACING",):
            if type(arg) not in [float, None, ]:
                raise TypeError


        old = self.vasp_params[token]
        self.vasp_params[token] = arg
        if old == arg:
            print_and_log("Warning! You did not change  "+token+"  in "+self.ise+" set\n")
            return
        self.history += " "+token+"  was changed from "+str(old)+" to "+str(arg) + "\n"
        log.write(token+"  was changed from "+str(old)+" to "+str(arg) +" - "+ des+" in set "+self.ise+" \n")
        self.update()



def clean_run():
    type = 'PBS'
    with open('run','w') as f:
    
        if type == 'SGE':
            f.write("#!/bin/tcsh\n")
            f.write("module load sge\n")
            f.write("module load vasp/parallel/5.2.12\n")
        if type == 'PBS':
            f.write("#!/bin/bash\n")

    # f.close()
    return


def complete_run(close_run = True):
    
    if close_run:

        f = open('run','a')
        f.write("qstat\n")
        f.write("sleep 2\n")
        f.write("mv run last_run\n")
        f.close()
        runBash('chmod +x run')

        log.write( runBash("rsync -zave ssh run "+cluster_address+":"+project_path_cluster) +"\n" )
        print 'run sent'
        clean_run()
    
    return


def plot_conv(list_of_calculations, calc, type_of_plot, conv_ext = [], labelnames = None):
    """
    Allows to fit and plot different properties;
    Input:
    'type_of_plot' - ("fit_gb_volume"-fits gb energies and volume and plot dependencies without relaxation and after it,
     )


    """

            
    def fit_and_plot(x1, y1, x2, y2, power, name = "", xlabel = "", ylabel = "", image_name = "test", lines = None):
        """Should be used in two below sections!
        Creates one plot with two dependecies and fit them;
        return minimum fitted value of x2 and corresponding valume of y2; 
        if name == "" image will not be plotted
        power - the power of polynom

        lines - add lines at x = 0 and y = 0

        """
        coeffs1 = np.polyfit(x1, y1, power)        
        coeffs2 = np.polyfit(x2, y2, power)
        
        fit_func1 = np.poly1d(coeffs1)
        fit_func2 = np.poly1d(coeffs2)
        
        #x_min  = fit_func2.deriv().r[power-2] #derivative of function and the second cooffecient is minimum value of x.
        #y_min  = fit_func2(x_min)
        
        if name:

            x_range = np.linspace(min(x2), max(x2))
            fit_y1 = fit_func1(x_range); 
            fit_y2 = fit_func2(x_range); 
            
            plt.figure(figsize=(8,6.1))
            # plt.title(name)
            plt.ylabel(ylabel)
            plt.xlabel(xlabel)
            plt.xlim(min(x2)-0.1*abs(min(x2) ), max(x2)+0.1*abs(min(x2)))

            plt.plot(x1, y1, 'ro', label = 'initial')
            plt.plot(x2, y2, 'bo', label = 'relaxed'   )
            plt.plot(x_range, fit_y1, 'r-',) #label = 'init_fit')
            plt.plot(x_range, fit_y2, 'b-',) #label = 'r_fit'   )
            plt.legend(loc =9)
            
            if lines == 'xy':
                plt.axvline(color='k')
                plt.axhline(color='k')



            plt.tight_layout()
            #plt.savefig('images/'+image_name)
            print 'Saving file ...',path_to_images+str(image_name)+'.png'
            plt.savefig(path_to_images+str(image_name)+'.png',format='png', dpi = 300)
        return fit_func2  




    conv = list_of_calculations
    name = []; n = conv[0]
    name.append( n[0] )
    image_name = n[0]+'_'+n[1]+'_'+str(n[2])

    energies = []; init_energies = []
    volumes = []
    gb_volumes = []
    pressures = []
    pressures_init = []
    sigma_xx = []
    sigma_yy = []
    sigma_zz = []
    e_gbs = [] 
    e_gbs_init = []

    if type_of_plot == "e_imp":
        e_imps = []
        v_imps = []
        lengths = []
        for id in conv:        
            e_imps.append(calc[id].e_imp*1000)
            v_imps.append(calc[id].v_imp)
            l = calc[id].vlength
            lengths.append( "%s\n%.1f\n%.1f\n%.1f"%(id[0],l[0], l[1], l[2]) )
        #l = lengths[0]
        #print str(l[0])+'\n'+str(l[1])+'\n'+str(l[2])


        xlabel = "Sizes, $\AA$"
        ylabel = "Impurity energy, meV"
        ylabel2 = "Impurity volume, $\AA^3$"
        plt.figure()
        plt.title(str(name)+' other cells')
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        x = range( len(e_imps) )
        plt.xticks(x, lengths)
        plt.plot(x, e_imps, 'ro-', label = 'energy')
        plt.legend()
        plt.twinx()
        plt.ylabel(ylabel2)
        plt.plot(x, v_imps, 'bo-', label = 'volume')

        plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None,
                wspace=None, hspace=None)
        #plt.ticker.formatter.set_scientific(True)
        plt.legend(loc =9)
        plt.savefig('images/e_imp_'+str(image_name)+'.png',format='png')#+str(image_name))#+'e_imp')


    if type_of_plot == "e_2imp":

        def dist_between_imp(cl):
            """Only for two impurities"""

            return np.linalg.norm(cl.end.xcart[-1] - cl.end.xcart[-2]) #assuming that impurities are at the end of xcart list.

        e_imps = [] # binding energy
        dist = [] #dist between impurities
        
        e_imps_ex = []
        dist_ex = []
        name_ex = []

        for id in conv:        
            cl = calc[id]
            e_imps.append(cl.e_imp*1000)
            #dist.append( "%s\n%.1f"%(id[0],dist_between_imp(cl) ) )
            dist.append( dist_between_imp(cl)  )







        xlabel = "Distance between atoms, $\AA$"
        ylabel = "Interaction energy, meV"
        plt.figure()
        
        # plt.title(str(name)+' v1-15')

        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        #x = range( len(e_imps) )
        #plt.xticks(x, dist)
        # plt.yscale('log')
        # plt.yscale('semilog')
        if labelnames:
            label = labelnames
        else:
            label = []
            label[0] = str(name)
            label[1] = name_ex[0]
            label[2] = name_ex[1]


        plt.plot(dist, e_imps, 'ro-', label = label[0], linewidth = 2 )
        

        if conv_ext: #manually add 
            for conv in conv_ext:
                e_imps_ex.append([])
                dist_ex.append([])
                for id in conv:        
                    cl = calc[id]
                    e_imps_ex[-1].append(cl.e_imp*1000)
                    #dist.append( "%s\n%.1f"%(id[0],dist_between_imp(cl) ) )
                    dist_ex[-1].append( dist_between_imp(cl)  )
                name_ex.append(id[0])
            plt.plot(dist_ex[0], e_imps_ex[0], 'go-', label = label[1], linewidth = 2)
            plt.plot(dist_ex[1], e_imps_ex[1], 'bo-', label = label[2], linewidth = 2)






        plt.axhline(color = 'k') #horizontal line

        plt.tight_layout()

        # plt.subplots_adjust(left=None, bottom=0.2, right=None, top=None,
        #         wspace=None, hspace=None)
        # #plt.ticker.formatter.set_scientific(True)
        plt.legend(loc =9)
        plt.savefig(path_to_images+'e_2imp_'+str(image_name)+'.png',format='png', dpi = 300)#+str(image_name))#+'e_imp')


    if type_of_plot == "fit_gb_volume_pressure":

        for id in conv:
            #energies.append(calc[id].energy_sigma0)
            #init_energies.append( calc[id].list_e_sigma0[0] ) 
            gb_volumes.append(calc[id].v_gb)
            #volumes.append(calc[id].end.vol)
            pressures.append(calc[id].extpress/1000. )
            pressures_init.append(calc[id].extpress_init/1000. )
            sigma_xx.append( calc[id].stress[0]  )
            sigma_yy.append( calc[id].stress[1]  )
            sigma_zz.append( calc[id].stress[2]  )
            #pressures_init = pressures
            e_gbs.append(calc[id].e_gb)
            e_gbs_init.append(calc[id].e_gb_init )           
            # print calc[id].bulk_extpress

        power = 3

        fit_ve = fit_and_plot(gb_volumes, e_gbs_init,  gb_volumes, e_gbs, power, 
            name, "Grain boundary expansion (m$\AA$)", "Grain boundary energy (mJ/m$^2$)", 
            image_name+"_fit_ve")


        fit = fit_and_plot(pressures_init, e_gbs_init,  pressures, e_gbs, power, 
            name, "External pressure (GPa)", "Grain boundary  energy (mJ/m$^2$)", 
            image_name+"_pe")
        #print fit
        ext_p_min  = fit.deriv().r[power-2] #external pressure in the minimum; derivative of function and the value of x in minimum

        fit_sxe = fit_and_plot(sigma_xx, e_gbs_init,  sigma_xx, e_gbs, power, 
            name, "Sigma xx (MPa)", "Grain boundary energy (mJ/m$^2$)", 
            image_name+"_sxe")
        sxe_min = fit_sxe.deriv().r[power-2] #sigma xx at the minimum of energy
        print "sigma xx at the minimum of energy is", sxe_min," MPa"


        fit1 = fit_and_plot(pressures_init, gb_volumes,  pressures, gb_volumes, 1,
            name, "External pressure (GPa)", "Grain boundary expansion (m$\AA$)", 
            image_name+"_pv", lines = 'xy')
        #print fit1
        pulay = - calc[id].bulk_extpress
        #print " At external pressure of %.0f MPa; Pulay correction is %.0f MPa." % (ext_p_min+pulay, pulay)       
        #print " Egb = %.1f mJ m-2; Vgb = %.0f mA;"%(fit(ext_p_min), fit1(ext_p_min)  )
        print ("%s.fit.pe_pv & %.0f & %.0f & %0.f & %0.f \\\\" %
            (n[0]+'.'+n[1], fit(ext_p_min), fit1(ext_p_min),ext_p_min, ext_p_min+pulay   ))


        #print "\n At zero pressure with Pullay correction:"
        #print " Egb = %.1f mJ m-2; Vgb = %.0f mA; " % (fit(-pulay), fit1(-pulay))
        outstring =  ("%s.fit.pe_pv & %.0f & %.0f & %0.f & %0.f\\\\" %(n[0]+'.'+n[1], fit(-pulay), fit1(-pulay),-pulay,0    ))
        # print outstring
        calc[conv[0]].egb = fit(-pulay)
        calc[conv[0]].vgb = fit1(-pulay)

        return outstring #fit(ext_p_min), fit1(ext_p_min) 


    if type_of_plot == "fit_gb_volume":
        """
        should be rewritten using fit_and_plot() function
        """

        for id in conv:
            #energies.append(calc[id].energy_sigma0)
            #init_energies.append( calc[id].list_e_sigma0[0] ) 
            gb_volumes.append(calc[id].v_gb)
            e_gbs.append(calc[id].e_gb)
            e_gbs_init.append(calc[id].e_gb_init )           


        power = 3
        fit_ve = fit_and_plot(gb_volumes, e_gbs_init,  gb_volumes, e_gbs, power, 
            name, "Excess volume ($m\AA$)", "Twin energy ($mJ/m^2$)", 
            image_name+"_fit_ve")

        vgb_min  = fit_ve.deriv().r[power-2]


        #print "Fit of excess volume against energy. Pressure is uknown:"
        #print "Test Egb_min = %.1f mJ m-2; v_min = %.0f mA;"%(fit_ve(vgb_min), vgb_min)
        print ("%s.fit.ve & %.0f & %.0f & - & - \\\\" %
            (n[0]+'.'+n[1], fit_ve(vgb_min), vgb_min,   ))

    if type_of_plot == "fit_gb_volume2":

        for id in conv:
            energies.append(calc[id].energy_sigma0)
            init_energies.append( calc[id].list_e_sigma0[0] ) 
            volumes.append(calc[id].end.vol)
            pressures.append(calc[id].extpress )
            pressures_init.append(calc[id].extpress_init )

        power = 3
        pulay = 500

        fit_ve = fit_and_plot(volumes, init_energies,  volumes, energies, power, 
            name, "Volume ($\AA^3$)", "Energy  sigma->0 ($eV$)", 
            image_name+"_fit_ve")
        
        Vmin  = fit_ve.deriv().r[power-2] # minimum volume at the minimum energy
        Emin  = fit_ve(Vmin)

        fit_pe = fit_and_plot(pressures_init, init_energies,  pressures, energies, power, 
            name, "External pressure ($MPa$)", "Energy  sigma->0 ($eV$)", 
            image_name+"_fit_pe")

        ext_p_min  = fit_pe.deriv().r[power-2] #external pressure in the minimum; derivative of function and the value of x in minimum
        

        fit_pv = fit_and_plot(pressures_init, volumes,  pressures, volumes, 1,
            name, "External pressure ($MPa$)", "Volume of cell ($\AA^3$)", 
            image_name+"_fit_pv")


             
        atP = (" Emin = %.3f meV;  Vmin = %.0f A^3; "%( fit_pe(ext_p_min), fit_pv(ext_p_min)  )  ) + \
              (" for the minimum of energy relative to external pressure. The value of pressure is %.0f MPa; Pulay correction is %.0f MPa." % (ext_p_min+pulay, pulay) )
        
        at_zeroP = (" Emin = %.3f meV;  Vmin = %.0f A^3; " % (fit_pe(-pulay), fit_pv(-pulay) )  ) + \
                   (" the value of energy and volume at zero pressure with Pullay correction" )
        
        #print " Emin = %.3f meV;  Vmin = %.0f A^3;  for the minimum of energy relative to volume at some external pressure"%(fit_ve(Vmin), Vmin )
        #print atP
        #print at_zeroP

        print "Compare V at -pulay and V for energy minimum", fit_pv(-pulay), Vmin

        return fit_pe(-pulay), fit_pv(-pulay), Emin, Vmin








    if type_of_plot == "kpoint_conv":
        energies = []
        kpoints = []
        times = []

        for id in list_of_calculations:
            if "4" not in calc[id].state:
                continue
            energies.append(calc[id].potenergy)
            kpoints.append(calc[id].kspacing[2])
            times.append(calc[id].time)

            name.append( id[1] )

        plt.figure()
        plt.title(name)
        plt.plot(kpoints, energies,'bo-')
        plt.ylabel("Total energy (eV)")
        plt.xlabel("KSPACING along 3rd recip. vector ($\AA ^{-1}$)")
        plt.twinx()
        plt.plot(kpoints,times,'ro-')
        plt.ylabel("Elapsed time (min)")
        plt.savefig('images/'+str(conv[0])+'kconv')


    if type_of_plot == "contour":
        alist = [] ;        clist = []
        nn = str(calc[conv[0]].id[0])+"."+str(calc[conv[0]].id[1])
        f = open("a_c_convergence/"+nn+"/"+nn+".out","w")
        f.write("END DATASET(S)\n")
        k = 1
        for id in conv: #Find lattice parameters and corresponding energies
            a = calc[id].a
            c = calc[id].c
            if a not in alist: alist.append(a); 
            if c not in clist: clist.append(c);
            f.write( "acell%i %f %f %f Bohr\n"%(k, calc[id].a/to_ang, calc[id].a/to_ang, calc[id].c/to_ang )   )
            #print "etotal%i %f\n"%(k, calc[id].energy_sigma0/to_eV ),
            k+=1;
        X,Y = np.meshgrid(alist, clist)
        Z = np.zeros(X.shape)
        Zinv = np.zeros(X.shape)

        
        k=1
        for i in range(len(alist)):
            for j in range(len(clist)):
                for id in conv:
                    if calc[id].a == alist[i] and calc[id].c == clist[j]:
                        Z[i][j] = calc[id].energy_sigma0
                        Zinv[j][i] = calc[id].energy_sigma0
                        f.write( "etotal%i %f\n"%(k, calc[id].energy_sigma0/to_eV )   )
                        k+=1
        f.write("+Overall time at end (sec) : cpu=     976300.2  wall=     976512.8")
        f.close

        #Make two plots for different a and c
        plt.figure()
        plt.title(name)
        for i in range(len(alist)):
            plt.plot(clist, Z[i],'o-',label='a='+str(alist[i]))
        plt.legend()
        plt.ylabel("Total energy (eV)")
        plt.xlabel("c parameter ($\AA$)")
        plt.savefig('images/'+str(conv[0])+'c')

        plt.figure()
        plt.title(name)
        for j in range(len(clist)):
            plt.plot(alist, Zinv[j],'o-',label='c='+str(clist[j]))
        plt.legend()
        plt.ylabel("Total energy (eV)")
        plt.xlabel("a parameter ($\AA$)")
        plt.savefig('images/'+str(conv[0])+'a')

        #Make contour
        plt.figure()
        cf = plt.contourf(X, Y, Z, 20,cmap=plt.cm.jet)
        cbar = plt.colorbar(cf)
        cbar.ax.set_ylabel('Energy (eV)')

        plt.xlabel('$a$ ($\AA$)')
        plt.ylabel('$c/a$')

        plt.legend()
        plt.savefig('images/ru-contourf.png')
        #plt.show()


        #Make equation of state
        eos = EquationOfState(clist,Z[2])
        v0, e0, B = eos.fit()
        #print "a = ", alist[2]
        print '''
v0 = {0} A^3
E0 = {1} eV
B  = {2} eV/A^3'''.format(v0, e0, B)
        eos.plot('images/a[2]-eos.png')

        eos = EquationOfState(alist,Zinv[2])
        v0, e0, B = eos.fit()
        #print "c = ", clist[2]
        print '''
v0 = {0} A^3
E0 = {1} eV
B  = {2} eV/A^3'''.format(v0, e0, B)
        eos.plot('images/c[2]-eos.png')


    if type_of_plot == "dimer":


        x1 = [] #list of distances

        cl =  calc[list_of_calculations[0]]
        if cl.end.natom > 2:
            raise RuntimeError


        for xcart in cl.end.list_xcart:
            # x = xcart[1]
            # d = (x[0]**2 + x[1]**2 + x[2]**2)**0.5
            d = np.linalg.norm(xcart[1]-xcart[0]) #assuming there are only two atoms


            x1.append(d)

        y1 = cl.list_e_without_entr
        power = 4
        name = 'dimer'
        xlabel = 'Bond length'
        ylabel = 'Full energy'
        coeffs1 = np.polyfit(x1, y1, power)        
      
        fit_func1 = np.poly1d(coeffs1)

        x_range = np.linspace(min(x1), max(x1))
        fit_y1 = fit_func1(x_range); 
        f = fit_func1.deriv()
        min_e = fit_func1(f.r[2]).real
        print "The minimum energy per atom and optimal length of dimer are {:.3f} eV and {:.3f} A".format( min_e/2., f.r[2].real)
        print "The atomization energy for dimer is {:.3f} eV ; The energy of atom in box is taken from the provided b_id".format(min_e - 2*cl.e_ref)

        plt.figure()
        plt.title(name)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.plot(x1, y1, 'ro', label = 'init')
        plt.plot(x_range, fit_y1, 'r-', label = 'init_fit')
        plt.show()
    return















def add_calculation(structure_name,inputset,version,first_version,last_version,input_folder,blockdir,calc,varset,update = "no"
    ,inherit_option = None, prevcalcver = None, coord = 'direct', savefile = "all", input_geo_format = 'abinit'):
    """Adds new or updates 
    up2 - allows to update only unfinished

    input_geo_format - abinit, vasp
    """
    id = (structure_name,inputset,version)
    try:
        calc[id]
        status = "exist"
        print_and_log(str(calc[id].name)+" has been already created and has state: "+str(calc[id].state)+"\n\n")
        if "4" in calc[id].state: 
            complete_state = calc[id].state
            status = "compl"
            if update == 'up2': 
                print 'Calculation', calc[id].name, 'is finished, continue'
                return
            pass
            if update != "up1": return #completed calculations updated only for "up1"

    except KeyError:
        #update = "up"
        status = "new"
        print_and_log( "There is no calculation with id "+ str(id)+". I create new with set "+str(inputset)+"\n" )        
    if "up" in update:
        if status in ["exist","compl"]: print_and_log("You asked to update existing calculation with id "+ str(id)+" Warning! I update creating new class \n" )         

        calc[id] = CalculationVasp( varset[id[1]] )
        calc[id].id = id 
        calc[id].name = str(id[0])+'.'+str(id[1])+'.'+str(id[2])

        calc[id].dir = blockdir+"/"+ str(id[0]) +'.'+ str(id[1])+'/'
        
        if update in ['up1', 'up2']:
            if not os.path.exists(calc[id].dir):
                log.write( runBash("mkdir -p "+calc[id].dir) )         #Create directory if does not exist
                log.write( runBash("ssh "+cluster_address+" ' mkdir -p "+calc[id].dir+" ' ") )
            if id[2] == first_version:
                calc[id].write_sge_script() #write header only once

        print_and_log("I am searching for geofiles in folder "+input_folder+"\n" )
        if input_geo_format == 'abinit': 
            searchinputtemplate = input_folder+'/*.geo*'
        elif input_geo_format == 'vasp': 
            searchinputtemplate = input_folder+'/POSCAR*'
        # print  input_geo_format
        geofilelist = glob.glob(searchinputtemplate) #Find input_geofile
        # print geofilelist

        #additional search in target folder if no files in root
        if not geofilelist:
            print_and_log("Attention! trying to find here "+input_folder+"/target\n" )
            geofilelist = glob.glob(input_folder+'/target/*.geo*') #Find input_geofile            


        for input_geofile in geofilelist:
            #print runBash("grep version "+str(input_geofile) )
            if input_geo_format == 'abinit':
                curv = int( runBash("grep version "+str(input_geofile) ).split()[1] )
            elif input_geo_format == 'vasp': 
                print input_geofile, 'is input_geofile'
                print input_geofile.split('-')[-1], ' is ver'

                curv = int(input_geofile.split('-')[-1] ) #!Applied only for phonopy POSCAR-n naming convention
                print curv, 'is curv'



            if curv == id[2]:
                
                if input_geo_format == 'abinit':
                    calc[id].read_geometry(input_geofile)
                
                elif input_geo_format == 'vasp':
                    calc[id].read_poscar(input_geofile)
                
                else:
                    raise RuntimeError

                

                break
        if calc[id].path["input_geo"] == None: 
            print_and_log("Could not find geofile in this list: "+ str(geofilelist)+  "\n")
            raise NameError #
        #print  calc[id].des
        calc[id].des += ' '+struct_des[id[0]].des + '; ' + varset[id[1]].des
        #print  calc[id].des

        calc[id].check_kpoints()    





        if update in ['up1', 'up2']:
            calc[id].write_structure(str(id[2])+".POSCAR", coord, inherit_option, prevcalcver)
            calc[id].write_sge_script(str(version)+".POSCAR", version, inherit_option, prevcalcver, savefile)
            #if calc[id].path["output"] == None:
            calc[id].path["output"] = calc[id].dir+str(version)+".OUTCAR" #set path to output

            if id[2] == first_version:
                calc[id].add_potcar()
            
            calc[id].calculate_nbands()

            if id[2] == last_version:        
                calc[id].make_incar_and_copy_all(update)
                calc[id].make_run()


        if status == "compl": calc[id].state = complete_state #Even if completed state was updated, the state does not change


        print_and_log("\nCalculation "+str(id)+" added or updated\n\n")

    return













        
def headers():
    j = (7,12,14,7,8,9,9,5,5,20,5,20,8,12,20,8,5,8,8)
    d="&"
    header_for_bands= "Set".ljust(j[0])+d+"Etot".center(j[1])+d+"a1,a2".center(j[2])+d+"c".center(j[3])\
                +d+"time, m".center(j[4])+d+"ittime, s".center(j[5])+d+"Nmd,Avr.".rjust(j[6])+d\
                +"Warn!"+d+"nband"+d+"Added, \%"+"\\\\"

    header_for_ecut= "Set".ljust(j[0])+d+"Etot".center(j[1])+d+"a1,a2".center(j[2])+d+"c".center(j[3])\
                +d+"time, m".center(j[4])+d+"ittime, s".center(j[5])+d+"Nmd,Avr.".rjust(j[6])+d\
                +"Warn!"+d+"Ecut,eV"+"\\\\"

    header_for_npar= "Set".ljust(j[0])+d+"Etot".center(j[1])+d+"a1,a2".center(j[2])+d+"c".center(j[3])\
                +d+"time, m".center(j[4])+d+"ittime, s".center(j[5])+d+"Nmd,Avr.".rjust(j[6])+d\
                +"Warn!"+d+"NPAR".center(j[16])+d+"LPLANE".center(j[17])+"\\\\"

    header_for_kpoints= "Set".ljust(j[0])+d+"Etot".center(j[1])+d+"a1,a2".center(j[2])+d+"c".center(j[3])\
                +d+"time, m".center(j[4])+d+"ittime, s".center(j[5])+d+"Nmd,Avr.".rjust(j[6])+d\
                +"Warn!"+d+"k-mesh".center(j[8])+d+"k-spacings".center(j[9])+d+"nkpt".center(j[10])+"\\\\"
    header_for_tsmear= "Set".ljust(j[0])+d+"Etot".center(j[1])+d+"a1,a2".center(j[2])+d+"c".center(j[3])\
                +d+"time, m".center(j[4])+d+"ittime, s".center(j[5])+d+"Nmd,Avr.".rjust(j[6])+d\
                +"Warn!"+d+"k-mesh".center(j[8])+d+"tsmear, meV".center(j[13])+d+"Smearing error, meV/atom".center(j[14])+"\\\\"

    header_for_stress= "Set".ljust(j[0])+d+"Etot".center(j[1])+d+"a1,a2".center(j[2])+d+"c".center(j[3])\
                +d+"time, m".center(j[4])+d+"ittime, s".center(j[5])+d+"Nmd,Avr.".rjust(j[6])+d\
                +"Warn!"+d+"Stress, intr u.*1000".center(j[11])+d+"Pressure, MPa".center(j[12])
    #print "\\hline"
    return header_for_kpoints

def gb_energy_volume(gb,bulk):
    if (gb.end.rprimd[1] != bulk.end.rprimd[1]).any() or (gb.end.rprimd[2] != bulk.end.rprimd[2]).any():
        print_and_log("Warning! You are trying to calculate gb_energy from cells with different lateral sizes:"+str(gb.end.rprimd)+" "+str(bulk.end.rprimd)+"\n")
    #print bulk.vol
    V_1at = bulk.vol / bulk.natom #* to_ang**3

    E_1at = bulk.energy_sigma0 / bulk.natom 
    A = np.linalg.norm( np.cross(gb.end.rprimd[1], gb.end.rprimd[2])  ) #surface area of gb
    #print A
    gb.v_gb =      ( gb.vol              - V_1at * gb.natom) / A / 2. * 1000
    gb.e_gb =      ( gb.energy_sigma0    - E_1at * gb.natom) / A / 2. * eV_A_to_J_m * 1000
    gb.e_gb_init = ( gb.list_e_sigma0[0] - E_1at * gb.natom) / A / 2. * eV_A_to_J_m * 1000
    gb.bulk_extpress = bulk.extpress     
    #print "Calc %s; e_gb_init = %.3f J/m^2; e_gb = %.3f J/m; v_gb = %.3f angstrom "%(gb.name, gb.e_gb_init, gb.e_gb, gb.v_gb )
    outst = "%15s&%7.0f&%7.0f"%(gb.name, gb.e_gb, gb.v_gb)
    return outst




def add_loop(it,ise,verlist, calc = None, conv = None, varset = None, up = 'test',typconv="",from_geoise = '', inherit_option = None, 
    coord = 'direct', savefile = "allx", show = [], comment = None, input_geo_format = 'abinit'):
    """To create folders and add calculations add_flag should have value 'add' 

    'from_geoise' - part of folder name with geometry input files. allows to use geomery files from different sets.
    please find how it was used

    'up' == "no_base": only relevant for typconv
     is the same as "up1", but the base set is ommited
     up2 - update only unfinished calculations

    savefile - controls files to be saved

    """
    if not calc:
        calc = header.calc
        conv = header.conv
        varset = header.varset

    header.close_run = True
    if type(verlist) == int: #not in [tuple, list]:
        #print "verlist"
        verlist = [verlist]; #transform to list
    #global history;
    #hstring = ( "add_loop('%s','%s',%s,calc,conv,varset, '%s', '%s', '%s' ) #at %s" %
    #    (it,ise,verlist,up, typconv, from_geoise, 
    # datetime.date.today() ) )
    
    #if hstring not in header.history:   header.history.append( hstring  )
    hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
    try:
        if hstring != header.history[-1]: header.history.append( hstring  )
    except:
        header.history.append( hstring  )

    #clean_run()

    #ifolder = "geo/"+it #the default place for input geometries

    #print struct_des['hs221O'] 
    nc = it+'.'+ise+typconv
    fv = verlist[0]; lv = verlist[-1];
    if typconv == "": setlist = (ise,)#
    elif up == "no_base": setlist = varset[ise].conv[typconv][1:]; up = "up1"
    else: setlist = varset[ise].conv[typconv] #
    try: conv[nc]; 
    except KeyError: 
        if typconv: conv[nc] = []    
    for inputset in setlist:

        if from_geoise:
            from_geoise = from_geoise[0]+inputset[1:] #it is supposed that the difference can be only in first digit
            ifolder = "geo/"+it+"/" + it+"."+from_geoise #+ "/" + "grainA_s" #geo used for fitted
        else: ifolder = "geo/"+struct_des[it].sfolder+"/"+it


        prevcalcver = None #version of previous calculation in verlist

        for v in verlist:
            id = (it,inputset,v)
            if typconv and id not in conv[nc]: conv[nc].append(id)
            try: blockfolder = varset[inputset].blockfolder
            except AttributeError: blockfolder = varset[ise].blockfolder
            cfolder = struct_des[it].sfolder+"/"+blockfolder #calculation folder
            add_calculation(it,inputset,v, fv, lv, ifolder, cfolder, calc, varset, up, inherit_option, prevcalcver, coord, savefile, input_geo_format)
            prevcalcver = v
    #complete_run(close_run)
    if up not in ('up1','up2'): print_and_log("You are in the test mode, please change up to up1; "); raise RuntimeError
    return



def res_loop(it, setlist, verlist,  calc = None, conv = {}, varset = {}, analys_type = 'no', b_id = (), 
    typconv='', up = 0 , imp1 = None, imp2 = None, matr = None, voronoi = False, r_id = None, readfiles = True, plot = True, show = [], 
    comment = None, input_geo_format = None, savefile = None, energy_ref = 0 ):
    """Read results
    INPUT:
    'analys_type' - ('gbe' - calculate gb energy and volume and plot it. b_id should be appropriete cell with 
        bulk material,
        'e_imp' ('e_imp_kp', 'e_imp_ecut') - calculate impurity energy - just difference between cells with impurity and without.
        'fit_ac' - fit a and c lattice constants using 2-dimensianal spline
        'clusters' - allows to calculate formation energies of clusters
        'diff' - difference of energies in meV, and volumes A^3; E(id) - E(b_id)
        )
    voronoi - True of False - allows to calculate voronoi volume of impurities and provide them in output. only if lammps is installed
    b_id - key of base calculation (for example bulk cell), used in several regimes; 
    r_id - key of reference calculation; defines additional calculation (for example atom in vacuum or graphite to calculate formation energies); can contain directly the energy per one atom

    up - controls if to download files from server
    readfiles - define if files to be readed at all or only additional analysis is required

    The next three used for 'clusters' regime:    
    imp1 - key of bulk cell with one imp1
    imp2 - key of bulk cell with one imp2
    matr - key of bulk cell with pure matrix.


    show - list, allows to show additional information (force)

    energy_ref - energy in eV; substracted from energy diffs

    """
    if type(show) == str:
        show = [show]
    if type(verlist) == int: #not in [tuple, list]:
        #print "verlist"
        verlist = [verlist]; 
    if type(setlist) == str: #not in [tuple, list,]:
        setlist = [setlist]
    if not calc:
        calc = header.calc

    """Copying files from server"""
    #to make
    # for inputset in setlist:
    #     for v in verlist:
    #         cl = calc[(it,inputset,v)]
    #         if not exists: 





    #print verlist
    if typconv == '': pass
    else: setlist = varset[setlist[0]].conv[typconv] #


    n = 'temp'; conv[n] = []
    base = 'base'; conv[base] = []
    conv[it] = []
    #print calc[b_id]

    emin = 0
    if len(b_id) == 3: # for all cases besides e_seg and coseg for wich b_id is determined every iteration
        # print "Start to read ", b_id
        if '4' not in calc[b_id].state:
            calc[b_id].read_results(1)
        e_b = calc[b_id].energy_sigma0
        v_b = calc[b_id].end.vol

    #define reference values
    e1_r = 0
    if type(r_id) in (float, int):
        e1_r = r_id
    elif type(r_id) == tuple:
        if '4' not in calc[r_id].state:
            print "start to read reference:"
            print calc[r_id].read_results(1)
            # print calc[r_id].read_results()
        e_r = calc[r_id].energy_sigma0 #reference calc
        nat_r = calc[r_id].end.natom # reference calc
        e1_r = e_r/nat_r # energy per one atom
        # print e1_r


    """Main loop"""
    final_outstring = 'no calculation found'
    for inputset in setlist:
        for v in verlist:
            # print 'Starting loops'

            id = (it,inputset,v)

            if id not in calc:
                print "Calculation does not exist!!!" 
                continue #pass non existing calculations

            cl = calc[id]
            if not hasattr(cl,'version'):
                calc[id].version = v



            path_to_contcar = cl.dir+str(v)+".CONTCAR"
            path_to_outcar = cl.dir+str(v)+".OUTCAR"
            path_to_xml   = cl.dir+str(v)+".vasprun.xml"
            # print path_to_contcar, path_to_outcar
            # print os.path.exists(path_to_contcar)
            # print os.path.exists(cl.path["output"])

            outst = ' File was not read '
            if readfiles:
                if  os.path.exists(path_to_outcar) and os.path.exists(path_to_xml):
                
                    outst = calc[id].read_results(up, analys_type, voronoi, show)
                
                else:
                    print "Trying to download OUTCAR and CONTCAR from server\n\n"
                    outst = calc[id].read_results(1, analys_type, voronoi, show)




            if analys_type in ('e_seg', 'coseg'): b_id = (b_id[0], id[1], id[2])
            if not hasattr(cl,'energy_sigma0'):
                print cl.name, 'is not finished!, continue; file renamed to _unfinished'
                outcar = cl.dir+str(v)+".OUTCAR"
                outunf = outcar+"_unfinished"
                runBash("mv "+outcar+" "+outunf)

                continue

            e = calc[id].energy_sigma0
            try:
                v = calc[id].end.vol
            except:
                v = 0
            #print e
            if e < emin: emin = e; id_min = id
            conv[n].append(id)
            # print base
            conv[base].append(b_id)


            outst2 = ("%s"%calc[id].name).ljust(22)
            outst2+='&'
            # print outst2+'&'            
            outst_end = '' 
            
            if   b_id :

                if "4" not in calc[b_id].state:    
                    calc[b_id].read_results(1)

                if calc[id].set.ngkpt != calc[b_id].set.ngkpt:
                    print_and_log("Warning! you are trying to compare calcs with "+str(calc[id].set.ngkpt)+" and "+str(calc[b_id].set.ngkpt)+"\n")
                    pass

                if calc[id].NKPTS != calc[b_id].NKPTS:
                    print_and_log("Warning! you are trying to compare calcs with "+str(calc[id].NKPTS)+" and "+str(calc[b_id].NKPTS)+" kpoints \n")

                if 'gbe' in analys_type:      
                    outst2 = gb_energy_volume(calc[id], calc[b_id])
            
                elif 'e_imp' in analys_type:
                    calc[id].e_imp = e - e_b
                    calc[id].v_imp = v - v_b

                    #calc[id].v_imp = e - e_b
                    outst2 += ("%.3f & %.2f  & %.2f &"% (e - e_b, v - v_b, (v - v_b)/v_b*100 ) )
                    conv['e_imp'].append(id)
                    a    = calc[id].hex_a;                     a_b  = calc[b_id].hex_a
                    c    = calc[id].hex_c;                     c_b  = calc[b_id].hex_c
                    ca = c/a;                                   ca_b = c_b/a_b
                    outst_end = " & {0:.1f} & {1:.1f} & {2:.1f} & {3:.3f}".format((a - a_b)/a_b*100,  (c - c_b)/c_b*100, (ca - ca_b)/ca_b*100, (e - e_b - e1_r) )

                elif analys_type == 'e_2imp': #"""For calculation of energies of two impurities in big cell"""
                    calc[id].e_imp = e - e_b            
                    outst2 += ("%.0f "% ( (e - e_b)*1000 ) )  
                    conv[it].append(id)

                elif analys_type in ('e_seg', 'coseg', 'diff'): #"""For calculation of segregation and cosegregation energies"""
                    

                    
                    e_b = calc[b_id].energy_sigma0
                    v_b = calc[b_id].end.vol


                    # outst2 += ("%.0f & %.2f "% ( (e - e_b)*1000, v - v_b ) )
                    
                    outst2 += " {:.0f} & {:.2f} ".format( (e - e_b - energy_ref)*1000, (v - v_b) ).center(6)
                    outst2 +='&'
                    # write_xyz(calc[id].end)
                    # write_xyz(calc[b_id].end)
                

            if analys_type == 'clusters':
                e1  = calc[imp1].energy_sigma0
                e2  = calc[imp2].energy_sigma0
                e_m = calc[matr].energy_sigma0
                n1 = calc[id].init.nznucl[1]
                if len(calc[id].init.nznucl) == 3:
                    n2 = calc[id].init.nznucl[2]
                else:
                    n2 = 0
                # print n1,n2
                outst2 += ("%.0f "% ( (e - n1*e1 - n2*e2 + (n1+n2-1)*e_m )*1000 / (n1+n2) ) )
            


            final_outstring = outst2+outst + outst_end              
            print final_outstring

        emin = 0
        





        """Aditional analysis, plotting"""
        if '4' not in calc[id].state:
            print "Calculation ",id, 'is unfinished;return'
            return
        final_list = () #if some part fill this list it will be returned instead of final_outstring
        if analys_type == 'gbe':
            print("\nGrain boundary energy and excess volume fit:")
            plot_conv( conv[n], calc, "fit_gb_volume")

        elif analys_type == 'gbep':
            print("\nGrain boundary energy and excess volume fit:")
            # plot_conv( conv[n], calc, "fit_gb_volume")
            final_outstring = plot_conv( conv[n], calc, "fit_gb_volume_pressure")

        elif analys_type in ('e_seg', 'coseg') and len(verlist) > 3:

            #Test lateral sizes
            A   = calc[id].end.yzarea
            A_b = calc[b_id].end.yzarea
            
            if A != A_b: 
                print_and_log("Warning! you are trying to compare calcs with different lateral sizes: "+str(A)+" "+str(A_b))
                print "Areas are ", A, A_b," A^3"
            
            #Show results 
            id1 = (it,inputset,verlist[0]) #choosen to save calculated values at first version of version set
            
            if readfiles and plot:           
                #print " \n\nImpurity at the interface :"
                e, v, emin, vmin       = plot_conv( conv[n], calc,  "fit_gb_volume2")
                #print " \n\nImpurity in the volume    :"
                e_b, v_b, e_bmin, v_bmin = plot_conv( conv[base], calc, "fit_gb_volume2")
                e_segmin = (emin - e_bmin) * 1000
                v_segmin =  vmin - v_bmin


                e_seg = (e - e_b) * 1000
                v_seg =  v - v_b

                calc[id1].e_seg = e_seg
                calc[id1].v_seg = v_seg
            
            if not hasattr(calc[id1], 'e_seg'): 
                print "Warning! Calculation ", id1, 'does not have e_seg and v_seg. Try to run with readfiles = True to calculate it.'
                calc[id1].e_seg = 0; calc[id1].v_seg = 0
            


            natom = calc[id1].natom
            calc[id1].X = 1./natom
            v1 = v / natom
            calc[id1].Xgb = v1 / A # for grain boundary with 1 A width. For other boundaries should be divided by width. 
            #print ("__________________________________________________________________________")
            
            #print (" At zero pressure: segregation energy is %.0f  meV; Seg. volume is %.1f A^3; excess seg. vol. is %.2f A" %(e_seg, v_seg, v_seg/A ) )
            # print ("%s.fit.pe & %.0f & %.1f & %.2f & %.3f & %.1f" %(id[0]+'.'+id[1], e_seg, v_seg, v_seg/A, 1./A, 1./calc[id].natom * 100  ) )
            
            #Calculate distance from impurity to boundary and number of neighbours for version 2!
            id2 =(it,inputset, 2)
            st = calc[id2].end
            gbpos2 = calc[id2].gbpos 
            gbpos1 = gbpos2 - st.rprimd[0][0]/2.
            d1 = abs(st.xcart[-2][0] - gbpos2)
            d2 = abs(st.xcart[-1][0] - gbpos2)
            dgb = d1; 
            iimp = -2
            if d2 < d1: 
                dgb = d2
                iimp = -1
            t = st.typat[iimp]
            z = st.znucl[t-1]
            segimp = element_name_inv(z) #Type of impurity closest to gb
            # print segimp, d

            id_m2   = (it+'.m',      '8'+inputset[1:], 2)
            
            if analys_type == 'e_seg':

                #calc e_seg2 and decomposition to mechanical and chemical contributions

                if id_m2 in calc: #additional analysis
                    b_id2 = (b_id[0],inputset, 2)
                    b_id_m2 = (b_id[0]+'.m', '8'+inputset[1:], 2)

                    e_seg2 = (calc[id2].energy_sigma0 - calc[b_id2].energy_sigma0) * 1000
                    e_m2   = (calc[id_m2].energy_sigma0 - calc[b_id_m2].energy_sigma0) * 1000
                    e_ch2  = e_seg2 - e_m2
                else:
                    e_seg2 = 0
                    e_m2   =0
                    e_ch2  =0


                #calculate number of close neibours around closest to gb imp
                x_central = st.xcart[iimp]
                st_r  = replic(st,   mul = (1,2,2), inv =  1 )
                st_rr = replic(st_r, mul = (1,2,2), inv = -1 ) # to be sure that impurity is surrounded by atoms

                dmax = 3
                list = [ x  for x, t  in zip(st_rr.xcart, st_rr.typat) if np.linalg.norm(x_central - x) < dmax and t == 1]
                nneigbours =  len(list)


                final_outstring = ("%s.fit.pe & %.0f & %.1f & %.2f & %.d & %4.0f & %4.0f & %4.0f & %s " %(
                    id2[0]+'.'+id2[1], calc[id1].e_seg, calc[id1].v_seg, dgb, nneigbours, e_seg2, e_ch2, e_m2, segimp  ))

                # final_outstring = ("%s.fit.pe & %.0f & %.0f & %.1f & %.1f & %.2f & %.d & %4.0f & %4.0f & %4.0f & %s " %(
                #     id2[0]+'.'+id2[1], calc[id1].e_seg, e_segmin, calc[id1].v_seg, v_segmin , dgb, nneigbours, e_seg2, e_ch2, e_m2, segimp  )) #e_segmin and v_segmin are minimum energy (but at some pressure) and corresponing volume



                final_list = [id2[0]+'.'+id2[1], calc[id1].e_seg, calc[id1].v_seg, dgb, nneigbours, e_seg2, e_ch2, e_m2, segimp]
            
            elif analys_type == 'coseg' :
                calc[id2].e_seg = calc[id1].e_seg #save in version 2
                calc[id2].v_seg = calc[id1].v_seg
                final_outstring = ("%s.fit.pe & %.0f & %.1f & %.1f & %.1f" %(id[0]+'.'+id[1], calc[id2].e_seg, calc[id2].v_seg, d1, d2 ))




            print  final_outstring
            print '\\hline'


        elif analys_type == 'e_2imp':
            # plot_conv( conv[it], calc,  analys_type, conv_ext) #instead use plot_conv( conv['hs443OO'], calc,  'e_2imp', [conv['hs443CO'],conv['hs443CC']]) 
            pass



        elif analys_type == 'fit_ac':
            #for x in calc[id_min].end.xred:
            #    print x[0], x[1], x[2]
            #print ( outst2 + " min_e & "+calc[id_min].read_results(v,) )
            #print ("name %s_template          acell  %.4f  %.4f  %.4f # fit parameters are &%.4f &%.4f &%i &%i"  % (fit_hex(0.0002,0.0003,400,600, it, inputset, verlist, calc) )  )    
            print ("name %s_template          acell  %.5f  %.5f  %.5f # fit parameters are &%.5f &%.5f &%i &%i"  % (fit_hex(0.00002,0.00003,4000,6000, it, inputset, verlist, calc) )  )    
            #name, e_min, a_min, c_min,a,b,c,d = fit_hex(0.0002,0.0003,400,600, it, inputset, verlist, calc)
            #calc[(name,'93',1)].energy_sigma0 = e_min
            #calc[(name,'93',1)].state = "4"

        elif analys_type == 'dimer':
            """Fit of md steps obtained from constant speed; see vasp description for dimer"""
            # print calc[id].list_e_sigma0
            # calc[id].list_dE = []
            # for E in calc[id].list_e_without_entr:

            #     calc[id].list_dE.append(E - 2*calc[b_id].e_without_entr)

            # print calc[id].list_e_without_entr
            # print calc[id].list_dE
            calc[id].e_ref = calc[b_id].e_without_entr

            plot_conv( [id], calc,  "dimer")

    if final_list:
        return final_list
    else:
        return final_outstring.split('&') # only for last version or fit depending on type of analysis

def inherit_icalc(inherit_type, it_new, ver_new, id_base, calc, 
    id_from = None,
    atom_new = None, atom_to_replace = None,  used_cell = 'end'):
    """
    Function for creating new geo files in geo folder based on different types of inheritance
    Input: it_new - name of new structure,
    id_base - new structure will be based on the final structure of this calculation,
    it_from - can be additionally used to adopt for example rprimd from id_from to it_new

    inherit_type = ('')
    remove_imp - removes all atoms with typat > 1
    used_cell - use init or end of id_base. Now realized only for replace atoms 



    Atoms of type 'atom_to_replace' in 'id_base' will be replaced by 'atom_new' type.

    Result: new geo file in the input geo folder


    Comments: changes len_units of new to Angstrom!!!
    """


    hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
    if hstring != header.history[-1]: header.history.append( hstring  )

    #if inherit_type not in header.history[-1] or \
    #it_new not in header.history[-1]:   header.history.append( hstring  )

    if id_from:
        if type(id_from) == str: # if string - treated like file name
            calc_from = CalculationVasp();
            calc_from.read_geometry(id_from)
            calc_from.end = calc_from.init
            calc_from_name = id_from
        else:
            print "id_from", id_from
            calc_from = calc[id_from]
            calc_from_name = calc_from.name

        if calc[id_base].len_units != calc_from.len_units:
            print_and_log("Calculations have different len_units"); raise RuntimeError



        if it_new == id_from[0] and ver_new == id_from[2]:
            print_and_log("Warning! check your versions, you are trying to overwrite existing from structures, nothing done")
            raise RuntimeError 


    if it_new == id_base[0] and ver_new == id_base[2]:
        print_and_log("Warning! check your versions, you are trying to overwrite existing base structures, nothing done")
        raise RuntimeError
  


    new = copy.deepcopy(  calc[id_base]  )

    new.len_units = 'Angstrom' #! Because from VASP


    new.path["input_geo"] = geo_folder + struct_des[it_new].sfolder + '/' + \
        it_new+"/"+it_new+'.inherit.'+inherit_type+'.'+str(ver_new)+'.'+'geo'
    print new.path["input_geo"]
    new.version = ver_new

    if inherit_type == "r2r3":
        des = ' Partly inherited from the final state of '+calc[id_base].name+'; r2 and r3 from '+calc_from_name
        new.des = struct_des[it_new].des + des
        new.end.rprimd[1] = calc_from.end.rprimd[1].copy()
        new.end.rprimd[2] = calc_from.end.rprimd[2].copy()       
        new.write_geometry("end",des)

    elif inherit_type == "r1r2r3":
        des = ' Partly inherited from the final state of '+calc[id_base].name+'; r1, r2, r3 from '+calc_from_name
        new.des = struct_des[it_new].des + des
        new.end.rprimd = copy.deepcopy( calc_from.end.rprimd )
        new.hex_a = calc_from.hex_a
        new.hex_c = calc_from.hex_c
        new.end.xcart = xred2xcart(new.end.xred, new.end.rprimd) #calculate new xcart from xred, because rprimd was changed
        new.write_geometry("end",des)


    elif inherit_type == "full":
        print_and_log("Warning! final xred and xcart was used from OUTCAR and have low precision. Please use CONTCAR file \n");
        des = 'Fully inherited from the final state of '+calc[id_base].name
        new.des = des + struct_des[it_new].des
        new.write_geometry("end",des)

    elif inherit_type == "remove_imp":
        """
        Assumed that typat == 1 is matrix atoms
        """
        des = 'All impurities removed from the final state of '+calc[id_base].name
        new.des = des + struct_des[it_new].des
     
        st = calc[id_base].end
        new.end.typat = []
        new.end.xred = []
        new.end.xcart = []
        new.end.ntypat = 1
        new.end.znucl = new.end.znucl[0:1]
        for i, t in enumerate(st.typat):
            if t == 1:
                new.end.typat.append(t)
                new.end.xred.append(st.xred[i])
                new.end.xcart.append(st.xcart[i])
        new.end.natom = len(new.end.xred)
        new.init = new.end #just for sure
        new.write_geometry("end",des)        


    elif inherit_type == "replace_atoms":
        """Simply replace in calculation one atoms by another """
        z_new     = element_name_inv(atom_new)
        z_replace = element_name_inv(atom_to_replace)


        if used_cell == 'init':
            st = new.init
        elif used_cell == 'end':
            st = new.end


        znucl = st.znucl
        
        if z_replace not in znucl: 
            print "Error! Calc "+new.name+" does not have atoms of this type\n"
            raise RuntimeError

        if atom_to_replace not in id_base[0] or atom_new not in it_new:
            print "Error! Something wrong with names of atom types\n"
            raise RuntimeError            
        print "replace ", z_replace, "by", z_new

        znucl = [int(z) for z in znucl] # convert to int
        i_r = znucl.index(z_replace)
        print "index ", i_r
        znucl[i_r] = z_new

        if znucl[-2] == znucl[-1]: #just for special case
            del znucl[-1]

            for i, t in enumerate(st.typat):
                if t == st.ntypat:
                    print "found ntypat" ,t
                    st.typat[i]-=1
                    #print t
            st.ntypat-=1
            #print st.typat
            #print new.end.typat
        
        st.znucl = znucl



        des = 'Fully inherited from the '+ used_cell +' state of '+calc[id_base].name+\
        ' by simple replacing of '+atom_to_replace+' by '+atom_new

        new.des = des + struct_des[it_new].des
        new.write_geometry(used_cell, des)

    else:
        print_and_log("Error! Unknown type of Calculation inheritance"); raise RuntimeError
    return







def for_phonopy(new_id, from_id = None, calctype = 'read', mp = [10, 10, 10], additional = None):
    #creates file for phonopy, run phonopy
    #new_id - will add this calculation or read; if string then interpreted as filename of thermal_properties.yaml
    #from_id - tuple - than will take coordinates from the end; or path to input poscar file
    #type - 'create', 'read'
    #additional - list of calculation names, if calculation was splited into several parts

    mpstr = " ".join(map(str, mp))

    def read_phonopy_data(filename, key = "free_energy" ):
        # with open(filename, 'r') as f:
        #     f.readline()
        F = []
        T = []
        #     for line in f:
        #         T.append(float(line.split()[0]))
        #         F.append(float(line.split()[1]))
        #     # print T, F
        import yaml
        f = open(filename)
        # use safe_load instead load
        dataMap = yaml.safe_load(f)
        f.close()
        prop = dataMap['thermal_properties']
        for i in range(len(prop)):
            T.append( prop[i]['temperature'] )
            F.append( prop[i][key]     )

        coeffs1 = np.polyfit(T, F, 8)
        fit_func = np.poly1d(coeffs1)
        T_range = np.linspace(min(T), max(T))
        print 'I return', key

        return T_range, fit_func




    if calctype == 'create':
        work_path = 'geo/'+struct_des[new_id[0]].sfolder+'/'+new_id[0]
    
        log_history(  "{:}    #on {:}".format( traceback.extract_stack(None, 2)[0][3],   datetime.date.today() )  )

        print type(from_id)
        if type(from_id) == str:
            from_cl = CalculationVasp()
            from_cl.read_poscar(from_id)
            state = 'init'
            from_st = from_cl.init
        else:
            from_cl = header.calc[from_id]
            state = 'end'
            from_st = from_cl.end

            # print  header.varset
            #create POSCAR
        posname = 'phonopy_input_poscar'
        from_cl.write_structure(name_of_output_file = posname, path = work_path, state = state)

        savedPath = os.getcwd()
        os.chdir(work_path)
        #create conf 
        confname = new_id[0]+'.conf'
        with open(confname, 'w') as f:
            f.write("DIM = 1 1 1\n")
            f.write("ATOM_NAME = ")
            for z in from_st.znucl:
                el = element_name_inv(z)
                f.write(el+' ')
            f.write("\nDIAG = .TRUE.\n")
            f.write("DISPLACEMENT_DISTANCE = 0.03\n")

         
        #run phonopy
        print runBash('export PYTHONPATH=~/installed/phonopy-1.9.5/lib/python:$PYTHONPATH; rm POSCAR-*; /home/dim/installed/phonopy-1.9.5/bin/phonopy '
            +confname+' -c '+posname+' -d --tolerance=0.01')

        ndis = len( glob.glob('POSCAR-*') )
        print ndis, ' displacement files was created'

        os.chdir(savedPath)


        #add
        add_loop(new_id[0],   new_id[1], range(1,ndis+1), up = 'up1', input_geo_format = 'vasp', savefile = 'allx')

        #copy SPOSCAR - an ideal cell
        src =  'geo/'+struct_des[new_id[0]].sfolder+'/'+new_id[0]
        dst = struct_des[new_id[0]].sfolder+'/'+new_id[0]+'.'+new_id[1]+'/'
        shutil.copy(src+'/SPOSCAR', dst)
        shutil.copy(src+'/disp.yaml', dst)

    if calctype == 'read':

        if type(new_id) == tuple:
            new_cl = header.calc[new_id]

            work_path = struct_des[new_id[0]].sfolder+'/'+new_id[0]+'.'+new_id[1]
            work_path_geo = 'geo/'+struct_des[new_id[0]].sfolder+'/'+new_id[0]



            npos = len( glob.glob(work_path+'/*.POSCAR') )
            print range(1,npos+1), 'range'
            if not os.path.exists(work_path+"/1.POSCAR"):
                res_loop(new_id[0],   new_id[1], range(1,npos+1), up = 'up1', input_geo_format = 'vasp', savefile = 'allx')

            if additional:
                for name in additional:
                    npos_new = len( glob.glob(struct_des[new_id[0]].sfolder+'/'+name+'.'+new_id[1]+'/*.POSCAR') )

                    if not os.path.exists(struct_des[new_id[0]].sfolder+'/'+name+'.'+new_id[1]+'/'+str(npos+1)+".POSCAR"):
                        res_loop(name,   new_id[1], range(npos+1, npos+npos_new+1), up = 'up1', input_geo_format = 'vasp', savefile = 'allx')
                    
                    npos = npos+npos_new

                    runBash("rsync "+struct_des[new_id[0]].sfolder+'/'+name+'.'+new_id[1]+'/*.vasprun.xml '+work_path)
                    # print 'Additional vasprun.xml files were copied to ', work_path

            savedPath = os.getcwd()
            os.chdir(work_path)


            #create conf 
            confname = new_id[0]+'_mesh.conf'
            with open(confname, 'w') as f:
                f.write("DIM = 1 1 1\n")
                f.write("ATOM_NAME = ")
                for z in new_cl.end.znucl:
                    el = element_name_inv(z)
                    f.write(el+' ')
                f.write("\nMP = {:}\n".format( mpstr ))
                f.write("\nTSTEP = {:}\n".format( 1 ))
                f.write("\nTMAX = {:}\n".format( 1155 ))

            if not os.path.exists("FORCE_SETS"):
                #run phonopy; read forces
                ndis = len( glob.glob('*.vasprun.xml') )
                print ndis, ' displacement files was Found'
                print runBash('export PYTHONPATH=~/installed/phonopy-1.9.5/lib/python:$PYTHONPATH; /home/dim/installed/phonopy-1.9.5/bin/phonopy '
                    +'  -f {1..'+str(ndis)+'}.vasprun.xml --tolerance=0.01')

            #calculate thermal prop
            result = 'thermal_properties_'+mpstr.replace(" ", "_")+'.yaml'
            if not os.path.exists(result):

                posname = 'SPOSCAR'
                print runBash('export PYTHONPATH=~/installed/phonopy-1.9.5/lib/python:$PYTHONPATH; /home/dim/installed/phonopy-1.9.5/bin/phonopy '
                    +confname+' -c '+posname+' -t -p -s --tolerance=0.01')

                shutil.copyfile('thermal_properties.yaml', result)
    

            T_range, fit_func = read_phonopy_data(result)

            os.chdir(savedPath)

        if type(new_id) == str:
            result = new_id
            T_range, fit_func = read_phonopy_data(result)



    return T_range, fit_func









#Deprecated classes, needed only for compatibility
class DataStuff:
    """Parent Class For reading data and manipulate"""
    def do(self,command):
        self.s = runBash(command)
        ch = re.compile('[a-zA-Z]+') # make object , matched anything other than digit
        #find first entrance of any word
        try:
            w = ch.match(self.s).group()
            self.sdig = self.s.replace(w, '')#delete
            self.sv =  self.sdig.split() #values but like strings
            self.v = []
            for i in range(len(self.sv)):
                self.v.append(float(self.sv[i]))  
        except AttributeError:
            self.sdig = self.s
                
       
class CartArray(DataStuff):
    "Objects of this class will storage lists of 3 cartesian components, like list of atoms coordinates "
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []
        self.r = []
    def do(self,command):
        DataStuff.do(self,command)
        if not self.sdig == '':
            lines = self.sdig.splitlines() #divide xcart string by lines
            for i in range(len(lines)):
                w = lines[i].split() #divide each line on components 
                if not len(w) == 3: 
                    print "Too many components, check sdig attribute - allowed only 3"
#                print self.s
                
                
                x = float(w[0]); y = float(w[1]); z = float(w[2]);
                self.x.append(x)
                self.y.append(y)
                self.z.append(z)
                self.r.append( math.sqrt(x*x + y*y + z*z) )









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
        printlog('self.path["input_geo"] ', self.path["input_geo"], imp='y')
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
>>>>>>> master
        """
        Convinient wrapper for add_loop()

        TODO: to be finished; It is more reasonable to move this method to class; Then a new calculation is created which requires a structure
        and set; But probably it can be here as well. It will be very simple to use. 

        please make it in such a way that a res command is suggested at the end, required to read the results
        """
        return









<<<<<<< HEAD
=======




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
            elif schedule_system == 'simple':
                pass
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


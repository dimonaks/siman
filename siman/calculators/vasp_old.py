
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
            path_to_ibzkpt  = path_to_outcar.replace('OUTCAR', "IBZKPT")
            path_to_wavecar = path_to_outcar.replace('OUTCAR', "WAVECAR")
            path_to_doscar = path_to_outcar.replace('OUTCAR', "DOSCAR")
            path_to_eigenval = path_to_outcar.replace('OUTCAR', "EIGENVAL")
            path_to_procar = path_to_outcar.replace('OUTCAR', "PROCAR")
            path_to_locpot = path_to_outcar.replace('OUTCAR', "LOCPOT")
            path_to_kpoints = path_to_outcar.replace('OUTCAR', "KPOINTS")
            path_to_waveder = path_to_outcar.replace('OUTCAR', "WAVEDER")  
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

        if 'w' in load:

            self.get_file(os.path.basename(path_to_wavecar), up = load)
            self.get_file(os.path.basename(self.dir+'/WAVECAR'), up = load)
        if 'b' in load:

            self.get_file(os.path.basename(path_to_ibzkpt), up = load)

        if 'e' in load:
            self.get_file(os.path.basename(path_to_eigenval), up = load)
            self.get_file(os.path.basename(self.dir+'/EIGENVAL'), up = load)

        if 'd' in load:
            self.get_file(os.path.basename(path_to_doscar), up = load)
            self.get_file(os.path.basename(self.dir+'/DOSCAR'), up = load)

        if 'l' in load:
            self.get_file(os.path.basename(path_to_locpot), up = load)
            self.get_file(os.path.basename(self.dir+'/LOCPOT'), up = load)

        if 'pr' in load:

            # self.get_file(os.path.basename(path_to_procar), up = load)
            self.get_file(os.path.basename(self.dir+'/PROCAR'), up = load)
        if 'k' in load:

            self.get_file(os.path.basename(self.dir+'/KPOINTS'), up = load)

        if 'wd' in load:
            self.get_file(os.path.basename(path_to_waveder), up = load)
            self.get_file(os.path.basename(self.dir+'/WAVEDER'), up = load)



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

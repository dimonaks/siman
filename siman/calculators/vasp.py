# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
import os, math, copy, glob, shutil, sys
import numpy as np
from siman import header
from siman.core.calculation import Calculation
from siman.core.structure import Structure
from siman.header import runBash

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




class CalculationVasp(Calculation):
    """Methods for calculations made using VASP DFT code"""
    def __init__(self, inset = None, iid = None, output = None):
        self.calculator = 'vasp'
        super(CalculationVasp, self).__init__(inset, iid, output)
        self.len_units = 'Angstrom'
        self.init = Structure()
        self.end = Structure()



    def set_output_filenames(self, out_name, version):
        cl = self

        if out_name:
            cl.path["output"] = cl.dir+out_name
        else:
            name_mod = ''
            cl.path["output"] = cl.dir+str(version)+name_mod+".OUTCAR" #set path to output
        
        #paths to other files
        cl.path["charge"] = cl.path["output"].replace('OUTCAR', 'CHGCAR')

        return cl.path["output"]


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



    def update_incar(self, parameter = None, value = None, u_ramp_step = None, write = True, f = None, run = False,):    
        """Modifications of INCAR. Take attention that *parameter* will be changed to new *value*
        if it only already exist in INCAR.  *u_ramp_step*-current step to determine u,
        *write*-sometimes just the return value is needed. 
        Returns U value corresponding to *u_ramp_step*.
        """


        # self = st 
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

        # elif parameter in ['IMAGES', 'ISPIN']:
        else:

            new_incar_string = parameter + ' = ' + str(value)
            command = "sed -i.bak '/"+parameter+"/c\\" + new_incar_string + "' INCAR\n"




        if write and f:
            f.write(command)

        if run:
            runBash(command)

        return  u_step #for last element





    
    def make_kpoints_file(self):

        struct_des = header.struct_des
        #Generate KPOINTS
        kspacing = self.set.vasp_params.get('KSPACING')

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
            if self.set.kpoints_file is True:

                print_and_log( "You said to generate KPOINTS file. set.kpoints_file =",self.set.kpoints_file," \n")
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

            self.path['kpoints'] = filename 


        else:
            print_and_log( "This set is without KPOINTS file.\n")
            filename = ''



        return [filename]



    def actualize_set(self, curset = None, params = None):
        """
        Makes additional processing of set parameters, which also depends on calculation
    
        adding parameters for atat


        TODO: 
            
            if number is used after element, then it should be considered as separate type

        RETURN:
            
            None

        AUTHOR:

            Aksyonov DA
        """
        vp = curset.vasp_params
        if not curset:
            curset = self.set


        #check if some parameters should be filled according to number of species
        #make element list
        el_list = [element_name_inv(el) for el in self.init.znucl] # length is equal to number of elements
        ldau = False
        if 'LDAUL' in vp and vp['LDAUL'] is not None: 
            ldau = True

            if 1:
                'This block checks if more types for one element should be added to structure'
                uels_set = vp['LDAUL'].keys()
                # print(uels_set)
                uels = [] # list of elements provided by set relevant for the given structure
                uelntypat = {} # number of types for each element from the structure, for which U should be used
                for el in set(el_list):
                    for el_set in uels_set:
                        if el == el_set.split('/')[0]: # after slash the required coordination is given.
                            uels.append(el_set)
                            if el not in uelntypat:
                                uelntypat[el] = 0
                            uelntypat[el] +=1
                # print(uels, uelntypat)
                anions_set = []
                for el in uelntypat.keys(): #
                    if uelntypat[el] > 1: # this condition shows that multitype regime was asked for
                        printlog('LDAUL are given for the following multitype elements: ',uels, imp = 'y')
                        for el_set in uels:
                            if el in el_set and '/' not in el_set: # check that correct format is used
                                printlog('Error! Element', el, f'has several types. LDAUL values should be given in format {el}/A but mixture of {el}/A and {el} was detected. Please correct')
                            #check that required anions are present 
                            A = el_set.split('/')[1]
                            anions_set.append(A)
                            # print(A)
                            if A not in el_list:
                                printlog(f'Warning! Anion {A} is absent in your structure but present in your set')

                        printlog(f'Checking if more types for {el} should be added ... ', imp = 'y')
                        self.init, cords = self.init.add_types_for_el(el)
                        if len(cords) == 1:
                            printlog(f'Error! This structure has only one symmetry non-equivalent position for {el} and incompatible with the chosen set' )
                        
                        if len(cords) != uelntypat[el]:
                            # print( len(cords), uelntypat[el] )
                            printlog(f'Error! Number of non-equivalent position for {el} is {len(cords)}, which is incompatible  with {uelntypat[el]} coordinations provided the chosen set' )

                        inter = list(set([item for sublist in cords for item in sublist]).intersection(set(anions_set))  ) 
                        if len(anions_set) != len(inter):
                            printlog('Warning! Coordinations in your structure and provided LDAUL are incompatible ')
                        # sys.exit()
                        # for A in anions_set:
                        #     for cord in cords:
                        #         print(A, cord)
                        #         if A in cords:
                        #             break
                        #     else:


            aniels = [invert(a) for a in header.ANION_ELEMENTS]
            # set(aniels)
            nintersections = len(set(aniels).intersection(el_list))
            # print(nintersections, aniels, el_list)

            el_list = self.init.get_unique_type_els(True) # new list after adding 
            # sys.exit()

            for key in ['LDAUL', 'LDAUU', 'LDAUJ']:
                # print( vp[key])
                try:
                    if set(vp[key].keys()).isdisjoint(set(el_list)): #no common elements at all
                        printlog('\n\n\nAttention! The '+str(key)+' doesnt not contain values for your elements! Setting to zero\n\n\n')
                        # raise RuntimeError

                    new = []
                    for el in el_list:
                        
                        if el in vp[key]:
                            val = vp[key][el]
                            
                            for A in aniels:
                                if A in el_list:  # use another U value for other anions, provided like Fe/S 
                                    kk = el+'/'+A 
                                    if kk in vp[key]:
                                        if nintersections == 1:
                                            val = vp[key][kk]
                                        else:
                                            printlog(f'Error! The chosen structure has more than one anion. The given U values {el}/{A} will be used for all {el} atoms')





                        else:

                            if key == 'LDAUL':
                                val = -1
                            else:
                                val =  0
                            if '/' in el:
                                printlog(f'Warning! no value is given in {key} for element {el}')

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
            for el in curset.magnetic_moments.keys():
                if '/' in el and ldau == False:
                    # if ldau is true and multitype regime is true then everything is updated above
                    # if ldau is true and multitype is false than we cant use multitype here because LDAUL will become incompatible with POSCAR 
                    self.init, cords = self.init.add_types_for_el(el.split('/')[0])
                    el_list = self.init.get_unique_type_els(True) # new list after adding 
                else:
                    ''

            for iat in range(self.init.natom):
                typ = self.init.typat[iat]
                el  = el_list[typ-1]
                if el in curset.magnetic_moments:
                    magmom.append(curset.magnetic_moments[el])
                else:
                    if '/' in el:
                        printlog(f'Warning! {el} was not found in your set, I use {mag_mom_other} for it')
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
        """
        Plot energy versus force for each ionic relaxation step.
        
        INPUT:
            
            - force_type (str) 
                - 'max' - maximum force
                - 'av'  - average force

        RETURN:
            
            None

        AUTHOR: 

            Aksyonov D.A.
        """

        import numpy as np
        from siman.picture_functions import fit_and_plot

        if 'max' in force_type:
            forces = [m[1] for m in self.maxforce_list ]
            lab = 'Max.'

        elif 'av' in force_type:
            forces = [m for m in self.average_list ]
            lab = 'Av.'

        energies = 1000*(np.array(self.list_e_sigma0)-self.energy_sigma0) # relative energies in meV
        numbers = list(range(len(forces)))
        annotates = []
        for n in numbers:
            if n%5 == 0 or n == len(numbers)-1:
                if n == 0:
                    annotates.append(str(n)+' first')
                elif n == len(numbers)-1:
                    annotates.append(str(n)+' last')
                else:
                    annotates.append(n)
            else:
                annotates.append('')



        fit_and_plot(data = {'x': forces, 'y':energies, 'fmt':'-o', 
            'annotates': annotates, 'annotate_fontsize':10, 'annotate_arrowprops':None}, annotate = 1,
            xlabel = lab+' force on atom (meV/$\AA$)', ylabel = 'Energy per cell relative to min (meV)',
            show = 1)


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
        """just wrapper to get chgcar files 

        - nametype (str)
            - 'asoutcar', e.g. 1.CHGCAR
            - '', CHGCAR

        """
        if 'nametype' not in kwargs:
            kwargs['nametype'] = 'asoutcar'
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

        P2A = header.PATH2ARCHIVE

        if P2A:
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
        
        if P2A:
        
            restore_CHG = "rsync "+CHG_scratch_gz+' '+path+' ; gunzip '+CHG+'.gz ' # on cluster
            restore_AEC = "rsync "+P2A+'/'+self.path['aeccar0']+' '+ P2A+'/'+self.path['aeccar2']+' '+path # on cluster


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
            
            if (remote(no_AECCAR0) or remote(no_AECCAR2)) and P2A:
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

    def band(self, ylim = None):
        """
        Plot band structure using pymatgen. Can be applied only to band structure calculations with correct KPOINTS and vasprun.xml files
        INPUT:
            - ylim (2*tuple of float)
        """

        from pymatgen.io.vasp import Vasprun, BSVasprun
        from pymatgen.electronic_structure.plotter import BSPlotter

        xml_file = self.get_file('vasprun.xml', nametype = 'asoutcar')
        # print(xml_file)
        # print(self.path['kpoints']) #os.getcwd()+'/'+
        if ylim is None:
            ylim = (-12,6)

        # sys.exit()

        if not self.path.get('kpoints'):
            self.path['kpoints'] = self.dir+'/KPOINTS'

        v = BSVasprun(xml_file)
        bs = v.get_band_structure(kpoints_filename = self.path['kpoints'], line_mode = True)
        plt = BSPlotter(bs,)
        ax = plt.get_plot(vbm_cbm_marker=True, ylim = ylim,)
        # ax.legend = None
        ax.legend().remove()
        ax.savefig('figs/png/'+str(self.name)+'_band.png', dpi = 300)
        ax.savefig('figs/'+str(self.name)+'_band.pdf')
        # ax.show()
        return
    def get_band_info(self):
        """
        Vasprun.xml is expected

        TODO
        can be done with EIGENVAL, please make

        """
        from pymatgen.io.vasp import Vasprun
        xml_file = self.get_file('vasprun.xml', nametype = 'asoutcar')
        # print('xml_file',xml_file)
        v = Vasprun(xml_file)
        ll = list(v.eigenvalue_band_properties)
        print('band gap = {:.1f} eV, cbm = {:.1f} eV, vbm = {:.1f} eV, is_band_gap_direct = {:}'.format(*ll))









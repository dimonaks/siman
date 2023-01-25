# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
from __future__ import division, unicode_literals, absolute_import, print_function
import json, sys
import copy

"""
Oграничения режима sequence_set:
1) OCCMATRIX не копируется для дочерних сетов
2) Режим U-ramping выключен для дочерних сетов
3) Есть еще, режим afm_ordering, возможно neb 
4) kpoints file только для первого сета
5) u-ramping inherit_xred - могут быть проблемы более чем для двух сетов


TODO:
ngkpt_dict_for_kspacings - when ngkpt is used could be problems, please test.

"""



from siman import header
from siman.header import print_and_log, printlog;
from siman.small_functions import is_list_like, red_prec, is_string_like
from siman.functions import invert

#Vasp keys
vasp_electronic_keys = [
'ALGO',
'PREC',
'LREAL',
'ENCUT',
'ENAUG',
'ISMEAR',
'SIGMA',
'EDIFF',
'NELM',
'NELMIN',
'NELMDL',
'MAXMIX',
'NELECT'
]

vasp_ionic_keys = [
'IBRION',
'ISIF',
'NSW',
'EDIFFG',
'POTIM',
'POMASS',
'ZVAL',
'SMASS'

]


vasp_other_keys = [
'AGGAC',
'LUSE_VDW',
'PARAM1',
'PARAM2',
'LVDW',
'LVHAR',
'LCALCPOL',
'EFIELD',
'VDW_RADIUS',
'VDW_SCALING',
'VDW_CNRADIUS',
'LWANNIER90_RUN',
'IVDW',
'VDW_D',
'MDALGO',
'TEBEG',
'TEEND',
'SYSTEM',
'ISTART',
'ICHARG',
'KGAMMA',
'KSPACING',
'EFIELD_PEAD',
'LPLANE',
'LSEPC',
'LSEPB',
'OMEGAMAX',
'ENCUTGW',
'NBANDSGW',
'NBANDSO',
'NBANDSV',
'ANTIRES',
'NOMEGA',
'OMEGATL',
'NCORE',
'NPAR',
'LSCALU',
'NSIM',
'ISYM',
'SYMPREC',
'LORBIT',
'EMIN',
'EMAX',
'NEDOS',
'LAECHG',
'LSORBIT',
'SAXIS',
'ISPIN',
'NBANDS',
'PSTRESS',
'ADDGRID',
'MAGMOM',
'GGA_COMPAT',
'IMAGES',
'LDAU',
'LDAUTYPE',
'LDAUL',
'LDAUU',
'LDAUJ',
'LDAUPRINT',
'LASPH',
'LMAXMIX',
'NFREE',
'AMIX',
'BMIX',
'AMIX_MAG',
'BMIX_MAG',
'WC',
'MAXMIX',
'OCCDIR1',
'OCCEXT',
'LHFCALC',
'HFSCREEN',
'TIME',
'PRECFOCK',
'NKRED',
'NGX',
'NGY',
'NGZ',
'NBMOD',
'LPARD',
'EINT',
'LWAVE',
'GGA',
'IALGO',
'LSCALAPACK',
'AMIN',
'IDIPOL',
'LDIPOL',
'DIPOL',
'LVTOT',
'AEXX',
'LDIAG',
'METAGGA',
'CMBJB',
'CMBJA',
'IMIX',
'LPEAD',
'LEPSILON',
'LCALCEPS',
'CSHIFT',
'LOPTICS',
'LRPA',
'LSPECTRAL',
'LCHARG',
'LELF',
'RWIGS',
'NUPDOWN',
'ALDAC',
'LAMBDA',
'SUBATOM',
'KPPRA',
'LAMBDA_D_K',
'USEPOT',
'M_CONSTR',
'I_CONSTRAINED_M',
'RWIGS',
'LSOL',
'EB_k',
'TAU',
'CORE_C',
'EB_K',
'LVDW_EWALD',
'NWRITE',
'MAGATOM',
'DOSTATIC',
'IWAVPRE'
]
vasp_keys = vasp_electronic_keys+vasp_ionic_keys+vasp_other_keys

siman_keys = [
'universal', # universal paramater with any content
'u_ramping_region', #deprecated
'u_ramping_nstep', #number of u ramping steps
'magnetic_moments',
'afm_ordering',
'set_sequence',# sequence of sets
'savefile', #additional keys pointing which files should be saved
'k_band_structure', # list, first position is number of points, then high-symmetry k-points in the form ['G', 0, 0, 0] in reciprocal space for calculating band structure 
'path2pot', # path to folder with potentials - used with potdir; if not provided that header.path2potentials is used
'path_to_potcar', # explicit path to potential - depreacated
'periodic', # 1 or 0, periodic boundary conditions or not; by default considered periodic
]

aims_keys = [
'k_grid',
'default_initial_moment',
'spin',
]

'put here quantum espresso keys'
qe_keys = [
]

'gaussian keys'
gaussian_keys = [
'functional',
'basis_set',
'job_type',
'optional',
'multiplicity',
'charge',
'SCRF',
]


def read_vasp_sets(user_vasp_sets, override_global = False):
    """
    Read user sets for different calculators and add them to project database
    Works not only for VASP but other codes as well, such as Gaussian
    
    INPUT:
        - varset (dict) - database dict with all sets of a project
        - user_vasp_sets (list) - list of user sets that describes creation of new sets based on inheritance 
        - override - allows to recreate all sets; can be usefull than you want to add some new property to all your sets - very dangerous to do!

    
    RETURN:
        - user_vasp_sets (list)

    """


    varset = header.varset

    bfolder = '' #by default no blockfolder
    for l in user_vasp_sets:
        if override_global or 'over' in l[-1]: 
            override = True
        else:
            override = False
        
        if override or l[0] not in varset:
            # print override, 'override'
            param = l[2]

            if 'bfolder' in param:
                bfolder = param['bfolder']
            else:
                bfolder = None


            s = inherit_iset(l[0], l[1], varset, override = override, newblockfolder = bfolder) 
            # print ('param', param,)
            s.load(param, inplace = True)
            
        header.varset = varset

    return varset




class InputSet():
    """docstring for InputSet
    The second important class which is used to store 
    parameters of calculation

    For VASP parameters *self.vasp_params* dict is used;
    usually it contains the parameters in the same format as INCAR file.
    However, several exceptions are:
    for 'LDAUU', 'LDAUJ', 'LDAUL' you should provide
    dictionaries with correponding values for each element in the form: {'Co':3.4,}.
    
    self.potdir (dict) - name of POTCAR folder for each element, for example {3:'Li', 8:'O'}

    self.blockfolder (str) - additional subfolder will be created calculation with this set

    self.save_last_wave (bool) - set True to save last WAVECAR in u-ramping mode

    self.kpoints_file - if True, k-points file is created, if string then it is considered as path to external kpoints file

    self.path_to_potcar (str) - explicit path to potcar, can be used instead of self.potdir

    self.set_sequence (list)  - list of InputSet() objects to make multiset runs. The current set is used as a first one.


    TODO
    Describe the difference between update() and load() methods !


    """
    def __init__(self, ise = None, path_to_potcar = None, calculator = "vasp"):
        #super(InputSet, self).__init__()
        self.ise = ise
        self.name = ise
        self.des   = "" # description
        self.potdir = {}
        self.units = "vasp"
        self.calculator = calculator # 'vasp', 

        if self.calculator is None:
            printlog('Error! Please provide calculator type!')

        self.params = {}
        self.vasp_params = self.params # params for any code!
        self.mul_enaug = 1
        self.history = "Here is my uneasy history( :\n"    
        self.tsmear = None 
        self.tolmxf = None
        self.ngkpt  = None
        self.blockfolder = ''
        self.set_sequence = None
        self.kpoints_file = None # can be path to external file
        self.save_last_wave = None #if True than do no remove last wavefunction
        self.periodic = 1 # PBC
        # self.use_ngkpt = False
        
        if path_to_potcar:
            self.path_to_potcar = path_to_potcar
        else:
            self.path_to_potcar = None


        #Code scpecific parameters, now only for Vasp
        for key in vasp_electronic_keys: 
            self.vasp_params[key] = None 
        for key in vasp_ionic_keys: 
            self.vasp_params[key] = None 
        for key in vasp_other_keys: 
            self.vasp_params[key] = None 

        for key in aims_keys: 
            self.params[key] = None 

        for key in gaussian_keys: 
            self.params[key] = None 


        #add to varset
        # if ise not in header.varset:
        header.varset[ise] = self




    def printme(self):
        """
        Print set

        """
        if hasattr(self, 'vasp_params') and self.vasp_params:
            self.params = self.vasp_params
            self.calculator = 'vasp'
        for key in self.params:
            if self.params[key] == None: 
                continue
            printlog( "{:30s} = {:s} ".format("s.params['"+key+"']", str(self.params[key]) ), imp = 'Y', end = '\n' )


        printlog('ngkpt:', self.ngkpt, imp = 'Y')
        # print(self.calculator)

        if hasattr(self, 'calculator') and self.calculator == 'vasp':
            printlog('POTDIR:', self.potdir, imp = 'Y', end = '\n' )

        return

    def update(self):
        #deprecated, but still can be usefull
        # print_and_log('Updating set ...\n')
        # c1 = 1; c2 = 1
        # if self.units == "abinit":
        #     c1 = to_eV
        #     c2 = Ha_Bohr_to_eV_A
        # #Update Vasp parameters
        # if self.units == "vasp":
        #     c1 = 1
        #     c2 = 1
        # if self.ecut == None:
        #     self.vasp_params['ENCUT'] = None
        #     self.vasp_params['ENAUG'] = None
        # else:           
        #     self.vasp_params['ENCUT'] = self.ecut * c1* self.dilatmx * self.dilatmx
        #     self.vasp_params['ENAUG'] = self.mul_enaug * self.vasp_params['ENCUT']
        # self.vasp_params['SIGMA'] = self.tsmear * c1
        vp = self.vasp_params

        self.tsmear = vp.get('SIGMA')
        self.tolmxf = vp.get('EDIFFG')

        if self.tolmxf and self.tolmxf < 0:
            self.tolmxf*=-1

        self.toldfe = vp.get('EDIFF')
        # self.vasp_params['EDIFF'] = self.toldfe * c1
        # self.vasp_params['NELM'] = self.nstep
        # self.vasp_params['NSW'] = self.ntime
        # self.vasp_params['EDIFFG'] = -self.tolmxf * c2
        self.kspacing = vp.get('KSPACING')
        self.ecut     = vp.get('ENCUT') 

        # print (self.vasp_params)
        if 'LDAUU' in self.vasp_params and self.vasp_params['LDAUU']:
            self.dftu = True
        else:
            self.dftu = False

        if 'ISPIN' in self.vasp_params and self.vasp_params['ISPIN'] == 2:
            self.spin_polarized = True
        else:
            self.spin_polarized = False


    def load(self,param, inplace = False):
        """
        Update parameters of set from dict param
        """
        # print(param)

        if inplace:
            s = self
        else:
            s = copy.deepcopy(self)
            
        for key in param:
            # print(key)
            if key in vasp_keys:
                s.set_params_dict(key, param[key])

            elif key == 'set_potential':
                for key2 in param[key]:
                    s.set_potential(key2, param[key][key2])

            elif key == 'add_nbands':
                s.set_add_nbands(param[key])

            elif key == 'ngkpt':
                s.set_ngkpt(param[key])

            elif key == 'kpoints_file':
                if param[key] in [1, True]:
                    s.kpoints_file = True 
                elif is_string_like(param[key]):
                    s.kpoints_file = param[key]
                else:
                    s.kpoints_file = False

            elif key == 'bfolder':
                print_and_log( 'New blockfolder', param[key])

            elif key in siman_keys:
                s.set_attrp(key, param[key] )

            elif key in aims_keys:
                s.set_params_dict(key, param[key] )
            
            elif key in gaussian_keys:
                # print(key )
                s.set_params_dict(key, param[key] )

            
            else:
                printlog('Error! Unknown key: '+key)
                raise RuntimeError
         

            if key == 'set_sequence':
                sets = []
                for se in s.set_sequence:
                    sets.append(copy.deepcopy(header.varset[se]))

                s.set_sequence = sets  #put objects instead of names

            # if hasattr(s, 'set_sequence') and s.set_sequence:
            #     sets = []
            #     for se in s.set_sequence:
            #         if type(se) == str:
            #             sets.append(copy.deepcopy(varset[se]))
            #         else:
            #             sets.append(copy.deepcopy(se))

            #     s.set_sequence = sets  #put objects instead of names
        return s






    def read_universal(self, filename):
        #read any file to univeral parameter
        with open(filename, 'r') as f:
            fil = f.read()
            self.params['universal'] = fil        



    def read_incar(self, filename):
        with open(filename, 'r') as f:
            fil = f.read()
            fil = fil.replace(';','\n').splitlines()
            for l in fil:
                if '=' in l:
                    (token, value) = l.split('=')
                    value = value.strip()
                    try:
                        if '.' in value:
                            value = float(value)
                        else:
                            value = int(value)
                    except:
                        pass

                    self.vasp_params[token.strip()] = value
        # self.update()


    def add_conv_kpoint(self,arg):
        if type(arg) is not str:
            sys.exit("\nadd_conv_kpoint error\n")
        if arg in self.conv_kpoint:
            print_and_log( "Warning! You already have this name in list")
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
            print_and_log( "Error! Set "+self.ise+" does not have conv_tsmear, I create new\n")
            self.conv_tsmear = []
        if arg in self.conv_tsmear:
            print_and_log( "Warning! You already have this name in list", imp = 'y')
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
            print_and_log( "Warning! Set "+self.ise+" does not have conv, I create new\n")
            self.conv = {}
        except KeyError:
            print_and_log( "Warning! Set "+self.ise+" does not have list for this key in conv, I add new\n")
            self.conv[type_of_conv] = []
        except IndexError:
            pass
        if arg in self.conv[type_of_conv]:
            print_and_log( "Warning! You already have name %s in list of conv %s. Nothing done.\n" % \
            (str(arg), str(self.conv[type_of_conv])    ) )
            return    
        self.conv[type_of_conv].append(arg)
        self.history += "Name "+arg+" was added to self.conv["+type_of_conv+"]\n"
        print_and_log( "Name "+arg+" was added to self.conv["+type_of_conv+"] of set "+self.ise+" \n")
        self.update()


    def set_compare_with(self,arg):
        if type(arg) is not str:
            raise TypeError ("\nset_compare_with error\n")
        self.compare_with += arg+" "


    def set_potential(self,znucl, arg = ''):
        # print arg
        
        if not arg:
            arg = header.PATH2POTENTIALS+'/'+invert(znucl)
            printlog('Attention!, Default potentials is chosen from ',header.PATH2POTENTIALS, 'for',invert(znucl) , imp ='Y')

        if type(arg) not in (str,):
            # sys.exit("\nset_potential error\n")
            raise RuntimeError



        if znucl in self.potdir:
            if arg == self.potdir[znucl]:
                print_and_log( "Warning! You already have the same potential for "+str(znucl)+" element\n" )
        # print type(self.potdir)
        self.potdir[znucl] = arg
        self.history += "Potential for "+str(znucl)+" was changed to "+arg+"\n"
        print_and_log( "Potential for "+str(znucl)+" was changed to "+arg+"\n" )

        # self.update()
        return



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
        print_and_log(name+"  was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        self.update()                
        #print self.history
        return

    def set_add_nbands(self,arg):
        name = "add_nbands"  
        # print(type(arg))
        if type(arg) not in [float, int, type(None) ]:
            raise TypeError
        try: 
            self.add_nbands
        except AttributeError: 
            self.add_nbands = None
        
        old = self.add_nbands    
        
        self.add_nbands = arg
        
        if old == arg:
            print_and_log("Warning! You did not change  "+name+"  in "+self.ise+" set\n")
            return
        self.history += " "+name+"  was changed from "+str(old)+" to "+str(arg) + "\n"
        print_and_log(" "+name+"  was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        return #ISTAR

    def set_ngkpt(self,arg):
        if not is_list_like(arg):
            printlog("Error! set_ngkpt type error")
        old = copy.copy(self.ngkpt)     
        self.ngkpt = copy.copy(arg)
        self.kpoints_file = True
        self.vasp_params['KSPACING'] = None
        if old == arg:
            print_and_log( "Warning! You did not change one of your parameters in new set", imp = 'y')
            return
        self.history += "ngkpt was changed from "+str(old)+" to "+str(arg) + " and KPOINTS file was swithed on\n"
        return


    def set_params_dict(self, token, arg, des = "see manual"):
        """
        Used for setting parameters for different calculators, such as VASP, Gaussian, etc

        """

        # print(token, arg)
        if token in ("ISMEAR",):
            if type(arg) not in [int, type(None), ]:
                raise TypeError
        if token in ("KSPACING",):
            # print(type(arg))
            if type(arg) not in [float, type(None), ]:
                raise TypeError

        if not hasattr(self, 'params'):
            self.params = self.vasp_params


        old = self.params.get(token)
        self.params[token] = arg
        
        if old == arg:
            print_and_log("Warning! You did not change  "+token+"  in "+self.ise+" set\n")
        else:
            self.history += " "+token+"  was changed from "+str(old)+" to "+str(arg) + "\n"
            print_and_log(token+"  was changed from "+str(old)+" to "+str(arg) +" - "+ des+" in set "+self.ise+" \n")
        
        self.update()                

        return

    def set_attrp(self, token, arg, des = "see manual"):
        """
        set any attribute.

        """
        # print (token)
        if hasattr(self, token):
            old = getattr(self, token)
            if old == arg:
                printlog("Warning! You did not change  "+token+"  in "+self.ise+" set\n")
            else:
                setattr(self, token, arg)

                self.history += " "+token+"  was changed from "+str(old)+" to "+str(arg) + "\n"
                printlog(token+"  was changed from "+str(old)+" to "+str(arg) +" - "+ des+" in set "+self.ise+" \n")
        
        else:
            setattr(self, token, arg)
            printlog("New attribute  "+token+"  added to "+self.ise+" set\n")
            self.history += " "+token+"  was added as a new attr with "+str(arg) + " value \n"

        return


    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, 
                sort_keys=True, indent=4)



    def toabinit(self, st):
        """
        Convert from VASP (add more codes in the future) to Abinit
        """

        def special_convert(vasp_param, vasp_dic, abi_dic):
            ''

            return dic

        special = {'EDIFFG', 'IBRION', 'ISIF', 'KSPACING', 'KGAMMA', 'ISMEAR', 'LDAU', 'LDAUL','LDAUU','LDAUJ',}
        
        skip = {'PREC', 'ALGO', 'POTIM'}

        VASP2Abi = {
        'ENCUT':'ecut',
        # 'ENAUG':'pawecutdg',
        'EDIFF':'toldfe', 
        'EDIFFG':'tolmxf',
        'NELM':'nstep', 
        'NSW':'ntime',
        # 'IBRION':'ionmov',
        # 'ISIF':'optcell',
        # 'PREC':['ngfft', 'boxcutmin',
        # 'ALGO':'iscf',
        # 'KSPACING':'ngkpt',
        # 'KGAMMA':'shiftk', #nshiftk
        'LREAL':None,
        'ISMEAR':'occopt',
        'SIGMA':'tsmear',
        'LPLANE':None,
        # 'POTIM':'dtion',
        'LORBIT':None,
        'ISPIN':'nsppol',
        'LDAU':'usepawu',
        'LDAUTYPE':None,
        'LDAUL':'lpawu',
        'LDAUU':'upawu',
        'LDAUJ':'jpawu',
        'LASPH':None,
        'LMAXMIX':None,
        }


        abi_dic = {}
        vp = self.vasp_params
        en = 1/header.to_eV
        fo = 1/header.Ha_Bohr_to_eV_A
        le = 1/header.to_ang
        for p in vp:
            ''
            if p in skip or p not in VASP2Abi:
                continue
            if VASP2Abi[p] is None:
                continue

            v = vp[p]
            abinam = VASP2Abi[p]

            if p == 'EDIFFG':
                aval = red_prec(v*-1*fo)
            elif p in ['ENCUT', 'EDIFF', 'ENAUG', 'SIGMA']:
                aval = red_prec(v*en )
            elif p in ['LDAU']:
                if 'T' in v:
                    aval = 1
                else:
                    aval = 0 
            elif p == 'LDAUL':
                aval = 2 # d metals


            elif p == 'ISMEAR':
                if v == 0:
                    #Gaussian
                    aval =7
                elif v == -5:
                    aval = 7 # still gauss !
            else:
                aval = vp[p]
            

            abi_dic[abinam] = aval

        for p in abi_dic:
            print(p, abi_dic[p])
        print('autoparal 1')
        print('boxcutmin 1.5') # prec normal
        print('pawecutdg', abi_dic['ecut']*2) # fine mesh
        print('ngkpt ','put here' )
        from textwrap import wrap
        import numpy as np 
        mag_str = '0 0 '+' 0 0  '.join(np.array(st.magmom).astype(str))

        print('spinat', '\n'.join(wrap(mag_str)) )


        return

























def inherit_iset(ise_new, ise_from, varset, override = False, newblockfolder = None):
    """ Create new set copying from existing and update some fields. If ise_from does not exist create new"""

    ise_new = ise_new.strip()
    ise_from = ise_from.strip()

    if ise_from not in varset:
        printlog( "\nWarning! Set "+ise_from+" does not exist. I return new empty set\n")
        return InputSet(ise_new)

    old = varset[ise_from]

    if hasattr(old, 'vasp_params') and old.vasp_params: # for compatability after renaming vasp_params to params
        old.params = old.vasp_params
        old.calculator = 'vasp'



    all_keys = vasp_electronic_keys+vasp_ionic_keys+vasp_other_keys+aims_keys+gaussian_keys

    for key in all_keys: #check if new keys was added
        if key not in old.params: 
            old.params[key] = None 

    if override:
        printlog( "\nAttention! You have chosen to override set "+ise_new+"\n")
    elif ise_new in varset:
        printlog( "\nSet "+ise_new+" already exists. I return it without changes. Be carefull not to spoil it\n")
        return varset[ise_new]           



    new = copy.deepcopy( old )
    new.ise = ise_new
    new.compare_with = ise_from+" "
    new.des = "no description for these set, see history"
    new.conv = {}

    printlog( "New set "+ise_new+" was inherited from set "+ise_from+"\n")
    new.history = old.history + "\nSet "+ise_new+" was inherited from: "+ ise_from +"\n"

    if newblockfolder:
        new.history += 'blockfolder changed from '+new.blockfolder+' to '+newblockfolder+'\n'
        new.blockfolder = newblockfolder
    
    varset[ise_new] = new
    

    return new





def make_sets_for_conv(isefrom,conv,list_of_parameters,varset):

    varset[isefrom].add_conv( isefrom, conv ); i = len(varset[isefrom].conv[conv])
    #print varset[isefrom].conv[conv]
    for param in list_of_parameters:
        newise = isefrom+conv[0:2]+str(i) ; i+=1
        if newise in varset:
            print_and_log("Set %s already in varset; continue\n" %( str(newise) ) ) 
            continue
           
        if conv == "kpoint_conv":
            for key in varset[isefrom].conv[conv]:
                if varset[key].ngkpt == param:
                    print_and_log( "Set %s already contains param %s; please check; return; \n" %( str(key), str(param) ) )
                    return
            #print newise
            s = inherit_iset(newise, isefrom, varset,newblockfolder = conv)
            s.set_ngkpt(param)
            #print s

        elif conv == "tsmear_conv":
            for key in varset[isefrom].conv[conv]:
                if varset[key].tsmear == param:
                    print_and_log( "Set %s already contains param %s; please check; return; \n" %( str(key), str(param) ) )
                    return
            s = inherit_iset(newise, isefrom, varset,newblockfolder = conv)
            s.set_tsmear(param)
        elif conv == "ecut_conv":
            #automatically set dilatmx == 1
            for key in varset[isefrom].conv[conv]:
                if varset[key].vasp_params["ENCUT"] == param:
                    print_and_log( "Set %s already contains param %s; please check; return; \n" %( str(key), str(param) ) )
                    return
            s = inherit_iset(newise, isefrom, varset,newblockfolder = conv)
            s.set_dilatmx(1.)
            s.set_ecut(param)           
        else:
            print_and_log( "Warning! Unknown type of conv; return\n")
            return


        varset[isefrom].add_conv( newise, conv )

    print_and_log( "The following sets are in varset[%s].conv %s \n"%(str(isefrom),str(varset[isefrom].conv)   ) ) 


    return



def init_default_sets(init = 0):
    """
    Pre-defined sets for Vasp
    Initialized in read_database(scratch = False, init_sets = 0) and can be updated with init_sets arg
    """
    varset = header.varset
    #print('init_default_sets():, init ', init)
    if init:
        printlog("Initializing defaults sets")



    setname = 'aks'
    if init or setname not in varset: #init only once
        s = InputSet(setname, calculator = 'vasp') #default starting set without relaxation
        s.kpoints_file = True
        s.add_nbands = 1.25
        s.params = {
            'NELM'      : 50,
            'IBRION'    : 1,
            'KGAMMA'    : ".TRUE.",
            'ENCUT'     : 441.0,
            'EDIFFG'    : 0,
            'SIGMA'     : 0.2,
            'NELMIN'    : 4,
            'ISTART'    : 0,
            'LSCALU'    : ".FALSE.",
            'MAXMIX'    : 40,
            'NSIM'      : 4,
            'ISIF'      : 2,
            'EDIFF'     : 6e-06,
            'ENAUG'     : 776.16,
            'NSW'       : 0,
            'LPLANE'    : ".TRUE.",
            'LREAL'     : "Auto",
            'ISMEAR'    : 2,
            'NPAR'      : 1,
            'ALGO'      : "Normal",
            'PREC'      : "Normal",
            'KSPACING'  : 0.235,
            }
        s.vasp_params = s.params
        s.potdir = copy.deepcopy(header.nu_dict)

        s.update()
        header.varset[setname] = copy.deepcopy(s)
    
    setname = 'static'
    if init or setname not in varset: #init only once
        s = InputSet(setname, calculator = 'vasp') #default starting set without relaxation
        s.kpoints_file = True
        s.add_nbands = 1.5
        s.params = {
            'ISTART'    : 0,
            'NELM'      : 50,
            'EDIFF'     : 1e-05,
            'NSW'       : 0,
            'EDIFFG'    : 0,
            'IBRION'    : 1,
            'ISIF'      : 2,
            'PREC'      : "Normal",
            'ALGO'      : "Normal",
            'ENCUT'     : 400,
            'ENAUG'     : 400*1.75,
            'KSPACING'  : 0.2,
            'KGAMMA'    : ".TRUE.",
            'LREAL'     : "Auto",
            'ISMEAR'    : 0,
            'SIGMA'     : 0.1,
            'LPLANE'    : ".TRUE.",
            'NPAR'      : 1,
            }
        s.vasp_params = s.params
        s.potdir = copy.deepcopy(header.nu_dict)

        s.update()
        header.varset[setname] = copy.deepcopy(s)
    



    setname = 'opt'
    if init or setname not in varset: #init only once
        # sys.exit()
        s = InputSet(setname, calculator = 'vasp') 
        s.kpoints_file = True
        s.add_nbands = 1.5
        s.params = {
            'IBRION'    : 1,
            'ENCUT'     : 150,
            'EDIFFG'    : -0.05,
            'SIGMA'     : 0.2,
            'ISIF'      : 2,
            'EDIFF'     : 1e-05,
            'NSW'       : 20,
            'ISMEAR'    : 2,
            'KSPACING'  : 0.2,
            }
        s.vasp_params = s.params
        s.potdir = copy.deepcopy(header.nu_dict)
        
        s.update()
        header.varset[setname] = copy.deepcopy(s)
    # print(header.varset[setname], setname)



    setname = 'gaus_sp'
    if init or setname not in varset: #init only once
        s = InputSet(setname, calculator = 'gaussian') #default starting set without relaxation
        # print('Init_sets', s.calculator)
        s.params = {
            'functional'    : 'B3LYP',
            'basis_set'     : '6-31G(d)',
            'job_type'      : 'SP',
            'optional'      : None,
            }
        # s.potdir = copy.deepcopy(header.nu_dict)

        # s.update()
        header.varset[setname] = copy.deepcopy(s)
    





    return

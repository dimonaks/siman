from __future__ import division, unicode_literals, absolute_import 

from header import log, print_and_log;
import copy


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
'MAXMIX'
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
'SYSTEM',
'ISTART',
'ICHARG',
'KGAMMA',
'KSPACING',
'LPLANE',
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
]

siman_keys = [
'u_ramping_region', #deprecated
'u_ramping_nstep', #number of u ramping steps
'magnetic_moments',
'afm_ordering',
]

def read_vasp_sets(varset, user_vasp_sets, override_global = False):
    """
    Read user sets and add them to project database
    Now for VASP
    ###INPUT:
        - varset (dict) - database dict with all sets of a project
        - user_vasp_sets (list) - list of user sets that describes creation of new sets based on inheritance 
        - override - allows to recreate all sets; can be usefull than you want to add some new property to all your sets - very dangerous to do!

    ###RETURN:
        - user_vasp_sets (list)

    """
    # print varset.keys()
    # varset['9'].printme()
    # print varset['9ml'].history
    # print varset['9'].history
    vasp_keys = vasp_electronic_keys+vasp_ionic_keys+vasp_other_keys
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
            # print param
            for key in param:
                
                if key in vasp_keys:
                    s.set_vaspp(key, param[key])

                elif key == 'set_potential':
                    for key2 in param[key]:
                        # print key2, 'key2'
                        s.set_potential(key2, param[key][key2])

                elif key == 'add_nbands':
                    # print param[key]

                    s.set_add_nbands(param[key])


                elif key == 'bfolder':
                    print 'New blockfolder', param[key]



                elif key in siman_keys:
                    s.set_attrp(key, param[key] )
                
                else:
                    print_and_log('Error! Uknown key: '+key)
                    raise RuntimeError




    return varset




class InputSet():
    """docstring for InputSet"""
    def __init__(self,ise):
        #super(InputSet, self).__init__()
        self.ise = ise
        self.des   = "" # description
        self.potdir = {}
        self.units = "vasp"
        self.vasp_params = {}
        self.mul_enaug = 1
        self.history = "Here is my uneasy history( :\n"    
        self.tsmear = None 
        self.tolmxf = None
        self.ngkpt  = None
        self.blockfolder = ''
        # self.kpoints_file = False
        # self.use_ngkpt = False
        #Code scpecific parameters, now only for Vasp
        for key in vasp_electronic_keys: 
            self.vasp_params[key] = None 
        for key in vasp_ionic_keys: 
            self.vasp_params[key] = None 
        for key in vasp_other_keys: 
            self.vasp_params[key] = None 

    def printme(self):
        for key in self.vasp_params:
            if self.vasp_params[key] == None: continue
            print "{:30s} = {:s} ".format("s.vasp_params['"+key+"']", str(self.vasp_params[key]) )

        print self.potdir

    def update(self):
        #deprecated, but still can be usefull
        # print_and_log('Updating set ...\n')
        c1 = 1; c2 = 1
        if self.units == "abinit":
            c1 = to_eV
            c2 = Ha_Bohr_to_eV_A
        #Update Vasp parameters
        if self.units == "vasp":
            c1 = 1
            c2 = 1
        # if self.ecut == None:
        #     self.vasp_params['ENCUT'] = None
        #     self.vasp_params['ENAUG'] = None
        # else:           
        #     self.vasp_params['ENCUT'] = self.ecut * c1* self.dilatmx * self.dilatmx
        #     self.vasp_params['ENAUG'] = self.mul_enaug * self.vasp_params['ENCUT']
        # self.vasp_params['SIGMA'] = self.tsmear * c1

        self.tsmear = self.vasp_params['SIGMA'] / c1
        self.tolmxf = - self.vasp_params['EDIFFG'] / c2
        # self.vasp_params['EDIFF'] = self.toldfe * c1
        # self.vasp_params['NELM'] = self.nstep
        # self.vasp_params['NSW'] = self.ntime
        # self.vasp_params['EDIFFG'] = -self.tolmxf * c2




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


    def set_compare_with(self,arg):
        if type(arg) is not str:
            raise TypeError ("\nset_compare_with error\n")
        self.compare_with += arg+" "


    def set_potential(self,znucl, arg):
        # print arg
        if type(arg) not in (str, unicode):
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
        log.write(name+"  was changed from "+str(old)+" to "+str(arg) + " in set "+self.ise+" \n")
        self.update()                
        #print self.history
        return

    def set_add_nbands(self,arg):
        name = "add_nbands"  
        if type(arg) not in [float, int ]:
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
        return

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
        return


    def set_vaspp(self, token, arg, des = "see manual"):
        """
        Used for setting vasp parameters.

        """


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
        else:
            self.history += " "+token+"  was changed from "+str(old)+" to "+str(arg) + "\n"
            print_and_log(token+"  was changed from "+str(old)+" to "+str(arg) +" - "+ des+" in set "+self.ise+" \n")
        
        self.update()                

        return

    def set_attrp(self, token, arg, des = "see manual"):
        """
        set any attribute.

        """
        # print token
        if hasattr(self, token):
            old = getattr(self, token)
            if old == arg:
                print_and_log("Warning! You did not change  "+token+"  in "+self.ise+" set\n")
            else:
                setattr(self, token, arg)

                self.history += " "+token+"  was changed from "+str(old)+" to "+str(arg) + "\n"
                print_and_log(token+"  was changed from "+str(old)+" to "+str(arg) +" - "+ des+" in set "+self.ise+" \n")
        
        else:
            setattr(self, token, arg)
            print_and_log("New attribute  "+token+"  added to "+self.ise+" set\n")
            self.history += " "+token+"  was added as a new attr with "+str(arg) + " value \n"

        return













def inherit_iset(ise_new,ise_from,varset,override = False, newblockfolder = None):
    """ Create new set copying from existing and update some fields. If ise_from does not exist create new"""

    ise_new = ise_new.strip()
    ise_from = ise_from.strip()

    if ise_from not in varset:
        log.write( "\nError! Set "+ise_from+" does not exist. I return new empty set\n")
        return InputSet(ise_new)

    old = varset[ise_from]

    for key in vasp_electronic_keys+vasp_ionic_keys+vasp_other_keys: #check if new keys was added
        if key not in old.vasp_params: 
            old.vasp_params[key] = None 

    if override:
        print_and_log( "\nAttention! You have chosen to override set "+ise_new+"\n")
    elif ise_new in varset:
        print_and_log( "\nSet "+ise_new+" already exists. I return it without changes. Be carefull not to spoil it\n")
        return varset[ise_new]           



    new = copy.deepcopy( old )
    new.ise = ise_new
    new.compare_with = ise_from+" "
    new.des = "no description for these set, see history"
    new.conv = {}

    print_and_log( "New set "+ise_new+" was inherited from set "+ise_from+"\n")
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

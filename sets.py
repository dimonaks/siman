from header import * 
from classes import InputSet
#1 - only volume; v  isif = 5
#2 - full relax;    vsa
#8 - no relaxation; 0
#9 - only atoms;    a


#Vasp keys
vasp_electronic_keys = [
'ALGO',
'PREC',
'LREAL',
'ENCUT',
'ENAUG',
'ISMEAR',
'SIGMA',
'EDIFF ',
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
]


def inherit_iset(ise_new,ise_from,varset,override = False, newblockfolder = ""):
    """ Create new set copying from existing and update some fields. If ise_from does not exist create new"""

    if ise_from not in varset:
        log.write( "\nError! Set "+ise_from+" does not exist. I return new empty set\n")
        return InputSet(ise_new)

    old = varset[ise_from]

    for key in vasp_electronic_keys+vasp_ionic_keys+vasp_other_keys: #check if new keys was added
        if key not in old.vasp_params: 
            old.vasp_params[key] = None 

    if override:
        print_and_log( "\nAttention! You have choosen to override set "+ise_new+"\n")
    elif ise_new in varset:
        print_and_log( "\nSet "+ise_new+" already exists. I return it without changes. Be carefull not to spoil it\n")
        return varset[ise_new]           



    new = copy.deepcopy( old )
    new.ise = ise_new
    new.compare_with = ise_from+" "
    new.des = "no description for these set, see history"
    new.conv = {}
    new.blockfolder = newblockfolder
    print_and_log( "New set "+ise_new+" was inherited from set "+ise_from+"\n")
    new.history = old.history + "\nSet "+ise_new+" was inherited from: "+ ise_from +"\n"
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



def read_archived_sets(varset):
    ise = '201'
    varset[ise] = InputSet(ise)
    s = varset[ise]
    s.units = "vasp"
    s.des = " relax shape and volume"
    # s.potdir =[]
    # s.potdir.append("potpaw_PBE_BELSU/Ti_sv")

    #Common most important parameters for majority of codes in Abinit format
    #s.description = '#Fine optimisation of volume and shape '
    s.np = 24*1 #Number of processors, cmmd
    s.ecut = 400 #eV
    s.dilatmx = 1.05
    s.tsmear = 0.2 #eV    0.005
    s.nband = 352
    #Control of relaxation
    s.ntime = 200
    s.tolmxf = 2.5e-3 # eV/A
    #Control of scf convergence
    s.nstep = 50
    s.toldff = s.tolmxf * 0.1
    s.toldfe = 1e-5 # eV

    #Code scpecific parameters, now only for Vasp
    for key in vasp_electronic_keys: 
        s.vasp_params[key] = None 
    for key in vasp_ionic_keys: 
        s.vasp_params[key] = None 
    for key in vasp_other_keys: 
        s.vasp_params[key] = None 


    s.vasp_params['ALGO'] = 'Fast'
    s.vasp_params['PREC'] = 'Normal'
    s.vasp_params['LREAL'] = 'Auto'
    s.vasp_params['ISMEAR'] = [2, " Methfessel-Paxton, the same as occopt 6 in Abinit"]
    s.vasp_params['NELMIN'] = 8 #Minimum number of scf loops
    s.vasp_params['MAXMIX'] = 20 #Number of vectors in Pulaymixer

    s.vasp_params['IBRION'] = 2 #CG algorithm for difficult relax problems, but 3 is more effective
    s.vasp_params['ISIF'] = 3 #Relax volume and shape

    s.vasp_params['KGAMMA'] = '.TRUE.'
    s.vasp_params['KSPACING'] = 0.35

    s.vasp_params['LPLANE'] = '.TRUE.'
    s.vasp_params['NPAR'] = 4 #?
    s.vasp_params['LSCALU'] = '.FALSE.'
    s.vasp_params['NSIM'] = 4
    s.update() #Calculate common parameters for VASPs



    ise = '202'
    varset[ise] = copy.deepcopy(varset['201']); s = varset[ise]
    s.vasp_params['IBRION'] = 1 #Relaxation of structures which are close to local minimum
    s.vasp_params['MAXMIX'] = 40 #Number of vectors in Pulaymixer
    s.update()

    ise = '203'
    varset[ise] = copy.deepcopy(varset['202']); s = varset[ise]
    s.mul_enaug = 1.76 #ENCUT is multiplied by this number to get ENAUG - augmented kinetic energy cutoff.
    s.kpoints_file = True
    s.vasp_params['KSPACING'] = None
    s.ngkpt = (4,4,4)
    s.toldfe = 1e-6
    s.update()



    #Make sets for kpoint convergence
    ise = '203'
    s = varset[ise]
    s.conv_kpoint = [ise]
    i = 0
    for ngkpt in (4,5,5), (5,6,6), (8,9,10): #numpy.linspace(ksi * 0.8, ksi * 1.2, num=4):
        i+=1
        varset[ise+'k'+str(i)] = copy.deepcopy(s)
        varset[ise+'k'+str(i)].ngkpt = copy.copy(ngkpt)
        s.conv_kpoint.append(ise+'k'+str(i))


    ise = '204'
    varset[ise] = copy.deepcopy(varset['203']); s = varset[ise]
    s.potdir[0] = "potpaw_PBE_MPIE/Ti_sv_new2"
    s.ngkpt = (5,6,6) #the same as 203k2
    s.update()


    #Make sets for kpoint convergence
    ise = '204'
    s = varset[ise]
    s.conv_kpoint = [ise]
    i = 0
    for ngkpt in (6,7,8), (7,8,9), (8,9,10): #numpy.linspace(ksi * 0.8, ksi * 1.2, num=4):
        i+=1
        varset[ise+'k'+str(i)] = copy.deepcopy(s)
        varset[ise+'k'+str(i)].ngkpt = copy.copy(ngkpt)
        s.conv_kpoint.append(ise+'k'+str(i))
  






def update_sets(varset):
    #print history
    if 0:
        print "This block is hand-made protocol of created sets, which can be found in database"

    
        s = inherit_iset('205','204',varset)
        s.set_KGAMMA(".FALSE.")




   







        s = inherit_iset('206','204',varset)
        s.set_potential(22,"potpaw_PBE_MPIE/Ti_sv_new")
        s.conv_kpoint = []
        ise = '206'; bs = varset[ise]
        bs.add_conv_kpoint(ise)
        i = 1 #len(bs.conv_kpoint)
        for ngkpt in (6,7,8), (7,8,9), (8,9,10):
            s = inherit_iset(ise+'k'+str(i), ise, varset)
            s.set_ngkpt( ngkpt )
            bs.add_conv_kpoint( ise+'k'+str(i) )
            i+=1

        #The next two should be compared with 203k2
        s = inherit_iset('207','203',varset)
        s.set_KGAMMA(".FALSE.")
        s.set_ngkpt( (5,6,6) )
        s.set_compare_with('203k2')

        s = inherit_iset('208','203',varset)
        s.set_LREAL(".FALSE.")
        s.set_ngkpt( (5,6,6) )
        s.set_compare_with('203k2')


        ise = '203'; bs = varset[ise]
        i = 7 #len(bs.conv_kpoint)
        for ngkpt in (5,6,6), (6,7,8), (7,8,9), (8,9,10):
            s = inherit_iset(ise+'k'+str(i), ise, varset)
            s.set_ngkpt( ngkpt )
            bs.add_conv_kpoint( ise+'k'+str(i) )
            i+=1


        #tsmear convergence 203k7, 203k8, 203k10
        ise = '203k7'; bs = varset[ise]
        bs.add_conv_tsmear(ise)
        i = 1 #len(bs.conv_kpoint)
        for tsmear in 0.1, 0.05, 0.025:
            s = inherit_iset(ise+'t'+str(i), ise, varset)
            s.set_tsmear( tsmear )
            bs.add_conv_tsmear( ise+'t'+str(i) )
            i+=1 



        s = inherit_iset('209','203k2',varset)
        s.set_PREC("Accurate")



        s = inherit_iset('215','203k2',varset)
        s.set_dilatmx(1.)
        s.set_ecut(275)
        #Ecut convergence
        ise = '215'; bs = varset[ise]
        bs.add_conv(ise,"ecut_conv")
        i = 1 #len(bs.conv_kpoint)
        for ecut in 350, 400, 550, 700:
            s = inherit_iset(ise+'e'+str(i), ise, varset)
            s.set_ecut( ecut )
            bs.add_conv( ise+'e'+str(i), "ecut_conv" )
            i+=1 

        s = inherit_iset('211','203k2',varset)
        s.set_add_nbands(1.0) 
        #nband convergence
        ise = '211'; bs = varset[ise]
        bs.add_conv(ise,"nband_conv")
        i = 1 #len(bs.conv_kpoint)
        for anb in 1.10, 1.25, 1.40, 1.55:
            s = inherit_iset(ise+'b'+str(i), ise, varset)
            s.set_add_nbands(anb)
            bs.add_conv( ise+'b'+str(i), "nband_conv" )
            i+=1   

        s = inherit_iset('800','203k2',varset)
        s.set_nmdsteps(1)

        s = inherit_iset('900','203k2',varset,"over")
        s.set_relaxation_type("ions")

        s = inherit_iset('212','203k2',varset)
        s.set_IBRION("damped")
        s.set_vaspp('POTIM',None)
        s.set_vaspp('SMASS',None)

        s = inherit_iset('210','203k2',varset)
        s.set_ALGO("Normal")

        #check npar
        s = inherit_iset('213','203k2',varset)
        s.set_NPAR(1)
        ise = '213'; bs = varset[ise]
        bs.add_conv(ise,"npar_conv")
        i = 1 #len(bs.conv_kpoint)
        for npar in 2, 3, 4, 5, 6:
            s = inherit_iset(ise+'np'+str(i), ise, varset)
            s.set_NPAR(npar)
            bs.add_conv( ise+'np'+str(i), "npar_conv" )
            i+=1 

        s = inherit_iset('214','213',varset)
        s.set_LPLANE(".FALSE.")
        ise = '214'; bs = varset[ise]
        bs.add_conv(ise,"npar_conv")
        print bs.history
        i = 1 #len(bs.conv_kpoint)
        for npar in 2, 3, 4, 5, 6:
            s = inherit_iset(ise+'np'+str(i), ise, varset)
            s.set_NPAR(npar)
            bs.add_conv( ise+'np'+str(i), "npar_conv" )
            i+=1

        ise = '210'; bs = varset[ise]
        bs.add_conv(ise,"npar_conv")
        bs.blockfolder = "other"
        i = 1 #len(bs.conv_kpoint)
        for npar in 3, 6, 1:
            s = inherit_iset(ise+'np'+str(i), ise, varset,"over")
            s.set_NPAR(npar)
            s.blockfolder = "npar_conv"
            bs.add_conv( ise+'np'+str(i), "npar_conv" )
            i+=1 



        ise = '203'; bs = varset[ise]
        i = 4 #len(bs.conv_kpoint)
        for ngkpt in (6,6,6), (8,8,8), (10,10,10):
            s = inherit_iset(ise+'k'+str(i), ise, varset)
            s.set_ngkpt( ngkpt )
            bs.add_conv_kpoint( ise+'k'+str(i) )
            i+=1 



        ise = '204'; bs = varset[ise]
        i = 4 #len(bs.conv_kpoint)
        for ngkpt in (5,6,7), (5,7,7), (6,7,7):
            s = inherit_iset(ise+'k'+str(i), ise, varset)
            s.set_ngkpt( ngkpt )
            bs.add_conv_kpoint( ise+'k'+str(i) )
            i+=1


        #
        s = inherit_iset('216','210',varset,newblockfolder = "other")
        s.set_potential(22,"potpaw_PBE_MPIE/Ti_sv_new2")
        s.set_NELMIN(4)
        s = inherit_iset('217','216',varset,newblockfolder = "other")
        s.set_KGAMMA(".FALSE.")

        paramlist = (6,7,8), (7,8,9), (8,9,10)
        make_sets_for_conv("216","kpoint_conv",paramlist,varset)
        make_sets_for_conv("217","kpoint_conv",paramlist,varset)


        s = inherit_iset('801','216',varset,"over",newblockfolder = "other")
        s.set_nmdsteps(1)
        #print varset.keys()
        paramlist = (6,7,8), (7,8,9), (8,9,10)
        s = inherit_iset('218','210',varset,newblockfolder = "other")
        s.set_NELMIN(4)
        make_sets_for_conv("218","kpoint_conv",paramlist,varset)

        s = inherit_iset('219','216',varset,newblockfolder = "other")
        s.set_tolmxf(2.5e-4)
        s.set_toldfe(1e-8)


        paramlist = (6,7,8), (7,8,9)
        make_sets_for_conv("801","kpoint_conv",paramlist,varset)

        make_sets_for_conv("801","ecut_conv",(700,),varset) #dilatmx is set to 1 automatically.

        s = inherit_iset('802','801',varset,newblockfolder = "other")
        s.set_vaspp("ISMEAR",-5,"Tetrahedron with Blochl corrections")
        make_sets_for_conv("802","ecut_conv",(700,),varset) #dilatmx is set to 1 automatically.


        make_sets_for_conv("216","tsmear_conv",(0.1, 0.05),varset)





        #constructing set for real cells
        s = inherit_iset('230','216',varset,"override",newblockfolder = "")
        s.kpoints_file = False
        s.vasp_params['KSPACING'] = 0.215
        s.set_add_nbands(1.25)
        s.set_toldfe(6e-6) #My big cell contains at least 6*8 = 48 atoms, that is why I can reduce this relative to toldfe for 8 atom cell


        s = inherit_iset('830','230',varset,newblockfolder = "")
        s.vasp_params['KSPACING'] = 0.215 #1 6 6 for C1
        s.set_nmdsteps(1)

        s = inherit_iset('929','230',varset,newblockfolder = "")
        s.vasp_params['KSPACING'] = 0.215
        s.set_relaxation_type("ions")
        s.set_tolmxf(2.5e-2)
        s.set_toldfe(1e-4)
        s = inherit_iset('8301','830',varset,newblockfolder = "")
        s.vasp_params['KSPACING'] = 0.2 #2 5 6 for T111b_ml

        s = inherit_iset('8302','830',varset,newblockfolder = "")
        s.set_vaspp('KSPACING',0.25) #'2 4 5 [0.15 0.22 0.21]' for t111b_r; '1 4 5 [0.22 0.22 0.21]' for t111g
        s = inherit_iset('9292','929',varset,newblockfolder = "")
        s.set_vaspp('KSPACING',0.25) 

        s = inherit_iset('2302','230',varset,newblockfolder = "")
        s.set_vaspp('KSPACING',0.25) 

        s = inherit_iset('9291','929',varset,newblockfolder = "")
        s.set_vaspp('KSPACING',0.2)
        s.set_NPAR(1)
        s = inherit_iset('9282','9292',varset,newblockfolder = "")
        s.set_vaspp('ISYM',0)
        s = inherit_iset('9292f','9292',varset,newblockfolder = "") #fine relaxation
        s.set_toldfe(6e-6)
        s.set_tolmxf(2.5e-3)

        s = inherit_iset('83','830',varset,newblockfolder = "")
        s.vasp_params['KSPACING'] = 0.235 # 2 3 6 for t21b_ml
        s.set_NPAR(1)

        s = inherit_iset('93','83',varset,newblockfolder = "")
        s.set_relaxation_type("ions") #200 md_steps automatically
        s.set_tolmxf(2.5e-2)
        s.set_toldfe(1e-4)
    
        s = inherit_iset('935','93',varset,newblockfolder = "") # for comparison with Lane2011
        s.set_potential(22,"potpaw_PBE_MPIE/Ti")
        s.set_ecut("default")

        s = inherit_iset('835','83',varset,newblockfolder = "")
        s.set_potential(22,"potpaw_PBE_MPIE/Ti")
        s.set_ecut("default")

        #sets for H cell tests
        s = inherit_iset('93','93',varset,newblockfolder = "")
        s.set_potential(6,"potpaw_PBE_MPIE/C")


        paramlist = (8,8,8), (10,10,10)
        #varset['93'].conv["kpoint_conv"] = []
        #del varset['93kp1']
        make_sets_for_conv("93","kpoint_conv",paramlist,varset)

        
        make_sets_for_conv("93","ecut_conv",(550, 700),varset)
        s = inherit_iset('83','83',varset,newblockfolder = "")
        s.set_potential(6,"potpaw_PBE_MPIE/C")
        paramlist = (8,8,8), (10,10,10)        
        make_sets_for_conv("83","kpoint_conv",paramlist,varset)
        make_sets_for_conv("83","ecut_conv",(550, 700),varset)

        s = inherit_iset('83is1','83',varset,newblockfolder = "") #ib1
        s.set_vaspp("ISMEAR",-5,"Tetrahedron with Blochl corrections")
        s = inherit_iset('83kp8','83',varset,newblockfolder = "")
        s.vasp_params['KSPACING'] = 0.225 # 2 4 6 for csl71b_r

        s = inherit_iset('93kp7','93',varset,newblockfolder = "") #93kp7 = 929 but with carbon for c1
        s.vasp_params['KSPACING'] = 0.215 # the same as 929

        s = inherit_iset('93kp6','93',varset,newblockfolder = "") #93kp6 = '9291' but with carbon
        s.vasp_params['KSPACING'] = 0.2 # the same as '9291'

        s = inherit_iset('93kp9','93',varset,newblockfolder = "") #93kp9 = '9292' '8302' but with carbon for t111
        s.vasp_params['KSPACING'] = 0.25 # the same as '9292

        #changing istart from auto to 0; 16.01.2014
        varset['93'].vasp_params['ISTART'] = 0
        varset['93kp9'].vasp_params['ISTART'] = 0
        varset['93kp7'].vasp_params['ISTART'] = 0
        #print varset['93kp7'].vasp_params['ISTART']
        varset['93'].set_potential(8,"potpaw_PBE_MPIE/O") 

        for key in varset: #935, 93kp9, 93kp6, 93kp7, 93kp2, 93kp1, 93, 93ec1, 93ec2
            if "93" in key:
                print key
                varset[key].set_potential(8,"potpaw_PBE_MPIE/O") 
                varset[key].vasp_params['ISTART'] = 0
        for key in "83is1",:
            varset[key].set_potential(8,"potpaw_PBE_MPIE/O") 
            varset[key].vasp_params['ISTART'] = 0
        
        s = inherit_iset('93kp9acc','93kp9',varset,newblockfolder = "") 
        s.set_PREC("Accurate")

        s = inherit_iset('83hsp','83',varset,newblockfolder = "")#high symprec
        s.set_vaspp('SYMPREC',0.05)
        s.set_vaspp('KSPACING',0.5)
        s = inherit_iset('83hsp_nosym','83hsp',varset,True,newblockfolder = "")
        s.set_vaspp('ISYM',0)

        varset['93'].set_nmdsteps(25)
        s = inherit_iset('93fe1','93',varset,newblockfolder = "")
        s.set_toldfe(1e-6) 

        s = inherit_iset('93acc','93',varset,newblockfolder = "")
        s.set_PREC("Accurate")
        
        s = inherit_iset('93PT1','93',varset,newblockfolder = "")
        s.set_vaspp('POTIM',0.2)

        s = inherit_iset('93CG','93',varset,newblockfolder = "")
        s.set_IBRION("CG")

        s = inherit_iset('93fe1MD','93',varset,newblockfolder = "")
        s.set_toldfe(1e-6)
        s.set_IBRION("MD")
        s.set_vaspp('POTIM',1.0)

        s = inherit_iset('93fe1D','93',varset,newblockfolder = "")
        s.set_toldfe(1e-6)
        s.set_IBRION("damped")
        s.set_vaspp('POTIM',0.8)

        s = inherit_iset('94','93fe1',varset,newblockfolder = "")
        s.set_vaspp('KSPACING',0.15) # 2 5 9 for t21g
        s.set_PREC("Accurate")

        s = varset['93kp9']
        s.set_vaspp('ISYM',0)
        s.set_nmdsteps(25)

        s = inherit_iset('93ns','93',varset,newblockfolder = "") #no sym
        s.set_vaspp('ISYM',0)

        s = inherit_iset('93hd','93ns',varset,newblockfolder = "") #hard potentials, high ecut, no sym
        s.set_potential(22,"potpaw_PBE_MPIE/Ti_sv_h")
        s.set_potential(6, "potpaw_PBE_MPIE/C_h")
        s.set_potential(8, "potpaw_PBE_MPIE/O_h")
        s.set_ecut(700)
    
        s = varset['93kp7']
        s.set_vaspp('ISYM',0)
        s.set_nmdsteps(25)
        s.set_toldfe(2e-5)

        s = inherit_iset('83kp9','83',varset,newblockfolder = "")
        s.vasp_params['KSPACING'] = 0.25 # the same as 8302; just new name

        varset['83'].set_vaspp('ISYM',0)


        #set for atoms
            
        s = inherit_iset('8at','83',varset) #atom
        s.des = "set for one atom in big box"
        s.kpoints_file = True
        s.vasp_params['KSPACING'] = None
        s.set_KGAMMA(".FALSE.")
        s.ngkpt = (1,1,1)
        s.set_LREAL(".FALSE.")
        s.set_vaspp('ISMEAR',0)
        s.set_tsmear(0.1)
        s.set_potential(8,"potpaw_PBE_MPIE/O") 

        s = inherit_iset('9dim','8at',varset) #dimer
        s.set_IBRION("MD")
        s.set_relaxation_type('ions')
        s.set_nmdsteps(10)
        s.set_vaspp('POTIM',1)
        s.set_vaspp("SMASS", -2)
        
        s = inherit_iset('83lf','83',varset) #atom
        s.set_LREAL(".FALSE.")

        s = inherit_iset('83is1lf','83is1',varset) #atom
        s.set_LREAL(".FALSE.")
        s.set_add_nbands(1.50) #just to check; the problem was not in bands but in parallelization

        s = inherit_iset('83kp7','93kp7',varset,newblockfolder = "")
        s.set_nmdsteps(1)
        s.set_toldfe(6e-6)
        s = varset['83kp9'] #was inherited from 83 before 83 had ISYM = 0 
        s.set_vaspp('ISYM',0) #as a result all 83* sets have small toldfe = 6e-06 eV, while 93* sets have toldfe = 1e-04 eV

        s = inherit_iset('93kp9p4','93kp9',varset,newblockfolder = "")
        s.set_potential(22,"potpaw_PBE_MPIE/Ti") #For comparison with Ghazisaeidi2014 (Trinkle) - 4 electron potential



        s = inherit_iset('93f','93',varset,newblockfolder = "") #fine relaxation
        s.set_toldfe(6e-6)
        s.set_tolmxf(2.5e-3)







        s = inherit_iset('83dos','83',varset, override = True) # The same k-points as 83
        s.set_potential(8,"potpaw_PBE_MPIE/O") 
        s.set_vaspp('LORBIT',12)
        s.set_vaspp('ISMEAR',-5) #Tetrahedron
        s.set_vaspp('LAECHG', '.TRUE.')
        s.set_vaspp('EMIN', -2) #was changed hs443C.r with old
        s.set_vaspp('EMAX',  14)
        s.set_vaspp('NEDOS', 2000) #was changed
        s.set_toldfe(1e-4) #was changed


        s = inherit_iset('83acdos','83dos',varset, override = True) # accurate and to test fft grid for bader
        s.set_PREC("Accurate")

        s = inherit_iset('dos','83dos',varset, override = True)
        s.set_vaspp('KSPACING', 0.14)# gives 5 5 4 for all hs443 cells

        s = inherit_iset('dosh','83dos',varset, override = True)#dos high
        s.set_vaspp('KSPACING', 0.1)# 

        s = inherit_iset('dosf','dos',varset, override = True) #fine energy convergence 
        s.set_toldfe(6e-6) #was changed


        s = inherit_iset('dosm','83dos',varset, override = True, newblockfolder = "dosm")#dos middle
        s.set_vaspp('KSPACING', 0.18)# 

        s = inherit_iset('lobster','83',varset,override = True,  newblockfolder = "lobster") # The same k-points as 83
        s.set_potential(8,"potpaw_PBE_MPIE/O") 
        # s.set_vaspp('LSORBIT', 'TRUE')# 
        s.kpoints_file = 'kfiles/t111sg.2'
        # s.set_add_nbands(1.5) 

        s = inherit_iset('23','93',varset)
        s.set_relaxation_type("full")

        s = inherit_iset('93kp7','93kp7',varset)
        s.set_potential(5,"potpaw_PBE_MPIE/B") 



        s = inherit_iset('8at_sv','8at',varset) #hard
        s.set_potential(8,"potpaw_PBE_MPIE/O_sv") #ecut 1066 eV 
        s.set_ecut(1000)



        s = inherit_iset('9dim_spin','9dim',varset)
        s.nstep = 100;s.update()
        s.set_vaspp('ISPIN', 2)
        s.set_tsmear(0.01)

        s = inherit_iset('8at_spin','8at',varset) #hard
        s.nstep = 100;s.update()
        s.set_vaspp('ISPIN', 2)
        s.set_tsmear(0.01)
        s.set_add_nbands(2.) #was not enough bands 

        s = inherit_iset('83nb15','83',varset) #hard
        s.set_add_nbands(1.5)

    # s = inherit_iset('85','83',varset) #hard
    # s.set_tsmear(0.05)
    # s.set_add_nbands(1.1)

        s = inherit_iset('95','93',varset) #hard
        s.set_tsmear(0.05)
        s.set_add_nbands(1.1)
        s.set_nmdsteps(10)

        """Start of Ti-H changes"""

        s = varset['83']
        s.kpoints_file = 1
        s.ngkpt = ()
        s.set_vaspp('ISYM',None)
        s.set_potential(1,"potpaw_PBE_MPIE/H") 

        s = inherit_iset('83ns','83',varset) #checking the influence of total drift
        s.set_vaspp('ISYM',0)



        s = varset['93']
        s.kpoints_file = 1
        s.ngkpt = ()
        s.set_vaspp('ISYM',None)
        s.set_potential(1,"potpaw_PBE_MPIE/H") 

        s = inherit_iset('23','23',varset)
        s.kpoints_file = 1
        s.ngkpt = ()
        s.set_vaspp('ISYM',None)
        s.set_vaspp('PSTRESS',-5) #kbar
        s.set_potential(1,"potpaw_PBE_MPIE/H") 

        s = inherit_iset('23f','23',varset) #fine relaxation
        s.set_toldfe(6e-6)
        s.set_tolmxf(2.5e-3)

        s = inherit_iset('83p','83',varset, override = True) #phonons, high k-points
        s.set_vaspp('KSPACING', 0.1)


        s = inherit_iset('83p1','83',varset, override = True) #phonons, high k-points
        s.set_vaspp('KSPACING', 0.08)

        s = inherit_iset('83p15','83',varset, override = True) #phonons, high k-points
        s.set_vaspp('KSPACING', 0.15)
        
        s = inherit_iset('23pfns','23f',varset, override = True)
        s.set_vaspp('KSPACING', 0.1)
        s.set_toldfe(6e-6)
        s.set_tolmxf(5e-3) # two times larger than previous. 2.5e-3 was two difficult to relax
        s.set_vaspp('ISYM',0)

        s = inherit_iset('83pns','83',varset, override = True) #phonons, high k-points
        s.set_vaspp('KSPACING', 0.1)
        s.set_vaspp('ISYM',0)
        s.set_toldfe(1e-6)


        s = inherit_iset('23ac','23pfns',varset, override = 0) #phonons, high k-points
        s.set_vaspp('KSPACING', 0.145)
        s.set_toldfe(1e-6)
        s.set_vaspp('LREAL', False)
        s.set_vaspp('PREC', "Accurate")
        
        s = inherit_iset('23ac2','23ac',varset, override = 1)
        s.set_vaspp('ADDGRID', True)

        s = inherit_iset('83ac2','23ac2',varset, override = 1)
        # s.set_vaspp('NSW', 0)
        s.set_nmdsteps(0)
        s.set_vaspp('IBRION', -1)
        s.set_vaspp('PSTRESS',None) #kbar

        s = inherit_iset('23ac3','23ac2',varset, override = 1)
        s.set_toldfe(1e-8)
        s.set_tolmxf(2e-3) #


# for phonons you should use:
# PREC = Accurate avoid wrap around errors
# LREAL = .FALSE. reciprocal space projection technique
# EDIFF = 1E-6 high accuracy required

    #May be useful
#For very accurate energy calculations:
#s.vasp_params['ISMEAR'] = [-5, " Tetrahedron with Blochl corrections"]
#However this approach needs at least 3 kpoints and produce errors in stress tensor and forces

        #





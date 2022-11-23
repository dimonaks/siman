from __future__ import division, unicode_literals, absolute_import 


import ctypes
from tabulate import tabulate
#from ctypes import *
from ctypes import cdll
from ctypes import c_float, byref

import numpy as np
import traceback, os, sys, datetime, glob, copy


from siman import header
from siman.header import print_and_log, printlog, geo_folder, runBash
from siman.classes import CalculationVasp, Structure
from siman.set_functions import InputSet
from siman.functions import  return_atoms_to_cell, element_name_inv
from siman.inout import write_xyz
from siman.geo import local_surrounding, local_surrounding2, xcart2xred, xred2xcart
from siman.small_functions import block_print, enable_print

libfile = os.path.dirname(__file__)+'/libfindpores.so'
if os.path.exists(libfile):
    lib = cdll.LoadLibrary(libfile)
else:
    lib = None
    printlog('Warning! libfindpores.so is not available, search of voids is not possible. Please compile it inside siman folder using makefile')


def create_c_array(pylist, ctype):
    if ctype == float:
        c_array = (ctypes.c_float * len(pylist))(*pylist)
    return c_array


def find_pores(st_in, r_matrix=1.4, r_impurity = 0.6, step_dec = 0.05, fine = 0.3, prec = 0.1, calctype = 'central', gbpos = 0,
    find_close_to = (), check_pore_vol = 0):
    """
    st_in - input Structure() object
    r_impurity (A)- all pores smaller than this radius will be found
    r_matrix (A) - radius of matrix atoms disregarding to their type
    step_dec - scanning step of the cell in Angstroms
    fine - allows to change density of local points; local_step = step_dec/fine
    prec - precicion of pore center determination
    check_pore_vol - allows to estimate volume of pores; has problems for big cells


    'find_close_to' - works in the cases of gb and grain_vol; allows to ignore all and find pore close to provided three reduced coordinates
    return - instance of Structure() class with coordinates of pores. Number and type of included pores depend on the argument of 'calctype'.
    """
  
    xred   = st_in.xred
    natom  = len(xred)
    rprimd = st_in.rprimd
    name = st_in.name
    #print xred

    """Additional values"""
    # check_pore_vol = 1
    #if calctype in ("pore_vol","gb","grain_vol","all_local" ): check_pore_vol = 1 #something wrong with this function, especially for big cells

    """----Conversions of types for C++"""
    r1 = create_c_array(rprimd[0], float)
    r2 = create_c_array(rprimd[1], float)
    r3 = create_c_array(rprimd[2], float)

    xred1 = (c_float * len(xred))(*[x[0] for x in xred])
    xred2 = (c_float * len(xred))(*[x[1] for x in xred])
    xred3 = (c_float * len(xred))(*[x[2] for x in xred])

    max_npores = 10000;
    ntot = ctypes.c_int32(); npores = ctypes.c_int32()
    l_pxred1 = (c_float * max_npores)(0) #make static arrays fol local points
    l_pxred2 = (c_float * max_npores)(0)
    l_pxred3 = (c_float * max_npores)(0)
    l_npores = (ctypes.c_int32 * max_npores)(0)

    pxred1 = (c_float * max_npores)(0) #make static arrays natoms + npore
    pxred2 = (c_float * max_npores)(0)
    pxred3 = (c_float * max_npores)(0)

    """----Run C++ function"""
    print_and_log("Starting C++ function lib.findpores()...\n")
    # print(r_matrix, r_impurity, step_dec, fine, prec)
    lib.findpores ( check_pore_vol, \
                    max_npores, \
                    byref(ntot),   l_pxred1, l_pxred2, l_pxred3, l_npores, \
                    byref(npores), pxred1,   pxred2,   pxred3,  \
                    natom,         xred1,    xred2,    xred3, \
                    c_float(r_matrix), c_float(r_impurity), c_float(step_dec), c_float(fine), c_float(prec), \
                    r1, r2, r3 )

    print_and_log( "ntot is ", ntot.value)
    print_and_log( "l_npores[0] is ",l_npores[0])

    v = np.zeros((3))
    l_pxred = []
    shift1 = 0; shift2 = 0
    for i_por in range(npores.value):
        l_pxred.append( [] )
        shift2+=l_npores[i_por]    

        for i in range(shift1, shift2):
            v[0] = l_pxred1[i];     v[1] = l_pxred2[i];     v[2] = l_pxred3[i]

            l_pxred[i_por].append( v.copy() )

        shift1 = shift2



    if shift2 != ntot.value: 
        print_and_log( "Error! final shift2 not equal to ntot")

    #print l_pxred[0]

    pxred = [] # only coordinates of pores
    #print pxred1[natom]
    for i in range(npores.value):
        v[0] = pxred1[i+natom]; v[1]= pxred2[i+natom]; v[2] = pxred3[i+natom] #with shift, because first natom elements are coordinates of atoms
        pxred.append( v.copy() )

    #print pxred
    """----End of C++; result is two lists: lpxred - local geometry of all pores, pxred - coordinates of all pores"""


    """ Analyse of pores """
    # st_result = Structure()
    st_result = st_in.new()


    st_result.rprimd = rprimd
 
    targetp = np.array((0.,0.,0.))
    if find_close_to: 
        targetp = np.asarray(find_close_to) #targer point
        print_and_log( "Target point is                        ",targetp)
 
    a = step_dec/fine #the side of little cube formed by the mesh which is used to find spheres inside the pore.
    aaa = a*a*a


    #find most central pore
    if calctype == 'central': #return coordinates of the most central pore
        st_result.name = "central_pore_from "+name 
        center = np.array((0.5,0.5,0.5))#center of cell
        d_min = 100
        for x in pxred:
            d = np.linalg.norm(x - center)
            #print x, x-center, d
            if d < d_min and x[0] <= 0.5 and x[1] <= 0.5 and x[2] <= 0.5:
                d_min = d
                x_min = x
        print_and_log( "The closest pore to the center has coordinates",x_min)
        st_result.xred.append( x_min )


    elif calctype == 'gb': #add impurity at gb
        st_result.name = "gb_pore_from "+name
        d_min = 100; #d2_min = 100
        dt_min =100
        i_min = 0; x_min = np.zeros((3))
        for i, x in enumerate(pxred):
            #print "l_npores ",l_npores[i]
            d = abs(x[0] - gbpos/rprimd[0][0]) #
            #print x[0], d
            if find_close_to:   closer = (np.linalg.norm(targetp - x) < dt_min)
            else:               closer = ( d < d_min ) # and x[1]>0.3 and x[2]>0.3:

            if closer:
                x_pre = x_min
                i_pre = i_min
                d_min = d
                dt_min = np.linalg.norm(targetp - x)
                x_min = x
                i_min = i

            #find and add impurity in bulk
            #d2 = abs( x[0] - (gbpos/rprimd[0][0] - 0.25) )
            #if d2 < d2_min: 
            #    d2_min = d2
            #    x2_min = x
            #    i2_min = i

        #print "rprimd[0][0]", rprimd[0][0]
        print_and_log( "Position of boundary is ",gbpos/rprimd[0][0])
     
        #x_min[0] = gbpos/rprimd[0][0]
        if find_close_to: print_and_log( "The closest pore to the target point is [ %.2f  %.2f  %.2f ]"%(x_min[0], x_min[1], x_min[2]))
        else: print_and_log( "The closest pore to the gb has coordinates",x_min)
        st_result.xred.append( x_min )
        #st_result.xred.append( x_pre )           
        #Calculate volume of the pore using local balls:
        print_and_log( "The number of pore is ",i_min," ; It has ",l_npores[i_min], "local balls")
        print_and_log( "Volume of pore is ", l_npores[i_min] * a*a*a, " A^3")
        #st_result.xred.extend( l_pxred[i_min] )
        #st_result.xred.extend( l_pxred[i_pre] )   

        #print "The closest pore to the center of bulk has coordinates",x2_min
        #st_result.xred.append( x2_min )           
        #Calculate volume of the pore using local balls:
        #print "The number of bulk pore is ",i2_min," ; It has ",l_npores[i2_min], "local balls"
        #print "Volume of pore is ", l_npores[i2_min] * a*a*a, " A^3";
        #st_result.xred.extend( l_pxred[i2_min] )  
    
    elif calctype == 'grain_vol': #add impurity to volume of grain
        st_result.name = "grain_volume_pore_from "+name
        d2_min = 100
        dt_min = 100
        i_min = 0; x_min = np.zeros((3))
        for i, x in enumerate(pxred):
            #find and add impurity to the  bulk
            d2 = abs( x[0] - (gbpos/rprimd[0][0] - 0.25) )
            
            if find_close_to:   closer = (np.linalg.norm(targetp - x) < dt_min)
            else:               closer = ( d2 < d2_min ) # and x[1]>0.3 and x[2]>0.3:

            if closer:
                dt_min = np.linalg.norm(targetp - x)
                d2_min = d2
                x2_min = x
                i2_min = i

        if find_close_to: print_and_log( "The closest pore to the target point is [ %.2f  %.2f  %.2f ]"%(x2_min[0], x2_min[1], x2_min[2]))
        else:             print_and_log( "The closest pore to the center of bulk has coordinates",x2_min)
        st_result.xred.append( x2_min )           
        #Calculate volume of the pore using local balls:
        print_and_log( "The number of bulk pore is ",i2_min," ; It has ",l_npores[i2_min], "local balls")
        print_and_log( "Volume of pore is ", l_npores[i2_min] * a*a*a, " A^3")
        st_result.xred.extend( l_pxred[i2_min] )  

    elif calctype == 'all_local':
        st_result.name = "all_local_points_from "+name
        v_max = 0
        i_max = 0
        for i in range(npores.value):
            v_pore = l_npores[i] * aaa
            print_and_log(  "Volume of pore is ", l_npores[i] * aaa, " A^3")
            if v_pore > v_max: v_max = v_pore; i_max = i
        print_and_log( "Pore number ", i_max,"has the largest volume ", v_max," A^3")
        # st_result.xred = l_pxred[i_max] # here coordinates of all local points to show geometry of pore with largerst volume
        st_result.xred = [x for group in l_pxred for x in group ] # all pores


    elif calctype == 'all_pores':
        st_result.name = "all_local_pores_from "+name
        st_result.xred = pxred


    st_result.rprimd = rprimd
    st_result.xred2xcart()
    st_result.typat = [1 for x in st_result.xred]
    st_result.ntypat = 1
    st_result.natom = len(st_result.typat)
    st_result.znucl = [200]
    st_ntypat = 1

    return st_result




def add_impurity(it_new, impurity_type = None, addtype = 'central', calc = [], r_pore = 0.5,
    it_to = '', ise_to = '', verlist_to = [], copy_geo_from = "", find_close_to = (),add_to_version = 0,
    write_geo = True, only_version = None, fine = 4, put_exactly_to = None, check_pore_vol = 0, replace_atom = None, override = False):

    """
    Add impurities in pores.

    Input:
    it_new - name of new structure with impurity
    
    impurity_type - name of impurity from Mendeley table, for example 'C'
    
    addtype - type of adding: ['central',]; 'central' means that impurity 
    will be placed as close to the geometrical center of cell as possible.
    
    it_to , ise_to , verlist_to - completed calculations in which impurity 
    will be added
    
    if 'verlist_to' is empty, function will try to find geometry files in 'geo_folder + struct_des[it_to].sfolder' folder;
    even if 'it_to' is empty it will try to find files in 'geo_folder + struct_des[it_new].sfolder+'/from' ' folder.
    'ise_to' also can be empty

    if 'copy_geo_from' is not empty, then programm copy all files from folder 'copy_geo_from' to 
    folder 'geo_folder + struct_des[it_to].sfolder+"/"+it_to' or  'geo_folder + struct_des[it_new].sfolder+"/from" '  

    'find_close_to' is tuple of three reduced coordinates of point close to which you want to find impurity. If empty - ignored; 

    'add_to_version' is integer number added to each 'verlist_to' number to produce ver_new.
    
    'only_version' - if == [v,], then instertion will be provided only for v. If None insertion will be made in all found versions

    If you want to add impurity to relaxed structure ...

    'fine' - integer number; allows to reduce number of small steps for defining center


    Possible addtype's:
    'central' - add one atom to the pore which is most close to  the center of the cell but with reduced coordinates less than 0.5 0.5 0.5
    'all_pore'  - add atoms in every found pore
    'all_local' - add atoms to every local point which allows to visualise topology of pores.
    'gb' - uses self.gbpos and places atom close to this value assuming that it will be at gb
    'grain_vol' - uses self.gbpos and assuming that cell contains two gb and two equal grains, places atom close to the centre of grain; y and z can be arbiratry


    put_exactly_to - will add impurity to this point 
    find_close_to - will try to find closest void and insert pore here.

    check_pore_vol - allows to estimate volume of pores; has problems for big cells

    replace_atom - if not None, than the specified atom is substituted


    Side effects: creates new geometry folder with input structures; 

    """

    struct_des = header.struct_des

    def test_adding_of_impurities(added, init, v):
        """
        Can be used only inside add_impurity()
        Replicates the structure and find again pores
        """
        global natoms_v1
        if added == None: return
        if v == 1: #TEST
            
            natoms_v1 = len(added.init.xcart) # for test 
            st_rep_after  = added.init.replic( (1,2,1) )

            rep = copy.deepcopy(init)

            rep.init = rep.init.replic( (1,2,1) );   
            #print rep
            rep = add(znucl, "", rep, write_geo = False)
            #print rep
            #print "xcart of replic after adding ", st_rep_after.xcart
            #print "xcart of adding to    replic ", rep.init.xcart
            if len(st_rep_after.xcart) != len(rep.init.xcart): raise RuntimeError
            p = 0
            #for x2 in st_rep_after.xcart:
            #    print x2
            for x in rep.init.xcart:
                a = any(  ( np.around(x2, p) == np.around(x, p) ).all() for x2 in st_rep_after.xcart   )
                #b = any(  ( np.ceil(x2, p)   == np.ceil(x, p)  ).all()  for x2 in st_rep_after.xcart   )
                #c = any(  ( np.floor(x2, p)  == np.floor(x, p) ).all()  for x2 in st_rep_after.xcart   )
                #print a, b, c
                #np.concatenate(a, b, c):
                if not a:
                    print_and_log( "Error! Can't find ", np.around(x,3), "in replic  ")
                    raise RuntimeError

            #assert all([ all( np.around(v1, 8) == np.around(v2, 8) ) for (v1, v2) in zip(st_rep_after.xcart, rep.init.xcart) ])
            print_and_log( "add_impurity: test succesfully done")

        if natoms_v1 != len(added.init.xcart): print_and_log("You have different number of pores in different versions\n");  raise RuntimeError
        return
    


    def add(znucl, xyzpath = "", new = None, write_geo = True, put_exactly_to = None):
        "if put_exactly_to is True, then atom just added and nothing are searched"


        if write_geo and os.path.exists(new.path["input_geo"]) and not override:
            print_and_log("add: File '"+new.path["input_geo"]+"' already exists; continue\n", imp = 'Y');
            return new

        #new.init = return_atoms_to_cell(new.init)
        if replace_atom:
            #atom substitution
            if znucl not in new.init.znucl:
                new.init.znucl.append(znucl)
                new.init.ntypat+=1
                new.init.typat[replace_atom] = new.init.ntypat
            else:
                ind = new.init.znucl.index(znucl)
                new.init.typat[replace_atom] = ind + 1
            new.init.nznucl = []
            for typ in range(1, new.init.ntypat+1):
                new.init.nznucl.append(new.init.typat.count(typ) )



        else:
            new_before = copy.deepcopy(new)
            
            # new.init.xcart[-2][0]-=0.9 #was made once manually for c1gCOi10.1
            # new.init.xcart[-2][2]+=0.2
            # new.init.xred = xcart2xred(new.init.xcart, new.init.rprimd)
            write_xyz(new.init)
            #step = 0.042
            step = 0.06
            #r_pore = 0.56
            #fine = 0.3 # for visualisation of pores
            #fine = 4   #controls small steps; the steps are smaller for larger numbers
            #r_pore = 0.54
            prec = 0.004 # precision of center Angs
            if new.hex_a == None:
                r_mat = 1.48 -step
            else:
                r_mat = new.hex_a / 2 - step

            if put_exactly_to:
                pores_xred = [np.array(put_exactly_to),]
                print_and_log( 'Inmpurity just put in ', pores_xred, imp = 'Y')
            else:
                pores = find_pores(new.init, r_mat, r_pore, step, fine, prec,  addtype, new.gbpos, find_close_to, check_pore_vol) #octahedral
                pores_xred = pores.xred
            


            npores = len(pores_xred)
            
            st = new.init

            #delete last oxygen; was made once manually for c1gCOi10.1
            # st.natom-=1
            # del st.xred[-1]
            # del st.typat[-1]




            st.natom += npores
            st.xred.extend( pores_xred )

            if znucl in st.znucl:
                print_and_log( "znucl of added impurity is already in cell")
                ind = st.znucl.index(znucl)
                typat = ind+1
                st.nznucl[ind]+=npores
            else:
                st.ntypat +=1
                typat = st.ntypat
                st.znucl.append( znucl )
                st.nznucl.append( npores )



            for i in range( npores  ):
                st.typat.append( typat )



            st.xred2xcart()

            new.init = st

            #print "Add impurity: len(xred ", len(new.init.xred)
            #print "natom", new.init.natom


            #For automatisation of fit
            try: 
                #new.build
                if new.build.nadded == None:      new.build.nadded=npores
                else: new.build.nadded+=npores
                if new.build.listadded == [None]: new.build.listadded = range(new.natom - npores, new.natom) #list of atoms which were added
                else: new.build.listadded.extend( range(new.natom - npores, new.natom) )
                #print "Warning!!! Information about added impurities rewritten"
            except AttributeError: 
                pass

            #new.init.znucl = new.znucl
            #new.init.typat = new.typat
            
            #write_xyz(replic(new.init, (2,1,2))  , xyzpath)

            #test_adding_of_impurities(new, new_before, v)

            print_and_log("Impurity with Z="+str(znucl)+" has been added to the found pore in "+new.name+"\n\n")
            




        if write_geo:
            write_xyz(new.init , xyzpath)
            new.write_geometry("init",new.des, override = override)

        print_and_log( "\n")


        return new

    """0.Begin----------------------------------------------------------------------------"""

    znucl = element_name_inv(impurity_type)


    if impurity_type != 'octa' and impurity_type not in it_new:
        print_and_log("add_impurity: Your name 'it_new' is incorrect!\n\n")
        raise RuntimeError
    #del header.history[-2]
    #
    
    #hstring = ("add_impurity('%s', '%s', '%s', calc, %.3f, '%s', '%s', %s, '%s')  #at %s" %
    #    (it_new, impurity_type, addtype, r_pore,
    #        it_to, ise_to, verlist_to, copy_geo_from,
    #     datetime.date.today() ) )
    
    hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
    if hstring != header.history[-1]: header.history.append( hstring  )


    #geo_exists = 
    





    """1. The case of insertion to existing calculations--------------------------------------------------"""

    if verlist_to:         

        for v in verlist_to:
            if only_version and v not in only_version: continue # only_version = None works for all versions 
            id = (it_to, ise_to, v)
            new = copy.deepcopy(calc[id])

            new.init = new.end #replace init structure by the end structure

            new.version = v+add_to_version
            new.name = it_new#+'.'+id[1]+'.'+str(id[2])
            new.des = 'Obtained from '+str(id)+' by adding '+impurity_type+' impurity '
            path_new_geo = struct_des[it_new].sfolder+"/"+it_new+"/"+it_new+'.imp.'+addtype+'.'+str(new.version)+'.'+'geo'
            new.init.name = it_new+".init."+str(new.version)
            xyzpath = struct_des[it_new].sfolder+"/"+it_new
            
            new.path["input_geo"] = geo_folder+path_new_geo
            
            print_and_log("File '"+new.path["input_geo"] +"' with impurity will be created\n");
            #new.init.name = 'test_before_add_impurity'

            new = add(znucl, xyzpath, new, write_geo, put_exactly_to = put_exactly_to)

            



        """2. The case of insertion to geo files------------------------------------------------------------"""
    else:     

        """ Please rewrite using new functions """

        print_and_log("You does not set 'id' of relaxed calculation. I try to find geometry files in "+it_new+" folder\n")

        if it_to:       geo_path = geo_folder + struct_des[it_to].sfolder  + "/" + it_to
        else:           geo_path = geo_folder + struct_des[it_new].sfolder + "/" + it_new+'/from'
        if copy_geo_from:
            print_and_log("You asked to copy geo files from "+copy_geo_from+" to " +geo_path+ " folder\n") 
            #if not os.path.exists(os.path.dirname(geo_path)): 
            runBash( "mkdir -p "+geo_path ) 
            runBash( "cp "+copy_geo_from+"/* "+geo_path  )




        if os.path.exists(geo_path):
            print_and_log("Folder '"+geo_path +"' was found. Trying to add impurity\n");
        else:
            print_and_log("Error! Folder "+geo_path+" does not exist\n"); raise RuntimeError

        #geofilelist = glob.glob(geo_path+'/*.geo*') #Find input_geofile
        #geofilelist = runBash('find '+geo_path+' -name "*grainA*.geo*" ').splitlines()
        #geofilelist = runBash('find '+geo_path+' -name "*.geo*" ').splitlines()
        geofilelist = glob.glob(geo_path+'/*.geo*')
        print_and_log( "There are several files here already: ", geofilelist, imp = 'y' )
        #print 'find '+geo_path+' -name "*.geo*" ',geofilelist
        #return


        for input_geofile in geofilelist:

            v = int( runBash("grep version "+str(input_geofile) ).split()[1] )
            
            if only_version and v not in only_version: continue # only_version = None works for all versions 
            
            new = CalculationVasp()
            new.version = v
            new.name = input_geofile
            
            

            new.read_geometry(input_geofile)
            init = copy.deepcopy(new)
            
            igl = input_geofile.split("/")
            #new.name = igl[-3]+'/'+igl[-3] #+input_geofile
            new.name = struct_des[it_new].sfolder+"/"+it_new+"/"+it_new
            print_and_log( "New path and part of name of file is ", new.name, imp = 'Y')
            #return
            new.des = 'Obtained from '+input_geofile+' by adding '+impurity_type+' impurity '
            #new.init.xred   = new.xred
            #new.init.rprimd = new.rprimd

            #print new.rprimd
            new.init.name = new.name+'.imp.'+addtype+'.'+str(new.version)
            #new.path["input_geo"] = geo_folder+it_new+"/"+new.end.name+'.'+'geo'
            new.path["input_geo"] = geo_folder+"/"+new.init.name+'.'+'geo'
            #new.init.name = 'test_before_add_impurity'
            
            new = add(znucl, "", new, write_geo, put_exactly_to = put_exactly_to)



    return new.path["input_geo"] #return for last version






def insert_cluster(insertion, i_center, matrix, m_center):
    """
    Take care of orientation; typat should be consistent
    Input:
    insertion -  object of class Structure(), which is supposed to be inserted in matrix
    in such a way that i_center will be combined with m_center.
    matrix - object of class Structure().
    i_center, m_center - numpy arrays (3) cartesian coordinates
    """
    ins = copy.deepcopy(insertion)
    mat = copy.deepcopy(matrix)
    r = mat.rprimd



    hproj = [ (r[0][i]+r[1][i]+r[2][i]) * 0.5 for i in (0,1,2) ] #projection of vectors on three axis
    if 1:
        for i, x in enumerate(ins.xcart):
            ins.xcart[i] = x - i_center

        for i, x in enumerate(mat.xcart):
            mat.xcart[i] = x - m_center

    max_dis = 1
    for i_x, ix in enumerate(ins.xcart):
        dv_min = max_dis
        print_and_log( "Insertion atom ",ix,)
        if 1:
            for j, mx in enumerate(mat.xcart):
                dv = mx - ix
                for i in 0,1,2:
                    if dv[i] >  hproj[i]: dv = dv - mat.rprimd[i] #periodic boundary conditions - can be not correct (in the sense that closest image can lie not 100 % in the neighbourhood image cell ) for oblique cells and large absolute values of dv 
                    if dv[i] < -hproj[i]: dv = dv + mat.rprimd[i]
                
                len1 = np.linalg.norm(dv)
                len2, second_len2 = mat.image_distance(mx, ix, r, 1) #check len1
                

                #print "Lengths calculated with two methods ", len1, len2
                len1 = len2 #just use second method
                #assert np.around(len1,1) == np.around(len2,1)

                if len1 < dv_min: 
                    dv_min = len1;   
                    j_r = j # number of matrix atom to replace





        if 1:
            #Modify to replace overlapping atoms 
            if dv_min == max_dis:
                print_and_log( " is more far away from any matrix atom than ",dv_min," A; I insert it")
                # mat.xcart.append( ix )
                # print_and_log( 'type of added atom is ', ins.typat[i_x])
                # mat.typat.append( ins.typat[i_x]   )
                mat = mat.add_atom(xc = ix, element = ins.get_elements()[i_x] )



            else:        
                print_and_log( "will replace martix atom", mat.xcart[j_r] )
                mat.xcart[j_r] = ix.copy()
    






    mat.rprimd = r
    mat.xcart2xred()
    mat.natom = len(mat.xcart)
    mat.name = 'test_of_insert'
    st = mat
    # print(st.natom, len(st.xcart), len(st.typat), len(st.znucl), max(st.typat) )
    # write_xyz(mat)
    mat = mat.return_atoms_to_cell()
    mat.write_poscar()
    return mat
    #write_xyz(mat)


def make_interface(main_slab, m_xc, second_slab, s_xc):
    """
    Make interfaces
    Both slabs should have close sizes along x and y and should be oriented correctly
    
    Input:
    main_slab (Structure) - slab
    second_slab (Structure) - slab, scaled to coincide with the main slab
    m_xc, s_xc (array(3)) - cartesian coordinates of pointis in main_slab and secondary slab to be combined
    

    Return Slab with interface and scaled second slab
    """
    ins = copy.deepcopy(second_slab)
    mat = copy.deepcopy(main_slab)


    if 1:
        #scale insertion
        mr = mat.rprimd_len()
        ir = ins.rprimd_len()
        print('Matrix vlength', mr)
        print('Insert vlength', ir)
        x_scale = mr[0]/ ir[0]
        y_scale = mr[1]/ ir[1]

        print('Scaling factors', x_scale, y_scale)

        # print('i_center', i_center)

        ins.rprimd[0] = ins.rprimd[0]*x_scale
        ins.rprimd[1] = ins.rprimd[1]*y_scale
        ir = ins.rprimd_len()
        s_xred = xcart2xred([s_xc], ins.rprimd)[0]
        print('Insert vlength after scaling', ir)

        ins.update_xcart()
        # ins.xcart2xred()
        ins_sc = ins.copy()
        ins_sc.name+='_scaled'

        s_xc = xred2xcart([s_xred], ins.rprimd)[0]
     
        # print('i_center', i_center)



    if 1:
        for i, x in enumerate(ins.xcart):
            ins.xcart[i] = x - s_xc

        for i, x in enumerate(mat.xcart):
            mat.xcart[i] = x - m_xc

    for i_x, ix in enumerate(ins.xcart):

        mat = mat.add_atom(xc = ix, element = ins.get_elements()[i_x] ) 


    mat.xcart2xred()
    mat.natom = len(mat.xcart)
    mat.name += 'inteface'
    mat = mat.return_atoms_to_cell()
    mat = mat.shift_atoms([0,0,0.5])

    return mat, ins_sc





def make_interface2(st1, st2, shift, mesh = [5,5], oxi = {}, at_fixed = [], mode = "top", tol = 0.05):
    """
    Creates interface from two input structures

    INPUT:
    st1 (Structure) - main input structure 
    st2 (Structure) - appending input structure 
    thickness (float) - remaining thickness of vacuum  
    shift (float array) - shift atoms along axis 
    tol (float) - relative difference between two steps with different meshes 
    mode (str) - type of interfaces or interfaces that will be returned 
        'top' - returns one structure that was made from highest atoms from st1 and lowest atom from st2 
        'min' - returns

    'min' mode requires these arguments:
        mesh (int array) - mesh in AB plane to find a structure with lowest Ewald energy
        oxi (float dic) - dictionary oxidation states. Look siman.geo.make_neutral() for more details 
            E.g {"Li": 1, "La": 2, "Zr":4, "O": -1}
        at_fixed (string array) - atoms with fixed oxidation states 

    RETURN:
        interface and new structure st2. 

    author - A. Burov 

    """
    from pymatgen.analysis.ewald import EwaldSummation

    st1 = copy.deepcopy(st1)
    st2 = copy.deepcopy(st2)
    st1_tmp = copy.deepcopy(st1)
    st2_tmp = copy.deepcopy(st2)

    st1_tmp = st1_tmp.remove_vacuum()
    st2_tmp = st2_tmp.remove_vacuum()
    z1 = st1_tmp.rprimd[2][2]
    z2 = st2_tmp.rprimd[2][2]
    st2 = st2.add_vacuum(vector_i = 2, thickness = z1 - z2)

    z_coord1 = [x[-1] for x in st1.xcart]
    z_coord2 = [x[-1] for x in st2.xcart]
    z1_max = [x for _, x in sorted(zip(z_coord1, range(st1.natom)))][-1]
    z2_min = [x for _, x in sorted(zip(z_coord2, range(st2.natom)))][0]

    block_print()
    interface, li_super_new = make_interface(st1, st1.xcart[z1_max], st2, st2.xcart[z2_min]+shift)
    interface = interface.remove_vacuum(thickness = abs(shift[-1]))
    z3 = interface.rprimd[2][2]
    enable_print()

    structures = []
    if (mode == "min"):
        print("Finding structure with minimal Ewald energy")
        shift_x = interface.rprimd[0][0] / mesh[0]
        shift_y = interface.rprimd[1][1] / mesh[1]

        en_min = 1e20
        diff = 1.0
        itr = 0

        while (diff > tol):
            en_prev = en_min
            en_min = 1e20
            itr += 1 
            for x in range(mesh[0]):
                for y in range(mesh[1]):
                    block_print()
                    interface_c = copy.deepcopy(interface)

                    for xc in interface_c.xcart:
                        if (xc[2] > z3 - z2  + (z1 - z2)):
                            xc[0] += shift_x * x
                            xc[1] += shift_y * y

                    interface_c.update_xred()
                    interface_c = interface_c.return_atoms_to_cell()
                    interface_c = interface_c.make_neutral(oxidation = oxi, at_fixed = at_fixed, mode = 'equal',
                                                            silent = 1, return_oxidation = 0)

                    en = EwaldSummation(interface_c)
                    en = en.total_energy
                    if (en < en_min) and (x != 0) and (y != 0):
                        en_min = en
                        move_x = x
                        move_y = y
                    enable_print()
                    print(en)
                 
            diff = abs((en_min - en_prev) / max(abs(en_min), abs(en_prev)))
            shift_x /= 2
            shift_y /= 2
            print("Interation: {}, difference: {}".format(itr, diff))
        
        block_print()
        interface_c = copy.deepcopy(interface)

        for xc in interface_c.xcart:
            if (xc[2] > z3 - z2  + (z1 - z2)):
                xc[0] += shift_x * x
                xc[1] += shift_y * y

        interface_c.update_xred()
        interface_c = interface_c.return_atoms_to_cell()
        enable_print()
        
    else:
        pass

    if (mode == "top") or (mode == "min"):
        return interface_c, li_super_new
    else:
        raise ValueError("Wrong mode, check function's description.")





def insert(it_ins, ise_ins, mat_path, it_new, calc, type_of_insertion = "xcart" ):
    """For insertion of atoms to cells with changed lateral sizes
    Input:
    'type_of_insertion = xred' used to add xred coordinates  
    mat_path - path to geo files which are supposed to be changed
    it_ins - already existed calculation; xred will be used from this calculation.
    it_new - new folder in geo folder for obtained structure
    
    This function finds version of calculation in folder mat_path and tries to use the same version of it_ins

    """
    if not os.path.exists(mat_path):
        print_and_log("Error! Path "+mat_path+" does not exist\n\n")
        raise RuntimeError

    if it_ins not in mat_path and it_ins not in it_new: 
        print_and_log('Cells are', it_ins, mat_path, it_new)
        print_and_log("Error! you are trying to insert coordinates from cell with different name\n\n")
        #raise RuntimeError       

    hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
    if hstring != header.history[-1]: header.history.append( hstring  )

    geofilelist = runBash('find '+mat_path+'/target -name "*.geo*" ').splitlines()
    
    if geofilelist == []:
        print_and_log("Warning! Target folder is empty. Trying to find in root folder ...")
        geofilelist = runBash('find '+mat_path+'/ -name "*.geo*" ').splitlines()

    ins = None
    for mat_geofile in geofilelist:
        mat = CalculationVasp()
        mat.name = mat_geofile
        mat.read_geometry(mat_geofile)
        #step = 0.27
        #r_pore = 0.56
        #r_mat = mat.hex_a / 2 - step
        #pores = find_pores(mat.init, r_mat, r_pore, step, 0.3, 'central') #octahedral
        #mat.xcart.append ( pores.xcart[0] )
        #mat.typat.append(1)
        try:
            ins_working = ins
            ins = calc[(it_ins, ise_ins, mat.version)]
        except KeyError: 
            print_and_log( "No key", (it_ins, ise_ins, mat.version), "I use previous working version !!!", imp = 'y' )
            ins = ins_working
            #return
        #ins.end.znucl = ins.znucl
        #ins.end.nznucl = ins.nznucl
        #ins.end.ntypat = ins.ntypat
        #ins.end.typat = ins.typat
        #print ins.xcart[-1]
        mat_geopath = geo_folder+struct_des[it_new].sfolder + '/'

        if type_of_insertion == "xcart":
            #Please update here!
            mat_filename = '/'+it_new+"."+"inserted."+str(mat.version)+'.'+'geo'
            
            v = np.zeros(3)
            result = insert_cluster(ins.end, v, mat.init, v )
            mat.end = result
            mat.init = result
            # mat.znucl  =   mat.end.znucl
            # mat.nznucl =   mat.end.nznucl
            # mat.ntypat =   mat.end.ntypat
            # mat.typat  =   mat.end.typat
            # mat.natom = len(mat.end.xred)    
            #mat.version = ins.version
            des = ins.name+" was inserted to "+mat_geofile
        
        elif type_of_insertion == "xred":

            mat_filename = '/from/'+it_new+".xred."+str(mat.version)+'.'+'geo'
          
            #mat.end.rprimd = mat.rprimd
            #mat.init.xred  = copy.deepcopy(ins.end.xred)
            #mat.init.typat = copy.deepcopy(ins.end.)
            #print ins.end.xcart
            rprimd   = copy.deepcopy(mat.init.rprimd)
            #build    = mat.build
            mat.init = copy.deepcopy(ins.end)
            #mat.build = build
            mat.init.rprimd = rprimd #return initial rprimd
            mat.init.xred2xcart() #calculate xcart with new rprimd
          
            des = "atoms with reduced coord. from "+ins.name+" was fully copied to "+mat_geofile
            mat.init.name = 'test_insert_xred'+str(mat.version)
            write_xyz(mat.init)


        mat.path["input_geo"] = mat_geopath + it_new + mat_filename
        if not mat.write_geometry("init",des): continue
        print_and_log("Xred from "+it_ins+" was inserted in "+mat_geofile+" and saved as "+mat_filename+" \n\n")

    return
    





def determine_voids(st, r_impurity, fine = 1, step_dec = 0.05):

    if not r_impurity:
        printlog('add_neb(): Error!, Please provide *r_impurity* (1.6 A?)')


    sums = []
    avds = []
    printlog('Searching for voids', important = 'y')
    st_pores = find_pores(st, r_matrix = 0.5, r_impurity = r_impurity, step_dec = step_dec, fine = fine, calctype = 'all_pores')

    printlog('List of found voids:\n', np.array(st_pores.xcart) )
    write_xyz(st.add_atoms(st_pores.xcart, 'H'), file_name = st.name+'_possible_positions')
    write_xyz(st.add_atoms(st_pores.xcart, 'H'), replications = (2,2,2), file_name = st.name+'_possible_positions_replicated')

    for x in st_pores.xcart:
        # summ = local_surrounding(x, st, n_neighbours = 6, control = 'sum', periodic  = True)
        # avd = local_surrounding(x, st, n_neighbours = 6, control = 'av_dev', periodic  = True)
        summ, avd = local_surrounding2(x, st, n_neighbours = 6, control = 'sum_av_dev', periodic  = True)
        # print (summ, avd)
        
        sums.append(summ)
        avds.append(avd[0])
    # print
    sums = np.array(sums)
    avds  = np.array(avds).round(0)

    print_and_log('Sum of distances to 6 neighboring atoms for each void (A):\n', sums, imp ='y')
    print_and_log('Distortion of voids (0 - is symmetrical):\n', avds, imp ='y')
    
    return st_pores, sums, avds

def determine_unique_voids(st_pores, sums, avds):
    crude_prec = 1 # number of signs after 0


    sums_crude = np.unique(sums.round(crude_prec))



    print_and_log('The unique voids based on the sums:', 
        '\nwith 0.01 A prec:',np.unique(sums.round(2)),
        '\nwith 0.1  A prec:',sums_crude,
        imp ='y')
    print_and_log('Based on crude criteria only', len(sums_crude),'types of void are relevant', imp = 'y') 


    insert_positions = []
    start_table = []
    for i, s in enumerate(sums_crude):
        index_of_first =  np.where(sums.round(crude_prec)==s)[0][0]

        start_table.append([i,  st_pores.xcart[index_of_first].round(2), index_of_first,
        avds[index_of_first], sums[index_of_first]     ])

        insert_positions.append( st_pores.xcart[index_of_first] )

    print_and_log( tabulate(start_table, headers = ['void #', 'Cart.', 'Index', 'Dev.', 'Sum'], tablefmt='psql'), imp = 'Y' )
    
    return insert_positions

def insert_atom(st, el, i_void = None, i_void_list = None, r_imp = 1.6, ):
    """Simple Wrapper for inserting atoms 

    i_void (int) has higher priority than i_void_list
    
    return st_new, i_add, sts_by_one
        st_new - all positions are filled 
        i_add - the number of last inserted atom
        sts_by_one - list of structures with only one inserted atom in all found positions


    """


    r_impurity = r_imp
    st_pores, sums, avds = determine_voids(st, r_impurity)

    insert_positions = determine_unique_voids(st_pores, sums, avds)

    printlog('To continue please choose *i_void* from the list above', imp = 'y')


    # st.name = st.name.split('+')[0]

    if i_void:
        i_void_list = [i_void]

    if i_void_list is None:
        i_void_list = list(range(len(insert_positions)))
        printlog('No i_void was provided, I insert all', imp = 'y')

    st_new = st.copy()
    sts_by_one = []
    for i in i_void_list:
        xc = insert_positions[i]
        
        st_new, i_add = st_new.add_atoms([xc], el, return_ins = True)
        st_one, _ = st.add_atoms([xc], el, return_ins = True)
        st_one.name+='+'+el+str(i)
        sts_by_one.append(st_one)

        st_new.name+='+'+el+str(i)
        st_new.des+=';Atom '+el+' added to '+ str(xc)
    
    printlog(st.des, imp = 'y')

    st_new.write_poscar()
    st_new.magmom = [None]

    return st_new, i_add, sts_by_one




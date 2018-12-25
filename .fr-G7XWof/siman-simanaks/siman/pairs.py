# -*- coding: utf-8 -*- 
#Copyright Aksyonov D.A
 
from operator import itemgetter

from siman import header
from siman.classes import CalculationVasp, Structure
from siman.functions import  element_name_inv
from siman.geo import image_distance, replic, local_surrounding, xcart2xred, xred2xcart
from siman.inout import write_xyz
from siman.impurity import insert, add_impurity



def PBC(dx, r):
    """
    can be incorrect
    Realisation of periodic boundary conditions in common case
    dx - vector[3] difference between coordinates of two atoms
    r - rprimd of cell
    return dx - the smallest distance between atoms

    """
    
    hproj = [ (r[0][i]+r[1][i]+r[2][i]) * 0.5 for i in (0,1,2) ] #projection of vectors on three axis
    
    dxl = dx.copy()
    for i in 0, 1, 2:
        #the smallest distance: can be incorrect for oblique cells
        if dx[i] >  hproj[i]: dx = dx - r[i] #periodic boundary conditions
        if dx[i] < -hproj[i]: dx = dx + r[i]
        #the largest distance: not correct
        #if dx[i] <  0: dxl = dx + r[i] 
        #f dx[i] >  0: dxl = dx - r[i]

    return math.sqrt(dx[0]**2 + dx[1]**2 + dx[2]**2) #, math.sqrt(dxl[0]**2 + dxl[1]**2 + dxl[2]**2)





def find_pairs(base_name, segtyp, in_calc, central_atoms = [], xcart1imp = None, input_dlist_coseg = None, prec = 2, gvolume_config_num = None,
    gbpos = None,
    take_final_rprimd_from = None, main_path = None, based_on = None, target_znucl = [22, 6, 8],
    max_dist_between_atoms = 4.8, add_typat = [2, 3]                 ):
    """
    Find uniq pairs of atoms and analyse them
    Input:
    
    segtyp - 
    three regimes for cells with grain boundaries:
    'segreg' assumes that in_calc contains carbon atom in grain volume, and creates all cases; 
    'coseg' assumes pure cell and creates only coseg cases.
    cosegregation cases of course should be the same for two regimes, however co-segregation configuations after 'coseg' is more easy to relax.
    'grainvol' - searching for pairs in grain volume

    two regimes for bulk cells:
    'bulk_triple' - used for bulk cells without grain boundaries; first step is searching for pairs, second step for triples.
    'bulk_pairs' - used for bulk cells without grain boundaries; searching for pairs.



    new_name - name of created structures; at first should be added to struct_des[]
    in_calc - Calculation() type or path to geo file
    region - list of numbers which determine region
    central_atoms - list of atoms for which pairs are constructed (Warinig! numbers in new array xcart_pores!);
    - parameter to change mode;

    xcart1imp - coordinates of first interstitial in the  grain interior 

    input_dlist_coseg - list of configurations with cosegregation cases. Needed to construct corresponding segregation cases. the format is quiet tricky

    prec - precision of lengths used to determine unique positions.

    gvolume_config_num - number of configuration with two atoms in grain volume choosen by user (usually should be the most favourable) 

    gbpos - position of grain boundary

    take_final_rprimd_from - path to geo file from which rprimd will be used


    target_znucl - numbers of target atoms

    max_dist_between_atoms - now at least used for 'bulk_pairs' and 'bulk_triple'; maximum length of found pairs.

    add_typat - mannualy set please update

    """

    def write_geometry_files(dlist, in_calc, xcart_pores, segtyp, take_final_rprimd_from = None, add_name_before = '', tlist = [], configver = False,
        add_typat = None):
        """Creating files
        dlist - list of pairs with distances and numbers
        in_calc - base calculation without pores

        tlist - list of additional atoms in the case of triples; list of structures
        configver - if True each configuration saved as a new version
        add_typat - manually determined; please automatize!
        """

        print(( "Warning! add_typat", add_typat))
        if tlist == []: #convert dlist to tlist - to do earlier
            for el in dlist:
                config = Structure()
                config.name = el[2]
                config.length = el[0]
                config.typat = add_typat
                config.xcart = [el[7], el[8] ]
                tlist.append(config)


        for i, el in enumerate( tlist): #by all found structures
            print(( "el name is ", el.name))
            stn = copy.deepcopy(in_calc.init)
            calc = copy.deepcopy(in_calc)

            stn.typat.extend(el.typat)
            stn.xcart.extend(el.xcart)


            stn.xred = xcart2xred(stn.xcart, stn.rprimd)

            xcart_check = xred2xcart(stn.xred, stn.rprimd)        
            assert len(xcart_check) == len(stn.xcart) #test
                  
            assert all([ all( np.around(v1, 8) == np.around(v2, 8) ) for (v1, v2) in zip(stn.xcart, xcart_check) ]) #check if xcart2xred(stn.xcart,r) and xred2xcart(stn.xred,r) are working correctly up to the eight digits after
            
            stn.natom = len(stn.xcart)


            

            """Adapt new rprimd"""
            print("take final rprimd is ", take_final_rprimd_from)
            if take_final_rprimd_from: #read final rprimd and version
                print("Start to read rprimd and version from "+take_final_rprimd_from)
                in_calc_rprimd = CalculationVasp()
                in_calc_rprimd.name = 'temp'
                in_calc_rprimd.read_geometry(take_final_rprimd_from)
                
                stn.rprimd = in_calc_rprimd.init.rprimd
                stn.xcart = xred2xcart(stn.xred, stn.rprimd)            

                calc.version = in_calc_rprimd.version
            elif configver:
                calc.version = i+1

            
            calc.init = stn
            
            des = ' was obtained by the insertion of C-O pair into '+in_calc_name+'; final vectors taken from corresponding ver'
            
            calc.build.ipairlen = el.length # Initial length of pair
            if    not hasattr(calc.build, 'nadded') or calc.build.nadded == None:      calc.build.nadded = 2
            else: calc.build.nadded+= 2
            if    not hasattr(calc.build, 'listadded') or calc.build.listadded == [None]: calc.build.listadded = list(range(stn.natom - 2, stn.natom)) #list of atoms which were added
            else: calc.build.listadded.extend(                           list(range(stn.natom - 2, stn.natom)) )

            structure_name = calc.name+el.name.split('.')[0]
            #calc.name = add_name_before+calc.name+ '.' +el[2]+'.'+str(calc.version)
            print('Structure_name', structure_name)
            if   structure_name in      struct_des:     

                if configver:
                    fname = structure_name# calc.name+'C2O2'
                    calc.path["input_geo"] = geo_folder + struct_des[fname].sfolder + '/' + fname + '/' + structure_name + '.' + segtyp+ '.' +str(calc.version) + '.geo'
                else:
                    calc.path["input_geo"] = geo_folder + struct_des[structure_name].sfolder + '/' + structure_name + '/' + structure_name + '.' + segtyp+ '.' +str(calc.version) + '.geo'
                print("write geo to ", calc.path["input_geo"])
                calc.write_geometry('init', des )


            print("write_geometry_files(): name ", el.name)
            stn.name = add_name_before+calc.name+ '' +str(el.name)+'.'+str(calc.version)
            #stn = replic(stn, (1,2,2))
            write_xyz(stn)
            print("__________________________\n\n\n")
        return

    def min_diff(f, list, diffprec):
        """
        calculates difference between one number and list of other numbers. return the index of number with smallest difference.
        if difference is smaller than diffprec returns true as the second argument.  
        """
        #print list
        if list == []: return 0, False
        mind  = min([abs(f - l) for l in list])
        with_i = np.asarray( [abs(f - l) for l in list] ).argmin()
        return with_i, (mind < diffprec)





    def pairs(in_calc, xcart_pores, central_atoms, prec = 2, max_dist = 20, max_dist_from_gb = 4 , pairtyp = 'gvol'):
        """
        Searhing for pairs and make list of distances and numbers of atoms
        prec - precision, allows to control which distances can be related to the same configurations
        max_dist - maximum distance between atoms in pair
        max_dist_from_gb - 
        pairtyp - 'gvol' assumes that central_atoms are in the grain volume, 'gb' assumes that central_atoms are in the grain boundary region

        """
        st = in_calc.init
        st_replic = replic(st, (2,2,2))
        st_replic = replic(st_replic, (2,2,2), -1) #replic in negative direction also
        r1x   = in_calc.rprimd[0][0]
        r3z   = in_calc.rprimd[2][2]
        print("Half length of r1x is", r1x/2)

        if segtyp in ['segreg', 'coseg', 'grainvol']:
            gbpos2 = in_calc.gbpos
            gbpos1 = gbpos2 - r1x/2.
            print("\n\nPositions of boundaries gb1 and gb2",gbpos1, gbpos2)
            print("Maximum possible distance between boundary and impurity", r1x/4)
        else:
            gbpos2 = 0
            gbpos1 = 0      


        dlist = []
        d1list = []
        d2list = []
        dgb2list = []

        
        n_neighbours = 8 # number of atoms to calculate sums

        sumrulist = [] #list of sums (sumr1 or sumr2) of unique pores
        unique_pores = [] #the same list but also with coordinates of pores        
        
        sumrlist = [] #list of sumr1+sumr2
        
        k = 1
        d2diff = 0
        d1diff = 0
        #z1 = 6 #charge of added impurity
        #z2 = 8

        diffprec = 0.02

        
        



        # print xcart_pores
        



        for i, x1 in enumerate(xcart_pores):
            if i not in central_atoms: continue
            #iz = z1
            for j, x2 in enumerate(xcart_pores):
                if all(x1 == x2): continue
                
                d = abs(x2[0]-in_calc.gbpos)
                if pairtyp == 'gb' and d > max_dist_from_gb: continue #second atom is too far from grain boundary 
                
                d1, d2 = image_distance(x1, x2, st.rprimd, 2) # the minimum distance and the next minimum dist
                if d1 > max_dist: continue
                if (d1, d2) != image_distance(x1, x2, st.rprimd, 3): raise RuntimeError  #test, searching in father images
                
                #d1 = round(d1,prec)
                #d2 = round(d2,prec)
                dgb1 = round(x2[0]-gbpos1,prec)
                dgb2 = round(gbpos2-x2[0],prec)
                
                sumr1 = local_surrounding(x1, st_replic, n_neighbours) # sum of distances to surrounding atoms
                sumr2 = local_surrounding(x2, st_replic, n_neighbours)
                sumr  = sumr2 + sumr1

                if sumr1 not in sumrulist:
                    sumrulist.append(sumr1) 
                    unique_pores.append((sumr1, x1) ) #determine unique pores 

                if sumr2 not in sumrulist:
                    sumrulist.append(sumr2) 
                    unique_pores.append((sumr2, x2) ) #determine unique pores 

                #if d1 in d1list: continue
                if sumr in sumrlist:# new condition based on sumr
                    ind = sumrlist.index(sumr)
                    i_min, smaller = min_diff(d1, d1list, diffprec)
                    if smaller: continue
                    
                # if 0:#d1list: 
                #     i_min, smaller = min_diff(d1, d1list, diffprec)# d1 has the smallest difference with di
                #     #print "exist"
                #     d2diff = abs(d2list[i_min]-d2)
                #     #print abs(d2list[i_min]-d2)
                #     #print central_atoms
                #     if smaller and abs(d2list[i_min]-d2) < diffprec*2  : continue #and abs(dgb2list[i_min]-dgb2) < diffprec
                    
                #     i_min, smaller = min_diff(d2, d2list, diffprec)# d1 has the smallest difference with di
                #     d1diff = abs(d1list[i_min]-d1)
                #     if smaller and abs(d1list[i_min]-d1) < diffprec*2  : continue

                    #print "skiped"
                #di, smaller = min_diff(d2, d2list, diffprec)
                #if di != None and smaller: continue
                #if min_diff(d2, d2list, diffprec): continue # be carefull here. this condition can pass some unique configrations; should make additional check like below
                #if d2 in d2list and dgb2list[d2list.index(d2)] == dgb2: continue
                #jz = z2

                sumrlist.append(sumr)
                d1list.append(d1)
                # d2list.append(d2)
                # dgb2list.append(dgb2)

                
                sym = ''
                if 0: #mannualy switched off
                    if abs(x1[1]-x2[1]) < diffprec: #Find symmetry
                        if abs(x1[2]-x2[2]) < diffprec:
                            sym = 'ms' # if y and z are the same, than mirror symmetry
                        elif abs(x1[2]-x2[2]) - r3z < diffprec:
                            sym = 'is' # inverse symmtry
                        elif abs(x1[2]+x2[2]) - 0.5*r3z < diffprec: # only for t111g; should be extended for general case of existing periods along y or z
                            sym = 'is'


                dlist.append([round(d1,prec), round(d2, prec), sym, sumr1, sumr2 , dgb1, dgb2, x1, x2, sumr1, sumr2] ) #the first sumr1, sumr2 below replaced by their types
                                



                k+=1

        dlist.sort(key = itemgetter(0))

        unique_pores.sort(key = itemgetter(0))
        sumrulist.sort()
        print('Number of unique pores     is', len(unique_pores))
        print('Pores have the following sums: ',unique_pores)



        print("Searching for similar pairs but with different distances ...")
        print("number, d1, d2, name,  sumr1, sumr2, dgb1, dgb2; parallel pair with larger distances")
        
        bname = element_name_inv(target_znucl[1]) + element_name_inv(target_znucl[2])
        for i, el1 in enumerate(dlist):

            typ1 = sumrulist.index(el1[3])+1 #typ of pore of the first atom
            typ2 = sumrulist.index(el1[4])+1 
            el1[3] = typ1
            el1[4] = typ2
            
            if pairtyp == 'gb':
                dlist[i][2] = bname+'i'+str(i+1)+'.'+str(el1[3])+'-'+str(el1[4]) + dlist[i][2]
            
            elif pairtyp == 'gvol':
                dlist[i][2] = bname+'.v'+str(i+1) + dlist[i][2]

            print(i+1, el1[:3], el1[-2:] , el1[-6], el1[-5], end=' ') #number, d1, d2, name,  sumr1, sumr2, dgb1, dgb2

            for el2 in dlist: #this loop looks for pairs which are parallel to the same direction as el1 but have larger interdistances
                
                if el1 == el2: continue
                
                mod = el2[0]/el1[0] % 1
                
                if ( mod < 0.005 or mod > 0.995 ) and abs(el1[0]-el2[0]) > dlist[0][0]: #only multiple distances and if difference is larger than smallest distance
                    #if round(el1[2],prec-1) != round(el2[2],prec-1): continue #In either case the sum the distances should be the same for the same direction
                    if el1[0] == el2[1]: continue
                    print(el2[0]/el1[0], end=' ') # el2, this pair of atoms is analogus to el1 but have larger interdistance
            print()
        print('Total number of structures is', len(dlist))


      
        if 0:
            print("\n\nSearching for pairs with equal distances by periodic boundary conditions:")
            for el1 in dlist:
                if el1[0] == el1[1]:
                    print(el1)        

            print("\nSearching for pairs with not equal distances by periodic boundary conditions:")
            for el1 in dlist:
                if el1[0] != el1[1]:
                    print(el1)  


            print("\nSearching for pairs with d2/d1>2:")
            for el1 in dlist:
                if el1[1]/el1[0] > 2:
                    print(el1) 


        dlist[0].append(unique_pores) # last element of dlist[0] is sum and coordinates of unique pores

        return dlist










    """0. BEGIN-------------------------------------------------------------------------------"""

    hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
    if hstring != header.history[-1]: header.history.append( hstring  )

    print_and_log("\n\n------Starting find_pairs()-----------...\n")


    if type(central_atoms) == int: #not in [tuple, list]:

        central_atoms = [central_atoms]; #transform to list


    if type(in_calc) == str:
        in_calc_name = in_calc
        in_calc = CalculationVasp()
        #in_calc.name = str(in_calc_name)
        
        in_calc.name = base_name 
        print("in_calc name is ", in_calc.name)
        

        in_calc.read_geometry(in_calc_name)
        if gbpos: in_calc.gbpos = gbpos    #rewrite gbpos

        st = in_calc.init
    else:
        """End relaxed structure is used!!!"""
        st = copy.deepcopy(in_calc.end)
        in_calc_name = str(in_calc.id)






    

    """1. Create separate list of pores and remove them from xcart--------------------------------------------------------"""
    
    if "hcp_octa_xred":
        in_calc.init.name = segtyp+'_all_pores'
        rep = replic(in_calc.init, (2,2,2), -1); write_xyz(rep) #just to check
        
        """Coordinates of octapores provided in xcart; znucl = 200;"""
        xcart = st.xcart; typat = st.typat; st.typat = []; st.xcart = []
        xcart_pores = []

        #clean structure from pores with z == 200 and construct xcart_pores
        for i, x in enumerate(xcart):
            z = st.znucl[ typat[i]-1  ]
            if z == 200:
                xcart_pores.append( x )
                #print "Found pore"
            else:
                st.xcart.append( x )
                st.typat.append( typat[i] )
        st.natom = len(st.xcart)
        print('Number of found pores with znucl = 200 is ', len(xcart_pores))
        for n in central_atoms:
            if n >= len(xcart_pores):
                raise RuntimeError



    """2. Preprocess segreg and grainvol cases--------------------------------------------------------"""
    # in_calc can be of two types: pure and with C in grain volume; using pure we construct co-segregation cases; using carbon in volume we can construct segregation cases

    if segtyp in ('segreg', 'grainvol'):

        if 2 in st.typat:# impurity in grain volume; (now assume that Carbon)

            iimp = st.typat.index(2); 

            xcart1imp = st.xcart[iimp]  #save coordinates of carbon atom 
            
            del st.xcart[iimp]; del st.typat[iimp]; st.natom-=1   #and remove it
            #del st.xred[iimp]
            st.ntypat -= 1
            del st.znucl[1]

            print("Impurity atom was removed from cell")  



        if xcart1imp: #for compatibility with previous cases; better not to use

            imp1 = len(xcart_pores); 
            xcart_pores.append( xcart1imp );


            xcart2imp = xcart1imp - 0.5 * st.rprimd[0]  #determine coordinates of second impurity assuming that we have mirror symmetry
            if xcart2imp[0] < 0 : xcart2imp = xcart1imp + 0.5 * st.rprimd[0]        

            imp2=imp1+1 
            xcart_pores.append( xcart2imp )


        else: #new version; both central pores are found in the pure cell!!!

            #We have pure cell here; Find central pore in 1st grain and 2nd grain:
            xcen1 = in_calc.gbpos - 0.25 * st.rprimd[0][0] #x center of the first grain
            xcen2 = in_calc.gbpos - 0.75 * st.rprimd[0][0] #z center of the second grain
            # print "xcen", xcen1, xcen2

            d1l = []; d2l = []; rpxx05 = st.rprimd[0][0]*0.5
            for x in xcart_pores:
                d1 = xcen1 - x[0]
                d2 = xcen2 - x[0] 
                if d2 < -rpxx05: d2 += st.rprimd[0][0] # assuming that periodic boundary conditions needed only here
                
                d1l.append(abs(d1) )
                d2l.append(abs(d2) )
                # print d1,d2
            imp1 = np.argmin(d1l) #needed numbers of pores
            imp2 = np.argmin(d2l)
            # print imp1, imp2
            xcart1imp = xcart_pores[imp1]
            xcart2imp = xcart_pores[imp2]
            # print "xcartimp", xcart1imp, xcart2imp









    """3. Define central atoms for segregation and co-segregation cases--------------------------------------------------------"""
    max_dist_from_gb = 100
    if segtyp in ('segreg', 'coseg'):

        # central_atoms = []

        max_dist_between_atoms = 4.8    
        max_dist_from_gb       = 3   #main controls
        
        for i, x in enumerate(xcart_pores): #generate central atoms
            d = abs(x[0]-in_calc.gbpos)
            if d < max_dist_from_gb: 
                central_atoms.append(i)

    """4. Assume that we always have target_znucl, but only three !!!--------------------------------------------------------"""

    st.znucl = target_znucl #Please make this part more general
    st.ntypat = len(set(st.znucl))
    print('Warning! Found only ', st.ntypat, 'of unique atoms in target_znucl')

    st.xred = xcart2xred(st.xcart, st.rprimd)

    """5. Find segreg and co-segreg cases--------------------------------------------------------"""
    in_calc.init = st

    dlist_coseg = []
    if segtyp == 'coseg':
        print("\nStart searching pairs in  gb")
        # main_path = 'T2/CO/' #! please make more general

        
        dlist_coseg = pairs(in_calc, xcart_pores, central_atoms, prec, max_dist_between_atoms, max_dist_from_gb, pairtyp = 'gb' ) 
        

        dlist_coseg_exc = []
        for el in copy.deepcopy(dlist_coseg):         #Exchange C and O only for unsymmetrical cases
            if 's' in el[2]: continue # not needed for symmetrical cases
            el[2] = el[2].replace('C','x'); el[2] = el[2].replace('O','C'); el[2] = el[2].replace('x','O')
            el[7], el[8] = el[8], el[7]
            el[3], el[4] = el[4], el[3]
            dlist_coseg_exc.append(el)

        for el in dlist_coseg+dlist_coseg_exc: # Helper
            stname = base_name+el[2]
            path   = main_path+base_name+'_coseg'
            print(( ("struct_des['{0:s}'] = des('{1:s}', 'co-segregation configurations; made from "+based_on+"'   )").format(stname, path)   )) 

        for el in dlist_coseg+dlist_coseg_exc: # Helper
            stname = base_name+el[2]
            path   = main_path+base_name+'_coseg'
            print("add_loop('"+stname+"','"+based_on.split('.')[1]+"',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')")




        write_geometry_files(dlist_coseg+dlist_coseg_exc, in_calc, xcart_pores, segtyp, take_final_rprimd_from, add_typat = add_typat ) 
    
    elif segtyp == 'segreg':

        """Produce segregation cases only in the case of segreg"""
        print("\nStart producing segragation cases")

        dlist_segreg = []
        # dlist_segreg1 = copy.deepcopy(input_dlist_coseg) #in this case we use input_dlist with co-segragation cases from pure cell. 
        # dlist_segreg2 = copy.deepcopy(input_dlist_coseg) #There is small error, because positions of pores at grain boundary 
        #                                            #differs in pure cell and cell with impurity in grain volume
        # for i, el in enumerate(input_dlist_coseg):
        #     sym = ''
        #     if   'is' in el[2]: sym = 'is'
        #     elif 'ms' in el[2]: sym = 'ms'
        #     dlist_segreg1[i][7] =  xcart1imp
        #     dlist_segreg1[i][2] = 'CvOi'+str(i+1)+sym

        #     dlist_segreg2[i][8] =  xcart1imp
            
        #     dlist_segreg2[i][2] = 'CiOv'+str(i+1)+sym
            
        """new determination based on input_dlist_coseg[0][-1]"""
        el = copy.deepcopy(input_dlist_coseg[0])
        unique = el[-1] #sums and coordinates of unique pores
        print("unique", unique)
        for i, sx in enumerate(unique):
            el[2] = 'Ci'+str(i+1)+'Ov'
            el[7] = sx[1]
            d1, dnext = image_distance(sx[1], xcart1imp, st.rprimd, 2)
            d2, dnext = image_distance(sx[1], xcart2imp, st.rprimd, 2) 
            if d1 > d2:
                el[8] = xcart1imp
            else:
                el[8] = xcart2imp # the farthest impurity in grain volume is used
            

            dlist_segreg.append( copy.deepcopy(el) )



        #dlist_segreg = dlist_segreg1 + dlist_segreg2 #Segregation of the first impurity and of the second


        dlist_segreg_exc = []
        for el in copy.deepcopy(dlist_segreg):         #Exchange C and O only for unsymmetrical cases
            if 's' in el[2]: continue # not needed for symmetrical cases
            el[2] = el[2].replace('C','x'); el[2] = el[2].replace('O','C'); el[2] = el[2].replace('x','O')
            el[7], el[8] = el[8], el[7]
            el[3], el[4] = el[4], el[3]
            dlist_segreg_exc.append(el)    


        #helper
        for el in dlist_segreg+dlist_segreg_exc: 
            stname = base_name+el[2]
            path   = main_path+base_name+'_segreg'
            print(( ("struct_des['{0:s}'] = des('{1:s}', 'co-segregation configurations; made from "+based_on+"'   )").format(stname, path)   )) 

        for el in dlist_segreg+dlist_segreg_exc: 
            stname = base_name+el[2]
            path   = main_path+base_name+'_segreg'

            print("add_loop('"+stname+"','"+based_on.split('.')[1]+"',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')")


        write_geometry_files(dlist_segreg+dlist_segreg_exc, in_calc, xcart_pores, segtyp, take_final_rprimd_from, add_typat = add_typat  ) 

    

    """6. Find volume cases--------------------------------------------------------"""
    # this part for construction volume cases    
    if segtyp == "grainvol": #take care that you have only one carbon atom in the grain
        print("\nStart searching pairs in the volume")
        central_atoms = [imp1]
        max_dist_between_atoms = 4.
        #gvolume_config_num = 0 #please choose manually    
       
        dlist = pairs(in_calc, xcart_pores, central_atoms, prec, max_dist_between_atoms, max_dist_from_gb, pairtyp = 'gvol' ) 

        #dlist = [dlist[0], copy.deepcopy(dlist[gvolume_config_num-1])  ]
        
        #dlist[0][4] = imp2 #no matter wht was dlist[0]; used for vv case
        dlist.append(copy.deepcopy(dlist[0]) )
        dlist[-1][8] = xcart2imp #last element for both atoms in grain volumes
        dlist[-1][2] = 'CvOvms'
        
        #helper
        for el in dlist: 
            stname = base_name+el[2]
            path   = main_path+base_name+'_gvol' #grain volume
            print("add_loop('"+stname+"','"+based_on.split('.')[1]+"',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')")

        for el in dlist: 
            stname = base_name+el[2]
            path   = main_path+base_name+'_gvol' #grain volume
            print(( ("struct_des['{0:s}'] = des('{1:s}', 'co-segregation configurations; made from "+based_on+"'   )").format(stname, path)   )) 



        write_geometry_files(dlist, in_calc, xcart_pores, segtyp, take_final_rprimd_from, add_typat = add_typat  )

    """. Triple cases--------------------------------------------------------"""

    def triples(addatom = ('O', 3), dlist = [], tlist = [], in_calc = None, xcart_pores = [], max_dist_to_next_atom = 3):
        """
        Add addatom to all configurations in tlist; 

        addatom[1] - type of atom in typat
        dlist - list of configurations with two impurity atoms; Used if tlist == []; the format of dlist is quiet special
        tlist - list of configurations with arbirtary number of atoms;

        RETURN:
        tlist - list of configurations with add atoms
        """
        st = in_calc.init
        
        if dlist and tlist == []:
            for el in dlist:
                par = el
                print('pair 1', par, end=' ') 
                x1 = par[7]; x2 = par[8]
                print('x1 = ', x1)
                print('x2 = ', x2)
                config = Structure()
                config.xcart = [x1,x2]
                config.typat = [2, 3]
                config.name  = el[2]
                tlist.append(config)



        tlist_new = []
        for config in tlist:
            xcart = config.xcart
            typat = config.typat
            name  = config.name        
            print('\n\n\nStart to adding atom to ',name)

            i = 1
            dlistlist = []

            diffprec = 0.001


            [dlistlist.append([]) for x in xcart]
            print(len(dlistlist))

            for xpor in xcart_pores:
                
                skip = True
                for j, x in enumerate(xcart): # list of 2 or 3 initial atoms to which additional atom will be added
                    if all(np.around(xpor,5) == np.around(x, 5) ): skip = True; break 
                
                    d1, d2 = image_distance(x, xpor, st.rprimd, 2) # the minimum distance and the next minimum dist
                    if d1 > max_dist_to_next_atom: skip = True; break #if only one pore is larger from atom than limit, the pore is skiped




                    # suml = d11+d21+par[0]
                    # for dl in dlistlist:
                    # print 'j is ', j
                    

                    i_min, smaller = min_diff(d1, dlistlist[j], diffprec) #old condition - bad - removes unique configurations
                    #if smaller: skip = True; continue # symmetrical pores deleted

                    dlistlist[j].append(d1)
                    skip = False # all conditions are fullfilled - this configuration is unique
                # else:
                # print 'List of distances to atoms'
                if skip: continue #

                #print "Pore can be used", xpor #sum of distances for triple
                
                new = Structure()
                new.name = name+addatom[0]+str(i)
                new.xcart = copy.deepcopy(xcart)
                new.xcart.append(xpor)
                new.typat = copy.deepcopy(typat)
                new.typat.append(addatom[1])
                # print 'new.typat  =', new.typat
                
                #calculate sum of lengths
                new.length = 0
                new.lengthCO = 0
                new.lengthCC = 0
                new.lengthOO = 0
                new.xcartC = []
                new.xcartO = []
                for m, x1 in enumerate(new.xcart):
                    if new.typat[m] == 2: new.xcartC.append(x1)
                    if new.typat[m] == 3: new.xcartO.append(x1)

                    for x2 in new.xcart:
                        d1, d2 = image_distance(x1, x2, st.rprimd, 2)
                        new.length += d1

                for xC in new.xcartC:
                    for xO in new.xcartO:
                        d1, d2 = image_distance(xC, xO, st.rprimd, 2)
                        new.lengthCO += d1                        

                for xC1 in new.xcartC:
                    for xC2 in new.xcartC:

                        d1, d2 = image_distance(xC1, xC2, st.rprimd, 2)
                        new.lengthCC += d1                     

                for xO1 in new.xcartO:
                    for xO2 in new.xcartO:

                        d1, d2 = image_distance(xO1, xO2, st.rprimd, 2)
                        new.lengthOO += d1  

                skip = False
                n = len(new.xcart)
                
                """additional conditions to leave only unique configurations"""
                for config in tlist_new:
                    if 1:
                        nr = 0
                        for (v1, t1) in zip(new.xcart, new.typat):
                            for (v2, t2) in zip(config.xcart, config.typat):
                                if all( np.around(v1, 8) == np.around(v2, 8) ) and t1 == t2:
                                    nr+=1;
                        if nr == n:
                            print("The configurations", new.name, 'and', config.name, 'consist of the same atoms, continue')
                            skip = True
                            break

                    # print  all([ all( np.around(v1, 8) == np.around(v2, 8) ) for (v1, v2) in zip(new.xcart, config.xcart) ])
                
                    #check identity using sum of distances
                    # i_min, smaller = min_diff(config.length, [new.length], diffprec)
                    # if smaller: 
                    #     print "Configuration ", new.name, "has the same sum of lengths as", config.name
                    i_min, smaller1 = min_diff(config.lengthCO, [new.lengthCO], diffprec)
                    i_min, smaller2 = min_diff(config.lengthCC, [new.lengthCC], diffprec)
                    i_min, smaller3 = min_diff(config.lengthOO, [new.lengthOO], diffprec)
                    # print 'Compare', new.name, config.name, smaller1, smaller2, smaller3
                    if smaller1 and smaller2 and smaller3: 
                        print("\nConfiguration ", new.name, "has the same sum of C-O, C-C  and O-O lengths as", config.name)
                        print() 
                        skip = True; break



                if skip:     continue
                print('\nSum of CO lengths in :',new.name, new.lengthCC, new.lengthOO, new.lengthCO)






                tlist_new.append(new)

                i+=1

        return tlist_new


    if segtyp == "bulk_triple" or segtyp == "bulk_pairs":

        # max_dist_between_atoms = 4.8

        print("\nSearching pairs ...")
        dlist = pairs(in_calc, xcart_pores, central_atoms, prec, max_dist_between_atoms, pairtyp = 'gvol' ) 

    if segtyp == "bulk_pairs":
        write_geometry_files(dlist, in_calc, xcart_pores, segtyp, take_final_rprimd_from, configver = True, add_typat = add_typat  )


        
    if segtyp == "bulk_triple":

        max_dist_to_next_atom = 5.5    
        print("\nSearching triples ...")#, tlist
        tlist = []
        tlist = triples(('O',3), dlist, tlist, in_calc, xcart_pores, max_dist_to_next_atom)

        tlist = triples(('C',2), dlist, tlist, in_calc, xcart_pores, max_dist_to_next_atom)

        write_geometry_files(dlist, in_calc, xcart_pores, segtyp, take_final_rprimd_from, tlist = tlist, configver = True, add_typat = add_typat  )






    return dlist_coseg, xcart_pores




























def create_coseg_samples(base_name, it_b, ise_b, ver, calc, path_template, imp_size = 0.44, fine_mul = 4, gbpos = None, 
    segtyp = None, input_dlist_coseg = None, xcart1imp = None , xcart2imp = None, main_path = None ):
    """
    Create all cells for co-segregation calculations.

    Input:
    base_name - the names of seg and coseg cases will be started from this name

    input_dlist_coseg - needed to construct corresponding segregation cases (segtyp = 'segreg'), but in relaxed cell, where
    one atom is in grain volume. Must be initialy constructed using segtyp = 'coseg';

    1. it_b, ise_b, base cell with relaxed grain boundaries (self.gbpos - position of one gb) 
    ver - version to create
    and one carbon atom in the volume of one grain 
    (can be several versions). 

    2. path_template - path to geo files with template cells with correct rprimd for cell with CO
    (can be several versions).     

    imp_size - size of pores, which are being found trying to insert impurity

    """
    # id = (it_b, ise_b, verlist_b[0])

    #cl = copy.deepcopy(calc[id])
    #print "gbpos is ", calc[id].gbpos



    temp    = path_template.split('/')[-1]
    temp    = temp+'_'+segtyp

    print_and_log("Please add these strings:\n")
    #print ( ("struct_des['{0:s}'] = des('path', 'template without impurities'   )\n").format(temp)   )
    print(( ("struct_des['{0:s}'] = des('path', 'template with correct rprimd; segreg - carbon in grain interior; coseg - pure'   )\n").format(temp)   ))

    """0. remove carbon but save its position and pass it to find_pairs"""
    #please remove both atoms if cell with carbon and oxygen
    #not neccessary to remove at all?
     #not needed 

       


    """1. Insert end.xred from base cell to template cell; for all found versions"""
    insert(it_b, ise_b, path_template, temp, calc, "xred" ) # resulted cell is still template with one atom in grain iterior but correct rprimd
    #Take care of versions!!!
    #return    
    """2. Find all octapores and place in them atoms with znucl = 200"""
    #v = 1
    path2temp = add_impurity(temp, 'octa', 'all_pores', calc, imp_size, write_geo = True, only_version = 1, fine = fine_mul)
    #print path2temp 
    
    #"""2.5. Inherit rprimd from other versions""" #pores generated only for first version. than 
    #for v in 1,3,4,5:
    #ver_new = v + 1
    #inherit_icalc('r1r2r3', it_new = temp, ver_new = ver_new, id_base = path2temp, calc = calc, id_from = path_template+'/from/'+temp+'.xred.'+str(ver_new)+'.geo')

    """3. Find pairs"""
    ver_rprimd_from = ver # this can be manually changed to obtain needed version
    dlist_coseg, xcart_pores = find_pairs(base_name, segtyp, path2temp,[], xcart1imp, input_dlist_coseg, 
                                          prec = 2, gvolume_config_num = 1,
                                          gbpos = gbpos,
                                          take_final_rprimd_from = geo_folder+ struct_des[temp].sfolder+'/'+temp+'/from/'+temp+'.xred.'+str(ver_rprimd_from)+'.geo',
                                          main_path = main_path, based_on = it_b+'.'+ise_b  )


    #print "Warning! Please check that dlist_coseg are the same:", dlist_coseg
    return dlist_coseg




def wrapper_create_coseg():
    basename =  't111g'
    verlist  =   list(range(1,6))
    templatepath = gb4_geo_folder+'T1/t111gCO_template'
    fine_mul = 4 # for center defining
    imp_size = 0.44

    gbpos = 20.66 #???

    #in_calc.gbpos -= 0.3   # manually improve the value using jmol
    #create_coseg_samples(               basename, 't111gCv', '93kp9', verlist, calc, templatepath, imp_size, fine_mul, gbpos, 'grainvol' ) #to construct volume cases

    dlist_coseg  = create_coseg_samples(basename, 't111g',   '9292',  verlist, calc, templatepath, imp_size, fine_mul, gbpos, 'coseg')  #cell without carbon used for cosegragation cases
    
    create_coseg_samples(               basename, 't111gCv', '93kp9', verlist, calc, templatepath, imp_size, fine_mul, gbpos, 'segreg', dlist_coseg ) #relaxed positions of cell with impurities in grain volume used for segregation cases

    # #assert len(dlist_coseg_p) == len(dlist_coseg_C)
    # print dlist_coseg_C
    # print dlist_coseg_p
    # print len(xcart_pores_C), len(xcart_pores_p)

    # print xcart_pores_C
    # print xcart_pores_p
    return

def wrapper_create_coseg_C1():
    basename =  'c1g'
    verlist  =   list(range(1,6))
    templatepath = gb4_geo_folder+'C1/c1gCO_template' #will search in target folder
    imp_size = 0.44 # for finding pores

    fine_mul = 5 # for center defining
    gbpos = 23.64   

    create_coseg_samples(               basename, 'c1gCv', '93kp7', verlist, calc, templatepath, imp_size, fine_mul, gbpos, 'grainvol' ) #to construct volume cases


    dlist_coseg  = create_coseg_samples(basename, 'c1g',   '929',  verlist, calc, templatepath, imp_size, fine_mul, gbpos, 'coseg')  #cell without carbon used for cosegragation cases
    

    xcart1imp = calc[('c1gCv', '93kp7', 1)].end.xcart[-1] 
    create_coseg_samples(               basename, 'c1g', '929', verlist, calc, templatepath,imp_size, fine_mul, gbpos, 'segreg', dlist_coseg, xcart1imp ) #relaxed positions of cell with impurities in grain volume used for segregation cases

    # #assert len(dlist_coseg_p) == len(dlist_coseg_C)
    # print dlist_coseg_C
    # print dlist_coseg_p
    # print len(xcart_pores_C), len(xcart_pores_p)

    # print xcart_pores_C
    # print xcart_pores_p
    return


def wrapper_create_pairs_H(basename = None, add_typat = [2, 3],  path2temp = None):
    # basename =  'hs443CO_test'
    # path2temp = gb4_geo_folder+'H/hps443CO2_template/grainA_s/hps443CO2_templategrainA_s.1.in.geo' #will search in target folder
    # path2temp = gb4_geo_folder+'H/hps443C2O2_template/grainA_s/hps443C2O2_templategrainA_s.1.in.geo' #will search in target folder
    # path2temp = gb4_geo_folder+'H/hps443CO_template/grainA_s/hps443CO_templategrainA_s.1.in.geo' #will search in target folder
    

    segtyp = 'bulk_pairs'
    find_pairs(basename, segtyp, path2temp, central_atoms = [10], prec = 2, max_dist_between_atoms = 8, target_znucl  = [22, 1, 1], add_typat = add_typat )

    # segtyp = 'bulk_triple'
    # find_pairs(basename, segtyp, path2temp, central_atoms = [57], prec = 2)



    return





def wrapper_create_coseg_T2(calc):
    basename =  't21g'
    main_path = 'T2/CO/'
    verlist  =   list(range(1,6))
    templatepath = gb4_geo_folder+'T2/t21gCO_template' #will search in target folder
    imp_size = 0.44 # for finding pores

    fine_mul = 5 # for center defining
    gbpos = 15.86   

    create_coseg_samples(               basename, 't21g', '93', verlist, calc, templatepath, imp_size, fine_mul, gbpos, 'grainvol', main_path = main_path ) #to construct volume cases


    dlist_coseg  = create_coseg_samples(basename, 't21g',   '93',  verlist, calc, templatepath, imp_size, fine_mul, gbpos, 'coseg', main_path = main_path)  #cell without carbon used for cosegragation cases
    

    # xcart1imp = calc[('t21gCv', '93', 1)].end.xcart[-1] 
    create_coseg_samples(               basename, 't21g', '93', verlist, calc, templatepath,imp_size, fine_mul, gbpos, 'segreg', dlist_coseg, main_path = main_path) #relaxed positions of cell with impurities in grain volume used for segregation cases

    # #assert len(dlist_coseg_p) == len(dlist_coseg_C)
    # print dlist_coseg_C
    # print dlist_coseg_p
    # print len(xcart_pores_C), len(xcart_pores_p)

    # print xcart_pores_C
    # print xcart_pores_p
    return



def wrapper_create_coseg_T1s(calc):
    basename =  't111sg'
    main_path = 'T1/CO/shift/'
    verlist  =   list(range(1,6))
    templatepath = gb4_geo_folder+'T1/t111sgCO_template' #will search in target folder
    imp_size = 0.44 # for finding pores

    fine_mul = 5 # for center defining
    gbpos = 20.97   

    for ver in verlist:
        create_coseg_samples(basename, 't111sg', '93kp9', ver, calc, templatepath, imp_size, fine_mul, gbpos, 'grainvol', main_path = main_path ) #to construct volume cases


        dlist_coseg  = create_coseg_samples(basename, 't111sg',   '93kp9',  ver, calc, templatepath, imp_size, fine_mul, gbpos, 'coseg', main_path = main_path )  #cell without carbon used for cosegragation cases
        

        create_coseg_samples(basename, 't111sg', '93kp9', ver, calc, templatepath,imp_size, fine_mul, gbpos, 'segreg', dlist_coseg, main_path = main_path  ) #relaxed positions of cell with impurities in grain volume used for segregation cases

    return

def wrapper_create_coseg_CSL7s(calc):
    basename =  'csl71sg'
    main_path = 'CSL7/CO/'
    verlist  =   list(range(1,6))
    templatepath = gb4_geo_folder+'CSL7/csl71sg15CO_template' #will search in target folder
    imp_size = 0.44 # for finding pores

    fine_mul = 5 # for center defining
    gbpos = 17.55   
    for ver in verlist:
        create_coseg_samples(basename, 'csl71sg15', '93', ver, calc, templatepath, imp_size, fine_mul, gbpos, 'grainvol', main_path = main_path ) #to construct volume cases


        dlist_coseg  = create_coseg_samples(basename, 'csl71sg15',   '93',  ver, calc, templatepath, imp_size, fine_mul, gbpos, 'coseg', main_path = main_path )  #cell without carbon used for cosegragation cases
        

        # create_coseg_samples(basename, 'csl71sg15', '93', ver, calc, templatepath,imp_size, fine_mul, gbpos, 'segreg', dlist_coseg, main_path = main_path  ) #relaxed positions of cell with impurities in grain volume used for segregation cases

    return
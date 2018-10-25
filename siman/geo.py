
import sys, copy, itertools, math
from operator import itemgetter

import numpy as np
try:
    from tabulate import tabulate
except:
    print('geo.py:tabulate is not avail')

from siman import header
from siman.header import printlog
from siman.small_functions import red_prec
# from impurity import find_pores

# sys.path.append('/home/aksenov/Simulation_wrapper/') 
# sys.path.append('/home/aksenov/Simulation_wrapper/savelyev') 



def image_distance(x1, x2, r, order = 1, sort_flag = True, return_n_distances = False):
    """
    Calculate smallest distance and the next smallest distance between two atoms 
    correctly treating periodic boundary conditions and oblique cells.
    x1, x2 - vector[3] xcart coordinates of two atoms
    r - rprimd of cell
    order - the order of periodic images which are accounted in the calcualtion of distances between atoms.
    for cubic cells, order = 1 always provide correct result.
    For highly oblique cell you should test and find the needed value of 'order' after which results are the same.
    sort_flag (bool) - use False if you do not need sorting of distances 

    return_n_distances(bool) - returns required number of smallest distances, depending on order

    return d1, d2 - the smallest and next smallest distances between atoms

    """
    d = [] # list of distances between 1st atom and images of 2nd atom
    for i in range(-order, order+1):
        for j in range(-order, order+1):
            for k in range(-order, order+1):
                x2i = x2 + (r[0] * i + r[1] * j + r[2] * k) #determine coordinates of image of atom 2 in corresponding image cells
                d.append(   np.linalg.norm(x1 - x2i)   )
    
    if sort_flag:
        d.sort()
    #print d
    # assert d[0] == min(d)

    if return_n_distances:
        return d[0:return_n_distances]
    else:
        return d[0], d[1] # old behaviour


def scale_cell_uniformly(st, scale_region = (-4,4), n_scale_images = 7, parent_calc_name = None, ):
    """
    Scale uniformly rprimd and xcart of structure() object *st* from *scale_region[0]* to *scale_region[1]*  (%) using *n_scale_images* images.
    *parent_calc_name* is added to st.des
    Return:
    list of scaled Structure() objects
    
    TODO: Take care of vol, recip and so on - the best is to create some method st.actual() that update all information 
    """
    # print scale_region
    scales = np.linspace(scale_region[0], scale_region[1], n_scale_images)
    printlog('Scales are', scales, imp = 'y')

    # print scales
    scaled_sts = []
    for j, s in enumerate(scales):
        st_s = copy.deepcopy(st)
        for i in (0,1,2):
            st_s.rprimd[i] *= (1 + s/100.)
        # print st_s.rprimd

        st_s.xred2xcart()
        st_s.des = 'obtained from '+str(parent_calc_name)+' by uniform scaling by '+str(s)+' %'
        st_s.name = str(j+1)
        scaled_sts.append(st_s)
        # print st_s.rprimd

    # plt.plot([np.linalg.norm(st.rprimd) for st in scaled_sts])
    # plt.show()
    return scaled_sts

def scale_cell_by_matrix(st, scale_region = (-4,4), n_scale_images = 7, parent_calc_name = None, mul_matrix = None ):
    """
    Scale  rprimd and xcart of structure() object *st* from *scale_region[0]* to *scale_region[1]*  (%) using *n_scale_images* images
    and mul_matrix.
    *parent_calc_name* is added to st.des
    Return:
    list of scaled Structure() objects
    
    TODO: Take care of vol, recip and so on - the best is to create some method st.actual() that update all information 
    """
    scales = np.linspace(scale_region[0], scale_region[1], n_scale_images)

    printlog('Scales are', scales, imp = 'y')
    # print(np.asarray(st.rprimd))

    scaled_sts = []
    for j, s in enumerate(scales):
        st_s = copy.deepcopy(st)
        # print(s)
        mul_matrix_f = s*(np.asarray(mul_matrix)-np.identity(3))+np.identity(3)
        st_s.rprimd = np.dot(mul_matrix_f, st_s.rprimd)
        # st_s.rprimd = np.dot(s/100*np.asarray(mul_matrix)+np.identity(3), st_s.rprimd)
        print(mul_matrix_f)
        print(np.asarray(st_s.rprimd))
        alpha, beta, gamma = st_s.get_angles()

        print(alpha, beta, gamma)    


        st_s.xred2xcart()
        st_s.des = 'obtained from '+str(parent_calc_name)+' by scaling by '+str(s)+' % '+str(mul_matrix)
        st_s.name = str(j+1)
        scaled_sts.append(st_s)
        # print st_s.rprimd

    # plt.plot([np.linalg.norm(st.rprimd) for st in scaled_sts])
    # plt.show()
    # sys.exit()


    return scaled_sts




def find_moving_atom(st1, st2):
    """
    find moving atom

    The cells should have the same rprimd!



    return number of atom which moves between two cell
    """

    for r1, r2 in zip(st1.rprimd, st2.rprimd):
        if np.linalg.norm(r1-r2)>0.001:
            printlog('Attention! find_moving_atom(): st1 and st2 have different rprimd')

    st1 = st1.return_atoms_to_cell()
    st2 = st2.return_atoms_to_cell()

    # diffv = np.array(st1.xcart) - np.array(st2.xcart)
    # diffn = np.linalg.norm(diffv, axis = 1)
    # st1.write_poscar()
    # st2.write_poscar()

    diffn = []
    for x1, x2 in zip(st1.xcart, st2.xcart):
        d1, d2 = image_distance(x1, x2, st1.rprimd)
        diffn.append(d1)

    # print('max', max(diffn))

    return np.argmax(diffn) # number of atom moving along the path







def calc_recip_vectors(rprimd):
    #Determine reciprocal vectors 
    #physics" definition
    recip = []
    vol = np.dot( rprimd[0], np.cross(rprimd[1], rprimd[2])  ); #volume
    #print vol
    recip.append(   np.cross( rprimd[1], rprimd[2] )   )
    recip.append(   np.cross( rprimd[2], rprimd[0] )   )
    recip.append(   np.cross( rprimd[0], rprimd[1] )   )
    for i in 0,1,2:
        recip[i] =  recip[i] * 2 * math.pi / vol;
    return recip



def calc_kspacings(ngkpt, rprimd):
    """Calculate kspacing from ngkpt and rprimd (A)
        ngkpt (list of int) - k-point mesh

    """
    kspacing = []

    recip = calc_recip_vectors(rprimd)

    for i in 0, 1, 2:
        a = np.linalg.norm( recip[i] ) / ngkpt[i]
        kspacing.append(red_prec(a))

    return  kspacing





def xcart2xred(xcart, rprimd):
    """Convert from cartesian coordinates xcart to
        dimensionless reduced coordinates 
        Input: xcart - list of numpy arrays, rprimd - list of numpy arrays
        Output: xred - list of numpy arrays"""
    xred = []
    gprimd = np.asarray( np.matrix(rprimd).I.T ) #Transpose of the inverse matrix of rprimd
    #print gprimd
    for xc in xcart:
        xred.append(  np.dot( gprimd , xc)  ) #dot product
    #print xred
    return xred

def xred2xcart(xred, rprimd):
    """Convert from dimensionless reduced coordinates to
    cartesian coordinates xcart;
        Input: xred - list of numpy arrays, rprimd - list of numpy arrays
        Output: xcart - list of numpy arrays"""
    xcart = []
    #print "rprimd ", rprimd
    for xr in xred:
        #for j in 0,1,2:
        #    print xr[0] * rprimd[0][j] + xr[1] * rprimd[1][j] + xr[2] * rprimd[2][j],
        #print ""
        #print np.dot( xr, rprimd)
        xcart.append(  np.dot( xr, rprimd)  ) #dot product

    #print xred
    return xcart





def replic(structure, mul = (1,1,1), inv = 1, only_atoms = None, cut_one_cell = None, include_boundary = (1,1) ):
    """
    Replicate structure() according to: mul[i]*rprimd[i]
    
    Input:
    structure - structure() type
    mul[] - is tuple of three integer numbers
    Use from structure:
    xcart, typat, rprimd, natom, xred
    inv - 1 or -1 allows to replicate in different directions
    
    inv = 0 - cell is replicated in both directions by mul[i];  2 still gives -1 0 1 but 3 gives -2 -1 0 1 2; for 'only_matrix' may work not correctly


    only_atoms - allows to replicate only specific atoms; now 
        'only_matrix'
            Warning - st.select is not working here

    cut_one_cell - allows to cut only one cell with replicated edge atoms
    include_boundary (A) - the width of region to include additional edge atoms (bottom, up)


    Return:
    replicated structure
    """
    st = copy.deepcopy(structure)
    # print 'Structure has before replication', st.natom,'atoms' 

    if hasattr(st, 'init_numbers') and st.init_numbers:
        numbers = st.init_numbers
    else:
        numbers = range(st.natom)


    #determine maximum and minimum values before replication
    xmax = [-1000000]*3
    xmin = [+1000000]*3
    for x in st.xcart:
        for j in 0,1,2:
            if xmax[j] < x[j]: xmax[j] = x[j] 
            if xmin[j] > x[j]: xmin[j] = x[j] 
    # print 'xmax, xmin', xmax, xmin



    inv_loc = inv
    for i in 0, 1, 2:
        
        axis_mul = range(mul[i])
        
        if inv == 0: # does not work propely; mul below is also should be accounted
            axis_mul = range(-mul[i]+1, mul[i])
            print_and_log('axis_mul = ', axis_mul)
            inv_loc = 1
        
        if only_atoms == 'only_matrix':
            st.xcart += [ x + inv_loc*k*st.rprimd[i] for x, t in zip(st.xcart[:], st.typat[:]) for k in axis_mul[1:] if t == 1] # fill by axis i by blocks
            
            st.typat += [ t                  for t in st.typat[:] for k in axis_mul[1:] if t == 1]
            st.magmom += [ t                  for t in st.magmom[:] for k in axis_mul[1:] if t == 1]
        else:
            st.xcart = [ x + inv_loc*k*st.rprimd[i] for x in st.xcart[:] for k in axis_mul ] # fill by axis i by blocks
            st.typat = [ t                  for t in st.typat[:] for k in axis_mul ]
            st.magmom = [ t                  for t in st.magmom[:] for k in axis_mul ]
            st.select = [ t                  for t in st.select[:] for k in axis_mul ]
            numbers = [n for n in numbers[:] for k in axis_mul]
            # print numbers
            # assert len(st.xcart) == abs(st.natom * reduce(lambda x, y: x*y, mul) )


        # print 'before', st.rprimd[i]
        
    if cut_one_cell:
        pass        
    else:

        for i in 0, 1, 2:
            if inv == 0:
                k = 2 * mul[i] - 1
            else:
                k = mul[i]
            st.rprimd[i] = st.rprimd[i] * k
        

    st.init_numbers = numbers

        # print st.init_numbers 
        # print 'after', st.rprimd[i]
    # print len(st.xcart)
    # print st.natom * reduce(lambda x, y: x*y, mul)
    #print st.xcart, 

    st.xred = xcart2xred(st.xcart, st.rprimd)


    if cut_one_cell:
        new_xred = []
        new_xcart = []
        new_typat = []
        new_mgmom = []

        precb = include_boundary[0]#/max(st.rprimd[2])
        precu = include_boundary[1]#/max(st.rprimd[2])
        


        bob = 0 - precb; upb = 1. + precu;
        # print bob, upb
        

        n = 0 
        # print st.xred
        for t, xr, m in  zip(st.typat, st.xcart, st.magmom):
            for j in 0,1,2:
                if (xr[j]  < xmin[j] - precb) or (xr[j]  > xmax[j] + precu): break  
            else:
                new_xcart.append(xr)
                new_typat.append(t)
                new_magmom.append(m)

        st.typat = new_typat
        st.xcart = new_xcart
        st.magmom = new_magmom
        print_and_log('After removing, cell has ', len(st.xred) )
        # print st.xred
        # st.xcart = xred2xcart(st.xred, st.rprimd)
        st.xred = xcart2xred(st.xcart, st.rprimd)

    st.get_nznucl()
    st.natom = len(st.xcart)
    # print 'Structure is replicated; now', st.natom,'atoms' 
    return st



def local_surrounding(x_central, st, n_neighbours, control = 'sum', periodic = False, only_elements = None, only_numbers = None, round_flag = 1):
    """
    Return list of distances to n closest atoms around central atom. (By defauld sum of distances)
    
    Input:
    - x_central - cartesian coordinates of central atom; vector
    - st - structure with xcart list of coordinates of all atoms in system
    - n_neighbours - number of needed closest neighbours

    - control - type of output; 
              sum - sum of distances, 
              av - average distance, 
              avsq - average squared dist
              'mavm': #min, av, max, av excluding min and max
              av_dev - return (average deviation, maximum deviation) from average distance in mA.
              list - list of distances; 
              atoms  - coordinates of neighbours

    - periodic - if True, then cell is additionaly replicated; needed for small cells
    Only for control = atoms
        - *only_elements* - list of z of elements to which only the distances are needed; 
        - only_numbers  (list of int) - calc dist only to this atoms 

    round_flag (bool) - if 1 than reduce distance prec to 2 points


    #TODO:
    the periodic boundary conditions realized very stupid by replicating the cell!

    """
    # round_orig = round
    if not round_flag:
        # overwrite round function with wrapper that do nothing
        def my_round(a, b):
            return a
    else:
        my_round = round


    def av_dev(n_neighbours):
        # nonlocal n_neighbours
        n_neighbours = float(n_neighbours)
        dav = sum(dlistnn)/n_neighbours
        av_dev = sum( [abs(d-dav) for d in dlistnn] ) / n_neighbours
        max_dev = max([abs(d-dav) for d in dlistnn])
        
            
        return my_round(av_dev*1000, 0), my_round(max_dev*1000, 0)

    st_original = copy.deepcopy(st)
    st.init_numbers = None
    if periodic:
        st = replic(st, mul = (2,2,2), inv = 1 ) # to be sure that impurity is surrounded by atoms
        st = replic(st, mul = (2,2,2), inv = -1 )

    xcart = st.xcart
    typat = st.typat
    natom = st.natom
    # print x_central

    #print len(xcart)
    if only_elements:
        only_elements = list(set(only_elements))
        # print(only_elements)
        # sys.exit()


    zlist = [int(st.znucl[t-1]) for t in st.typat]
    

    dlist_unsort = [np.linalg.norm(x_central - x) for x in xcart ]# if all (x != x_central)] # list of all distances

    if only_elements:
        dlist = [np.linalg.norm(x_central - x) for x, z in zip(xcart, zlist) if z in only_elements]
    else:
        dlist = copy.deepcopy(dlist_unsort)
    dlist.sort()
    # print('local_surrounding(): dlist', dlist)


    if len(dlist) > 0 and abs(dlist[0]) < 0.01:
        dlistnn = dlist[1:n_neighbours+1] #without first impurity which is x_central
    else:
        dlistnn = dlist[:n_neighbours]

    # print('dlistnn', dlistnn)
    # os._exit(1)

    if control == 'list':
        output = dlistnn

    elif control == 'sum':
        
        output = my_round(sum(dlistnn), 2)
    
    elif control == 'av':
        n_neighbours = float(n_neighbours)
        dav = sum(dlistnn)/n_neighbours
        output = my_round(dav, 2)

    elif control == 'avsq':
        n_neighbours = float(n_neighbours)
        # print(dlistnn)
        davsq = sum([d*d for d in dlistnn])/n_neighbours
        davsq = davsq**(0.5)
        output = my_round(davsq, 2)


    elif control == 'mavm': #min, av, max
        dsort = sorted(dlistnn)
        if n_neighbours > 2:
            output = (my_round(dsort[0], 2), sum(dsort[1:-1])/(n_neighbours-2), my_round(dsort[-1], 2) ) #min, av excluding min and max, max
        else:
            output = (my_round(dsort[0], 2), 0, my_round(dsort[-1], 2) ) #min, av excluding min and max, max

       
    elif control == 'av_dev':
        output = av_dev(n_neighbours)

    elif control == 'sum_av_dev':
        output = (my_round(sum(dlistnn), 2), av_dev(n_neighbours))



    elif control == 'atoms':
        # print dlist_unsort
        if hasattr(st, 'init_numbers') and st.init_numbers:
            numbers = st.init_numbers
        else:
            numbers = range(natom)
        temp = list(zip(dlist_unsort, xcart, typat, numbers, zlist) )
        
        temp.sort(key = itemgetter(0))


        if only_elements:
            centr_type = temp[0][4]
            if centr_type in only_elements:
                first = []
            else:
                first = temp[0:1]
            temp = first+[t for t in temp if t[4] in only_elements] #including central; included ionce even if only elements are and central are the same

        if only_numbers:
            temp = temp[0:1]+[t for t in temp if t[3] in only_numbers]



        temp2 = list( zip(*temp) )
        dlist       = temp2[0][:n_neighbours+1]
        xcart_local = temp2[1][:n_neighbours+1]
        typat_local = temp2[2][:n_neighbours+1]
        numbers     = temp2[3][:n_neighbours+1]
        # print temp2[0][:n_neighbours]
        # print xcart_local[:n_neighbours]
        






        #check if atoms in output are from neighboring cells
        if 0:
            xred_local = xcart2xred(xcart_local, st_original.rprimd)
            # print 'xred_local', xred_local
            for x_l in xred_local:
                for i, x in enumerate(x_l):
                    if x > 1: 
                        x_l[i]-=1
                        # print 'returning to prim cell', x,x_l[i]
                    if x < 0: 
                        x_l[i]+=1
                        # print 'returning to prim cell', x,x_l[i]
            xcart_local = xred2xcart(xred_local, st_original.rprimd)

        # print 'Warning! local_surrounding() can return several atoms in one position due to incomplete PBC implementation; Improve please\n'

        output =  (xcart_local, typat_local, numbers, dlist )

    return output


def local_surrounding2(x_central, st, n_neighbours, control = 'sum', periodic = False, only_elements = None, only_numbers = None, round_flag = 1):
    """
    !!! Attempt to improve speed of periodic conditions!
    #control = 'atoms' could work wrong!!! check

    Return list of distances to n closest atoms around central atom. (By defauld sum of distances)
    
    Input:
    - x_central - cartesian coordinates of central atom; vector
    - st - structure with xcart list of coordinates of all atoms in system
    - n_neighbours - number of needed closest neighbours

    - control - type of output; 
              sum - sum of distances, 
              av - average distance, 
              avsq - average squared dist
              'mavm': #min, av, max, av excluding min and max
              av_dev - return (average deviation, maximum deviation) from average distance in mA.
              list - list of distances; 
              atoms  - coordinates of neighbours

    - periodic - if True, then cell is additionaly replicated; needed for small cells
    Only for control = atoms
        - *only_elements* - list of z of elements to which only the distances are needed; 
        - only_numbers  (list of int) - calc dist only to this atoms 

    round_flag (bool) - if 1 than reduce distance prec to 2 points


    #TODO:
    the periodic boundary conditions realized very stupid by replicating the cell!

    """
    # round_orig = round
    if not round_flag:
        # overwrite round function with wrapper that do nothing
        def my_round(a, b):
            return a
    else:
        my_round = round


    def av_dev(n_neighbours):
        # nonlocal n_neighbours
        n_neighbours = float(n_neighbours)
        dav = sum(dlistnn)/n_neighbours
        av_dev = sum( [abs(d-dav) for d in dlistnn] ) / n_neighbours
        max_dev = max([abs(d-dav) for d in dlistnn])
        
            
        return my_round(av_dev*1000, 0), my_round(max_dev*1000, 0)

    st_original = copy.deepcopy(st)
    st.init_numbers = None
    if periodic:
        ''
        # not needed anymore, since image_distance is used,
        # however for 'atoms' regime more actions can be needed
        # st = replic(st, mul = (2,2,2), inv = 1 ) # to be sure that impurity is surrounded by atoms
        # st = replic(st, mul = (2,2,2), inv = -1 )

    xcart = st.xcart
    typat = st.typat
    natom = st.natom
    # print x_central

    #print len(xcart)
    if only_elements:
        only_elements = list(set(only_elements))
        # print(only_elements)
        # sys.exit()


    zlist = [int(st.znucl[t-1]) for t in st.typat]
    

    dlist_unsort = [image_distance(x_central, x, st.rprimd)[0] for x in xcart ]# if all (x != x_central)] # list of all distances

    if only_elements:
        dlist = [image_distance(x_central, x, st.rprimd)[0]  for x, z in zip(xcart, zlist) if z in only_elements]
    else:
        dlist = copy.deepcopy(dlist_unsort)
    dlist.sort()
    # print('local_surrounding(): dlist', dlist)


    if len(dlist) > 0 and abs(dlist[0]) < 0.01:
        dlistnn = dlist[1:n_neighbours+1] #without first impurity which is x_central
    else:
        dlistnn = dlist[:n_neighbours]

    # print('dlistnn', dlistnn)
    # os._exit(1)

    if control == 'list':
        output = dlistnn

    elif control == 'sum':
        
        output = my_round(sum(dlistnn), 2)
    
    elif control == 'av':
        n_neighbours = float(n_neighbours)
        dav = sum(dlistnn)/n_neighbours
        output = my_round(dav, 2)

    elif control == 'avsq':
        n_neighbours = float(n_neighbours)
        # print(dlistnn)
        davsq = sum([d*d for d in dlistnn])/n_neighbours
        davsq = davsq**(0.5)
        output = my_round(davsq, 2)


    elif control == 'mavm': #min, av, max
        dsort = sorted(dlistnn)
        if n_neighbours > 2:
            output = (my_round(dsort[0], 2), sum(dsort[1:-1])/(n_neighbours-2), my_round(dsort[-1], 2) ) #min, av excluding min and max, max
        else:
            output = (my_round(dsort[0], 2), 0, my_round(dsort[-1], 2) ) #min, av excluding min and max, max

       
    elif control == 'av_dev':
        output = av_dev(n_neighbours)

    elif control == 'sum_av_dev':
        output = (my_round(sum(dlistnn), 2), av_dev(n_neighbours))



    elif control == 'atoms':
        # print dlist_unsort
        if hasattr(st, 'init_numbers') and st.init_numbers:
            numbers = st.init_numbers
        else:
            numbers = range(natom)
        temp = list(zip(dlist_unsort, xcart, typat, numbers, zlist) )
        
        temp.sort(key = itemgetter(0))


        if only_elements:
            centr_type = temp[0][4]
            if centr_type in only_elements:
                first = []
            else:
                first = temp[0:1]
            temp = first+[t for t in temp if t[4] in only_elements] #including central; included ionce even if only elements are and central are the same

        if only_numbers:
            temp = temp[0:1]+[t for t in temp if t[3] in only_numbers]



        temp2 = list( zip(*temp) )
        dlist       = temp2[0][:n_neighbours+1]
        xcart_local = temp2[1][:n_neighbours+1]
        typat_local = temp2[2][:n_neighbours+1]
        numbers     = temp2[3][:n_neighbours+1]
        # print temp2[0][:n_neighbours]
        # print xcart_local[:n_neighbours]
        






        #check if atoms in output are from neighboring cells
        if 0:
            xred_local = xcart2xred(xcart_local, st_original.rprimd)
            # print 'xred_local', xred_local
            for x_l in xred_local:
                for i, x in enumerate(x_l):
                    if x > 1: 
                        x_l[i]-=1
                        # print 'returning to prim cell', x,x_l[i]
                    if x < 0: 
                        x_l[i]+=1
                        # print 'returning to prim cell', x,x_l[i]
            xcart_local = xred2xcart(xred_local, st_original.rprimd)

        # print 'Warning! local_surrounding() can return several atoms in one position due to incomplete PBC implementation; Improve please\n'

        output =  (xcart_local, typat_local, numbers, dlist )

    return output


def ortho_vec_old(rprim, ortho_sizes = None):
    """
    old function
    Function returns mul_mat - 3 vectors of integer numbers (ndarray)
    By calculating np.dot(mul_matrix, st.rprimd) you will get rprim of orthogonal supercell (actually as close as possible to it) 
    """
    from savelyev import vector_i 

    a = vector_i.Vector()
    a.vec_new_in_vec_old(vec_new = np.diag(ortho_sizes), vec_old = rprim)
    mul_matrix = np.array(a.vec_new_in_old)
    mul_matrix = mul_matrix.round(0)
    mul_matrix = mul_matrix.astype(int)

    for i in [0,1,2]:
        if mul_matrix[i][i] == 0:
            mul_matrix[i][i] = 1

    return mul_matrix

def ortho_vec(rprim, ortho_sizes = None):
    """
    Function returns mul_mat - 3 vectors of integer numbers (ndarray)
    By calculating np.dot(mul_matrix, rprim) you will get rprim of orthogonal supercell (actually as close as possible to it) 
    """

    printlog('Calculating mul_matrix for ortho:',ortho_sizes, imp = 'y',)
    printlog('rprim is;', rprim)
    vec_new = np.diag(ortho_sizes)

    # print(rprim)
    # t = rprim[1]
    # rprim[1] = rprim[0]
    # rprim[0] = t


    mul_matrix_float = np.dot( vec_new,  np.linalg.inv(rprim) )

    # ortho_test = np.dot(mul_matrix_float, rprim )

    # print(ortho_test)
    # print(mul_matrix_float)
    printlog('mul_matrix_float:\n',mul_matrix_float, imp = 'y', end = '\n')

    mul_matrix = np.array(mul_matrix_float)
    mul_matrix = mul_matrix.round(0)
    mul_matrix = mul_matrix.astype(int)


    for i in [0,1,2]:
        if mul_matrix[i][i] == 0:
            # mul_matrix[i][i] = 1
            ''

    printlog('mul_matrix:\n',mul_matrix, imp = 'y', end = '\n')


    return mul_matrix

# def mul_matrix(rprimd1, rprimd2):
#     """
#     Determines mul matrix needed to obtain rprimd2 from rprimd1
#     """




def create_supercell(st, mul_matrix, test_overlap = False, mp = 4, bound = 0.01): 
    """ 
    st (Structure) -  
    mul_matrix (3x3 ndarray of int) - for example created by *ortho_    vec()* 


    bound (float) - shift (A) allows to correctly account atoms on boundaries
    mp    (int)  include additionall atoms before cutting supecell
    test_overlap (bool) - check if atoms are overlapping -  quite slow
    """ 
    sc = st.new() 
    # st = st.return_atoms_to_cell()
    sc.name = st.name+'_supercell'
    sc.rprimd = list(np.dot(mul_matrix, st.rprimd  ))
    printlog('Old vectors (rprimd):\n',np.round(st.rprimd,1), imp = 'y', end = '\n')
    # printlog('Mul_matrix:\n',mul_matrix, imp = 'y', end = '\n')

    printlog('New vectors (rprimd) of supercell:\n',np.round(sc.rprimd,1), imp = 'y', end = '\n')
    sc.vol = np.dot( sc.rprimd[0], np.cross(sc.rprimd[1], sc.rprimd[2])  )
    st.vol = np.dot( st.rprimd[0], np.cross(st.rprimd[1], st.rprimd[2])  )
    # sc_natom_i = int(sc.vol/st.vol*st.natom) # test
    # print(st.natom)

    if hasattr(st, 'magmom'):
        if len(st.typat) != len(st.magmom):
            st.magmom = [None]*st.natom
            mag_flag = False
        else:
            mag_flag = True
    else:
        st.magmom = [None]*st.natom
        mag_flag = False        

    sc_natom = sc.vol/st.vol*st.natom # test
    printlog('The supercell should contain', sc_natom, 'atoms ... \n', imp = 'y', end = ' ')
    sc.xcart = []
    sc.typat = []
    sc.xred  = []
    sc.magmom  = []
    #find range of multiplication
    mi = np.min(mul_matrix, axis = 0)
    ma = np.max(mul_matrix, axis = 0)
    mi[mi>0] = 0  # 

    # print(mi, ma)




    # find bound values
    lengths = np.linalg.norm(sc.rprimd, axis = 1)
    bounds = bound/lengths # in reduced coordinates
    # print(bounds)
    # print(st.xcart)
    # print([range(*z) for z in zip(mi-mp, ma+mp)])
    # print(st.rprimd)
    # print(sc.rprimd)
    for uvw in itertools.product(*[range(*z) for z in zip(mi-mp, ma+mp)]): #loop over all ness uvw
        # print(uvw)
        xcart_mul = st.xcart + np.dot(uvw, st.rprimd) # coordinates of basis for each uvw
        # print(xcart_mul)
        xred_mul  = xcart2xred(xcart_mul, sc.rprimd)

        # print(len(xred_mul), len(xcart_mul), len(st.typat), len(st.magmom) )
        for xr, xc,  t, m in zip(xred_mul, xcart_mul, st.typat, st.magmom):
            # if 0<xr[0]<1 and 0<xr[1]<1 and 0<xr[2]<1:
                # print (xr)
            if all([0-b <= r < 1-b for r, b in zip(xr, bounds)]): #only that in sc.rprimd box are needed
                sc.xcart.append( xc )
                sc.xred.append ( xr )
                sc.typat.append( t  )
                sc.magmom.append(m)
    

    sc.natom = len(sc.xcart)


    if abs(sc.natom - sc_natom)>1e-5: #test 1, number of atoms
        printlog('Error! Supercell contains wrong number of atoms:', sc.natom  , 'instead of', sc_natom, 
            'try to increase *mp* of change *bound* ')

    else:
        printlog('OK', imp = 'y')
    
    if test_overlap: #test 2: overlapping of atoms
        enx = list(enumerate(sc.xcart))
        for (i1, x1), (i2, x2) in itertools.product(enx, enx):
            # print image_distance(x1, x2, sc.rprimd)[0]
            if i1 != i2 and image_distance(x1, x2, sc.rprimd)[0] < 0.1: #less than 0.1 Angstrom
                printlog('Error! Atoms in supercell are overlapping. Play with *bound*')

    sc.recip = sc.get_recip()
    sc.znucl = copy.copy(st.znucl)
    sc.ntypat = st.ntypat
    sc.nznucl = sc.get_nznucl()

    if mag_flag is False:
        sc.magmom = [None]


    return sc


def supercell(st, ortho_sizes):
    """
    wrapper
    """
    mul_matrix = ortho_vec(st.rprimd, ortho_sizes)
    return create_supercell(st, mul_matrix)

def cubic_supercell(st, ortho_sizes):
    """
    wrapper
    """
    mul_matrix = ortho_vec(st.rprimd, ortho_sizes)
    return create_supercell(st, mul_matrix)


def determine_symmetry_positions(st, element, silent = 0):
    """
    determine non-equivalent positions for atoms of type *element*

    element (str) - name of element, for example Li

    return list of lists -  atom numbers for each non-equivalent position
    """

    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


    stp = st.convert2pymatgen()

    spg = SpacegroupAnalyzer(stp)

    info = spg.get_symmetry_dataset()

    positions = {}
    for i, (el, pos) in enumerate(zip(st.get_elements(), info['equivalent_atoms'])):
        
        if el == element and pos not in positions:
            positions[pos] = []

        if el == element:
            positions[pos].append(i)

    if not silent:
        printlog('I have found ', len(positions), 'non-equivalent positions for', element, ':',positions.keys(), imp = 'y', end = '\n')
    positions_for_print = {}
    for key in positions:
        positions_for_print[key] = [p+1 for p in positions[key]]

    if not silent:
        printlog('Atom numbers: ', positions_for_print, imp = 'y')
    
    sorted_keys = sorted(list(positions.keys()))
    pos_lists = [positions[key] for key in sorted_keys ]

    return pos_lists



def remove_atoms(st, atoms_to_remove):
    """
    remove atoms either of types provided in *atoms_to_remove* or having numbers provided in *atoms_to_remove*
    st (Structure)
    atoms_to_remove (list) - list of element names or numbers

    """
    st = st.remove_atoms(atoms_to_remove)

    return st


def remove_one_atom(st, element, del_pos = None, iat = 0):
    """
    removes one atom of element type from position del_pos
    iat - number of atom inside subset
    """
    # if not del_pos:
    #     del_pos = 1
    positions = determine_symmetry_positions(st, element)
    
    if not del_pos and len(positions) > 1:
        printlog('Error! More than one symmetry position is found, please choose del position starting from 1')
    elif len(positions) == 1:
        del_pos = 1
    else:
        printlog('Position', del_pos, 'was chosen', imp = 'y')

    pos = positions[ del_pos - 1 ]
    i_del = pos[iat]
    st = st.del_atom(i_del) # remove just iat atom
    st.name += '.'+element+str(i_del)+'del'
    st.magmom = [None]
    return st, i_del

def create_deintercalated_structure(st, element, del_pos = 1):

    """
    returns deintercalated structures

    del_pos(int) - number of position starting from 1
    """
    positions = determine_symmetry_positions(st, element)
    # position_list = sorted(list(positions.keys()))
    printlog('Choose from the following list using *del_pos*:', end = '\n', imp = 'y')
    
    for i, pos in enumerate(positions):
        printlog('     ', i+1,'--->' , pos[0], end = '\n', imp = 'y')

    # pos = position_list[ del_pos - 1 ]
    pos = positions[ del_pos - 1 ]

    printlog('You have chosen position:', pos[0], imp = 'y')

    # print(st.get_elements())

    st1 = remove_atoms(st, atoms_to_remove = pos)
    st1.name += '.'+element+str(pos)+'del'
    # print(st1.get_elements())
    # sys.exit()

    return st1


def create_replaced_structure(st, el1, el2, rep_pos = 1, only_one = False):

    """
    allow to replace symmetry non-equivalent positions structures

    rep_pos(int) - number of position starting from 1
    only_one - replace only one first atom
    """
    positions = determine_symmetry_positions(st, el1)
    # position_list = sorted(list(positions.keys()))
    printlog('Choose from the following list using *del_pos*:', end = '\n', imp = 'y')
    
    for i, pos in enumerate(positions):
        printlog('     ', i+1,'--->' , pos[0], end = '\n', imp = 'y')

    # pos = position_list[ del_pos - 1 ]
    pos = positions[ rep_pos - 1 ]

    printlog('You have chosen position:', pos[0], imp = 'y')

    # print(st.get_elements())
    if only_one:
        pos = pos[0:1]


    st1 = st.replace_atoms(atoms_to_replace = pos, el_new = el2)
    st1.name += '.'+el1+str(pos)+el2+'rep'
    # print(st1.get_elements())
    # sys.exit()

    return st1








def create_antisite_defect_old(st, cation_positions = None):
    """
    exchange cation and transition metal
    st (Structure)

    cation_positions (list of numpy arrays) - reduced coordinates of deintercalated cation positions

    """

    #1. Find first alkali ion
    def find_alkali_ion(st, j_need = 1):
        # return the number of found alk ion of *j_need* occurrence 
        elements = st.get_elements_z()
        # print (elements)
        j = 0
        for i, el in enumerate(elements):
            if el in header.ALKALI_ION_ELEMENTS:
                # print (i,el)
                j+=1
                if j == j_need:
                    return i



    i_alk = find_alkali_ion(st, 3)
    x_alk = st.xcart[i_alk]


    #2. Find closest transition metal
    sur = local_surrounding(x_alk, st, n_neighbours = 1, 
        control = 'atoms', only_elements = header.TRANSITION_ELEMENTS, periodic  = True)

    i_tr = sur[2][0]
    x_tr = st.xcart[i_tr]


    
    #3. Swap atoms
    st.write_xyz(file_name = st.name+'_antisite_start')
    st = st.mov_atoms(i_alk, x_tr)
    st = st.mov_atoms(i_tr,  x_alk)

    st.write_xyz(file_name = st.name+'_antisite_final')

    printlog('Atom ',i_alk+1,'and', i_tr+1,'were swapped')
    printlog('The distance between them is ', sur[3][0])

    st.magmom = [None]

    return st



def create_antisite_defect2(st_base, st_from, cation = None, trans = None, trans_pos = 1,  mode = None):
    """
    exchange cation and transition metal
    st_base (Structure) - basic structure in which defects are created
    st_from (Structure) - structure from which the positions of *cation* are chosen;  st_from should be consistent with st_base 
    cation (str) - element, position of which is extracted from st_from and added to st_base

    trans (str) - element name transition metal for exchange
    trans_pos (int) - number of non-equiv position of trans starting from 1

    mode - 
        'add_alk' or 'a1' - add alkali cation
        'mov_trs' or 'a2' - mov trans to alkali pos
        'add_swp' or 'a3' - add alk and swap with trans
    """

    printlog('create_antisite_defect2(): mode = ', mode, imp = 'y')


    st = st_base

    cation_xred = st_from.get_element_xred(cation)
    printlog('For cation ', cation, 'reduced coordinates:',cation_xred, ' were chosen', imp = 'y')


    positions = determine_symmetry_positions(st_from, trans)
    printlog('Transition atom ', trans, 'has ', len(positions), 'non-equiv positions', imp = 'y')
    transition_numbers = positions[trans_pos -1]



    #1. Insert cation
    cation_xcart = xred2xcart([cation_xred], st.rprimd)[0]
    st_i = st.add_atoms([cation_xcart], cation)
    st_i.write_xyz(filename = st.name+'_'+cation+'_added')


    #2. Find transition metal atoms close to cation_xcart and move it here
    sur = local_surrounding(cation_xcart, st, n_neighbours = 1, 
        control = 'atoms', only_numbers = transition_numbers, periodic  = True)

    i_trans = sur[2][0]
    x_trans = sur[0][0]

    st_m = st.mov_atoms(sur[2][0], cation_xcart)
    st_m.write_xyz(filename = st.name+'_trans_moved')


    #3. Put cation to empty trans metal pos
    st_s = st_m.add_atoms([x_trans], cation)
    st_s.write_xyz(filename = st.name+'_swapped')


    if 'add_alk' in mode or 'a1' in mode:
        st = st_i
    elif 'mov_trs' in mode or 'a2' in mode:
        st = st_m
    elif 'add_swp' in mode or 'a3' in mode:
        st = st_s

    st.magmom = [None]

    return st



def create_antisite_defect3(st, el1, el2, tol = 0.1, max_sep = 4, iatom = None):
    """
    Looks for all unique antisites for el1 and el2
    
    iatom (int) - create antistes only using this atom number

    Todo
    #check that distances through  PBC are two small
    """
    # tol = 0.1 #tolerance for distinguishing antisites within one group
    # max_sep = 4 # maximum separation of antisite

    r = st.rprimd
    pos1 = determine_symmetry_positions(st, el1)
    pos2 = determine_symmetry_positions(st, el2)

    anti = {}

    for eqv_atoms1 in pos1:
        for eqv_atoms2 in pos2:
            uniq1 = eqv_atoms1[0]
            uniq2 = eqv_atoms2[0]
            lab = (uniq1, uniq2)
            if lab not in anti:
                anti[lab] = []

            for at1 in eqv_atoms1:
                for at2 in eqv_atoms2:

                    if iatom != None:
                        if at1 != iatom:
                            continue

                    x1 = st.xcart[at1]
                    x2 = st.xcart[at2]
                    d = image_distance(x1, x2, r)[0]
                    
                    if d > max_sep:
                        continue # skip larger than asked

                    for tup in anti[lab]:
                        if abs(d-tup[2]) < tol:  #antisite already included 
                            break
                    else:
                        anti[lab].append([at1, at2, round(d,3)])
                        # print(lab, at1, at2, d)
    
    structures = []
    table = []
    i = 0
    for k in anti:
        anti[k].sort(key=itemgetter(2))
        # print([k]+anti[k])
        for a in anti[k]:
            table.append([i, k]+[a[0]+1,a[1]+1,a[2]])
            # st_as = copy.deepcopy(st)
            st_as = st.swap_atoms(a[0], a[1])
            suf = 'as'+str(a[0])+'-'+str(a[1])
            # st_as.name+='_as_'+str(k)+'_with_atoms_'+str(a[0]+1)+'_and_'+str(a[1]+1)+'_swapped'
            st_as.name+='_'+suf
            st_as.magmom = None
            
            st_as.i_el1 = a[0]
            st_as.i_el2 = a[1]             

            i+=1

            structures.append(st_as)
            st_as.write_poscar()
    st.write_xyz()

    printlog('List of antisites:', imp  = 'y')

    printlog( tabulate(table, headers = ['No.', 'Antisite type', 'at1', 'at2', 'Separation, A'], tablefmt='psql'), imp = 'Y' )



    return structures



create_antisite_defect = create_antisite_defect3


def calc_k_point_mesh(rprimd, kspacing):
    """
    rprimd (list of lists 3x3 of floats) - vectors of cell (Angstroms)
    kspacing (float) - required spacing between k-points in reciprocal space (A-1); paramter KSPACING in VASP

    the provided optimal k-mesh has the smallest sum of squared deviations of kspacings

    returns k-point mesh (list of int)
    """
    N = []
    recip = calc_recip_vectors(rprimd)
    # print(recip)


    for i in 0, 1, 2:
        n = (np.linalg.norm(recip[i])) / kspacing
        N.append( math.ceil(n) )

    N_options = [ng for ng in itertools.product( *[(n-1, n, n+1) for n in N] ) ]

    errors = [  np.sum( np.square( np.array(calc_kspacings(N, rprimd) ) - kspacing ) ) for N in N_options] # sum of squared deviation from kspacing for each option
    i_min = np.argmin(errors)

    N_opt = N_options[i_min] # k-mesh with smallest error



    printlog('I recommend k-point mesh:', N_opt, 'with k-spacings:', np.array( calc_kspacings(N_opt, rprimd) ).round(2), end = '\n', imp = 'y' )
    printlog('Other options are:', end = '\n', imp = 'y' )
    printlog('{:13s} |    {:20s}'.format('Mesh', 'k-spacings'), end = '\n', imp = 'y'  )

    for ngkpt in itertools.product( *[(n-1, n, n+1) for n in N_opt] ):
        
        printlog('{:13s} |    {:26s}'.format(str(ngkpt), str(np.array(calc_kspacings(ngkpt, rprimd) ).round(2))), end = '\n', imp = 'y'  )


    return N_opt






def remove_half_based_on_symmetry(st, sg = None, info_mode = 0):
    """
    Generate all possible configurations by removing half of atoms
    sg (int) - give back structure with specific space group

    info_mode (bool) if 1 then return list of possible space groups
    
    return list of structures with sg space groups


    """
    from collections import Counter

    def spin(ls, i):
        """
        Find recursivly all possible orderings
        ls - initial list of atoms 
        i - index in ls  

        """
        for s in 1,-1:
            
            ls[i] = s
            
            if i < len(ls)-1:
            
                spin(ls, i+1)
            
            else:
                if sum(ls) == 0:
                    orderings.append(copy.deepcopy(ls) )  
        return


    structures = []
    orderings = []
    ls = [0]*st.natom
    spin(ls, 0)
    symmetries = []


    for order in orderings:
        atoms_to_remove = [i for i, s in enumerate(order) if s < 0]
        # print(atoms_to_remove)
        st_rem = st.remove_atoms(atoms_to_remove)
        nm = st_rem.sg(silent = True)[1]
        # if nm > 50:
            # print(nm)
        symmetries.append(nm)
        
        if nm == sg:
            # st_rem.jmol()
            # sc = supercell(st_rem, [14,14,14])
            # sc.jmol()
            # sc.write_poscar('xyz/POSCAR_SC2_half')
            # sc.write_cif('xyz/POSCAR_SC2_half')
            # sys.exit()
            structures.append(st_rem)

    # print(len(orderings))
    print('The following space groups were found', Counter(symmetries))
    if info_mode:
        return list(set(symmetries))

    return structures







def remove_half(st, el, sg = None, info_mode = 0):
    """
    # works only for 

    sg - required space group

    TODO
    1. Take care about matching the initial cell and supercell from primitive
    Now the manual shift is done

    2. Make full conversion from pymat structure to mine

    """

    prim = 0

    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.io.vasp.inputs import Poscar
    from calc_manage import smart_structure_read
    st_ohne_el = st.remove_atoms([el])



    st_only_el = st.leave_only(el)

    st_mp = st_only_el.convert2pymatgen()

    if prim:
        sf = SpacegroupAnalyzer(st_mp)
        st_mp_prim = sf.find_primitive() # find primitive based only on element el
    else:
        st_mp_prim = st_mp

    #convert back to my format! please improve!!!
    p = Poscar(st_mp_prim)
    p.write_file('xyz/POSCAR')
    st_prim = smart_structure_read('xyz/POSCAR')

    if info_mode:
        return remove_half_based_on_symmetry(st_prim, info_mode = 1)

    sts = remove_half_based_on_symmetry(st_prim, sg)




    st_only_el_half = sts[0]   # now only first configuration is taken, they could be different


    if prim:
        mul_matrix_float = np.dot( st.rprimd,  np.linalg.inv(st_prim.rprimd) )
        mul_matrix = np.array(mul_matrix_float)
        mul_matrix = mul_matrix.round(0)
        mul_matrix = mul_matrix.astype(int)



        sc_only_el_half = create_supercell(st_only_el_half, mul_matrix = mul_matrix)

        sc_only_el_half = sc_only_el_half.shift_atoms([0.125,0.125,0.125])
        sc_only_el_half = sc_only_el_half.return_atoms_to_cell()

    else:
        sc_only_el_half = st_only_el_half
        # sc_only_el_half


    # st_only_el.write_poscar('xyz/POSCAR1')
    # sc_only_el_half.write_poscar('xyz/POSCAR2')


    st_half = st_ohne_el.add_atoms(sc_only_el_half.xcart, el)

    st_half.name+='_half'+str(sg)
    
    return st_half




def remove_x_based_on_symmetry(st, sg = None, info_mode = 0, x = None):
    """
    Generate all possible configurations by removing half of atoms
    sg (int) - give back structure with specific space group

    info_mode (bool) if 1 then return list of possible space groups
    
    return list of structures with sg space groups


    """
    from collections import Counter

    def spin(ls, i):
        """
        Find recursivly all possible orderings
        ls - initial list of atoms 
        i - index in ls  

        """
        for s in 1,-1:
            
            ls[i] = s
            
            if i < len(ls)-1:
            
                spin(ls, i+1)
            
            else:
                if abs(ls.count(-1)/st.natom - x ) < 0.001:
                    orderings.append(copy.deepcopy(ls) )  
        return


    structures = []
    orderings = []
    ls = [0]*st.natom
    spin(ls, 0)
    symmetries = []


    for order in orderings:
        atoms_to_remove = [i for i, s in enumerate(order) if s < 0]
        # print(atoms_to_remove)
        st_rem = st.remove_atoms(atoms_to_remove)
        nm = st_rem.sg(silent = True)[1]
        # if nm > 50:
            # print(nm)
        symmetries.append(nm)
        
        if nm == sg:
            # st_rem.jmol()
            # sc = supercell(st_rem, [14,14,14])
            # sc.jmol()
            # sc.write_poscar('xyz/POSCAR_SC2_half')
            # sc.write_cif('xyz/POSCAR_SC2_half')
            # sys.exit()
            structures.append(st_rem)

    # print(len(orderings))
    print('The following space groups were found', Counter(symmetries))
    if info_mode:
        return list(set(symmetries))

    return structures



def remove_x(st, el, sg = None, info_mode = 0, x = None):
    """
    # works only for 
    x - remove x of atoms, for example 0.25 of atoms

    sg - required space group

    TODO
    1. Take care about matching the initial cell and supercell from primitive
    Now the manual shift is done

    2. Make full conversion from pymat structure to mine

    """

    prim = 0

    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.io.vasp.inputs import Poscar
    from calc_manage import smart_structure_read
    st_ohne_el = st.remove_atoms([el])



    st_only_el = st.leave_only(el)

    st_mp = st_only_el.convert2pymatgen()

    if prim:
        sf = SpacegroupAnalyzer(st_mp)
        st_mp_prim = sf.find_primitive() # find primitive based only on element el
    else:
        st_mp_prim = st_mp

    #convert back to my format! please improve!!!
    p = Poscar(st_mp_prim)
    p.write_file('xyz/POSCAR')
    st_prim = smart_structure_read('xyz/POSCAR')

    if info_mode:
        return remove_x_based_on_symmetry(st_prim, info_mode = 1, x = x)

    sts = remove_x_based_on_symmetry(st_prim, sg, x = x )




    st_only_el_half = sts[0]   # now only first configuration is taken, they could be different


    if prim:
        mul_matrix_float = np.dot( st.rprimd,  np.linalg.inv(st_prim.rprimd) )
        mul_matrix = np.array(mul_matrix_float)
        mul_matrix = mul_matrix.round(0)
        mul_matrix = mul_matrix.astype(int)



        sc_only_el_half = create_supercell(st_only_el_half, mul_matrix = mul_matrix)

        sc_only_el_half = sc_only_el_half.shift_atoms([0.125,0.125,0.125])
        sc_only_el_half = sc_only_el_half.return_atoms_to_cell()

    else:
        sc_only_el_half = st_only_el_half
        # sc_only_el_half


    # st_only_el.write_poscar('xyz/POSCAR1')
    # sc_only_el_half.write_poscar('xyz/POSCAR2')


    st_half = st_ohne_el.add_atoms(sc_only_el_half.xcart, el)

    st_half.name+='_half'+str(sg)
    
    return st_half
















def primitive(st):
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    
    st.sg()
    # st.jmol()
    st_mp = st.convert2pymatgen()
    # print(st_mp)


    sf = SpacegroupAnalyzer(st_mp)

    st_mp_prim = sf.find_primitive()

    # print(st_mp_prim)

    st  =st.update_from_pymatgen(st_mp_prim)
    # st.sg()
    return st





def create_surface(st, miller_index, min_slab_size = 10, min_vacuum_size = 10, surface_i = 0, oxidation = None, ):
    """
    INPUT:
        st (Structure) - Initial input structure. Note that to
                ensure that the miller indices correspond to usual
                crystallographic definitions, you should supply a conventional
                unit cell structure.

        miller_index ([h, k, l]): Miller index of plane parallel to
                        surface. Note that this is referenced to the input structure. If
                        you need this to be based on the conventional cell,
                        you should supply the conventional structure.


        oxidation (dic) - dictionary of effective oxidation states, e. g. {'Y':'Y3+', 'Ba':'Ba2+', 'Co':'Co2.25+', 'O':'O2-'}
                          allows to calculate dipole moment

        surface_i (int) - choose particular surface 

        min_slab_size (float) - minimum slab size

        min_vacuum_size (float) - vacuum thicknes in A

    """

    from pymatgen.core.surface import SlabGenerator
    from pymatgen.io.vasp.inputs import Poscar
    from siman.geo import replic


    pm = st.convert2pymatgen(oxidation = oxidation)
    # pm = st.convert2pymatgen()


    slabgen = SlabGenerator(pm, miller_index, min_slab_size, min_vacuum_size)
    # print(slabgen.oriented_unit_cell)
    slabs = slabgen.get_slabs()

    printlog(len(slabs), 'surfaces were generated, choose required surface using *surface_i* argument\nWriting POSCARs to xyz', imp = 'y')

    for i, slab in enumerate(slabs):
        pos = Poscar(slab)
        pos.write_file('xyz/POSCAR_suf'+str(i))

    return slabs[surface_i]



def create_surface2(st, miller_index, shift = None, min_slab_size = 10, min_vacuum_size = 10, surface_i = 0, oxidation = None, suf = '', 
    primitive = None, symmetrize = False, cut_thickness = None, return_one = False, write_poscar = 1):
    """
    INPUT:
        st (Structure) - Initial input structure. Note that to
                ensure that the miller indices correspond to usual
                crystallographic definitions, you should supply a conventional
                unit cell structure.


        pymatgen-related:
            miller_index ([h, k, l]): Miller index of plane parallel to
                            surface. Note that this is referenced to the input structure. If
                            you need this to be based on the conventional cell,
                            you should supply the conventional structure.


            oxidation (dic) - dictionary of effective oxidation states, e. g. {'Y':'Y3+', 'Ba':'Ba2+', 'Co':'Co2.25+', 'O':'O2-'}
                              allows to calculate dipole moment

            surface_i (int) - choose particular surface 

            min_slab_size (float) - minimum slab size

            min_vacuum_size (float) - vacuum thicknes in A

            symmetrize - try to make both surfaces exact


        my_paramters:
        shift (float) - shift along z 
        cut_thickness (float) - in A - allow to remove more layers from top
        return_one (bool) - allows to return only one Structure, otherwise list of pymatgen slabs is returned 
        write_poscar (bool) -self-explained
    """

    from pymatgen.core.surface import SlabGenerator
    from pymatgen.io.vasp.inputs import Poscar
    from siman.geo import replic

    if shift:
        st = st.shift_atoms([0,0,shift])

    pm = st.convert2pymatgen(oxidation = oxidation)
    # pm = st.convert2pymatgen()

    # print(min_vacuum_size)
    # sys.exit()
    slabgen = SlabGenerator(pm, miller_index, min_slab_size, min_vacuum_size, in_unit_planes= False ,  primitive = primitive )
    # print(slabgen.oriented_unit_cell)
    slabs = slabgen.get_slabs(symmetrize = symmetrize)

    printlog(len(slabs), 'surfaces were generated, choose required surface using *surface_i* argument', imp = 'y')
    st = st.update_from_pymatgen(slabs[surface_i])

    if write_poscar:
        for i, slab in enumerate(slabs):
            pos = Poscar(slab)
            # \nWriting POSCARs to xyz
            pos.write_file('xyz/POSCAR_suf'+str(i)+str(suf))

    if cut_thickness:
        return_one = True
        # print(slabs[surface_i])
        
        z = st.get_surface_pos()[1]
        # st.printme()

        st = st.del_layers(xcart_range = [z-cut_thickness+0.1, z+0.1])

        # print(st.rprimd[2])
        st.rprimd[2][2]-=cut_thickness
        # print(st.rprimd[2])

        st.update_xred()

        st.name+='cutted'
        if write_poscar:
            st.write_poscar()


    if return_one:
        return st
    else:
        return slabs

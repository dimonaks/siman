

import sys, re, copy, itertools, math
from operator import itemgetter

import itertools
flatten = itertools.chain.from_iterable

import numpy as np
try:
    from tabulate import tabulate
except:
    print('geo.py:tabulate is not avail')

from siman import header

if header.pymatgen_flag:
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.io.vasp.inputs import Poscar

from siman.header import printlog
from siman.small_functions import red_prec
from siman.functions import invert
# from impurity import find_pores

# sys.path.append('/home/aksenov/Simulation_wrapper/') 
# sys.path.append('/home/aksenov/Simulation_wrapper/savelyev') 



def image_distance(x1, x2, r, order = 1, sort_flag = True, return_n_distances = False, coord_type = 'xcart'):
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

    coord_type (str) 
        - 'xred'
        - 'xcart'


    return d1, d2 - the smallest and next smallest distances between atoms

    """
    d = [] # list of distances between 1st atom and images of 2nd atom
    
    if coord_type == 'xcart':
        def dr(i,j,k):
            return (r[0] * i + r[1] * j + r[2] * k)
    
    if coord_type == 'xred':
        a1=np.array([1,0,0])
        a2=np.array([0,1,0])
        a3=np.array([0,0,1])
        def dr(i,j,k):
            return (a1 * i + a2 * j + a3 * k)

    for i in range(-order, order+1):
        for j in range(-order, order+1):
            for k in range(-order, order+1):
                x2i = x2 +  dr(i,j,k) #determine coordinates of image of atom 2 in corresponding image cells
                d.append(   np.linalg.norm(x1 - x2i)   )
    
    if sort_flag:
        d.sort()
    #print d
    # assert d[0] == min(d)

    if return_n_distances:
        return d[0:return_n_distances]
    else:
        return d[0], d[1] # old behaviour

def image_vector(st, x1, x2,  coord_type = 'xcart'):
    """
    Calculate smallest vector between two atoms 
    correctly treating periodic boundary conditions and oblique cells.
    x1, x2 - vector[3] xcart coordinates of two atoms
    r - rprimd of cell
   
    coord_type (str) 
        - 'xred'
        - 'xcart'

    """
    # import numpy as np

    d = [] # list of distances between 1st atom and images of 2nd atom
    x2i_list = []
    r= st.rprimd
    order = 1
    # x1 = st.xcart[x1]
    # x2 = st.xcart[x2]
    if coord_type == 'xcart':
        def dr(i,j,k):
            return (r[0] * i + r[1] * j + r[2] * k)
    if coord_type == 'xred':
        a1=np.array([1,0,0])
        a2=np.array([0,1,0])
        a3=np.array([0,0,1])
        def dr(i,j,k):
            return (a1 * i + a2 * j + a3 * k)

    for i in range(-order, order+1):
        for j in range(-order, order+1):
            for k in range(-order, order+1):
                x2i = x2 +  dr(i,j,k) #determine coordinates of image of atom 2 in corresponding image cells
                d.append(   np.linalg.norm(x1 - x2i)   )
                x2i_list.append(x2i)
    
    dmin = min(d)
    pos = d.index(dmin)
    x2i = x2i_list[pos]

    vec = (x1-x2i)/2
    return  vec

def scale_cell_uniformly(st, scale_region = None, n_scale_images = 7, parent_calc_name = None, ):
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

def scale_cell_by_matrix(st, scale_region = None, n_scale_images = 7, parent_calc_name = None, mul_matrix = None ):
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
    # print(mul_matrix)
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



    return number of atom which moves between two cell (from zero)
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
        if ngkpt[i] != 0:
            a = np.linalg.norm( recip[i] ) / ngkpt[i]
            kspacing.append(red_prec(a))
        else:
            printlog('Warning! ngkpt = 0 in geo.calc_kspacings', imp='y')

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

def rms_between_structures(st1, st2):
    """
    Compare atoms position and calculated RMS between them
    Useful to see the effect of relaxation, or atoms moves during NEB
    st1 and st2 should be imposed on each other

    INPUT
        - st1 (Structure) - structure 1
        - st2 (Structure) - structure 2

    RETURN:

    """
    els1 = st1.get_elements()
    els2 = st2.get_elements()
    sumv = 0
    sumd = 0
    sumdsqr = 0
    
    for j, x1 in enumerate(st2.xcart):
        i, d, v = st1.find_closest_atom(x1)
        if len(els1) == len(els2):
            print(i, '{:.2f}'.format(d), v, els1[i], els2[j] )
        else:
            print(i, '{:.2f}'.format(d), v, els1[i] )
        sumv= sumv + v
        sumd+=d
        sumdsqr+=d**2

    print('Average deviation {:.2f} A'.format(sumd/st2.natom) )
    print('Average squared deviation {:.2f} A'.format(np.sqrt(sumdsqr/st2.natom)) )
    print('Shift of first cell relative to second cell', sumv/st2.natom, np.linalg.norm(sumv/st2.natom))

    return 

def rms_between_structures2(st1, st2, el = None):
    """
    Compare atoms position and calculated RMS between two structures with the same order of atoms
    Useful to see the effect of relaxation, or atoms moves during NEB
    st1 and st2 should be imposed on each other

    INPUT
        - st1 (Structure) - structure 1
        - st2 (Structure) - structure 2
        - el (str) - calculate only for element *el*, if None then for all elments

    RETURN:
        - d_list (list) - list of shifts for each atom
        - rms (float) - rms over all atoms

    Aksyonov
    """
    els1 = st1.get_elements()
    els2 = st2.get_elements()

    st = st1.combine([st2], add2end = 1)
    # print(st.natom, st1.natom)
    # print(st1.xcart[1], st.xcart[1], )
    # print(st2.xcart[1], st.xcart[st1.natom+1])
    # st.printme()
    if 1:
        d_list = []
        els = st1.get_elements()
        for i in range(st1.natom):
            x1 = st.xcart[i]
            x2 = st.xcart[i+st1.natom]
            d = image_distance(x1, x2, st1.rprimd, )[0] # shift
            d_list.append(d)
            print('Atom {:4d} {:2s} d = {:5.2f} A'.format(i, els[i], d))


        n_list = range(len(d_list))
        if el:
            dn_list = [[d, n] for d, ele, n in zip(d_list, els, n_list) if ele == el]

            # print(dn_list)
            # sys.exit()

            d_list = np.array(dn_list).T[0]
            n_list = np.array(dn_list).T[1].astype(int)
            print('For element', el, ':')
        else:
            print('For all elements:')

        rms = np.sqrt( np.mean(np.array(d_list)**2)  )
        av = np.mean(d_list)
        i_max = np.array(d_list).argmax()
        print('Rms = {:5.2f}, Av = {:5.2f}, Max d = {:5.2f} A for i = {:4d} '.format(rms, av, max(d_list), n_list[i_max]))

    return d_list, rms




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


    TODO:
        oxi_states are not added yet
    Return:
    replicated structure


    """
    st = copy.deepcopy(structure)
    # print 'Structure has before replication', st.natom,'atoms' 
    if not hasattr(st, 'select'):
        st.select = [None]


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
    st.oxi_state = None
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
              avharm - average harmonic - it minimal average
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
        # output = my_round(dav, 2)
        # print(dlistnn)
        output = dav

    elif control == 'avsq':
        n_neighbours = float(n_neighbours)
        # print(dlistnn)
        davsq = sum([d*d for d in dlistnn])/n_neighbours
        davsq = davsq**(0.5)
        # output = my_round(davsq, 2)
        output = davsq

    elif control == 'avharm':
        n_neighbours = float(n_neighbours)
        davharm = n_neighbours/sum([1./d for d in dlistnn])
        output = davharm


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
    #In case of small cell also works wrong with PBC. Does not take into account the several atoms should be counted more
    than once

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
        # print('only')
        # print(xcart)
        # print(zlist)
        dlist = [image_distance(x_central, x, st.rprimd)[0]  for x, z in zip(xcart, zlist) if z in only_elements]
        # for i, x, z in zip(list(range(natom)), xcart, zlist):
        #     if z in only_elements:
        #         print(i, x, z, image_distance(x_central, x, st.rprimd)[0] )

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
        print(dlistnn)
        print(n_neighbours)
        output = my_round(dav, 2)

    elif control == 'avsq':
        n_neighbours = float(n_neighbours)
        # print(dlistnn)
        davsq = sum([d*d for d in dlistnn])/n_neighbours
        davsq = davsq**(0.5)
        output = my_round(davsq, 2)


    elif control == 'avharm':
        n_neighbours = float(n_neighbours)
        davharm = n_neighbours/sum([1./d for d in dlistnn])
        output = davharm


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

def ortho_vec(rprim, ortho_sizes = None, silent = 0):
    """
    Function returns mul_mat - 3 vectors of integer numbers (ndarray)
    By calculating np.dot(mul_matrix, rprim) you will get rprim of orthogonal supercell (actually as close as possible to it) 
    """
    if not silent:
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
    if not silent:
        printlog('mul_matrix_float:\n',mul_matrix_float, imp = 'y', end = '\n')

    mul_matrix = np.array(mul_matrix_float)
    mul_matrix = mul_matrix.round(0)
    mul_matrix = mul_matrix.astype(int)


    for i in [0,1,2]:
        if mul_matrix[i][i] == 0:
            # mul_matrix[i][i] = 1
            ''
    if not silent:
        printlog('mul_matrix:\n',mul_matrix, imp = 'y', end = '\n')


    return mul_matrix


def find_mul_mat(rprimd1, rprimd2, silent = 0):
    #find mul_matrix to convert from rprimd1 to rprimd2

    mul_matrix_float = np.dot( rprimd2,  np.linalg.inv(rprimd1) )
    if not silent:
        printlog('mul_matrix_float:\n',mul_matrix_float, imp = 'y', end = '\n')

    mul_matrix = np.array(mul_matrix_float)
    mul_matrix = mul_matrix.round(0)
    mul_matrix = mul_matrix.astype(int)
    if not silent:

        printlog('mul_matrix:\n',mul_matrix, imp = 'y', end = '\n')

    return mul_matrix_float, mul_matrix


# def mul_matrix(rprimd1, rprimd2):
#     """
#     Determines mul matrix needed to obtain rprimd2 from rprimd1
#     """





def create_supercell(st, mul_matrix, test_overlap = False, mp = 4, bound = 0.01, mul = (1,1,1), silent = 0): 

    """ 
    st (Structure) -  
    mul_matrix (3x3 ndarray of int) - for example created by *ortho_    vec()* 

    mul - multiply mul matrix - allows to choose fractions of new vectors

    bound (float) - shift (A) allows to correctly account atoms on boundaries
    mp    (int)  include additionall atoms before cutting supecell
    test_overlap (bool) - check if atoms are overlapping -  quite slow
    """ 
    sc = st.new() 
    # st = st.return_atoms_to_cell()
    sc.name = st.name+'_supercell'
    # sc.rprimd = list(np.dot(mul_matrix, st.rprimd))
    # print(sc.rprimd)
    sc.rprimd = list( np.dot(mul_matrix, st.rprimd)*np.array(mul)[:, np.newaxis]  )
    
    # print(sc.rprimd)
    if not silent:
        printlog('Old vectors (rprimd):\n',np.round(st.rprimd,1), imp = 'y', end = '\n')

        # printlog('Mul_matrix:\n',mul_matrix, imp = 'y', end = '\n')

        # printlog('Mul_matrix:\n',mul_matrix, imp = 'y', end = '\n')
    if not silent:


        printlog('New vectors (rprimd) of supercell:\n',np.round(sc.rprimd,2), imp = 'y', end = '\n')
    
    # print(sc.rprimd)
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

    if not silent:
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
        if not silent:
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
        positions_for_print[key] = [p for p in positions[key]]

    if not silent:
        printlog('Atom numbers (from zero!): ', positions_for_print, imp = 'y')
    
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
    printlog('Choose from the following list using *del_pos* starting from 1:', end = '\n', imp = 'y')
    
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

def create_single_antisite(st, el1, el2, i_el1, i_el2_list = None,
    tol = 0.1, max_sep = 4, iatom = None, 
    return_with_table = False, 
    disp_AS1 = None, mag_AS1 = None, disp_AS2 = None,
    AP_on = False, i_AP = None, mag_AP = None, disp_AP = None, confs = None ):
    
    
    """
    Looks for all unique single antisites for el1 and el2
    takes into account formation of polaron and change of oxidation state


    confs (dict)
        i_el1 - choose specific atom, from 0
    """

    if i_AP is None:
        AP_on = False


    "Start determining unique positions"
    r = st.rprimd
    
    pos1 = determine_symmetry_positions(st, el1, silent = 0)
    if i_el1 and i_el1 not in list(flatten(pos1)):
        printlog('Error!', el1, 'and', i_el1, 'are not compatible.')


    # print()
    printlog('Use confs to chose Non-equivalent sites for AS; Only first atom is taken',imp = 'y')
    sts = []
    i_el1s = []
    for conf, p in enumerate(pos1):
        if confs is not None and conf not in confs:
            continue
        if i_el1:
            i = i_el1
        else:
            i = p[0] # only first atom is taken

        st_as = st.replace_atoms([i], el2, silent = 0, mode = 2)
        # st_as.jmol(r=2)

        # print(st.get_elements()[i])
        # print(st_as.get_elements()[i])
        # print(el2)
        # print( st_as.get_el_z(i)  )
        # sys.exit()

        smag_i = ''
        if mag_AS1 is not None:
            smag_i = 'm'+str(mag_AS1)
            if st_as.get_el_z(i) not in header.TRANSITION_ELEMENTS:
                printlog('Warning! Your element in antisite is ', st.get_el_name(i), ' which is not a TM'  )
            if disp_AS1 is None:
                printlog('Error! Provide disp_i')
        
        suf = 'as'+str(i)+smag_i



        if AP_on:

            'Determine possible atom candidates near AS1 for changing oxidation state'
            z2 = st_as.get_el_z(i) # e.g. Ni_Li
            el2 = invert(z2)
            el_pol = st_as.get_elements()[i_AP]
            z_pol = invert(el_pol)
            out = st_as.nn(i, n= 40,only  = [z_pol], from_one = 0, silent = 1)
            
            d1 = 'd({0}-{1}_AP), A'.format(el2, el_pol)
            tabheader = ['No of AP '+el_pol, d1]
            tab_ap = []
            kts = []
            for d, kt in zip(out['dist'], out['numbers']):
                if kt not in kts:
                    kts.append(kt)
                    # print(kt, i, st_as.distance(kt, i))
                    tab_ap.append([kt, d, ])
                # print('AP ',d, st_as.distance(kt, i) , 'has k=', kt)
            printlog('Possible positions for additional polaron:', imp = 'Y')
            printlog( tabulate(tab_ap[1:], headers = tabheader, tablefmt='psql', floatfmt=".2f"), imp = 'Y')

            if i_AP is None:
                printlog('Error! Youve chosen AP_on, Provide i_AP based on suggestions above')

            if st_as.get_el_z(i_AP) not in header.TRANSITION_ELEMENTS:
                printlog('Warning! You want to change oxidation state on ', st_as.get_el_name(i), ' which is not a TM'  )

            if mag_AP is None:
                printlog('Error! Provide mag_AP')
            if disp_AP is None:
                printlog('Error! Provide disp_AP')                

            suf+='-'+str(i_AP)+'m'+str(mag_AP)


        if mag_AS1  is None and mag_AP is None:
            st_as.magmom = [None]

        if mag_AS1 is not None:
            st_as.magmom[i] = mag_AS1

        if disp_AS1:

            if st_as.if_surface_atom(i):
                nn = 5
            else:
                nn = 6

            st_as = st_as.localize_polaron(i, disp_AS1, nn)


        if i_AP and mag_AP is not None:
            st_as.magmom[i_AP] = mag_AP


        if i_AP and disp_AP is not None:
            # av1 = st.nn(atTM,          6, only = [8], from_one = 0, silent = 1)['av(A-O,F)']

            if st_as.if_surface_atom(i_AP):
                nn = 5
            else:
                nn = 6

            st_as = st_as.localize_polaron(i_AP, disp_AP, nn = nn)


        st_as.name+='_'+suf
        sts.append(st_as)
        i_el1s.append(i)

    return sts, i_el1s




def create_antisite_defect3(st, el1, el2, i_el2_list = None,
    tol = 0.1, max_sep = 4, iatom = None, 
    return_with_table = False, 
    disp_AS1 = None, mag_AS1 = None, disp_AS2 = None,
    AP_on = False, i_AP = None, mag_AP = None, disp_AP = None, confs = None ):
    
    
    """
    Looks for all unique antisite pairs for el1 and el2
    takes into account formation of polaron and change of oxidation state


    Antisite consisits of three parts:
        AS1 - el2_el1 (e.g. Ni_Li)
        AS2 - el1_el2 (Li_Ni)
        AP - additional polaron. if AS1 changes its oxidation state (e.g. from +3 to +2 in oxide then
        additional polaron should compensate this by oxidizing from +3 to +4)

    INPUT:
        el1 - first element name from periodic table for exchange
        el2 - second element name from periodic table for exchange
        i_el2_list - use only this specific atom numbers searching AS (duplicates iatom)

        tol - tolerance for determining unique antisite configurations (A)
        max_sep - maximum separation between antisite components (A)
        iatom (int) - create antistes using this atom number, from 0
        return_with_table (bool) - in addition to structures return table with basic information

        disp_AS1 - polaronic displacement around first component (-0.2 for hole, +0.2 for electron)
            transition metal is assumed here
        mag_AS1 - magnetic moment of TM in AS1
        
        AP_on (bool) - turn on Additional polaron suggestion and creation
        i_AP - number of AP TM atom. Positions are suggested by code depending on their position relative
        to AS1 and AS2
        mag_AP - magnetic moment of AP_nn atom
        disp_AP - polaronic displacement around AP

        confs (list) - create only this configurations, others are skipped

    RETURN:
        sts (list) - list of structures
        if return_with_table:
            table (list) - see code



    Todo
    #check that distances through  PBC could be two small
    """


    def make_antisite(st, i, j, disp_AS1, mag_AS1, disp_AS2, AP_on, i_AP, disp_AP, mag_AP):
        """
        Sub-function for making antisite, 
        taking into account change of transition metal oxidation state i.e.
        sets magnetic moments and create small polarons by displacements
        

        i, j - numbers of atoms to swap; i is AS2, j is AS1

        
        AP_on - turn on AP search and creation
        i_AP - number of atom that change oxidation state
        disp_AS1, disp_AP - displacement of surrounding oxygen for creating polaron (relevant only for TM)
        mag_AS1, mag_AP - new magnetic moments; i - should be a TM to create small polaron

        """


        # print(st.get_elements()[43])
        st_as = st.swap_atoms(i, j)

        smag_j = ''
        if mag_AS1 is not None:
            smag_j = 'm'+str(mag_AS1)
            if st.get_el_z(j) not in header.TRANSITION_ELEMENTS:
                printlog('Warning! Your second element in antisite is ',j, st.get_el_name(j), ' which is not a TM'  )


            if disp_AS1 is None:
                printlog('Error! Provide disp_i')
        
        suf = 'as'+str(i)+'-'+str(j)+smag_j

        if AP_on:
           
            'Determine possible atom candidates near AS1 and AS2 for changing oxidation state'
            z1 = st_as.get_el_z(i) # e.g. Li
            el1 = st_as.get_el_name(i) # e.g. Li
            z2 = st_as.get_el_z(j) # e.g. Ni
            el2 = st_as.get_el_name(j) # e.g. Ni
            out = st_as.nn(j, n= 40,only  = [z2], from_one = 0, silent = 1)
            
            d1 = 'd({0}_{1}-{0}_AP), A'.format(el2, el1)
            d2 = 'd({1}_{0}-{0}_AP), A'.format(el2, el1)
            tabheader = ['No of AP '+el2, d1, d2 , 'd Sum, A ' ]
            tab_ap = []
            kts = []
            for d, kt in zip(out['dist'], out['numbers']):
                if kt not in kts:
                    kts.append(kt)
                    # print(kt, i, st_as.distance(kt, i))
                    tab_ap.append([kt, d, st_as.distance(kt, i), ])
                # print('AP ',d, st_as.distance(kt, i) , 'has k=', kt)
            printlog('Possible positions for additional polaron:', imp = 'Y')
            printlog( tabulate(tab_ap[1:], headers = tabheader, tablefmt='psql', floatfmt=".2f"), imp = 'Y')




            if i_AP is None:
                printlog('Error! Youve chosen AP_on, Provide i_AP based on suggestions above')


            # print(header.TRANSITION_ELEMENTS)
            # sys.exit()
            if st_as.get_el_z(i_AP) not in header.TRANSITION_ELEMENTS:
                printlog('Warning! You want to change oxidation state on ', st_as.get_el_name(i), ' which is not a TM'  )


            "section for determining parameters for AP"
            # z = st_as.get_el_z(i_AP)



            "end of section"



            if mag_AP is None:
                printlog('Error! Provide mag_AP')
            if disp_AP is None:
                printlog('Error! Provide disp_AP')                



            suf+='-'+str(i_AP)+'m'+str(mag_AP)


        st_as.i_el1 = i
        st_as.i_el2 = j            


        if mag_AS1  is None and mag_AP is None:
            st_as.magmom = [None]

        if mag_AS1:
            st_as.magmom[j] = mag_AS1

        if disp_AS1:
            st_as = st_as.localize_polaron(j, disp_AS1)

        # print(disp_AS2)
        if disp_AS2:
            st_as = st_as.localize_polaron(i, disp_AS2)        
        # sys.exit()

        if mag_AP is not None:
            st_as.magmom[i_AP] = mag_AP


        if disp_AP is not None:
            # st_as.magmom[i_AP] = mag_AP
            st_as = st_as.localize_polaron(i_AP, disp_AP)


        st_as.name+='_'+suf

        return st_as



    # if confs == None:
        # confs = []

    "Start determining unique positions"
    r = st.rprimd
    pos1 = determine_symmetry_positions(st, el1)

    if i_el2_list:
        pos2 = [i_el2_list]
    else:
        pos2 = determine_symmetry_positions(st, el2)

    anti = {}


    'Create dictionary with all possible antisite exchanges below max_sep'

    #Loop over unique positions
    for eqv_atoms1 in pos1:
        for eqv_atoms2 in pos2:
            uniq1 = eqv_atoms1[0]
            uniq2 = eqv_atoms2[0]
            lab = (uniq1, uniq2) #label
            if lab not in anti:
                anti[lab] = []

            #Loop over equivalent atoms to scan all possible distances
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
                        # print(d-tup[2], tol)
                        if abs(d-tup[2]) < tol:  #antisite already included 
                            break
                    else:
                        anti[lab].append([at1, at2, round(d,3)])
                        # print(lab, at1, at2, d)
    



    structures = []
    numbers = []
    table = []
    i = 0
    for k in anti:
        anti[k].sort(key=itemgetter(2))
        # print([k]+anti[k])
        for a in anti[k]:
            
            if confs is None or i in confs:
                st_as = make_antisite(st, i = a[0], j = a[1], disp_AS1 = disp_AS1, mag_AS1 = mag_AS1, disp_AS2 = disp_AS2,
                    AP_on = AP_on, i_AP = i_AP, disp_AP = disp_AP, mag_AP = mag_AP )
                st_as.write_poscar()
                structures.append(st_as)
                table.append([i, k]+['', a[0],  a[1], a[2]])
                numbers.append(i)
            i+=1

    st.write_xyz()


    printlog('List of antisites:', imp  = 'y')

    printlog( tabulate(table, headers = ['No.', 'Antisite type', 'it', 'at1', 'at2', 'Separation, A'], tablefmt='psql'), imp = 'Y' )


    if return_with_table:
        return structures, table, numbers
    else:
        return structures



create_antisite_defect = create_antisite_defect3


def calc_k_point_mesh(rprimd, kspacing, silent = 0):
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


    if not silent:
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

    from siman.inout import read_poscar

    st_new = st.copy()
    st_prim = read_poscar(st_new, 'xyz/POSCAR')

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




def remove_x_based_on_symmetry(st, sg = None, info_mode = 0, x = None, silent = None):
    """
    Generate all possible configurations by removing x of atoms
    
    st (Structure) - structure with only one element!

    sg (int) - give back structure with specific space group

    info_mode (bool) if 1 then return list of possible space groups (and structures)
    
    return list of structures (only first ten) with sg space groups


    """
    from collections import Counter

    def order(ls, i):
        """
        Find recursivly all possible orderings for the given x
        ls - initial list of atoms 
        i - index in ls  

        """
        for s in 1,-1:
            
            ls[i] = s
            
            if i < len(ls)-1:
            
                order(ls, i+1)
            
            else:
                if abs(ls.count(-1)/st.natom - x ) < 0.001:
                    orderings.append(copy.deepcopy(ls) )  
        return


    structures = {}
    orderings = []
    ls = [0]*st.natom
    order(ls, 0)
    symmetries = []


    for order in orderings:
        atoms_to_remove = [i for i, s in enumerate(order) if s < 0]
        # print(atoms_to_remove)
        st_rem = st.remove_atoms(atoms_to_remove)
        nm = st_rem.sg(silent = True)[1]
        # if nm > 50:
            # print(nm)
        symmetries.append(nm)
        # print(nm)
        if nm not in structures:
            structures[nm] = []
        if len(structures[nm]) < 10:
            structures[nm].append(st_rem)
        else:
            continue  
        # if nm == sg:
            # st_rem.jmol()
            # sc = supercell(st_rem, [14,14,14])
            # sc.jmol()
            # sc.write_poscar('xyz/POSCAR_SC2_half')
            # sc.write_cif('xyz/POSCAR_SC2_half')
            # sys.exit()
            # structures.append(st_rem)

    # print(len(orderings))
    print('The following space groups were found', Counter(symmetries))
    if info_mode:
        return list(set(symmetries)), structures

    return structures[sg]



def remove_x(st, el, sg = None, info_mode = 0, x = None, silent = 0, return_sts = None):
    """

    Allows to remove x of element el from the structure.
    You should know which space group you want to get.
    If you don't know the space group, first use info_mode = 1

    st (Structure) - input structure
    el (str) - element name, e.g. Li

    x - remove x of atoms, for example 0.25 of atoms
    
    info_mode (bool) - more information

    sg - number of required space group obtained with info_mode = 1

    return_sts - return all structure, otherwise only first one is returned

    TODO
    1. Take care about matching the initial cell and supercell from primitive
    Now the manual shift is done

    2. Make full conversion from pymat structure to mine

    """

    prim = 0 # find primitive based only on element el


    st_ohne_el = st.remove_atoms([el])



    st_only_el = st.leave_only(el)

    st_mp = st_only_el.convert2pymatgen()

    if prim:
        sf = SpacegroupAnalyzer(st_mp)
        st_mp_prim = sf.find_primitive() # find primitive based only on element el
    else:
        st_mp_prim = st_mp

    if 0:
        #convert back to my format! please improve!!!
        p = Poscar(st_mp_prim)
        p.write_file('xyz/POSCAR')

        from siman.inout import read_poscar
        st_new = st.copy()
        st_prim = read_poscar(st_new, 'xyz/POSCAR')
    else:
        #new way of converting
        st_prim = st_only_el.update_from_pymatgen(st_mp_prim)


    if info_mode:
        syms, sts_dic = remove_x_based_on_symmetry(st_prim, info_mode = 1, x = x)
        sts_dic2 = {}
    
    else:
        sts = remove_x_based_on_symmetry(st_prim, sg, x = x )
        syms = [sg]
        sts_dic = {}
        sts_dic2 = {}
        sts_dic[sg] = sts
    # st_prim.jmol()
    # print(sts)


    sts_dic_one = {} # only first st for each sg
    for sg in syms:
        sts = sts_dic[sg]
        if len(sts) == 0:
            printlog('Warning! number of structures for sg',sg,'is zero')

        sts2 = []

        if not return_sts:
            sts = sts[0:1]
        for st_only_el_x in sts:
            # st_only_el_x = sts[0]   # now only first configuration is taken, they could be different


            if prim:
                mul_matrix_float = np.dot( st.rprimd,  np.linalg.inv(st_prim.rprimd) )
                mul_matrix = np.array(mul_matrix_float)
                mul_matrix = mul_matrix.round(0)
                mul_matrix = mul_matrix.astype(int)
                sc_only_el_half = create_supercell(st_only_el_x, mul_matrix = mul_matrix)
                sc_only_el_half = sc_only_el_half.shift_atoms([0.125,0.125,0.125])
                sc_only_el_half = sc_only_el_half.return_atoms_to_cell()
            else:
                sc_only_el_x = st_only_el_x


            # st_only_el.write_poscar('xyz/POSCAR1')
            # sc_only_el_half.write_poscar('xyz/POSCAR2')


            st_x = st_ohne_el.add_atoms(sc_only_el_x.xcart, el)

            st_x.name+='_'+str(x)+'_'+str(sg)
            # st_x.write_poscar()
            sts2.append(st_x)
        
        sts_dic_one[sg] = sts2[0]
        sts_dic2[sg] = sts2
    if info_mode:
        return syms, sts_dic_one
    elif return_sts:
        return sts_dic2[sg]
    else:
        return sts_dic_one[sg] # only one structure is returned





def replace_x_based_on_symmetry(st, el1, el2, x = None, sg = None, info_mode = 0, silent  = 0, mag = 0.6, mode = 'rep'):
    """
    Generate all possible configurations by replacing x of element el1 by el2 from the structure.
    You should know which space group you want to get.
    If you don't know the space group, first use info_mode = 1

    st (Structure) - input structure
    el1 (str) - element name to replace, e.g. Li
    el2 (str) - replace by
    mag (float) - magnetic moment of new element
    x - replace x of atoms, for example 0.25 of atoms

    info_mode (bool) - print all possible configurations

    mode 
        - 'rep' - replace atoms 
        - 'pol' - create polarons

    sg - number of required space group obtained with info_mode = 1
    return list of structures with sg space groups


    """

    from collections import Counter

    def order(ls, i):
        """
        Find recursivly all possible orderings for the given x
        ls - initial list of atoms 
        i - index in ls  
        
        """
        for s in 1,-1:
            
            ls[i] = s
            
            if (i < len(ls)-1):
            
                order(ls, i+1)
            
            else:
                if abs(ls.count(-1)/ntot - x ) < 0.001:
                    orderings.append(copy.deepcopy(ls) )  
        return


    if silent:
        warn = 'n'
    else:
        warn = 'Y'

    structures = []
    orderings = []
    
    req = st.get_specific_elements([invert(el1)])
    # print(req)
    # sys.exit()
    ntot = len(req)
    ls = [0]*ntot
    # print(ls)
    # sys.exit()
    if x is None:
        printlog('Error! Please provide x' )

    order(ls, 0)
    symmetries = []
    if not silent:
    
        print('Total number of orderings is', len(orderings))
    at_replace = []
    els = st.get_elements()
    for order in orderings:
        atoms_to_replace = [req[i] for i, s in enumerate(order) if s < 0]
        printlog('Atoms to replace:', list(els[i] for i in atoms_to_replace), atoms_to_replace, imp = warn)
        
        if 'rep' in mode:
            # print(st.magmom)
            st_rep = st.replace_atoms(atoms_to_replace, el2, silent = silent, mag_new = mag, mode = 1)
        if 'pol' in mode:
            # print(mag)
            # sys.exit()
            st_rep = st.make_polarons(atoms_to_replace, silent = silent, mag = mag)


        # printlog('magmom:', st_rep.magmom, imp = warn)

        nm = st_rep.sg(silent = silent)[1]
        symmetries.append(nm)
        if nm == sg:
            structures.append(st_rep)
            at_replace.append(atoms_to_replace)
    
    if not silent:
        print('The following space groups were found', Counter(symmetries))
    if info_mode:
        return list(set(symmetries))

    return structures, at_replace












def two_cell_to_one(st1, st2):
    """Join two cells 
    st1 - first cell 
    st2 - second cell
    """
    # xcart = []
    # sorts = []

    n_at = st1.natom + st2.natom

    #print(n_at, dir(st1), st2.typat)
    els2 = st2.get_elements()
    for i in range(0, st2.natom):
       st1 = st1.add_atom( xc = st2.xcart[i], element = els2[i])


    return st1




















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


def stoichiometry_criteria(st1,st2):

    natom1 = st1.get_natom()
    natom2 = st2.get_natom()

    tra1 = st1.get_transition_elements()
    tra2 = st2.get_transition_elements()
    ntra1 = len(tra1)
    if ntra1 == 0: 
        ntra1 = natom1
    ntra2 = len(tra2)
    if ntra2 == 0: 
        ntra2 = natom2
    rat1 = natom1/ntra1
    rat2 = natom2/ntra2
    mul = ntra1/ntra2

    if rat1 == rat2:
        return 1
    else:
        return 0

def stoichiometry_criteria2(st1,st2, silent = 1):
    atoms1 = st1.get_elements()
    atoms2 = st2.get_elements()

    from collections import Counter
    el_dict1 = Counter(atoms1)
    el_dict2 = Counter(atoms2)
    el1 = list(el_dict1.keys())[0]
    el2 = list(el_dict1.keys())[1]
    # print(el_dict1)
    # print(el_dict2)
    ratio1 = el_dict1[el1]/el_dict1[el2]
    ratio2 = el_dict2[el1]/el_dict2[el2]

    if ratio1 == ratio2:
        if not silent:
            print('Stoichiometric')
        return 1
    else:
        if not silent:
            print('Non-stoichiometric')
            print(round(ratio1,2), round(ratio2,2))
        return 0

def symmetry_criteria(st):
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    st = st.convert2pymatgen()
    sym_criteria = SpacegroupAnalyzer(st).is_laue()
    if sym_criteria == True:
        print('Symmetric')
        return 1
    else:
        print('Non-symmetric')
        return 0

def symmetry_criteria_at(st):
    from collections import Counter
    
    els = Counter(st.get_elements())
    sym_criteria = 0
    for el in els:
        suf_at1 = st.get_surface_atoms(el, surface = 0, surface_width=1.5)
        suf_at2 = st.get_surface_atoms(el, surface = 1, surface_width=1.5)
        print(el, suf_at1, suf_at2)
        if len(suf_at1) == len(suf_at2):
            sym_criteria += 0
        else:
            sym_criteria += 1

    if sym_criteria == 0:
        print('Symmetric')
        return 1
    else:
        print('Non-symmetric')
        return 0

def create_surface2(st, miller_index, shift = None, min_slab_size = 10, min_vacuum_size = 10, surface_i = 0, oxidation = None, suf = '', 
    primitive = None, symmetrize = False, cut_thickness = None, return_one = False, write_poscar = 1, lll_reduce  = 0, silent = 0 ):
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

            lll_reduce - try to find smaller basis vectors

        my_paramters:
        shift (float) - shift along z 
        cut_thickness (float) - in A - allow to remove more layers from top
        return_one (bool) - allows to return only one Structure, otherwise list of pymatgen slabs is returned 
        write_poscar (bool) -self-explained
    """

    from pymatgen.core.surface import SlabGenerator
    from pymatgen.io.vasp.inputs import Poscar
    from pymatgen.core.composition import Composition

    from siman.geo import replic

    void_param = 0
    st_bulk = copy.deepcopy(st)

    if 'void' in st.get_elements():
        print('\nAttention! Voids are found in the structure!\nChange it by Po atoms\n')
        void_n = st.get_numbers('void')
        st = st.replace_atoms(void_n, 'Po')
        oxidation = {**oxidation, 'Po':'Po2+'}
        void_param = 1



    if shift:
        st = st.shift_atoms([0,0,shift])

    pm = st.convert2pymatgen(oxidation = oxidation)
    # pm = st.convert2pymatgen()

    # print(min_vacuum_size)
    # sys.exit()
    slabgen = SlabGenerator(pm, miller_index, min_slab_size, min_vacuum_size, primitive = primitive, lll_reduce = lll_reduce )
    # print(slabgen.oriented_unit_cell)
    slabs = slabgen.get_slabs(symmetrize = symmetrize)

    for i, slab in enumerate(slabs):
        sl = st.update_from_pymatgen(slab)
        if stoichiometry_criteria(sl, st_bulk):
            stoi = 'Stoichiometric slab'
        else:
            stoi = 'Non-stoichiometric slab'
        

        # printlog( 
        #     ';Polar:',       slab.is_polar(),
        #     ';Eqvuivalent:', slab.have_equivalent_surfaces(), 
        #     ';symmetric:',   slab.is_symmetric(), imp = 'n')
        if not silent:
            cm = Composition(slab.formula)
            print(i, cm.reduced_composition, stoi)
        

    if not silent:
        printlog(len(slabs), 'surfaces were generated, choose required surface using *surface_i* argument', imp = 'y')

    if len(slabs) >= 1:
        # surface_i =0
        st = st.update_from_pymatgen(slabs[surface_i])

    if write_poscar:
        for i, slab in enumerate(slabs):
            pos = Poscar(slab)
            # \nWriting POSCARs to xyz
            pos.write_file('xyz/POSCAR_suf'+str(i)+str(suf))

    if cut_thickness:
        return_one = True
        # print(slabs[surface_i])
        
        z = st.get_surface_pos(reduced = True)[1]
        # st.printme()
        print('surface position ', z )
        red_thick = cut_thickness/np.linalg.norm(st.rprimd[2])
        st = st.del_layers(xred_range = [z-red_thick+0.001, z+0.001], )
        # st = st.del_layers(xcart_range = [z-cut_thickness+0.1, z+0.1], )

        # print(st.rprimd[2])
        st.rprimd[2][2]-=cut_thickness
        # print(st.rprimd[2])

        st.update_xred()
        print(oxidation)
        slab = st.convert2pymatgen(oxidation = oxidation, slab = 1)
        # print('Eqvuivalent after cutting:', slab.have_equivalent_surfaces() )


        st.name+='cutted'+str(cut_thickness)
        if write_poscar:
            st.write_poscar()


    

    if return_one:

        if void_param:
            print('Return voids in the structure\n')
            void_n = st.get_numbers('Po')
            st = st.replace_atoms(void_n, 'void')
            void_param = 0

        print('Final structure contains ', st.natom, 'atoms')
        return st
    else:
        if void_param:
            print('Don\'t forget replace Po atoms by voids in the chosen structure\n')
        return slabs


def interpolate(st1, st2, images, write_poscar = 0, poscar_folder = '', omit_edges = 1):
    """
    Linear interpolation between two structures.
    The number of atoms and order should be the same

    INPUT:
    images (int) - number of intermediate images
    write_poscar (int) - starting from given number

    omit_edges (bool) - first and last corresponding to st1 and st2 are omitted, 

    """


    xl = np.linspace(0, 1, images+2)
    if omit_edges:
        xl = xl[1:-1]

    # print(xl)
    # st1.printme()
    # st2.printme()
    R = st1.rprimd
    nl = range(st1.natom)
    sts = []
    for j, x in enumerate(xl):
        st_inter = copy.deepcopy(st1)
        for i, x1, x2, xc1,xc2 in zip(nl, st1.xred, st2.xred, st1.xcart, st2.xcart):
            # d1,d2 = image_distance(xc1, xc2, st1.rprimd)
            d = np.linalg.norm(xc1-xc2)
            # if d> 10:
                # print('d =',d)
                # print(x1, x2)
                # print(xc1, xc2)
                # print((1-x) * x1 + x * x2)
            for k in 0,1,2: #periodic boundary conditions
                # print('j = ',k, x1[k] - x2[k])

                if x1[k] - x2[k] > 0.5:
                    x1[k] -= 1
                    # print(x1)
                if x1[k] - x2[k] <= -0.5:
                    x1[k] += 1
                    # print(x1)


            st_inter.xred[i] = (1-x) * x1 + x * x2
        st_inter.update_xcart()
        sts.append(st_inter)
        if write_poscar:
            st_inter.name+='.'+str(j)
            st_inter.write_poscar(poscar_folder+str(write_poscar+j)+'.POSCAR')

    return sts


def rms_pos_diff(st1, st2):

    """
    Calculate rms difference of atomic positions, excluding moving atom 
    """

    atom_num = find_moving_atom(st1, st2)

    atoms = range(st1.natom)
    summa = 0
    for i, x1, x2 in zip(atoms, st1.xcart, st2.xcart):
        if i == atom_num:
            continue
        dx = image_distance(x1, x2, r = st1.rprimd)[0]
        summa+=dx**2
    rms = (summa/st1.natom)**0.5

    return rms 


def removed_atoms(st1, st2, tol = 1e-2):
    """
    This function finds voids by comparing ideal structure and structure with removed atoms
    Input: st1 - ideal, st2 - with removed atoms
    Return list with atomic numbers of removed atoms 
    """

    removed_atoms = []
    for i in range(0, st1.natom):
        n = 0
        for j in range(0,st2.natom):
            d = np.linalg.norm(st1.xred[i] - st2.xred[j])
            # if round(st1.xred[i][0],2) == round(st2.xred[j][0],2) and round(st1.xred[i][1],2) == round(st2.xred[j][1],2) and round(st1.xred[i][2],2) == round(st2.xred[j][2],2): 
            if d < tol: 
                None
            else: 
                n+=1

        if n == st2.natom:
            removed_atoms.append(i)
    print('Removed atoms: ',removed_atoms)
    return removed_atoms

def find_voids(st1, st2):
    """
    Function returns structure with voids in the position of removed atoms
    """
    removed_at = removed_atoms(st1, st2)
    st = st1.replace_atoms(removed_at, 'void')
    return st

def hex2rhombo(h,k,l):
    #https://chem.libretexts.org/Bookshelves/Inorganic_Chemistry/Supplemental_Modules_(Inorganic_Chemistry)/Crystallography/Fundamental_Crystallography/Miller_Indices#Rhombohedral_crystals
    i = -h - k
    hr = int(1/3*(-k + i + l))
    kr = int(1/3*( h - i + l))
    lr = int(1/3*(-h + k + l))
    print(hr,kr,lr)
    return hr,kr,lr

def rhombo2hex(h,k,l):
    #https://chem.libretexts.org/Bookshelves/Inorganic_Chemistry/Supplemental_Modules_(Inorganic_Chemistry)/Crystallography/Fundamental_Crystallography/Miller_Indices#Rhombohedral_crystals
    hh = k - l 
    kh = l - h 
    lh = h + k + l 
    print(hh,kh,lh)
    return hh,kh,lh



def create_ads_molecule(st, molecule = ['O'], mol_xc = [[0,0,0]], conf_i = [0], fix_layers = False, fix_xc_range = None, 
    under_atom = 0, find_args = {'distance':0.5, 'positions' : ['ontop']}):
    """
    The function uses special module AdsorbateSiteFinder  from pymatgen


    https://static-content.springer.com/esm/art%3A10.1038%2Fs41524-017-0017-z/MediaObjects/41524_2017_17_MOESM1_ESM.pdf
    @article{montoya2017high,
          title={A high-throughput framework for determining adsorption energies on solid surfaces},
          author={Montoya, Joseph H and Persson, Kristin A},
          journal={npj Computational Materials},
          volume={3},
          number={1},
          pages={14},
          year={2017},
          publisher={Nature Publishing Group}
        }

    molecule -  'H', 'CO' ...
    mol_xc - list with xcart of atoms in molecule: [[0,0,0]], [[0,0,0],[0,0,1.23]]
    return structure with adsorbed molecule on the surface
    conf_i - [0,1,2] - list of ads configuration numbers
            key 'all' means all constructed configurations
    under_atom return configuration with ads atom strongly under me and neme atoms in surface
    """



    from pymatgen import Structure, Lattice, Molecule
    from pymatgen.analysis.adsorption import AdsorbateSiteFinder
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.io.vasp.inputs import Poscar

    pm = st.convert2pymatgen()
    # pm = st
    asf_pm = AdsorbateSiteFinder(pm)
    st_ads_pack = []


    ads_sites = asf_pm.find_adsorption_sites()
    # print(dir(ads_sites))
    # print(dir(asf_pm))
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    sym_criteria = SpacegroupAnalyzer(pm).is_laue()
    # print(sym_criteria)



    adsorbate = Molecule(molecule, mol_xc)
    ads_structs = asf_pm.generate_adsorption_structures(adsorbate,repeat=None, min_lw=8.0,  find_args=find_args)
    print('\nI\'ve found ',len(ads_structs), ' configurations for ', molecule, ' on the surface\n')

    st_ads_pack = []
    st_ads_pack_under = {'me': None, 'nme': None}
    closest_neighbor = {'me': None, 'nme': None}


    if conf_i == 'all':
        conf_i = np.arange(len(ads_structs))
        print(conf_i)

    for i in conf_i:
        p = st.update_from_pymatgen(ads_structs[i])

        if fix_layers:
            p = p.fix_layers(xcart_range = fix_xc_range)

        st_ads_pack.append(p)

        # find closest atom for ads atom
        if under_atom:
            i_close, dist_close, delta_close = p.find_closest_neighbor(p.natom-1)
            neighbor_el = p.get_elements()[i_close]
            # print(p.find_closest_neighbor(p.natom-1))

            if neighbor_el in header.nme_list:
                neighbor_type = 'nme'
            else:
                neighbor_type = 'me'

            try:
                if dist_close < closest_neighbor[neighbor_type][1] and abs(delta_close[0]) < 0.25 and  abs(delta_close[1]) < 0.25:

                    closest_neighbor[neighbor_type] = [i, dist_close]
            except (KeyError, TypeError):
                closest_neighbor[neighbor_type] = [i, dist_close]




    if under_atom:
        print(closest_neighbor)

        if closest_neighbor['me']:
            i_m = closest_neighbor['me'][0]
            st_ads_pack_under['me'] = st_ads_pack[i_m]
        else:
            st_ads_pack_under['me'] = None

        if closest_neighbor['nme']:
            i_nm = closest_neighbor['nme'][0]
            st_ads_pack_under['nme'] = st_ads_pack[i_nm]
        else:
            st_ads_pack_under['nme'] = None


        # st_ads_pack_f = [st_ads_pack[i_m], st_ads_pack[i_nm]]

        return st_ads_pack_under
         
    else:
        if len(st_ads_pack) == 1:
            st_ads = st_ads_pack[0]

            return st_ads
        else:
            return st_ads_pack



def best_miller(hkl):
    #find best representation of hkl
    #returns float
    if min(hkl) == 0: 
        n = abs(max(hkl))
    else:
        n = abs(min(hkl)) 
    hkl = hkl/n

    d_m = 100
    for mul in range(1,10):
        hklm = hkl*mul
        hkli = hklm.round(0).astype(int)
        d = np.linalg.norm(hkli-hklm)
        # print(d, d < d_m, hklm, hkli, )
        if d < d_m:
            d_m = d
            hkli_opt = hkli
            hklm_opt = hklm
        if d < 0.001: # the obtained multiplier is nice
            break

    return hklm_opt


def hkl2uvw(hkl, rprimd):
    #convert hkl to uvw
    # print(rprimd)
    # print('rprimd', rprimd)

    recip = calc_recip_vectors(rprimd)
    # print('recip', recip)
    ghkl = np.dot(hkl, recip) # convert to cartesian
    # print('hkl', hkl)
    # print('ghkl', ghkl)
    grprimd = np.asarray( np.matrix(rprimd).I.T ) #Transpose of the inverse matrix of rprimd
    uvw = np.dot(grprimd, ghkl )
    # print('uvw', uvw)
    m = np.linalg.norm(uvw)
    uvw = uvw/m # normalize
    # uvw = uvw.round(0).astype(int)
    uvwo = best_miller(uvw)


    return uvwo

def uvw2hkl(uvw, rprimd):
    #convert,
    #tested vice versa
    recip = calc_recip_vectors(rprimd)
    ruvw = np.dot(uvw, rprimd) # convert to cartesian
    grecip = np.asarray( np.matrix(recip).I.T ) #Transpose of the inverse matrix of rprimd
    hkl = np.dot(grecip, ruvw )
    # .round(0).astype(int)
    # print(hkl)
    m = np.linalg.norm(hkl)
    hkl = hkl/m # normalaze
    # print(hkl)
    hklo = best_miller(hkl)

    return hklo


def transform_miller(rprimd1, rprimd2, hkl, silent = 1):
    """
    Convert miller indicies between two choices of primitive vectors for the same lattice.
    defined in rprimd1 to miller indicies defined in rprimd2. 
    rprimd1 and rprimd2 are two primitive vectors chosen 

    Can be used for different homomorphic lattices, but first be 
    sure that they are correctly oriented in space.

    rprimd1 (list of arrays) - first set of vectors 
    rprimd2 (list of arrays) - second set of vectors
    hkl (list of int) - hkl - miller index for first set of vectors

    RETURN
    hkl2 - miller index for second set of vectors corresponding to hkl


    """

    recip1 = calc_recip_vectors(rprimd1)
    recip2 = calc_recip_vectors(rprimd2)
    grecip2 = np.asarray( np.matrix(recip2).I.T ) # transpose of the inverse matrix 

    ghkl1 = np.dot(hkl, recip1) # convert reciprocal vector to cartesian space

    hkl2 = np.dot(grecip2, ghkl1) # get miller indices of this vector for new set of vectors


    hkl2o = best_miller(hkl2)
    hkl2i = hkl2o.round(0).astype(int)

    printlog('Check my conversion of hkl2 from float to integer', hkl2o, '->' , hkl2i)

    if not silent:
        printlog('new Miller is', hkl2i, imp = 'y')

    return hkl2i



def test_transform_miller(rprimd1, rprimd2, hkl, silent = 1):
    """
    Different attempts to find correct way of index transformstion
    rprimd1 (list of arrays) - primitive vectors 1
    rprimd2 (list of arrays) - primitive vectors 2
    hkl (list of int) - miller index for rprimd1

    """




    mul_mat, _ = find_mul_mat(rprimd1, rprimd2, silent = 1)
    norm = np.linalg.norm
    # print('Lengths R1=', norm(rprimd1[0]),norm(rprimd1[1]),norm(rprimd1[2]) )
    
    test = 4

    if test == 1:
        #transform cartesian normals,  wrong
        uvw = hkl2uvw(hkl, rprimd1)
        if not silent:
            print('Converting hkl to uvw', hkl, '->',uvw)

        #transform using cartesian - correct
        ruvw = np.dot(uvw, rprimd1) # convert to cartesian

        ruvw2 = np.dot(mul_mat, ruvw)

        grprimd2 = np.asarray( np.matrix(rprimd2).I.T ) #Transpose of the inverse matrix of rprimd
        uvw2 = np.dot(grprimd2, ruvw2 )
        if not silent:
            print('Converting uvw to ruvw', uvw, '->',ruvw)
            print('Transforming ruvw to ruvw2', ruvw, '->',ruvw2)
            print('Converting ruvw2 to uvw2', ruvw2, '->',uvw2)
            
        if 0:
            #transform using direct - gives wrong result
            uvw2 = np.dot(mul_mat, uvw)
            print('Transforming uvw to uvw2', uvw, '->',uvw2)


        hkl2 = uvw2hkl(uvw2, rprimd2)
        print('Converting uvw2 to hkl2', uvw2, '->',hkl2)

    if test == 2:
        #transform miller indexes

        #transformation matrix between two lattices allow to obtain new millre indexes
        #! only  works in case, when lattices coincide with each other 

        hkl2 = np.dot(mul_mat, hkl)
        if not silent:
            print('Transforming hkl to hkl2', hkl, '->',hkl2)


    recip1 = calc_recip_vectors(rprimd1)
    recip2 = calc_recip_vectors(rprimd2)
    grecip2 = np.asarray( np.matrix(recip2).I.T ) #Transpose of the inverse matrix of rprimd


    if test == 3:
        #use g, equivalent to test2,
        #in fact it gives some new plane, as vector is transformed!!!
        ghkl1 = np.dot(hkl, recip1) # convert to cartesian

        ghkl2 = np.dot(mul_mat, ghkl1 ) # transform

        hkl2 = np.dot(grecip2, ghkl2) # convert to index

    if test == 4:
        #correct if orientation of two phases in cartesian space coincide!
        #seems equivalent to test=2
        #more clear to understand!
        ghkl1 = np.dot(hkl, recip1) # convert to cartesian

        hkl2 = np.dot(grecip2, ghkl1) # convert to index





    d_m = 100
    for mul in range(1,10):
        hkl2m = hkl2*mul
        hkl2i = hkl2m.round(0).astype(int)
        # print(hkl2m, hkl2i)
        d = np.linalg.norm(hkl2i-hkl2m)
        # print(mul, d)
        if d < d_m:
            d_m = d
            hkl2i_opt = hkl2i
    if d > 0.1 or not silent:
        printlog('Attention! Check my conversion of hkl2 from float to integer', hkl2, '->' , hkl2i_opt)


    # print(uvw2_int_opt)
    # hkl = uvw2hkl(uvw2_int_opt, rprimd2)
    if not silent:
        # print('Converting uvw to hkl', uvw2_int_opt, '->',hkl)
        printlog('new Miller is', hkl2i_opt, imp = 'y')


    return hkl2i_opt



def calc_volume(v1, v2, v3):
    return np.dot( v1, np.cross(v2, v3)  )

def triangle_area_points(v1, v2, v3):
    # if one vector is zero, then return difference of two non zero vectors
    v1v2  = v1 - v2
    v1v3  = v1 - v3
    
    if np.linalg.norm(v1) == 0:
        a = np.linalg.norm(v2-v3)
    elif np.linalg.norm(v2) == 0:
        a = np.linalg.norm(v1-v3)
    elif np.linalg.norm(v3) == 0:
        # print(v1,v2, v1-v2)
        a = np.linalg.norm(v1-v2)
    else:
        a = np.linalg.norm(np.cross(v1v2, v1v3) ) / 2


    return a


def sl_misfit(st1,st2, silent = 0):
    size1 = st1.rprimd_len()
    size2 = st2.rprimd_len()
    misfit = [(j-i)*100/j for i,j in zip(size1,size2)]
    # print('\n\nSize 1: {},{},{} A'.format(round(size1[0],2),round(size1[1],2),round(size1[2],2)))
    # print('Size 2: {},{},{} A'.format(round(size2[0],2),round(size2[1],2),round(size2[2],2)))
    if silent == 0:
        print('Misfit: {},{} % \n\n'.format(round(misfit[0],2),round(misfit[1],2)))
    return misfit

def fit2host(st_host, st_oxide):
    replic = [1,1,1]
    misf = sl_misfit(st_host,st_oxide, silent = 1)
    for m in (0,1):
        if 60 < abs(misf[m]) < 150:
            replic[m] +=1
        elif 150 < abs(misf[m]) < 250:
            replic[m] +=2
        elif 250 < abs(misf[m]) < 350:
            replic[m] +=3
        elif 350 < abs(misf[m]) < 450:
            replic[m] +=3
    st_oxide = st_oxide.replic(replic)

    return st_oxide


def hkl_slab(st, st_host, hkl, i_suf = None):

    os.remove('/home/anton/media/vasp/log_best1')
    f = open('/home/anton/media/vasp/log_best1', 'a')
    if slabs2:
        if not i_suf:
            for sl_i in range(0,len(slabs2)):
                # print(hkl)
                st2_new = st.update_from_pymatgen(slabs2[sl_i])
                misf = sl_misfit(st_host,st2_new, silent = 0)

                replic = [1,1,1]
                for m in (0,1):
                    if 80 < abs(misf[m]) < 110:
                        replic[m] +=1
                st2_new = st2_new.replic(replic)
                # print(replic)
                misf = sl_misfit(st_host,st2_new, silent = 1)
                print(hkl, sl_i)
                string = str(hkl) + ' ' + str(sl_i) + '  Misfit: {},{} % \n\n'.format(round(misf[0],2),round(misf[1],2))
                if abs(misf[0]) < 20 and abs(misf[1])<20:
                    f.write(string)
    # else:
    f.close()
    return misf, slabs2


def create_interface_solid(st_host, st_oxide, suf_host, i_suf_host = 0, 
    seek_mode = 0, seek_range = [0,2], check_shift = None, 
    hkl_lio = None, i_suf_lio = None, size = [5,5], ads_pos = None, z_shift = 1.5, lio_thick = 8):


    st1_init = st_host.copy()
    st2_init = st_oxide.copy()
    
    if suf_host:
        sc1 = st1_init.get_conventional_cell()#.replic([2,2,1])
        slabs1 = create_surface2(sc1, suf_host, shift = None, min_slab_size = 10, min_vacuum_size = 15, 
                 surface_i = 0, oxidation = None, suf = '', 
                 symmetrize = 1, cut_thickness = 0, return_one = 0, lll_reduce = 1, primitive = 1)
        slab1 = sc1.update_from_pymatgen(slabs1[i_suf_host])
        # slab1 = slabs1
        mul_matrix = ortho_vec(slab1.rprimd, [5,5,15], silent = 1) # matrix which allows to obtain supercell close to 10x10x10 A cube 
        st1 = create_supercell(slab1, mul_matrix, silent = 1)
    else:
        st1 = st_host


    
    if seek_mode:
        for h in range(seek_range[0],seek_range[1]):
            for k in range(seek_range[0],seek_range[1]):
                for l in range(seek_range[0],seek_range[1]):
                    hkl = [h,k,l]

                    if hkl != [0,0,0]:
                        slabs2 = create_surface2(st2_init, hkl, shift = None, min_slab_size = 10, min_vacuum_size = 25, surface_i = 0, oxidation = None, suf = '', 
                        primitive = 1, symmetrize = 0, cut_thickness = None, return_one = 0, write_poscar = 0, lll_reduce  = 1)
                        for sl in range(0,len(slabs2)):
                            print(hkl, sl)
                            st2 = st2_init.update_from_pymatgen(slabs2[sl])
                            slab_2 = fit2host(st1, st2)
                            sl_misfit(st1, slab_2, silent = 0)
    else:

        slabs2 = create_surface2(st2_init, hkl_lio, shift = 0, min_slab_size = lio_thick, min_vacuum_size = 15, surface_i = i_suf_lio, oxidation = None, suf = '', 
                primitive = 0, symmetrize = 0, cut_thickness = 0, return_one = 1, write_poscar = 0, lll_reduce  = 1)
        # st2 = st2_init.update_from_pymatgen(slabs2[i_suf_lio])
        st2 = slabs2
        # st2.jmol()
        slab_2 = fit2host(st1, st2)
        st2_init = slab_2

        from siman.impurity import make_interface
        z_max1 = 0
        z_max2 = 50
        st1 = st1.add_vacuum(2,40)
        for r in st1.xcart:
            if r[2] > z_max1:
                z_max1 = r[2]
        for r in slab_2.xcart:
            if r[2] < z_max2:
                z_max2 = r[2]
        

        if 0:
            elements1 = list(set(st1.get_elements()))
            elements2 = list(set(slab_2.get_elements()))
            # st1.jmol()
            suf_ats1 = []
            suf_ats2 = []
            for el in elements1:
                try:
                    suf_ats1.extend(st1.get_surface_atoms(el, surface = 1, surface_width=1.5))
                    print(el, suf_ats1)
                except TypeError:
                    None
            for el in elements2:
                try:
                    suf_ats2.extend(slab_2.get_surface_atoms(el, surface = 0, surface_width=1.5))
                    print(el, suf_ats2)
                except TypeError:
                    None


            for el1 in suf_ats1[0:1]:
                for el2 in suf_ats2[1:3]:
                    print('\n\nStarting configuration with {} and {} atoms\n'.format(el1, el2))
                    # xc1 = st1.xcart[23]
                    # xc2 = slab_2.xcart[8]
                    xc1 = st1.xcart[el1]
                    xc2 = slab_2.xcart[el2]
                    xc2[2] -= 1.8
                    # print(xc1, xc2)
                    mat0, pas_scaled2 = make_interface(st1, xc1, slab_2, xc2)
                    # mat0.replic([2,2,1]).jmol(r=2)
                    
                    interface_z_position = min([m[2] for m in mat0.xcart])/2 + max([m[2] for m in mat0.xcart])/2
                    print(interface_z_position)
                    av_pack = []
                    
                    for st in [st1, slab_2, mat0]:
                        av_dist = 0
                        n_atoms_i = 0
                        for i in range(0,st.natom):
                            if st == mat0: 
                                if (interface_z_position - 4) < st.xcart[i][2]< (interface_z_position+4):
                                    n_atoms_i += 1
                                    # print('ok')
                                    xxx = st.nn(i, n = 6, ndict = None, silent = 1, 
                                    from_one = 0, more_info = 0, oxi_state = 0, print_average = 0)
                                    d = 0
                                    for n in xxx['dist']:
                                        d+=n
                                    d = d/6
                                    av_dist += d
                                    # print(d, n_atoms_i)
                            else:
                                n_atoms_i = st.natom
                                xxx = st.nn(i, n = 6, ndict = None, silent = 1, 
                                from_one = 0, more_info = 0, oxi_state = 0, print_average = 0)
                                d = 0
                                for n in xxx['dist']:
                                    d+=n
                                d = d/6
                                av_dist += d

                        try:
                            av = round(av_dist/n_atoms_i,2)
                        except ZeroDivisionError:
                            av = 999
                        av_pack.append(av)
                        print('Average bond lengths is {} A'.format(av))

                    # mat0.replic([2,2,1]).jmol()
                    print(av_pack)
                    if abs(av_pack[2] - av_pack[0]) < 0.03:
                        print('\nGood interface!\n\n')


        print('\n\n\n', st1.rprimd, slab_2.rprimd, '\n\n\n',st1.get_angles(), slab_2.get_angles(),)
        if ads_pos:
            interface_list = []
            suf_ats2 = (slab_2.get_surface_atoms('OLiNa', surface = 0, surface_width=0.5))
            # print(suf_ats2)
            xc1 = ads_pos
            for at in suf_ats2[0:]:
                xc2 = slab_2.xcart[at]
                mat0, pas_scaled2 = make_interface(st1, xc1, slab_2, xc2)
                # mat0.jmol(rep=[3,3,1])
                interface_list.append(mat0)
                interface_list.append(pas_scaled2)
            return interface_list







        mat2, pas_scaled2 = make_interface(st1, [0, 4, z_max1+z_shift], slab_2, [0, 0.626, z_max2],)
       
        # st_wide = st1.add_vacuum(2,15)
        # mat3, pas_scaled3 = make_interface(st_wide, [0, 4, z_max1+10.5], slab_2, [0, 0.626, z_max2],)
        return [mat2, pas_scaled2]






def symmetry_multiply(st_ideal, st, el, ops = None, rm_ovrlp = None, name = ''):
    """
    Allows to multiply atomic positions of atoms *el* in *st* according to symmetry of *st_ideal*
    Usually st_ideal and st are commensurate crystal structures, but st has some defects, which are required to 
    multiply according the symmetry of the ideal structure.


    st_ideal (Structure) - the structure with required symmetry  
    st (Structure) - the structure to be modified
    el (str) - name of element, atoms of which should be multiplied; several elements can be given separated by space
    rm_ovrlp (float) - atoms closer than rm_ovrlp (in A) will be removed, if None - nothing is done

    ops (list ) - pymatgen type list of symmetry operations used in addition to that obtained from st_ideal

    name (str)

    TODO:

    RETURN:
        st (Structure) - the structure in which all atoms *el* are multiplied according to the given symmetry

    Author: DA

    """

    st  = st.copy()

    ops_full = st_ideal.get_symmetry_operations() # the first one is identity
    
    if ops:
        ops_full.extend(ops)
    
    sym_xr   = []

    els = el.split()

    my_op = 0
    if my_op:
        from pymatgen.core.operations import SymmOp
        op = SymmOp(np.arange(1, 17).reshape(4,4)) # just any 4x4 matrix
        my_op = op.from_rotation_and_translation([[1,0,0],[0,-1,0],[0,0,1]], [0.5, 0.5, 0.5])

        # print(my_op)
        # sys.exit()
        ops_full.append(my_op)



    for v, elem in zip(st.xred, st.get_elements() ):
        if elem in els:
            for op in ops_full:
                R = op.rotation_matrix
                T = op.translation_vector
                # print(R, T)
                # print(v)
                rv = np.dot(R,v)+T
                # print(rv)
                sym_xr.append(rv)
                # print('\n')
                st = st.add_atom(xr = rv, element = elem)

    # print(elem)
    st = st.return_atoms_to_cell()
    st.name+='_multiplied'+name
    
    if rm_ovrlp:
        # st, _, _ = st.remove_close_lying(tol = rm_ovrlp)
        st = st.remove_close_lying2(tol = rm_ovrlp)
    
    st.write_poscar()
    st.write_cif(symprec = None)

    return st 


def rotate_align_with_vector(st1, at1, at2):
    '''
    The function orients the given structure by aligning the z-axis with a vector between two atoms.
    
    cl1 (str)  - Structure() object
    at1 (int)  - atomic number since 0
    at2 (int)  - atomic number since 0


    returns a structure oriented along a vector 

    Author: A.Boev
    
    '''
    from siman.functions import rotation_matrix_from_vectors
    from siman.small_functions import vec_l, angle, normal
    from siman.geo import image_vector

    # cl = cl1.copy()
    st = st1.copy()
    def get_vector(st, at1, at2):
        xc1=st.xcart[at1]
        xc2=st.xcart[at2]

        vector = [i-j for i,j in zip(xc1,xc2)]
        print(xc1,xc2, vector)
        return vector

    # vec1 = st.rprimd

    vec1 = st.rprimd[2]
    vec2 = get_vector(st, at1, at2)
    # vec2 = image_vector(st, st.xcart[at1], st.xcart[at2])
    # print(vec2,vec22,vec_l(vec2),vec_l(vec22))
    m=rotation_matrix_from_vectors(vec1, vec2)
    ang = angle(vec1,vec2)
    if ang > 90:
        ang = 180 - ang
    if vec2[0] < 0:
        vec2 = [-i for i in vec2]
    print(vec2,normal(vec1,vec2), ang)

    st1 = st.rotate(normal(vec1, vec2), ang)
    # st1.jmol(r=3)
    # print(st.rprimd,st1.rprimd)

    # print(dir(cl),cl.id)
    return st1


def remove_closest(self, el, nn = 6, n = 0, x = 0.0):
    """
    Removes closest lying atoms of type el  

    INPUT:
        st (Structure) - input structure 
        el (int array) - list of elements to remove
        nn (int) - number of closest atoms 
        n (int array) - number of removing atoms 
        x (float array) - relative number of removing atoms 

    RETURN:
        st (Structure) - modified structure 

    author - A. Burov

    """
    st = copy.deepcopy(self)
    
    atoms = {} 
    if (n != 0):
        natoms = sum(n)
        for idx, el_c in enumerate(el):
            if (n[idx] != 0):
                atoms_c = st.get_specific_elements(required_elements = [el_c], fmt = 'n', z_range = None, zr_range = None)
                if (len(atoms_c) == 0):
                    n[idx] = 0  
                atoms[el_c] = atoms_c
    elif (x != 0):
        natoms = 0
        n = []
        for item in atoms.items():
            n.append(int(len(item)*x))
            natoms += n[-1]
        print("Atoms of each species will be removed: {}".format(n))
        for idx, el_c in enumerate(el):    
            if (x[idx] != 0):
                atoms_c = st.get_specific_elements(required_elements = [el_c], fmt = 'n', z_range = None, zr_range = None)
                if (len(atoms_c) == 0):
                    n[idx] = 0 
                atoms[el_c] = atoms_c    
    else:
        natoms = 0

    atoms_removed = {} 
    for el_c in el:
        atoms_removed[el_c] = []
    for i in range(natoms):
        dist_min = 1e3 
        idx_min = -1  
        for el_idx, el_c in enumerate(el):
            if (n[el_idx] == 0):
                continue 
            for atom_idx in atoms[el_c]:
                dist = st.nn(atom_idx, nn, from_one = 0, silent = 1)['dist'][1:]
                dist_cur = sum(dist) / len(dist)
                if (dist_cur <= dist_min):
                    dist_min, idx_min = dist_cur, atom_idx
                    el_min = el_idx
        st = st.remove_atoms([idx_min], from_one = 0, clear_magmom  = 1)
        del atoms[el[el_min]][idx_min]  
        if (n[el_min] == 0):
            del atoms[el[el_min]]   
        for el_c in el:
            for idx_c, atom_c in enumerate(atoms[el_c]):
                if (atom_c > idx_min):
                    atoms[el_c][idx_c] -= 1
        atoms_removed[el[el_min]].append(idx_min)
        print("Atoms were removed: {} / {}".format(i+1, natoms))
        n[el_min] -= 1
    for el_c in el:
        print("For element {}, atoms with indicies {} were removed".format(el_c, atoms_removed[el_c]))
    print("The final reduced formula is {}".format(st.get_reduced_formula()))
    return st


def remove_vacuum(self, thickness = 0.0):
    """
    Removes vacuum from a structure  
    
    INPUT:
        st (Structure) - input structure 
        thickness (float) - remaining thickness of vacuum  
    
    RETURN:
        slab with vacuum specified thickness 

    author - A. Burov 

    """

    if (thickness < 0):
        raise ValueError('The thickness of remaining vacuum should not be negative')

    st = copy.deepcopy(self)
    xyz = list(map(list, zip(st.xcart)))
    z_coord = [i[-1][-1] for i in xyz]
    z_length = st.rprimd[2][2]
    max_z = max(z_coord)-min(z_coord)
    st_new = st.shift_atoms(vector_cart = [0, 0, -min(z_coord)], return2cell = 1)
    st_new = st_new.add_vacuum(vector_i = 2, thickness = -(z_length - max_z - thickness))

    return st_new


def make_neutral(self, oxidation = None, type = 'element', at_fixed = None, mode = 'equal', 
                return_oxidation = 1, silent = 1):
    """
    Makes slab with total a charge equal to 0 

    INPUT:
        st (Structure) - input structure 
        oxidation (dir integer) - list of oxidation states  
            E.g oxi_state = {"Li": 1, "La": 2, "Zr":4, "O": -2}
        type (dir integer) - assign oxidation states based on algorythm
            'element' - by chemical element, requires oxidation in format: 
                E.g oxi_state = {"Li": 1, "La": 2, "Zr":4, "O": -2}
            'position' - by position of chemical element, requires oxidation in format:
                E.g oxi_state = {"Li": 1, "La": 2, "Zr":4, "O": -2}
        at_fixed (dir string) - list of atoms with fixed oxidation states
        mode (string) - how uncompensated charge will be redistributed between unfixed atoms    
            'equal' - equally between unfixed atoms 
            'propotional' - proportionally to oxidation state
    
    RETURN:
        if (return_oxidation == True)
            returns a new structure with a neutral charge and new oxidation states
        else
            returns only a new structure 

    author - A. Burov 

    """

    st = copy.deepcopy(self)
    st = st.convert2pymatgen()

    st.add_oxidation_state_by_element(oxidation)    
    diff_chr = st.charge

    if (silent  == 1):
        print("Uncompensated charge is {}".format(st.charge))

    atoms = st.formula.split()
    at_init = {}

    for atom in atoms:
        at_type = re.sub(r'[0-9]+', '', atom) 
        at_number = re.sub(r'[^0-9]+', '', atom)
        at_init[at_type] = int(at_number)

    at_sum = 0
    if (mode == 'equal'):
        for key, item in at_init.items():
            if (diff_chr > 0):
                if (key not in at_fixed):
                    if (oxidation[key] > 0):
                        at_sum += item
            else:
                if (key not in at_fixed):
                    if (oxidation[key] < 0):
                        at_sum += item
            
        rel_chr = diff_chr / at_sum
        for key, item in oxidation.items():
            if (diff_chr > 0):
                if (key not in at_fixed):
                    if (oxidation[key] > 0):
                        oxidation[key] = item - rel_chr
            else:
                if (key not in at_fixed):
                    if (oxidation[key] < 0):
                        oxidation[key] = item - rel_chr

    elif (mode == 'propotional'):
        ox_sum = 0
        for key, item in at_init.items():
            if (diff_chr > 0):
                if (key not in at_fixed):
                    if (oxidation[key] > 0):
                        ox_sum += oxidation[key]
            else:
                if (key not in at_fixed):
                    if (oxidation[key] < 0):
                        ox_sum += oxidation[key]

        for key, item in oxidation.items():
            rel_chr = item / ox_sum * diff_chr / at_init[key]
            if (diff_chr > 0):
                if (key not in at_fixed):
                    if (oxidation[key] > 0):
                        oxidation[key] = item - rel_chr
            else:
                if (key not in at_fixed):
                    if (oxidation[key] < 0):
                        oxidation[key] = item - rel_chr 

        st.add_oxidation_state_by_element(oxidation)    
        print(st.charge)

    else:
        print("Wrong mode, check the function's description")

    if (silent == 1):
        print("New oxidation states are {}".format(oxidation))

    if (return_oxidation == 1):
        return st, oxidation
    else:
        return st


def move_edge(self, mode="bottom", tol=0.0):
    """
    Removes closest lying atoms of type el  

    INPUT:
        st (Structure) - input structure 
        mode (str) - to which edge of the cell structure shall be shifted
        tol (float) - tolerance factor, the distance from the edge. May be used to preserve 
                      the structure as a whole on one cell side.  
    RETURN:
        st (Structure) - modified structure 

    COMMENT: make for all vectors
    author - A. Burov
    """
    st = copy.deepcopy(self)
    st1 = copy.deepcopy(self)
    
    coords = st1.xred
    coords.sort(key=lambda x: x[2])
    # print(coords[0][-1], coords[-1][-1])
    coord_bot = coords[0][-1]
    coord_top = coords[-1][-1]

    if mode == "bottom":
        st = st.shift_atoms(vector_red = -coord_bot+0.01, return2cell = 1)
        # st = st.return_atoms_to_cell()
        print()
    elif mode == "top":
        st = st.shift_atoms(vector_red = -coord_top, return2cell = 1)
    else:
        raise ValueError("Unknown mode, check the description")

    return st



def find_slab_width(self, vacuum="no"):
    """
        Calculate width of the sample without vacuum.
        INPUT:
            st - input strucutre
            vacuum - width of the vacuum or crystal structure. Defaultly, the width of the crystal structure
        RETURN:
            width of the necessary structure in A 
    """
    st = copy.deepcopy(self)
    c_old = st.rprimd_len()[-1]
    
    st_no_vacuum = st.remove_vacuum()
    c_new = st_no_vacuum.rprimd_len()[-1]
    
    if (vacuum == "yes"):
        width = c_old - c_new
    else:
        width = c_new
        
    return round(width, 1)







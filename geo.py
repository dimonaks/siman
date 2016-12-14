
import sys, copy, itertools, math
import numpy as np

from header import printlog

from small_functions import red_prec
# sys.path.append('/home/aksenov/Simulation_wrapper/') 
# sys.path.append('/home/aksenov/Simulation_wrapper/savelyev') 


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

    vec_new = np.diag(ortho_sizes)

    # print(rprim)
    mul_matrix_float = np.dot( vec_new,  np.linalg.inv(rprim) )

    # ortho_test = np.dot(mul_matrix_float, rprim )

    # print(ortho_test)
    # print(mul_matrix_float)

    mul_matrix = np.array(mul_matrix_float)
    mul_matrix = mul_matrix.round(0)
    mul_matrix = mul_matrix.astype(int)

    for i in [0,1,2]:
        if mul_matrix[i][i] == 0:
            mul_matrix[i][i] = 1

    return mul_matrix

# def mul_matrix(rprimd1, rprimd2):
#     """
#     Determines mul matrix needed to obtain rprimd2 from rprimd1
#     """




def create_supercell(st, mul_matrix, test_overlap = False, mp = 2, bound = 0.01): 
    """ 
    st (Structure) -  
    mul_matrix (3x3 ndarray of int) - for example created by *ortho_    vec()* 


    bound (float) - shift (A) allows to correctly account atoms on boundaries
    mp    (int)  include additionall atoms before cutting supecell
    test_overlap (bool) - check if atoms are overlapping -  quite slow
    """ 
    sc = copy.deepcopy(st) 

    sc.rprimd = list(np.dot(mul_matrix, st.rprimd  ))
    printlog('New vectors (rprimd) of supercell:\n',np.round(sc.rprimd,1), imp = 'y', end = '\n')
    sc.vol = np.dot( sc.rprimd[0], np.cross(sc.rprimd[1], sc.rprimd[2])  )
    st.vol = np.dot( st.rprimd[0], np.cross(st.rprimd[1], st.rprimd[2])  )
    # sc_natom_i = int(sc.vol/st.vol*st.natom) # test
    sc_natom = sc.vol/st.vol*st.natom # test
    printlog('The supercell should contain', sc_natom, 'atoms ... ', imp = 'y', end = ' ')
    sc.xcart = []
    sc.typat = []
    sc.xred  = []
    #find range of multiplication
    mi = np.min(mul_matrix, axis = 0)
    ma = np.max(mul_matrix, axis = 0)
    mi[mi>0] = 0  # 

    # print(mi, ma)


    # find bound values
    lengths = np.linalg.norm(sc.rprimd, axis = 1)
    bounds = bound/lengths # in reduced coordinates


    for uvw in itertools.product(*[range(*z) for z in zip(mi-mp, ma+mp)]): #loop over all ness uvw
        # print(uvw)
        xcart_mul = st.xcart + np.dot(uvw, st.rprimd) # coordinates of basis for each uvw
        # print(xcart_mul)
        xred_mul  = xcart2xred(xcart_mul, sc.rprimd)
        
        for xr, xc,  t in zip(xred_mul, xcart_mul, st.typat):
            
            if all([0-b <= r < 1-b for r, b in zip(xr, bounds)]): #only that in sc.rprimd box are needed
                sc.xcart.append( xc )
                sc.xred.append ( xr )
                sc.typat.append( t  )
    

    sc.natom = len(sc.xcart)
    sc.magmom = [None]


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



    return sc





def determine_symmetry_positions(st, element):
    """
    determine non-equivalent positions for atoms of type *element*

    element (str) - name of element, for example Li

    return dictionary of atom numbers for each non-equivalent position
    """

    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


    stp = st.convert2pymatgen()

    spg = SpacegroupAnalyzer(stp)

    info = spg.get_symmetry_dataset()

    positions = {}
    for i, (el, pos) in enumerate(zip(st.get_elements(), info['equivalent_atoms'])):
        
        if el in element and pos not in positions:
            positions[pos] = []

        if el in element:
            positions[pos].append(i)


    printlog('I have found ', len(positions), 'non-equivalent positions for', element, ':',positions.keys(), imp = 'y', end = '\n')
    printlog('Atom numbers: ', positions, imp = 'y')
    return positions



def remove_atoms(st, atoms_to_remove):
    """
    remove atoms either of types provided in *atoms_to_remove* or having numbers provided in *atoms_to_remove*
    st (Structure)
    atoms_to_remove (list) - list of str or int

    """
    st = copy.deepcopy(st)
    numbers = list(range(st.natom))


    atom_exsist = True

    while atom_exsist:

        for i, (n, el) in enumerate(  zip(numbers, st.get_elements()) ):
            # print(i)
            
            if el in atoms_to_remove or n in atoms_to_remove:
                # print(n)
                # atoms_to_remove.remove(i)
                st = st.del_atom(i)
                del numbers[i]

                break
        else:
            atom_exsist = False
    printlog('remove_atoms(): Atoms', atoms_to_remove, 'were removed')

    # print(st.get_elements())
    return st


def create_deintercalated_structure(st, element, del_pos = 1):

    """
    returns deintercalated structures

    del_pos(int) - number of position starting from 1
    """
    positions = determine_symmetry_positions(st, element)
    position_list = sorted(list(positions.keys()))
    printlog('Choose from the following list using *del_pos*:', end = '\n', imp = 'y')
    
    for i, pos in enumerate(position_list):
        printlog('     ', i+1,'--->' , pos, end = '\n', imp = 'y')

    pos = position_list[ del_pos - 1 ]

    printlog('You have chosen position:', pos, imp = 'y')

    # print(st.get_elements())

    st1 = remove_atoms(st, atoms_to_remove = positions[pos])
    st1.name += '.'+element+str(pos)+'del'
    # print(st1.get_elements())
    # sys.exit()

    return st1






def create_antisite_defect(st, cation_positions = None):
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
    write_xyz(st, file_name = st.name+'_antisite_start')
    st = st.mov_atoms(i_alk, x_tr)
    st = st.mov_atoms(i_tr,  x_alk)

    write_xyz(st, file_name = st.name+'_antisite_final')

    printlog('Atom ',i_alk+1,'and', i_tr+1,'were swapped')
    printlog('The distance between them is ', sur[3][0])

    st.magmom = [None]

    return st
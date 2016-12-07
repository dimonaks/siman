
import sys,copy,itertools
import numpy as np

from header import printlog

from functions import xcart2xred
# sys.path.append('/home/aksenov/Simulation_wrapper/') 
# sys.path.append('/home/aksenov/Simulation_wrapper/savelyev') 


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


    printlog('I have found ', len(positions), 'non-equivalent positions for', element, ':',positions.keys(), imp = 'y')
    # print(positions)
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



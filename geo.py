
import sys,copy,itertools
import numpy as np
from functions import xcart2xred
sys.path.append('/home/aksenov/Simulation_wrapper/') 
sys.path.append('/home/aksenov/Simulation_wrapper/savelyev') 


def ortho_vec(rprim, ortho_sizes = None):
    """
    Function returns mul_mat - 3 vectors of integer numbers (ndarray)
    By calculating rprim*mul_mat.T you will get rprim of orthogonal supercell (actually as close as possible to it) 
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



def create_supercell(st, mul_matrix, test_overlap = False): 
    """ 
    st (Structure) -  
    mul_matrix (3x3 ndarray of int) - for example created by *ortho_    vec()* 
    test_overlap - quite slow
    """ 
    sc = copy.deepcopy(st) 

    sc.rprimd = list(np.dot(mul_matrix, st.rprimd  ))
    print('rprimd of supercell\n',np.round(sc.rprimd,1))
    sc.vol = np.dot( sc.rprimd[0], np.cross(sc.rprimd[1], sc.rprimd[2])  )
    st.vol = np.dot( st.rprimd[0], np.cross(st.rprimd[1], st.rprimd[2])  )
    sc_natom_i = int(sc.vol/st.vol*st.natom) # test
    sc_natom = sc.vol/st.vol*st.natom # test
    print('The supercell should contain ', sc_natom,sc_natom_i, 'atoms')
    sc.xcart = []
    sc.typat = []
    sc.xred  = []
    #find range of multiplication
    mi = np.min(mul_matrix, axis = 0)
    ma = np.max(mul_matrix, axis = 0)
    mi[mi>0] = 0  # 


    bound = 1e-7 # to overcome numerical errors
    mp    = 2    # include additionall space for searching atoms
    for uvw in itertools.product(*[range(*z) for z in zip(mi-mp, ma+mp)]): #loop over all ness uvw

        xcart_mul = st.xcart + np.dot(uvw, st.rprimd) # coordinates of basis for each uvw
        xred_mul  = xcart2xred(xcart_mul, sc.rprimd)
        
        for xr, xc,  t in zip(xred_mul, xcart_mul, st.typat):
            
            if all([0-bound <= r < 1-bound for r in xr]): #only that in sc.rprimd box are needed
                sc.xcart.append( xc )
                sc.xred.append ( xr )
                sc.typat.append( t  )
    


    sc.natom = len(sc.xcart)
    sc.magmom = [None]
    # nznucl&
    # print np.array(sc.xred).round(5)

    if abs(sc.natom - sc_natom)>1e-5: #test 1, number of atoms
        print('Error! Supercell contains wrong number of atoms:', sc_natom, 'instead of', sc.natom , 
            'try to increase *mp* or even play with *bound*')
        raise RuntimeError
    else:
        print('Number of atoms ... OK')
    
    if test_overlap: #test 2: overlapping of atoms
        enx = list(enumerate(sc.xcart))
        for (i1, x1), (i2, x2) in itertools.product(enx, enx):
            # print image_distance(x1, x2, sc.rprimd)[0]
            if i1 != i2 and image_distance(x1, x2, sc.rprimd)[0] < 0.1: #less than 0.1 Angstrom
                print('Error! Atoms in supercell are overlapping. Play with *bound*')
                raise RuntimeError



    return sc

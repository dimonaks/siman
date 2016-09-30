
import sys
import numpy as np
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
    mul_matrix = np.array(a.vec_new_in_old).astype(int)

    return mul_matrix
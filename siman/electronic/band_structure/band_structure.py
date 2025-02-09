import os
import sys
from siman.header import db
from siman.calc_manage import smart_structure_read, add_loop, res_loop
from siman.set_functions import read_vasp_sets
import project_sets
import time
import itertools
from mp_api.client import MPRester
import seekpath
import re
from project_sets import sp_pack, band_pack
from siman.header import printlog
from siman import header

atoms_dict ={'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64, 'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100, 'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108, 'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Uuq': 114, 'Uuh': 116}

def create_POSCAR(st, st_name, home_path, debug=False):
    """
    This function creates the POSCAR file for particular structure

    INPUT:
        - st <class 'pymatgen.core.structure.Structure'> - structure
        - st_name (str) - name of structure
        - home_path (str) - path to folder with calculations of band structure; here the folder with st_name will be created and collect all calculated data
        - debug (bool) - if True, displays logs
    RETURN:
        - path_to_POSCAR (str) - path to created POSCAR
        - path_to_main_folder (str) - path to folder with all future calculations
    """
    folder_name = st_name
    path_to_main_folder = f'{home_path}/{folder_name}'

    os.makedirs(path_to_main_folder, exist_ok=True)
    st.to(f"{st_name}.POSCAR", f'{path_to_main_folder}')

    path_to_POSCAR = f'{path_to_main_folder}/{st_name}.POSCAR'

    if debug:
        print(f'\n\nFrom create_POSCAR: the file {path_to_POSCAR} was crested')

    return path_to_POSCAR, path_to_main_folder


def find_kpath(st=None, method="AFLOW", path_to_POSCAR=None, debug=False):
    """
    The function finds a highly symmetric k-path for a given structure and creates an input KPOINTS file for calculations
    Currently, two path detection methods are supported: using libraries AFLOW and seekpath
    Also you can use your own KPOINTS file by using method='reference'

    INPUT:
        - st <class 'pymatgen.core.structure.Structure'> - structure
        - method (str) - the method of path detection ("AFLOW", "seekpath" or "reference")
        - path_to_POSCAR (str) - path to POSCAR file with the structure
        - debug (bool) - if True, shows the obtained k-path to created KPOINTS
    RETURN:
        - k_path (list of tuples) - list with coordinates and names of k-points in the obtained k-path in the format for siman set
                                    Example: [128 */number of points per line/*), ('G', 0.0, 0.0, 0.0), ('X', 0.5, 0.0, 0.5), ('W', 0.5, 0.25, 0.75), ...]
    """
    current_folder = os.getcwd()

    if method == "AFLOW":
        if not path_to_POSCAR:
            print("\n\n! ERROR --- From find_kpath with AFLOW: there is no path_to_POSCAR\n\n")
            exit()
        if not header.PATH2AFLOW:
            print("\n\n! ERROR --- From find_kpath with AFLOW: there is no path to AFLOW software. Please, add PATH2AFLOW in the simanrc.py\n\n")
            exit()

        os.environ["PATH"] = f"{os.path.dirname(header.PATH2AFLOW)}:{os.environ.get('PATH', '')}"
        os.system(f"{header.PATH2AFLOW} --kpath < {path_to_POSCAR} > KPOINTS")

        k_path = band_set_from_KPOINTS()
        os.remove('KPOINTS')

    elif method == "seekpath":
        if not st:
            print("\n\n! ERROR --- From find_kpath with seekpath: there is no structure info\n\n")
            exit()
        k_path = seekpath_method(st)

    elif method == 'reference':
        k_path = band_set_from_KPOINTS()

    else:
        print("\n\nERROR: The method of k-path detection is set incorrectly!\n\n")
        exit()

    if debug:
        print("\n\nFrom find_kpath: generated k_path is\n\n", k_path)

    os.chdir(f'{current_folder}')

    return k_path


def seekpath_method(st):
    """
    The function finds a highly symmetric k-path by seekpath library

    INPUT:
        - st <class 'pymatgen.core.structure.Structure'> - structure
    RETURN:
        - k_path (list of tuples) - list with coordinates and names of k-points in the obtained k-path in the format for siman set
                                    Example: [128 */number of points per line/*, ('G', 0.0, 0.0, 0.0), ('X', 0.5, 0.0, 0.5), ('W', 0.5, 0.25, 0.75), ...]
    TO DO:
        Eliminate the atoms_dict
    """
    lattice = list(list(elem) for elem in st.lattice.matrix)
    coords = []
    atoms_num = []
    for site in st.sites:
        coords.append(list(site.frac_coords))
        atoms_num.append(atoms_dict[site.specie.name])
    seekpath_out = seekpath.get_path((lattice, coords, atoms_num))
    points_coord, path_list = seekpath_out['point_coords'], seekpath_out['path']

    k_path = [40]
    prev_point = None
    for points in path_list:
        for point in points:
            coord = points_coord[point]
            point = point[0]
            if prev_point != point:
                prev_point = point
                elem = (point, *coord)
                k_path.append(elem)
    return k_path


def band_set_from_KPOINTS(path_to_KPOINTS=None):
    """
    The function gets k-path in format for siman set from the ready KPOINTS file

    INPUT:
        - path_to_KPOINTS (str) - path to KPOINTS file; if None, the file will be searched in the current directory
    RETURN:
        - k_path (list of tuples) - list with coordinates and names of k-points in the obtained k-path in the format for siman set
                                    Example: [128 */number of points per line/*, ('G', 0.0, 0.0, 0.0), ('X', 0.5, 0.0, 0.5), ('W', 0.5, 0.25, 0.75), ...]
    """
    if path_to_KPOINTS:
        os.chdir(path_to_KPOINTS)

    k_path = []
    pattern = r'(?i)Gamma'
    with open('KPOINTS', 'r') as file:
        for line in file:
            if "!" in line:
                data = line.replace("\\", "").split()
                data = [re.sub(pattern, 'G', elem) for elem in data]
                if 'grids' in data:
                    k_path.append(int(data[0])) # number of points per line
                else:
                    point_name = data[-1]
                    if k_path and (len(k_path) == 1 or (len(k_path) >= 1 and k_path[-1][0] != point_name)):
                        k_path.append((point_name, float(data[0]), float(data[1]), float(data[2])))
    return k_path
















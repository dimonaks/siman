"""
Author: Kobernik Tatiana
"""
import os
import re
from siman.calc_manage import smart_structure_read, add_loop, res_loop
from siman.set_functions import read_vasp_sets
from siman.header import printlog
from siman import header


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

    st.to(filename=f"{path_to_main_folder}/{st_name}.POSCAR", fmt='poscar')

    path_to_POSCAR = f'{path_to_main_folder}/{st_name}.POSCAR'

    if debug:
        printlog(f'\n\nFrom create_POSCAR: the file {path_to_POSCAR} was crested', imp = 'y')

    return path_to_POSCAR, path_to_main_folder


def find_kpath(st=None, siman_st=None, method="AFLOW", path_to_POSCAR=None, refer_path_to_KPOINTS=None, debug=False):
    """
    The function finds a highly symmetric k-path for a given structure and creates an input KPOINTS file for calculations
    Currently, two path detection methods are supported: using libraries AFLOW and seekpath
    Also you can use your own KPOINTS file by using method='reference'

    INPUT:
        - st <class 'pymatgen.core.structure.Structure'> - structure
        - siman_st <class 'siman.core.structure.Structure'> - structure
        - method (str) - the method of path detection ("AFLOW", "seekpath" or "reference")
        - path_to_POSCAR (str) - path to POSCAR file with the structure
        - refer_path_to_KPOINTS (str) - path to reference KPOINTS file for band structure
        - debug (bool) - if True, shows the obtained k-path to created KPOINTS
    RETURN:
        - k_path (list of tuples) - list with coordinates and names of k-points in the obtained k-path in the format for siman set
                                    Example: [128 */number of points per line/*), ('G', 0.0, 0.0, 0.0), ('X', 0.5, 0.0, 0.5), ('W', 0.5, 0.25, 0.75), ...]
    """
    
    if method == "AFLOW":
        if not path_to_POSCAR:
            printlog("ERROR --- From find_kpath with AFLOW: there is no path_to_POSCAR", imp = 'y')
            exit()
        if not header.PATH2AFLOW:
            printlog("ERROR --- From find_kpath with AFLOW: there is no path to AFLOW software. Please, add PATH2AFLOW in the simanrc.py", imp = 'y')
            exit()

        os.environ["PATH"] = f"{os.path.dirname(header.PATH2AFLOW)}:{os.environ.get('PATH', '')}"
        path_to_dir = '/'.join(path_to_POSCAR.split('/')[:-1])
        path_to_KPOINTS = os.path.join(path_to_dir, 'KPOINTS.band')

        os.system(f"{header.PATH2AFLOW} --kpath < {path_to_POSCAR} > {path_to_KPOINTS}")

        k_path = band_set_from_KPOINTS(path_to_KPOINTS=path_to_KPOINTS)
        os.remove(f"{path_to_KPOINTS}")

    elif method == "seekpath":
        if not st and not siman_st:
            printlog("ERROR --- From find_kpath with seekpath: there is no structure info", imp = 'y')
            exit()
        k_path = seekpath_method(st=st, siman_st=siman_st)

    elif method == 'reference':
        if not refer_path_to_KPOINTS:
            printlog("ERROR --- From find_kpath with reference: there is no refer_path_to_KPOINTS", imp = 'y')
            exit()
        k_path = band_set_from_KPOINTS(path_to_KPOINTS=refer_path_to_KPOINTS)

    else:
        printlog("ERROR: The method of k-path detection is set incorrectly!", imp = 'y')
        exit()

    if debug:
        printlog(f"From find_kpath: generated k_path is {k_path}", imp = 'y')

    return k_path


def seekpath_method(st=None, siman_st=None):
    """
    The function finds a highly symmetric k-path by seekpath library

    INPUT:
        ** you need to provide one of them **
        - st <class 'pymatgen.core.structure.Structure'> - structure 
        - siman_st <class 'siman.core.structure.Structure'> - structure
    RETURN:
        - k_path (list of tuples) - list with coordinates and names of k-points in the obtained k-path in the format for siman set
                                    Example: [128 */number of points per line/*, ('G', 0.0, 0.0, 0.0), ('X', 0.5, 0.0, 0.5), ('W', 0.5, 0.25, 0.75), ...]
    """
    import seekpath
    from siman.functions import element_name_inv

    seekpath_input = []
    if st:
        lattice = list(list(elem) for elem in st.lattice.matrix)
        coords = []
        atoms_num = []
        for site in st.sites:
            coords.append(list(site.frac_coords))
            atoms_num.append(element_name_inv(site.specie.name))
        seekpath_input = (lattice, coords, atoms_num)
    elif siman_st:
        seekpath_input.append(siman_st.rprimd)
        seekpath_input.append(siman_st.xred)
        seekpath_input.append([element_name_inv(elem) for elem in siman_st.get_elements()])
        seekpath_input = (seekpath_input)
    else:
        printlog('ERROR: You did not provide the structure for seekpath method', imp = 'y')

    seekpath_out = seekpath.get_path(seekpath_input)
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


def band_set_from_KPOINTS(path_to_KPOINTS):
    """
    The function gets k-path in format for siman set from the ready KPOINTS file

    INPUT:
        - path_to_KPOINTS (str) - path to KPOINTS file
    RETURN:
        - k_path (list of tuples) - list with coordinates and names of k-points in the obtained k-path in the format for siman set
                                    Example: [128 */number of points per line/*, ('G', 0.0, 0.0, 0.0), ('X', 0.5, 0.0, 0.5), ('W', 0.5, 0.25, 0.75), ...]
    """

    k_path = []
    pattern = r'(?i)Gamma'
    with open(path_to_KPOINTS, 'r') as file:
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
















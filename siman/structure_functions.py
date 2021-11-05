#!/usr/bin/env python3
""" 
Author: Kartamyshev A.I. (Darth Feiwante)
"""

def inherit_icalc_isotropic(new_structure = '', start_new_version = None,  base_calculation = (None, None, None), database = None, min_mult = 1, max_mult = 1, num_points = 2, geo_folder = '', it_folder = '', override = False):
    """
    This function makes set of structures uniformly scaled from the initial one within the range of deformation

    INPUT:
        - new_structure (str)       - arbitary name for your crystal structure 
        - start_new_version (int)   - start version for newly built structures
        - base_calculation (tuple)  - tuple describing initial Calculation object in form ('structure', 'set', 'version')
        - database (dict)           - dictionary with the project's results
        - min_mult (int)            - minimal deformation of the initial structure
        - max_mult (int)            - maximal deformation of the initial structure
        - num_points (int)          - number of newly built structures
        - geo_folder (str)          - path to the folder to save *.geo files of newly built structures
        - it folder (str)           - section folder
        - override (boolean)        - if True then the old structures with the same names will be overwritten
        
    RETURN:
        None
    SOURCE:
        None
    TODO:
        Some improvements
    """

    from calc_manage import inherit_icalc
    min_mult = min_mult
    max_mult = max_mult
    num_points = num_points
    step = (max_mult - min_mult)/(num_points - 1)
    mult_list = [min_mult+step*i for i in range(num_points)]
    version = start_new_version
    for j in mult_list:
        inherit_icalc('isotropic',   new_structure,  version, base_calculation, database, mult_rprimd = j, geo_folder=geo_folder, override=override)
        version += 1

def inherit_icalc_c_a(new_structure = '', start_new_version = None,  base_calculation = (None, None, None), database = None, min_mult_a = 1, max_mult_a = 1, num_points_a = 2, 
                      min_mult_c = 1, max_mult_c = 1,num_points_c = 2, geo_folder='', it_folder =''):
    """
    This function makes set of structures deformed uniformly in the plane presented by the vectors 1 and 2 and separately deformed along the vector 3 of the lattice

    INPUT:
        - new_structure (str)       - arbitary name for your crystal structure 
        - start_new_version (int)   - start version for newly built structures
        - base_calculation (tuple)  - tuple describing initial Calculation object in form ('structure', 'set', 'version')
        - database (dict)           - dictionary with the project's results
        - min_mult_a (float)        - minimal simultaneous deformation of the vector 1 and 2 of the final structure from "base_calculation"
        - max_mult_a (float)        - maximal simultaneous deformation of the vector 1 and 2 of the final structure from "base_calculation"
        - num_points_a (int)        - number of different simultaneous deformations of the vectors 1 and 2
        - min_mult_c (float)        - minimal deformation of the vector 3 of the structure from "base_calculation"
        - max_mult_c (float)        - maximal deformation of the vector 3 of the structure from "base_calculation"
        - num_points_c (int)        - number of different deformations of the vector 3 from "base_calculation"
        - geo_folder (str)          - path to the folder to save *.geo files of newly built structures
        - it folder (str)           - section folder
        - override (boolean)        - if True then the old structures with the same names will be overwritten

    RETURN:
        None
    SOURCE:
        None
    TODO:
        Some improvements
    """
    from classes import inherit_icalc
    
    if num_points_a > 1:    
        # Lattice parameter a
        min_mult_a = min_mult_a
        max_mult_a = max_mult_a
        num_points_a = num_points_a
        step_a = (max_mult_a - min_mult_a)/(num_points_a - 1)
        mult_list_a = [min_mult_a+step_a*i for i in range(num_points_a)]

    if num_points_c > 1:  
        # Lattice parameter c
        min_mult_c = min_mult_c
        max_mult_c = max_mult_c
        num_points_c = num_points_c
        step_c = (max_mult_c - min_mult_c)/(num_points_c - 1)
        mult_list_c = [min_mult_c+step_c*i for i in range(num_points_c)]

    print('database', database)

    version = start_new_version

    if num_points_a > 1 and num_points_c > 1:    
        for j in mult_list_a:
            for k in mult_list_c:
                inherit_icalc('c_a',   new_structure,  version, base_calculation, database, mult_a = j, mult_c = k, geo_folder=geo_folder)
                version += 1

    elif num_points_c == 1:
        for j in mult_list_a:
            inherit_icalc('c_a',   new_structure,  version, base_calculation, database, mult_a = j, mult_c = 1, geo_folder=geo_folder, override=override)
            version += 1

    elif num_points_a == 1:
        for j in mult_list_c:
            inherit_icalc('c_a',   new_structure,  version, base_calculation, database, mult_a = 1, mult_c = j, geo_folder=geo_folder, override=override)
            version += 1
            
def inherit_icalc_x_y(new_structure = '', start_new_version = None,  base_calculation = (None, None, None), database = None, 
                      min_mult_a = 1, max_mult_a = 1, num_points_a = 2, min_mult_b = 1, max_mult_b = 1,num_points_b = 2, geo_folder='', it_folder ='',
                      override = False):
    """
    This function makes set of structures separately deformed along the vectors 1 and 2 of the lattice

    INPUT:
        - new_structure (str)       - arbitary name for your crystal structure 
        - start_new_version (int)   - start version for newly built structures
        - base_calculation (tuple)  - tuple describing initial Calculation object in form ('structure', 'set', version)
        - database (dict)           - dictionary with the project's results
        - min_mult_a (float)        - minimal deformation of the vector 1 of the structure from "base_calculation"
        - max_mult_a (float)        - maximal deformation of the vector 1 of the structure from "base_calculation"
        - num_points_a (int)        - number of different deformations of the vector 2
        - min_mult_b (float)        - minimal deformation of the vector 2 of the structure from "base_calculation"
        - max_mult_b (float)        - maximal deformation of the vector 2 of the structure from "base_calculation"
        - num_points_b (int)        - number of different deformations of the vector 2
        - geo_folder (str)          - path to the folder to save *.geo files of newly built structures
        - it folder (str)           - section folder
        - override (boolean)        - if True then the old structures with the same names will be overwritten

    RETURN:
        None
    SOURCE:
        None
    TODO:
        Some improvements
    """
    from calc_manage import inherit_icalc


    if num_points_a > 1:
        # Coordinate x in rprimd
        step_a = (max_mult_a - min_mult_a)/(num_points_a - 1)
        mult_list_a = [min_mult_a+step_a*i for i in range(num_points_a)]

    if num_points_b > 1:
        # Coordinate y in rprimd
        step_b = (max_mult_b - min_mult_b)/(num_points_b - 1)
        mult_list_b = [min_mult_b+step_b*i for i in range(num_points_b)]


    version = start_new_version

    if num_points_a > 1 and num_points_b > 1:
        for j in mult_list_a:
            for k in mult_list_b:
                inherit_icalc('xy',   new_structure,  version, base_calculation, database, mult_a = j, mult_b = k, geo_folder=geo_folder, it_folder = it_folder, override=override)
                version += 1
    elif num_points_b == 1:
        for j in mult_list_a:
            inherit_icalc('xy',   new_structure,  version, base_calculation, database, mult_a = j, mult_b = 1, geo_folder=geo_folder, it_folder = it_folder, override=override)
            version += 1

    elif num_points_a == 1:
        for j in mult_list_b:
            inherit_icalc('xy',   new_structure,  version, base_calculation, database, mult_a = 1, mult_b = j, geo_folder=geo_folder, it_folder = it_folder, override=override)
            version += 1



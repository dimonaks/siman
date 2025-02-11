"""
Author: Kobernik Tatiana
"""

from siman.set_functions import read_vasp_sets
import project_sets
from siman.header import printlog


def create_new_set(new_set, old_set, parameters, k_path=None, path_key=None, debug=False):
    """
    The function creates a new set in siman based on an existing one.
    It can add or change calculation parameters,
    as well as add k-paths to calculate the band structure or effective mass

    INPUT:
        - new_set (str) - name of the old set defined in the varset
        - old_set (str) - name of the new set, which will be created
        - parameters (dict) - dictionary with new parameters of calculations, like {'LORBIT': 11, 'ICHARG': 11}
        -  k_path (list of tuples) - list with coordinates and names of k-points in the k-path
                                    for band structure or effective mass calculations in the format for siman set
                                    Example: [128 */number of points per line/*, ('G', 0.0, 0.0, 0.0), ('X', 0.5, 0.0, 0.5), ('W', 0.5, 0.25, 0.75), ...]
        - path_key (str) - set's key that determines the type of calculations when adding a k-path: 'k_band_structure' or 'effective_mass'
        - debug - if True, displays the parameters of the new set and all available sets in varset
    RETURN:
        None
    """

    varset = read_vasp_sets(project_sets.user_vasp_sets, override_global=0)

    if path_key and k_path:
        parameters[path_key] = k_path
        if debug:
            print('\n\nFrom create_new_set: the k-path was added in parameters of new set')

    read_vasp_sets([(new_set, old_set, parameters, 'override')])
    printlog(f"New set [{new_set}] has been created\n\n")

    if debug:
        varset[new_set].printme()
        print("Available sets: ", varset.keys(), '\n')
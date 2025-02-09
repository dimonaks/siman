import os
from siman.header import db
from siman.calc_manage import add_loop, res_loop
from siman.set_functions import read_vasp_sets
import project_sets
import time
import itertools
from siman.header import printlog


def charge_calc(st_name, st, set):
    """
    This function manages single-point calculations via VASP; it includes a part with waiting for the end of calculations.
    It uses for obtaining CHGCAR for band structure calculations

    INPUT:
        - st <class 'pymatgen.core.structure.Structure'> - structure
        - st_name (str) - name of structure
        - set (str) - name of siman's set for calculation; must be defined in the varset
    RETURN:
        None
    """

    add_loop(st_name, setlist=set, verlist=1, input_st=st, it_folder=st_name, run=2) #charge_density

    spinner = itertools.cycle(['-', '\\', '|', '/'])  # Animation of the spinner

    while res_loop(st_name, setlist=set, verlist=1, up='x', log_flag=False) == ({}, []):
        print(f"\rWaiting for the calculation to be completed {next(spinner)}", end="", flush=True)
        time.sleep(0.4)


def non_self_consist_calc(st_name, new_set, old_set, path_to_POSCAR, path_to_vasprun=None):
    """
    This function manages non-self-consistent calculations via VASP; it includes a part with waiting for the end of calculations.
    It uses for band structure and effective mass calculations

    INPUT:
        - st_name (str) - name of structure
        - new_set (str) - name of the siman's set which was used for the previous single-point calculations; must be defined in the varset
        - old_set (str) - name of the siman's set for current non-self-consistent calculation; must be defined in the varset
        - path_to_POSCAR (str) - path to the POSCAR file of input structure
        - path_to_vasprun (str) - it is necessary if the calculation is restarted again;
                                  this variable specifies the path to an outdated previously received vasprun.xml file that is being deleted.
                                  If you do not specify this parameter, but the new vasprun.xml will not be pulled into the folder.
    RETURN:
        None
    """

    if path_to_vasprun:
        if os.path.exists(path_to_vasprun):
            try:
                os.remove(path_to_vasprun)
                print(f"The old file {path_to_vasprun} is deleted.")
            except:
                print(f'The defined file {path_to_vasprun} was not found')


    add_loop(st_name, setlist=old_set, verlist=1, ise_new=new_set, inherit_option='full', savefile='ocxe',
             input_geo_file = path_to_POSCAR, it_folder = st_name, override=1, run=2)

    spinner = itertools.cycle(['-', '\\', '|', '/'])  # Animation of the spinner

    while res_loop(st_name + '.if', setlist=[new_set], verlist=[1], show = 'fo', log_flag=False) == ({}, []):
        print(f"\rWaiting for the calculation to be completed {next(spinner)}", end="", flush=True)
        time.sleep(0.4)

    db[st_name + '.if', new_set, 1].get_file('1.vasprun.xml')


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
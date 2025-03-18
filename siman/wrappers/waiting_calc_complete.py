"""
Author: Kobernik Tatiana
"""

import os
from siman.header import db
from siman.calc_manage import add_loop, res_loop
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
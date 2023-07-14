#!/usr/bin/env python3
""" 
Author: Kartamyshev A.I. (Darth Feiwante)
"""



def formation_energy(database = None, calc_def = (), calc_id = ()):
    """
    This function is used to find minimal distance in the lattice presented 
    as a part of the 'Calculation' object

    INPUT:
        - database (.gbdm3) - database dictionary; could be provided; if not then are taken from header
        - calc_def (tuple) - tuple describing the Calculation object for the 
                             lattice containing a defect in form ('structure', 'set', version)
        - calc_id (tuple) - tuple describing the Calculation object for the 
                             lattice without a defect in form ('structure', 'set', version)
    RETURN:
        None
    SOURCE:
        None
    TODO:
        - Add different type of defects
    """
    defect = database[calc_def]
    ideal = database[calc_id]
    n_at_def = defect.natom
    e_def = defect.energy_free
    e_id_at = ideal.energy_free/ideal.natom
    E_f = e_def - n_at_def*e_id_at
    print('Formation energy for defect '+calc_def[0]+' = '+str(E_f)+' eV')


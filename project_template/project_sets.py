# -*- coding: utf-8 -*-
"""
List of user sets (e.g. for VASP) obtained on inheritance principle, see user_vasp_sets
Syntax:  ("set_new", "set_old", {"param1":value1, "param2":value2, ...}, over)
        - set_new (str) - name of new set
        - set_old (str) - name of base set used for creating new set
        - {} (dict)     - dictionary of parameters to be updated in set_new
        - over (str)   - by default is '', if 'over' then set_new is reinitializied (override), otherwise it is protected
Naming conventions of sets:
    0 - no relaxation, single point calculation
    1 - relax only atoms
    2 - full relax
    4 - relax only cell shape at fixed volume
"""
from vasp_params import *

user_vasp_sets = [ # this a description of user sets for VASP, based on inheritance principle
('0', 'static', pot_pack), # static is a predefined set; choose required potentials
('1',    '0', ion_relax_packet),
('0m',   '0', mag_packet ),
('1m',   '1', mag_packet ),
('0u',    '0m', dftu_packet ),
('1u',    '1m', dftu_packet ),
('4' ,'1', {'ISIF':4} ), #relax everything except volume
('4m' ,'1m', {'ISIF':4} ), #relax everything except volume
('4u' ,'1u', {'ISIF':4} ), #relax everything except volume
('dos' ,'0', dos_pack,  ), # set for dos calculations
]


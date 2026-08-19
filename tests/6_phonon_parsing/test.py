"""
Test parsing phonon calculations and plot phonons 

phonopy should be installed pip install phonopy

Author: Dmitry Aksyonov
"""

from siman.header import db, _update_configuration
from siman.calc_manage import add, res
from siman.set_functions import init_default_sets, read_vasp_sets
from siman.calc_manage import smart_structure_read
from siman.database import read_database, write_database


if 0:
    init_default_sets(1)
    _update_configuration('/home/d.aksenov/simanrc.py') # read configuration, required to run job

    read_vasp_sets([('phdos', 'static', {'ISIF':2, 'IBRION':6, 'POTIM':0.015, 'NFREE':2, 'NSW':1, # create new set 'phdos' from 'static' for normal mode calculations
    'PREC':'Accurate', 'ADDGRID':'.TRUE.', 'EDIFF':1e-08, 'NPAR':None, 'savefile':'x'})]) #  finite difference with two displacements by 0.015 A
    
    st = smart_structure_read('LiF.cif')
    add('LiF', 'phdos', 1, input_st = st, it_folder = 'LCO/phonon', run = 2)
    write_database()
else:
    read_database()
    cl = db['LiF', 'phdos', 1]
    cl.res(up = 'up2', show='freq')
    # cl.read_pdos_using_phonopy()
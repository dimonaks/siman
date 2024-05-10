#!/usr/bin/python3
# -*- coding: utf-8 -*- 
"""
To use this template, please install Siman package https://github.com/dimonaks/siman/wiki
"""
from __future__ import division, unicode_literals, absolute_import, print_function

if 1:
    import sys
    # sys.path.append('/home/aksenov/Simulation_wrapper/siman') #path to siman package
    from siman import header
    from siman.header import printlog, runBash

    from siman.SSHTools import SSHTools
    from siman.calc_manage   import (smart_structure_read, update_des, add_loop, res_loop, add, res, complete_run)
    from siman.database      import read_database, write_database
    from siman.set_functions import read_vasp_sets

    
    if 0:
        #run this once to make migration from old database
        from siman.header import pickle_module_migration_script
        pickle_module_migration_script()


    header.conv, header.varset, size_on_start = read_database()
    header.struct_des = update_des(header.struct_des, header.MANUALLY_ADDED); #read manually added calculations from project_conf.py file
    db                = header.db # main database dictionary

    import project_sets # should be after read_database
    varset = read_vasp_sets(project_sets.user_vasp_sets, override_global = 0) #read user sets




"""Control"""
save = 1
header.warnings = 'neyY'
header.warnings = 'yY'




"""Start working"""









"""End working"""









complete_run(header.close_run)



if save:
    write_database(db, header.conv, header.varset, size_on_start)








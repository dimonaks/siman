#!/usr/bin/python3
# -*- coding: utf-8 -*- 
"""
To use this template, please install Siman package https://github.com/dimonaks/siman/wiki
"""
if 1:
    import sys
    from siman import header
    from siman.header import printlog, runBash

    from siman.SSHTools import SSHTools
    from siman.calc_manage   import (smart_structure_read, update_des, add_loop, res_loop, add, res, complete_run)
    from siman.database      import read_database, write_database
    from siman.set_functions import read_vasp_sets


    header.conv, header.varset, size_on_start = read_database()
    header.struct_des = update_des(header.struct_des, header.MANUALLY_ADDED); #read manually added calculations from project_conf.py file
    db                = header.db # main database dictionary

    import project_sets # should be after read_database
    # header._update_configuration('/home/username/simanrc.py') #here you can provide global control of siman for all projects
    header._update_configuration('./project_conf.py') # here you can put control specific for this project

    varset = read_vasp_sets(project_sets.user_vasp_sets, override_global = 0) #read user sets, set override_global=1 if you want to reread all sets


"""Control: Below is some of the frequently used parameters. For more check project_conf.py"""
save = 1
# header.warnings = 'neyY' # show more warnings
header.warnings = 'yY' # show less  warnings
# header.check_job = 1 # check job status in the queue on cluster
# header.siman_run = 0 # special simplified regime for siman, see documentation
# header.copy_to_cluster_flag = 0 # allows to prevent copying files to cluster
# header.corenum = 4 # overwrite  number of cores used for calculations.






"""Start working"""



#put here your commands, such as add() and res()




"""End working"""
complete_run(header.close_run)



if save:
    write_database(db, header.conv, header.varset, size_on_start)








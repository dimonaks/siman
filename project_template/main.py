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
# header.check_job = 1
# header.siman_run = 0
# header.copy_to_cluster_flag = 0
# header.corenum = 5
# header.schedule_system = 'SLURM'



"""Start working"""
# header.ssh_object = SSHTools()
# header.ssh_object.setup(user="aksenov",host="fri.cnm.utexas.edu",pkey="/home/aksenov/.ssh/id_rsa")

# print(runBash('ssh aksenov@fri.cnm.utexas.edu "pwd"'))








"""End working"""









complete_run(header.close_run)



if save:
    write_database(db, header.conv, header.varset, size_on_start)




























"""
TODO:
исправить calc[id].path["output"] для U_ramping - вроде сделано


0) read_geometry(), Если переменная не найдена, то подставляется [None] - это не очень удобно и лучше сделать просто None.
1) read project_conf explicitly from here and not from header - maybe not needed. The values from project_conf can be needed everywhere in siman and header is universal file for siman 
2) перенести настройки matplotlib из header в конкретные функции, которые строят графики
3) project_sets.update_sets(varset) нужно удалить, она остается пока для этого проекта

4) inherit_option = continue и sequence_set совместно не тестировались!

5) sequence_set and self.associated_outcars ????

!How to make tables and pictures more straightforward?
!How to make inheritance of last relaxed configuration more straightforward? - добавить возможность продолжения расчёта

!Для нового проекта нужно подумать об объединении папки geo с исходными структурами и выходной папки; Лучше все что касается отдельного расчета хранить в одной папке, просто использовать разные имена для файлов или подпапки.


!Добавление нового атома подразумевает набор стандартных действий. Написать маленькую функцию для этого. Сейчас код подобной функции используется
в двух местах: внутри create_segregation_cases() и add_impurity().add()


!gbpos в самом старте определяется вручную для первой версии и просто копируется для других версий.


!make_incar_and_copy_all проверить magmom



! В классе Structure() добавить методы: 
удалить атом, добавить атом, заменить атом; 
наследовать rprimd; растянуть;
потом с помощью этих методов упростить функции inherit_icalc и add_impurity


Changes to siman2; please move this section to siman2 folder.
1. latex_table() moved to functions.py


"""


# -*- coding: utf-8 -*-
"""
Control of project
"""
from __future__ import division, unicode_literals, absolute_import 

"""General parameters"""
PATH2POTENTIALS = '~/PAW_PBE_VASP' # put potential into this directory
pmgkey = "" # please get your own key from materials project for pymatgen: 16 digits (https://legacy.materialsproject.org/dashboard) 
mpkey = "" # please get your own key from materials project api: 32 digits (https://next-gen.materialsproject.org/api)
PATH2DATABASE        = './database'
AUTO_UPDATE_DB = False # if True, then execute write_database() at the end of add() and res()
PATH2JMOL     = 'jmol'
PATH2PHONOPY  = 'phonopy'
PATH2NEBMAKE = '~/vts/nebmake.pl' # http://theory.cm.utexas.edu/vasp/scripts.html
geo_folder           = './' # duplicate structures using .geo format files in this folder

"""Cluster parameters"""
PATH2PROJECT = 'my_project' # path to project on cluster relative to home folder
DEFAULT_CLUSTER = 'mag'

#description of user clusters:
username = 'username' # user name at cluster
CLUSTERS = {}
CLUSTERS['mag'] = {
'address':username+'@10.30.16.168',
'vasp_com':'mpirun vasp_std',
'gaussian_command':'g16',
'homepath':'/home/'+username+'/',
'schedule':'SLURM',
'corenum':4,
'memory':24, #GB, 
'partition':'AMG-medium',
# 'any_command':['--nodelist=node-amg01,node-amg02,node-amg03,node-amg09,node-amg12'],
'modules':'module load Compiler/Intel/16u4 Q-Ch/VASP/5.4.4 ScriptLang/python/3.6i_2018u3; \nulimit -s unlimited\n' 
}


CLUSTERS['ar'] = {
'address':'username@ip',
'vasp_com':'mpirun  vasp_std',
'homepath':'/home/username/',
'schedule':'PBS',
'walltime':'72:00:00',
'corenum':16,
'pythonpath':'/usr/lib64/python2.7/site-packages/numpy',
'modules':'module load Compilers/Intel/psxe_2015.6; module load MPI/intel/5.1.3.258/intel; module load QCh/VASP/5.4.1p1/psxe2015.6; module load ScriptLang/python/2.7',
}




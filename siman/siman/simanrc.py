# -*- coding: utf-8 -*-
"""
User-related parameters for siman, file is installed to home folder

"""
from __future__ import division, unicode_literals, absolute_import 
from siman.header import CLUSTERS

"""Cluster constants"""
DEFAULT_CLUSTER = 'pardus'
PATH2ARCHIVE = '' # path to archive; if no files are found at home folder, siman will check here; relative paths should be same

CLUSTERS['CEE'] = {
'address':'username@10.30.16.62', # command for ssh
'schedule':'SLURM',
'corenum':4,
'pythonpath':'/usr/lib64/python2.7/site-packages/numpy',
'vasp_com':'prun /opt/vasp/bin/vasp5.4.1MPI_aksenov',
'modules':'module add prun/1.0; module add intel/16.0.2.181; module add impi/5.1.3.181\n'
}

CLUSTERS['pardus'] = {
'address':'name.surname@10.30.17.12',
'vasp_com':'mpirun  vasp_std',
'schedule':'PBS',
'walltime':'72:00:00',
'corenum':4,
'pythonpath':'/usr/lib64/python2.7/site-packages/numpy',
'modules':'module load Compilers/Intel/psxe_2015.6; module load MPI/intel/5.1.3.258/intel; module load QCh/VASP/5.4.1p1/psxe2015.6; module load ScriptLang/python/2.7\n',
}


"""Local constants"""
PATH2NEBMAKE = '/hdd/home/aksenov/Simulation_wrapper/vts/nebmake.pl'
PATH2POTENTIALS = '~/potcars/'
PATH2JMOL = 'jmol'
PATH2PHONOPY = 'phonopy'
pmgkey = "AWqKPyV8EmTRlf1t" #MAPI_KEY


cluster_tools = 'tools'
show_head = None # show header for res_loop()







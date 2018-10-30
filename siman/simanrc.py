# -*- coding: utf-8 -*-
"""
User-related parameters for siman, file is installed to home folder

"""
from __future__ import division, unicode_literals, absolute_import 


"""Cluster constants"""
DEFAULT_CLUSTER = 'default'
PATH2ARCHIVE = '' # path to archive; if no files are found at home folder, siman will check here; relative paths should be same


CLUSTERS['default'] = {
'address':'aksenov@10.30.16.62', # command for ssh
# 'homepath':'/home/aksenov/',  # deprecated, determined automatically
'schedule':'SLURM',
'corenum':16,
'pythonpath':'/usr/lib64/python2.7/site-packages/numpy',
'vasp_com':'prun /opt/vasp/bin/vasp5.4.1MPI', # path to vasp binary
'modules':'module add prun/1.0; module add intel/16.0.2.181; module add impi/5.1.3.181\n'
}

"""Local constants"""
PATH2POTENTIALS = '~/potcars/'
PATH2JMOL = 'jmol'
PATH2PHONOPY = 'phonopy'
pmgkey = "AWqKPyV8EmTRlf1t" #MAPI_KEY


cluster_tools = 'tools'
show_head = None # show header for res_loop()







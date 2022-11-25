# -*- coding: utf-8 -*-
"""
Default control parameters for siman

"""
from __future__ import division, unicode_literals, absolute_import 


"""Cluster constants"""
CLUSTERS = {}
DEFAULT_CLUSTER = 'slurm'
PATH2PROJECT = '' # path to project on cluster relative to home folder
PATH2ARCHIVE = '' # path to archive; if no files are found at home folder, siman will check here; relative paths should be same


CLUSTERS['slurm'] = {
'address':'man@10.30.16.62', # command for ssh
# 'homepath':'/home/aksenov/',  # deprecated, determined automatically
'schedule':'SLURM',
'corenum':16,
'pythonpath':'/usr/lib64/python2.7/site-packages/numpy',
'vasp_com':'prun /opt/vasp/bin/vasp5.4.1MPI', # path to vasp binary
'modules':'module add prun/1.0; module add intel/16.0.2.181; module add impi/5.1.3.181\n'
}

CLUSTERS['pbs'] = {
'address':'man@10.30.17.10',
# 'homepath':'/home/Dmitry.Aksenov/',
'schedule':'PBS',
'corenum':16,
'pythonpath':'/usr/lib64/python2.7/site-packages/numpy'

}

CLUSTERS['sge'] = {'address':'ut1',
'schedule':'SGE',
'vasp_com':'mpirun  -n $NSLOTS vasp',
'pe':'mpi24',
'shell':'/bin/bash',
'corenum':24,
}


"""Local constants"""
PATH2POTENTIALS = './potcars/'
PATH2NEBMAKE = 'nebmaker.pl'
PATH2JMOL = 'jmol'
PATH2PHONOPY = 'phonopy'
PATH2VASPKIT = 'vaspkit'
pmgkey = "" #PMG_MAPI_KEY

path_to_paper        = './'
PATH2DATABASE        = './'

cluster_tools = 'tools'
geo_folder           = './'
path_to_images       = './'
path_to_wrapper      = '~/Simulation_wrapper/'
show_head = None # show header for res_loop()


"""List of constants determined during installation"""
CIF2CELL = True 






"""List of manually added calculations:"""
MANUALLY_ADDED = [# calc name, calc folder, calc des  
    ( 'Li111'        ,"Li",        "2 Li"                                  ),
    ]




















"""
Naming conventions:

endings:
'_ml' - was used to show that this calculation uses manual equilibrium lattice determination and 
contains several versions of identical structures with different
lattice constants. Now not in use, because I always use this method. Usually 16 versions for hcp;

'_r' - calculation with structure constructed for fitted lattice constants; 
Now was replaced with '.f'; Usually one version.
'.ur' - unrelaxed
.r - relaxed atomic positions
.o - optimised cell and volume and atomic positions automatically
'.f'  - fitted
'.fr' - means that current calculation based on the structure for which lattice constants were fitted and
positions of atoms were relaxed. However see description to know for wich set they were fitted and relaxed.
Calculations with '.f' and '.fr' can have different versions which are correspondig to different sets.

.m - only matrix, all impurities were removed and matrix was freezed


letters in name, wich are usually between didgits and element's names:
b - stands for bulk, which denote ideal cells without boundaries.
g - cells with grain boundary;
v - means that impurity is in the volume of grain; far away from boundaries;
i - means that impurity is close to interface plane (grain boundary)

Versions:
20 - usually means that lattice constatns was used from other calculation and this is very good assumtion.


"""




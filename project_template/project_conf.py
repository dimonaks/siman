# -*- coding: utf-8 -*-
"""
Control of project

TODO:
cif2cell installation check or add to siman
"""
from __future__ import division, unicode_literals, absolute_import 

NEW_BATCH = True   # for testing new batch system based on set sequences 

PBS_PROCS = True # if true than #PBS -l procs="+str(number_cores) is used
WALLTIME_LIMIT = True

"""Cluster constants"""
# cluster_address = 'aksenov@10.30.16.62' #
# CLUSTER_ADDRESS = cluster_address
# cluster_home    = '/home/aksenov/' # needed only for SLURM std out and std err
# CLUSTER_PYTHONPATH = '/usr/lib64/python2.7/site-packages/numpy'

# SCHEDULE_SYSTEM = 'SLURM' #see write_batch_header()
# corenum = 16; #queue = ' -l cmmd '
# CORENUM = corenum

CLUSTERS = {}
DEFAULT_CLUSTER = 'cee'
PATH2PROJECT = '' # path to project on cluster relative to home folder
# project_path_cluster = '' 
# PATH_TO_PROJECT_ON_CLUSTER = project_path_cluster


CLUSTERS['cee'] = {'address':'aksenov@10.30.16.62',
'vasp_com':'prun /opt/vasp/bin/vasp5.4.1MPI_aksenov',
# 'vasp_com':'prun /opt/vasp/bin/vasp5.4.1MPI_ncl',

'homepath':'/home/aksenov/',
'schedule':'SLURM',
'corenum':16,
'pythonpath':'/usr/lib64/python2.7/site-packages/numpy'
}

CLUSTERS['skol'] = {'address':'Dmitry.Aksenov@10.30.17.12',
'vasp_com':'mpirun  vasp_std',
'homepath':'/home/Dmitry.Aksenov/',
'schedule':'PBS',
'walltime':'72:00:00',
'corenum':16,
'pythonpath':'/usr/lib64/python2.7/site-packages/numpy',
'modules':'module load Compilers/Intel/psxe_2015.6; module load MPI/intel/5.1.3.258/intel; module load QCh/VASP/5.4.1p1/psxe2015.6; module load ScriptLang/python/2.7',
}

CLUSTERS['bsu'] = {'address':'aleksenov_d@95.167.109.79',
'vasp_com':'mpiexec --prefix /home/aleksenov_d/mpi/openmpi-1.6.3/installed vasp',
'homepath':'/home/aleksenov_d',
'schedule':'PBS',
'corenum':16,
'pythonpath':'/usr/lib64/python2.7/site-packages/numpy'
}






"""Local constants"""
PATH2POTENTIALS = '/home/aksenov/scientific_projects/PAW_PBE_VASP'
pmgkey = "AWqKPyV8EmTRlf1t"

path_to_paper        = '/home/aksenov/Research/CEStorage/aksenov_report/'
# PATH2DATABASE        = '/home/aksenov/Data/CEStorage/_aksenov'
PATH2DATABASE        = '/home/aksenov/Data/CEStorage/'
# gb4_geo_folder       = '/home/dim/Simulation_wrapper/gb4/out/'
#we have gb5!
PATH2JMOL  = 'jmol'
PATH2NEBMAKE = '~/Simulation_wrapper/vts/nebmake.pl'

geo_folder           = './'
# path_to_images       = '/home/aksenov/ydisk/cathode_report/images/'
path_to_images       = path_to_paper+'/fig/'
# path_to_jmol         = '/home/dim/installed/jmol-14.2.12_2015.02.11/jmol.sh '
path_to_wrapper      = '/home/aksenov/Simulation_wrapper/'

# RAMDISK              = '/mnt/ramdisk/'
RAMDISK              = None
EXCLUDE_NODES  = 1


"""List of constants determined during installation"""
CIF2CELL = True 






"""List of manually added calculations:"""
MANUALLY_ADDED = [# calc name, calc folder, calc des  
    ( 'Li111'        ,"Li",        "2 Li"                                  ),
    ( 'Rb111'        ,"Rb/bcc",        "2 Rb"                                  ),
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




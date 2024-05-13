"""
User-related settings for siman
"""

local_path = '/home/anton/jupyter/siman/'
PATH2POTENTIALS = local_path+'potcars'
PATH2JMOL = 'java -jar Jmol.jar'
AUTO_UPDATE_DB = True
pmgkey = "AWqKPyV8EmTRlf1t" #API_KEY can be generated in the following webpage: https://materialsproject.org/dashboard 


"""Cluster settings"""
DEFAULT_CLUSTER = 'geo' #short name of cluster
user = 'a.boev'

from siman.header import CLUSTERS

CLUSTERS['geo'] = {'address':user+'@95.71.121.195', #cluster address
'vasp_com':'mpirun /opt/vasp/bin/vasp_std', #command for VASP perfoming on cluster
'homepath':'/home/'+user, #path to home directory on cluster
'schedule':'SLURM', #type of schedule system using on cluster
'walltime':'2:00:00', #maximum time for job execution, hours:minutes:seconds, after this time since job was started process will be killed by system
'corenum':2, #number of cores for perfoming of one job on cluster
}

CLUSTERS['raz'] = {'address':user+'@10.16.77.19', #cluster address
'vasp_com':'mpirun vasp_std', #command for VASP perfoming on cluster
'homepath':'/home/'+user, #path to home directory on cluster
'schedule':'SLURM', #type of schedule system using on cluster
'walltime':'2:00:00', #maximum time for job execution, hours:minutes:seconds, after this time since job was started process will be killed by system
'corenum':2, #number of cores for perfoming of one job on cluster
'modules':'module load devtools/compiler/nvhpc/20.11; \
module load q-ch/vasp/5.4.4_OPT; \
\nulimit -s unlimited\n\
'
}



CLUSTERS['cee'] = {'vasp_com':'mpirun vasp_std',
'address':'sk',
'homepath':'/home/a.boev/',
'schedule':'SLURM',
'walltime':'24:00:00',
'corenum':16,
'modules':'module load Compiler/Intel/17u8; module load Q-Ch/VASP/5.4.4_OPT; module load ScriptLang/python/3.6i_2018u3\n ulimit -s unlimited\n',
'partition':'AMG-medium',
}

CLUSTERS['cee-omc'] = {'vasp_com':'mpirun vasp_std',
'address':'sk',
'homepath':'/home/a.boev/',
'schedule':'SLURM',
'walltime':'72:00:00',
'corenum':16,
'modules':'module load Compiler/Intel/17u8; module load Q-Ch/VASP/5.4.4_OMC; module load ScriptLang/python/3.6i_2018u3\n ulimit -s unlimited\n',
'partition':'AMG-medium',
}


CLUSTERS['cee-pcm'] = {'vasp_com':'mpirun vasp_std',
'address':'sk',
'homepath':'/home/a.boev/',
'schedule':'SLURM',
'walltime':'72:00:00',
'corenum':16,
'modules':'module load Compiler/Intel/17u8; module load Q-Ch/VASP/5.4.4_SOL; module load ScriptLang/python/3.6i_2018u3\n ulimit -s unlimited\n',
'partition':'AMG-medium',
}





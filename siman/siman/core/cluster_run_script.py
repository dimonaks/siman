# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
from siman import header
from siman.header import printlog
from siman.functions import push_to_server, run_on_server

def prepare_run():
    """
    INPUT:
        schedule_system - type of job scheduling system:'PBS', 'SGE', 'SLURM'
    """
    schedule_system = header.schedule_system
    with open('run','w', newline = '') as f:
    
        if schedule_system == 'SGE':
            if 'shell' in header.cluster:
                f.write("#!"+header.cluster['shell']+'\n')

            else:
                f.write("#!/bin/tcsh\n")
            # if 'modules' in header.cluster:
                # f.write(header.cluster['modules']+'\n')
            # f.write("module load sge\n")
            # f.write("module load vasp/parallel/5.2.12\n")
        elif schedule_system in ('PBS', 'SLURM', 'none'):
            f.write("#!/bin/bash\n")
        else:
            ''
            # print_and_log('Please provide schedule_system!')
            # raise RuntimeError

    # f.close()
    return




def make_run(cl_dir, schedule_system, run_name):
    """
    Generate run file

    INPUT:
        - cl_dir (str) - calculation directory cl.dir
        - schedule_system (str) - ['SGE', 'PBS', 'SLURM', 'none']
    """

    with open('run','a', newline = '') as f:

        if schedule_system == 'SGE':
            #'qsub -pe 'mpi*' NCORES -l CLUSTER_TAG script.parallel.sh' for mpi-jobs which should run on CLUSTER_TAG (cmmd or cmdft)
            #IMPORTANT: NCORES must be a multiple of 8(24) on cmdft(cmmd). 
            # f.write("qsub -pe 'mpi*' "+str(header.corenum)+" "+header.queue+" "+run_name+"\n") #str(self.set.np) #-l cmmd; on MPIE
            
            f.write("qsub "+" "+run_name+"\n") 
        
            # f.write('sleep 5\n')
            # runBash('chmod +x run')
        
        elif schedule_system in ['PBS']:
            if header.PATH2PROJECT == '':
                header.PATH2PROJECT = '.'

            f.write("cd "+header.PATH2PROJECT+'/'+cl_dir+"\n")
            f.write("qsub "+run_name.split('/')[-1]+"\n") 
            f.write("cd -\n")
            f.write('sleep 1\n')                        
        
        
        elif schedule_system == 'SLURM':
            f.write("squeue\n") 
            f.write("sbatch " + run_name+"\n") 
            # f.write("sbatch -p AMG " + run_name+"\n") 

        elif schedule_system in ['none']:
            if header.PATH2PROJECT == '':
                header.PATH2PROJECT = '.'

            f.write("cd "+header.PATH2PROJECT+'/'+cl_dir+"\n")
            f.write('./'+run_name.split('/')[-1]+"\n") 
            f.write("cd -\n")
            # f.write('sleep 1\n')     


        else:
            printlog('Error! Unknown schedule_system', schedule_system)
            



    printlog("\nRun file created\n")     
    return










def complete_run(close_run = True):
    header.first_run = False
    if close_run:

        with open('run','a', newline = '') as f:
            if header.schedule_system in ["PBS", 'PBS_bsu']:
                f.write("qstat\n")
                f.write("sleep 2\n")
            elif header.schedule_system == "SLURM":
                f.write("squeue\n")
            elif header.schedule_system == "SGE":
                f.write("qstat\n")
            elif header.schedule_system == "none":
                f.write("\n")
            f.write("mv run last_run\n")


        if header.copy_to_cluster_flag:
            push_to_server('run',  header.cluster_home, header.cluster_address)
            run_on_server('chmod +x run', header.cluster_address)
            printlog('run sent')
    
    return




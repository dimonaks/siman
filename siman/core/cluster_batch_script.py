# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
import os
from siman.small_functions import list2string
from siman import header
from siman.header import printlog

def write_batch_header(cl, batch_script_filename = None,
    schedule_system = None, path_to_job = None, job_name = 'SuperJob', corenum = None):
    """
    Write header for different schedule systems
    path_to_job (str) - absolute path to job folder 
    """
    NC = str(corenum)
    with open(batch_script_filename,'w', newline = '') as f:


        if schedule_system == 'SGE':
            if 'shell' in header.cluster:
                def_shell = header.cluster['shell']
            else:
                def_shell = '/bin/tcsh'
            
            f.write("#!"+def_shell+"\n")
            # f.write("#$ -M aksenov@mpie.de\n")
            f.write("#$ -m be\n")
            f.write("#$ -S "+def_shell+"\n")
            f.write("#$ -cwd \n")
            f.write("#$ -R y \n")
            f.write("#$ -V  \n") # use variables
            f.write("#$ -o "+path_to_job+" -j y\n")
            if 'pe' in header.cluster:
                f.write("#$ -pe "+header.cluster['pe']+' '+NC+'\n')
            f.write("\n")


            f.write("cd "+path_to_job+"\n")

            if 'modules' in header.cluster:
                f.write(header.cluster['modules']+'\n')
            # f.write("module load sge\n")
            # f.write("module load vasp/parallel/5.2.12\n\n")


        if schedule_system == 'PBS':

            prefix = '#PBS'
            f.write("#!/bin/bash   \n")
            f.write("#PBS -N "+job_name+"\n")
            
            if 'walltime' in header.cluster:
                f.write("#PBS -l walltime="+str(header.cluster['walltime'])+'\n')
            else:
                if header.WALLTIME_LIMIT: #deprecated remove
                    f.write("#PBS -l walltime=72:00:00 \n")
            

            # f.write("#PBS -l nodes=1:ppn="+str(corenum)+"\n")
            nodes = 1
            if 'nodes' in header.cluster:
                nodes = header.cluster['nodes']



            if header.PBS_PROCS or nodes == 0: # parameter in header
                f.write("#PBS -l procs="+str(corenum)+"\n")
            else: #  node option 

                if type(nodes) is str:
                    f.write("#PBS -l nodes="+nodes+"\n")
                else:
                    f.write("#PBS -l nodes="+str(nodes)+":ppn="+str(corenum)+"\n")
            


            if 'queue' in header.cluster:
                f.write("#PBS -q "+str(header.cluster['queue'])+'\n')


            
            if 'procmemgb' in header.cluster:
                f.write("#PBS -l pmem="+str(header.cluster['procmemgb'])+'gb\n')

            if 'feature' in header.cluster:
                f.write("#PBS -l feature="+header.cluster['feature']+"\n")



            f.write("#PBS -r n\n")
            f.write("#PBS -j eo\n")
            f.write("#PBS -m bea\n")

            if 'any_commands' in header.cluster:
                lines = header.cluster['any_commands']
                if not is_list_like(lines):
                    printlog('Error! Please use list for sbatch key in cluster description')
                for line in lines:
                    f.write(prefix+' '+line+'\n')



            f.write("cd $PBS_O_WORKDIR\n")
            

            if 'modules' in header.cluster:
                f.write(header.cluster['modules']+'\n')



        if schedule_system == 'PBS_bsu':
            f.write("#!/bin/bash   \n")
            f.write("#PBS -N "+job_name+"\n")
            if header.WALLTIME_LIMIT:
                f.write("#PBS -l walltime=72:00:00 \n")
            # f.write("#PBS -l nodes=1:ppn="+str(corenum)+"\n")
            if header.PBS_PROCS:
                f.write("#PBS -l nodes=node07:ppn="+str(corenum)+"\n")
            else: # 1 node option 
                f.write("#PBS -l nodes=node07:ppn="+str(corenum)+"\n")
            # f.write("#PBS -l pmem=16gb\n") #memory per processor, Skoltech
            f.write("#PBS -r n\n")
            f.write("#PBS -j eo\n")
            f.write("#PBS -m bea\n")
            f.write("#PBS -M boev.anton.olegovich@gmail.com\n")
            f.write("cd $PBS_O_WORKDIR\n")
            f.write("echo $LD_LIBRARY_PATH \n")




        if schedule_system == 'none':
            f.write("#!/bin/bash   \n")
            if 'modules' in header.cluster:
                f.write(header.cluster['modules']+'\n')



        if schedule_system == 'SLURM':
            
            hc = header.cluster
            if '~' in path_to_job:
                print_and_log('Error! For slurm std err and out you need full paths')

            f.write("#!/bin/bash   \n")
            f.write("#SBATCH -J "+job_name+"\n")
            if 'walltime' in header.cluster:
                f.write("#SBATCH -t "+str(header.cluster['walltime'])+'\n')
            else:
                ''
                # f.write("#SBATCH -t 250:00:00 \n")

            f.write("#SBATCH -N 1\n")
            f.write("#SBATCH -n "+str(corenum)+"\n")
            f.write("#SBATCH -o "+path_to_job+"sbatch.out\n")
            f.write("#SBATCH -e "+path_to_job+"sbatch.err\n")
            if header.MEM_CPU:
                f.write("#SBATCH --mem-per-cpu=7675\n") # this is mem per core
            
            # print(header.cluster)
            # sys.exit()
            if 'partition' in hc:
                f.write('#SBATCH -p '+hc['partition']+'\n')

            if 'any_commands' in header.cluster:
                lines = header.cluster['any_commands']
                if not is_list_like(lines):
                    printlog('Error! Please use list for sbatch key in cluster description')
                for line in lines:
                    f.write('#SBATCH '+line+'\n')

            # f.write("#SBATCH -I other=avx\n") # AVX2 instructions for new node to improve speed by 18% 

            # f.write("#SBATCH --nodelist=node-amg03\n")
            if header.siman_run: #only for me
                if header.EXCLUDE_NODES:
                    f.write("#SBATCH --exclude=node-amg11\n")
                # f.write("#SBATCH --mail-user=\n")
                # f.write("#SBATCH --mail-type=END\n")
            f.write("cd "+path_to_job+"\n")
            # f.write("export OMP_NUM_THREADS=1\n")

            if 'modules' in header.cluster:
                f.write(header.cluster['modules']+'\n')
            
           
            f.write("export PATH=$PATH:"+header.cluster['homepath'] +"/tools/\n")

        
        if schedule_system == 'simple':
            f.write("#!/bin/bash   \n")
            # f.write("mpirun -np "+str(corenum)+" /home/hieuvatly/vasp.5.4.4/bin/vasp_std \n")
            print('"write_batch_header" was launched successfully!')  


        if cl.calculator == 'gaussian':
            f.write('export GAUSS_SCRDIR=/scr/$SLURM_JOB_USER/$SLURM_JOB_ID\n')


        f.write("touch RUNNING\n")




    return









def prepare_input(cl, prevcalcver = None, option = None, input_geofile = None, name_mod_prev = '', write = True, 
    curver = None, copy_poscar_flag = True, f = None ):
    """1. Input files preparation

        curver - current version
    """  


    if write:
        # if not 'only_neb' in cl.calc_method:
        precont = str(prevcalcver)+name_mod_prev+'.CONTCAR' #previous contcar
        if option == 'inherit_xred' and prevcalcver:
            if copy_poscar_flag:

                f.write('grep -A '+str(cl.init.natom)+ ' "Direct" '+precont+' >> '+input_geofile+ ' \n')

        if copy_poscar_flag: # only for first set 
            if option == 'continue': #test for the case of sequence set - OK
                ''
                precont = str(curver)+name_mod_prev+'.CONTCAR ' #previous contcar
                preout  = str(curver)+name_mod_prev+'.OUTCAR ' #previous outcar
                f.write("cp "+precont+" POSCAR  # inherit_option = continue\n")
                f.write("cp "+preout+'prev.'+preout+" # inherit_option = continue\n")
                f.write('mv CHGCAR prev.CHGCAR   # inherit_option = continue\n')
            
            else:
                if cl.calculator == 'vasp':
                    f.write("cp "+input_geofile+" POSCAR\n")
                elif cl.calculator == 'gaussian':
                    pass




    return





def run_command(cl, option, name, parrallel_run_command,
    write = True, mpi = None, corenum = None, f = None ):
    """
    Write commands used for running calculator. 

    INPUT:

        - cl (Calculation) - 
        - mpi (bool) - write mpirun command with explicit number of cores 

    TODO:
        improve mpi, cores 

    """ 

    if write:

        if mpi is False:

            if option == 'master':
                f.write("vasp >"+name+".log\n")

            elif 'monte' in cl.calc_method:
                f.write("python "+header.cluster_home+'/'+ header.cluster_tools+'/siman/monte.py > monte.log\n')

            elif 'polaron' in cl.calc_method:
                f.write("python "+header.cluster_home+'/'+ header.cluster_tools+'/siman/polaron.py > polaron.log\n')
            elif 'polaron2' in cl.calc_method:
                f.write("python "+header.cluster_home+'/'+ header.cluster_tools+'/siman/polaron_mod.py > polaron.log\n')

            elif 'atat' in  cl.calc_method:
                f.write('maps -d&\npollmach runstruct_vasp mpirun\n')

            else:
                if cl.calculator == 'vasp':
                    f.write(parrallel_run_command +" >"+name+".log\n")
                elif cl.calculator == 'gaussian':
                    f.write(parrallel_run_command +" < input.gau > "+name+".out\n")
                else:
                    printlog('Error! Calculator ', cl.calculator, 'is unknown!')

        elif mpi == True:
            f.write('mpirun -np '+str(cores)+' '+parrallel_run_command+" >"+name+".log\n")


        f.write("sleep 20\n")
    return


def mv_files_according_versions(cl, savefile, v, name_mod = '', write = True, 
    rm_chg_wav = 'cw', f = None    ):    
    """3. Out files saving block
        
        rm_chg_wav - if True than CHGCAR and WAVECAR are removed

        savefile (str) - key, which determines what files should be saved
            'o' - OUTCAR
            'i' - INCAR
            'v' - CHG
            'c' - CHGCAR
            'p' - PARCHG
            'l' - LOCPOT
            'd' - DOSCAR
            'a' - AECCAR0, AECCAR2
            'x' - vasprun.xml
            't' - XDATCAR
            'z' - OSZICAR
            'w' - WAVECAR

    """   
    printlog('The value of savefile is', savefile)
    


    if 'polaron' in cl.calc_method:
        write = 0 # not needed, since files are automatically saved by python script

    pre = v + name_mod
    
    if cl.calculator == 'vasp':
        final_structure_file = pre+'.CONTCAR'

        try:
            printlog('The files to be removed are ', header.clean_vasp_files)
        except AttributeError:
            raise RuntimeError('The variable "clean_vasp_files" is absent! It should be initially stated in project_conf.py file as a list of file names ["NAME"]!!!')

        try:
            printlog('The files to be removed are ', header.clean_vasp_files_ignore)
        except AttributeError:
            raise RuntimeError('The variable "clean_vasp_files_ignore" is absent! It should be initially stated in project_conf.py file as a list of file names ["NAME"]!!!')

        if write:

            files_key_dict = {"o": "OUTCAR",
                              "s": "CONTCAR",
                              # "i": "INCAR", #bad idea for INCAR
                              "e": "EIGENVAL",
                              "v": "CHG",
                              "c": "CHGCAR",
                              "p": "PARCHG",
                              "r": "PROCAR",
                              "l": "LOCPOT",
                              "d": "DOSCAR",
                              "a0": "AECCAR0",
                              "a2": "AECCAR2",
                              "x": "vasprun.xml",
                              "t": "XDATCAR",
                              "z": "OSZICAR",
                              "w": "WAVECAR",
                              "f": "WAVEDER"}

            if 'a' in savefile:
                savefile+='a0a2' # for compatibility with previous behavior

            if header.SAVE_CONTCAR:
                savefile+='s'


            for i in files_key_dict.keys():
                if i in savefile:
                    fln = files_key_dict[i]
                    prefln = pre+'.'+fln
                    if i+'!' in savefile:
                        f.write('cp '+fln+' '+prefln+'\n')
                    else:
                        f.write('mv '+fln+' '+prefln+'\n')
                    if (i+'!@' in savefile) or (i+'@' in savefile): 
                        f.write("gzip -f "+prefln+"\n")
                    else:
                        pass
                else:
                    if files_key_dict[i] in header.clean_vasp_files_ignore: 
                        pass
                    else: 
                        header.clean_vasp_files.append(files_key_dict[i])                                                   



            if 0:
                'obsolete'
                if 'c' in rm_chg_wav:
                    f.write("rm CHGCAR   # rm_chg_wav flag\n") #file is important for continuation
                if 'w' in rm_chg_wav:
                    ''
                    f.write("rm WAVECAR  # rm_chg_wav flag\n") #
                if 'v' in rm_chg_wav: #chgcar for visualization
                    ''
                    f.write("rm CHG   # rm_chg_wav flag\n") #

    elif cl.calculator == 'gaussian':
        final_structure_file = None
        pass




    return final_structure_file


def analysis_script(cl, write = True, f = None ):
    #now only for u-ramping
    if write:
        f.write("touch ENERGIES\n")                

        for outcar in cl.associated_outcars:
            f.write("grep 'energy  without entropy' "+outcar+" | awk '{print $7}' >> ENERGIES\n")

    return


def name_mod_U_last(cl):
    name_mod_last = 'U'+str(
                cl.update_incar(parameter = 'LDAUU', 
                    u_ramp_step = cl.set.u_ramping_nstep-1, write = False, f = f)).replace('.','') #used to det last U

    return name_mod_last



 
def write_body(cl, version = None, savefile = None, set_mod = '', copy_poscar_flag = True,
    final_analysis_flag = True, penult_set_name = None, 
    curset = None, f = None, prevcalcver = None, option = None,
    input_geofile = None, parrallel_run_command = None, mpi = False, corenum = 1):
    """
    cl (Calculation) - calculation object
    v (int) - version
    set_mod (str) - additional modification of names needed for *set_sequence* regime, should be '.setname'
    
    f (file ) - file object
    """
    v = version
    if 'only_neb' in cl.calc_method:
        write = False
        write_poscar = False
    else:
        write = True
        write_poscar = True                

    #neb
    if 'neb' in cl.calc_method: 
        
        if write: 
            f.write("#NEB run, start and final configurations, then IMAGES:\n") 
        
        cl.update_incar(parameter = 'IMAGES', value = 0, write = write, f = f) # start and final runs


    if 'u_ramping' in cl.calc_method:

        if write: 
            f.write("#U-ramping run:\n")  

        # name_mod_last = '.'+name_mod_U_last()
        name_mod_last = '.U00' #since u_ramp starts from u = 00, it is more correct to continue from 00
        if penult_set_name:
            name_mod_last += '.'+penult_set_name #however, for multiset run, the structure for u=00 exists only
                                              #for penult set or maybe even for the first only set
            # print (name_mod_last, penult_set_name)
            # sys.exit()
            
        # print 'prevcalver', prevcalcver

        if write and copy_poscar_flag: 
            f.write("rm CHGCAR    #u-ramp init from scratch\n")                

        prepare_input(cl, prevcalcver = prevcalcver, option = option,
         input_geofile = input_geofile, name_mod_prev = name_mod_last, write = write_poscar, curver = version,
         copy_poscar_flag = copy_poscar_flag, f = f)

        if copy_poscar_flag:
            usteps = range(cl.set.u_ramping_nstep)
        else:
            usteps = [cl.set.u_ramping_nstep-1]  # now it the case of sequence_set for contin sets only the last U is used

        u_last = 100
        for i_u in usteps:

            u = cl.update_incar(parameter = 'LDAUU', u_ramp_step = i_u, write = write, f = f)
            if u == u_last:
                continue
            name_mod   = '.U'+str(u).replace('.', '')+set_mod
           
            run_command(cl, option = option, name = cl.name+name_mod, 
                parrallel_run_command = parrallel_run_command, write = write, mpi = mpi, corenum = corenum, f = f)

            if write: 
                if copy_poscar_flag:

                    f.write("cp CONTCAR POSCAR   #u-ramp preparation\n")                

    # print(savefile)
            contcar_file = mv_files_according_versions(cl, savefile, v, 
                name_mod = name_mod, write = write, rm_chg_wav = '', f = f)
        


            cl.associated_outcars.append( v + name_mod +  ".OUTCAR"  )
            u_last = u
        
        if final_analysis_flag:
            rm_chg_wav = 'w' #The wavcar is removed for the sake of harddrive space
        
        else:
            rm_chg_wav = ''

        if curset.save_last_wave:
            save_last = 'cw'
        else:
            save_last = 'c'


        mv_files_according_versions(cl, savefile = save_last, v=version, 
        name_mod = name_mod, rm_chg_wav = rm_chg_wav, f = f) #save more files for last U
        

        analysis_script(cl, write = write, f=f)
        # print cl.associated





    elif 'afm_ordering' in cl.calc_method:

        #Comment - inherit_xred option is not available here
        f.write("rm CHGCAR\n")                
        if not savefile: 
            savefile = 'o'

        for i, magmom in enumerate(cl.magnetic_orderings):

            name_mod   = '.AFM'+str(i)+set_mod

            cl.update_incar(parameter = 'MAGMOM', value = magmom, write = write, f = f)

            prepare_input(cl, prevcalcver = prevcalcver, option = option, input_geofile = input_geofile,
                copy_poscar_flag = copy_poscar_flag, f = f)
            
            run_command(cl, option = option, name = cl.name+name_mod, 
                parrallel_run_command = parrallel_run_command, mpi = mpi, corenum = corenum,  f = f)

            contcar_file = mv_files_according_versions(cl, savefile, version, name_mod = name_mod, f = f)
        
            cl.associated_outcars.append( v + name_mod +  ".OUTCAR"  )

        analysis_script(cl, write = write, f=f)
    

    else: #simple run
        
        if not savefile: 
            savefile = 'vco'

        if write: 
            f.write("#Basic run:\n")  

        name_mod   = set_mod
        name_mod_last = ''

        prepare_input(cl, prevcalcver = prevcalcver, option = option, name_mod_prev = name_mod_last,
            input_geofile = input_geofile, write = write_poscar, curver = version,
            copy_poscar_flag = copy_poscar_flag, f = f)

        run_command(cl = cl, option = option, name = cl.name+name_mod, 
            parrallel_run_command = parrallel_run_command, write = write, mpi = mpi, corenum = corenum, f = f)

        if final_analysis_flag:
            rm_chg_wav = 'w' #The wavcar is removed for the sake of harddrive space
        
        else:
            rm_chg_wav = ''


        contcar_file = mv_files_according_versions(cl, savefile, version, 
            write = write, name_mod = name_mod, rm_chg_wav = rm_chg_wav, f = f)

        cl.associated_outcars.append( version + name_mod +  ".OUTCAR"  )

    return contcar_file


def u_ramp_prepare(cl):
    if 'u_ramping' in cl.calc_method:
        u = cl.update_incar(parameter = 'LDAUU', u_ramp_step = cl.set.u_ramping_nstep-1, write = False, f = f)
        name_mod   = '.U'+str(u).replace('.', '')
        # name_mod_last = name_mod_U_last()+'.'
        name_mod_last = '.'+'U00' #since u_ramp starts from u = 00, it is more correct to continue from 00
    
    else:
        name_mod_last = ''
        name_mod   = ''                

    return name_mod, name_mod_last

def u_ramp_loop(cl, ver_prefix = '', subfolders = None, run_name_prefix = None, set_mod = '', f = None):

    if not subfolders:
        subfolders = ['.']

    

    if run_tool_flag:
        usteps = range(cl.set.u_ramping_nstep)
    else:
        usteps = [cl.set.u_ramping_nstep-1]  # now it the case of sequence_set for contin sets only the last U is used


    u_last = 100

    for i_u in usteps:


        u = cl.update_incar(parameter = 'LDAUU', u_ramp_step = i_u, write = 1,  f = f)
        if u == u_last:
            continue

        name_mod   = ver_prefix+'U'+str(u).replace('.', '')+set_mod

        
        run_command(cl, option = option, name = run_name_prefix+'.'+name_mod, 
            parrallel_run_command = parrallel_run_command, write = True, mpi = mpi, corenum = corenum)
        
        u_last = u


        for n_st in subfolders:

            f.write('cp '+n_st+'/CONTCAR '+n_st+'/POSCAR'+'               #u_ramp_loop()\n' )
            f.write('cp '+n_st+'/OUTCAR  '+n_st+'/'+name_mod+'.OUTCAR'+'  #u_ramp_loop()\n' )
            contcar = name_mod+'.CONTCAR'
            f.write('cp '+n_st+'/CONTCAR  '+n_st+'/'+contcar+'            #u_ramp_loop()\n' )

            # cl.associated_outcars.append( v + name_mod +  ".OUTCAR"  )


    return contcar



def write_footer(cl, set_mod = '', run_tool_flag = True, 
    savefile = None, final_analysis_flag = True, neb_flag = None, f = None, mpi = False, corenum = 1, option = None,parrallel_run_command=None, output_files_names = None):
    """footer"""
    

    subfolders = None
    contcar_file = None
    
    if neb_flag:
        printlog('Writing scripts for NEB method', important = 'n')
        nim = cl.set.vasp_params['IMAGES']
        nim_str = str(nim)

        subfolders = []
        for n in range(1, nim+1):
            if n < 10:
                n_st = '0'+str(n)
            elif n < 100:
                n_st = str(n)
            subfolders.append(n_st)


        name_mod, name_mod_last = u_ramp_prepare(cl)


        start = '1'+name_mod+'.OUTCAR '
        final = '2'+name_mod+'.OUTCAR '
        startC = start.replace('OUT','CONT')
        finalC = final.replace('OUT','CONT')
        start_folder = '00'
        if nim+1 < 10: 
            final_folder = '0'+str(nim+1)
        else:
            final_folder = str(nim+1)



        f.write("\n\n#Starting NEB script \n")

        if option and 'continue' in option:
            prevout = name_mod_last+'OUTCAR '

            for n_st in subfolders:
                f.write('cp '+n_st+'/'+prevout+n_st+'/'+'prev.'+prevout+'  # inherit_option = continue\n' )
                f.write('cp '+n_st+'/CONTCAR '+n_st+'/POSCAR  # inherit_option = continue\n')
                f.write('mv '+n_st+'/CHGCAR '+n_st+'/prev.CHGCAR   # inherit_option = continue\n')

        else:
            
            if run_tool_flag:
                if hasattr(cl, 'init_neb_geo_fld') and cl.init_neb_geo_fld:
                    #only 1.CONTCAR and 2.CONTCAR should be copied to start and end images
                    ff = start_folder 
                    lf = final_folder
                    f.write('cp '+ff+'/POSCAR '+ff+'/POSCAR_init\n')
                    f.write('cp '+lf+'/POSCAR '+lf+'/POSCAR_init\n')
                    f.write('cp '+startC+ff+'/POSCAR\n')
                    f.write('cp '+finalC+lf+'/POSCAR\n')
                else:
                    f.write('export PATH=$PATH:'+header.cluster_home+'/tools/vts/\n') #header.project_path_cluster

                    f.write('nebmake.pl '+ startC + finalC + nim_str +' \n')


        if run_tool_flag:
            f.write('cp '+start +  '00/OUTCAR\n')
            f.write('cp '+final +  final_folder + '/OUTCAR\n' )


        cl.update_incar(parameter = 'IMAGES', value = nim, write  =1, f  = f)


        if 'u_ramping' in cl.calc_method:


            contcar_file = u_ramp_loop(cl, subfolders = subfolders, 
                run_name_prefix = cl.name+'.n_'+nim_str, 
                set_mod = set_mod, f = f)
      

        else:

            run_command(cl, option = option, name = cl.name+set_mod+'.n_'+nim_str+name_mod, 
            parrallel_run_command = parrallel_run_command, write = True, mpi = mpi, corenum = corenum, f=f)
            # print(set_mod)
            # sys.exit()
            if '.' in set_mod and set_mod[0] == '.':
                set_mod_loc = set_mod[1:]
            else:
                set_mod_loc = set_mod

            name_mod   = set_mod_loc
            if name_mod:
                contcar = name_mod+'.CONTCAR'
                outcar  = name_mod+'.OUTCAR'
                for n_st in subfolders:
                    f.write('cp '+n_st+'/OUTCAR  '+n_st+'/'+outcar  +'  #sequence set: save file\n' )
                    f.write('cp '+n_st+'/CONTCAR  '+n_st+'/'+contcar+'  #sequence set: save file\n' )
            else:
                contcar = 'CONTCAR'

            contcar_file = contcar






        if final_analysis_flag:
            # f.write('export PATH=$PATH:'+header.cluster_home+'/tools/gnuplot/bin/ \n')
            # f.write(header.cluster_home+'/tools/vts/nebresults.pl  \n')
            f.write('find . -name WAVECAR -delete\n')
            f.write('find . -name PROCAR -delete\n')
        # for n in range



    # print (calc[id].calc_method )
    # sys.exit()
    if 'uniform_scale' in cl.calc_method or 'c_scale' in cl.calc_method:
        # print (input_geofile)
        name_mod = set_mod
        
        if run_tool_flag:
            f.write("\n\n#Starting fitting tool \n")
            outputs = [ os.path.basename(out) for out in output_files_names ]
            # f.write('export PYTHONPATH=$PYTHONPATH:'+CLUSTER_PYTHONPATH+'\n')
            # f.write('/home/aksenov/tools/fit_tool.py '+list2string(outputs)+'\n' )
            f.write('python '+header.cluster_home+'/tools/fit_tool.py '+list2string(outputs)+'\n' )
            

            f.write('cp 100.POSCAR POSCAR \n')
        
        if 'u_ramping' in cl.calc_method:
            

            contcar_file = u_ramp_loop(cl, ver_prefix = '100.', 
                run_name_prefix = cl.id[0]+'.fitted', set_mod = set_mod, f = f)

        else:
            if final_analysis_flag:
                rm_chg_wav = 'w' #The wavcar is removed for the sake of harddrive space
            
            else:
                rm_chg_wav = ''




            run_command(cl, option = option, name = cl.id[0]+'.'+cl.id[1]+'.100'+name_mod+'.fitted', 
                parrallel_run_command = parrallel_run_command, write = True, mpi = mpi, corenum = corenum, f = f)

            # print(final_analysis_flag)
            # sys.exit()

            contcar_file = mv_files_according_versions(cl, savefile, '100', 
                name_mod = name_mod, write = True, rm_chg_wav = rm_chg_wav, f = f)

        # sys.exit()


    #clean at the end
    if final_analysis_flag: 
        if header.final_vasp_clean:
            # print(header.clean_vasp_files)
            uni_files = list(set(header.clean_vasp_files))
            uni_files_string = ''
            for i in uni_files:
                uni_files_string += i+' '
            # print(uni_files_string)
            if uni_files_string:
                f.write('rm '+uni_files_string+'\n')
            # f.write('rm LOCPOT CHGCAR CHG PROCAR DOSCAR OSZICAR PCDAT REPORT XDATCAR vasprun.xml\n')
        



        f.write('rm RUNNING\n')



    return contcar_file, subfolders







def write_batch_body(cl, input_geofile = "header", version = 1, option = None, 
    prevcalcver = None, savefile = None, schedule_system = None,
    output_files_names = None,
    mode = None,
    batch_script_filename = None, mpi = False, corenum = None):
    """
    Create job script for cluster for different calculation types, such as volume scan, neb, and so on.
    without arguments writes header, otherwise appends sequence of calculations.

    INPUT:
        - option (str) - the same as inherit_option: 
            - 'inherit_xred' - control inheritance, or 
            - 'master'       - run serial on master 
        - prevcalcver (int) - version of previous calc; for the first version equal to None
        - savefile (str) - 'cdawx', where c-charge, d-dos, a- AECCAR, w-wavefile, x-xml
        - schedule_system (str) - type of job scheduling system:
            - 'PBS' 
            - 'SGE' 
            - 'SLURM' 
            - 'none' - just run without any system
        - mode (str)
            - 'body'
            - 'footer'
    RETURN:
        Name of output file
    """

    self = cl
    varset = header.varset
    
    f = open(batch_script_filename, 'a', newline = '')

    printlog("The following output files will be saved: savefile =", savefile,)

    if self.calculator == 'vasp':
        parrallel_run_command = header.vasp_command

    elif self.calculator == 'gaussian':
        parrallel_run_command = self.cluster.get('gaussian_command') or 'please provide command for gaussian in *cluster* dict in simanrc.py'

    else:
        printlog('Error! Unknown calculator', self.calculator )


    run_name = batch_script_filename     
    job_name = self.id[0]+"."+self.id[1]
    neb_flag = ('neb' in self.calc_method or 'only_neb' in self.calc_method)

    if hasattr(self.set, 'set_sequence') and self.set.set_sequence and any(self.set.set_sequence):
        sets = [self.set]+[se for se in self.set.set_sequence]
    else:
        sets = [self.set]


    nsets = len(sets)
    footer_flag = not set(self.calc_method).isdisjoint(['uniform_scale', 'neb', 'only_neb' ])


    if mode == "body": #control part of script
        self.associated_outcars = []

    penult_set_name = None

    for k, curset in enumerate(sets):
        
        if nsets > 1: #the incar name is modified during creation only if more than 1 set is detected
            if mode == 'body' or footer_flag:
                f.write('\n#sequence set: '+curset.ise+' \n')
                f.write('cp '+curset.ise+'.INCAR  INCAR\n')
                if hasattr(curset, 'savefile') and len(curset.savefile) > 0:
                    savefile = curset.savefile 


            penult_set_name = sets[-2].ise
        

        if k < nsets-1:
            set_mod = '.'+curset.ise
            final_analysis_flag = False
        else: #last set
            set_mod = '' # the last step do not use modifications of names 
            final_analysis_flag = True #for footer



        if k == 0: # additional control of prepare_input routine and footer
            copy_poscar_flag = True # the flag is also used to detect first set
            run_tool_flag = True
        else:
            copy_poscar_flag = False
            run_tool_flag  = False

        if mode == "body":
            
            contcar_file = write_body(self, version = str(version), savefile = savefile, 
                set_mod = set_mod, copy_poscar_flag = copy_poscar_flag, 
                final_analysis_flag = final_analysis_flag, 
                penult_set_name = penult_set_name, curset = curset, 
                f = f, prevcalcver = prevcalcver, option = option,
                input_geofile = input_geofile, 
                parrallel_run_command = parrallel_run_command, mpi = mpi, corenum = corenum)

            

        elif mode == 'footer':
            if copy_poscar_flag: 
                f.write('\n#Footer section: \n')


            # print(savefile)
            # sys.exit()
            contcar_file, subfolders = write_footer(self, set_mod = set_mod, run_tool_flag = run_tool_flag, savefile = savefile,
             final_analysis_flag = final_analysis_flag, neb_flag = neb_flag, f = f, mpi = mpi, corenum = corenum, option = option, parrallel_run_command=parrallel_run_command, output_files_names = output_files_names)

        
        if k < nsets-1 and contcar_file:
            if 'o' in savefile:
                if neb_flag and mode == 'footer':
                    for n_st in subfolders:
                        f.write('cp '+n_st+'/'+contcar_file+' '+n_st+'/POSCAR  # sequence_set: preparation of input geo for next neb set\n' )
                else:
                    f.write('cp '+contcar_file+' POSCAR  #sequence_set: preparation of input geo for next set\n')


    if hasattr(self, 'associated_outcars') and  self.associated_outcars:
        out = self.associated_outcars[-1]
    else:
        out = None
    # print 'write_sge() out=', out
    f.close()
    
    return  out#return OUTCAR name
    

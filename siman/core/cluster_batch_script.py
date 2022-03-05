# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
from siman import header
from siman.header import printlog

def write_batch_header(batch_script_filename = None,
    schedule_system = None, path_to_job = None, job_name = 'SuperJob', number_cores = 1):
    """
    cl-explanatory)
    path_to_job (str) - absolute path to job folder 
    """
    NC = str(number_cores)
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
            

            # f.write("#PBS -l nodes=1:ppn="+str(number_cores)+"\n")
            nodes = 1
            if 'nodes' in header.cluster:
                nodes = header.cluster['nodes']



            if header.PBS_PROCS or nodes == 0: # parameter in header
                f.write("#PBS -l procs="+str(number_cores)+"\n")
            else: #  node option 

                if type(nodes) is str:
                    f.write("#PBS -l nodes="+nodes+"\n")
                else:
                    f.write("#PBS -l nodes="+str(nodes)+":ppn="+str(number_cores)+"\n")
            


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
            # f.write("#PBS -l nodes=1:ppn="+str(number_cores)+"\n")
            if header.PBS_PROCS:
                f.write("#PBS -l nodes=node07:ppn="+str(number_cores)+"\n")
            else: # 1 node option 
                f.write("#PBS -l nodes=node07:ppn="+str(number_cores)+"\n")
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
            f.write("#SBATCH -n "+str(number_cores)+"\n")
            f.write("#SBATCH -o "+path_to_job+"sbatch.out\n")
            f.write("#SBATCH -e "+path_to_job+"sbatch.err\n")
            if header.MEM_CPU:
                f.write("#SBATCH --mem-per-cpu=7675\n")
            
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
                # f.write("#SBATCH --mail-user=d.aksenov@skoltech.ru\n")
                # f.write("#SBATCH --mail-type=END\n")
            f.write("cd "+path_to_job+"\n")
            # f.write("export OMP_NUM_THREADS=1\n")

            if 'modules' in header.cluster:
                f.write(header.cluster['modules']+'\n')
            # f.write("module add prun/1.0\n")
            # f.write("module add intel/16.0.2.181\n")
            # f.write("module add impi/5.1.3.181\n")
            
            # if header.siman_run: #only for me
            lib64 = header.cluster['homepath'] + '/tools/lib64'
            atlas = header.cluster['homepath'] + '/tools/atlas'
            # f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"+lib64+'\n')
            # f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"+atlas+'\n')
            
            f.write("export PATH=$PATH:"+header.cluster['homepath'] +"/tools/\n")
            

            f.write("touch RUNNING\n")


    return









def prepare_input(prevcalcver = None, option = None, input_geofile = None, name_mod_prev = '', write = True, 
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
                f.write("cp "+input_geofile+" POSCAR\n")




    return





def run_command(cl, option, name, parrallel_run_command, condition = False, 
    write = True, f = None ):
    """
    Write commands used for running calculator. 

    INPUT:

        - cl (Calculation) - 
        - condition (bool) - True - allows override additional conditions (obsolete)

    """ 

    if write:

        if option == 'master':
            f.write("vasp >"+name+".log\n")

        elif 'monte' in cl.calc_method:
            f.write("python "+header.cluster_home+'/'+ header.cluster_tools+'/siman/monte.py > monte.log\n')

        elif 'polaron' in cl.calc_method:
            f.write("python "+header.cluster_home+'/'+ header.cluster_tools+'/siman/polaron.py > polaron.log\n')

        elif 'atat' in  cl.calc_method:
            f.write('maps -d&\npollmach runstruct_vasp mpirun\n')

        else:
            if cl.calculator == 'vasp':
                f.write(parrallel_run_command +" >"+name+".log\n")
            elif cl.calculator == 'gaussian':
                f.write(parrallel_run_command +" < input.gau > "+name+".out\n")
            else:
                printlog('Error! Calculator ', cl.calculator, 'is unknown!')

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
    contcar = pre+'.CONTCAR'

    if write:

        if "o" in savefile:

            f.write("mv OUTCAR "          + v + name_mod +  ".OUTCAR\n")
            f.write("mv CONTCAR "         + contcar+'\n')

        if "i" in savefile:
            f.write("cp INCAR "           + v + name_mod +  ".INCAR\n")
        
        if "v" in savefile: # v means visualization chgcar
            chg  = pre + '.CHG'
            f.write("mv CHG "+chg+"\n")
            f.write("gzip -f "+chg+"\n")
        
        if 'c' in savefile: # 
            fln = 'CHGCAR'
            chgcar  = pre +'.'+fln
            f.write('cp '+fln+' '+chgcar+'\n') #use cp, cause it may be needed for other calcs in run
            f.write('gzip -f '+chgcar+'\n')                

        if 'p' in savefile: # 
            fln = 'PARCHG'
            parchg  = pre +'.'+fln
            f.write('cp '+fln+' '+parchg+'\n') #use cp, cause it may be needed for other calcs in run
            f.write('gzip -f '+parchg+'\n') 

        if 'l' in savefile: # 
            fln = 'LOCPOT'
            locpot  = pre +'.'+fln
            f.write('cp '+fln+' '+locpot+'\n') #use cp, cause it may be needed for other calcs in run
            f.write('gzip -f '+locpot+'\n') 


        # else:
        #     f.write("rm CHG \n") #file can be used only for visualization


        if "d" in savefile:
            fln = 'DOSCAR'
            doscar  = pre +'.'+fln
            f.write('mv '+fln+' '+doscar+'\n')
            f.write('gzip -f '+doscar+'\n')                
        


        if "a" in savefile:
            f.write("mv AECCAR0 "     + v + name_mod + ".AECCAR0\n")
            f.write("mv AECCAR2 "     + v + name_mod + ".AECCAR2\n")
        
        if 'x' in savefile:
            f.write("mv vasprun.xml " + v + name_mod + ".vasprun.xml\n")

        if 't' in savefile:
            f.write("mv XDATCAR " + v + name_mod + ".XDATCAR\n")

        if 'z' in savefile:
            f.write("mv OSZICAR " + v + name_mod + ".OSZICAR\n")
       
       
       
        if 'w' in savefile:
            fln = 'WAVECAR'
            wavecar  = pre +'.'+fln
            # f.write("mv WAVECAR "     + v + name_mod + ".WAVECAR\n")
            f.write('mv '+fln+' '+wavecar+'\n') #
            f.write('gzip -f '+wavecar+'\n')  
            rm_chg_wav = rm_chg_wav.replace('w','')
        # else:
        #     f.write("rm WAVECAR\n")


        if 'c' in rm_chg_wav:
            f.write("rm CHGCAR   # rm_chg_wav flag\n") #file is important for continuation
        if 'w' in rm_chg_wav:
            ''
            f.write("rm WAVECAR  # rm_chg_wav flag\n") #
        if 'v' in rm_chg_wav: #chgcar for visualization
            ''
            f.write("rm CHG   # rm_chg_wav flag\n") #


    return contcar


def analysis_script(write = True, f = None ):
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
    input_geofile = None, parrallel_run_command = None):
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

        prepare_input(prevcalcver = prevcalcver, option = option,
         input_geofile = input_geofile, name_mod_prev = name_mod_last, write = write_poscar, curver = version,
         copy_poscar_flag = copy_poscar_flag)

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
           
            run_command(option = option, name = cl.name+name_mod, 
                parrallel_run_command = parrallel_run_command, write = write)

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
        

        analysis_script(write = write, f=f)
        # print cl.associated





    elif 'afm_ordering' in cl.calc_method:

        #Comment - inherit_xred option is not available here
        f.write("rm CHGCAR\n")                
        if not savefile: 
            savefile = 'o'

        for i, magmom in enumerate(cl.magnetic_orderings):

            name_mod   = '.AFM'+str(i)+set_mod

            cl.update_incar(parameter = 'MAGMOM', value = magmom, write = write, f = f)

            prepare_input(prevcalcver = prevcalcver, option = option, input_geofile = input_geofile,
                copy_poscar_flag = copy_poscar_flag)
            
            run_command(option = option, name = cl.name+name_mod, parrallel_run_command = parrallel_run_command)

            contcar_file = mv_files_according_versions(cl, savefile, version, name_mod = name_mod, f = f)
        
            cl.associated_outcars.append( v + name_mod +  ".OUTCAR"  )

        analysis_script()
    

    else: #simple run
        
        if not savefile: 
            savefile = 'vco'

        if write: 
            f.write("#Basic run:\n")  

        name_mod   = set_mod
        name_mod_last = ''

        prepare_input(prevcalcver = prevcalcver, option = option, name_mod_prev = name_mod_last,
            input_geofile = input_geofile, write = write_poscar, curver = version,
            copy_poscar_flag = copy_poscar_flag, f = f)

        run_command(cl = cl, option = option, name = cl.name+name_mod, 
            parrallel_run_command = parrallel_run_command, write = write, f = f)

        if final_analysis_flag:
            rm_chg_wav = 'w' #The wavcar is removed for the sake of harddrive space
        
        else:
            rm_chg_wav = ''


        contcar_file = mv_files_according_versions(cl, savefile, version, 
            write = write, name_mod = name_mod, rm_chg_wav = rm_chg_wav, f = f)

        cl.associated_outcars.append( version + name_mod +  ".OUTCAR"  )

    return contcar_file



def write_footer(cl, set_mod = '', run_tool_flag = True, 
    savefile = None, final_analysis_flag = True, neb_flag = None, f = None):
    """footer"""
    
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

            
            run_command(option = option, name = run_name_prefix+'.'+name_mod, 
                parrallel_run_command = parrallel_run_command, write = True)
            
            u_last = u


            for n_st in subfolders:

                f.write('cp '+n_st+'/CONTCAR '+n_st+'/POSCAR'+'               #u_ramp_loop()\n' )
                f.write('cp '+n_st+'/OUTCAR  '+n_st+'/'+name_mod+'.OUTCAR'+'  #u_ramp_loop()\n' )
                contcar = name_mod+'.CONTCAR'
                f.write('cp '+n_st+'/CONTCAR  '+n_st+'/'+contcar+'            #u_ramp_loop()\n' )

                # cl.associated_outcars.append( v + name_mod +  ".OUTCAR"  )


        return contcar

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



        f.write("\n\n#Starting NEB script \n")

        if option and 'continue' in option:
            prevout = name_mod_last+'OUTCAR '

            for n_st in subfolders:
                f.write('cp '+n_st+'/'+prevout+n_st+'/'+'prev.'+prevout+'  # inherit_option = continue\n' )
                f.write('cp '+n_st+'/CONTCAR '+n_st+'/POSCAR  # inherit_option = continue\n')
                f.write('mv '+n_st+'/CHGCAR '+n_st+'/prev.CHGCAR   # inherit_option = continue\n')

        else:
            
            if run_tool_flag:
                f.write('export PATH=$PATH:'+header.cluster_home+'/tools/vts/\n') #header.project_path_cluster

                f.write('nebmake.pl '+ start.replace('OUT','CONT') + final.replace('OUT','CONT') + nim_str +' \n')


        if nim+1 < 10: 
            nim_plus_one_str = '0'+str(nim+1)

        if run_tool_flag:
            f.write('cp '+start +  '00/OUTCAR\n')
            f.write('cp '+final +  nim_plus_one_str + '/OUTCAR\n' )


        cl.update_incar(parameter = 'IMAGES', value = nim, write  =1, f  = f)


        if 'u_ramping' in cl.calc_method:


            contcar_file = u_ramp_loop(cl, subfolders = subfolders, 
                run_name_prefix = cl.name+'.n_'+nim_str, 
                set_mod = set_mod, f = f)
      

        else:

            run_command(option = option, name = cl.name+set_mod+'.n_'+nim_str+name_mod, 
            parrallel_run_command = parrallel_run_command, write = True)
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




            run_command(option = option, name = cl.id[0]+'.'+cl.id[1]+'.100'+name_mod+'.fitted', 
                parrallel_run_command = parrallel_run_command, write = True)

            # print(final_analysis_flag)
            # sys.exit()

            contcar_file = mv_files_according_versions(cl, savefile, '100', 
                name_mod = name_mod, write = True, rm_chg_wav = rm_chg_wav, f = f)

        # sys.exit()


    #clean at the end
    if final_analysis_flag: 
        if header.final_vasp_clean:
            f.write('rm LOCPOT CHGCAR CHG PROCAR DOSCAR OSZICAR PCDAT REPORT XDATCAR vasprun.xml\n')
        f.write('rm RUNNING\n')



    return contcar_file, subfolders

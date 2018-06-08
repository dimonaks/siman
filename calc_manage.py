#Copyright Aksyonov D.A
from __future__ import division, unicode_literals, absolute_import 
from operator import itemgetter
import copy, traceback, datetime, sys, os, glob, shutil, re
from itertools import product

import numpy as np

try:
    # pmg config --add VASP_PSP_DIR $VASP_PSP_DIR MAPI_KEY $MAPI_KEY
    from pymatgen.ext.matproj import MPRester
    from pymatgen.io.vasp.inputs import Poscar
    from pymatgen.io.cif import CifParser
    pymatgen_flag = True 
except:
    print('pymatgen is not available')
    pymatgen_flag = False 





import header
from header import print_and_log, runBash, mpl, plt

from small_functions import is_list_like, makedir, list2string
from classes import Calculation, CalculationVasp, Description
from functions import (gb_energy_volume, element_name_inv 
     , get_from_server,  run_on_server, push_to_server, wrapper_cp_on_server)


from inout import write_xyz, read_xyz, write_occmatrix

from picture_functions import plot_mep, fit_and_plot, plot_conv
from analysis import calc_redox, matrix_diff, fit_a
from geo import replic, image_distance, scale_cell_uniformly, scale_cell_by_matrix, remove_atoms, create_deintercalated_structure, create_antisite_defect, create_antisite_defect2, local_surrounding, find_moving_atom


from set_functions import init_default_sets

from database import push_figure_to_archive
from analysis import find_polaron, neb_analysis


printlog = print_and_log

init_default_sets()




def log_func_exec():
    """
    save to history the executed form
    """
    ''





def write_batch_header(batch_script_filename = None,
    schedule_system = None, path_to_job = None, job_name = 'SuperJob', number_cores = 1  ):
    """
    self-explanatory)
    path_to_job (str) - absolute path to job folder 
    """
    with open(batch_script_filename,'w', newline = '') as f:


        if schedule_system == 'SGE':
            f.write("#!/bin/tcsh   \n")
            f.write("#$ -M aksenov@mpie.de\n")
            f.write("#$ -m be\n")
            f.write("#$ -S /bin/tcsh\n")
            f.write("#$ -cwd \n")
            f.write("#$ -R y \n")
            f.write("#$ -o "+path_to_job+" -j y\n\n")

            f.write("cd "+path_to_job+"\n")
            f.write("module load sge\n")
            f.write("module load vasp/parallel/5.2.12\n\n")


        if schedule_system == 'PBS':
            f.write("#!/bin/bash   \n")
            f.write("#PBS -N "+job_name+"\n")
            if header.WALLTIME_LIMIT:
                f.write("#PBS -l walltime=72:00:00 \n")
            # f.write("#PBS -l nodes=1:ppn="+str(number_cores)+"\n")
            if header.PBS_PROCS:
                f.write("#PBS -l procs="+str(number_cores)+"\n")
            else: # 1 node option 
                f.write("#PBS -l nodes=1:ppn="+str(number_cores)+"\n")
            # f.write("#PBS -l pmem=16gb\n") #memory per processor, Skoltech
            f.write("#PBS -r n\n")
            f.write("#PBS -j eo\n")
            f.write("#PBS -m bea\n")
            f.write("#PBS -M dimonaks@gmail.com\n")
            f.write("cd $PBS_O_WORKDIR\n")
            f.write("PATH=/share/apps/vasp/bin:/home/aleksenov_d/mpi/openmpi-1.6.3/installed/bin:/usr/bin:$PATH \n")
            f.write("LD_LIBRARY_PATH=/home/aleksenov_d/lib64:$LD_LIBRARY_PATH \n")
            f.write("module load Compilers/Intel/psxe_2015.6\n")
            f.write("module load MPI/intel/5.1.3.258/intel \n")
            f.write("module load QCh/VASP/5.4.1p1/psxe2015.6\n")
            f.write("module load ScriptLang/python/2.7\n\n")



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
     

        if schedule_system == 'SLURM':
            if '~' in path_to_job:
                print_and_log('Error! For slurm std err and out you need full paths')

            f.write("#!/bin/bash   \n")
            f.write("#SBATCH -J "+job_name+"\n")
            f.write("#SBATCH -t 250:00:00 \n")
            f.write("#SBATCH -N 1\n")
            f.write("#SBATCH -n "+str(number_cores)+"\n")
            f.write("#SBATCH -o "+path_to_job+"sbatch.out\n")
            f.write("#SBATCH -e "+path_to_job+"sbatch.err\n")
            f.write("#SBATCH --mem-per-cpu=7675\n")
            # f.write("#SBATCH -I other=avx\n") # AVX2 instructions for new node to improve speed by 18% 

            # f.write("#SBATCH --nodelist=node-amg03\n")
            if header.siman_run: #only for me
                if header.EXCLUDE_NODES:
                    f.write("#SBATCH --exclude=node-amg13\n")
                # f.write("#SBATCH --mail-user=d.aksenov@skoltech.ru\n")
                # f.write("#SBATCH --mail-type=END\n")
            f.write("cd "+path_to_job+"\n")
            f.write("export OMP_NUM_THREADS=1\n")

            f.write("module add prun/1.0\n")
            f.write("module add intel/16.0.2.181\n")
            f.write("module add impi/5.1.3.181\n")
            
            # if header.siman_run: #only for me
            f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/aksenov/tools/lib64:/home/aksenov/tools/atlas\n")
            f.write("export PATH=$PATH:/home/aksenov/tools\n")
            

            f.write("touch RUNNING\n")


    return

























def clean_history_file(history_list):
    seen = set()
    seen_add = seen.add
    return [x for x in history_list if not (x in seen or seen_add(x))]





def prepare_run():
    """
    INPUT:
        schedule_system - type of job scheduling system:'PBS', 'SGE', 'SLURM'
    """
    schedule_system = header.schedule_system
    with open('run','w', newline = '') as f:
    
        if schedule_system == 'SGE':
            f.write("#!/bin/tcsh\n")
            f.write("module load sge\n")
            f.write("module load vasp/parallel/5.2.12\n")
        elif schedule_system in ('PBS', 'SLURM'):
            f.write("#!/bin/bash\n")
        else:
            ''
            # print_and_log('Please provide schedule_system!')
            # raise RuntimeError

    # f.close()
    return


def complete_run(close_run = True):
    
    if close_run:

        with open('run','a', newline = '') as f:
            if header.schedule_system == "PBS":
                f.write("qstat\n")
                f.write("sleep 2\n")
            elif header.schedule_system == "SLURM":
                f.write("squeue\n")

            f.write("mv run last_run\n")


        if header.copy_to_cluster_flag:
            push_to_server('run',  header.cluster_home, header.cluster_address)
            run_on_server('chmod +x run', header.cluster_address)
            printlog('run sent')
    
    return












def create_additional(struct_des):
    """
    Automatically make objects in struct_des
    with .f and .fvac index


    """
    for key in copy.deepcopy(struct_des):
        if 'auto-created' in struct_des[key].des: continue
        new = copy.deepcopy(struct_des[key])
        new.des+=' fitted; des auto-created;'
        struct_des[key+'.f'] = new

        new = copy.deepcopy(struct_des[key])
        new.des+=' fitted and relaxed; des auto-created;'
        struct_des[key+'.fr'] = new


        new = copy.deepcopy(struct_des[key])
        new.des+=' with vacancy; des auto-created;'
        struct_des[key+'.fvac'] = new


        new = copy.deepcopy(struct_des[key])
        new.des+=' with vacancy; des auto-created;'
        struct_des[key+'.r'] = new


    return struct_des

def add_des(struct_des, it, it_folder, des = 'Lazy author has not provided description for me :( ', override = False):
    """
    Function adds description to struct_des dictionary;


    ###INPUT:
        - struct_des (dict)         - dict from project database
        * it (str)        - name of calculation
        * it_folder (str) - path and name of folder used for calculation in current project both on local and remote machines
        * des (str)       - description of calculation
        * override (bool) - allows to override existing field


    ###RETURN:
        None
    """

    if it not in struct_des or override:
        struct_des[it] = Description(it_folder, des)
        # hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
        hstring = 'add_des("{:s}", "{:s}", "{:s}")    #on {:})'.format(it, it_folder, des, datetime.date.today())

        try:
            if hstring != header.history[-1]: header.history.append( hstring  )
        except:
            header.history.append( hstring  )
        print_and_log("New structure name "+it+ " added to struct_des dict"+"\n")


    else:
        print_and_log("Attention! "+it+' already exist in struct_des, skipping; use override = True if you really need it; or first remove using manually_remove_from_struct_des()')
        # raise RuntimeError


    return




def update_des(struct_des, des_list):
    """
    Manuall adding of information to struct_des dictionary
    

    ###INPUT:
        - struct_des (dict)         - dict from project database
        - des_list (list of tuples) - list of new calculations to be added to database


    ###RETURN:
        - struct_des (dict) 
    """



    for des in des_list:
        if des[0] not in struct_des:
            add_des(struct_des, *des)

    return create_additional(struct_des)


def cif2poscar(cif_file, poscar_file):




    if pymatgen_flag:
        parser = CifParser(cif_file)
        s = parser.get_structures(primitive=0)[0]
        
        Poscar(s).write_file(poscar_file)
        printlog('File',poscar_file, 'created.')
    
    elif header.CIF2CELL: #using cif2cell for conversion

        print_and_log( runBash("cif2cell "+cif_file+"  -p vasp -o "+poscar_file)  )
        printlog('File',poscar_file, 'created.')

        #check
        if not os.path.exists(poscar_file):
            print_and_log("Error! cif2cell failed")

    else:
        printlog('Error! Support of cif files requires pymatgen or cif2cell; install it with "pip install pymatgen" or provide POSCAR or Abinit input file')




def determine_file_format(input_geo_file):

    supported_file_formats = {'abinit':'.geo',   'vasp':'POSCAR',   'cif':'.cif', 'xyz':'.xyz'} #format name:format specifier

    if 'POSCAR' in input_geo_file or 'CONTCAR' in input_geo_file:
        input_geo_format = 'vasp'
    elif '.cif' in input_geo_file:
        input_geo_format = 'cif'
    elif '.geo' in input_geo_file:
        input_geo_format = 'abinit'
    elif '.xyz' in input_geo_file:
        input_geo_format = 'xyz'
    else:
        printlog("Attention! File format of",input_geo_file ,"is unknown, should be from ", supported_file_formats, 'skipping')
        input_geo_format = None
    return input_geo_format



def get_file_by_version(geofilelist, version):

    """
    Find file with needed version from filelist according to certain rules
    """
    curv = None
    
    matched_files = []
    # print(geofilelist)

    for input_geofile in geofilelist: 
        
        input_geofile = os.path.normpath(input_geofile)


        input_geo_format = determine_file_format(input_geofile)
        # sys.exit()

        printlog('For file', input_geofile, input_geo_format, ' format was detected', imp = 'n')

        if input_geo_format in ['abinit',]: #version determined from token 
            # curv = int( runBash("grep version "+str(input_geofile) ).split()[1] )
            
            with open(input_geofile, 'r') as f:
                for line in f:
                    if 'version' in line:
                        curv = int(line.split()[1])


        elif input_geo_format == 'vasp': #from filename
            if '-' in input_geofile:
                curv = int(input_geofile.split('-')[-1] ) #!Applied only for phonopy POSCAR-n naming convention
            # try: 
            #     curv = int(input_geofile.split('-')[-1] ) #!Applied only for phonopy POSCAR-n naming convention
            # except:
            #     printlog('Error! Could not determine version of poscar file')

        elif input_geo_format == 'cif': #from filename
            printlog('I found cif file ', input_geofile)
            try:
                curv = int(os.path.basename(input_geofile).split('.')[0] )
            except:
                curv = None
                printlog('Failed to determine version, skipping')
        else:
            curv = None

        if curv == version:
            printlog('match', curv, version)
            matched_files.append(input_geofile)

    if len(matched_files) > 1:
        printlog('get_file_by_version(): Error! Several files have same versions:', matched_files)
    elif len(matched_files) == 0:
        input_geofile = None
    else:
        input_geofile = matched_files[0]


    return input_geofile





def smart_structure_read(filename = None, curver = 1, calcul = None, input_folder = None, input_geo_format = None, input_geo_file = None, ):
    """
    Wrapper for reading geometry files
    calcul (Calculation()) - object to which the path and version read

    curver (int) - version of file to be read
    input_geo_file or filename (str) - explicitly provided input file, has higher priority
    input_folder (str)   - folder with several input files, the names doesnot matter only versions
    input_geo_format (str) - explicitly provided format of input file 



    returns Structure()
    """

    search_templates =       {'abinit':'*.geo*', 'vasp':'*POSCAR*', 'cif':'*.cif'}

    if filename:
        input_geo_file = filename

    if input_geo_file:

        printlog("You provided the following geo file explicitly ", input_geo_file, 
            '; Version of file does not matter, I use *curver*=',curver, 'as a new version' )
        
    elif input_folder:

        print_and_log("I am searching for geofiles in folder "+input_folder+"\n" )

        if input_geo_format: 
            geofilelist = glob.glob(input_folder+'/'+search_templates[input_geo_format]) #Find input_geofile of specific format
        else:
            geofilelist = glob.glob(input_folder+'/*') 

        geofilelist = [file for file in geofilelist if os.path.basename(file)[0] != '.'   ]  #skip hidden files

        printlog('List of files:', geofilelist)

        input_geo_file = get_file_by_version(geofilelist, curver)
        # sys.exit()

        printlog('The result of getting by version', input_geo_file)

        if input_geo_file:
            printlog('File ', input_geo_file, 'was found')
        else:
            printlog('Error! No input file with version ', curver, 'was found in', input_folder)
    
    else:
        printlog('Neither *input_geo_file* nor *input_folder* were provided')


    input_geo_format = determine_file_format(input_geo_file)
    printlog(input_geo_format,' format is detected')


    if calcul:
        cl = calcul
    else:
        cl = CalculationVasp()


    if input_geo_format   == 'abinit':
        cl.read_geometry(input_geo_file)

    
    elif input_geo_format == 'vasp':
        cl.read_poscar(input_geo_file, version = curver)


    elif input_geo_format == 'cif':
        
        cif2poscar(input_geo_file, input_geo_file.replace('.cif', '.POSCAR'))
        input_geo_file = input_geo_file.replace('.cif', '.POSCAR')
        cl.read_poscar(input_geo_file)

    elif input_geo_format == 'xyz':
        #version = 1
        st = cl.init
        st = read_xyz(st, input_geo_file)
        cl.init = st
        cl.path["input_geo"] = input_geo_file
        cl.version = 1

    else:
        print_and_log("Error! smart_structure_read(): File format", input_geo_format, "is unknown")


    if cl.path["input_geo"] == None: 
        printlog("Error! Input file was not properly read for some reason")
    

    return cl.init








def name_mod_supercell(ortho = None, mul_matrix = None):

    if ortho:

        if len(set(ortho))==1:
            mod = '.s'+str(ortho[0])
        else:
            mod =  '.s'+list2string(ortho).replace(' ','')
    else:
        mod = '.s'+str(mul_matrix[0][0])+str(mul_matrix[1][1])+str(mul_matrix[2][2])

    return mod



def inherit_ngkpt(it_to, it_from, inputset):
    """
    inherit ngkpt from it_from to it_to

    """
    # print(it_to, it_from)
    
    struct_des = header.struct_des
    
    if it_to and it_from and it_from in struct_des:
        ks = inputset.vasp_params['KSPACING']

        for it in [it_to, it_from]: #just create ngkpt_dict_for_kspacings property
            if not hasattr(struct_des[it], 'ngkpt_dict_for_kspacings' ):
                struct_des[it].ngkpt_dict_for_kspacings = {}

        k_dict1 = struct_des[it_from].ngkpt_dict_for_kspacings
        k_dict2 = struct_des[it_to].ngkpt_dict_for_kspacings

        if ks in k_dict1:
            k_dict2[ks] = k_dict1[ks]
            printlog('inherit_ngkpt(): the k-grid from', it_from, 'was inherited to', it_to, imp = 'Y')
        else:
            printlog('no ngkpt for k-spacing', ks, 'in ngkpt_dict_for_kspacings of', it_from, 'ngkpt will determined from inputset k-spacing', imp = 'Y')
    
    return


def choose_cluster(cluster_name, cluster_home, corenum):
    """
    *cluster_name* should be in header.project_conf.CLUSTERS dict
    """
    if cluster_name in header.CLUSTERS:
        printlog('We use', cluster_name,'cluster')
        clust = header.CLUSTERS[cluster_name]



    else:
        printlog('Cluster', cluster_name, 'is not found, using default', header.DEFAULT_CLUSTER)
        clust = header.CLUSTERS[header.DEFAULT_CLUSTER]


    header.cluster_address = clust['address']
    header.CLUSTER_ADDRESS = clust['address']
    

    if cluster_home is None:
        header.cluster_home    = clust['homepath']
    else:
        header.cluster_home    = cluster_home
    

    #Determine cluster home using ssh
    # run_on_server('touch ~/.hushlogin', header.cluster_address)
    if header.copy_to_cluster_flag:
        header.cluster_home = run_on_server('pwd', header.cluster_address)
    else:
        header.cluster_home = ''

    printlog('The home folder on cluster is ', header.cluster_home)



    header.CLUSTER_PYTHONPATH    = clust['pythonpath']
    # header.SCHEDULE_SYSTEM    = clust['schedule']
    header.schedule_system    = clust['schedule']
    # header.CORENUM    = clust['corenum']
    
    if corenum:
        header.corenum    = corenum

    else:
        header.corenum    = clust['corenum']

    header.project_path_cluster = header.cluster_home +'/'+ header.PATH2PROJECT

    try:
        header.vasp_command = clust['vasp_com']
    except:
        header.vasp_command = None


    return





def add_loop(it, setlist, verlist, calc = None, varset = None, 
    up = 'up2', inherit_option = None, id_from = None, inherit_args = None, confdic = None,
    i_atom_to_remove = None,
    coord = 'direct', savefile = 'oc', show = '', comment = '', 
    input_geo_format = None, ifolder = None, input_geo_file = None, input_st = None,
    corenum = None,
    calc_method = None, u_ramping_region = None, it_folder = None, 
    mat_proj_cell = '',
    mat_proj_id = None, 
    cee_args = None,
    ise_new = None, it_suffix = None,
    scale_region = None, n_scale_images = 7,
    n_neb_images = None, occ_atom_coressp = None,ortho = None,
    mul_matrix = None,
    ngkpt = None,
    cluster = None, cluster_home = None,
    override = None,
    ssh_object = None,
    run = False, check_job  = 1, params = None
    ):
    """
    Main subroutine for creation of calculations, saving them to database and sending to server.

    Input:
        - it - arbitary name for your crystal structure 
        - setlist (list of str or str) - names of sets with vasp parameters from *varset* dictionary
        - verlist - list of versions of new calculations
        - calc, varset - database dictionaries; could be provided; if not then are taken from header

        - input_geo_format - format of files in input geo folder 
            'abinit' - the version is determined from the value inside the file
            'vasp', 'cif' -   the version is determined from the name of file; the names should be like POSCAR-1 for vasp files and 1.name.cif for cif
            'mat_proj' - take structure from materialsproject.org; use it_folder, len(verlist) = 1


        - up - string, possible values are: 'up1', 'up2', 'no_base'; if empty then test run is performed without saving and sending 
            'up1' - needed for normal creation of calculation and copy to server all files
            'no_base': only relevant for typconv
            is the same as "up1", but the base set is ommited
            'up2' - update only unfinished calculations
            'up3' - run only if id does not exist
        - coord - type of cooridnates written in POSCAR:
            'direct'
            'cart'

        - savefile - controls which files are saved during VASP run on server; check
            'ocvdawx' - outcar, chgcar, chg, dos, AECCAR,WAVECAR, xml

        - ifolder - explicit path to folder where to search for input geo file.

        - input_geo_file - explicit file name of input file

        - input_st - see in add_calculation()

        - it_folder - section folder (sfolder) used in struct_des; here needed with input_geo_format = mat_proj

        - show - only for read_results() ?.

        - comment - arbitrary comment for history.



        #inherit flags:
        inherit_args (dict) - to pass parameters to inherit_icalc; 
        confdic (dicts) - additional configuration parameters to inherit_icalc in dict, used for antisites
        - inherit_option (str):
            - 'continue'     - copy last contcar to poscar, outcar to prev.outcar and run again; on the next launch prev.outcar
                will be rewritten, please improve the code to save all previous outcars 
            - 'inherit_xred' - if verlist is provided, xred are copied from previous version to the next
            - all options available for inherit_icalc() subroutine (now only 'full' is tested)
        - *id_from* - see inherit_icalc()
        - ise_new (str) - name of new set for inherited calculation  ('uniform_scale')

        - occ_atom_coressp (dict) see inherit_icalc()
        - ortho, mul_matrix - transfered to inherit_icalc
        
        - ngkpt (list) - the list of k-points provided explicitly added to struct_des

        - corenum - number of cores used for calculation; overwrites header.corenum

        - calc_method - provides additional functionality:
            - 'u_ramping'    - realizes U ramping approach #Phys Rev B 82, 195128
            - 'afm_ordering' - 
            - 'uniform_scale' - creates uniformly scaled copies of the provided calculations
            - 'scale' - arbitrary scale according to mul_matrix
            using *scale_region*  and *n_scale_images* (see *scale_cell_uniformly()*)
            The copies are available as versions from 1 to *n_scale_images* and
            suffix .su appended to *it* name
            Copies to cluster *fit* utility that finds volume corresp. to energy minimum, creates 100.POSCAR and continues run 


        - u_ramping_region - used with 'u_ramping'=tuple(u_start, u_end, u_step)


        - cluster_home - override value of header.CLUSTERS

        cee_args - arguments for taking files from cee database; see get_structure_from_cee_database
            - cee_file (str) - name of file to be taken from cee database
            - section (str) - CEStorage, Catalysts, ets

        - run (bool) - complete the run file copy to server and run

        - params (dic) - dictionary of additional parameters, please move here numerous arguments
            - 'occmatrix' - explicit path to occmatrix file
            - 'update_set_dic' (dict) - additional parameters to override the existing set



    Comments:
        !Check To create folders and add calculations add_flag should have value 'add' 


    TODO:
    Now number of images is taken from self.set.vasp_params['IMAGES']; 
    In the case of generalization to other codes, set.nimages should be added and used

    Make class for cluster to save all information about cluster in one object, 
    like schedule_system, cluster address, corenum and so on 
    
    read structure in add_loop, to add calculation provied only structure, mat_proj_st_id pass in structure

    - occmatrix in params better to rename to occfile or something like this, 


    """



    def add_loop_prepare():

        nonlocal calc, it, it_folder, verlist, setlist, varset, calc_method, inherit_args, params

        if not params:
            params = {}


        params['show'] = show
        # if header.copy_to_cluster_flag:
        choose_cluster(cluster, cluster_home, corenum)
        
        if run:
            prepare_run()

        if header.first_run and header.copy_to_cluster_flag:
            prepare_run()
            header.first_run = False


        if not calc:
            calc = header.calc
            varset = header.varset



        it = it.strip()
        
        if it_folder: 
            it_folder = it_folder.strip()

        if not is_list_like(verlist):
            verlist = [verlist]

        if not is_list_like(setlist):
            setlist = [setlist]
        # print(setlist)
        if None in setlist:
            printlog('Error! None detected in setlist:', setlist)
        setlist = [s.strip() for s in setlist]


        if not is_list_like(calc_method):
            calc_method = [calc_method]


        if inherit_args is None:
            inherit_args = {}






        if ifolder: 
            if it not in ifolder: # just to be consistent with names
                print_and_log('Check ifolder !!! it is not in ifolder')
                raise RuntimeError


        return


    def add_loop_inherit():
        """

        inherit options:
        full_chg - including chg file
        """
        nonlocal  it, setlist, id_base, it_suffix
        struct_des = header.struct_des

        printlog('add_loop: starting add_loop_inherit ...', imp ='n')
        #inherit option
        inh_opt_ngkpt = ['full', 'full_chg', 'full_nomag', 'occ', 'r1r2r3', 'remove_imp', 'replace_atoms', 'make_vacancy', 'antisite'] #inherit also ngkpt
        inh_opt_other = ['supercell', 'r2r3'] # do not inherit ngkpt
        # if inherit_option in inh_opt_ngkpt+inh_opt_other:
        omit_inh_opt = ['inherit_xred', 'continue']
        if inherit_option and inherit_option not in omit_inh_opt:
            
            if 'id_base_st_type' in inherit_args and inherit_args['id_base_st_type'] == 'init':
                iti = it+'.init'
            else:
                iti = it

            if inherit_option == 'full':
                it_new = iti+'.if'
            
            if inherit_option == 'full_chg':
                it_new = iti+'.ifc'
            
            elif inherit_option == 'full_nomag':
                it_new = iti+'.ifn'

            elif inherit_option == 'occ':
                #please add additional vars to control for which atoms the inheritance should take place (added)
                it_new = iti+'.ifo' #full inheritence + triggering OMC from some other source        

            elif inherit_option == 'supercell':
               mod = name_mod_supercell(ortho, mul_matrix)
               it_new = iti+mod

            elif 'antisite' in inherit_option:
                suf = inherit_option.split('.')[-1]
                it_new = iti+'.'+suf
                # print (it_new)
                # sys.exit()

            elif inherit_option == 'make_vacancy':
                it_new = iti+'.vac'



            if it_suffix: # override default behaviour
                it_new = it+'.'+it_suffix
                it_suffix = None



            if 'up' in up:

                if it_folder:
                    section_folder = it_folder
                else:
                    section_folder = struct_des[it].sfolder

                for inputset in setlist:
                    for v in verlist:
                        id_base = (it,inputset,v)

                        inherit_icalc(inherit_option, it_new, v, id_base, calc, st_base = input_st, id_from = id_from, confdic = confdic,
                            it_folder = section_folder, occ_atom_coressp = occ_atom_coressp, i_atom_to_remove = i_atom_to_remove,
                            ortho = ortho, mul_matrix = mul_matrix, override =override, **inherit_args)
                
                    if inherit_option in inh_opt_ngkpt:
                        inherit_ngkpt(it_new, it, varset[inputset]) # 



            if ise_new:
                print_and_log('Inherited calculation uses set', ise_new)

                setlist = [ise_new,]

            else:
                print_and_log('Inherited calculation uses the same sets', setlist)


            it = it_new
        return



    def add_loop_scale():

        struct_des = header.struct_des
        nonlocal it, verlist, setlist, input_st,it_suffix
        u_scale_flag = False
        fitted_v100_id = None


        if calc_method and ('scale' in calc_method or 'uniform_scale' in calc_method):
            # print('sdfsdf')
            if 'uniform_scale' in calc_method:
                u_scale_flag = True

                it_new = it+'.su' #scale uniformly 
            else:
                it_new = it+'.sm' #scale according to mul_matrix
            

            if it_suffix:
                it_new = it+'.'+it_suffix
                it_suffix = None



            v = verlist[0]
            # printlog('add_loop_scale(): version is ', v)
            # sys.exit()


            # if up != 'up3':
            print_and_log('add_loop_scale(): Preparing   scale  calculation ... ', imp = 'Y')

            if len(verlist) > 1:
                print_and_log('Error! Currently   scale  is allowed only for one version')
            


            if it_new not in struct_des:
                if it_folder:
                    section_folder = it_folder
                else:
                    section_folder = struct_des[it].sfolder

                add_des(struct_des, it_new, section_folder, 'scale: scaled "images" for '+it+'.'+str(setlist)+'.'+str(v)   )




            verlist_new = []

            if ise_new and len(setlist) > 1:
                printlog('Error, ise_new and setlist > 1 detected!')


            for inputset in setlist:
                if inputset in varset:
                    inherit_ngkpt(it_new, it, varset[inputset])

                id_s = (it,inputset,v)




                if input_st:
                    st = input_st
                    pname = st.name
                    printlog('add_loop_scale():using input_st', pname)
                    input_st = None
                elif id_s in calc:
                    try:
                        st = calc[id_s].end
                        assert len(st.xcart) == st.natom
                        printlog('add_loop_scale(): end state of ', id_s, 'is used', st.name)
                    except:
                        st = calc[id_s].init
                        printlog('add_loop_scale(): init state of ', id_s, 'is used', st.name)

                    pname = str(id_s)
                else:
                    printlog('add_loop_scale(): starting to read input file')
                    st = smart_structure_read(curver = v, input_folder = struct_des[it].sfolder+'/'+it, 
                        input_geo_format = input_geo_format, input_geo_file = input_geo_file)
                    pname = st.name

                st = copy.deepcopy(st)
                st.magmom = [None] # added on 24.06.2017

                write_xyz(st, file_name = st.name+'_used_for_scaling')
                printlog('Scale_region is', scale_region, imp = 'y')
                
                if 'uniform_scale' in calc_method:

                    sts = scale_cell_uniformly(st, scale_region = scale_region, n_scale_images = n_scale_images, parent_calc_name = pname)
                
                else:
                    sts = scale_cell_by_matrix(st, scale_region = scale_region, n_scale_images = n_scale_images, parent_calc_name = pname, mul_matrix = mul_matrix)

                if ise_new:
                    inputset = ise_new
                    id_s = (it,inputset,v)

                cl_temp = CalculationVasp(varset[inputset], id_s)

                for i, s in enumerate(sts):
                    ver_new = i+ 1 # start from 1; before it was v+i
                    s.name = it_new+'.'+s.name
                    cl_temp.init = s
                    cl_temp.version = ver_new
                    cl_temp.path["input_geo"] = header.geo_folder + struct_des[it_new].sfolder + '/' + \
                                                it_new+"/"+it_new+'.auto_created_scaled_image'+'.'+str(ver_new)+'.'+'geo'

                    cl_temp.write_siman_geo(geotype = "init", 
                        description = s.des, override = True)
                    write_xyz(s)
                    verlist_new.append(ver_new)
                
                if 'uniform_scale' in calc_method:
                    #make version 100
                    cl_temp.version = 100
                    cl_temp.des = 'fitted with fit_tool.py on cluster, init is incorrect'
                    cl_temp.id = (it_new, inputset, 100)
                    cl_temp.state = '2. separately prepared'
                    blockdir = struct_des[it_new].sfolder+"/"+varset[inputset].blockfolder #calculation folder
                    # iid = cl_temp.id          
                    cl_temp.name = cl_temp.id[0]+'.'+cl_temp.id[1]+'.'+str(cl_temp.id[2])
                    cl_temp.dir = blockdir+"/"+ str(cl_temp.id[0]) +'.'+ str(cl_temp.id[1])+'/'
                    cl_temp.path["output"] = cl_temp.dir+str(cl_temp.version)+'.OUTCAR'
                    cl_temp.cluster_address      = header.cluster_address
                    cl_temp.project_path_cluster = header.project_path_cluster
                    calc[cl_temp.id] = cl_temp
                    # cl_temp.init = None
                    fitted_v100_id = cl_temp.id


                verlist = verlist_new

                print_and_log(len(sts), 'scale images have been created.', imp = 'y')
            




            it      = it_new
            if ise_new:
                setlist = [ise_new]
            # sys.exit()
        return u_scale_flag, fitted_v100_id






    def add_loop_take_from_database():

        nonlocal input_geo_format, it

        mat_proj_st_id = None
        if mat_proj_id:
            input_geo_format = 'mat_proj'

        if input_geo_format == 'mat_proj':
            print_and_log("Taking structure "+it+" from materialsproject.org ...", imp = 'Y')
            if it_folder == None:
                print_and_log('Error! Please provide local folder for new ', it, 'structure using *it_folder* argument! ', imp = 'Y')
            
            st = get_structure_from_matproj(it, it_folder, verlist[0], mat_proj_cell, mat_proj_id)
            mat_proj_st_id = st.mat_proj_st_id
            input_geo_file = st.input_geo_file
            input_geo_format = 'vasp'

        elif input_geo_format == 'cee_database':
            
            if it_folder == None:
                print_and_log('Error! Please provide local folder for new ', it, 'structure using *it_folder* argument! ', imp = 'Y')

            get_structure_from_cee_database(it, it_folder, verlist[0], **cee_args) #will transform it to vasp
            input_geo_format = 'vasp'
        return mat_proj_st_id


    def add_loop_neb():
        
        nonlocal n_neb_images
        nebsets = []
        neb_flag = calc_method and not set(['neb', 'only_neb']).isdisjoint(calc_method)
        if neb_flag: #put nimage values for set_sequence
            curset = varset[ setlist[0] ]
            if not n_neb_images:
                n_neb_images = varset[curset.vasp_params['IMAGES']]

            if not n_neb_images:
                print_and_log('Error! You did not provide number of NEB images nor in *n_neb_images* nor in your set!')
                raise RuntimeError


            if header.corenum % n_neb_images > 0:
                print_and_log('Error! Number of cores should be dividable by number of IMAGES')
                raise RuntimeError

            nebsets = [curset]
            
            if hasattr(curset, 'set_sequence') and curset.set_sequence:
                for s in curset.set_sequence:
                    nebsets.append(s)

            for s in nebsets:
                s.init_images_value = copy.deepcopy(s.vasp_params['IMAGES'])
                s.vasp_params['IMAGES'] = n_neb_images
            #     print s.vasp_params['IMAGES']
            # sys.exit()
            print_and_log('Attention, I update number of images in the set to', n_neb_images, 'for this calculation; ', imp = 'y')
        return neb_flag, nebsets


    def add_loop_neb2(neb_flag, nebsets):
        if neb_flag:

            if len(setlist) > 1:
                print_and_log('In "neb" mode only one set is allowed' )
                raise RuntimeError
            print_and_log('Preparing   neb  calculation ... ')

            #create necessary calculations without
            nimages = varset[setlist[0]].vasp_params['IMAGES']
            # verlist+=[ 3+v for v in range(nimages)  ] #list of images starts from 3 (1 and 2 are final and start)

            # write_batch_list+=[False for v in range(nimages)] #not used now

            #probably the add_calculation() should be used instead of the duplicating for code below, but then 
            #the creation of footer should be taken out
            #from add_calculation and put in add_loop() in the end.
            cl = calc[it, setlist[0], 2]
            
            for i in range(nimages):
                i+=3
                cl_i = copy.deepcopy(cl)
                cl_i.version = i
                cl_i.id = (cl.id[0], cl.id[1], cl_i.version)
                cl_i.name = str(cl_i.id[0])+'.'+str(cl_i.id[1])+'.'+str(cl_i.id[2])
                
                n = i - 2
                if n < 10:
                    n_st = '0'+str(n)
                elif n < 100:
                    n_st = str(n)


                cl_i.path["output"] = cl_i.dir + n_st + "/OUTCAR"
                print_and_log(i , cl_i.path["output"], 'overwritten in database')

                cl_i.associated_outcars = list([a.replace('2.', '', 1) for a in cl.associated_outcars])




                cl_i.state = '2. Ready to read outcar'

                if cl_i.id in calc: # for 'continue' mode the add_calculation() takes care, but here it should be
                                    # repeated. Again think about using only add_calculation and  write_batch_list
                    ''
                    # print_and_log('Please test code below this message to save prev calcs')
                    # if cl_i != calc[cl_i.id]
                    #     if hasattr(calc[cl_i.id], 'prev') and calc[cl_i.id].prev:
                    #         prevlist = calc[cl_i.id].prev
                    #     else:
                    #         prevlist = [calc[cl_i.id]]
                    #     cl_i.prev = prevlist
                    #     calc[cl_i.id] = cl_i
                else:
                    ''
                calc[cl_i.id] = cl_i

            #return back images values
            for s in nebsets:
                s.vasp_params['IMAGES'] = s.init_images_value #return back


    def add_loop_prepare2():

        nonlocal it
        struct_des = header.struct_des

        # print(it)
        # sys.exit()
        if it_suffix:
            it = it+'.'+it_suffix

        if it not in struct_des:
            if not it_folder:
                printlog('Error! Structure',it,'is not in struct_des, Please provide *it_folder*')
            else:
                add_des(struct_des, it, it_folder, 'auto add_des '  )


        if ngkpt: # add to struct_des
            # print (setlist)
            if len(setlist) > 1:
                printlog('Error! add_loop(): *ngkpt* parameter is allowed only with one inputset')
            else:
                kspacing = varset[setlist[0]].kspacing
                if not kspacing:
                    printlog('Error! add_loop(): no kspacing. In order to use inheritence of ngkpt I should know corresponding approximate kspacing, please provide')

                printlog('add_loop(), you provided *ngkpt*, I add', ngkpt,'to description of',it,'for kspacing',kspacing, imp = 'Y')
                struct_des[it].ngkpt_dict_for_kspacings[kspacing] = ngkpt



    def add_loop_choose_input_folder():
        
        if ifolder:
            input_folder = ifolder
        else:
            input_folder = header.geo_folder+header.struct_des[it].sfolder+"/"+it

        return input_folder

    def add_loop_try_to_read():
        # check if it is possible to read as is; please improve that section
        #not this functional is in add_calculation, probably move here
        if id in header.calc:
            clc = header.calc[id]
            clc.res()
            if '3' in clc.state or '4' in clc.state:
                printlog(id, 'is running or already read', imp = 'y')
                # continue


    def add_loop_finalize(u_scale_flag, fitted_v100_id):
        nonlocal id_base
        if u_scale_flag:
            #modify output names for fitted version 100, since it is created manually above and 
            #by add_calculation; for u-ramping names are different
            cl = calc[it, setlist[0], 1]

            calc[fitted_v100_id].path["output"] = cl.path["output"].replace('/1.', '/100.')
            
            calc[fitted_v100_id].associated_outcars = [out.replace('1.', '100.', 1) for out in cl.associated_outcars]


            # print (fitted_v100_id, calc[fitted_v100_id].associated_outcars)
            # sys.exit()

        if ise_new and hasattr(varset[ise_new], 'k_band_structure') and varset[ise_new].k_band_structure: #copy chgcar
            
            # calc[id_base].path["charge"]





            printlog('Coping CHGCAR for band structure', imp = 'y')
            wrapper_cp_on_server(calc[id_base].path["charge"], header.project_path_cluster + '/' + calc[id].dir + '/')




        if inherit_option  == 'full_chg':


            # cl.path["charge"] = cl.path["output"].replace('OUTCAR', 'CHGCAR')

            # print(calc[id_base].path)

            printlog('Coping CHGCAR for band structure', imp = 'y')

            wrapper_cp_on_server(calc[id_base].path["charge"], header.project_path_cluster + '/' + calc[id].dir + '/')







        hstring = "res_loop('{:s}', {:s}, {:s}, show = 'fo'  )     # {:s}, on {:s}  ".format(
            it, str(setlist), str(verlist), comment, str(datetime.date.today() )  )

        if hstring != header.history[-1]: 
            header.history.append( hstring  )



        if up not in ('up1','up2','up3'): 
            print_and_log("Warning! You are in the test mode, to add please change up to up1; "); 
            sys.exit()
        
        if run: #
            complete_run() # for IPython notebook
            printlog(run_on_server('./run', header.CLUSTER_ADDRESS), imp= 'Y' )
            printlog('To read results use ', hstring, '; possible options for show: fit, fo, fop, en, mag, magp, smag, maga, occ, occ1, mep, mepp', imp = 'Y')

        return u_scale_flag

    id_base = None

    # id1 = (it, setlist[0], verlist[0])

    add_loop_prepare()

    # print(verlist)

    mat_proj_st_id = add_loop_take_from_database()
    
    neb_flag, nebsets     = add_loop_neb()


    # if 
    u_scale_flag, fitted_v100_id = add_loop_scale()
    
    add_loop_inherit()
    
    add_loop_prepare2()

    """Main Loop by setlist and verlist"""
    output_files_names = []
    input_folder = add_loop_choose_input_folder()


    for inputset in setlist:

        prevcalcver = None # version of previous calculation in verlist

        for v in verlist:
            id = (it,inputset,v)
            
           
            blockdir = header.struct_des[it].sfolder+"/"+varset[inputset].blockfolder #calculation folder

            add_calculation(it,inputset,v, verlist[0], verlist[-1], 
                input_folder, blockdir, calc, varset, up, 
                inherit_option, prevcalcver, coord, savefile, input_geo_format, input_geo_file, 
                calc_method = calc_method, 
                u_ramping_region = u_ramping_region,
                mat_proj_st_id = mat_proj_st_id,
                output_files_names = output_files_names,
                run = run, input_st = input_st, check_job = check_job, params = params)
            
            prevcalcver = v


 
    add_loop_neb2(neb_flag, nebsets)

    add_loop_finalize(u_scale_flag, fitted_v100_id)
    

    return it












def add_calculation(structure_name, inputset, version, first_version, last_version, input_folder, blockdir, 
    calc, varset, up = "no",
    inherit_option = None, prevcalcver = None, coord = 'direct', savefile = None, input_geo_format = 'abinit', 
    input_geo_file = None, calc_method = None, u_ramping_region = None,
    mat_proj_st_id = None, output_files_names = None, run = None, input_st = None, check_job = 1, params = None):
    """

    schedule_system - type of job scheduling system:'PBS', 'SGE', 'SLURM'

    prevcalcver - version of previous calculation in verlist

    output_files_names - the list is updated on every call

    if inherit_option == 'continue' the previous completed calculation is saved in cl.prev list

    input_st (Structure) - Structure object can be provided instead of input_folder and input_geo_file, has highest priority


    TODO:
    make init of seqset inside calculate_nbands(), actualize_set, check_kpoints 

    """
    struct_des = header.struct_des


    id = (structure_name, inputset, version)
    id_first = (structure_name, inputset, first_version)


    cl_prev = None
    if not params:
        params = {}

    if 'show' not in params:
        params['show'] = ''

    if 'update_set_dic' not in params:
        params['update_set_dic'] = {}



    if id in calc: 
        cl = calc[id]
        status = "exist"
        printlog('add_calculation():',str(calc[id].name), " has been already created and has state: ", str(calc[id].state),)# imp = 'y')

        # print(cl.state)

        if check_job:
            if '2' in cl.state or '5' in cl.state:
                status = "ready"
                if up != 'up2':
                    cl.res(check_job = check_job, show = params['show']) 
                    return

            if "3" in cl.state: #attention, should be protected from running the same calculation once again
                status = "running"
                cl.res(check_job = check_job, show = params['show'])
                # print(check_job)
                # sys.exit()
                if '3' in cl.state and check_job: 
                    return

            elif "4" in cl.state: 
                status = "compl"
                if up == 'up2':
                    cl.init.select = None
                cl.res(check_job = check_job, show = params['show']) 
                # sys.exit()

                if up != 'up2':
                    return


    else:
        status = "new"
        print_and_log( "There is no calculation with id "+ str(id)+". I create new with set "+str(inputset)+"\n" )        




    if "up" in up:

        header.close_run = True


        if status in ["exist","compl"]: 
            print_and_log("You asked to update existing calculation with id "+ str(id)+"; results are overwritten" )         

        if status == 'compl' and inherit_option == 'continue':
            print_and_log(id, 'is completed, I will make its copy in self.prev[]', imp = 'Y' )         

            cl_prev = copy.deepcopy(calc[id])

        calc[id] = CalculationVasp( varset[id[1]] )
        
        cl = calc[id]

        cl.id = id 
        cl.name = str(id[0])+'.'+str(id[1])+'.'+str(id[2])
        cl.dir = blockdir+'/'+ str(id[0]) +'.'+ str(id[1])+'/'
        

        batch_script_filename = cl.dir +  cl.id[0]+'.'+cl.id[1]+'.run'        


        # all additional properties:

        if hasattr(cl.set, 'savefile'):
            for s in cl.set.savefile:
                if s not in savefile:
                    savefile+=s



        cl.calc_method = calc_method


        if hasattr(cl.set, 'u_ramping_nstep') and cl.set.u_ramping_nstep:
            print_and_log("Attention! U ramping method is detected from set\n\n")
            cl.calc_method.append('u_ramping')

        if hasattr(cl.set, 'afm_ordering'):
            print_and_log("Attention! afm_ordering method is detected from set\n\n")
            cl.calc_method.append('afm_ordering')



        #pass using object
        # if header.copy_to_cluster_flag:
        cl.cluster_address      = header.cluster_address
        cl.project_path_cluster = header.project_path_cluster
        cl.cluster_home = header.cluster_home
        cl.corenum = header.corenum 
        cl.schedule_system = header.schedule_system


        if mat_proj_st_id:
            cl.mat_proj_st_id = mat_proj_st_id




        if inherit_option == 'continue' and cl_prev:
            if hasattr(cl_prev, 'prev') and cl_prev.prev:
                cl.prev.extend(cl_prev.prev) #if 'continue' flag is used several times, we do not need cl.prev.prev.prev.... but cl.prev = [cl1, cl2 ..]
            else:
                cl.prev.append(cl_prev)


        if up in ['up1', 'up2', 'up3']:
            if not os.path.exists(cl.dir):
                os.makedirs(cl.dir)
                if header.copy_to_cluster_flag:
                    run_on_server("mkdir -p "+cl.dir, addr = cl.cluster_address)

            if id[2] == first_version:
                write_batch_header(batch_script_filename = batch_script_filename,
                    schedule_system = cl.schedule_system, 
                    path_to_job = header.project_path_cluster+'/'+cl.dir, 
                    job_name = cl.id[0]+"."+cl.id[1], number_cores = cl.corenum)



        if input_st:
            cl.init  = input_st
            # sys.exit()
        else:
            cl.init = smart_structure_read(curver = cl.id[2], calcul = cl, input_folder = input_folder, 
                input_geo_format = input_geo_format, input_geo_file = input_geo_file)


        if cl.path["input_geo"]:

            calc_geofile_path = os.path.join(cl.dir, os.path.basename(cl.path["input_geo"]) )
            
            # sys.exit()
            if cl.path["input_geo"] != calc_geofile_path: # copy initial geo file and other files to calc folder

                makedir(calc_geofile_path)

                shutil.copyfile(cl.path["input_geo"] , calc_geofile_path)

                #copy OCCMATRIX file as well             
                dir_1 = os.path.dirname(cl.path["input_geo"] )
                dir_2 = cl.dir
                # sys.exit()
        if 'OCCEXT' in cl.set.vasp_params and cl.set.vasp_params['OCCEXT'] == 1: #copy occfile
            if 'occmatrix' in params:
                shutil.copyfile(params['occmatrix'], cl.dir+'/OCCMATRIX' ) # file is provided explicitly

            else:
                try:
                    shutil.copyfile(dir_1+'/OCCMATRIX', cl.dir+'/OCCMATRIX' )
                except:
                    printlog('Attention! no OCCMATRIX file was found!!!')


        # if cl.des:
        cl.des = ' '+struct_des[id[0]].des + '; ' + varset[id[1]].des


        setseq = [cl.set]                                                                                                    
        if hasattr(cl.set, 'set_sequence') and cl.set.set_sequence:
            for s in cl.set.set_sequence:
                setseq.append(s)
    

        cl.check_kpoints()    
        

        for curset in setseq:
            if len(setseq) > 1:
                printlog('sequence set mode: set', curset.ise,':', end = '\n')
            curset.load(params['update_set_dic'], inplace = True)
            cl.actualize_set(curset)




        if up in ['up1', 'up2', 'up3']:
            
            cl.write_structure(str(id[2])+".POSCAR", coord, inherit_option, prevcalcver)
            
        
            out_name = cl.write_sge_script(str(version)+".POSCAR", version, 
                inherit_option, prevcalcver, savefile, 
                schedule_system = cl.schedule_system, mode = 'body',
                batch_script_filename = batch_script_filename)
            
                        
            if out_name:
                cl.path["output"] = cl.dir+out_name
            else:
                name_mod = ''
                cl.path["output"] = cl.dir+str(version)+name_mod+".OUTCAR" #set path to output
            
            #paths to other files
            cl.path["charge"] = cl.path["output"].replace('OUTCAR', 'CHGCAR')



            output_files_names.append( cl.path["output"] )




            if id == id_first:
                path_to_potcar = cl.add_potcar()
            
            for curset in setseq: #for each set
                cl.calculate_nbands(curset, calc[id_first].path['potcar'])



            if id[2] == last_version:
                list_to_copy = []
                
                cl.write_sge_script(mode = 'footer', schedule_system = cl.schedule_system, option = inherit_option, 
                    output_files_names = output_files_names, batch_script_filename = batch_script_filename, savefile = savefile )
                
                list_to_copy.extend( cl.make_incar() )
                
                list_to_copy.extend( cl.make_kpoints_file() )

                list_to_copy.append(batch_script_filename)
                
                if header.copy_to_cluster_flag: 
                    cl.copy_to_cluster(list_to_copy, up)

                    batch_on_server = cl.project_path_cluster+'/'+batch_script_filename 

                    printlog('Setting executable rights for batch script on server', batch_on_server)
                    run_on_server('chmod +x '+batch_on_server, header.cluster_address)

                    if header.siman_run or run: #for IPython should be only for run = 1 
                        cl.make_run(cl.schedule_system, batch_on_server)

            cl.state = "2. Ready for start"



        if status == "compl": 
            cl.state = '2. Can be completed but was reinitialized' #new behavior 30.08.2016


        print_and_log("\nCalculation "+str(id)+" successfully created\n\n", imp = 'Y')




    return













def inherit_icalc(inherit_type, it_new, ver_new, id_base, calc = None, st_base = None,
    id_from = None, confdic = None,
    atom_new = None, atom_to_replace = None,  
    id_base_st_type = 'end', 
    atoms_to_remove = None, del_pos = None,
    i_atom_to_remove = None, 
    id_from_st_type = 'end',
    atom_to_shift = None, shift_vector = None,
    it_folder = None, occ_atom_coressp = None, ortho = None, mul_matrix = None, override = None, use_init = None,
    ):
    """
    Function for creating new geo files in geo folder based on different types of inheritance
    Input args: 
        it_new, ver_new - name of new structure,
        id_base - new structure will be based on the final structure of this calculation;     (can be either Calculation() object or path to geo file)
        id_from - can be additionally used to adopt for example rprimd from id_from to it_new; (can be either Calculation() object or path to geo file)

        st_base - if not None then used instead of id_base; not implemented yet
        confdic (dict) - to pass more parameters


        
        inherit_type = '':
            full          - full inheritance of final state
            full_chg      - full + chg file, works only if chg file is on the same cluster
            'full_nomag'  - full except magmom which are set to None
            r2r3          - use r2 and r3 from id_from
            r1r2r3        - use r1, r2 and r3 from id_from
            remove_atoms  - removes atoms specified with *atoms_to_remove* (list of element names or list of atom numbers)
                del_pos (int) - choose specific position if several positions exist for the same ion

            replace_atoms - atoms of type 'atom_to_replace' in 'id_base' will be replaced by 'atom_new' type.
            make_vacancy  - produce vacancy by removing 'i_atom_to_remove' starting from 0
            occ           - take occ from *id_from* and create file OCCMATRIX for 
                            OMC [https://github.com/WatsonGroupTCD/Occupation-matrix-control-in-VASP]
                            - occ_atom_coressp (dict) {iatom_calc_from:iatom_calc_base, ... } (atomno starting from 0!!!)
            supercell - create orthogonal supercel using *ortho* list [a,b,c] or *mul_matrix* (3x3) ( higher priority)
            antisite  - create anitsite defect:
                        curent implimintation takes the first alkali cation and the closest to it transition metal and swap them
                confdic
                    - st_from
                    - cation
                    - trans
                    - mode 



        id_base_st_type - use init or end structure of id_base calculation.
        id_from_st_type  - init or end for id_from

        atom_to_shift - number of atom to be shifted; starting from 1.
        shift_vector - vector in decart cooridinates (Angstrom!!) by which the atom will be shifted

        - it_folder - section folder
        
        - use_init (bool) use init structure if end is empty

    Result: 
        new geo file in the input geo folder

    Output:
        no

    Depends from:
        header.struct_des
        header.calc


    Comments: 
        changes len_units of new to Angstrom!!!
        !nznucl is not calculated, since only geo is created here!
        make use of new methods for atom manipulation
        add to des which type of st is used: 'end', 'init'

    """


    # hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
    hstring = "inherit_icalc(it_new = '{:s}', ver_new = {:s}, id_base = {:s}, id_from = {:s})   # on {:s}".format(
        it_new, str(ver_new), str(id_base), str(id_from), str( datetime.date.today())   )
    if hstring != header.history[-1]: 
        header.history.append( hstring  )

    #if inherit_type not in header.history[-1] or \
    #it_new not in header.history[-1]:   header.history.append( hstring  )
    calc = header.calc
    struct_des = header.struct_des
    # override  = False
    if type(id_base) == str: #use function below for other formats
        printlog('Reading id_base from file', id_base)
        cl_base = CalculationVasp()
        cl_base.read_geometry(id_base)
        cl_base.id = ('from_file', 'from_file', cl_base.version)
        cl_base.name = id_base
        cl_base.end = cl_base.init

    else:
        # print(id_base)
            # id_base = (bytes(id_base[0]),bytes(id_base[1]),id_base[2])

        if id_base in calc and hasattr(calc[id_base], 'id'):
            ''
            valid_calc = False
            cl_base = calc[id_base]
            printlog('Taking id_base from calc:', id_base)
        

        else:

            cl_temp = CalculationVasp()
            input_folder = struct_des[id_base[0]].sfolder+'/'+id_base[0]
            printlog('Searching for id_base', id_base, 'in ',input_folder)

            cl_temp.init = smart_structure_read( curver = id_base[2], input_folder = input_folder)
            cl_temp.end = copy.deepcopy(cl_temp.init)
            cl_temp.name = id_base[0]+'from_file'
            cl_temp.id = ('temp','temp',id_base[2])
            cl_base = cl_temp





    if id_from:
        if type(id_from) == str: # if string - treated like file name
            print_and_log("I detect *id_from* path provided; taking some information from:", id_from)

            calc_from = CalculationVasp();
            calc_from.read_geometry(id_from)
            calc_from.end = calc_from.init
            calc_from_name = id_from
        else:
            print_and_log("I detect *id_from* Calculation(); taking some information from:", id_from, id_from_st_type)
            

            calc_from = calc[id_from]
            calc_from_name = calc_from.name

        if id_from_st_type == 'end':
            st_from = calc_from.end
        elif id_from_st_type == 'init':
            st_from = calc_from.init



        if cl_base.len_units != calc_from.len_units:
            print_and_log("Calculations have different len_units"); raise RuntimeError

        if it_new == id_from[0] and ver_new == id_from[2]:
            print_and_log("Warning! check your versions, you are trying to overwrite existing from structures, nothing done")
            raise RuntimeError 


    if it_new == cl_base.id[0] and ver_new == cl_base.id[2]:
        print_and_log("Warning! check your versions, you are trying to overwrite existing base structures, nothing done")
        raise RuntimeError
  



    new = copy.deepcopy(  cl_base  )

    new.len_units = 'Angstrom' #! Because from VASP

    # print(id_base_st_type)
    # sys.exit()

    new.version = ver_new

    if id_base_st_type == 'init':
        st = new.init
    elif id_base_st_type == 'end':
        st = new.end
        if not hasattr(st, 'znucl'):
            if use_init:
                st = new.init
                printlog('Attention! *use_init* flag detected, init is used instead of end')

            else:

                printlog('Error! end structure of', new.id, 'is empty! Use either init or finish calculation, check *use_init* flag!')

    # print(st.typat)


    #path to new calc
    if it_folder:
        # add_des(struct_des, it_new, it_folder, des = 'auto by inherit_icalc '+inherit_type) see below
        section_folder = it_folder

    else:
        if it_new not in struct_des:
            printlog('Error! please provide *it_folder*')
        section_folder = struct_des[it_new].sfolder


    it_new_folder = header.geo_folder + section_folder + '/' + it_new
    new.path["input_geo"] = it_new_folder + '/' +it_new+'.inherit.'+inherit_type+'.'+str(ver_new)+'.'+'geo'

    makedir(new.path["input_geo"])
    print_and_log('Path for inherited calc =', it_new_folder)




    if inherit_type == "r2r3":
        des = ' Partly inherited from the final state of '+cl_base.name+'; r2 and r3 from '+calc_from_name
        st.rprimd[1] = st_from.rprimd[1].copy()
        st.rprimd[2] = st_from.rprimd[2].copy()       
        st.update_xcart() #calculate new xcart from xred, because rprimd was changed


    elif inherit_type == "r1r2r3":
        des = ' Partly inherited from the final state of '+cl_base.name+'; r1, r2, r3 from '+calc_from_name
        st.rprimd = copy.deepcopy( st_from.rprimd )
        try:
            new.hex_a = calc_from.hex_a
            new.hex_c = calc_from.hex_c
        except:
            printlog('Attention! hex_a and hex_c were not found')

        st.update_xcart() #calculate new xcart from xred, because rprimd was changed


    elif inherit_type in ["full", ]:
        # print_and_log("Warning! final xred and xcart was used from OUTCAR and have low precision. Please use CONTCAR file \n");
        des = 'Fully inherited from the final state of '+cl_base.name


    elif inherit_type  == 'full_chg':

        des = 'Fully inherited (including chg file ) from the final state of '+cl_base.name

        # The file is copied in add_loop_finalize !!!

    elif inherit_type == "full_nomag":
        # print_and_log("Warning! final xred and xcart was used from OUTCAR and have low precision. Please use CONTCAR file \n");
        des = 'Fully inherited from the final state of '+cl_base.name+'; "magmom" set to [None]'
        st.magmom = [None]

    elif inherit_type == "occ":
        des = 'Fully inherited from the final state of '+cl_base.name+'; occupation matrix is taken from '+calc_from_name

        print_and_log('Inherit option: "occ", reading occupation matrices from',calc_from_name)
        
        if not calc_from.occ_matrices:
            print_and_log('Error! calc_from.occ_matrices is empty')
            raise RuntimeError


        #additional control to which atoms should be applied
        #if cells are different
        print_and_log('You can use *occ_atom_coressp* to control for which atoms you inherit occupations')



        if calc_from.end.natom != new.natom or occ_atom_coressp:
            
            if calc_from.end.natom != new.natom:
                print_and_log('Attention! Numbers of atoms are different. Please use *occ_atom_coressp* or I will try to use most closest to alkali ion d-atoms ')

            if not occ_atom_coressp:
                # raise RuntimeError
                print_and_log('Please run res_loop(show = "occ") for *id_base*=',id_base, ' and *id_from*=',id_from, 'to save self.dist_numb')


                occ_atom_coressp = {}
                occ_atom_coressp[calc_from.dist_numb[0][1] ] = new.dist_numb[0][1] 
            
            
            print_and_log('The occ matrices will be inherited from atom # in id_from to atom # in id_base:', end = '\n')
            occs = {}
                
            for iat_from in occ_atom_coressp: # key is atom number in id_from 
                iat_new = occ_atom_coressp[iat_from] #new is based on id_base
                occs[ iat_new ] = calc_from.occ_matrices[iat_from]
                print_and_log('        occ:',iat_from+1,'-->', iat_new+1)
        else:
            print_and_log('The cells seems to be consistent; full inheritence of occ_matrices')
            occs = calc_from.occ_matrices
        
        # sys.exit()

        write_occmatrix(occs, it_new_folder)


        # st.magmom = [None]
        



        # sys.exit()

    elif inherit_type == 'supercell':
        from geo import ortho_vec, create_supercell
        print_and_log('       inherit_icalc(): starting supercell mode ...', imp = 'Y')
        # print_and_log('rprimd is \n', st.rprimd)
        
        if mul_matrix is None: #if *mul_matrix* is not provided, try to use *ortho*
            mul_matrix = ortho_vec(st.rprimd, ortho_sizes = ortho)
            printlog('*mul_matrix* was calculated from *ortho*')
        else:
            printlog('*mul_matrix* was explicitly provided')

        print_and_log('Mul matrix is\n',mul_matrix)
        st = create_supercell(st, mul_matrix)
        # sc.mul_matrix = mul_matrix.copy()
        # new.init = sc
        # new.end  = sc
        des = 'obtained from '+cl_base.name+' by creating supercell '+str(ortho)
        override = True
    

    elif inherit_type == "atom_shift":
        des = 'obtainded from final state of '+cl_base.name+' by shifting atom '+ str(atom_to_shift) +' by '+ str(shift_vector)
        
        st.xcart[atom_to_shift-1] += np.asarray(shift_vector) 
        st.xred = xcart2xred(new.end.xcart, new.end.rprimd)



    elif inherit_type == "remove_atoms":
        """
        remove atoms either of types provided in *atoms_to_remove* or having numbers provided in *atoms_to_remove*
        """

        des = 'All atoms of type ' + str(atoms_to_remove)+' removed from the final state of '+cl_base.name


        if del_pos:
            if len(atoms_to_remove) > 1 and not is_string_like(atoms_to_remove[0]):
                printlog('Error! When *del_pos* is given,  *atoms_to_remove* should be list with one element, but it is', atoms_to_remove)

            st = create_deintercalated_structure(st, atoms_to_remove[0], del_pos = del_pos)
        else:        
            st = remove_atoms(st, atoms_to_remove)




        new.init = st
        new.end  = st
        st.name = it_new+'_from_'+new.name
        override = True









    elif inherit_type == "make_vacancy":
        """Remove  atom 'i_atom_to_remove' from final state of id_base"""


        print_and_log('Warning! Please check inherit_type == "make_vacancy", typat can be wrong  if more than one element present in the system\n ',
            'Use del_atom() method ')
        # raise RuntimeError

        des = 'Atom '+str(i_atom_to_remove)+' removed from  '+cl_base.name

        del st.typat[i_atom_to_remove]
        del st.xcart[i_atom_to_remove]
        del st.xred[i_atom_to_remove]
        ntypat = len(set(st.typat))
        znucl_new = []
        for t in sorted(set(st.typat)):
            znucl_new.append( st.znucl[t-1]  )
        st.znucl = znucl_new
        st.natom -=1
        # new.write_geometry("end", des)      

        #make visualization by adding to vacancy hydrogen atom
        st_b_copy = copy.deepcopy(cl_base.end)
        st_b_copy.typat[i_atom_to_remove] = max(st_b_copy.typat)+1
        st_b_copy.ntypat+=1
        st_b_copy.znucl.append(1)
        write_xyz(st_b_copy, file_name = "test_of_vacancy_creation."+str(new.version)+"."+st_b_copy.name)






    elif inherit_type == "replace_atoms":
        """Simply replace in calculation one atoms by another """
        z_new     = element_name_inv(atom_new)
        z_replace = element_name_inv(atom_to_replace)




        znucl = st.znucl
        
        if z_replace not in znucl: 
            print_and_log("Error! Calc "+new.name+" does not have atoms of this type")
            raise RuntimeError

        if atom_to_replace not in id_base[0] or atom_new not in it_new:
            print_and_log("Error! inherit_icalc(): Something wrong with names of atom types")
            raise RuntimeError            
        print_and_log( "replace ", z_replace, "by", z_new)

        znucl = [int(z) for z in znucl] # convert to int
        i_r = znucl.index(z_replace)
        print_and_log( "index ", i_r)
        znucl[i_r] = z_new

        if znucl[-2] == znucl[-1]: #just for special case
            del znucl[-1]

            for i, t in enumerate(st.typat):
                if t == st.ntypat:
                    print_and_log( "found ntypat" ,t)
                    st.typat[i]-=1
                    #print t
            st.ntypat-=1
            #print st.typat
            #print new.end.typat
        
        st.znucl = znucl



        des = 'Fully inherited from the '+ id_base_st_type +' state of '+cl_base.name+\
        ' by simple replacing of '+atom_to_replace+' by '+atom_new

        override = 1



    elif 'antisite' in inherit_type:

        if 'as' in inherit_type:
            st = create_antisite_defect(st)
            des = 'Fully inherited from the '+ id_base_st_type +' state of '+cl_base.name+\
            ' by simple swapping of alkali and transition atoms'

        elif 'a' in inherit_type:
            st = create_antisite_defect2(st, st_from = confdic['st_from'], cation = confdic['cation'], trans = confdic['trans'], mode = confdic['mode'])
            des = 'create_antisite_defect2 '

        override = True






    else:
        print_and_log("Error! Unknown type of Calculation inheritance")



    




    #auto addition of description
    if it_new not in struct_des: 
        add_des(struct_des, it = it_new, it_folder = it_folder, des = 'auto '+des)
        new.des =  struct_des[it_new].des
    else:
        new.des = des + struct_des[it_new].des
        if it_folder:
            struct_des[it_new].sfolder = it_folder #update itfolder,




    if mul_matrix is not None:
        struct_des[it_new].mul_matrix = mul_matrix
    #write files

    # print new.end.xcart

    #print(len(new.end.xred))
    # print (id_base_st_type)
    new.init = st
    new.write_geometry('init', des, override = override)
    

    write_xyz(st)





    return




























def res_loop(it, setlist, verlist,  calc = None, varset = None, analys_type = 'no', b_id = None, 
    typconv='', up = "", imp1 = None, imp2 = None, matr = None, voronoi = False, r_id = None, readfiles = True, plot = True, show = 'fo', 
    comment = None, input_geo_format = None, savefile = None, energy_ref = 0, ifolder = None, bulk_mul = 1, inherit_option = None,
    calc_method = None, u_ramping_region = None, input_geo_file = None, corenum = None, run = None, input_st= None,
    ortho = None, mat_proj_cell = None,
    it_folder = None, choose_outcar = None, choose_image = None, 
    cee_args = None, mat_proj_id = None, ise_new = None, push2archive = False,
    description_for_archive = None, old_behaviour  = False,
    alkali_ion_number = None, cluster = None, ret = None, override = None, check_job = 1, fitplot_args = None, style_dic = None, params = None):
    """Read results
    INPUT:
        'analys_type' - ('gbe' - calculate gb energy and volume and plot it. b_id should be appropriete cell with 
            bulk material,
            'e_imp' ('e_imp_kp', 'e_imp_ecut') - calculate impurity energy - just difference between cells with impurity and without.
            'fit_ac' - fit a and c lattice constants using 2-dimensianal spline
            'clusters' - allows to calculate formation energies of clusters
            'diff' - difference of energies in meV, and volumes A^3; E(id) - E(b_id)
            'matrix_diff' - difference normalized by matrix atoms

            'redox_pot' - calculate redox potential relative to b_id() (deintercalated cathode) and energy_ref ( energy per one ion atom Li, Na in metallic state or in graphite)
            'neb' - make neb path. The start and final configurations 
            should be versions 1 and 2, the intermidiate images are starting from 3 to 3+nimages
            )
        voronoi - True of False - allows to calculate voronoi volume of impurities and provide them in output. only if lammps is installed
        b_id - key of base calculation (for example bulk cell), used in several regimes; 
        r_id - key of reference calculation; defines additional calculation (for example atom in vacuum or graphite to calculate formation energies); can contain directly the energy per one atom

        up - 
            if equal to 'up2' the files are redownloaded; also can be used to download additional files can be 'xo' (deprecated?)
            - if 'un' is found in up then siman will try to read unfinished outcars


        readfiles (bool) - True - read from outcar, False - read from database; 



        The next three used for 'clusters' regime:    
        imp1 - key of bulk cell with one imp1
        imp2 - key of bulk cell with one imp2
        matr - key of bulk cell with pure matrix.


        - show - (str), allows to show additional information:
            - mag - magnetic moments on magnetic atoms
              maga
                *alkali_ion_number* (int) - number of atom around which to sort mag moments from 1
            - en  - convergence of total energy vs max force
            - mep - neb path
            - fo  - max force on each md step
            - polaron - determine polaron positon, write local surroundin
            - mig_path - write migration path xyz



        energy_ref - energy in eV; substracted from energy diffs
        
        bulk_mul - allows to scale energy and volume of bulk cell during calculation of segregation energies

        choose_outcar (int, starting from 1)- if calculation have associated outcars, you can check them as well, by default
        the last one is used during creation of calculation in write_sge_script()

        choose_image (int) - relative to NEB, allows to choose specific image for analysis, by default the middle image is used

        - push2archive (bool) - if True produced images are copied to header.project_conf.path_to_images
        - description_for_archive - caption for images

        ret (str) - return some more information in results_dic
            'energies' - just list of full energies
        check_job - (bool) check status on server, use 0 if no internet connection

        fitplot_args - additional arguments for fit_and_plot function

        style_dic - passed to plot_mep()

        - ise_new - dummy
        - inherit_option - dummy
        - savefile - dummy
        - cluster - dummy
        -override - dummy
        
        - params - dictionary of additional parameters to control internal, many arguments could me moved here
            mep_shift_vector - visualization of mep in xyz format

    RETURN:
        (results_dic,    result_list)
        or 
        (result_string, result_list)

        result_list - list of results; was used in previous versions, now left for compatibility
        results_dic - should be used in current version! actually once was used as list, now should be used as dict

    TODO:
        Make possible update of b_id and r_id with up = 'up2' flag; now only id works correctly


    """


    """Setup"""




    if not is_list_like(verlist):
        verlist = [verlist]

    if not is_list_like(setlist):
        setlist = [setlist]


    if not calc:
        calc = header.calc



    try:
        b_ver_shift = b_id[2] #add to version of base with respect to version of main
    except:
        b_ver_shift = 0



    if '2' in up:
        loadflag = up+'o'
    else:
        loadflag = up

    # if choose_outcar:


    name_field_length = 30

    header.show_head = 1 # head before the string of read_results()

    conv = {}
    base = 'base'; 
    n    = 'temp'
    conv[n] = []
    conv[base] = []
    conv[it] = []
    result_list = []
    energies = []



    emin = 0
    if b_id:
        if b_id in calc:
            if len(b_id) == 3: # for all cases besides e_seg and coseg for wich b_id is determined every iteration
                # print "Start to read ", b_id
                # if '4' not in calc[b_id].state:
                if readfiles:
                    calc[b_id].read_results(loadflag, choose_outcar = choose_outcar)
                
                e_b = 1e10; v_b = 1e10
                if '4' in calc[b_id].state:
                    e_b = calc[b_id].energy_sigma0
                    v_b = calc[b_id].end.vol
                else:
                    print_and_log('Warning! Calculation ',b_id, 'was not finished; please check, now skipping ...', important = 'y')
        else:
            printlog('Attention! res_loop(): b_id', b_id, 'does not exist. return {} []')
            # return {}, []

    #define reference values
    e1_r = 0
    if type(r_id) in (float, int):
        e1_r = r_id
    elif type(r_id) == tuple:
        # if '4' not in calc[r_id].state:
        #     print "Start to read reference:"
        if readfiles:
            print_and_log( calc[r_id].read_results(loadflag, choose_outcar = choose_outcar)  )
        e_r = calc[r_id].energy_sigma0 #reference calc
        nat_r = calc[r_id].end.natom # reference calc
        e1_r = e_r/nat_r # energy per one atom
        # print e1_r




    """Main loop"""
    final_outstring = 'no calculation found'
    for inputset in setlist:
        for v in verlist:
            id = (it,inputset,v)
            # print(id)
            if id not in calc:
                printlog('Key', id,  'not found in calc!', imp = 'Y')
                continue #pass non existing calculations
            cl = calc[id]
           
            if readfiles and check_job:
                if '3' in cl.check_job_state():
                    printlog( cl.name, 'has state:',cl.state,'; I will continue', cl.dir, imp = 'y')
                    # cl.res()
                    continue


            if 'path' in show:
                printlog(os.getcwd()+'/'+cl.path['output'], imp = 'Y')
                # sys.exit()
                return

            if not hasattr(cl,'version'):
                cl.version = v

            
            if readfiles:
                printlog('Starting self.read_results() ...')
                outst = calc[id].read_results(loadflag, analys_type, voronoi, show, 
                    choose_outcar = choose_outcar, alkali_ion_number = alkali_ion_number)
                
                if '5' in cl.state:
                    continue

            else:
                outst = ' output was not read '

            printlog('read_results() output', outst)



            if analys_type in ('e_seg', 'coseg'):
                try:
                    b_id[1] 
                    b_id = (b_id[0], b_id[1], id[2] + b_ver_shift)
                except:
                    b_id = (b_id[0], id[1], id[2] + b_ver_shift)
                printlog('b_id', b_id)


            # print(id)
            # sys.exit()
            e = calc[id].energy_sigma0
            n_m = calc[id].end.nznucl[0] # number of matrix atoms


            try:
                v = calc[id].end.vol
            except:
                v = 0
            #print e
            if e < emin: 
                emin = e; 
                id_min = id
            
            conv[n].append(id)
            # print base
            conv[base].append(b_id)


            outst2 = ("%s"%calc[id].name).ljust(name_field_length)
            outst2+='|'
            outst_end = '' 

            energies.append(cl.energy_sigma0)


            
            if   b_id :

                # if "4" not in calc[b_id].state:
                if readfiles:    
                    # print(b_id)
                    calc[b_id].read_results(loadflag, choose_outcar = choose_outcar)

                # print(b_id)
                if "4" in calc[b_id].state:    

                    if calc[id].set.ngkpt != calc[b_id].set.ngkpt:
                        print_and_log("Warning! you are trying to compare calcs with "+str(calc[id].set.ngkpt)+" and "+str(calc[b_id].set.ngkpt)+"\n")
                        pass

                    if calc[id].NKPTS != calc[b_id].NKPTS:
                        print_and_log("Warning! you are trying to compare calcs with "+str(calc[id].NKPTS)+" and "+str(calc[b_id].NKPTS)+" kpoints \n")

                    if 'gbe' in analys_type:      
                        outst2 = gb_energy_volume(calc[id], calc[b_id])
                
                    elif 'e_imp' in analys_type:
                        calc[id].e_imp = e - e_b
                        calc[id].v_imp = v - v_b

                        #calc[id].v_imp = e - e_b
                        outst2 += ("%.3f & %.2f  & %.2f &"% (e - e_b, v - v_b, (v - v_b)/v_b*100 ) )
                        conv['e_imp'].append(id)
                        a    = calc[id].hex_a;                     a_b  = calc[b_id].hex_a
                        c    = calc[id].hex_c;                     c_b  = calc[b_id].hex_c
                        ca = c/a;                                   ca_b = c_b/a_b
                        outst_end = " & {0:.1f} & {1:.1f} & {2:.1f} & {3:.3f}".format((a - a_b)/a_b*100,  (c - c_b)/c_b*100, (ca - ca_b)/ca_b*100, (e - e_b - e1_r) )

                    elif analys_type == 'e_2imp': #"""For calculation of energies of two impurities in big cell"""
                        calc[id].e_imp = e - e_b            
                        outst2 += ("%.0f "% ( (e - e_b)*1000 ) )  
                        conv[it].append(id)


                    elif analys_type in ('e_seg', 'coseg'): #"""For calculation of segregation and cosegregation energies"""
                        e_b = calc[b_id].energy_sigma0 * bulk_mul
                        n_m_b = calc[b_id].end.nznucl[0]
                        v_b = calc[b_id].end.vol * bulk_mul
                        # diffE = e - e_b/n_m_b*n_m
                        diffE = e - e_b

                        # outst2 += ("%.0f & %.2f "% ( (e - e_b)*1000, v - v_b ) )
                        
                        outst2 += " {:.3f} & {:.2f} ".format( (diffE - energy_ref), (v - v_b) ).center(6)
                        outst2 +='&'
                        # write_xyz(calc[id].end)
                        # write_xyz(calc[b_id].end)
                        result_list = [diffE - energy_ref, v - v_b]


                    elif analys_type == 'matrix_diff': #
                        print_and_log( 'Calculating matrix_diff...')
                        

                        diffE, diffV = matrix_diff(calc[id], calc[b_id], energy_ref)
                        outst2 += " {:.3f} & {:.2f} &".format( diffE, diffV ).center(6)
                        result_list = [diffE, diffV]


                    elif analys_type == 'diff': #
                        print_and_log( 'Calculating diff...')
                        e_b = calc[b_id].energy_sigma0
                        v_b = calc[b_id].end.vol
                        diffE = e - e_b
                        
                        outst2 += " {:.3f} & {:.2f} &".format( (diffE - energy_ref), (v - v_b) ).center(6)
                        result_list = [diffE - energy_ref, v - v_b]


            if analys_type == 'clusters':
                e1  = calc[imp1].energy_sigma0
                e2  = calc[imp2].energy_sigma0
                e_m = calc[matr].energy_sigma0
                n1 = calc[id].init.nznucl[1]
                if len(calc[id].init.nznucl) == 3:
                    n2 = calc[id].init.nznucl[2]
                else:
                    n2 = 0
                # print n1,n2
                outst2 += ("%.0f "% ( (e - n1*e1 - n2*e2 + (n1+n2-1)*e_m )*1000 / (n1+n2) ) )
            


            final_outstring = outst2+outst + outst_end     
            # print ([final_outstring])         
            print_and_log( final_outstring, end = '\n',  imp = 'Y')

        emin = 0
        





        """Aditional analysis, plotting"""
        results_dic = {} #if some part fill this list it will be returned instead of final_outstring
        if ret == 'energies':
            results_dic[ret] = energies

        cl = calc[id]


        if id not in calc or '4' not in calc[id].state:
            # printlog(calc[id].state, imp = 'Y')
            print_and_log( "res_loop(): Calculation ",id, 'is unfinished; return \{\} []',cl.dir, imp = 'Y')
            return {}, []
        

        if b_id: 
            bcl = calc[b_id]

        if analys_type == 'gbe':
            print_and_log("\nGrain boundary energy and excess volume fit:")
            plot_conv( conv[n], calc, "fit_gb_volume")

        elif analys_type == 'gbep':
            print_and_log("\nGrain boundary energy and excess volume fit:")
            # plot_conv( conv[n], calc, "fit_gb_volume")
            final_outstring = plot_conv( conv[n], calc, "fit_gb_volume_pressure")

        elif analys_type in ('e_seg', 'coseg') and len(verlist) > 3:

            #Test lateral sizes
            A   = calc[id].end.yzarea
            A_b = calc[b_id].end.yzarea
            
            if A != A_b: 
                print_and_log("Warning! you are trying to compare calcs with different lateral sizes: "+str(A)+" "+str(A_b))
                print_and_log( "Areas are ", A, A_b," A^3")
            
            #Show results 
            id1 = (it,inputset,verlist[0]) #choosen to save calculated values at first version of version set
            
            if readfiles and plot:           
                #print " \n\nImpurity at the interface :"
                e, v, emin, vmin       = plot_conv( conv[n], calc,  "fit_gb_volume2")
                #print " \n\nImpurity in the volume    :"
                e_b, v_b, e_bmin, v_bmin = plot_conv( conv[base], calc, "fit_gb_volume2")
                e_segmin = (emin - e_bmin) * 1000
                v_segmin =  vmin - v_bmin

                
                e_seg = (e - e_b * bulk_mul) * 1000
                v_seg =  v - v_b * bulk_mul


                # print(e_seg, e_segmin)

                calc[id1].e_seg = e_seg
                calc[id1].v_seg = v_seg
            
            if not hasattr(calc[id1], 'e_seg'): 
                print_and_log( "Warning! Calculation ", id1, 'does not have e_seg and v_seg. Try to run with readfiles = True to calculate it.')
                calc[id1].e_seg = 0; calc[id1].v_seg = 0
            


            natom = calc[id1].natom
            calc[id1].X = 1./natom
            v1 = v / natom
            calc[id1].Xgb = v1 / A # for grain boundary with 1 A width. For other boundaries should be divided by width. 
            #print ("__________________________________________________________________________")
            
            # print (" At zero pressure: segregation energy is %.0f  meV; Seg. volume is %.1f A^3; excess seg. vol. is %.2f A" %(e_seg, v_seg, v_seg/A ) )
            print ("At min: segregation energy is %.0f  meV; Seg. volume is %.1f A^3; excess seg. vol. is %.2f A" %(e_segmin, v_segmin, v_segmin/A ) )
            # print ("%s.fit.pe & %.0f & %.1f & %.2f & %.3f & %.1f" %(id[0]+'.'+id[1], e_seg, v_seg, v_seg/A, 1./A, 1./calc[id].natom * 100  ) )
            
            #Calculate distance from impurity to boundary and number of neighbours for version 2!
            id2 =(it,inputset, 2)
            st = calc[id2].end
            gbpos2 = calc[id2].gbpos 
            print_and_log( id2, 'is id2')
            # print gbpos2
            # print st.rprimd[0][0]/2.
            if gbpos2 == None:
                gbpos2 = 100
            gbpos1 = gbpos2 - st.rprimd[0][0]/2.
            d1 = abs(st.xcart[-2][0] - gbpos2)
            d2 = abs(st.xcart[-1][0] - gbpos2)
            dgb = d1; 
            iimp = -2
            if d2 < d1: 
                dgb = d2
                iimp = -1
            
            t = st.typat[iimp]
            z = st.znucl[t-1]
            segimp = element_name_inv(z) #Type of impurity closest to gb
            # print segimp, d

            id_m2   = (it+'.m',      '8'+inputset[1:], 2)
            
            if analys_type == 'e_seg':

                #calc e_seg2 and decomposition to mechanical and chemical contributions

                if id_m2 in calc: #additional analysis
                    b_id2 = (b_id[0],inputset, 2)
                    b_id_m2 = (b_id[0]+'.m', '8'+inputset[1:], 2)

                    e_seg2 = (calc[id2].energy_sigma0 - calc[b_id2].energy_sigma0) * 1000
                    e_m2   = (calc[id_m2].energy_sigma0 - calc[b_id_m2].energy_sigma0) * 1000
                    e_ch2  = e_seg2 - e_m2
                else:
                    e_seg2 = 0
                    e_m2   =0
                    e_ch2  =0


                #calculate number of close neibours around closest to gb imp
                x_central = st.xcart[iimp]
                st_r  = replic(st,   mul = (1,2,2), inv =  1 )
                st_rr = replic(st_r, mul = (1,2,2), inv = -1 ) # to be sure that impurity is surrounded by atoms

                dmax = 3
                nlist = [ x  for x, t  in zip(st_rr.xcart, st_rr.typat) if np.linalg.norm(x_central - x) < dmax and t == 1]
                nneigbours =  len(nlist)


                final_outstring = ("%s.fit.pe & %.0f & %.1f & %.2f & %.d & %4.0f & %4.0f & %4.0f & %s " %(
                    id2[0]+'.'+id2[1], calc[id1].e_seg, calc[id1].v_seg, dgb, nneigbours, e_seg2, e_ch2, e_m2, segimp  ))

                # final_outstring = ("%s.fit.pe & %.0f & %.0f & %.1f & %.1f & %.2f & %.d & %4.0f & %4.0f & %4.0f & %s " %(
                #     id2[0]+'.'+id2[1], calc[id1].e_seg, e_segmin, calc[id1].v_seg, v_segmin , dgb, nneigbours, e_seg2, e_ch2, e_m2, segimp  )) #e_segmin and v_segmin are minimum energy (but at some pressure) and corresponing volume



                results_dic = [id2[0]+'.'+id2[1], calc[id1].e_seg, calc[id1].v_seg, dgb, nneigbours, e_seg2, e_ch2, e_m2, segimp]
            
            elif analys_type == 'coseg' :
                calc[id2].e_seg = calc[id1].e_seg #save in version 2
                calc[id2].v_seg = calc[id1].v_seg
                final_outstring = ("%s.fit.pe & %.0f & %.1f & %.1f & %.1f" %(id[0]+'.'+id[1], calc[id2].e_seg, calc[id2].v_seg, d1, d2 ))




            print_and_log(  final_outstring)
            print_and_log( '\\hline')


        elif analys_type == 'e_2imp':
            # plot_conv( conv[it], calc,  analys_type, conv_ext) #instead use plot_conv( conv['hs443OO'], calc,  'e_2imp', [conv['hs443CO'],conv['hs443CC']]) 
            pass



        elif analys_type == 'fit_ac':

            print_and_log ("name %s_template          acell  %.5f  %.5f  %.5f # fit parameters are &%.5f &%.5f &%i &%i"  % (fit_hex(0.00002,0.00003,4000,6000, it, inputset, verlist, calc) )  )    

        elif 'fit_a' in analys_type:
            
            fit_a(conv, n, description_for_archive, analys_type, show, push2archive)


        # elif analys_type == 'fit_angle':




        elif analys_type == 'dimer':
            """Fit of md steps obtained from constant speed; see vasp description for dimer"""
            # print calc[id].list_e_sigma0
            # calc[id].list_dE = []
            # for E in calc[id].list_e_without_entr:

            #     calc[id].list_dE.append(E - 2*calc[b_id].e_without_entr)

            # print calc[id].list_e_without_entr
            # print calc[id].list_dE
            calc[id].e_ref = calc[b_id].e_without_entr

            plot_conv( [id], calc,  "dimer")


        elif analys_type == 'redox_pot':
            
            if '4' not in bcl.state:
                print_and_log("res_loop: Calculation ",bcl.id, 'is unfinished; return', imp = 'Y')
                return {}, []

            results_dic = calc_redox(cl, bcl, energy_ref)
            



        if analys_type == 'neb':
            results_dic = neb_analysis(cl, show, up, push2archive, old_behaviour, results_dic, fitplot_args, style_dic, params)








    if results_dic:
        return results_dic, result_list
    else:
        return final_outstring.split('&'), result_list # only for last version or fit depending on type of analysis








def create_phonopy_conf_file(st, path = '', mp = [10, 10, 10], filetype = 'mesh'):

    """filetype
            mesh - mesh.conf
            band - band.conf
    """


    mpstr = " ".join(map(str, mp))


    if filetype == 'band':
        filename = path+'/band.conf'

    else:
        filename = path+'/mesh.conf'

    with open(filename, 'w', newline = '') as f:
        f.write("DIM = 1 1 1\n")
        f.write("ATOM_NAME = ")
        for z in st.znucl:
            el = element_name_inv(z)
            f.write(el+' ')
        f.write('\n')

        # f.write("\nDIAG = .TRUE.\n")
        # f.write("DISPLACEMENT_DISTANCE = 0.03\n")    

        # f.write("TETRAHEDRON\n")
        # f.write("SIGMA = 0.1\n")
        if filetype == 'mesh':

            f.write("MP = {:}\n".format( mpstr ))
        
        if filetype == 'band':
            f.write("BAND = {:}\n".format( '0.5 0.5 0.5  0.0 0.0 0.0  0.5 0.5 0.0  0.0 0.5 0.0' ))




def read_phonopy_dat_file(filename):

    """
    read .dat file from phonopy
    for reading yaml see self.read_phonopy_data()
    should be probably combined

    freq - frequency in THz
    tot - total energy for freq
    """

    dos = {'tot':[], 'freq':[]}

    with open(filename, 'r', ) as f:
        f.readline()
        for line in f:
            val = line.split()
            dos['freq'].append(float(val[0]))
            dos['tot'].append(float(val[1]))

    return dos


def read_phonopy_data(filename, key = "free_energy", convert = False):
    """
    convert (bool) - convert kJ/mol to eV

    """
    # with open(filename, 'r') as f:
    #     f.readline()
    F = []
    T = []
    #     for line in f:
    #         T.append(float(line.split()[0]))
    #         F.append(float(line.split()[1]))
    #     # print T, F
    import yaml

    if convert:
        mul = header.kJ_mol2eV
    else:
        mul = 1
    # print(mul)

    f = open(filename)
    # use safe_load instead load
    dataMap = yaml.safe_load(f)
    f.close()
    prop = dataMap['thermal_properties']
    for i in range(len(prop)):
        T.append( prop[i]['temperature'] )
        F.append( prop[i][key]*mul     )

    coeffs1 = np.polyfit(T, F, 8)
    fit_func = np.poly1d(coeffs1)
    T_range = np.linspace(min(T), max(T))
    

    fit_and_plot(d1 = (T, F, 'b-', 'orig'), d2 = (T_range, fit_func(T_range), 'r--', 'fit'), show = 1)

    print_and_log( 'I return', key)

    return T_range, fit_func




def for_phonopy(new_id, from_id = None, calctype = 'read', mp = [10, 10, 10], additional = None):
    #creates file for phonopy, run phonopy
    #new_id - will add this calculation or read; if string then interpreted as filename of thermal_properties.yaml
    #from_id - tuple - than will take coordinates from the end; or path to input poscar file
    #type - 'create', 'read'
    #additional - list of calculation names, if calculation was splited into several parts

    mpstr = " ".join(map(str, mp))





    if calctype == 'create':
        work_path = 'geo/'+struct_des[new_id[0]].sfolder+'/'+new_id[0]
    
        log_history(  "{:}    #on {:}".format( traceback.extract_stack(None, 2)[0][3],   datetime.date.today() )  )

        # print type(from_id)
        if type(from_id) == str:
            from_cl = CalculationVasp()
            from_cl.read_poscar(from_id)
            state = 'init'
            from_st = from_cl.init
        else:
            from_cl = header.calc[from_id]
            state = 'end'
            from_st = from_cl.end

            # print  header.varset
            #create POSCAR
        posname = 'phonopy_input_poscar'
        from_cl.write_structure(name_of_output_file = posname, path = work_path, state = state)

        savedPath = os.getcwd()
        os.chdir(work_path)
        #create conf 
        confname = new_id[0]+'.conf'
        with open(confname, 'w', newline = '') as f:
            f.write("DIM = 1 1 1\n")
            f.write("ATOM_NAME = ")
            for z in from_st.znucl:
                el = element_name_inv(z)
                f.write(el+' ')
            f.write("\nDIAG = .TRUE.\n")
            f.write("DISPLACEMENT_DISTANCE = 0.03\n")

         
        #run phonopy
        print_and_log(
            runBash('export PYTHONPATH=~/installed/phonopy-1.9.5/lib/python:$PYTHONPATH; rm POSCAR-*; /home/dim/installed/phonopy-1.9.5/bin/phonopy '
            +confname+' -c '+posname+' -d --tolerance=0.01'), imp = 'y' )

        ndis = len( glob.glob('POSCAR-*') )
        print_and_log( ndis, ' displacement files was created', )

        os.chdir(savedPath)


        #add
        add_loop(new_id[0],   new_id[1], range(1,ndis+1), up = 'up1', input_geo_format = 'vasp', savefile = 'ocdx')

        #copy SPOSCAR - an ideal cell
        src =  'geo/'+struct_des[new_id[0]].sfolder+'/'+new_id[0]
        dst = struct_des[new_id[0]].sfolder+'/'+new_id[0]+'.'+new_id[1]+'/'
        shutil.copy(src+'/SPOSCAR', dst)
        shutil.copy(src+'/disp.yaml', dst)

    if calctype == 'read':

        if type(new_id) == tuple:
            new_cl = header.calc[new_id]

            work_path = struct_des[new_id[0]].sfolder+'/'+new_id[0]+'.'+new_id[1]
            work_path_geo = 'geo/'+struct_des[new_id[0]].sfolder+'/'+new_id[0]



            npos = len( glob.glob(work_path+'/*.POSCAR') )
            # print range(1,npos+1), 'range'
            if not os.path.exists(work_path+"/1.POSCAR"):
                res_loop(new_id[0],   new_id[1], range(1,npos+1), up = 'up1', input_geo_format = 'vasp', )

            if additional:
                for name in additional:
                    npos_new = len( glob.glob(struct_des[new_id[0]].sfolder+'/'+name+'.'+new_id[1]+'/*.POSCAR') )

                    if not os.path.exists(struct_des[new_id[0]].sfolder+'/'+name+'.'+new_id[1]+'/'+str(npos+1)+".POSCAR"):
                        res_loop(name,   new_id[1], range(npos+1, npos+npos_new+1), up = 'up1', input_geo_format = 'vasp', )
                    
                    npos = npos+npos_new

                    runBash("rsync "+struct_des[new_id[0]].sfolder+'/'+name+'.'+new_id[1]+'/*.vasprun.xml '+work_path)
                    # print 'Additional vasprun.xml files were copied to ', work_path

            savedPath = os.getcwd()
            os.chdir(work_path)


            #create conf 
            confname = new_id[0]+'_mesh.conf'
            with open(confname, 'w', newline = '') as f:
                f.write("DIM = 1 1 1\n")
                f.write("ATOM_NAME = ")
                for z in new_cl.end.znucl:
                    el = element_name_inv(z)
                    f.write(el+' ')
                f.write("\nMP = {:}\n".format( mpstr ))
                f.write("\nTSTEP = {:}\n".format( 1 ))
                f.write("\nTMAX = {:}\n".format( 1155 ))

            if not os.path.exists("FORCE_SETS"):
                #run phonopy; read forces
                ndis = len( glob.glob('*.vasprun.xml') )
                print_and_log( ndis, ' displacement files was Found')
                print_and_log( runBash('export PYTHONPATH=~/installed/phonopy-1.9.5/lib/python:$PYTHONPATH; /home/dim/installed/phonopy-1.9.5/bin/phonopy '
                    +'  -f {1..'+str(ndis)+'}.vasprun.xml --tolerance=0.01'), imp = 'Y' )

            #calculate thermal prop
            result = 'thermal_properties_'+mpstr.replace(" ", "_")+'.yaml'
            if not os.path.exists(result):

                posname = 'SPOSCAR'
                print_and_log( runBash('export PYTHONPATH=~/installed/phonopy-1.9.5/lib/python:$PYTHONPATH; /home/dim/installed/phonopy-1.9.5/bin/phonopy '
                    +confname+' -c '+posname+' -t -p -s --tolerance=0.01'), imp = 'y' )

                shutil.copyfile('thermal_properties.yaml', result)
    

            T_range, fit_func = read_phonopy_data(result)

            os.chdir(savedPath)

        if type(new_id) == str:
            result = new_id
            T_range, fit_func = read_phonopy_data(result)



    return T_range, fit_func






def get_structure_from_cee_database(it, it_folder, ver, section = 'CEStorage', cee_struct_type = 'exp', cee_file = None):
    """
    cee_struct_type (str) - 
        'exp' - experimental structures
        '' - all

    """
    it_base = it.split('.')[0]
    print_and_log("Taking structure "+it_base+" from CEE CREI database of Skoltech ...", imp = 'Y')

    # database_server = 'aksenov@10.30.100.28'
    database_server = 'sd'
    # database_path   = '/home/Data/CEStorage/'
    database_path   = '/home/Data/'+section+'/'

    if 'exp' in cee_struct_type:
        templ = '*exp*.cif'
    else:
        templ = '*.cif'

    local_folder = it_folder+'/'+it+'/'

    makedir(local_folder)

    out = get_from_server(database_path+'/'+it_base+'/'+templ, local_folder, addr = database_server)

    geofiles = glob.glob(local_folder+templ)
    printlog(out, 'The following files have been downloaded', geofiles, imp ='Y'  )
    



    if len(geofiles) > 1:
        if cee_file:
            for file in geofiles:
                if cee_file in file :
                    geofile_from_server = file

        else:
            printlog('Error! More than one file, check what you need using *cee_file* parameter')
    else:
    
        geofile_from_server = geofiles[0]
    
    printlog('You have chosen', geofile_from_server)



    local_poscar_file = local_folder+it+'.POSCAR-'+str(ver)

    cif2poscar(geofile_from_server, poscar_file = local_poscar_file)

    add_des(header.struct_des, it, it_folder, des = 'taken automatically from cee_database: '+geofile_from_server)

    
    return






def get_structure_from_matproj(it = None, it_folder = None, ver = None, mat_proj_cell = '', mat_proj_id = None):
    """
    Take structures from Mat. projects

    Find material with 'it' stoichiometry (lowest energy) from materialsproject.org, 
    download and create field in struct_des and input POSCAR file
    ###INPUT:
        - struct_des-  
        - it        - materials name, such as 'LiCoO2', .... By default the structure with minimum *e_above_hull* is taken
        - it_folder - section folder in which the Poscar will be placed
        - ver       - version of structure defined by user
        - mat_proj_id (str) - the id can be provided explicitly
        - mat_proj_cell (str)- 
                - 'conv' - conventional

    ###RETURN:
        - ?
        - ?


    """
    with MPRester(header.pmgkey) as m:
        # print m.get_materials_id_references('mp-24850')
        # print m.get_structures('mp-24850')
        # mp_entries = m.get_entries_in_chemsys(["Co", "O"])
        # for e in mp_entries:
        #     if not 'is_hubbard = True' in e: continue
        # energy_per_atom
        # print mp_entries
        # print m.supported_task_properties
        # print it, 'it'
        # print type(it)
        # it = "LiFePO4"
        # print m.get_data(it, data_type='vasp', prop='e_above_hull')

        if mat_proj_id:
            groundstate_st_id = mat_proj_id
        else: 
            prop_dic_list =  m.get_data(it, data_type='vasp', prop='e_above_hull')




            newlist = sorted(prop_dic_list, key=itemgetter('e_above_hull')) 

            groundstate_st_id = newlist[0]['material_id']
            # print groundstate_st_id
            # print m.get_data(groundstate_st_id, data_type='vasp', prop='hubbards')
        
        st_pmg =  m.get_structure_by_material_id(groundstate_st_id, final=True)
        
        if 'conv' in mat_proj_cell:
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            sf = SpacegroupAnalyzer(st_pmg, symprec = 0.01)
            st_pmg = sf.get_conventional_standard_structure()

        # print(st_pmg.lattice)
        # sys.exit()

    if it:
        add_des(header.struct_des, it, it_folder, des = 'taken automatically from materialsproject.org: '+groundstate_st_id,)
        path2poscar = it_folder+'/'+it+'/'+groundstate_st_id+".POSCAR-"+str(ver)
        makedir(path2poscar)
    else:
        path2poscar = groundstate_st_id+".POSCAR"

    
    Poscar(st_pmg).write_file(path2poscar, direct=True, vasp4_compatible=True, )
    print_and_log('Structure', groundstate_st_id, 'downloaded from materialsproject.org\n',
        'File '+path2poscar+" was written", imp = 'y')

    st = smart_structure_read(input_geo_file = path2poscar)
    st.groundstate_st_id = groundstate_st_id
    st.mat_proj_st_id    = groundstate_st_id
    st.input_geo_file    = path2poscar


    return st

    #pymatgen drafts, can be useful

    # with MPRester(pmgkey) as m:
    #     print dir(m)
    #     print m.supported_properties
    #     print m.get_data('mp-540111', data_type='vasp', prop='total_magnetization')
    #     # 'total_magnetization'

    #Pymatgen symmetry analyzer
    # from pymatgen.io.vasp.inputs import Poscar
    # from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as SA
    # from pymatgen import Lattice, Structure, Molecule
    # a = Structure.from_file('/home/aksenov/scientific_projects/cathode/Li2FePO4F/scaled/Li2FePO4F.ir.su.4uis/100.CONTCAR' )
    # b = SA(a)
    # print(b.get_space_group_symbol())
    # print(a.get_space_group_info())
    # # lattice = Lattice(st.rprimd)
    # # print (lattice)

    # struct = Structure(st.rprimd, st.get_elements(), st.xred)
    # print(struct)








def manually_remove_from_struct_des(struct_des, key):
    """
    
    """
    del struct_des[key]
    print_and_log('Attention! Entry '+key+' was removed from struct_des dict/\n')












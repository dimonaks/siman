# from header import *
from __future__ import division, unicode_literals, absolute_import 


from small_functions import is_list_like
import copy, traceback, datetime, sys, os, glob, shutil
import header
from header import print_and_log, runBash
from classes import Calculation, CalculationVasp, Description
from functions import list2string, gb_energy_volume, element_name_inv, write_xyz, makedir, get_from_server, scale_cell_uniformly, image_distance
from picture_functions import plot_mep

import matplotlib.pyplot as plt
import numpy as np

from pymatgen.matproj.rest import MPRester
from pymatgen.io.vasp.inputs import Poscar

from operator import itemgetter

# sys.path.append('/home/aksenov/Simulation_wrapper/ase') #path to siman library
from ase.utils.eos import EquationOfState



log             = header.log
geo_folder      = header.geo_folder
cluster_address = header.cluster_address
corenum         = header.corenum
project_path_cluster = header.project_path_cluster
pmgkey          = header.project_conf.pmgkey
path_to_images = header.path_to_images




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
    with open(batch_script_filename,'w') as f:


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
            f.write("#PBS -l walltime=99999999:00:00 \n")
            f.write("#PBS -l nodes=1:ppn="+str(number_cores)+"\n")
            f.write("#PBS -r n\n")
            f.write("#PBS -j eo\n")
            f.write("#PBS -m bea\n")
            f.write("#PBS -M dimonaks@gmail.com\n")
            f.write("cd $PBS_O_WORKDIR\n")
            f.write("PATH=/share/apps/vasp/bin:/home/aleksenov_d/mpi/openmpi-1.6.3/installed/bin:/usr/bin:$PATH \n")
            f.write("LD_LIBRARY_PATH=/home/aleksenov_d/lib64:$LD_LIBRARY_PATH \n")


        if schedule_system == 'SLURM':
            if '~' in path_to_job:
                print_and_log('Error! For slurm std err and out you need full paths')
                raise RuntimeError
            f.write("#!/bin/bash   \n")
            f.write("#SBATCH -J "+job_name+"\n")
            f.write("#SBATCH -t 250:00:00 \n")
            f.write("#SBATCH -N 1\n")
            f.write("#SBATCH -n "+str(number_cores)+"\n")
            f.write("#SBATCH -o "+path_to_job+"sbatch.out\n")
            f.write("#SBATCH -e "+path_to_job+"sbatch.err\n")
            f.write("#SBATCH --mem-per-cpu=7675\n")
            f.write("#SBATCH --mail-user=d.aksenov@skoltech.ru\n")
            f.write("#SBATCH --mail-type=END\n")
            f.write("cd "+path_to_job+"\n")
            f.write("export OMP_NUM_THREADS=1\n")

            f.write("module add prun/1.0\n")
            f.write("module add intel/16.0.2.181\n")
            f.write("module add impi/5.1.3.181\n")
            f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/tools/lib64:~/tools/atlas\n")
            f.write("export PATH=$PATH:~/tools\n")
            f.write("touch RUNNING\n")


    return























def push_figure_to_archive(local_figure_path, caption, figlabel = None ):
    shutil.copy(local_figure_path, header.project_conf.path_to_images)
    print_and_log('push_figure_to_archive():', local_figure_path, 'copied to', header.project_conf.path_to_images, imp = 'y')
    
    name_without_ext =   '.'.join( os.path.basename(local_figure_path).split('.')[:-1]) 
    figfile = '{{'+name_without_ext+'}}'

    if not figlabel:
        figlabel = '.'.join(name_without_ext.split('.')[:-1])
    

    tex_text = \
    ("\\begin{{figure}} \n\includegraphics[width=\columnwidth]{{{:s}}}\n"
    "\caption{{\label{{fig:{:s}}} {:s} }}\n"
    "\end{{figure}}\n").format(figfile, figlabel, caption+' for '+figlabel )
    # print (tex_text)
    with open(header.project_conf.path_to_paper+'/auto_fig.tex', 'a+') as f:
        if tex_text not in f.read():
            f.write(tex_text)
    return


def clean_history_file(history_list):
    seen = set()
    seen_add = seen.add
    return [x for x in history_list if not (x in seen or seen_add(x))]





def clean_run(schedule_system = None):
    """
    INPUT:
        schedule_system - type of job scheduling system:'PBS', 'SGE', 'SLURM'
    """
    with open('run','w') as f:
    
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

        with open('run','a') as f:
            if header.project_conf.SCHEDULE_SYSTEM == "PBS":
                f.write("qstat\n")
                f.write("sleep 2\n")
            elif header.project_conf.SCHEDULE_SYSTEM == "SLURM":
                f.write("squeue\n")

            f.write("mv run last_run\n")

        

        runBash('chmod +x run')

        log.write( runBash("rsync -zave ssh run "+cluster_address+":"+project_path_cluster) +"\n" )
        print_and_log('run sent')
        # clean_run(header.project_conf.SCHEDULE_SYSTEM)
    
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

    if it in struct_des and not override:
        print_and_log("Error! "+it+' already exist in struct_des; use override = True if you really need it; or first remove using manually_remove_from_struct_des()')
        raise RuntimeError
    else:
        struct_des[it] = Description(it_folder, des)
        # hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
        hstring = 'add_des("{:s}", "{:s}", "{:s}")    #on {:})'.format(it, it_folder, des, datetime.date.today())

        try:
            if hstring != header.history[-1]: header.history.append( hstring  )
        except:
            header.history.append( hstring  )
        print_and_log("New structure name "+it+ " added to struct_des dict"+"\n")






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




def smart_structure_read(curver, inputset = '', cl = None, input_folder = None, input_geo_format = None, input_geo_file = None):
    """
    Wrapper for reading geometry files
    Also copies geofile and OCCMATRIX
    returns Structure()
    """

    if input_geo_file:
        geofilelist = glob.glob(input_geo_file) 
        print_and_log("You provided the following geo file explicitly "+str(geofilelist)+"\n" )
    
    else:
        print_and_log("I am searching for geofiles in folder "+input_folder+"\n" )
        
        if input_geo_format == 'abinit': 
            searchinputtemplate = input_folder+'/*.geo*'
        
        elif input_geo_format == 'vasp': 
            searchinputtemplate = input_folder+'/*POSCAR*'



        elif input_geo_format == 'cif': 
            searchinputtemplate = input_folder+'/*.cif'

        # print 'searchinputtemplate = ', searchinputtemplate

        # print  input_geo_format
        geofilelist = glob.glob(searchinputtemplate) #Find input_geofile
        # print geofilelist
    # print os.path.basename(file[0])
    geofilelist = [file for file in geofilelist if os.path.basename(file)[0] != '.'   ]  #skip hidden files


    #additional search in target folder if no files in root # !!!Add for Vasp also 
    if not geofilelist:
        print_and_log("Attention! trying to find here "+input_folder+"/target\n" )
        geofilelist = glob.glob(input_folder+'/target/*.geo*') #Find input_geofile            

    if not geofilelist:
        input_folder += '.'+inputset
        print_and_log("Attention! trying to find here "+input_folder+"\n" )
        geofilelist = glob.glob(input_folder+'/*.geo*') #Find input_geofile    


    for input_geofile in geofilelist: #quite stupid to have this loop here - much better to move this to upper function, and the loop will not be needed
        
        #print runBash("grep version "+str(input_geofile) )
        input_geofile = os.path.normpath(input_geofile)
        if input_geo_format in ['abinit',]:
            curv = int( runBash("grep version "+str(input_geofile) ).split()[1] )

        elif input_geo_format == 'vasp': 
            curv = int(input_geofile.split('-')[-1] ) #!Applied only for phonopy POSCAR-n naming convention

        elif input_geo_format == 'cif': 
            curv = int(os.path.basename(input_geofile).split('.')[0] )


        if curv == curver:

            if cl:
                calc_geofile_path = os.path.normpath(cl.dir + input_geofile.split('/')[-1])
                
                if input_geofile != calc_geofile_path: # copy initial geo file and other files to calc folder

                    makedir(calc_geofile_path)

                    shutil.copyfile(input_geofile, calc_geofile_path)
                    dir_1 = os.path.dirname(input_geofile)
                    dir_2 = os.path.dirname(calc_geofile_path) 
                    
                    if 'OCCEXT' in cl.set.vasp_params and cl.set.vasp_params['OCCEXT'] == 1:

                        shutil.copyfile(dir_1+'/OCCMATRIX', dir_2+'/OCCMATRIX' )
            else:
                cl = Calculation()


            if input_geo_format == 'abinit':
                cl.read_geometry(input_geofile)
            
            elif input_geo_format == 'vasp':
                cl.read_poscar(input_geofile)

            elif input_geo_format == 'cif':
                if header.project_conf.CIF2CELL:
                    print_and_log( runBash("cif2cell "+input_geofile+"  -p vasp -o "+input_geofile.replace('.cif', '.POSCAR'))  )
                    input_geofile = input_geofile.replace('.cif', '.POSCAR')
                    
                    #check
                    if not os.path.exists(input_geofile):
                        print_and_log("Error! Something wrong with conversion of cif2cell: \n")
                        raise RuntimeError


                    cl.read_poscar(input_geofile)
                


                else:
                    print_and_log("Error! cif2cell is not available in your system")
                    raise RuntimeError

            
            else:
                raise RuntimeError

            

            break
    
    if cl and cl.path["input_geo"] == None: 
        print_and_log("Error! Could not find geofile in this list: "+ str(geofilelist)+  "\n")
        raise NameError #
    

    return cl.init





def add_loop(it, setlist, verlist, calc = None, conv = None, varset = None, 
    up = 'up1', typconv="", from_geoise = '', inherit_option = None, 
    coord = 'direct', savefile = 'ov', show = None, comment = '', 
    input_geo_format = 'abinit', ifolder = None, input_geo_file = None, corenum = None,
    calc_method = None, u_ramping_region = None, it_folder = None, mat_proj_id = None, ise_new = None,
    scale_region = None, n_scale_images = 7, id_from = None,
    n_neb_images = None, occ_atom_coressp = None,ortho = None,
    ):
    """
    Main subroutine for creation of calculations, saving them to database and sending to server.

    Input:
        - it - arbitary name for your crystal structure 
        - setlist (list of str or str) - names of sets with vasp parameters from *varset* dictionary
        - verlist - list of versions of new calculations
        - calc, conv, varset - database dictionaries; could be provided; if not then are taken from header

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

        - it_folder - section folder (sfolder) used in struct_des; here needed with input_geo_format = mat_proj

        - show - only for read_results() ?.

        - comment - arbitrary comment for history.


        #inherit flags:
        - inherit_option (str):
            - 'continue'     - copy last contcar to poscar, outcar to prev.outcar and run again; on the next launch prev.outcar
                will be rewritten, please improve the code to save all previous outcars 
            - 'inherit_xred' - if verlist is provided, xred are copied from previous version to the next
            - all options available for inherit_icalc() subroutine (now only 'full' is tested)
        - *id_from* - see inherit_icalc()
        - ise_new (str) - name of new set for inherited calculation  ('uniform_scale')

        - occ_atom_coressp (dict) see inherit_icalc()


        - typconv - ? to be described.
        
        - 'from_geoise' - part of folder name with geometry input files. allows to use geometry files from different sets.
        please find how it was used

        - corenum - number of cores used for calculation; overwrites header.corenum

        - calc_method - provides additional functionality:
            - 'u_ramping'    - realizes U ramping approach #Phys Rev B 82, 195128
            - 'afm_ordering' - 
            - 'uniform_scale' - creates uniformly scaled copies of the provided calculations
            using *scale_region* and *n_scale_images* (see *scale_cell_uniformly()*)
            The copies are available as versions from 1 to *n_scale_images* and
            suffix .su appended to *it* name
            Copies to cluster *fit* utility that finds volume corresp. to energy minimum, creates 100.POSCAR and continues run 


        - u_ramping_region - used with 'u_ramping'=tuple(u_start, u_end, u_step)


    Comments:
        !Check To create folders and add calculations add_flag should have value 'add' 


    TODO:
    Now number of images is taken from self.set.vasp_params['IMAGES']; 
    In the case of generalization to other codes, set.nimages should be added and used

    """

    header.close_run = True

    schedule_system = header.project_conf.SCHEDULE_SYSTEM

    if corenum:
        # corenum = ppn
        ''
    else:
        corenum = header.CORENUM

    struct_des = header.struct_des

    # print type('NaFePO4.pnma')
    # print struct_des['NaFePO4.pnma']

    if not calc:
        calc = header.calc
        conv = header.conv
        varset = header.varset



    it = it.strip()
    
    if it_folder: 
        it_folder = it_folder.strip()

    if not is_list_like(verlist):
        verlist = [verlist]

    if not is_list_like(setlist):
        setlist = [setlist]

    setlist = [s.strip() for s in setlist]


    if not is_list_like(calc_method):
        calc_method = [calc_method]






    if ifolder: 
        if it not in ifolder: # just to be consistent with names
            print_and_log('Check ifolder !!! it is not in ifolder')
            raise RuntimeError




    





    if typconv: 
        setlist = varset[ise].conv[typconv] #
        nc = it+'.'+ise[0]+typconv
        if nc not in conv: 
            conv[nc] = []    
    

    if up == "no_base": 
        setlist = varset[ise].conv[typconv][1:]; 
        up = "up1"
    







    mat_proj_st_id = None
    if input_geo_format == 'mat_proj':
        print_and_log("Taking structure "+it+" from materialsproject.org ...\n")
        if it_folder == None:
            raise RuntimeError



        mat_proj_st_id, input_geo_file = get_structure_from_matproj(struct_des, it, it_folder, fv, mat_proj_id)
        input_geo_format = 'vasp'

    neb_flag = calc_method and not set(['neb', 'only_neb']).isdisjoint(calc_method)
    # print neb_flag
    # print calc_method
    # sys.exit()

    if neb_flag: #put nimage values for set_sequence
        curset = varset[ setlist[0] ]
        if not n_neb_images:
            n_neb_images = varset[curset.vasp_params['IMAGES']]

        if not n_neb_images:
            print_and_log('Error! You did not provide number of NEB images nor in *n_neb_images* nor in your set!')
            raise RuntimeError


        if corenum % n_neb_images > 0:
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



    if calc_method and 'uniform_scale' in calc_method:


        it_new = it+'.su'
        v = verlist[0]

        # if up != 'up3':
        print_and_log('Preparing   uniform_scale  calculation ... ')

        if len(verlist) > 1:
            print_and_log('Error! Currently   uniform_scale  is allowed only for one version')
            raise RuntimeError
        


        if it_new not in struct_des:
            if it_folder:
                section_folder = it_folder
            else:
                section_folder = struct_des[it].sfolder

            add_des(struct_des, it_new, section_folder, 'uniform_scale: scaled "images" for '+it+'.'+str(setlist)+'.'+str(v)   )


        verlist_new = []

        if ise_new and len(setlist) > 1:
            raise RuntimeError

        for inputset in setlist:


            id_s = (it,inputset,v)
            if id_s in calc:
                st = calc[id_s].end
                pname = str(id_s)
            else:
                st = smart_structure_read(curver = v, inputset = inputset, cl = None, input_folder = struct_des[it].sfolder+'/'+it, 
                    input_geo_format = input_geo_format, input_geo_file = input_geo_file)
                pname = st

            write_xyz(st, file_name = st.name+'_used_for_scaling')
            sts = scale_cell_uniformly(st, scale_region = scale_region, n_scale_images = n_scale_images, parent_calc_name = pname)
            


            if ise_new:
                inputset = ise_new
                id_s = (it,inputset,v)

            cl_temp = CalculationVasp(varset[inputset], id_s)

            for i, s in enumerate(sts):
                ver_new = i+1
                s.name = it_new+'.'+s.name
                cl_temp.init = s
                cl_temp.version = ver_new
                cl_temp.path["input_geo"] = geo_folder + struct_des[it_new].sfolder + '/' + \
                                            it_new+"/"+it_new+'.auto_created_scaled_image'+'.'+str(ver_new)+'.'+'geo'

                cl_temp.write_siman_geo(geotype = "init", 
                    description = s.des, override = True)
                write_xyz(s)
                verlist_new.append(ver_new)
            #make version 100
            cl_temp.version = 100
            cl_temp.des = 'fitted with fit_tool.py on cluster, init is incorrect'
            cl_temp.id = (it_new, inputset, 100)
            cl_temp.state = '2.ready to read outcar'
            blockdir = struct_des[it_new].sfolder+"/"+varset[inputset].blockfolder #calculation folder
            # iid = cl_temp.id          
            cl_temp.name = cl_temp.id[0]+'.'+cl_temp.id[1]+'.'+str(cl_temp.id[2])
            cl_temp.dir = blockdir+"/"+ str(cl_temp.id[0]) +'.'+ str(cl_temp.id[1])+'/'
            cl_temp.path["output"] = cl_temp.dir+str(cl_temp.version)+'.OUTCAR'
            cl_temp.cluster_address      = header.cluster_address
            cl_temp.project_path_cluster = header.project_path_cluster
            calc[cl_temp.id] = cl_temp
            # cl_temp.init = None

            verlist = verlist_new

            print_and_log(len(sts), 'uniform images have been created.')
        




        it      = it_new
        if ise_new:
            setlist = [ise_new]
        # sys.exit()



    write_batch_list = [True for v in verlist] #which version should be written in batch_script - not used now





    #inherit option
    if inherit_option in ['supercell', 'occ', 'full', 'full_nomag', 'r2r3', 'r1r2r3', 'remove_imp', 'replace_atoms', 'make_vacancy',]:
        if inherit_option == 'full':
            it_new = it+'.if'
        if inherit_option == 'full_nomag':
            it_new = it+'.ifn'

        if inherit_option == 'occ':
            #please add additional vars to control for which atoms the inheritance should take place
            it_new = it+'.ifo' #full inheritence + triggering OMC        
        if inherit_option == 'supercell':
           if len(set(ortho))==1:
                mod = '.s'+str(ortho[0])
           else:
                mod =  '.s'+list2string(ortho).replace(' ','')
           it_new = it+mod






        if it_folder:
            section_folder = it_folder
        else:
            section_folder = struct_des[it_new].sfolder




        # if up != 'up3':
        if it_new not in struct_des:
            add_des(struct_des, it_new, section_folder, 'Inherited '+inherit_option+' from '+it+'.'+str(setlist)+'.'+str(verlist)   )

        for inputset in setlist:
            for v in verlist:
                id = (it,inputset,v)
                # print (calc[id].end.magmom)
                # sys.exit()
                
                inherit_icalc(inherit_option, it_new, v, id, calc, id_from = id_from, it_folder = it_folder, occ_atom_coressp = occ_atom_coressp,ortho = ortho)
        



        if ise_new:
            if up != 'up3':
                print_and_log('Inherited calculation uses set', ise_new)

            setlist = [ise_new,]

        else:
            if up != 'up3':
                print_and_log('Inherited calculation uses the same sets', setlist)




        it = it_new




    if 0:
        hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
        args = hstring.split('(')[1].split(',')
        hstring = hstring.replace(args[0], "'"+it+"'")
        hstring = hstring.replace(args[1], str(setlist))
    else: #more useful and convenient
        hstring = "res_loop('{:s}', {:s}, {:s}, show = 'fo'  )     # {:s}, on {:s}  ".format(
            it, str(setlist), str(verlist), comment, str(datetime.date.today() )  )
    # try:
    if hstring != header.history[-1]: 
        header.history.append( hstring  )
    # except:
    #     header.history.append( hstring  )








    """Main Loop by setlist and verlist"""
    fv = verlist[0]; #first version
    lv = verlist[-1];#last version

    output_files_names = []


    for inputset in setlist:

        if ifolder:
            input_folder = ifolder

        else:
            if from_geoise:
                from_geoise = from_geoise[0]+inputset[1:] #it is supposed that the difference can be only in first digit
                input_folder = geo_folder+it+"/" + it+"."+from_geoise #+ "/" + "grainA_s" #geo used for fitted
            else: 
                input_folder = geo_folder+struct_des[it].sfolder+"/"+it


        prevcalcver = None #version of previous calculation in verlist

        for v in verlist:
            id = (it,inputset,v)
            
            if typconv and id not in conv[nc]: 
                conv[nc].append(id)
            
            try: 
                blockfolder = varset[inputset].blockfolder
            except AttributeError: 
                blockfolder = varset[ise].blockfolder
            
            blockdir = struct_des[it].sfolder+"/"+blockfolder #calculation folder

            add_calculation(it,inputset,v, fv, lv, input_folder, blockdir, calc, varset, up, 
                inherit_option, prevcalcver, coord, savefile, input_geo_format, input_geo_file, 
                schedule_system = schedule_system, 
                calc_method = calc_method, u_ramping_region = u_ramping_region,
                mat_proj_st_id = mat_proj_st_id,
                output_files_names = output_files_names,
                corenum = corenum
                )
            
            prevcalcver = v





    # if calc_method and 'neb' in calc_method:
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










    if up not in ('up1','up2','up3'): 
        print_and_log("Warning! You are in the test mode, to add please change up to up1; "); 
        sys.exit()
    return it












def add_calculation(structure_name, inputset, version, first_version, last_version, input_folder, blockdir, 
    calc, varset, update = "no",
    inherit_option = None, prevcalcver = None, coord = 'direct', savefile = None, input_geo_format = 'abinit', 
    input_geo_file = None, schedule_system = None, calc_method = None, u_ramping_region = None,
    mat_proj_st_id = None, output_files_names = None, corenum = 1):
    """

    schedule_system - type of job scheduling system:'PBS', 'SGE', 'SLURM'

    prevcalcver - version of previous calculation in verlist

    output_files_names - the list is updated on every call

    if inherit_option == 'continue' the previous completed calculation is saved in cl.prev list

    """
    struct_des = header.struct_des


    id = (structure_name,inputset,version)

    cl_prev = None


    if id in calc: 
        cl = calc[id]
        status = "exist"
        if update != 'up3': #
            print_and_log(str(calc[id].name)+" has been already created and has state: "+str(calc[id].state)+"\n\n")

        if "4" in calc[id].state: 
            complete_state = calc[id].state
            status = "compl"

            if update == 'up2': 
                print_and_log( 'Calculation', calc[id].name, 'is finished, continue')

                return

            if update != "up1": 
                return #completed calculations updated only for "up1"
        
        if status == 'exist' and update == 'up3':
            return #

    else:
        #update = "up"
        status = "new"
        print_and_log( "There is no calculation with id "+ str(id)+". I create new with set "+str(inputset)+"\n" )        

    if "up" in update:

        if status in ["exist","compl"]: 
            print_and_log("You asked to update existing calculation with id "+ str(id)+" Warning! I update creating new class \n" )         

        if status == 'compl' and inherit_option == 'continue':
            print_and_log(id, 'is completed, I will make its copy in self.prev[]', imp = 'Y' )         

            cl_prev = copy.deepcopy(calc[id])

        calc[id] = CalculationVasp( varset[id[1]] )
        cl = calc[id]

        calc[id].id = id 
        calc[id].name = str(id[0])+'.'+str(id[1])+'.'+str(id[2])
        calc[id].dir = blockdir+"/"+ str(id[0]) +'.'+ str(id[1])+'/'
        

        batch_script_filename = cl.dir+cl.id[0]+"."+cl.id[1]+'.run'        

        # all additional properties:
        calc[id].calc_method = calc_method


        if hasattr(calc[id].set, 'u_ramping_nstep') and calc[id].set.u_ramping_nstep:
            print_and_log("Attention! U ramping method is detected from set\n\n")
            calc[id].calc_method.append('u_ramping')

        if hasattr(calc[id].set, 'afm_ordering'):
            print_and_log("Attention! afm_ordering method is detected from set\n\n")
            calc[id].calc_method.append('afm_ordering')




        calc[id].cluster_address = header.CLUSTER_ADDRESS
        calc[id].project_path_cluster = header.project_path_cluster
        
        calc[id].corenum = corenum #this is correct - corenum provided to functions
        calc[id].schedule_system = schedule_system


        if mat_proj_st_id:
            calc[id].mat_proj_st_id = mat_proj_st_id




        if inherit_option == 'continue' and cl_prev:
            if hasattr(cl_prev, 'prev') and cl_prev.prev:
                calc[id].prev.extend(cl_prev.prev) #if 'continue' flag is used several times, we do not need cl.prev.prev.prev.... but cl.prev = [cl1, cl2 ..]
            else:
                calc[id].prev.append(cl_prev)


        if update in ['up1', 'up2', 'up3']:
            if not os.path.exists(calc[id].dir):
                log.write( runBash("mkdir -p "+calc[id].dir) )         #Create directory if it does not exist
                log.write( runBash("ssh "+calc[id].cluster_address+" ' mkdir -p "+calc[id].dir+" ' ") )
            
            if id[2] == first_version:
                write_batch_header(batch_script_filename = batch_script_filename,
                    schedule_system = cl.schedule_system, 
                    path_to_job = header.PATH_TO_PROJECT_ON_CLUSTER+cl.dir, 
                    job_name = cl.id[0]+"."+cl.id[1], number_cores = cl.corenum  )




                
        cl.init = smart_structure_read(curver = cl.id[2], inputset = inputset, cl = cl, input_folder = input_folder, 
            input_geo_format = input_geo_format, input_geo_file = input_geo_file)






        calc[id].des += ' '+struct_des[id[0]].des + '; ' + varset[id[1]].des


        setlist = [cl.set]                                                                                                    
        if hasattr(cl.set, 'set_sequence') and cl.set.set_sequence:
            for s in cl.set.set_sequence:
                setlist.append(s)
        calc[id].check_kpoints()    
        for curset in setlist:
            calc[id].actualize_set(curset)




        if update in ['up1', 'up2', 'up3']:
            
            calc[id].write_structure(str(id[2])+".POSCAR", coord, inherit_option, prevcalcver)
            
        
            out_name = calc[id].write_sge_script(str(version)+".POSCAR", version, 
                inherit_option, prevcalcver, savefile, 
                schedule_system = schedule_system, mode = 'body',
                batch_script_filename = batch_script_filename)
            
                        
            if out_name:
                calc[id].path["output"] = calc[id].dir+out_name
            else:
                name_mod = ''
                calc[id].path["output"] = calc[id].dir+str(version)+name_mod+".OUTCAR" #set path to output
            
            output_files_names.append( calc[id].path["output"] )




            if id[2] == first_version:
                calc[id].add_potcar()
            for curset in setlist:
                calc[id].calculate_nbands(curset)



            if id[2] == last_version:
                list_to_copy = []
                
                calc[id].write_sge_script(mode = 'footer', schedule_system = schedule_system, option = inherit_option, 
                    output_files_names = output_files_names, batch_script_filename = batch_script_filename )
                
                runBash('chmod +x '+batch_script_filename)

                list_to_copy.extend( cl.make_incar() )
                
                list_to_copy.extend( cl.make_kpoints_file() )
                
                cl.copy_to_cluster(list_to_copy, update)

                calc[id].make_run(schedule_system = schedule_system)




        if status == "compl": 
            calc[id].state = '2. Can be completed but was reinitialized' #new behavior 30.08.2016


        print_and_log("\nCalculation "+str(id)+" added or updated\n\n")

    return













def inherit_icalc(inherit_type, it_new, ver_new, id_base, calc = None,
    id_from = None,
    atom_new = None, atom_to_replace = None,  id_base_st_type = 'end', atoms_to_remove = None, i_atom_to_remove = None, id_from_st_type = 'end',
    atom_to_shift = None, shift_vector = None,
    it_folder = None, occ_atom_coressp = None, ortho = None,
    ):
    """
    Function for creating new geo files in geo folder based on different types of inheritance
    Input args: 
        it_new, ver_new - name of new structure,
        id_base - new structure will be based on the final structure of this calculation;     (can be either Calculation() object or path to geo file)
        id_from - can be additionally used to adopt for example rprimd from id_from to it_new; (can be either Calculation() object or path to geo file)

        
        inherit_type = '':
            full          - full inheritance of final state
            'full_nomag'  - full except magmom which are set to None
            r2r3          - use r2 and r3 from id_from
            r1r2r3        - use r1, r2 and r3 from id_from
            remove_atoms  - removes atoms of type *atoms_to_remove (list of str)*
            replace_atoms - atoms of type 'atom_to_replace' in 'id_base' will be replaced by 'atom_new' type.
            make_vacancy  - produce vacancy by removing 'i_atom_to_remove' starting from 0
            occ           - take occ from *id_from* and create file OCCMATRIX for 
                            OMC [https://github.com/WatsonGroupTCD/Occupation-matrix-control-in-VASP]
                            - occ_atom_coressp (dict) {iatom_calc_from:iatom_calc_base, ... } (atomno starting from 0!!!)
            supercell - create orthogonal supercel using ortho list [a,b,c]
        id_base_st_type - use init or end structure of id_base calculation.
        id_from_st_type  - init or end for id_from

        atom_to_shift - number of atom to be shifted; starting from 1.
        shift_vector - vector in decart cooridinates (Angstrom!!) by which the atom will be shifted

        - it_folder - section folder



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
        take care of magmom, not implemented
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
    override  = False
    if type(id_base) == str:
        print_and_log('Reading id_base\n')
        cl_base = CalculationVasp()
        cl_base.read_geometry(id_base)
        cl_base.id = ('from_file', 'from_file', cl_base.version)
        cl_base.name = id_base
        cl_base.end = cl_base.init

    else:
        if id_base not in calc:
            ''
            id_base = (bytes(id_base[0]),bytes(id_base[1]),id_base[2])
        cl_base = calc[id_base]


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



    new.version = ver_new

    if id_base_st_type == 'init':
        st = new.init
    elif id_base_st_type == 'end':
        st = new.end



    #path to new calc
    if it_folder:
        section_folder = it_folder
    else:
        section_folder = struct_des[it_new].sfolder


    it_new_folder = geo_folder + section_folder + '/' + it_new
    new.path["input_geo"] = it_new_folder + '/' +it_new+'.inherit.'+inherit_type+'.'+str(ver_new)+'.'+'geo'

    makedir(new.path["input_geo"])
    print_and_log('Path for inherited calc =', it_new_folder)





    if inherit_type == "r2r3":
        des = ' Partly inherited from the final state of '+cl_base.name+'; r2 and r3 from '+calc_from_name
        # new.des = struct_des[it_new].des + des
        st.rprimd[1] = st_from.rprimd[1].copy()
        st.rprimd[2] = st_from.rprimd[2].copy()       
        # new.write_geometry("end",des)

    elif inherit_type == "r1r2r3":
        des = ' Partly inherited from the final state of '+cl_base.name+'; r1, r2, r3 from '+calc_from_name
        # new.des = struct_des[it_new].des + des
        st.rprimd = copy.deepcopy( st_from.rprimd )
        new.hex_a = calc_from.hex_a
        new.hex_c = calc_from.hex_c
        st.xcart = xred2xcart(new.end.xred, new.end.rprimd) #calculate new xcart from xred, because rprimd was changed
        # new.write_geometry("end",des)


    elif inherit_type == "full":
        print_and_log("Warning! final xred and xcart was used from OUTCAR and have low precision. Please use CONTCAR file \n");
        des = 'Fully inherited from the final state of '+cl_base.name
        # new.des = des + struct_des[it_new].des
        # new.write_geometry("end",des)

    elif inherit_type == "full_nomag":
        # print_and_log("Warning! final xred and xcart was used from OUTCAR and have low precision. Please use CONTCAR file \n");
        des = 'Fully inherited from the final state of '+cl_base.name+'; "magmom" set to [None]'
        # new.des = des + struct_des[it_new].des
        st.magmom = [None]
        # new.write_geometry("end",des)

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

        #create OCCMATRIX 
        print_and_log('I create OCCMATRIX in ', it_new_folder)
        with open(it_new_folder+'/OCCMATRIX', 'w') as f:
            numat = len(occs)
            f.write(str(numat)+'  #num of atoms to be specified\n')
            
            at_nums = occs.keys()
            at_spin = [] # # 2 or 1
            at_ltyp = [] # l - orbital type, 1 - s, 2 - d, 3 - f
            for key in occs: 
                occ = occs[key]
                if len(occ) == 10: # spin polarized, d orbital
                    at_spin.append(2)
                    at_ltyp.append(2)
                else:
                    raise RuntimeError # please write by yourself for other cases


            for i, l, s in zip(at_nums, at_spin, at_ltyp):

                f.write(list2string([i+1, l, s])+'    #i, l, s\n')
                # for sp in range(s):
                f.write('spin 1\n')
                for row in occs[i][ 0:len(occs[i])//s ]:
                    f.write(list2string(row)+'\n')
                if s == 2:
                    f.write('spin 2\n')
                    for row in occs[i][ len(occs[i])//s: ]:
                        f.write(list2string(row)+'\n')
                f.write('\n')

        # st.magmom = [None]
        



        # sys.exit()

    elif inherit_type == 'supercell':
        from geo import ortho_vec, create_supercell
        mul_matrix = ortho_vec(st.rprimd, ortho_sizes = ortho)
        print_and_log('Mul matrix is\n',mul_matrix)
        sc = create_supercell(st, mul_matrix)
        new.init = sc
        new.end  = sc
        des = 'obtained from'+cl_base.name+'by creating supercell'+str(ortho)
        override = True
    elif inherit_type == "atom_shift":
        des = 'obtainded from final state of '+cl_base.name+' by shifting atom '+ str(atom_to_shift) +' by '+ str(shift_vector)
        # new.des = des + struct_des[it_new].des
        
        st.xcart[atom_to_shift-1] += np.asarray(shift_vector) 
        st.xred = xcart2xred(new.end.xcart, new.end.rprimd)
        # new.write_geometry("end",des)
        # write_xyz(new.end)



    elif inherit_type == "remove_atoms":
        """
        
        """
        des = 'All atoms of type' + str(atoms_to_remove)+' removed from the final state of '+cl_base.name
        
        atoms = [ element_name_inv(st.znucl[t-1])    for t in st.typat ]

        # print (atoms)
        atom_exsist = True
        
        while atom_exsist:
            atoms = [ element_name_inv(st.znucl[t-1])    for t in st.typat ]

            for i, at in enumerate(atoms):
                
                if at in atoms_to_remove:
                    
                    st = st.del_atoms(i)

                    break
            else:
                atom_exsist = False
            # print (atoms)

        
        st.name = it_new+'_from_'+new.name
        override = True

        # sys.exit()
     
        # st_copy = copy.deepcopy(st)
        # st.typat = []
        # st.xred = []
        # st.xcart = []
        # st.ntypat = 1
        # st.znucl = st.znucl[0:1]
        # for i, t in enumerate(st_copy.typat):
        #     if t == 1:
        #         st.typat.append(t)
        #         st.xred.append(st_copy.xred[i])
        #         st.xcart.append(st_copy.xcart[i])
        # st.natom = len(st.xred)
        



        # new.init = new.end #just for sure
        # new.write_geometry("end",des)        
        # write_xyz(new.end)



    elif inherit_type == "make_vacancy":
        """Remove  atom 'i_atom_to_remove' from final state of id_base"""


        print_and_log('Warning! Please check inherit_type == "make_vacancy", typat can be wrong  if more than one element present in the system\n ',
            'Use del_atoms() method ')
        raise RuntimeError

        des = 'Atom '+str(i_atom_to_remove)+' removed from  '+cl_base.name
        # new.des = des + struct_des[it_new].des

        # st = new.end
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
            print_and_log("Error! Something wrong with names of atom types")
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

        # new.des = des #+ struct_des[it_new].des
        # write_xyz(st)

    else:
        print_and_log("Error! Unknown type of Calculation inheritance"); raise RuntimeError



    


    #auto addition of description
    if it_new not in struct_des: 
        add_des(struct_des, it = it_new, it_folder = it_folder, des = 'auto '+des)
        new.des =  struct_des[it_new].des
    else:
        new.des = des + struct_des[it_new].des

    #write files

    # print new.end.xcart

    #print(len(new.end.xred))
    new.write_geometry(id_base_st_type, des, override = override)
    write_xyz(st)





    return




























def res_loop(it, setlist, verlist,  calc = None, conv = {}, varset = {}, analys_type = 'no', b_id = (), 
    typconv='', up = "", imp1 = None, imp2 = None, matr = None, voronoi = False, r_id = None, readfiles = True, plot = True, show = '', 
    comment = None, input_geo_format = None, savefile = None, energy_ref = 0, ifolder = None, bulk_mul = 1, inherit_option = None,
    calc_method = None, u_ramping_region = None, input_geo_file = None,
    it_folder = None, choose_outcar = None, choose_image = None, mat_proj_id = None, ise_new = None, push2archive = False,
    description_for_archive = None, old_behaviour  = False,
    alkali_ion_number = None):
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

            )
        voronoi - True of False - allows to calculate voronoi volume of impurities and provide them in output. only if lammps is installed
        b_id - key of base calculation (for example bulk cell), used in several regimes; 
        r_id - key of reference calculation; defines additional calculation (for example atom in vacuum or graphite to calculate formation energies); can contain directly the energy per one atom

        up - controls if to download files from server; can be 'xo'
        
        readfiles (bool) - True - read from outcar, False - read from database; 



        The next three used for 'clusters' regime:    
        imp1 - key of bulk cell with one imp1
        imp2 - key of bulk cell with one imp2
        matr - key of bulk cell with pure matrix.


        - show - (str), allows to show additional information:
            - mag - magnetic moments on magnetic atoms
                *alkali_ion_number* (int) - number of atom around which to sort mag moments 
            - en  - convergence of total energy vs max force
            - mep - neb path
            - fo  - max force on each md step

        energy_ref - energy in eV; substracted from energy diffs
        
        bulk_mul - allows to scale energy and volume of bulk cell during calculation of segregation energies

        choose_outcar (int, starting from 1)- if calculation have associated outcars, you can check them as well, by default
        the last one is used during creation of calculation in write_sge_script()

        choose_image (int) - relative to NEB, allows to choose specific image for analysis, by default the middle image is used

        - push2archive (bool) - if True produced images are copied to header.project_conf.path_to_images
        - description_for_archive - caption for images



        - ise_new - dummy
        - inherit_option - dummy
        - savefile - dummy


    RETURN:
        result_list - list of results

    TODO:
        Make possible update of b_id and r_id with up = 'up2' flag; now only id works correctly


    """


    """Setup"""
    if not is_list_like(verlist):
        verlist = [verlist]

    if not is_list_like(setlist):
        setlist = [setlist]
        # print (setlist)


    if not calc:
        calc = header.calc

    try:
        b_ver_shift = b_id[2] #add to version of base with respect to version of main
    except:
        b_ver_shift = 0

    if '2' in up:
        loadflag = 'o'
    else:
        loadflag = ''

    # if choose_outcar:


    name_field_length = 30





    if typconv == '': pass
    else: setlist = varset[setlist[0]].conv[typconv] #


    n = 'temp'; conv[n] = []
    base = 'base'; conv[base] = []
    conv[it] = []
    #print calc[b_id]

    result_list = []


    emin = 0
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
    # print (setlist)
    final_outstring = 'no calculation found'
    for inputset in setlist:
        for v in verlist:
            # print 'Starting loops'

            id = (it,inputset,v)
            # print(id)
            if id not in calc:
                id = (bytes(it, 'utf-8'), bytes(inputset, 'utf-8'), v) #try non-unicode for compatability with python2
                if id not in calc:
                    print_and_log('Key', id,  'not found!', imp = 'Y')
                    continue #pass non existing calculations

            cl = calc[id]


            if not hasattr(cl,'version'):
                calc[id].version = v


            outst = ' File was not read '
            
            if readfiles:

                    outst = calc[id].read_results(loadflag, analys_type, voronoi, show, 
                        choose_outcar = choose_outcar, alkali_ion_number = alkali_ion_number)



            if analys_type in ('e_seg', 'coseg'):
                try:
                    b_id[1] 
                    b_id = (b_id[0], b_id[1], id[2] + b_ver_shift)
                except:
                    b_id = (b_id[0], id[1], id[2] + b_ver_shift)
            
            if not hasattr(cl,'energy_sigma0'):
                
                #check if job in queue
                if 'SLURM' in header.SCHEDULE_SYSTEM:
                    job_in_queue = cl.id[0]+'.'+cl.id[1] in runBash('ssh '+header.CLUSTER_ADDRESS+""" squeue -o '%o' """)
                    # print(job_in_queue, cl.name)
                
                else:
                    print_and_log('Error! I do not know how to check job status with chosen SCHEDULE_SYSTEM; Please teach me here! ')

                if not get_from_server(cl.dir+'/RUNNING', addr = header.CLUSTER_ADDRESS, trygz = False): #if exist than '' is returned
                    cl.state = '3. Running'
                



                elif job_in_queue:
                    cl.state = '3. In queue'
                    # print_and_log('Job is in queue')
                    # sys.exit()

                else:

                    if '2' in cl.state:
                        ''
                    else:
                        cl.state = '5. Some fault most probably'

                # sys.exit()

                print_and_log( cl.name, 'has state = ,',cl.state,'; I will continue; outcar file renamed to _unfinished')
                outcar = cl.path['output']
                outunf = outcar+"_unfinished"
                runBash("mv "+outcar+" "+outunf)

                continue

            e = calc[id].energy_sigma0
            n_m = calc[id].end.nznucl[0] # number of matrix atoms


            try:
                v = calc[id].end.vol
            except:
                v = 0
            #print e
            if e < emin: emin = e; id_min = id
            conv[n].append(id)
            # print base
            conv[base].append(b_id)


            outst2 = ("%s"%calc[id].name).ljust(name_field_length)
            outst2+='|'
            outst_end = '' 
            
            if   b_id :

                # if "4" not in calc[b_id].state:
                if readfiles:    
                    calc[b_id].read_results(loadflag, choose_outcar = choose_outcar)


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
                        
                        e_b = calc[b_id].energy_sigma0
                        n_m_b = calc[b_id].end.nznucl[0]
                        v_b = calc[b_id].end.vol
                        diffE = e - e_b/n_m_b*n_m
                        
                        outst2 += " {:.3f} & {:.2f} &".format( (diffE - energy_ref), (v - v_b) ).center(6)
                        result_list = [diffE - energy_ref, v - v_b]


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
            print_and_log( final_outstring, end = '')

        emin = 0
        





        """Aditional analysis, plotting"""
        if '4' not in calc[id].state:
            print_and_log( "res_loop(): Calculation ",id, 'is unfinished; return')
            return
        final_list = () #if some part fill this list it will be returned instead of final_outstring
        

        cl = calc[id]
        if b_id: bcl = calc[b_id]

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
            
            #print (" At zero pressure: segregation energy is %.0f  meV; Seg. volume is %.1f A^3; excess seg. vol. is %.2f A" %(e_seg, v_seg, v_seg/A ) )
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



                final_list = [id2[0]+'.'+id2[1], calc[id1].e_seg, calc[id1].v_seg, dgb, nneigbours, e_seg2, e_ch2, e_m2, segimp]
            
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

        elif analys_type == 'fit_a':
            """Fit equation of state for bulk systems.

            The following equation is used::

               sjeos (default)
                   A third order inverse polynomial fit 10.1103/PhysRevB.67.026103

                                   2      3        -1/3
               E(V) = c + c t + c t  + c t ,  t = V
                       0   1     2      3

               taylor
                   A third order Taylor series expansion about the minimum volume

               murnaghan
                   PRB 28, 5480 (1983)

               birch
                   Intermetallic compounds: Principles and Practice,
                   Vol I: Principles. pages 195-210

               birchmurnaghan
                   PRB 70, 224107

               pouriertarantola
                   PRB 70, 224107

               vinet
                   PRB 70, 224107

               antonschmidt
                   Intermetallics 11, 23-32 (2003)

               p3
                   A third order polynomial fit

                Use::

                   eos = EquationOfState(volumes, energies, eos='sjeos')
                   v0, e0, B = eos.fit()
                   eos.plot()

            """
            # e, v, emin, vmin       = plot_conv( conv[n], calc,  "fit_gb_volume2")
            alist = []
            vlist = []
            etotlist  = []
            magn1 = []
            magn2 = []
            for id in conv[n]:
                cl = calc[id]
                alist.append(cl.end.rprimd[0][0])
                etotlist.append(cl.energy_sigma0)
                vlist.append(cl.end.vol)
                magn1.append(cl.magn1)
                magn2.append(cl.magn2)
            eos = EquationOfState(vlist, etotlist, eos = 'sjeos')
            v0, e0, B = eos.fit()
            #print "c = ", clist[2]
            print_and_log( '''
            v0 = {0} A^3
            a0 = {1} A
            E0 = {2} eV
            B  = {3} eV/A^3'''.format(v0, v0**(1./3), e0, B)  )
            savedpath = 'figs/'+cl.name+'.eps'
            # plt.close()
            # plt.clf()
            # plt.close('all')
            
            eos.plot(savedpath, show = True)
            # plt.clf()

            if push2archive:
                push_figure_to_archive(local_figure_path = savedpath, caption = description_for_archive)







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
                print_and_log("Calculation ",bcl.id, 'is unfinished; return')
                return [[], []]


            #normalize numbers of atoms by some element except Li and Na 
            iLi = None; jLi = None
            # print cl.end.znucl
            for i, z in enumerate(cl.end.znucl):
                # print i, z
                if z in [3, 11]: 
                    iLi = i
                    # print 'iLi is found'
                    continue
                # print i, z

                for j, zb in enumerate(bcl.end.znucl):
                    if zb in [3, 11]: 
                        jLi = j
                        continue

                    if z == zb:
                        # print "I use ", z, " to normalize"
                        i_n = i
                        j_n = j


            # print "i, j",i, j
            # print 'nznucl cl',  cl.end.nznucl
            # print 'znucl cl',  cl.end.znucl
            n  = cl.end.nznucl[i_n]
            bn = bcl.end.nznucl[j_n]
            if iLi != None:
                nLi  = cl.end.nznucl[iLi]
            else:
                raise RuntimeError

            if jLi != None:
                bnLi  = bcl.end.nznucl[jLi]
            else:
                bnLi  = 0

            # print n, bn, nLi

            # print nLi/n

            mul = 1. / (float(nLi) / n)             

            # print mul


            redox = -(  ( cl.energy_sigma0 / n - bcl.energy_sigma0 / bn ) * mul  -  energy_ref  )



            final_outstring = ("{:} | {:.2f} eV \n".format(id[0]+'.'+id[1], redox  ))
            
            print_and_log( final_outstring )

            final_list = {'is':id[0], 'redox_pot':redox, 'id_is':id, 'id_ds':b_id, 
            'kspacing':cl.set.vasp_params['KSPACING'], 'time':cl.time/3600.,
            'mdstep':cl.mdstep, 'ecut':cl.set.vasp_params['ENCUT'], 'niter':cl.iterat/cl.mdstep,
            'set_is':id[1] }





        if cl.calc_method and ('neb' in cl.calc_method or 'only_neb' in cl.calc_method):
            path2mep_s = cl.dir+'/mep.eps'
            itise = cl.id[0]+'.'+cl.id[1]
            name_without_ext = 'mep.'+itise
            path2mep_l = cl.dir+name_without_ext+'.eps'
            if not os.path.exists(path2mep_l) or '2' in up:
                ''
                get_from_server(files = path2mep_s, to = path2mep_l, addr = cluster_address, )
                get_from_server(files = cl.dir+'/movie.xyz', to = cl.dir+'/movie.xyz', addr = cluster_address, )
            



            # trying to get one image closest to the saddle point
            if old_behaviour and cl.version == 2: #old behaviour, now created automatically in add callc
                im = cl.set.vasp_params['IMAGES']
                # if im % 2 > 0: #odd
                #     i = im//2 + 1
                # else:
                #     i = im/2
                # if choose_image:
                #     i = choose_image

                for i in range(im):
                    i+=1
                    cl_i = copy.deepcopy(cl)
                    cl_i.version+=i
                    cl_i.id = (cl.id[0], cl.id[1], cl_i.version)
                    cl_i.name = str(cl_i.id[0])+'.'+str(cl_i.id[1])+'.'+str(cl_i.id[2])
                    # print cl_i.name
                    cl_i.path["output"] = cl_i.dir+'0'+str(i)+"/OUTCAR"
                    # for i in range():

                    cl_i.associated_outcars = [ aso[2:] for aso in cl_i.associated_outcars  ]

                    # print cl_i.path["output"] 
                    cl_i.state = '2. Ready to read outcar'
                    # if not os.path.exists(cl_i.path["output"]):
                    #     load = 'o'
                    outst2 = ("%s"%cl_i.name).ljust(name_field_length)
                    if readfiles:
                        print(outst2+'|'+cl_i.read_results(loadflag, show = show, choose_outcar = choose_outcar) )
                    else:
                        print_and_log(outst2+' | File was not read')
                    

                    if cl_i.id in calc: #move creation of calcs with images to add_neb
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
                        calc[cl_i.id] = cl_i






            # print path2mep_l
            if 0:
                if os.path.exists(path2mep_l):
                    # get_from_server(file = path2mep_s, to = path2mep_l, addr = cluster_address)

                    runBash('evince '+path2mep_l)
                else:
                    a =  glob.glob(cl.dir+'*mep*')
                    if a:
                        runBash('evince '+a[0])


            
            #find moving atom
            cl1 = calc[cl.id[0], cl.id[1], 1]
            cl2 = calc[cl.id[0], cl.id[1], 2]
            
            diffv = np.array(cl1.init.xcart) - np.array(cl2.init.xcart)
            diffn = np.linalg.norm(diffv, axis = 1)
            atom_num = np.argmax(diffn) # number of atom moving along the path

            #prepare lists
            ni = cl.set.vasp_params['IMAGES']
            vlist = [1]+list(range(3, ni+3) )+[2]
            # print vlist
            mep_energies = []
            atom_pos     = []
            for v in vlist:
                cli = calc[cl.id[0], cl.id[1], v]
                # print cli.id
                # cli.end = return_to_cell(cli.end)
                # mep_energies.append(  min(cli.list_e_sigma0)   ) #use minimum energy - not very good, sometimes unconverged energy could be lower! 
                mep_energies.append(  cli.energy_sigma0   ) #use last energy 
                atom_pos.append( cli.end.xcart[atom_num] )

            # print np.array(atom_pos)

            #test if the distances between points are not spoiled by PBC 
            nbc = range(-1, 2)
            jj=0
            for x in atom_pos:

                x2 = atom_pos[jj+1]
                r = cl.end.rprimd
                d1, _ = image_distance(x, x2, r, order = 1) #minimal distance
                x2_gen = (x2 + (r[0] * i  +  r[1] * j  +  r[2] * k) for i in nbc for j in nbc for k in nbc) #generator over PBC images
                x2c = copy.deepcopy(x2)
                ii = 0
                while  np.linalg.norm(x - x2c) > d1: #find the closest PBC image position
                    if ii > 100:
                        break
                    ii+=1
                    x2c = next(x2_gen)
                atom_pos[jj+1] = x2c
                jj+=1
                if jj == len(atom_pos)-1: # the last point is not needed, we could not use slice since we need to use changed atom_pos in place
                    break
                # print np.linalg.norm(x - x2c), d1




            
            if 'mep' in show:
                plot_mep(atom_pos, mep_energies)
                plot_mep(atom_pos, mep_energies, image_name = 'figs/'+name_without_ext+'_my.png')


            if push2archive:
                path2saved = plot_mep(atom_pos, mep_energies, image_name = 'figs/'+name_without_ext+'_my')
                push_figure_to_archive(local_figure_path = path2saved, caption = description_for_archive)








                if 0: #copy files according to chosen outcar to run nebresults locally 
                    wd = cl_i.dir
                    out_i = cl_i.associated_outcars[choose_outcar-1]
                    out_1 = calc[cl.id[0],cl.id[1], 1].associated_outcars[choose_outcar-1]
                    out_2 = calc[cl.id[0],cl.id[1], 2].associated_outcars[choose_outcar-1]
                    # print out_1
                    # print out_2 
                    shutil.copyfile(wd+out_1, wd+'00/OUTCAR')
                    shutil.copyfile(wd+out_2, wd+'04/OUTCAR')
                    for d in ['01/','02/','03/' ]:
                        shutil.copyfile(wd+d+out_i, wd+d+'OUTCAR')

                        # print wd+d+out_i





    if final_list:
        return final_list, result_list
    else:
        return final_outstring.split('&'), result_list # only for last version or fit depending on type of analysis














def for_phonopy(new_id, from_id = None, calctype = 'read', mp = [10, 10, 10], additional = None):
    #creates file for phonopy, run phonopy
    #new_id - will add this calculation or read; if string then interpreted as filename of thermal_properties.yaml
    #from_id - tuple - than will take coordinates from the end; or path to input poscar file
    #type - 'create', 'read'
    #additional - list of calculation names, if calculation was splited into several parts

    mpstr = " ".join(map(str, mp))

    def read_phonopy_data(filename, key = "free_energy" ):
        # with open(filename, 'r') as f:
        #     f.readline()
        F = []
        T = []
        #     for line in f:
        #         T.append(float(line.split()[0]))
        #         F.append(float(line.split()[1]))
        #     # print T, F
        import yaml
        f = open(filename)
        # use safe_load instead load
        dataMap = yaml.safe_load(f)
        f.close()
        prop = dataMap['thermal_properties']
        for i in range(len(prop)):
            T.append( prop[i]['temperature'] )
            F.append( prop[i][key]     )

        coeffs1 = np.polyfit(T, F, 8)
        fit_func = np.poly1d(coeffs1)
        T_range = np.linspace(min(T), max(T))
        print_and_log( 'I return', key)

        return T_range, fit_func




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
        with open(confname, 'w') as f:
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
            with open(confname, 'w') as f:
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







"""Take structures from Mat. projects"""


def get_structure_from_matproj(struct_des, it, it_folder, ver, mat_proj_id = None):
    """
    Find material with 'it' stoichiometry (lowest energy) from materialsproject.org, 
    download and create field in struct_des and input POSCAR file
    ###INPUT:
        - struct_des-  
        - it        - materials name, such as 'LiCoO2', .... By default the structure with minimum *e_above_hull* is taken
        - it_folder - section folder in which the Poscar will be placed
        - ver       - version of structure defined by user
        - mat_proj_id (str) - the id can be provided explicitly
    
    ###RETURN:
        - ?
        - ?


    """
    with MPRester(pmgkey) as m:
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


    add_des(struct_des, it, it_folder, des = 'taken automatically from materialsproject.org: '+groundstate_st_id,)
    path2poscar = it_folder+'/'+it+'/'+groundstate_st_id+".POSCAR-"+str(ver)
    makedir(path2poscar)
    Poscar(st_pmg).write_file(path2poscar, direct=True, vasp4_compatible=True, )
    print_and_log("File "+path2poscar+" was written\n")
    
    return groundstate_st_id, path2poscar


# with MPRester(pmgkey) as m:
#     print dir(m)
#     print m.supported_properties
#     print m.get_data('mp-540111', data_type='vasp', prop='total_magnetization')
#     # 'total_magnetization'


def manually_remove_from_struct_des(struct_des, key):
    """
    
    """
    del struct_des[key]
    print_and_log('Attention! Entry '+key+' was removed from struct_des dict/\n')

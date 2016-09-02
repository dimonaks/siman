from header import *
import header
from classes import CalculationVasp, Description
from functions import gb_energy_volume
# from ext_databases import get_structure_from_matproj

from functions import makedir

from pymatgen.matproj.rest import MPRester
from pymatgen.io.vasp.inputs import Poscar

from operator import itemgetter



pmgkey = header.project_conf.pmgkey



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
            print_and_log('Please provide schedule_system!')
            raise RuntimeError

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
        print 'run sent'
        clean_run(header.project_conf.SCHEDULE_SYSTEM)
    
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
        print_and_log("Error! "+it+' already exist in struct_des; use override = True')
        raise RuntimeError
    else:
        struct_des[it] = Description(it_folder, des)
        # hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
        hstring = 'add_des(struct_des, {:s}, {:s}, {:s})    #on {:})'.format(it, it_folder, des, datetime.date.today())

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









def add_loop(it, ise, verlist, calc = None, conv = None, varset = None, 
    up = 'test', typconv="", from_geoise = '', inherit_option = None, 
    coord = 'direct', savefile = None, show = [], comment = None, 
    input_geo_format = 'abinit', ifolder = None, input_geo_file = None, ppn = None,
    calc_method = None, u_ramping_region = None, it_folder = None ):
    """
    Main subroutine for creation of calculations, saving them to database and sending to server.

    Input:
        - it, ise - structure name and set name of new calculation
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

        - coord - type of cooridnates written in POSCAR:
            'direct'
            'cart'

        - savefile - controls which files are saved during VASP run on server; check
            'ocdawx' - outcar, chgcar, dos, AECCAR,WAVECAR, xml

        - ifolder - explicit path to folder where to search for geo file.

        - input_geo_file - explicit file name of input file

        - it_folder - section folder (sfolder) used in struct_des; here needed with input_geo_format = mat_proj

        - show - allows to choose type of formatting. See read_results ?.

        - comment - arbitrary comment for history.

        - inherit_option - ? to be described.

        - typconv - ? to be described.
        
        - 'from_geoise' - part of folder name with geometry input files. allows to use geometry files from different sets.
        please find how it was used

        - ppn - number of cores used for calculation; overwrites header.corenum

        - calc_method - provides additional functionality:
            'u_ramping' - realizes U ramping approach #Phys Rev B 82, 195128

        - u_ramping_region - used with 'u_ramping'=tuple(u_start, u_end, u_step)


    Comments:
        To create folders and add calculations add_flag should have value 'add' 


    """
    it = it.strip()
    if it_folder: it_folder = it_folder.strip()

    schedule_system = header.project_conf.SCHEDULE_SYSTEM

    if ppn:
        header.corenum = ppn

    struct_des = header.struct_des

    if not calc:
        calc = header.calc
        conv = header.conv
        varset = header.varset

    header.close_run = True

    if type(verlist) == int: #not in [tuple, list]:
        #print "verlist"
        verlist = [verlist]; #transform to list




    hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
    args = hstring.split('(')[1].split(',')
    # arg1 = args[0], arg2 = args[1]
    # print args[0], it
    hstring = hstring.replace(args[0], "'"+it+"'")
    hstring = hstring.replace(args[1], "'"+str(ise)+"'")
    # print hstring
    try:
        if hstring != header.history[-1]: header.history.append( hstring  )
    except:
        header.history.append( hstring  )






    if ifolder:
        if it not in ifolder:
            print_and_log('Check ifolder !!! it not in ifolder')
            raise RuntimeError






    fv = verlist[0]; lv = verlist[-1];
    

    if typconv == "": setlist = list(ise)#
    elif up == "no_base": setlist = varset[ise].conv[typconv][1:]; up = "up1"
    else: setlist = varset[ise].conv[typconv] #
    setlist = [s.strip() for s in setlist]
    # print setlist

    nc = it+'.'+ise[0]+typconv
    try: conv[nc]; 
    except KeyError: 
        if typconv: conv[nc] = []    
    






    mat_proj_st_id = None
    if input_geo_format == 'mat_proj':
        print_and_log("Taking structure "+it+" from materialsproject.org ...\n")
        if it_folder == None:
            raise RuntimeError



        mat_proj_st_id, input_geo_file = get_structure_from_matproj(struct_des, it, it_folder, fv)
        input_geo_format = 'vasp'










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
            
            if typconv and id not in conv[nc]: conv[nc].append(id)
            
            try: blockfolder = varset[inputset].blockfolder
            except AttributeError: blockfolder = varset[ise].blockfolder
            
            cfolder = struct_des[it].sfolder+"/"+blockfolder #calculation folder
            
            add_calculation(it,inputset,v, fv, lv, input_folder, cfolder, calc, varset, up, 
                inherit_option, prevcalcver, coord, savefile, input_geo_format, input_geo_file, 
                schedule_system = schedule_system, 
                calc_method = calc_method, u_ramping_region = u_ramping_region,
                mat_proj_st_id = mat_proj_st_id)
            
            prevcalcver = v

    if up not in ('up1','up2'): 
        print_and_log("Warning! You are in the test mode, please change up to up1; "); 
        raise RuntimeError
    return
















def add_calculation(structure_name, inputset, version, first_version, last_version, input_folder, blockdir, calc, varset, update = "no",
    inherit_option = None, prevcalcver = None, coord = 'direct', savefile = None, input_geo_format = 'abinit', 
    input_geo_file = None, schedule_system = None, calc_method = None, u_ramping_region = None,
    mat_proj_st_id = None):
    """

    schedule_system - type of job scheduling system:'PBS', 'SGE', 'SLURM'

    prevcalcver - version of previous calculation in verlist

    """
    struct_des = header.struct_des





    id = (structure_name,inputset,version)
    try:
        calc[id]
        status = "exist"
        print_and_log(str(calc[id].name)+" has been already created and has state: "+str(calc[id].state)+"\n\n")
        if "4" in calc[id].state: 
            complete_state = calc[id].state
            status = "compl"
            if update == 'up2': 
                print 'Calculation', calc[id].name, 'is finished, continue'
                return
            pass
            if update != "up1": return #completed calculations updated only for "up1"

    except KeyError:
        #update = "up"
        status = "new"
        print_and_log( "There is no calculation with id "+ str(id)+". I create new with set "+str(inputset)+"\n" )        
    if "up" in update:
        if status in ["exist","compl"]: print_and_log("You asked to update existing calculation with id "+ str(id)+" Warning! I update creating new class \n" )         
        # print 'State', calc[id].state
        calc[id] = CalculationVasp( varset[id[1]] )
        # print 'State', calc[id].state

        calc[id].id = id 
        calc[id].name = str(id[0])+'.'+str(id[1])+'.'+str(id[2])
        calc[id].dir = blockdir+"/"+ str(id[0]) +'.'+ str(id[1])+'/'
        

        # all additional properties:
        calc[id].calc_method = calc_method
        



        calc[id].u_ramping_region = u_ramping_region
        if hasattr(calc[id].set, 'u_ramping_region'):
            print_and_log("Attention! U ramping method is detected from set\n\n")
            calc[id].calc_method = 'u_ramping'
            calc[id].u_ramping_region =calc[id].set.u_ramping_region


        # print dir(calc[id].set)

        if hasattr(calc[id].set, 'afm_ordering'):
            print_and_log("Attention! afm_ordering method is detected from set\n\n")
            calc[id].calc_method = 'afm_ordering'




        calc[id].cluster_address = cluster_address
        calc[id].project_path_cluster = project_path_cluster
        if mat_proj_st_id:
            calc[id].mat_proj_st_id = mat_proj_st_id







        if update in ['up1', 'up2']:
            if not os.path.exists(calc[id].dir):
                log.write( runBash("mkdir -p "+calc[id].dir) )         #Create directory if does not exist
                log.write( runBash("ssh "+cluster_address+" ' mkdir -p "+calc[id].dir+" ' ") )
            if id[2] == first_version:
                calc[id].write_sge_script(schedule_system = schedule_system) #write header only once

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

            print 'searchinputtemplate = ', searchinputtemplate

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


        for input_geofile in geofilelist:
            #print runBash("grep version "+str(input_geofile) )
            input_geofile = os.path.normpath(input_geofile)
            if input_geo_format in ['abinit',]:
                curv = int( runBash("grep version "+str(input_geofile) ).split()[1] )

            elif input_geo_format == 'vasp': 
                curv = int(input_geofile.split('-')[-1] ) #!Applied only for phonopy POSCAR-n naming convention

            elif input_geo_format == 'cif': 
                curv = int(os.path.basename(input_geofile).split('.')[0] )


            if curv == id[2]:

                copy_to = os.path.normpath(calc[id].dir+input_geofile.split('/')[-1])
                
                if input_geofile != copy_to:
                    # print os.path.normpath(input_geofile)
                    # print os.path.normpath(copy_to)

                    dire = os.path.dirname(copy_to)
                    if not os.path.exists(dire):
                        os.makedirs(dire)

                    shutil.copyfile(input_geofile, copy_to)

                
                if input_geo_format == 'abinit':
                    calc[id].read_geometry(input_geofile)
                
                elif input_geo_format == 'vasp':
                    calc[id].read_poscar(input_geofile)

                elif input_geo_format == 'cif':
                    if header.project_conf.CIF2CELL:
                        print runBash("cif2cell "+input_geofile+"  -p vasp -o "+input_geofile.replace('.cif', '.POSCAR'))
                        input_geofile = input_geofile.replace('.cif', '.POSCAR')
                        
                        #check
                        if not os.path.exists(input_geofile):
                            print_and_log("Error! Something wrong with conversion of cif2cell: \n")
                            raise RuntimeError


                        calc[id].read_poscar(input_geofile)
                    


                    else:
                        print_and_log("Error! cif2cell is not available in your system")
                        raise RuntimeError

                
                else:
                    raise RuntimeError

                

                break
        
        if calc[id].path["input_geo"] == None: 
            print_and_log("Error! Could not find geofile in this list: "+ str(geofilelist)+  "\n")
            raise NameError #
        
        #print  calc[id].des
        calc[id].des += ' '+struct_des[id[0]].des + '; ' + varset[id[1]].des
        #print  calc[id].des

        calc[id].check_kpoints()    
        calc[id].actualize_set()





        if update in ['up1', 'up2']:
            calc[id].write_structure(str(id[2])+".POSCAR", coord, inherit_option, prevcalcver)
            calc[id].write_sge_script(str(version)+".POSCAR", version, inherit_option, prevcalcver, savefile, schedule_system = schedule_system)
            #if calc[id].path["output"] == None:
            calc[id].path["output"] = calc[id].dir+str(version)+".OUTCAR" #set path to output

            if id[2] == first_version:
                calc[id].add_potcar()
            
            calc[id].calculate_nbands()

            if id[2] == last_version:        
                calc[id].make_incar_and_copy_all(update)
                calc[id].make_run(schedule_system = schedule_system)


        # if status == "compl": calc[id].state = complete_state #Even if completed state was updated, the state does not change
        if status == "compl": calc[id].state = '3. Can be completed but was reinitialized' #new behavior 30.08.2016


        print_and_log("\nCalculation "+str(id)+" added or updated\n\n")

    return













def inherit_icalc(inherit_type, it_new, ver_new, id_base, calc = None,
    id_from = None,
    atom_new = None, atom_to_replace = None,  used_cell = 'end', atom_to_remove = None, id_from_used_cell = 'end',
    atom_to_shift = None, shift_vector = None):
    """
    Function for creating new geo files in geo folder based on different types of inheritance
    Input args: 
        it_new, ver_new - name of new structure,
        id_base - new structure will be based on the final structure of this calculation;     (can be either Calculation() object or path to geo file)
        id_from - can be additionally used to adopt for example rprimd from id_from to it_new; (can be either Calculation() object or path to geo file)

        
        inherit_type = '':
            full       - full inheritance of final state
            r2r3       - use r2 and r3 from id_from
            r1r2r3     - use r1, r2 and r3 from id_from
            remove_imp - removes all atoms with typat > 1
            replace_atoms - atoms of type 'atom_to_replace' in 'id_base' will be replaced by 'atom_new' type.

            make_vacancy - produce vacancy by removing 'atom_to_remove' starting from 0

        used_cell - use init or end structure of id_base calculation. Now realized only for 'replace_atoms' regime of inheritance


        atom_to_shift - number of atom to be shifted; starting from 1.
        shift_vector - vector in decart cooridinates (Angstrom!!) by which the atom will be shifted


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


    """


    hstring = ("%s    #on %s"% (traceback.extract_stack(None, 2)[0][3],   datetime.date.today() ) )
    if hstring != header.history[-1]: header.history.append( hstring  )

    #if inherit_type not in header.history[-1] or \
    #it_new not in header.history[-1]:   header.history.append( hstring  )
    calc = header.calc
    struct_des = header.struct_des

    if type(id_base) == str:
        print_and_log('Reading id_base\n')
        cl_base = CalculationVasp()
        cl_base.read_geometry(id_base)
        cl_base.id = ('from_file', 'from_file', cl_base.version)
        cl_base.name = id_base
        cl_base.end = cl_base.init

    else:
        cl_base = calc[id_base]


    if id_from:
        if type(id_from) == str: # if string - treated like file name
            calc_from = CalculationVasp();
            calc_from.read_geometry(id_from)
            calc_from.end = calc_from.init
            calc_from_name = id_from
        else:
            print "id_from", id_from
            calc_from = calc[id_from]
            calc_from_name = calc_from.name

        if id_from_used_cell == 'end':
            st_from = calc_from.end
        elif id_from_used_cell == 'init':
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


    new.path["input_geo"] = geo_folder + struct_des[it_new].sfolder + '/' + \
        it_new+"/"+it_new+'.inherit.'+inherit_type+'.'+str(ver_new)+'.'+'geo'
    print new.path["input_geo"]
    new.version = ver_new

    if inherit_type == "r2r3":
        des = ' Partly inherited from the final state of '+cl_base.name+'; r2 and r3 from '+calc_from_name
        new.des = struct_des[it_new].des + des
        new.end.rprimd[1] = st_from.rprimd[1].copy()
        new.end.rprimd[2] = st_from.rprimd[2].copy()       
        new.write_geometry("end",des)

    elif inherit_type == "r1r2r3":
        des = ' Partly inherited from the final state of '+cl_base.name+'; r1, r2, r3 from '+calc_from_name
        new.des = struct_des[it_new].des + des
        new.end.rprimd = copy.deepcopy( st_from.rprimd )
        new.hex_a = calc_from.hex_a
        new.hex_c = calc_from.hex_c
        new.end.xcart = xred2xcart(new.end.xred, new.end.rprimd) #calculate new xcart from xred, because rprimd was changed
        new.write_geometry("end",des)


    elif inherit_type == "full":
        print_and_log("Warning! final xred and xcart was used from OUTCAR and have low precision. Please use CONTCAR file \n");
        des = 'Fully inherited from the final state of '+cl_base.name
        new.des = des + struct_des[it_new].des
        new.write_geometry("end",des)

    elif inherit_type == "atom_shift":
        des = 'obtainded from final state of '+cl_base.name+' by shifting atom '+ str(atom_to_shift) +' by '+ str(shift_vector)
        new.des = des + struct_des[it_new].des
        
        new.end.xcart[atom_to_shift-1] += np.asarray(shift_vector) 
        new.end.xred = xcart2xred(new.end.xcart, new.end.rprimd)
        new.write_geometry("end",des)
        write_xyz(new.end)



    elif inherit_type == "remove_imp":
        """
        Assumed that typat == 1 is matrix atoms
        """
        des = 'All impurities removed from the final state of '+cl_base.name
        new.des = des + struct_des[it_new].des
     
        st = cl_base.end
        new.end.typat = []
        new.end.xred = []
        new.end.xcart = []
        new.end.ntypat = 1
        new.end.znucl = new.end.znucl[0:1]
        for i, t in enumerate(st.typat):
            if t == 1:
                new.end.typat.append(t)
                new.end.xred.append(st.xred[i])
                new.end.xcart.append(st.xcart[i])
        new.end.natom = len(new.end.xred)
        new.init = new.end #just for sure
        new.write_geometry("end",des)        
        write_xyz(new.end)

    elif inherit_type == "make_vacancy":
        """Remove  atom 'atom_to_remove' from final state of id_base"""
        des = 'Atom '+str(atom_to_remove)+' removed from  '+cl_base.name
        new.des = des + struct_des[it_new].des

        st = new.end
        del st.typat[atom_to_remove]
        del st.xcart[atom_to_remove]
        del st.xred[atom_to_remove]
        ntypat = len(set(st.typat))
        znucl_new = []
        for t in sorted(set(st.typat)):
            znucl_new.append( st.znucl[t-1]  )
        st.znucl = znucl_new
        st.natom -=1
        new.write_geometry("end", des)      

        #make visualization by adding to vacancy hydrogen atom
        st = copy.deepcopy(cl_base.end)
        st.typat[atom_to_remove] = max(st.typat)+1
        st.ntypat+=1
        st.znucl.append(1)
        write_xyz(st, file_name = "test_of_vacancy_creation."+str(new.version)+"."+st.name)






    elif inherit_type == "replace_atoms":
        """Simply replace in calculation one atoms by another """
        z_new     = element_name_inv(atom_new)
        z_replace = element_name_inv(atom_to_replace)


        if used_cell == 'init':
            st = new.init
        elif used_cell == 'end':
            st = new.end


        znucl = st.znucl
        
        if z_replace not in znucl: 
            print "Error! Calc "+new.name+" does not have atoms of this type\n"
            raise RuntimeError

        if atom_to_replace not in id_base[0] or atom_new not in it_new:
            print "Error! Something wrong with names of atom types\n"
            raise RuntimeError            
        print "replace ", z_replace, "by", z_new

        znucl = [int(z) for z in znucl] # convert to int
        i_r = znucl.index(z_replace)
        print "index ", i_r
        znucl[i_r] = z_new

        if znucl[-2] == znucl[-1]: #just for special case
            del znucl[-1]

            for i, t in enumerate(st.typat):
                if t == st.ntypat:
                    print "found ntypat" ,t
                    st.typat[i]-=1
                    #print t
            st.ntypat-=1
            #print st.typat
            #print new.end.typat
        
        st.znucl = znucl



        des = 'Fully inherited from the '+ used_cell +' state of '+cl_base.name+\
        ' by simple replacing of '+atom_to_replace+' by '+atom_new

        new.des = des + struct_des[it_new].des
        new.write_geometry(used_cell, des)
        # write_xyz(st)

    else:
        print_and_log("Error! Unknown type of Calculation inheritance"); raise RuntimeError
    return




























def res_loop(it, setlist, verlist,  calc = None, conv = {}, varset = {}, analys_type = 'no', b_id = (), 
    typconv='', up = "", imp1 = None, imp2 = None, matr = None, voronoi = False, r_id = None, readfiles = True, plot = True, show = [], 
    comment = None, input_geo_format = None, savefile = None, energy_ref = 0, ifolder = None, bulk_mul = 1, inherit_option = None,
    calc_method = None, u_ramping_region = None, input_geo_file = None,
    it_folder = None):
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
        
        readfiles - define if files to be read at all or only additional analysis is required; 



        The next three used for 'clusters' regime:    
        imp1 - key of bulk cell with one imp1
        imp2 - key of bulk cell with one imp2
        matr - key of bulk cell with pure matrix.


        show - list, allows to show additional information (force)

        energy_ref - energy in eV; substracted from energy diffs
        
        bulk_mul - allows to scale energy and volume of bulk cell during calculation of segregation energies

    RETURN:
        result_list - list of results

    """
















    if type(show) == str:
        show = [show]
    if type(verlist) == int: #not in [tuple, list]:
        #print "verlist"
        verlist = [verlist]; 
    if type(setlist) == str: #not in [tuple, list,]:
        setlist = [setlist]
    if not calc:
        calc = header.calc

    try:
        b_ver_shift = b_id[2] #add to version of base with respect to version of main
    except:
        b_ver_shift = 0


    """Copying files from server"""
    #to make
    # for inputset in setlist:
    #     for v in verlist:
    #         cl = calc[(it,inputset,v)]
    #         if not exists: 





    #print verlist
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
        if '4' not in calc[b_id].state:
            calc[b_id].read_results('o')
        e_b = calc[b_id].energy_sigma0
        v_b = calc[b_id].end.vol

    #define reference values
    e1_r = 0
    if type(r_id) in (float, int):
        e1_r = r_id
    elif type(r_id) == tuple:
        if '4' not in calc[r_id].state:
            print "start to read reference:"
            print calc[r_id].read_results('o')
            # print calc[r_id].read_results()
        e_r = calc[r_id].energy_sigma0 #reference calc
        nat_r = calc[r_id].end.natom # reference calc
        e1_r = e_r/nat_r # energy per one atom
        # print e1_r


    """Main loop"""
    final_outstring = 'no calculation found'
    for inputset in setlist:
        for v in verlist:
            # print 'Starting loops'

            id = (it,inputset,v)

            if id not in calc:
                print "Calculation does not exist!!!" 
                continue #pass non existing calculations

            cl = calc[id]

            # cl.set.update()
            

            if not hasattr(cl,'version'):
                calc[id].version = v



            path_to_outcar = cl.path["output"]
            # path_to_contcar = cl.dir+str(v)+".CONTCAR"
            # path_to_xml   = cl.dir+str(v)+".vasprun.xml"
            # print path_to_contcar, path_to_outcar
            # print os.path.exists(path_to_contcar)
            # print os.path.exists(cl.path["output"])

            outst = ' File was not read '
            
            if readfiles:
                load = up
                # print path_to_outcar
                flag1 = 1#os.path.exists(path_to_outcar)
                flag3 = 1; #ot load
                flag2 = 1
                # if 'x' in load:
                #     flag2 = os.path.exists(path_to_xml)
                # else:
                #     flag2 = True
                # print flag1
                if flag1 and flag2 and flag3:
                
                    outst = calc[id].read_results(load, analys_type, voronoi, show)
                
                else:
                    print "Trying to download OUTCAR and CONTCAR from server\n\n"
                    outst = calc[id].read_results('o', analys_type, voronoi, show)




            if analys_type in ('e_seg', 'coseg'):
                try:
                    b_id[1] 
                    b_id = (b_id[0], b_id[1], id[2] + b_ver_shift)
                except:
                    b_id = (b_id[0], id[1], id[2] + b_ver_shift)
            if not hasattr(cl,'energy_sigma0'):
                print cl.name, 'is not finished!, continue; file renamed to _unfinished'
                outcar = cl.dir+str(v)+".OUTCAR"
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


            outst2 = ("%s"%calc[id].name).ljust(22)
            outst2+='&'
            # print outst2+'&'            
            outst_end = '' 
            
            if   b_id :

                if "4" not in calc[b_id].state:    
                    calc[b_id].read_results(1)

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
                    print 'Calculating matrix_diff...'
                    
                    e_b = calc[b_id].energy_sigma0
                    n_m_b = calc[b_id].end.nznucl[0]
                    v_b = calc[b_id].end.vol
                    diffE = e - e_b/n_m_b*n_m
                    
                    outst2 += " {:.3f} & {:.2f} &".format( (diffE - energy_ref), (v - v_b) ).center(6)
                    result_list = [diffE - energy_ref, v - v_b]


                elif analys_type == 'diff': #
                    print 'Calculating diff...'
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
            print final_outstring

        emin = 0
        





        """Aditional analysis, plotting"""
        if '4' not in calc[id].state:
            print "Calculation ",id, 'is unfinished; return'
            return
        final_list = () #if some part fill this list it will be returned instead of final_outstring
        

        cl = calc[id]
        if b_id: bcl = calc[b_id]

        if analys_type == 'gbe':
            print("\nGrain boundary energy and excess volume fit:")
            plot_conv( conv[n], calc, "fit_gb_volume")

        elif analys_type == 'gbep':
            print("\nGrain boundary energy and excess volume fit:")
            # plot_conv( conv[n], calc, "fit_gb_volume")
            final_outstring = plot_conv( conv[n], calc, "fit_gb_volume_pressure")

        elif analys_type in ('e_seg', 'coseg') and len(verlist) > 3:

            #Test lateral sizes
            A   = calc[id].end.yzarea
            A_b = calc[b_id].end.yzarea
            
            if A != A_b: 
                print_and_log("Warning! you are trying to compare calcs with different lateral sizes: "+str(A)+" "+str(A_b))
                print "Areas are ", A, A_b," A^3"
            
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
                print "Warning! Calculation ", id1, 'does not have e_seg and v_seg. Try to run with readfiles = True to calculate it.'
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
            print id2, 'is id2'
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
                list = [ x  for x, t  in zip(st_rr.xcart, st_rr.typat) if np.linalg.norm(x_central - x) < dmax and t == 1]
                nneigbours =  len(list)


                final_outstring = ("%s.fit.pe & %.0f & %.1f & %.2f & %.d & %4.0f & %4.0f & %4.0f & %s " %(
                    id2[0]+'.'+id2[1], calc[id1].e_seg, calc[id1].v_seg, dgb, nneigbours, e_seg2, e_ch2, e_m2, segimp  ))

                # final_outstring = ("%s.fit.pe & %.0f & %.0f & %.1f & %.1f & %.2f & %.d & %4.0f & %4.0f & %4.0f & %s " %(
                #     id2[0]+'.'+id2[1], calc[id1].e_seg, e_segmin, calc[id1].v_seg, v_segmin , dgb, nneigbours, e_seg2, e_ch2, e_m2, segimp  )) #e_segmin and v_segmin are minimum energy (but at some pressure) and corresponing volume



                final_list = [id2[0]+'.'+id2[1], calc[id1].e_seg, calc[id1].v_seg, dgb, nneigbours, e_seg2, e_ch2, e_m2, segimp]
            
            elif analys_type == 'coseg' :
                calc[id2].e_seg = calc[id1].e_seg #save in version 2
                calc[id2].v_seg = calc[id1].v_seg
                final_outstring = ("%s.fit.pe & %.0f & %.1f & %.1f & %.1f" %(id[0]+'.'+id[1], calc[id2].e_seg, calc[id2].v_seg, d1, d2 ))




            print  final_outstring
            print '\\hline'


        elif analys_type == 'e_2imp':
            # plot_conv( conv[it], calc,  analys_type, conv_ext) #instead use plot_conv( conv['hs443OO'], calc,  'e_2imp', [conv['hs443CO'],conv['hs443CC']]) 
            pass



        elif analys_type == 'fit_ac':
            #for x in calc[id_min].end.xred:
            #    print x[0], x[1], x[2]
            #print ( outst2 + " min_e & "+calc[id_min].read_results(v,) )
            #print ("name %s_template          acell  %.4f  %.4f  %.4f # fit parameters are &%.4f &%.4f &%i &%i"  % (fit_hex(0.0002,0.0003,400,600, it, inputset, verlist, calc) )  )    
            print ("name %s_template          acell  %.5f  %.5f  %.5f # fit parameters are &%.5f &%.5f &%i &%i"  % (fit_hex(0.00002,0.00003,4000,6000, it, inputset, verlist, calc) )  )    
            #name, e_min, a_min, c_min,a,b,c,d = fit_hex(0.0002,0.0003,400,600, it, inputset, verlist, calc)
            #calc[(name,'93',1)].energy_sigma0 = e_min
            #calc[(name,'93',1)].state = "4"

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
            print '''
            v0 = {0} A^3
            a0 = {1} A
            E0 = {2} eV
            B  = {3} eV/A^3'''.format(v0, v0**(1./3), e0, B)
            eos.plot("images/fit_a.png", show = True)

            # my_plot(xlim = (2.79, 2.86),
            #     # etotl = (alist, etotlist, '-'),  
            #     # magn1 = (alist, magn1, '-'),
            #     magn2 = (alist, magn2, '-'),
            #     )



            pass


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



            final_outstring = ("{:} & {:.2f} eV \n".format(id[0]+'.'+id[1], redox  ))
            
            print final_outstring

            final_list = [id[0], redox, id, b_id, ]

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
        print 'I return', key

        return T_range, fit_func




    if calctype == 'create':
        work_path = 'geo/'+struct_des[new_id[0]].sfolder+'/'+new_id[0]
    
        log_history(  "{:}    #on {:}".format( traceback.extract_stack(None, 2)[0][3],   datetime.date.today() )  )

        print type(from_id)
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
        print runBash('export PYTHONPATH=~/installed/phonopy-1.9.5/lib/python:$PYTHONPATH; rm POSCAR-*; /home/dim/installed/phonopy-1.9.5/bin/phonopy '
            +confname+' -c '+posname+' -d --tolerance=0.01')

        ndis = len( glob.glob('POSCAR-*') )
        print ndis, ' displacement files was created'

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
            print range(1,npos+1), 'range'
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
                print ndis, ' displacement files was Found'
                print runBash('export PYTHONPATH=~/installed/phonopy-1.9.5/lib/python:$PYTHONPATH; /home/dim/installed/phonopy-1.9.5/bin/phonopy '
                    +'  -f {1..'+str(ndis)+'}.vasprun.xml --tolerance=0.01')

            #calculate thermal prop
            result = 'thermal_properties_'+mpstr.replace(" ", "_")+'.yaml'
            if not os.path.exists(result):

                posname = 'SPOSCAR'
                print runBash('export PYTHONPATH=~/installed/phonopy-1.9.5/lib/python:$PYTHONPATH; /home/dim/installed/phonopy-1.9.5/bin/phonopy '
                    +confname+' -c '+posname+' -t -p -s --tolerance=0.01')

                shutil.copyfile('thermal_properties.yaml', result)
    

            T_range, fit_func = read_phonopy_data(result)

            os.chdir(savedPath)

        if type(new_id) == str:
            result = new_id
            T_range, fit_func = read_phonopy_data(result)



    return T_range, fit_func







"""Take structures from Mat. projects"""


def get_structure_from_matproj(struct_des, it, it_folder, ver):
    """
    Find material with 'it' stoichiometry (lowest energy) from materialsproject.org, 
    download and create field in struct_des and input POSCAR file
    ###INPUT:
        - struct_des-  
        - it        - materials name, such as 'LiCoO2', ...
        - it_folder - section folder in which the Poscar will be placed
        - ver       - version of structure defined by user
    
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
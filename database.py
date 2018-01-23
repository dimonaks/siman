#Copyright Aksyonov D.A
from __future__ import division, unicode_literals, absolute_import 

import shelve, sys, datetime, shutil, tempfile, os, json, re, glob

import pandas as pd


import header
from header import runBash, print_and_log, printlog
from classes import CalculationVasp
from small_functions import makedir
from set_functions import init_default_sets
from functions import invert

"""
Module contains utilities for project database management
and working with external databases

TODO:
1) in write_database() make update of history file more clever
"""
calc_key       = 'calc'
conv_key       = 'conv'
varset_key     = 'varset'
history_key    = 'history'
struct_des_key = 'struct_des'


def read_database(scratch = False):
    """
    Read database of calculations

    INPUT:
        scratch - not used
    
    RETURN:
        calc   - dict, contains all calculations of the project
        conv   - dict, convergence sequences
        varset - dict, parameter sets of the project
        size_on_start, int - not used now
    """

    # databasefile = 'calc.s' #was used with python2
    # databasefile = 'calc.gdbm'
    databasefile3 = 'calc.gdbm3'
    # if header.RAMDISK:
    #     databasefile3 = header.RAMDISK+databasefile3




    # if scratch == True: databasefile =   '/scratch/aksenov/calc.s'
    
    printlog("\nLaunch at "+str( datetime.datetime.today() )+'\n')
    
    # mod = __import__("gdbm")
    # d = shelve.Shelf(mod.open(databasefile, protocol=1))

    # print(databasefile3)
    d = shelve.open(databasefile3, protocol = 3)


    try:
        # calc              = d[calc_key]; 
        # calc              = {}; 
        header.conv       = d[conv_key]; 
        header.varset     = d[varset_key]; 
        header.history    = d[history_key]
        header.struct_des = d[struct_des_key]

    except KeyError: 
        
        # try: calc = d[calc_key] #dictionary of calculations
        # except KeyError:
        #     printlog( "There is no database of calculations. I create new"); calc = {}
        


        try: header.conv = d[conv_key] #dictionary of convergence lists
        except KeyError:
            printlog( "There is no dictionary of convergence lists. I create new"); conv = {}   
        


        try: header.varset = d[varset_key] 
        except KeyError:
            printlog( "There is no dictionary of inputsets. I create new");  varset = {} 
        


        try: header.history = d[history_key] 
        except KeyError:
            header.history = ['Project started on '+ str( datetime.date.today() ) ]

            printlog( "There is still no history in database. The list is in header module ");
        


        try: header.struct_des = d[struct_des_key] 
        except KeyError:
            printlog( "There is no struct_des in database. The dict is in header module "); 

    d.close()
    #print history
    init_default_sets()


    return header.conv, header.varset, sys.getsizeof(d)


def write_database(calc = None, conv = None, varset = None, size_on_start = None):
    """
    The function writes main dictionaries to database file calc.s
    Also creates copy of calc.s

    INPUT:
        calc - dict, contains all calculations of the project
        conv - dict, convergence sequences
        varset - dict, parameter sets of the project
        size_on_start - not used now

    RETURN:
        None    

    """
    #size_on_finish = sys.getsizeof(dbase)
    #if size_on_finish != size_on_start:
    # runBash("cp calc.s calc_copy.s") #create copy before writing
    databasefile3 = 'calc.gdbm3'

    # if header.RAMDISK:
    #     databasefile3 = header.RAMDISK+databasefile3

    shutil.copyfile(databasefile3, 'calc_copy.gdbm3')
    if 0:
        d = shelve.open('calc.s', protocol=1) #Write database of calculations
        d[calc_key]       = calc
        d[conv_key]       = conv
        d[varset_key]     = varset
        d[history_key]    = header.history
        d[struct_des_key] = header.struct_des 
        d.close()

    python2 = False
    if python2:
        import gdbm# use in python2

        d = shelve.Shelf(gdbm.open('calc.gdbm', 'c'), protocol=1) #Write dbm database for python3
        d[unicode(calc_key)]       = calc
        d[unicode(conv_key)]       = conv
        d[unicode(varset_key)]     = varset
        d[unicode(history_key)]    = header.history
        d[unicode(struct_des_key)] = header.struct_des 
        d.close()
    else: #python3 
        import dbm
        # d = shelve.Shelf(dbm.open('calc.gdbm', 'c'), protocol=1) #Write dbm database for python3
        # d[calc_key]       = calc
        # d[conv_key]       = conv
        # d[varset_key]     = varset
        # d[history_key]    = header.history
        # d[struct_des_key] = header.struct_des 
        # d.close()        

        d = shelve.Shelf(dbm.open(databasefile3, 'n'), protocol = 3) #Write dbm database for python3
        # d[calc_key]       = calc
        d[conv_key]       = header.conv
        d[varset_key]     = header.varset
        d[history_key]    = header.history
        d[struct_des_key] = header.struct_des 
        d.close()   

        with shelve.Shelf(dbm.open(header.calc_database, 'w'), protocol = 3) as d:
            for key in header.calc:
                # print(key)
                d[str(key)] = header.calc[key]
        
        # with dbm.open(header.calc_database, 'w') as d:
        #     d.reorganize()


    printlog("\nEnd of work at "+str(datetime.datetime.now())+'\n')

    try:
        header.log.close()
    except:
        pass

    #Update history file
    with  open('history','w') as his:
        #print history
        for i in header.history:
            #print i
            his.write(i+"\n")
    
    print("\nDatabase has been successfully updated\n")
    
    return






def get_from_database(x1, x2, mat, inquiry_keys = None, silent = None, ssh_object = None):
    """
    inquiry_keys (list) - list of keys that should exist in filenames both for x1 and x2
    ssh_object (SSHTools) - ssh object based on paramiko with access details

    """

    def check(key, inquiry_keys):
        return all([k in key for k in inquiry_keys])


    path2database = '/home/Data/CEStorage/'

    hash_dict_file = 'hash_dict.json'

    cluster_path2hash = os.path.join(path2database, hash_dict_file)

    if inquiry_keys is None:
        inquiry_keys = []

    if ssh_object:
        # ssh_object.get()
        tempdir = tempfile.gettempdir()
        local_path2hash = os.path.join(tempdir, hash_dict_file)

        ssh_object.get(cluster_path2hash,  local_path2hash  )

        # sys.exit()

    with open(local_path2hash, 'r') as fp:
        hash_dict = json.load(fp)

    # print(hash_dict)
    x1s = []
    x2s = []
    # print(hash_dict)
    for key, val in hash_dict.items():
        if check(key, inquiry_keys+[x1, mat]):
            x1s.append(key)

        if check(key, inquiry_keys+[x2, mat]):
            x2s.append(key)

    x1s = sorted(x1s, key = lambda el: len(el) )
    x2s = sorted(x2s, key = lambda el: len(el) )


    for xi, xis in (x1, x1s), (x2, x2s):
        if not silent:
            print('\nFiles for',xi,':')
        for i, f in enumerate(xis):
            if not silent:
            
                print(i+1,f)


    if len(x1s) == 0 or len(x2s) == 0:
        print('No information in database for this inquire:', x1, x2, mat, str(inquiry_keys) )
        return None, None
    
    key1 = x1s[0]
    key2 = x2s[0]

    if not silent:

        print('\nI choose first entries for both concentrations:',key1, 'and', key2,'\n')
    # print('Use *inquiry_keys* arg to clarify the output results.\n')

    #get files
    loc1 = os.path.join(tempdir, hash_dict[key1])
    loc2 = os.path.join(tempdir, hash_dict[key2])
    makedir(loc1)
    makedir(loc2)
    # print()/

    ssh_object.get(os.path.join(path2database, hash_dict[key1]), loc1  )
    ssh_object.get(os.path.join(path2database, hash_dict[key2]), loc2  )



    cl1 = CalculationVasp().deserialize(loc1)
    cl2 = CalculationVasp().deserialize(loc2)

    return cl1, cl2




def push_figure_to_archive(local_figure_path, caption, figlabel = None, autocompl = True ):
    shutil.copy(local_figure_path, header.path_to_images)
    print_and_log('push_figure_to_archive():', local_figure_path, 'copied to', header.path_to_images, imp = 'y')
    
    name_without_ext =   '.'.join( os.path.basename(local_figure_path).split('.')[:-1]) 
    figfile = '{{'+name_without_ext+'}}'

    if not figlabel:
        figlabel = '.'.join(name_without_ext.split('.')[:-1])
    
    if autocompl:
        caption+=' for '+figlabel 

    tex_text = \
    ("\\begin{{figure}} \n\includegraphics[width=\columnwidth]{{{:s}}}\n"
    "\caption{{\label{{fig:{:s}}} {:s} }}\n"
    "\end{{figure}}\n").format(figfile, figlabel, caption )


    # print (tex_text)
    with open(header.project_conf.path_to_paper+'/auto_fig.tex', 'a+', newline = '') as f:
        f.seek(0)
        a = f.read()
        # print (a)
        if tex_text not in a:
            f.write(tex_text)
    return





def read_cvs_database(columns):
    """
    Allows to read cvs file with experimental results
    """

    dfs = pd.read_csv(r'database/literature.csv')[columns]
    dfs.drop(0, inplace=True)
    dfs.dropna(inplace=True)
    dfs.set_index('is', inplace = True)
    dfs = dfs.apply(lambda x: pd.to_numeric(x, errors='ignore'))
    # dfs['owner'] = 'world'    
    # print(dfs.sort_index())
    # print(dfs)
    # sys.exit()
    return dfs




def add_to_archive_database(cl, subgroup):
    """
    cl is Calculation which should be added to database
    subgroup (str) - subgroup folder

    """
    join = os.path.join

    save_format = 'azh'
    dbpath = header.PATH2DATABASE
    it = cl.id[0]
    material = header.struct_des[it].sfolder.split('/')[0]
    sfolder = os.path.join(dbpath, material)

    name = []

    if 'azh' in save_format:
        #1. Single point calculation of total energy
        # print(sfolder)
        makedir( join( sfolder, 'dummy')  )

        #determine x for alkali ion from structure name
        parsed = re.findall(r'([A-Z][a-z]*)(\d*)', it)
        parsed = [(el, x if x else '1') for (el, x) in parsed ]
        print('detected element is ', parsed[0][0])

        if parsed[0][0] in [invert(z) for z in header.ALKALI_ION_ELEMENTS]:
            x = parsed[0][0]

            if hasattr(cl, 'max_alk_ion_content'):
                x = float(x)/cl.max_alk_ion_content
            else:
                x = '1'

        else:
            x = '0'
        name.append('x'+x)

        cl.read_results()
        
        if material in ['LiCoO2', 'LiTiO2', 'LiFePO4', 'NaFePO4', 'LiMnPO4', 
        'LiNiO2', 'LiTiS2', 'LiMn2O4', 'LiVP2O7', 'LiVPO4F', 
        'NaMnAsO4', 'Na2FePO4F', 'Na2FeVF7', 'KFeSO4F', 'NaLiCoPO4F', 'KVPO4F' ]: 
            sfolder = join(sfolder, subgroup)
            makedir( join(sfolder,'dummy') )

        cl.set.update()




        (pot, func) = cl.potcar_lines[0][0].split('_')
        
        if cl.set.spin_polarized:
            func = 'U'+func #unrestricted


        if hasattr(cl.set, 'u_ramping_nstep') and cl.set.u_ramping_nstep:
            func += '-UR'
        elif cl.set.dftu:
            func += '-U'

        func+=pot.lower()
        ecut = str(round(cl.set.ecut ))

        func+=ecut
        # print(func)
        name.append(func)

        name.extend(it.split('.')[1:]+[cl.id[1]]+[str(cl.id[2])])

        name_str = '_'.join(name)
        # print('_'.join(name) )

        # sys.exit()



        outcar_name = name_str+'.out'

        shutil.copyfile(cl.path["output"], join(sfolder, outcar_name)  )

        cl.end.write_xyz(path = sfolder, filename =  name_str)

        pickle_file = cl.serialize(os.path.join(sfolder, 'dat', name_str) )
        # cl




        #write input, problem with fitted version 100, which does not have input geometry, since they are created on cluster
        # makedir(sfolder+'input/dummy')
        # shutil.copyfile(cl.path["input_geo"], sfolder+'input/'+name_str+'.geo')


        st_mp = cl.end.convert2pymatgen()
        sg_before =  st_mp.get_space_group_info() 
        # from pymatgen.symmetry.finder import SymmetryFinder
        # sf = SymmetryFinder(st_mp_prim)
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        symprec = 0.1
        sf = SpacegroupAnalyzer(st_mp, symprec = symprec)

        st_mp_prim = sf.find_primitive()
        # st_mp_prim = sf.get_primitive_standard_structure()
        # st_mp_prim = sf.get_conventional_standard_structure()


        # st_mp_conv = sf.get_conventional_standard_structure()
        # print(st_mp_conv)
        # print(st_mp_conv.lattice.matrix)
        # print(st_mp_prim)
        # print(st_mp_prim.lattice)

        sg_after = st_mp_prim.get_space_group_info()

        if sg_before[0] != sg_after[0]:
            printlog('Attention! the space group was changed after primitive cell searching', sg_before, sg_after)
            printlog('I will save supercell in cif and reduce symprec to 0.01')
            st_mp_prim = st_mp
            symprec = 0.01

        if st_mp_prim:
            from pymatgen.io.cif import CifWriter
            cif = CifWriter(st_mp_prim, symprec = symprec)
            cif_name =  name_str+'.cif'
            cif.write_file(  join(sfolder, cif_name)  )
            printlog('Writing cif', cif_name)

        if 0:
            #get multiplication matrix which allows to obtain the supercell from primitive cell.
            #however this matrix is not integer which is not convinient.
            print(st_mp.lattice.matrix.round(2))
            print(st_mp_prim.lattice.matrix.round(2))

            mul_matrix = np.dot(st_mp.lattice.matrix, np.linalg.inv(st_mp_prim.lattice.matrix) )

            print(mul_matrix.round(1))

            rprimd = np.dot(mul_matrix, st_mp_prim.lattice.matrix  )

            print(rprimd.round(2))

        #write chg
        if 1:
            path_to_chg = cl.get_chg_file('CHGCAR')
            if path_to_chg:
                makedir( join(sfolder,'bin','dummy') )
                printlog('path to chgcar',path_to_chg)
                gz = '.gz'
                if gz not in path_to_chg:
                    gz = ''
                shutil.copyfile(path_to_chg, join( sfolder, 'bin', name_str+'.chg'+gz)  )


        #make dat
        #incars
        makedir(  join(sfolder, 'dat','dummy')  )
        incars = glob.glob(  join(cl.dir, '*INCAR*')  )
        # print(incars)
        for inc in incars:
            shutil.copy(  inc, join(sfolder, 'dat')  )


        #kpoints
        import json
        with open(  join(sfolder, 'dat', 'kpoints_for_kspacings.json'), 'w', newline = '') as fp:
            json.dump(header.struct_des[it].ngkpt_dict_for_kspacings, fp,)

        # print(cl.set.toJSON())


        #prepare for neb
        # makedir(sfolder+'neb_'+name_str+'/dummy')


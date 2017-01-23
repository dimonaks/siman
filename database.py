#Copyright Aksyonov D.A
from __future__ import division, unicode_literals, absolute_import 

import shelve, sys, datetime, shutil, tempfile, os, json

import header
from header import runBash, print_and_log, printlog
from classes import CalculationVasp
from small_functions import makedir
from set_functions import init_default_sets


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

    # shutil.copyfile(databasefile3, 'calc_copy.gdbm3')
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









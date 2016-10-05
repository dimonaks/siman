from __future__ import print_function

import header
from header import log, runBash, print_and_log

import shelve, sys, datetime, shutil


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

    # if scratch == True: databasefile =   '/scratch/aksenov/calc.s'
    
    log.write("\nLaunch at "+str( datetime.datetime.today() )+'\n')
    
    # mod = __import__("gdbm")

    # d = shelve.open(databasefile, protocol = 1)
    d = shelve.open(databasefile3, protocol = 1)
    # d = shelve.Shelf(mod.open(databasefile, protocol=1))


    try:
        calc              = d[calc_key]; 
        conv              = d[conv_key]; 
        varset            = d[varset_key]; 
        header.history    = d[history_key]
        header.struct_des = d[struct_des_key]

    except KeyError: 
        
        try: calc = d[calc_key] #dictionary of calculations
        except KeyError:
            log.write( "There is no database of calculations. I create new"); calc = {}
        


        try: conv = d[conv_key] #dictionary of convergence lists
        except KeyError:
            log.write( "There is no dictionary of convergence lists. I create new"); conv = {}   
        


        try: varset = d[varset_key] 
        except KeyError:
            log.write( "There is no dictionary of inputsets. I create new");  varset = {} 
        


        try: header.history = d[history_key] 
        except KeyError:
            header.history = ['Project started on '+datetime.date.today()]
            log.write( "There is still no history in database. The list is in header module ");
        


        try: header.struct_des = d[struct_des_key] 
        except KeyError:
            log.write( "There is no struct_des in database. The dict is in header module "); 

    d.close()
    #print history

    return calc, conv, varset, sys.getsizeof(d)


def write_database(calc, conv, varset, size_on_start = None):
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
    shutil.copyfile('calc.gdbm3', 'calc_copy.gdbm3')
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

        d = shelve.Shelf(dbm.open('calc.gdbm3', 'n'), protocol = 3) #Write dbm database for python3
        d[calc_key]       = calc
        d[conv_key]       = conv
        d[varset_key]     = varset
        d[history_key]    = header.history
        d[struct_des_key] = header.struct_des 
        d.close()   


    log.write("\nEnd of work at "+str(datetime.datetime.now())+'\n')

    log.close()
    

    #Update history file
    with  open('history','w') as his:
        #print history
        for i in header.history:
            #print i
            his.write(i+"\n")
    
    print("\nDatabase has been successfully updated\n")
    
    return




import header
from header import log, runBash, print_and_log

import shelve, sys, datetime


"""
Module contains utilities for project database management
and working with external databases

TODO:
1) in write_database() make update of history file more clever
"""

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

    databasefile = 'calc.s'

    # if scratch == True: databasefile =   '/scratch/aksenov/calc.s'
    
    log.write("\nLaunch at "+str( datetime.datetime.today() )+'\n')
    

    d = shelve.open(databasefile, protocol=1)
    
    try:
        calc = d['calc']; 
        conv = d['conv']; 
        varset = d['varset']; 
        header.history    = d['history']
        header.struct_des = d['struct_des']
        #print struct_des 
       #
    except KeyError: 
        
        try: calc = d['calc'] #dictionary of calculations
        except KeyError:
            log.write( "There is no database of calculations. I create new"); calc = {}
        


        try: conv = d['conv'] #dictionary of convergence lists
        except KeyError:
            log.write( "There is no dictionary of convergence lists. I create new"); conv = {}   
        


        try: varset = d['varset'] 
        except KeyError:
            log.write( "There is no dictionary of inputsets. I create new");  varset = {} 
        


        try: header.history = d['history'] 
        except KeyError:
            header.history = ['Project started on '+datetime.date.today()]
            log.write( "There is still no history in database. The list is in header module ");
        


        try: header.struct_des = d['struct_des'] 
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
    runBash("cp calc.s calc_copy.s") #create copy before writing


    d = shelve.open('calc.s', protocol=1) #Write database of calculations
    d['calc']       = calc
    d['conv']       = conv
    d['varset']     = varset
    d['history']    = header.history
    d['struct_des'] = header.struct_des 
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




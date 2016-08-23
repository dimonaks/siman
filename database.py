import header
from header import log, runBash
import shelve, sys, datetime


def read_database(scratch = False):
    #2. Read database of calculations
    #global history; 
    #global struct_des;
    databasefile = 'calc.s'
    if scratch == True: databasefile =   '/scratch/aksenov/calc.s'
    
    log.write("\nLaunch at "+str( datetime.datetime.today() )+'\n')
    d = shelve.open(databasefile, protocol=1)
    try:
        calc = d['calc']; 
        conv = d['conv']; 
        varset = d['varset']; 
        header.history = d['history']
        #struct_des = d['struct_des']
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
            log.write( "There is still no history in database. The list is in header module ");  #history = [] 
        #try: struct_des = d['struct_des'] 
        #except KeyError:
            #log.write( "There is still no struct_des in database. The dict is global "); # struct_des = {} 

    d.close()
    #print history

    return calc,conv,varset,sys.getsizeof(d)


def write_database(calc,conv,varset,size_on_start):
    #size_on_finish = sys.getsizeof(dbase)
    #if size_on_finish != size_on_start:
    runBash("cp calc.s calc_copy.s") #create copy before writing
    d = shelve.open('calc.s', protocol=1) #Write database of calculations
    #print struct_des 
    d['calc'] = calc
    d['conv'] = conv
    d['varset'] = varset
    d['history'] = header.history
    #d['struct_des'] = struct_des 

    d.close()
    log.write("\nEnd of work at "+str(datetime.datetime.now())+'\n')
    log.close()
    with  open('history','w') as his:
        #print history
        for i in header.history:
            #print i
            his.write(i+"\n")
    print("\nDatabase was succesfully updated\n")
    return

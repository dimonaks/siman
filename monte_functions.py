import os, subprocess
import math
import numpy as np
import copy
import datetime
import shutil
import traceback
import glob
from random import randint, random

class Structure():

    def __init__(self):
        self.natom = 0
        self.acell = [0,0,0]
        self.rprimd = []
        self.xcart = []
        self.xred  = []
        self.xcart10 = []
        self.xred10  = []
        self.typat = []
        self.typat1 = []
        self.typat2 = []
        self.typat10 = [] #vacancies
        self.etot = 0



def runBash(cmd):
    """Input - string; Executes Bash commands and returns stdout
Need: import subprocess
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    out = p.stdout.read().strip()
    return out  #This is the stdout from the shell command

def print_and_log(mystring):
    print mystring,
    #log.write(mystring)



def words(fileobj):
    """Generator of words. However does not allow to use methods of list for returned"""
    for line in fileobj:
        for word in line.split():
            yield word

def read_vectors(token, number_of_vectors, list_of_words):
    """Returns the list of numpy vectors for the last match"""

    number_of_matches = list_of_words.count( token )
    if number_of_matches == 0: 
        print_and_log("Warning token '"+token+"' was not found! return empty\n")
        return [None]

    if number_of_matches > 1:
        print_and_log("Warning token '"+token+"' was found more than one times\n")
        raise RuntimeError


    index = list_of_words.index(token, number_of_matches - 1 )     #Return the index of the last match
    #print list_of_words[index]
    list_of_vectors = []
    vector = np.zeros((3))
    for i in range(number_of_vectors):
        vector[0] = float(list_of_words[index + 1])
        vector[1] = float(list_of_words[index + 2])
        vector[2] = float(list_of_words[index + 3])
        index+=3
        list_of_vectors.append(vector.copy())
    return list_of_vectors


def read_lammps_xcart(token, number_of_vectors, list_of_words):
    """Returns the list of numpy vectors for the last match"""

    number_of_matches = list_of_words.count( token )
    if number_of_matches == 0: 
        print_and_log("Warning token '"+token+"' was not found! return empty\n")
        return [None]

    if number_of_matches > 1:
        print_and_log("Warning token '"+token+"' was found more than one times\n")
        raise RuntimeError


    index = list_of_words.index(token, number_of_matches - 1 )     #Return the index of the last match
    #print list_of_words[index]
    list_of_vectors = []
    vector = np.zeros((3))
    for i in range(number_of_vectors):
        vector[0] = float(list_of_words[index + 3])
        vector[1] = float(list_of_words[index + 4])
        vector[2] = float(list_of_words[index + 5])
        index+=5
        list_of_vectors.append(vector.copy())
    return list_of_vectors


def xred2xcart(xred, acell):
    xcart = []
    for x in xred:
        xr = x.copy()
        for i in 0,1,2:
            xr[i] = x[i] * acell[i]
        xcart.append (xr)
    return xcart   


def read_lammps_poscar(name):
    st = Structure()
    #print st.natom
    with open(name, "r") as initposcar:

        for i, line in enumerate(initposcar):
            if "ITEM: NUMBER OF ATOMS" in line:
                st.natom = int( initposcar.next() ) 

            if "ITEM: BOX BOUNDS pp pp pp" in line:
                for j in 0,1,2:
                    st.acell[j] = float(initposcar.next().split()[1])


            if "zs" in line:
                for i_at in range(st.natom):
                    string = initposcar.next()
                    #print string
                    st.typat.append(   int(string.split()[1])    ) 
                    st.xred.append(   np.asarray( [ float(x) for x in string.split()[2:5]  ] )  )

                print st.xred


                st.xcart = xred2xcart(st.xred, st.acell)
    return st


def write_lammps_poscar(st, name):
    """
    Writes conf.dump file from Structure() type
    """


    with open(name,'w') as poscar:
        poscar.write("ITEM: TIMESTEP\n")

        poscar.write("0\n")

        poscar.write("ITEM: NUMBER OF ATOMS\n")

        poscar.write(str(st.natom)  )

        poscar.write("\nITEM: BOX BOUNDS pp pp pp")

        poscar.write("\n0 "+str(st.acell[0])+"\n0 "+str(st.acell[1])+"\n0 "+str(st.acell[2])   )         

        poscar.write("\nITEM: ATOMS id type xs ys zs\n"    ) 
        for i in range(st.natom):
            #if st.typat[i] == 10: continue
            poscar.write(str(i+1)+" "+str(st.typat[i])+" "+str(st.xred[i][0])+" " \
            +str(st.xred[i][1]) +" "+str(st.xred[i][2])+"\n"      )           


def calc_energy():
    lammpsout = runBash("/home/soezen/./lmp_serial -in md.in").splitlines()

    for i, line in enumerate(lammpsout):
        if "Energy initial" in line:
            endenergy = float( lammpsout[i+1].split()[2] )
    return endenergy





def metropolis(E1, E2, T = 1):
    """
    Metropolis algorithm
    """
    decrease = False # energy reduction
    
    kb = 1.3806488*10**-23  / 1.6 * 10 **19
    dE = E2 - E1
    
    print("metropolis(): dE is ", dE)
    if dE < -0.000001:
        print( "dE is ", dE, "Accept!")
        decrease = True
    elif  1 > math.exp(-dE/kb/T) > random():
        print ("Accepted due to the temperature; exponent is ", math.exp(-dE/kb/T) )
        decrease = True
    else:
        print('Not accepted')

    return decrease




def exchange_atom_vac(cur):
    ir1 = randint(0, len(cur.xred)-1 ) #random Al or Fe atom
    #print "The number of atom to exchange", ir1+1
    ir2 = randint(0, len(cur.xred10)-1 ) #random vacancy

    x_new_vac = cur.xred[ir1].copy() #coordinates of new vacancy

    atom_typ  = cur.typat[ir1] #type of removed atom

    del cur.xred[ir1] #remove atom
    del cur.typat[ir1] #remove atom

    cur.xred.append( cur.xred10[ir2] ) #insert atom in old vacancy by adding it to the end of list
    cur.typat.append( atom_typ )          #type of atom is also added to the end of the list

    cur.xred10[ir2] = x_new_vac          #save new vacancy
    
    cur.xcart= xred2xcart(cur.xred, cur.acell)
    cur.xcart10= xred2xcart(cur.xred10, cur.acell)
    return cur

def exchange_atom_atom(cur, type = None ):

    alnum, fenum = [], []
    for i, t in enumerate(cur.typat):
        if t == 1: alnum.append(i)
        if t == 2: fenum.append(i)


    if type == 11:
        pass
    elif type == 22:
        pass

    elif type == 12: #change Al Fe
        ir1 = randint(0, len(alnum)-1 ) #random Al 
        ir2 = randint(0, len(fenum)-1 ) #random  Fe atom
        ir1 = alnum[ir1]
        ir2 = fenum[ir2]

    else: # change 
        ir1 = randint(0, len(cur.xred)-1 ) #random Al or Fe atom
        #print "The number of atom to exchange", ir1+1
        ir2 = randint(0, len(cur.xred)-1 ) #random Al or Fe atom
        while ir1 == ir2:
            ir2 = randint(0, len(cur.xred)-1 ) #random Al or Fe atom
    
    typat1 = cur.typat[ir1]
    cur.typat[ir1] = cur.typat[ir2]
    cur.typat[ir2] = typat1

    #st.xcart = xred2xcart(st.xred, acell)

    return cur

def write_xyz(st, path = '', repeat = 1):
    """Writes st structure in xyz format in the folder xyz/pat"""
  
    # rprimd = st.rprimd

    xcart  = st.xcart
    xcart10  = st.xcart10
    # xred   = st.xred
    typat = st.typat
    # znucl = st.znucl
    # name = st.name
    natom = st.natom +len(xcart10)

    # if natom != len(xred) != len(xcart) != len(typat) or len(znucl) != max(typat): 
    #     print "Error! write_xyz: check your arrays.\n\n"    
    
    # if name == '': name = 'noname'
    # if xcart == [] or len(xcart) != len(xred):
    #     print "Warining! write_xyz: len(xcart) != len(xred) making xcart from xred.\n"
    #     xcart = xred2xcart(xred, rprimd)
    #     #print xcart[1]


    xyzfile = 'xyz/'+path+".xyz"
    if not os.path.exists(os.path.dirname(xyzfile)):
        os.makedirs(os.path.dirname(xyzfile))

    #print_and_log("Writing xyz "+name+" \n")

    with open(xyzfile,'w') as f:
        for i in range(repeat):
            
            f.write(str(natom)+"\n")
            f.write(path+"\n")
            for i in range(st.natom):
                # typ = typat[i] - 1
                
                # z = int ( znucl[ typ ] )
                #print "typ", znucl
                #print "z", z
                if   typat[i] == 1:  f.write( "Al " )
                elif typat[i] == 2:  f.write( "Fe " )
                else:                f.write( "Pu " )
                f.write( "%s %s %s \n"%( xcart[i][0], xcart[i][1], xcart[i][2] ) )
            print len(xcart10)
            for i in range( len(xcart10)  ):
                f.write( "Pu %s %s %s \n"%( xcart10[i][0], xcart10[i][1], xcart10[i][2] ) )

    return


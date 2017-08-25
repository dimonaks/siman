# -*- coding: utf-8 -*-
#Copyright Aksyonov D.A

from __future__ import division, unicode_literals, absolute_import, print_function

"""
Siman - management of VASP calculations
Author: Aksyonov D.A.


!Internal units are Angstroms!
Bohrs should be converted during reading!

TODO:


"""

import os, subprocess, sys, shelve

import matplotlib as mpl



"""Global matplotlib control"""
# size = 22 #for one coloumn figures
size = 16 #for DOS
# size = 16 #for two coloumn figures
mpl.rc('font',family='Serif')
# mpl.rc('xtick', labelsize= size) 
# mpl.rc('ytick', labelsize= size) 
# mpl.rc('axes', labelsize = size) 
# mpl.rc('legend', fontsize= size) 
# mpl.rc('Axes.annotate', fontsize= size) #does not work
mpl.rcParams.update({'font.size': size})
mpl.rcParams.update({'mathtext.fontset': "stix"})
# plt.rcParams['mathtext.fontset'] = "stix"

#paths to libraries needed by siman
# sys.path.append('/home/dim/Simulation_wrapper/ase') #


history = []

try:
    from project_conf import *
    import project_conf
    siman_run = True
    log = open('log','a')
    warnings = 'neyY'
    warnings = 'Y'


except:
    # print('Some module is used separately; default_project_conf.py is used')
    mpl.use('agg') #switch matplotlib on or off; for running script using ssh
    siman_run = False
    from default_project_conf import *
    history.append('separate run')
    warnings = 'yY'


import matplotlib.pyplot as plt


calc_database = 'only_calc.gdbm3'

class CalcDict(dict):
    def __getitem__(self, key):
        # print(self)

        if dict.__contains__(self, key):
            # print('key', key, 'is  in self')
            val = dict.__getitem__(self, key)
        
        else:
            with shelve.open(calc_database, protocol = 3) as d:
                try:
                    # print(type(key)==str)
                    if type(key) == str:
                        # print('String key detected', key)
                        l = key.split('.')
                        if len(l) > 2:
                            key = ('.'.join(l[0:-2]), l[-2], int(l[-1]))
                        # print(key)
                    val = d[str(key)]
                    # print(len(d))
                    dict.__setitem__(self, key, val)
                    # print('reading ',str(key), 'from db')
                # print(val)
                # else:


                except:
                    val = None
        

        return val

    def __contains__(self, key):

        if dict.__contains__(self, key):
            return True
        
        else:
            with shelve.open(calc_database, protocol = 3) as d:
                # print('checking if key',key, 'is in db:', str(key) in d)
                # print(self)
                return str(key) in d



#Global variables
final_vasp_clean     = True 
copy_to_cluster_flag = True
close_run = False # alows to control close run file automatically after each add_loop
first_run = True  # needed to write header of run script
ssh_object = None # paramiko ssh_object
show = None
corenum = 1

calc = CalcDict()
# global db
db = calc
conv = {};
varset = {};
struct_des = {};








#Constants
to_ang = 0.52917721092
to_eV = 27.21138386
Ha_Bohr_to_eV_A = 51.4220641868956
kB_to_GPa = 0.1
eV_A_to_J_m = 16.021765
THz2eV = 0.00413566553853599


kB = 8.617e-5 # eV/K
TRANSITION_ELEMENTS = [22, 23, 25, 26, 27, 28]
ALKALI_ION_ELEMENTS = [3, 11, 19, 37]
MAGNETIC_ELEMENTS = [26, 27, 28]
# EXCLUDE_NODES = False
el_dict = {'octa':200, 'n':0, 'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10, 'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15, 'S':16, 'Cl':17, 'Ar':18, 'K':19, 'Ca':20, 'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, 'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 'Kr':36, 'Rb':37, 'Sr':38, 'Y':39, 'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48, 'In':49, 'Sn':50, 'Sb':51, 'Te':52, 'I':53, 'Xe':54, 'Cs':55, 'Ba':56, 'La':57, 'Ce':58, 'Pr':59, 'Nd':60, 'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64, 'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 'Yb':70, 'Lu':71, 'Hf':72, 'Ta':73, 'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80, 'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86, 'Fr':87, 'Ra':88, 'Ac':89, 'Th':90, 'Pa':91, 'U':92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96, 'Bk':97, 'Cf':98, 'Es':99, 'Fm':100, 'Md':101, 'No':102, 'Lr':103, 'Rf':104, 'Db':105, 'Sg':106, 'Bh':107, 'Hs':108, 'Mt':109, 'Ds':110, 'Rg':111, 'Cn':112, 'Uuq':114, 'Uuh':116, }
nu_dict = { 200:'octa', 0:'n', 1:'H', 2:'He', 3:'Li', 4:'Be', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F', 10:'Ne', 11:'Na', 12:'Mg', 13:'Al', 14:'Si', 15:'P', 16:'S', 17:'Cl', 18:'Ar', 19:'K', 20:'Ca', 21:'Sc', 22:'Ti', 23:'V', 24:'Cr', 25:'Mn', 26:'Fe', 27:'Co', 28:'Ni', 29:'Cu', 30:'Zn', 31:'Ga', 32:'Ge', 33:'As', 34:'Se', 35:'Br', 36:'Kr', 37:'Rb', 38:'Sr', 39:'Y', 40:'Zr', 41:'Nb', 42:'Mo', 43:'Tc', 44:'Ru', 45:'Rh', 46:'Pd', 47:'Ag', 48:'Cd', 49:'In', 50:'Sn', 51:'Sb', 52:'Te', 53:'I', 54:'Xe', 55:'Cs', 56:'Ba', 57:'La', 58:'Ce', 59:'Pr', 60:'Nd', 61:'Pm', 62:'Sm', 63:'Eu', 64:'Gd', 65:'Tb', 66:'Dy', 67:'Ho', 68:'Er', 69:'Tm', 70:'Yb', 71:'Lu', 72:'Hf', 73:'Ta', 74:'W', 75:'Re', 76:'Os', 77:'Ir', 78:'Pt', 79:'Au', 80:'Hg', 81:'Tl', 82:'Pb', 83:'Bi', 84:'Po', 85:'At', 86:'Rn', 87:'Fr', 88:'Ra', 89:'Ac', 90:'Th', 91:'Pa', 92:'U', 93:'Np', 94:'Pu', 95:'Am', 96:'Cm', 97:'Bk', 98:'Cf', 99:'Es', 100:'Fm', 101:'Md', 102:'No', 103:'Lr', 104:'Rf', 105:'Db', 106:'Sg', 107:'Bh', 108:'Hs', 109:'Mt', 110:'Ds', 111:'Rg', 112:'Cn', 114:'Uuq', 116:'Uuh', }


def print_and_log(*logstrings, **argdic):
    """
    '' - silent
    e - errors and warnings
    a - attentions
    m - minimalistic output of scientific procedures - only obligatory mess are shown
    M - maximalistic output of scientific procedures  
    debug_level importance:
        'n' - not important at all - for debugging
        ''  - almost all actions, no flag is needed
        'y' - important - major actions
        'Y' - super important, or output asked by user
    """
    end = '\n\n'# no argument for end, make one separate line
    
    debug_level  = 'e' #empty
    for key in argdic:
        if 'imp' in key:
            debug_level = argdic[key]
        
        if 'end' in key:
            end = argdic[key]


    mystring = ''
    for m in logstrings:
        mystring+=str(m)+' '


    if len(mystring.splitlines()) == 1:
        mystring = '-- '+mystring
    else:
        ''
        # mystring = '    '+mystring.replace('\n', '\n    ') 
    
    mystring+=end


    if 'Error' in mystring:# or 'Warning' in mystring:
        mystring+='\n\n\n'
    
    if 'Warning' in mystring or 'Attention' in mystring:
        debug_level = 'Y'

    for level in 'neyY':
        # print(level, debug_level)
        if (level in warnings and level in debug_level):
            print (mystring,  end = "")

    # if warnings:
    #     ''
    #     # print(debug_level)
    #     if 'n' in debug_level and 'n' not in warnings:
    #         pass
    #     else:
    #         print (mystring,  end = "")

    if siman_run:
        log.write(mystring)
    
    if 'Error!' in mystring:
        print (mystring)
        print ('Error! keyword was detected in message; invoking RuntimeError ')
        # sys.exit()
        raise RuntimeError

    return

printlog = print_and_log








def runBash(cmd, env = None, detached = False):
    """Input - string; Executes Bash commands and returns stdout
Need: import subprocess
    """
    if detached:
        stdout = None
        stderr = None
    else:
        stdout = subprocess.PIPE
        stderr = subprocess.STDOUT

    my_env = os.environ.copy()
    # my_env["PATH"] = "/opt/local/bin:/opt/local/sbin:" + my_env["PATH"]
    p = subprocess.Popen(cmd, executable='/bin/bash', shell=True, stdout=stdout, stderr = stderr, stdin = None, env = my_env)
    # print (cmd)
    # print 'Bash output is\n'+out
    # print ( str(out, 'utf-8') ) 
    out = ''
    try:
        out = p.stdout.read().strip()

        out = str(out, 'utf-8')
    except:
        pass

    return out  #This is the stdout from the shell command





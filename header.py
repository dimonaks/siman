# -*- coding: utf-8 -*-
#Copyright Aksyonov D.A

from __future__ import division, unicode_literals, absolute_import, print_function

"""
Siman - management of VASP calculations
Author: Aksyonov D.A.
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
    print('Some module is used separately; default_project_conf.py is used')
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
                    val = d[str(key)]
                    dict.__setitem__(self, key, val)
                # print('reading ',key, 'from db')
                # print(val)
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

calc = CalcDict()
conv = {};
varset = {};
struct_des = {};








#Constants
to_ang = 0.52917721092
to_eV = 27.21138386
Ha_Bohr_to_eV_A = 51.4220641868956
kB_to_GPa = 0.1
eV_A_to_J_m = 16.021765
kB = 8.617e-5 # eV/K
TRANSITION_ELEMENTS = [22, 23, 25, 26, 27, 28]
ALKALI_ION_ELEMENTS = [3, 11, 19]
MAGNETIC_ELEMENTS = [26, 27, 28]
# EXCLUDE_NODES = False


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
        mystring = '    '+mystring.replace('\n', '\n    ') 
    
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








def runBash(cmd, env = None):
    """Input - string; Executes Bash commands and returns stdout
Need: import subprocess
    """
    my_env = os.environ.copy()
    # my_env["PATH"] = "/opt/local/bin:/opt/local/sbin:" + my_env["PATH"]
    p = subprocess.Popen(cmd, executable='/bin/bash', shell=True, stdout=subprocess.PIPE, stderr = subprocess.STDOUT, env = my_env)
    out = p.stdout.read().strip()
    # print (cmd)
    # print 'Bash output is\n'+out
    # print ( str(out, 'utf-8') ) 
    try:
        out = str(out, 'utf-8')
    except:
        pass

    return out  #This is the stdout from the shell command




def red_prec(value, precision = 100.):
    a = value * precision
    return round(a)/1./precision




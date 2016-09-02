# -*- coding: utf-8 -*-
"""
Siman - management of VASP calculations
Author: Aksyonov D.A.
TODO:
1. Use our makedir() function


"""


import os, subprocess
import math
import numpy as np
import copy
import datetime
import shutil
import traceback
import glob
from operator import itemgetter
import sys
from math import exp
import optparse
import re
import colorsys
import pandas as pd
import matplotlib.pyplot as plt

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy

plt.rcParams['mathtext.fontset'] = "stix"


"""Global matplotlib control"""

# size = 22 #for one coloumn figures
size = 16 #for DOS
# size = 16 #for two coloumn figures
# mpl.rc('font',family='Times New Roman')
mpl.rc('font',family='Serif')
# mpl.rc('xtick', labelsize= size) 
# mpl.rc('ytick', labelsize= size) 
# mpl.rc('axes', labelsize = size) 
# mpl.rc('legend', fontsize= size) 
# mpl.rc('Axes.annotate', fontsize= size) #does not work
mpl.rcParams.update({'font.size': size})
# mpl.use('agg') #switch matplotlib on or off; for running script using ssh





"""
TODO

"""
#paths to libraries needed by siman
sys.path.append('/home/dim/Simulation_wrapper/ase') #

# sys.path.append('../../Simulation_wrapper/')
siman_run = True
try:
    from project_conf import *
    import project_conf
except:
    print 'Some module is used separatly'
    siman_run = False

# print 'header, geo_folder = ', geo_folder


#Global names
close_run = False # alows to control close run file automatically after each add_loop
calc = {};
conv = {};
varset = {};

struct_des = {};
history = []

if siman_run:
    log = open('log','a')
#history = open('log','a')
#Constants
to_ang = 0.52917721092
to_eV = 27.21138386
Ha_Bohr_to_eV_A = 51.4220641868956
kB_to_GPa = 0.1
eV_A_to_J_m = 16.021765
kB = 8.617e-5 # eV/K

warnings = True

def print_and_log(*logstrings, **debug_level):
    """
    '' - silent
    e - errors and warnings
    a - attentions
    m - minimalistic output of scientific procedures - only obligatory mess are shown
    M - maximalistic output of scientific procedures  
    """
    try:
        debug_level  = debug_level.values()[0]
    except:
        debug_level  = 'Mea'

    print debug_level
    mystring = ''
    for m in logstrings:
        mystring+=m

    if '\n' not in mystring[-2:-1]:
        mystring+='\n'

    if 'Error' in mystring or 'Warning' in mystring:
        mystring+='\n\n\n'
    

    if warnings:
    
        print '--'+mystring,

    log.write(mystring)

def runBash(cmd, env = None):
    """Input - string; Executes Bash commands and returns stdout
Need: import subprocess
    """
    my_env = os.environ.copy()
    # my_env["PATH"] = "/opt/local/bin:/opt/local/sbin:" + my_env["PATH"]
    p = subprocess.Popen(cmd, executable='/bin/bash', shell=True, stdout=subprocess.PIPE, env = my_env)
    out = p.stdout.read().strip()
    # print cmd
    # print 'Bash output is\n'+out
    return out  #This is the stdout from the shell command




def red_prec(value, precision = 100.):
    a = value * precision
    return round(a)/1./precision




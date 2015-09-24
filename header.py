# -*- coding: utf-8 -*-
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

sys.path.append('../../Simulation_wrapper/')

#Global names
#corenum = 8; queue = ' ' # Number of cores
cluster_address = 's'
corenum = 24; queue = ' -l cmmd '
close_run = False # alows to control close run file automatically after each add_loop
calc = {}
conv = {}
varset = {}
struct_des = {};
history = []
gb4_geo_folder = '/home/dim/Simulation_wrapper/gb4/out/'
project_path_cluster = "~/"
geo_folder = 'geo/'
# path_to_images = 'images/'
path_to_images = '/home/dim/gdrive/Наука/paper5/fig/'

path_to_jmol = '/home/dim/installed/jmol-14.2.12_2015.02.11/jmol.sh '

log = open('log','a')
#history = open('log','a')
#Constants
to_ang = 0.52917721092
to_eV = 27.21138386
Ha_Bohr_to_eV_A = 51.4220641868956
kB_to_GPa = 0.1
eV_A_to_J_m = 16.021765

def print_and_log(mystring):
    print mystring,
    log.write(mystring)

def runBash(cmd):
    """Input - string; Executes Bash commands and returns stdout
Need: import subprocess
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    out = p.stdout.read().strip()
    # print out
    return out  #This is the stdout from the shell command

class des():
    def __init__(self, sectionfolder = "forgot_folder", description = "forgot_description"):
        self.des = description
        self.sfolder = sectionfolder

class empty_struct():
    def __init__(self):
        pass

def update_des():
    """
    Naming conventions:

    endings:
    '_ml' - was used to show that this calculation uses manual equilibrium lattice determination and 
    contains several versions of identical structures with different
    lattice constants. Now not in use, because I always use this method. Usually 16 versions for hcp;

    '_r' - calculation with structure constructed for fitted lattice constants; 
    Now was replaced with '.f'; Usually one version.
    '.ur' - unrelaxed
    '.f'  - fitted
    '.fr' - means that current calculation based on the structure for which lattice constants were fitted and
    positions of atoms were relaxed. However see description to know for wich set they were fitted and relaxed.
    Calculations with '.f' and '.fr' can have different versions which are correspondig to different sets.
    
    .m - only matrix, all impurities were removed and matrix was freezed


    letters in name, wich are usually between didgits and element's names:
    b - stands for bulk, which denote ideal cells without boundaries.
    g - cells with grain boundary;
    v - means that impurity is in the volume of grain; far away from boundaries;
    i - means that impurity is close to interface plane (grain boundary)

    Versions:
    20 - usually means that lattice constatns was used from other calculation and this is very good assumtion.

    

    """

    #1.Description of structures
    struct_des['bcc']       = des("test",   "test cells for testing total drift v1 - ideal and v2 with atom shift")
    struct_des['gr221']       = des("GR",   "graphite 16 atoms Baskin1955")
    struct_des['O']       = des("box",   "O in box")
    struct_des['OO']       = des("box",   "O-O dimer in box; with initial velocity for constant change of length")
    struct_des['TiC']       = des("Ti-C",   "TiC(225)")
    struct_des['Ti2C227']       = des("Ti-C",   "Ti2C(227)")
    struct_des['Ti2O227']       = des("Ti-O",   "Ti2O(227)")
    
    struct_des['Ti2C164']       = des("Ti-C",   "Ti2C(164)")
    struct_des['Ti2O164']       = des("Ti-O",   "Ti2O(164)")



    #C1 twin
    struct_des['c1b']       = des("bulk",   "Bulk cell without boundary correspondig to cell with C1 twin")
    struct_des['c1b_ml']    = des("bulk",  "The same as c1b, but with 16 versions"   )
    struct_des['c1bC']      = des("C1/C",  " c1b_ml with carbon v1-16"   )
    struct_des['c1bC.f']      = des("C1/C",  " fitted v1"   )

    struct_des['c1bO_template'] = des("C1/O",  " template without oxygen; lattice constants from t111bO.f"   )
    struct_des['c1bO'] = des("C1/O",  "v1-16 from c1bC; from template with added oxygen v20 "   )
    struct_des['c1bO.f']      = des("C1/O",  " fitted v1"   )

    
    struct_des['c1gCv']    = des("C1/C",  " carbon in grain volume; v1-v5 the same -0.3 0.5"   )
    struct_des['c1gCi1']   = des("C1/C",  " carbon in place 1; v1-v5 the same as -0.3 0.5"   )

    struct_des['c1gBv']    = des("C1/B",  " boron in grain volume; v1-v5 the same -0.3 0.5"   )
    struct_des['c1gBi1']   = des("C1/B",  " boron in place 1; v1-v5 the same as -0.3 0.5"   )


    struct_des['c1gO_template'] = des("C1/O",  " template without oxygen; lattice constants from t111bO.f"   )
    struct_des['c1gOv']     = des("C1/O",  "oxygen in grain volume; v1-v5  "   )
    struct_des['c1gOi1']    = des("C1/O",  " oxygen in place 1; v1-v5 "   )




    struct_des['c1b_r']     = des("C1",   "The same as c1b, but with fitted lattice constants acell 2.9380 2.9380 4.6445 #angstrom #obtained using fit_hex(0.0002,0.0003,400,700, 'c1b_ml', '830', range(1,17), calc)")
    struct_des['c1g']       = des("C1", "Cell with C1 twin with the same lateral sizes as c1b_r and 10 different specific volumes")
    struct_des['c1g112']       = des("C1", "Cell with C1 twin with the same lattice const as c1b_r and 5 different specific volumes")
    struct_des['c1g121']       = des("C1", "Cell with C1 twin with the same lattice const as c1b_r and 5 different specific volumes")

    struct_des['c1gCO_template'] = des('C1/CO', 'template without impurities; lattice constants from t111gCO_template'   )
    struct_des['c1gCO_template_grainvol'] = des('C1/CO', 'template for generation of grain volume impurity cases'   )
    struct_des['c1gCO_template_coseg'] = des('C1/CO', 'template for generation of coseg cases'   )
    struct_des['c1gCO_template_segreg'] = des('C1/CO', 'template for generation of segreg cases'   )

    struct_des['c1gCOv1'] = des('C1/CO/c1g_gvol', 'corresponding configurations with impurities in volume; made from c1gCv.93kp7'   )
    struct_des['c1gCOv2'] = des('C1/CO/c1g_gvol', 'corresponding configurations with impurities in volume; made from c1gCv.93kp7'   )
    struct_des['c1gCOv3'] = des('C1/CO/c1g_gvol', 'corresponding configurations with impurities in volume; made from c1gCv.93kp7'   )
    struct_des['c1gCOv4'] = des('C1/CO/c1g_gvol', 'corresponding configurations with impurities in volume; made from c1gCv.93kp7'   )
    struct_des['c1gCOv5'] = des('C1/CO/c1g_gvol', 'corresponding configurations with impurities in volume; made from c1gCv.93kp7'   )
    struct_des['c1gCOv6'] = des('C1/CO/c1g_gvol', 'corresponding configurations with impurities in volume; made from c1gCv.93kp7'   )
    struct_des['c1gCOv7'] = des('C1/CO/c1g_gvol', 'corresponding configurations with impurities in volume; made from c1gCv.93kp7'   )
    struct_des['c1gCvOvms'] = des('C1/CO/c1g_gvol', 'corresponding configurations with impurities in volume; made from c1gCv.93kp7'   )

    struct_des['c1gCOi1.2-1'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gCOi2.2-1'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gCOi3.1-2'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gCOi4.2-1'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gCOi5.2-1'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gCOi6.1-2'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gCOi7.1-1'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gCOi8.2-2'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gOCi1.2-1'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gOCi2.2-1'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gOCi3.1-2'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gOCi4.2-1'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gOCi5.2-1'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gOCi6.1-2'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gOCi7.1-1'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gOCi8.2-2'] = des('C1/CO/c1g_coseg', 'co-segregation configurations; made from c1g.929'   )
    struct_des['c1gCOi10.1' ] = des('C1/CO/c1g_coseg',  'C and O in one pore; made from c1gCi1Ov O was removed'   )


    struct_des['c1gCi1Ov'] = des('C1/CO/c1g_segreg', 'corresponding segregation configurations; made also from c1g.929'   )
    struct_des['c1gCi2Ov'] = des('C1/CO/c1g_segreg', 'corresponding segregation configurations; made also from c1g.929'   )
    struct_des['c1gOi1Cv'] = des('C1/CO/c1g_segreg', 'corresponding segregation configurations; made also from c1g.929'   )
    struct_des['c1gOi2Cv'] = des('C1/CO/c1g_segreg', 'corresponding segregation configurations; made also from c1g.929'   )
    #T1 twin
    struct_des['t111b_ml']  = des("T1", "Bulk cell (only grainA 48 atoms) without boundary correspondig to cell with T1.1 (1) twin with 64 atoms [Wang2012]")
    struct_des['t111b_r']   = des("T1", "11-15 versions for optimization along x; Bulk cell with fitted lattice constants 2.9360 2.9360 4.6499 #angstrom fit_hex(0.0002,0.0003,400,600, 't111b_ml', '8301', range(1,17), calc) )")
    struct_des['t111g']     = des("T1", "Cell with T1.1 (1) twin  3 periods along x; lateral sizes as in t111b_r; 5 versions of free volume 64 atoms [Wang2012]; ")
    struct_des['t112g']     = des("T1", "Cell with T1.1 (2) twin, 4 periods along x; lateral sizes as in t111b_r; 5 versions of free volume 96 atoms [Wang2012]; ")
    
    struct_des['t111b']     = des("T1", "Bulk cell with initial lattice 't111b_r', full relaxation")
    struct_des['t111bC']    = des("T1/C",  " t111b_ml with carbon v1-16"   )
    struct_des['t111bC.f']    = des("T1/C",  " fitted v1"   )

    struct_des['t111bO']    = des("T1/O",  " t111b_ml with oxygen v1-16; relaxed positions used from t111bC"   )
    struct_des['t111bO.f.93kp9_template'] = des("T1/O",  " template without oxygen"   )
    struct_des['t111bO.f'] = des("T1/O",  "from template with oxygen "   )

    struct_des['t111gCv']    = des("T1/C",  " carbon in grain volume; v1-v5 volumes -0.3 0.5"   )
    struct_des['t111gCi1']   = des("T1/C",  " carbon in place 1; v1-v5 volumes -0.3 0.5"   )
    struct_des['t113gC_template'] = des("T1/C",  " empty template 96 atoms without carbon for preliminary relaxation of boundaries; v1-v5 volumes -0.3 0.5"   )
    struct_des['t113gCi1']   = des("T1/C",  " carbon in the same place as t111gCi1; v2; vol -0.1"   )
    struct_des['t113gCv']    = des("T1/C",  " carbon in the same place as t111gCv;  v2; vol -0.1"   )

    struct_des['t113gCO_template'] = des("T1/CO",  " empty template 96 atoms without impurites; the same sizes as t111gCO; v1-v5 volumes -0.3 0.5"   )
    struct_des['t113gCOi6.1-1is'] = des('T1/CO/t113g_coseg', 'co-segregation configurations; the same sizes as t111gCO '   )
    struct_des['t113gCvOvms'] = des('T1/CO/t113g_gvol', 'corresponding configurations with impurities in volume; same sizes as t111gCO'   )

    struct_des['t114g'] = des("T1",  "increased two times by r3 comparing to t111g; v1-v5 volumes -0.3 0.5"   )


    struct_des['t112gCO_template'] = des("T1/CO/t112", 'empty template 96 atoms without impurites; sizes calculated with calc_ac; v1-v5 volumes -0.3 0.5'   )
    struct_des['t112gCOi6.1-1is']  = des('T1/CO/t112', 'co-segregation configurations;'   )
    struct_des['t112gCvOvms']      = des('T1/CO/t112', 'corresponding configurations with impurities in volume;'   )

    struct_des['t114gCi1']   = des("T1/C",  " double of t111 by r3; 128 at; carbon in the same place as t111gCi1; v1-5; "   )
    struct_des['t114gCv']    = des("T1/C",  " double of t111 by r3; 128 at; carbon in the same place as t111gCv;  v1-5; "   )
    struct_des['t115gCi1']   = des("T1/C",  " double of t111 by r2; 128 at; carbon in the same place as t111gCi1; v1-5; vol -0.1"   )
    struct_des['t115gCv']    = des("T1/C",  " double of t111 by r2; 128 at; carbon in the same place as t111gCv;  v1-5; vol -0.1"   )


    struct_des['t111gO_template'] = des("T1/O",  " template without oxygen; lattice constants from t111bO.f"   )
    struct_des['t111gOv']     = des("T1/O",  "oxygen in grain volume; v1-v5  "   )
    struct_des['t111gOi1']    = des("T1/O",  " oxygen in place 1; v1-v5 "   )

    struct_des['t111gCO_template']    = des('T1/CO', 'template without impurities'   )
    struct_des['t111gCO_template_Cv'] = des('T1/CO', 'template with carbon in grain interior'   )
    struct_des['t111gCO_template_segreg'] = des('T1/CO', 'template with correct rprimd; segreg -  carbon and oxygen in grain interior; '   )
    struct_des['t111gCO_template_coseg'] = des('T1/CO', 'template with correct rprimd; coseg - pure'   )
    struct_des['t111gCO_template_grainvol'] = des('T1/CO', 'template with correct rprimd; grainvol - carbon in grain interior'  )

    struct_des['t111gCOv1'] = des('T1/CO/t111g_gvol', 'corresponding configurations with impurities in volume; made from t111gCv.93kp9'   )
    struct_des['t111gCOv2'] = des('T1/CO/t111g_gvol', 'corresponding configurations with impurities in volume; made from t111gCv.93kp9'   )
    struct_des['t111gCOv3is'] = des('T1/CO/t111g_gvol', 'corresponding configurations with impurities in volume; made from t111gCv.93kp9'   )
    struct_des['t111gCOv4'] = des('T1/CO/t111g_gvol', 'corresponding configurations with impurities in volume; made from t111gCv.93kp9'   )
    struct_des['t111gCOv5'] = des('T1/CO/t111g_gvol', 'corresponding configurations with impurities in volume; made from t111gCv.93kp9'   )
    struct_des['t111gCOv6'] = des('T1/CO/t111g_gvol', 'corresponding configurations with impurities in volume; made from t111gCv.93kp9'   )
    struct_des['t111gCvOvms'] = des('T1/CO/t111g_gvol', 'corresponding configurations with impurities in volume; made from t111gCv.93kp9'   )

    struct_des['t111gCOi1.4-3'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi2.2-1'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi3.4-1'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi4.4-4is'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi5.2-2is'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi6.1-1is'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi7.3-3is'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi8.2-3'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi9.4-4ms'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi10.2-2ms'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi11.4-3'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi12.1-3'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi13.4-2'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi14.2-1'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi15.4-4is'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gCOi16.2-2is'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gOCi1.4-3'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gOCi2.2-1'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gOCi3.4-1'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gOCi8.2-3'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gOCi11.4-3'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gOCi12.1-3'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gOCi13.4-2'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )
    struct_des['t111gOCi14.2-1'] = des('T1/CO/t111g_coseg', 'co-segregation configurations; made from t111g.9292'   )

    struct_des['t111gCi1Ov'] = des('T1/CO/t111g_segreg', 'corresponding segregation configurations; made from t111gCv.93kp9'   )
    struct_des['t111gCi2Ov'] = des('T1/CO/t111g_segreg', 'corresponding segregation configurations; made from t111gCv.93kp9'   )
    struct_des['t111gCi3Ov'] = des('T1/CO/t111g_segreg', 'corresponding segregation configurations; made from t111gCv.93kp9'   )
    struct_des['t111gCi4Ov'] = des('T1/CO/t111g_segreg', 'corresponding segregation configurations; made from t111gCv.93kp9'   )
    struct_des['t111gOi1Cv'] = des('T1/CO/t111g_segreg', 'corresponding segregation configurations; made from t111gCv.93kp9'   )
    struct_des['t111gOi2Cv'] = des('T1/CO/t111g_segreg', 'corresponding segregation configurations; made from t111gCv.93kp9'   )
    struct_des['t111gOi3Cv'] = des('T1/CO/t111g_segreg', 'corresponding segregation configurations; made from t111gCv.93kp9'   )
    struct_des['t111gOi4Cv'] = des('T1/CO/t111g_segreg', 'corresponding segregation configurations; made from t111gCv.93kp9'   )


    struct_des['t111sg']     = des("T1", "Cell with T1.1 (1) twin, but with shift r2/2 along r2; lateral sizes as in t111b_r; v1-5; 64 atoms [Wang2012]; ")


    struct_des['t111sgCO_template_coseg'] = des('T1/CO/shift', 'template with correct rprimd; segreg - carbon in grain interior; coseg - pure'   )
    struct_des['t111sgCO_template_segreg'] = des('T1/CO/shift', 'template with correct rprimd; segreg - carbon in grain interior; coseg - pure'   )
    struct_des['t111sgCO_template_grainvol'] = des('T1/CO/shift', 'template with correct rprimd; segreg - carbon in grain interior; coseg - pure'   )

    struct_des['t111sgCvOvms'] = des('T1/CO/shift/t111sg_gvol', 'corresponding configurations with impurities in volume; made from t111sg.93kp9'   )
    struct_des['t111sgCOv6'  ] = des('T1/CO/shift/t111sg_gvol', 'corresponding configurations with impurities in volume; made from csl71sg15.93'   )
    struct_des['t111sgCOv2'  ] = des('T1/CO/shift/t111sg_gvol', 'corresponding configurations with impurities in volume; made from csl71sg15.93'   )


    struct_des['t111sgCi1Ov'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )
    struct_des['t111sgCi2Ov'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )
    struct_des['t111sgCi3Ov'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )
    struct_des['t111sgCi4Ov'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )
    struct_des['t111sgCi5Ov'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )
    struct_des['t111sgCi6Ov'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )
    struct_des['t111sgCi7Ov'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )
    struct_des['t111sgOi1Cv'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )
    struct_des['t111sgOi2Cv'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )
    struct_des['t111sgOi3Cv'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )
    struct_des['t111sgOi4Cv'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )
    struct_des['t111sgOi5Cv'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )
    struct_des['t111sgOi6Cv'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )
    struct_des['t111sgOi7Cv'] = des('T1/CO/shift/t111sg_segreg', 'corresponding segregation configurations; made also from t111sg.93kp9'   )



    struct_des['t111sgCOi17.7-4'] = des('T1/CO/shift/t111sg_coseg', 'co-segregation configurations; made from t111sg.93kp9'   )
    struct_des['t111sgCOi28.4-5'] = des('T1/CO/shift/t111sg_coseg', 'co-segregation configurations; made from t111sg.93kp9'   )
    struct_des['t111sgCOi29.4-6'] = des('T1/CO/shift/t111sg_coseg', 'co-segregation configurations; made from t111sg.93kp9'   )
    struct_des['t111sgCOi30.7-5'] = des('T1/CO/shift/t111sg_coseg', 'co-segregation configurations; made from t111sg.93kp9'   )
    struct_des['t111sgOCi8.4-7'  ] = des('T1/CO/shift/t111sg_coseg', 'co-segregation configurations; made from t111sg.93kp9'   )
    struct_des['t111sgOCi9.6-2' ] = des('T1/CO/shift/t111sg_coseg', 'co-segregation configurations; made from t111sg.93kp9'   )
    struct_des['t111sgOCi15.6-6'] = des('T1/CO/shift/t111sg_coseg', 'co-segregation configurations; made from t111sg.93kp9'   )
    struct_des['t111sgOCi22.6-3'] = des('T1/CO/shift/t111sg_coseg', 'co-segregation configurations; made from t111sg.93kp9'   )
    struct_des['t111sgOCi25.6-7'] = des('T1/CO/shift/t111sg_coseg', 'co-segregation configurations; made from t111sg.93kp9'   )
    struct_des['t111sgOCi27.6-2'] = des('T1/CO/shift/t111sg_coseg', 'co-segregation configurations; made from t111sg.93kp9'   )











    struct_des['t21b_ml']   = des("T2", "Bulk cell (only grainA 44 atoms) without boundary correspondig to cell with T2 (1) twin with 64 atoms [Lane2011]; versions 1-16")
    struct_des['t21bC']     = des("T2/C",  " t21b_ml with carbon v1-16"   )
    struct_des['t21bC.f']   = des("T2/C",  " fitted v1"   )

    struct_des['t21bO_template'] = des("T2/O",  " template without oxygen; lattice constants from t111bO.f"   )
    struct_des['t21bO']     = des("T2/O",  "v1-16 from t21bC; from template with added oxygen v20 "   )
    struct_des['t21bO.f']   = des("T2/O",  " fitted v1"   )


    struct_des['t21b_r']    = des("T2", "Bulk Ti acell 2.9380 2.9380 4.6427 #after fit_hex(0.0002,0.0003,400,600, 't21b_ml', '83', range(1,17), calc) correspondig to cell with T2 (1) twin with 64 atoms [Lane2011]")
    struct_des['t21g']      = des("T2", "Cell with T2 (1) twin with 64 atoms [Lane2011]; v 1-5 -0.3 +0.3")
    struct_des['t21gCv']    = des("T2/C",  " carbon in grain volume; v1-v5 the same as t21g"   )
    struct_des['t21gCi1']   = des("T2/C",  " carbon in place 1; v1-v5 the same as t21g"   )

    struct_des['t21gO_template'] = des("T2/O",  " template without oxygen; lattice constants from t111bO.f"   )
    struct_des['t21gOv']     = des("T2/O",  "oxygen in grain volume; v1-v5  "   )
    struct_des['t21gOi1']    = des("T2/O",  " oxygen in place 1; v1-v5 "   )

    struct_des['t21gCO_template_grainvol'] = des('T2/CO', 'template with correct rprimd; segreg - carbon in grain interior; coseg - pure'   )
    struct_des['t21gCO_template_coseg'] = des('T2/CO', 'template with correct rprimd; segreg - carbon in grain interior; coseg - pure'   )
    struct_des['t21gCO_template_segreg'] = des('T2/CO', 'template with correct rprimd; segreg - carbon in grain interior; coseg - pure'   )

    struct_des['t21gCOv2'] = des('T2/CO/t21g_gvol', 'corresponding configurations with impurities in volume; made from t21gCv.93'   )
    struct_des['t21gCOv6'] = des('T2/CO/t21g_gvol', 'corresponding configurations with impurities in volume; made from t21gCv.93'   )
    struct_des['t21gCvOvms'] = des('T2/CO/t21g_gvol', 'corresponding configurations with impurities in volume; made from t21gCv.93'   )


    struct_des['t21gCOi1.4-3' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi2.3-2' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi3.2-2' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi4.2-1' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi5.3-1' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi6.4-2' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi7.2-1' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi8.4-1' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi9.4-2' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi10.3-3'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi11.4-3'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi12.3-2'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi13.4-1'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi14.1-4'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi15.3-2'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gCOi16.4-2'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi1.4-3' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi2.3-2' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi3.2-2' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi4.2-1' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi5.3-1' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi6.4-2' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi7.2-1' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi8.4-1' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi9.4-2' ] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi10.3-3'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi11.4-3'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi12.3-2'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi13.4-1'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi14.1-4'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi15.3-2'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )
    struct_des['t21gOCi16.4-2'] = des('T2/CO/t21g_coseg', 'co-segregation configurations; made from t21g.93'   )

    struct_des['t21gCi1Ov'] = des('T2/CO/t21g_segreg', 'corresponding segregation configurations; made also from t21g.93'   )
    struct_des['t21gCi2Ov'] = des('T2/CO/t21g_segreg', 'corresponding segregation configurations; made also from t21g.93'   )
    struct_des['t21gCi3Ov'] = des('T2/CO/t21g_segreg', 'corresponding segregation configurations; made also from t21g.93'   )
    struct_des['t21gCi4Ov'] = des('T2/CO/t21g_segreg', 'corresponding segregation configurations; made also from t21g.93'   )
    struct_des['t21gOi1Cv'] = des('T2/CO/t21g_segreg', 'corresponding segregation configurations; made also from t21g.93'   )
    struct_des['t21gOi2Cv'] = des('T2/CO/t21g_segreg', 'corresponding segregation configurations; made also from t21g.93'   )
    struct_des['t21gOi3Cv'] = des('T2/CO/t21g_segreg', 'corresponding segregation configurations; made also from t21g.93'   )
    struct_des['t21gOi4Cv'] = des('T2/CO/t21g_segreg', 'corresponding segregation configurations; made also from t21g.93'   )











    struct_des['csl71b_r']  = des("CSL7", "bulk; 56 atoms; corresponds to csl71g; a from t111b, c from hs221")
    struct_des['csl71bC']  = des("CSL7/C", "bulk; 57 atoms; Lattice constants as csl71b_r; to check voronoi volume of ideal octapore ")

    struct_des['csl71g']    = des("CSL7", "twin; v2-6; 52 atoms; a from t111b, c from hs221")
    struct_des['csl71sg5']    = des("CSL7", "twin shifted; v5-1; 52 atoms; a from t111b, c from hs221")
    struct_des['csl71sg10']    = des("CSL7", "twin shifted; v5-1; 52 atoms; a from t111b, c from hs221")
    struct_des['csl71sg15']    = des("CSL7", "twin shifted; v5-1; 52 atoms; a from t111b, c from hs221")
    struct_des['csl71sg10z2y']    = des("CSL7", "twin shifted; v5-1; 52 atoms; a from t111b, c from hs221")


    struct_des['csl71sg15CO_template_grainvol'] = des('CSL7/CO', 'template with correct rprimd; segreg - carbon in grain interior; coseg - pure'   )
    struct_des['csl71sg15CO_template_coseg'] = des('CSL7/CO', 'template with correct rprimd; segreg - carbon in grain interior; coseg - pure'   )
    struct_des['csl71sg15CO_template_segreg'] = des('CSL7/CO', 'template with correct rprimd; segreg - carbon in grain interior; coseg - pure'   )


    struct_des['csl71sgCvOvms'] = des('CSL7/CO/csl71sg_gvol', 'corresponding configurations with impurities in volume; made from csl71sg15.93'   )

    struct_des['csl71sgCOv5'] = des('CSL7/CO/csl71sg_gvol', 'corresponding configurations with impurities in volume; made from csl71sg15.93'   )
    struct_des['csl71sgCOv9'] = des('CSL7/CO/csl71sg_gvol', 'corresponding configurations with impurities in volume; made from csl71sg15.93'   )

    struct_des['csl71sgCi1Ov'] = des('CSL7/CO/csl71sg_segreg', 'corresponding segregation configurations; made also from csl71sg15.93'   )
    struct_des['csl71sgCi2Ov'] = des('CSL7/CO/csl71sg_segreg', 'corresponding segregation configurations; made also from csl71sg15.93'   )
    struct_des['csl71sgCi3Ov'] = des('CSL7/CO/csl71sg_segreg', 'corresponding segregation configurations; made also from csl71sg15.93'   )
    struct_des['csl71sgCi4Ov'] = des('CSL7/CO/csl71sg_segreg', 'corresponding segregation configurations; made also from csl71sg15.93'   )
    struct_des['csl71sgCi5Ov'] = des('CSL7/CO/csl71sg_segreg', 'corresponding segregation configurations; made also from csl71sg15.93'   )
    struct_des['csl71sgCi6Ov'] = des('CSL7/CO/csl71sg_segreg', 'corresponding segregation configurations; made also from csl71sg15.93'   )
    struct_des['csl71sgOi1Cv'] = des('CSL7/CO/csl71sg_segreg', 'corresponding segregation configurations; made also from csl71sg15.93'   )
    struct_des['csl71sgOi2Cv'] = des('CSL7/CO/csl71sg_segreg', 'corresponding segregation configurations; made also from csl71sg15.93'   )
    struct_des['csl71sgOi3Cv'] = des('CSL7/CO/csl71sg_segreg', 'corresponding segregation configurations; made also from csl71sg15.93'   )
    struct_des['csl71sgOi4Cv'] = des('CSL7/CO/csl71sg_segreg', 'corresponding segregation configurations; made also from csl71sg15.93'   )
    struct_des['csl71sgOi5Cv'] = des('CSL7/CO/csl71sg_segreg', 'corresponding segregation configurations; made also from csl71sg15.93'   )
    struct_des['csl71sgOi6Cv'] = des('CSL7/CO/csl71sg_segreg', 'corresponding segregation configurations; made also from csl71sg15.93'   )



    struct_des['csl71sgCOi1.6-4' ] = des('CSL7/CO/csl71sg_coseg', 'co-segregation configurations; made from csl71sg15.93'   )
    struct_des['csl71sgCOi11.4-4'] = des('CSL7/CO/csl71sg_coseg', 'co-segregation configurations; made from csl71sg15.93'   )
    struct_des['csl71sgCOi22.1-2'] = des('CSL7/CO/csl71sg_coseg', 'co-segregation configurations; made from csl71sg15.93'   )
    struct_des['csl71sgOCi6.5-4'  ] = des('CSL7/CO/csl71sg_coseg', 'co-segregation configurations; made from csl71sg15.93'   )
    struct_des['csl71sgOCi7.1-6'  ] = des('CSL7/CO/csl71sg_coseg', 'co-segregation configurations; made from csl71sg15.93'   )


    struct_des['csl77gCv']  = des("CSL7/C", "carbon at GB; from abinit")
    struct_des['csl77gCi6']    = des("CSL7/C", "carbon in bulk; from abinit")
    struct_des['csl77gCi7']    = des("CSL7/C", "testing search of pores in csl77gCv;")
    struct_des['csl77gCv.m']  = des("CSL7/C",    "noimp; v1")
    struct_des['csl77gCi6.m']    = des("CSL7/C", "noimp; v1")


    #without imp:
    struct_des['t111gCvOvms.m'] = des('T1/CO/t111g_segreg_m', 'without imp; made from t111gCvOvms.93kp9; ver are the same'   )
    struct_des['t111gCi1Ov.m'] = des('T1/CO/t111g_segreg_m', 'without imp; made from t111gCi1Ov.93kp9; ver are the same'   )
    struct_des['t111gCi2Ov.m'] = des('T1/CO/t111g_segreg_m', 'without imp; made from t111gCi2Ov.93kp9; ver are the same'   )
    struct_des['t111gCi3Ov.m'] = des('T1/CO/t111g_segreg_m', 'without imp; made from t111gCi3Ov.93kp9; ver are the same'   )
    struct_des['t111gCi4Ov.m'] = des('T1/CO/t111g_segreg_m', 'without imp; made from t111gCi4Ov.93kp9; ver are the same'   )
    struct_des['t111gOi1Cv.m'] = des('T1/CO/t111g_segreg_m', 'without imp; made from t111gOi1Cv.93kp9; ver are the same'   )
    struct_des['t111gOi2Cv.m'] = des('T1/CO/t111g_segreg_m', 'without imp; made from t111gOi2Cv.93kp9; ver are the same'   )
    struct_des['t111gOi3Cv.m'] = des('T1/CO/t111g_segreg_m', 'without imp; made from t111gOi3Cv.93kp9; ver are the same'   )
    struct_des['t111gOi4Cv.m'] = des('T1/CO/t111g_segreg_m', 'without imp; made from t111gOi4Cv.93kp9; ver are the same'   )
    struct_des['t111sgCvOvms.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgCvOvms.93kp9; ver are the same'   )
    struct_des['t111sgCi1Ov.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgCi1Ov.93kp9; ver are the same'   )
    struct_des['t111sgCi2Ov.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgCi2Ov.93kp9; ver are the same'   )
    struct_des['t111sgCi3Ov.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgCi3Ov.93kp9; ver are the same'   )
    struct_des['t111sgCi4Ov.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgCi4Ov.93kp9; ver are the same'   )
    struct_des['t111sgCi5Ov.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgCi5Ov.93kp9; ver are the same'   )
    struct_des['t111sgCi6Ov.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgCi6Ov.93kp9; ver are the same'   )
    struct_des['t111sgCi7Ov.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgCi7Ov.93kp9; ver are the same'   )
    struct_des['t111sgOi1Cv.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgOi1Cv.93kp9; ver are the same'   )
    struct_des['t111sgOi2Cv.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgOi2Cv.93kp9; ver are the same'   )
    struct_des['t111sgOi3Cv.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgOi3Cv.93kp9; ver are the same'   )
    struct_des['t111sgOi4Cv.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgOi4Cv.93kp9; ver are the same'   )
    struct_des['t111sgOi5Cv.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgOi5Cv.93kp9; ver are the same'   )
    struct_des['t111sgOi6Cv.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgOi6Cv.93kp9; ver are the same'   )
    struct_des['t111sgOi7Cv.m'] = des('T1/CO/t111sg_segreg_m', 'without imp; made from t111sgOi7Cv.93kp9; ver are the same'   )
    struct_des['c1gCvOvms.m'] = des('C1/CO/c1g_segreg_m', 'without imp; made from c1gCvOvms.93kp7; ver are the same'   )
    struct_des['c1gCi1Ov.m'] = des('C1/CO/c1g_segreg_m', 'without imp; made from c1gCi1Ov.93kp7; ver are the same'   )
    struct_des['c1gCi2Ov.m'] = des('C1/CO/c1g_segreg_m', 'without imp; made from c1gCi2Ov.93kp7; ver are the same'   )
    struct_des['c1gOi1Cv.m'] = des('C1/CO/c1g_segreg_m', 'without imp; made from c1gOi1Cv.93kp7; ver are the same'   )
    struct_des['c1gOi2Cv.m'] = des('C1/CO/c1g_segreg_m', 'without imp; made from c1gOi2Cv.93kp7; ver are the same'   )
    struct_des['t21gCvOvms.m'] = des('T2/CO/t21g_segreg_m', 'without imp; made from t21gCvOvms.93; ver are the same'   )
    struct_des['t21gCi1Ov.m'] = des('T2/CO/t21g_segreg_m', 'without imp; made from t21gCi1Ov.93; ver are the same'   )
    struct_des['t21gCi2Ov.m'] = des('T2/CO/t21g_segreg_m', 'without imp; made from t21gCi2Ov.93; ver are the same'   )
    struct_des['t21gCi3Ov.m'] = des('T2/CO/t21g_segreg_m', 'without imp; made from t21gCi3Ov.93; ver are the same'   )
    struct_des['t21gCi4Ov.m'] = des('T2/CO/t21g_segreg_m', 'without imp; made from t21gCi4Ov.93; ver are the same'   )
    struct_des['t21gOi1Cv.m'] = des('T2/CO/t21g_segreg_m', 'without imp; made from t21gOi1Cv.93; ver are the same'   )
    struct_des['t21gOi2Cv.m'] = des('T2/CO/t21g_segreg_m', 'without imp; made from t21gOi2Cv.93; ver are the same'   )
    struct_des['t21gOi3Cv.m'] = des('T2/CO/t21g_segreg_m', 'without imp; made from t21gOi3Cv.93; ver are the same'   )
    struct_des['t21gOi4Cv.m'] = des('T2/CO/t21g_segreg_m', 'without imp; made from t21gOi4Cv.93; ver are the same'   )
    struct_des['csl71sgCvOvms.m'] = des('CSL7/CO/csl71sg_segreg_m', 'without imp; made from csl71sgCvOvms.93; ver are the same'   )
    struct_des['csl71sgCi1Ov.m'] = des('CSL7/CO/csl71sg_segreg_m', 'without imp; made from csl71sgCi1Ov.93; ver are the same'   )
    struct_des['csl71sgCi2Ov.m'] = des('CSL7/CO/csl71sg_segreg_m', 'without imp; made from csl71sgCi2Ov.93; ver are the same'   )
    struct_des['csl71sgCi3Ov.m'] = des('CSL7/CO/csl71sg_segreg_m', 'without imp; made from csl71sgCi3Ov.93; ver are the same'   )
    struct_des['csl71sgCi4Ov.m'] = des('CSL7/CO/csl71sg_segreg_m', 'without imp; made from csl71sgCi4Ov.93; ver are the same'   )
    struct_des['csl71sgCi5Ov.m'] = des('CSL7/CO/csl71sg_segreg_m', 'without imp; made from csl71sgCi5Ov.93; ver are the same'   )
    struct_des['csl71sgCi6Ov.m'] = des('CSL7/CO/csl71sg_segreg_m', 'without imp; made from csl71sgCi6Ov.93; ver are the same'   )
    struct_des['csl71sgOi1Cv.m'] = des('CSL7/CO/csl71sg_segreg_m', 'without imp; made from csl71sgOi1Cv.93; ver are the same'   )
    struct_des['csl71sgOi2Cv.m'] = des('CSL7/CO/csl71sg_segreg_m', 'without imp; made from csl71sgOi2Cv.93; ver are the same'   )
    struct_des['csl71sgOi3Cv.m'] = des('CSL7/CO/csl71sg_segreg_m', 'without imp; made from csl71sgOi3Cv.93; ver are the same'   )
    struct_des['csl71sgOi4Cv.m'] = des('CSL7/CO/csl71sg_segreg_m', 'without imp; made from csl71sgOi4Cv.93; ver are the same'   )
    struct_des['csl71sgOi5Cv.m'] = des('CSL7/CO/csl71sg_segreg_m', 'without imp; made from csl71sgOi5Cv.93; ver are the same'   )
    struct_des['csl71sgOi6Cv.m'] = des('CSL7/CO/csl71sg_segreg_m', 'without imp; made from csl71sgOi6Cv.93; ver are the same'   )


    #Additional for dos
    struct_des['c1gCi1Ov.r'] = des('C1/CO/c1g_segreg', 'relaxed'   )
    struct_des['c1gCi2Ov.r'] = des('C1/CO/c1g_segreg', 'relaxed'   )
    struct_des['c1gOi1Cv.r'] = des('C1/CO/c1g_segreg', 'relaxed'   )
    struct_des['c1gOi2Cv.r'] = des('C1/CO/c1g_segreg', 'relaxed'   )
    struct_des['t111gCi2Ov.r'] = des('T1/CO/t111g_segreg', 'relaxed'   )
    struct_des['t111gCi3Ov.r'] = des('T1/CO/t111g_segreg', 'relaxed'   )
    struct_des['t111gOi2Cv.r'] = des('T1/CO/t111g_segreg', 'relaxed'   )
    struct_des['t111gOi3Cv.r'] = des('T1/CO/t111g_segreg', 'relaxed'   )
    struct_des['t111sgCi6Ov.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['t111sgOi6Cv.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['csl71sgCi4Ov.r'] = des('CSL7/CO/csl71sg_segreg', 'relaxed'   )
    struct_des['csl71sgOi4Cv.r'] = des('CSL7/CO/csl71sg_segreg', 'relaxed'   )

    struct_des['t21gCi1Ov.r'] = des('T2/CO/t21g_segreg', 'relaxed'   )
    struct_des['t21gCi4Ov.r'] = des('T2/CO/t21g_segreg', 'relaxed'   )
    struct_des['t21gOi1Cv.r'] = des('T2/CO/t21g_segreg', 'relaxed'   )
    struct_des['t21gOi4Cv.r'] = des('T2/CO/t21g_segreg', 'relaxed'   )



    struct_des['t111gCi1Ov.r'] = des('T1/CO/t111g_segreg', 'relaxed'   )
    struct_des['t111gCi4Ov.r'] = des('T1/CO/t111g_segreg', 'relaxed'   )
    struct_des['t111gOi1Cv.r'] = des('T1/CO/t111g_segreg', 'relaxed'   )
    struct_des['t111gOi4Cv.r'] = des('T1/CO/t111g_segreg', 'relaxed'   )
    struct_des['t111sgCi1Ov.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['t111sgCi2Ov.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['t111sgCi3Ov.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['t111sgCi4Ov.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['t111sgCi5Ov.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['t111sgCi7Ov.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['t111sgOi1Cv.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['t111sgOi2Cv.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['t111sgOi3Cv.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['t111sgOi4Cv.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['t111sgOi5Cv.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['t111sgOi7Cv.r'] = des('T1/CO/shift/t111sg_segreg', 'relaxed'   )
    struct_des['t21gCi2Ov.r'] = des('T2/CO/t21g_segreg', 'relaxed'   )
    struct_des['t21gCi3Ov.r'] = des('T2/CO/t21g_segreg', 'relaxed'   )
    struct_des['t21gOi2Cv.r'] = des('T2/CO/t21g_segreg', 'relaxed'   )
    struct_des['t21gOi3Cv.r'] = des('T2/CO/t21g_segreg', 'relaxed'   )





    struct_des['c1gCvOvms.r'] = des('C1/CO/c1g_gvol', 'relaxed'   )
    struct_des['t111sgCvOvms.r'] = des('T1/CO/shift/t111sg_gvol', 'relaxed'   )
    struct_des['t111gCvOvms.r'] = des('T1/CO/t111g_gvol', 'relaxed'   )
    struct_des['csl71sgCvOvms.r'] = des('CSL7/CO/csl71sg_gvol', 'relaxed'   )
    struct_des['t21gCvOvms.r'] = des('T2/CO/t21g_gvol', 'relaxed'   )






    #scaled 
    struct_des['c1gCi1Ov.r2d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['c1gCi2Ov.r2d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['c1gOi2Cv.r2d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['c1gCvOvms.r2d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['t111sgCi6Ov.r2d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['t111sgOi6Cv.r2d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['t111sgCvOvms.r2d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['csl71sgCi4Ov.r2d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['csl71sgOi4Cv.r2d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['csl71sgCvOvms.r2d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['t21gCi4Ov.r2d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['t21gOi4Cv.r2d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['t21gCvOvms.r2d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )


    struct_des['c1gCi1Ov.r3d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['c1gCi2Ov.r3d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['c1gOi2Cv.r3d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['c1gCvOvms.r3d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['t111sgCi6Ov.r3d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['t111sgOi6Cv.r3d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['t111sgCvOvms.r3d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['csl71sgCi4Ov.r3d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['csl71sgOi4Cv.r3d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['csl71sgCvOvms.r3d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['t21gCi4Ov.r3d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['t21gOi4Cv.r3d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )
    struct_des['t21gCvOvms.r3d'] = des('scaled', 'scaled only matrix; the sizes are doubled without any manipulations'   )




    struct_des['o1b']       = des("pure", "Ortogonal bulk cell for convergence tests")
    struct_des['o1b_ml']    = des("pure",  "The same as o1b, but with 25 versions of lattice parameters \
        for manual lattice (ml) search")

    
    struct_des['hs221']     = des("H",   "hexagonal cell with bulk Ti hcp; 8 atoms; ")
    struct_des['hs221.f']   = des("H",   "fitted;")
    struct_des['hs221C']    = des("H",   "hexagonal cell with bulk Ti hcp  and one carbon atom in central octapore; 9 atoms; ")
    struct_des['hs221C.f']  = des("H",   "fitted by 16 cells; ")

    struct_des['hs221O']    = des("H/O",   "hexagonal cell with bulk Ti hcp  and one oxygen atom in central octapore; 9 atoms; ")
    struct_des['hs221O.f.93_template'] = des("H/O", "fitted template without oxygen")
    struct_des['hs221O.f'] = des("H/O", "obtained from template by adding oxygen")
    struct_des['hs221O.fr'] = des("H/O", "atomic positions were relaxed in previous calculations;")

    struct_des['hs332']     = des("H",   "hexagonal cell with bulk Ti hcp; v1-16; 36 atoms; ")
    struct_des['hs332.f']   = des("H",   "fitted 0.0002,0.0003,400,600; v1; 36 atoms; ")

    struct_des['hs332C']    = des("H",   "one C in hexagonal cell with bulk Ti hcp; v1-16; 37 atoms; ")
    struct_des['hs332C.f']  = des("H",   "fitted &0.0002 &0.0003 &400 &600; v1; 37 atoms; ")
    
    struct_des['hs332O']    = des("H/O",   "made from hs332C.93 by replacing impurity; v1-16; 37 atoms; ")
    struct_des['hs332O.f.93_template'] = des("H/O", "fitted template without oxygen 36 atoms")
    struct_des['hs332O.f']  = des("H/O",   "obtained from template by adding oxygen; v1; 37 atoms; ")


    struct_des['hs443']     = des("H",   "hex cell with bulk Ti hcp; v1-16; v20; 96 atoms; ")
    struct_des['hs443.f']     = des("H",   "fitted; 96 atoms; ")
    struct_des['hs443C']    = des("H",   "one C in hexagonal cell with bulk Ti hcp; v1-16; 97 atoms; ")
    struct_des['hs443C.r']    = des("H",   " relaxed; for DOS; one C in hexagonal cell with bulk Ti hcp; v1-16; 97 atoms; ")
    struct_des['hs443C.f']  = des("H",   "fitted; v1; 97 atoms; ")
    struct_des['hs443C.fm']  = des("H/C",   "fitted relaxed and C removed; v1; 96 atoms; ")
    struct_des['hs443C.fr']  = des("H/C",   "fitted and relaxed; v1; 97 atoms; ")
    struct_des['hs443C.m']  = des("H/C",   "without imp; v1-16 and v20 for fitted; 96 atoms; ")


    struct_des['hs443O']    = des("H/O",   "made from hs443C.93 by replacing impurity; v1-16; 97 atoms; ")
    struct_des['hs443O.r']    = des("H/O",   "relaxed; made from hs443C.93 by replacing impurity; v1-16; 97 atoms; ")
    struct_des['hs443O.ur']    = des("H/O",   "unrelaxed; impurity added to hs443 ; v1-16; 97 atoms; ")

    # struct_des['hs443O.f.93_template'] = des("H/O", "fitted template without oxygen 96 atoms")
    struct_des['hs443O.f']  = des("H/O",   "obtained from template by adding oxygen; v1; 97 atoms; ")
    struct_des['hs443O.fm']  = des("H/O",   "fitted relaxed and O removed; v1; 96 atoms; ")
    struct_des['hs443O.fr']  = des("H/O",   "fitted and relaxed; v1; 97 atoms; ")
    struct_des['hs443O.m']  = des("H/O",   "without imp; v1-16 and v20 for fitted; 96 atoms; ")
    


    struct_des['hs443CO_template']   = des("H/CO",   "template without impurities, used for charge density calc; v1; 96 atoms; ")
    struct_des['hs443CO']   = des("H/CO",   "obtained from template by adding carbon and oxygen; v1-15 - different topo configurations; 98 atoms; ")
    struct_des['hs443CO.m']   = des("H/CO",   "relaxed hs443CO but without impurities; v1-15  96 atoms; ")
    
    struct_des['hs443CC']   = des("H/CC",   "inherited from relaxed hs443CO by replacing; v1-15 - different topo configurations; 98 atoms; ")
    struct_des['hs443OO']   = des("H/OO",   "inherited from relaxed hs443CO by replacing; v1-15 - different topo configurations; 98 atoms; ")

    struct_des['hs443CO2']   = des("H/CO2",   "obtained from template by adding carbon and two oxygen; v1-15 - different topo configurations; 99 atoms; ")
    struct_des['hs443C2O2']   = des("H/C2O2",   "obtained from template by adding two carbon and two oxygen; v1-42 - different topo configurations; 100 atoms; ")

    struct_des['hs554']     = des("H",   "hex cell with bulk Ti hcp; v20; 200 atoms; ")
    struct_des['hs554.f']   = des("H",   "v1 hex cell with bulk Ti hcp; fitted from v1-16; 200 atoms; ")


    struct_des['hs554C']     = des("H/C",   "hex cell plus C; vol with calc_ac; v1; 201 atoms; ")
    struct_des['hs554O']     = des("H/O",   "hex cell plus O; vol with calc_ac; v1; 201 atoms; ")
    struct_des['hs554C.f']  = des("H/C",   "obtained from template by adding carbon; v1; 201 atoms; ")
    struct_des['hs554O.f']  = des("H/O",   "obtained from template by adding carbon; v1; 201 atoms; ")


    return


def red_prec(value, precision = 100.):
    a = value * precision
    return round(a)/1./precision




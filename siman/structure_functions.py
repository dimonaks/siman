#!/usr/bin/env python3
""" 
Include:
1. runBash(cmd)
2. CalcResults
3. interstitial()
4. out_for_paper()
5. shift_analys(st)
6. write_geo(st)
"""

import subprocess
import optparse
import re
import glob
import os
import math
import sys
import colorsys
from pylab import *
from scipy.optimize import leastsq
from sympy import solve, diff, sqrt, Matrix



def inherit_icalc_isotropic(new_structure = '', start_new_version = None,  base_calculation = (None, None, None), database = None, min_mult = 1, max_mult = 1, num_points = 2, geo_folder = ''):
    from calc_manage import inherit_icalc
    min_mult = min_mult
    max_mult = max_mult
    num_points = num_points
    step = (max_mult - min_mult)/(num_points - 1)
    mult_list = [min_mult+step*i for i in range(num_points)]
    version = start_new_version
    for j in mult_list:
        inherit_icalc('isotropic',   new_structure,  version, base_calculation, database, mult_rprimd = j, geo_folder=geo_folder)
        version += 1

def inherit_icalc_c_a(new_structure = '', start_new_version = None,  base_calculation = (None, None, None), database = None, min_mult_a = 1, max_mult_a = 1, num_points_a = 2, min_mult_c = 1, max_mult_c = 1,num_points_c = 2, geo_folder=''):
    from classes import inherit_icalc
    
    if num_points_a > 1:    
        # Lattice parameter a
        min_mult_a = min_mult_a
        max_mult_a = max_mult_a
        num_points_a = num_points_a
        step_a = (max_mult_a - min_mult_a)/(num_points_a - 1)
        mult_list_a = [min_mult_a+step_a*i for i in range(num_points_a)]

    if num_points_c > 1:  
        # Lattice parameter c
        min_mult_c = min_mult_c
        max_mult_c = max_mult_c
        num_points_c = num_points_c
        step_c = (max_mult_c - min_mult_c)/(num_points_c - 1)
        mult_list_c = [min_mult_c+step_c*i for i in range(num_points_c)]

    print('database', database)

    version = start_new_version

    if num_points_a > 1 and num_points_c > 1:    
        for j in mult_list_a:
            for k in mult_list_c:
                inherit_icalc('c_a',   new_structure,  version, base_calculation, database, mult_a = j, mult_c = k, geo_folder=geo_folder)
                version += 1

    elif num_points_c == 1:
        for j in mult_list_a:
            inherit_icalc('c_a',   new_structure,  version, base_calculation, database, mult_a = j, mult_c = 1, geo_folder=geo_folder, override=override)
            version += 1

    elif num_points_a == 1:
        for j in mult_list_c:
            inherit_icalc('c_a',   new_structure,  version, base_calculation, database, mult_a = 1, mult_c = j, geo_folder=geo_folder, override=override)
            version += 1
            
def inherit_icalc_x_y(new_structure = '', start_new_version = None,  base_calculation = (None, None, None), database = None, min_mult_a = 1, max_mult_a = 1, num_points_a = 2, min_mult_b = 1, max_mult_b = 1,num_points_b = 2, geo_folder=''):
    from calc_manage import inherit_icalc
    
    if num_points_a > 1:
        # Coordinate x in rprimd
        step_a = (max_mult_a - min_mult_a)/(num_points_a - 1)
        mult_list_a = [min_mult_a+step_a*i for i in range(num_points_a)]

    if num_points_b > 1:
        # Coordinate y in rprimd
        step_b = (max_mult_b - min_mult_b)/(num_points_b - 1)
        mult_list_b = [min_mult_b+step_b*i for i in range(num_points_b)]

    print('database', database)

    version = start_new_version

    if num_points_a > 1 and num_points_b > 1:
        for j in mult_list_a:
            for k in mult_list_b:
                inherit_icalc('xy',   new_structure,  version, base_calculation, database, mult_a = j, mult_b = k, geo_folder=geo_folder)
                version += 1
    elif num_points_b == 1:
        for j in mult_list_a:
            inherit_icalc('xy',   new_structure,  version, base_calculation, database, mult_a = j, mult_b = 1, geo_folder=geo_folder)
            version += 1

    elif num_points_a == 1:
        for j in mult_list_b:
            inherit_icalc('xy',   new_structure,  version, base_calculation, database, mult_a = 1, mult_b = j, geo_folder=geo_folder)
            version += 1
    



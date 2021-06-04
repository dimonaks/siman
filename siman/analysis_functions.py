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


def min_distance(database = None, calculation = ()):
    c = database[calculation] 
    min_dist = 1000000
    for i in range(c.natom):
        for j in range(c.natom):
            if j == i: continue
            at1 = c.xcart[i]
            at2 = c.xcart[j]
            dist = ((at1[0]-at2[0])**2 + (at1[1]-at2[1])**2 + (at1[2]-at2[2])**2)**0.5
            if min_dist>dist: 
                min_dist = dist
                atom1 = at1
                atom2 = at2
    print('Minimal distance = ', min_dist)
    print('Atom 1 = ', atom1)
    print('Atom 2 = ', atom2)

def formation_energy(database = None, calc_def = (), calc_id = ()):
    defect = database[calc_def]
    ideal = database[calc_id]
    n_at_def = defect.natom
    e_def = defect.energy_free
    e_id_at = ideal.energy_free/ideal.natom
    E_f = e_def - n_at_def*e_id_at
    print('Formation energy for defect '+calc_def[0]+' = '+str(E_f)+' eV')

# Fitting of the E(a,c) dependence for the equilibrium c/a searching

import xalglib

class ALGLIB:

    def build_2d_bicubic_spline(self, x, m, y, n, z, d): self.bicubicv2d = xalglib.spline2dbuildbicubicv(x, m, y, n, z, d)

    def calc(self, x, y, ind): 
        l = xalglib.spline2ddiff(self.bicubicv2d,x,y)
        if ind==0: return l[0]    # z
        elif ind==1: return l[1]  # dz/dx
        elif ind==2: return l[2]  # dz/dy
        elif ind==3: return l[3]  # d2z/dxdy
        else: raise RuntimeError ('Unknown ind = '+str(ind))



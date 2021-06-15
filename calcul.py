# -*- coding: utf-8 -*- 
#Copyright Aksyonov D.A
"""
module contains different calculators of energy, etc.
"""
from __future__ import division, unicode_literals, absolute_import 
import os, copy, shutil, sys
import numpy as np

from siman.geo import image_distance

def buck(r, el1, el2):
    #Repulsive energy according to Kunz1992
    A = {'Li':1.061, 'Na':1.432, 'K':1.794, 'Fe':1.335, 'V':1.276, 'Mn':1.366, 'O':1.853, 'F':1.492, 'P':0.887} # +2 oxidation states for metals (small diff from +3)
    B = {'Li':0.07, 'Na':0.082, 'K':0.11, 'Fe':0.066, 'V':0.086, 'Mn':0.066, 'O':0.168, 'F':0.133, 'P':0.017}

    kcal_mol2eV = 4.3363e-2 
    if r > 0.001:
        E = kcal_mol2eV* (B[el1]+B[el2]) * np.exp( (A[el1] + A[el2] - r)/ (B[el1] +B[el2])   )
    else:
        E = 0
    return E

def site_repulsive_e(st, i):
    """
    calculate repulsive energy for site i of structure st
    """

    xc = st.xcart[i]
    el = st.get_elements()
    E = 0
    for j, x in enumerate(st.xcart):
        d = image_distance(xc, x, st.rprimd)[0]
        E += buck(d, el[i], el[j])

    return E
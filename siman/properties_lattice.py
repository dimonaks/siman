#!/usr/bin/env python3
""" 
Author: Kartamyshev A.I. (Darth Feiwante)
"""

def min_distance(database = None, calculation = ()):
    """
    This function is used to find minimal distance in the lattice presented 
    as a part of the 'Calculation' object

    INPUT:
        - database (.gbdm3) - database dictionary; could be provided; if not then are taken from header
        - calculation (tuple) - tuple describing the Calculation object in form ('structure', 'set', version)
    RETURN:
        None
    SOURCE:
        None
    TODO:
        Some improvements
    """
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


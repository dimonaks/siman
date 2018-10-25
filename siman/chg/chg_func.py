#!/usr/bin/env python
import sys, os
import shutil
import numpy as np


from ase.calculators.vasp import VaspChargeDensity
from siman.chg.vasputil_chgarith_module import chgarith

from siman.header import runBash, printlog


def chg_at_point(chgfile, xred1, ):
    """
    Return the the value of charge density at coordinate xred1; Actually it provides charge density for the closest grid point
    Most probably the units are (el/A^3)

    chgfile - full path to the file with charge density
    xred1 - reduced coordinate;

    RETURN: 
    Charge density at given point
    """
    vasp_charge = VaspChargeDensity(chgfile)
    density = vasp_charge.chg[-1]
    atoms = vasp_charge.atoms[-1]
    del vasp_charge
    # print density[0][0][0]

    ngridpts = numpy.array(density.shape) # size of grid
    print ('Size of grid', ngridpts)
    # rprimd = atoms.get_cell()
    # print rprimd

    # xred1 = [0.5, 0.5, 0.5]
    # rprimd_lengths=numpy.sqrt(numpy.dot(rprimd,rprimd.transpose()).diagonal()) #length of cell vectors
    i,j,k =  [ int(round(x * (n-1) ) ) for x, n in zip(xred1, ngridpts)]# corresponding to xred1 point
    print (i,j,k)
    print ('Density at xred', xred1, 'is',  density[i][j][k])
    return density[i][j][k]



def cal_chg_diff(cl1, cl2, wcell, chg = 'CHGCAR'):
    """1. Calculate differences of charge densities
    Works on local computer
    wcell = 0 or 1 - which cell to use to show
    chg (str) - which file to use 
        CHGCAR - the name as outcar
            if not exist CHG is used
        PARCHG - partial charge, the name without any additions
    

    TO DO:
    instead of paths to files, work with objects
    d = d(cl1) - d(cl2)
    d is calculated on server


    """
    files = []
    for cl in cl1, cl2:
        if chg == 'CHGCAR':
            file = cl.get_chg_file(nametype  = 'asoutcar')
            if not file:
                printlog('No CHGCAR for cl',cl.id[0], 'trying CHG', imp = 'Y')
                file = cl.get_chg_file('CHG', nametype  = 'asoutcar')
        elif chg == 'PARCHG':
            file = cl.get_file('PARCHG')


        files.append(file)



    file1 = files[0]
    file2 = files[1]


    if file1 == None or file2 == None:
        printlog('Error!, chg not found for one cl:', file1, file2)


    working_dir = cl1.dir

    dendiff_filename = working_dir + (chg+'_'+str(cl1.id[0])+'-'+str(cl2.id[0])).replace('.', '_')

    printlog('Diff =', file1, '-', file2)
    chgarith(file1, file2, '-', dendiff_filename, wcell)

    printlog('Charge difference saved to', dendiff_filename, imp = 'Y')

    return dendiff_filename
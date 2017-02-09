#!/usr/bin/env python
import sys, os
import shutil
import numpy as np


from ase.calculators.vasp import VaspChargeDensity
from chg.vasputil_chgarith_module import chgarith
from header import runBash, printlog


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



def cal_chg_diff(cl1, cl2, wcell):
    """1. Calculate differences of charge densities
    Works on local computer

    TO DO:
    instead of paths to files, work with objects
    d = d(cl1) - d(cl2)
    d is calculated on server


    """

    file1 = cl1.get_chg_file()
    file2 = cl2.get_chg_file()

    working_dir = cl1.dir

    dendiff_filename = working_dir + ('CHGCAR_'+str(cl1.id[0])+'-'+str(cl2.id[0])).replace('.', '_')

    chgarith(file1, file2, '-', dendiff_filename, wcell)
    printlog('Charge difference saved to', dendiff_filename, imp = 'Y')

    return dendiff_filename
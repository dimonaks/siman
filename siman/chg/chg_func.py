#!/usr/bin/env python
import sys, os
import shutil
import numpy as np

try:
    from ase.calculators.vasp import VaspChargeDensity
except:
    print('No ase')
from siman.header import runBash, printlog
from siman.chg.vasputil_chgarith_module import chgarith


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



def cal_chg_diff(cl1, cl2, cl3=None, wcell=0, chg = 'CHGCAR'):
    """1. Calculate differences of charge densities
    Works on local computer
    wcell = 0 or 1 - which cell to use to show atom (cl1 or cl2)
    chg (str) - which file to use 
        CHGCAR - the name as outcar
            if not exist CHG is used
        PARCHG - partial charge, the name without any additions
    

    TO DO:
    instead of paths to files, work with objects
    if cl3 is None
        d = d(cl1) - d(cl2)
    else
        d = d(cl1) - [d(cl2)+d(cl3)]

    d is calculated on server


    """
    files = []
    for cl in cl1, cl2, cl3:
        if cl is None:
            continue
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
    if cl3:
        printlog('cl3 is detected Diff =', dendiff_filename, '-', files[2], imp = 'Y')

        chgarith(dendiff_filename, files[2], '-', dendiff_filename+'cl3', wcell=0)
        printlog('Charge difference saved to', dendiff_filename+'cl3', imp = 'Y')


    return dendiff_filename


def cal_chg_diff_files(chg_file1, chg_file2, cl3=None, wcell=0, chg = 'CHGCAR'):
    """1. Calculate differences of charge densities
    Works on local computer
    wcell = 0 or 1 - which cell to use to show atom (cl1 or cl2)
    chg (str) - which file to use 
        CHGCAR - the name as outcar
            if not exist CHG is used
        PARCHG - partial charge, the name without any additions
    

    TO DO:
    instead of paths to files, work with objects
    if cl3 is None
        d = d(cl1) - d(cl2)
    else
        d = d(cl1) - [d(cl2)+d(cl3)]

    d is calculated on server


    """
    files = []
    for cl in cl1, cl2, cl3:
        if cl is None:
            continue
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
    if cl3:
        printlog('cl3 is detected Diff =', dendiff_filename, '-', files[2], imp = 'Y')

        chgarith(dendiff_filename, files[2], '-', dendiff_filename+'cl3', wcell=0)
        printlog('Charge difference saved to', dendiff_filename+'cl3', imp = 'Y')


    return dendiff_filename


def chg_at_z_direct(cl, k_p = 20, plot = None, filetype = 'CHGCAR'):
    """
    Return the the value of charge density or electrostatic potential along z direction of slab; 

    chgfile - full path to the file with charge density
    cl - Calculation() with slab structure; 
    it is needed for the definition of the correct coordinates of points in a structure  in which  will be calculated of el/stat pot 

    RETURN: 
    List of z-coordinates and respective average value of electrostatic pot in the z slice.
    """
    from siman.picture_functions import fit_and_plot
    

    if filetype == 'CHGCAR':
        chgfile = cl.get_chg_file()       
    else: 
        chgfile = cl.get_file(filetype = filetype)

    st = cl.end
    vasp_charge = VaspChargeDensity(chgfile)
    density = vasp_charge.chg[-1]
    atoms = vasp_charge.atoms[-1]
    del vasp_charge


    ngridpts = np.array(density.shape) # size of grid
    # print ('Size of grid', ngridpts)

    z = int(cl.vlength[2]*10)

    elst = []
    z_coord = []
    xred1 = [0,0,0]
    for n3 in range(0,z):
        dens = 0
        # for more accurate calculation of electrostatic pot, it is needed to split the z-sliced plane into a grid of points k_p (20 x 20) 
        # and find the average value for the slice
        for n1 in range(0,k_p):
            for n2 in range(0,k_p):
                xred1[0] = n1/k_p
                xred1[1] = n2/k_p
                xred1[2] = n3/z
                i,j,k =  [ int(round(x * (n-1) ) ) for x, n in zip(xred1, ngridpts)]# corresponding to xred1 point
                dens += density[i][j][k]*st.vol

        elst.append(dens/k_p**2 * st.vol)
        z_coord.append(n3/10)


    if plot:
        # print(z_coord, elst)
        fit_and_plot(a=(z_coord, elst, '-b'), xlabel = 'Z coordinate, $\AA$', 
            ylabel = 'Potential, eV', filename = 'figs/'+st.id[0]+'_pot')

    return z_coord, elst
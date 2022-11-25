#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (c) 2008, 2009, 2010 Janne Blomqvist

# This source code file is subject to the terms of the MIT (Expat)
# License. See the file LICENSE for details.

"""Simple arithmetic operations (+,-,*,/) on CHG/CHGCAR files"""

from optparse import OptionParser
import sys
try:
    sys.path.append('/home/aksenov/Simulation_wrapper/ase')
    from ase.calculators.vasp import VaspChargeDensity
    import ase.calculators.vasp
    # print(ase.calculators.vasp.__file__) 
except:
    print('No ase')


usage = """%prog [options] chgfile0 op chgfile1

chgfile0 and chgfile1 are the charge density files in VASP CHG or
CHGCAR format. If the charge density files are not compatible, you
get what you deserve.

op is an elementary arithmetic operator (+,-,*,/), or avg to calculate
the average."""

# parser = OptionParser(usage)
# parser.add_option('-o', '--output', dest='outfile', help='Output to file \
# named by this argument. If omitted, defaults to a file chgfile0_op_chgfile1')
# parser.add_option('-w', '--which-cell', dest='wcell', help='Take the \
# embedded supercell information from which charge density file. Must be \
# either 0 or 1.')
# (options, args) = parser.parse_args()

def chgarith(chgf1, chgf2, op, filename, wcell):

    # chgf1 = args[0]
    # chgf2 = args[2]
    # op = args[1]

    chg1 = VaspChargeDensity(chgf1)
    chg2 = VaspChargeDensity(chgf2)

    # if options.wcell:
    #     wcell = int(options.wcell)
    if wcell == 0:
        chga = chg1
    elif wcell == 1:
        chga = chg2
    #     else:
    #         print ('Error, invalid argument to -w option')
    #         sys.exit(1)
    # else:
    # chga = chg1

    if len(chg1.chg) != len(chg2.chg):
        print ('Number of images in charge density files not equal. Using just   the final images in both files.')
        print ('len(chg.chg)', len(chg1.chg), len(chg2.chg) )
        chg1.chg = [chg1.chg[-1]]
        chg1.atoms = [chg1.atoms[-1]]
        if chg1.is_spin_polarized():
            chg1.chgdiff = [chg1.chgdiff[-1]]
            chg2.chgdiff = [chg2.chgdiff[-1]]
        chg2.chg = [chg2.chg[-1]]
        chg2.atoms = [chg2.atoms[-1]]

    newchg = VaspChargeDensity(None)

    print ('Start charge manipul')

    for i, atchg in enumerate(chg1.chg):
        c1 = atchg
        c2 = chg2.chg[i]
        newchg.atoms.append(chga.atoms[i].copy())
        if op == '+':
            nc = c1 + c2
            oplong = '_add_'
        elif op == '-':
            nc = c1 - c2
            oplong = '_sub_'
        elif op == '*':
            nc = c1 * c2
            oplong = '_mult_'
        elif op == '/':
            nc = c1 / c2
            oplong = '_div_'
        elif op == 'avg':
            nc = (c1 + c2) / 2
            oplong = '_avg_'
        newchg.chg.append(nc)


    if chg1.is_spin_polarized():
        print ('Spin polarized')

        for i, cd in enumerate(chg1.chgdiff):
            cd2 = chg2.chgdiff[i]
            if op == '+':
                nd = cd + cd2
            elif op == '-':
                nd = cd - cd2
            elif op == '*':
                nd = cd * cd2
            elif op == '/':
                nd = cd / cd2
            elif op == 'avg':
                nd = (cd + cd2) / 2
            newchg.chgdiff.append(nd)

    # Screw doing anything fancy with the augmentation charges
    # Just take them from the same file as the embedded atoms object.
    newchg.aug = chga.aug
    newchg.augdiff = chga.augdiff

    # if options.outfile:
    #     fname = options.outfile
    # else:
    #     from os.path import basename
    #     fname = basename(chgf1) + oplong + basename(chgf2)

    newchg.write(filename)

    return filename
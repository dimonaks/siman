#!/usr/bin/python3
from __future__ import division, unicode_literals, absolute_import, print_function
import sys, glob, os, json, pickle, re
from pathlib import Path
import pandas as pd
from pymatgen.symmetry.analyzer import PointGroupAnalyzer  

import numpy as np
# sys.path.append('/home/aksenov/Simulation_wrapper/siman') #path to siman library
# sys.path.append('siman') #path to siman library on sd cluster
from siman.classes import CalculationVasp

"""
NAME
    create .cvs files from .pickle for web visualization

SYNOPSIS
    make_database_cvs.py path_to_database pickle_type

    pickle_type:
    1 - siman vasp
    2 - pymatgen gaussian


DESCRIPTION

    Walk through path_to_database, read objects from pickle file and create
    required cvs files.

REQUIRE
        numpy, siman

AUTHOR
        Aksyonov Dmitry, Skoltech, Moscow

"""

__author__ = "Aksyonov Dmitry"
__copyright__ = "Copyright 2017, Skoltech"
__version__ = "0.1"
__maintainer__ = "Aksyonov Dmitry"
__email__ = "dimonaks@gmail.com"
__status__ = "alpha"
__date__ = ""


def read_siman_vasp(cl, row, filename):



    c = filename.split('_')[0]
    # print(c) 
    c = re.sub('\D', '', c) # replace any non-digit character with ''
    if c:
        c = int(c)
        if c == 1:
            c *= 100
    else:
        c = None
    row['conc'] = c


    row['name'] = cl.name
    row['energy'] = '{:.3f}'.format(cl.energy_sigma0)
    row['energy_at'] = '{:.3f}'.format(cl.energy_sigma0/cl.end.natom)
    row['status'] = 'OK'
    row['spg'] = cl.end.get_space_group_info()[0]
    row['nat'] = cl.end.natom
    row['vol'] = cl.end.vol

    return row

def read_pymatgen_gaussian(cl, row, filename):

    # print('cl ', type(cl), dir(cl))
    try:
        st = cl.final_structure # Molecule

        # print('st ', type(st), dir(st))
        pg = PointGroupAnalyzer(st)
        row['spg'] = pg.sch_symbol

        nat = st.num_sites
        row['name'] = st.formula
        row['energy'] = '{:.3f}'.format(cl.final_energy)
        row['energy_at'] = '{:.3f}'.format(cl.final_energy/nat)
        row['status'] = 'OK'
        row['nat'] = nat
        print(row['name'], row['energy'], row['spg'], row['path'])
        return row

    except:
        st = None
        return None





path2database = sys.argv[1]

if len(sys.argv) > 2:
    pickle_type   = int(sys.argv[2])
else:
    pickle_type = 1

print(pickle_type)

if pickle_type == 1:

    columns = ['name', 'spg', 'conc', 'energy_at', 'nat', 'vol',                      'path', 'status']
    col_rename = {'name':'Name', 'spg':'Group', 'conc':'X, %', 'energy_at':'E, eV/at', 'nat':'Nsites', 'vol':'Vol (A3)',   'path':'path', 'status':'status'}
    read_func = read_siman_vasp

elif pickle_type == 2:

    columns = ['Formula', 'spg', 'energy_at', 'nat', 'path', 'status']
    col_rename = {'name':'Name', 'spg':'Group', 'energy_at':'E, eV/at', 'nat':'Nsites', 'path':'path', 'status':'status'}

    read_func = read_pymatgen_gaussian



os.chdir(path2database)


csv = []

if 1:

    with open('hash_dict.json') as fp:
       hash_dict =  json.load(fp)
    # print(hash_dict)

    for pickle_file in hash_dict:
        # pickle_fileP = Path(pickle_file+'.pickle')
        pickle_fileP = Path(pickle_file)
        row = {'path':pickle_file}

        if pickle_fileP.is_file():
            with pickle_fileP.open('rb') as f:
                cl = pickle.load(f)
            # a = cl.read_results()
            # print(a    )
            # try:
            
            row = read_func(cl, row, pickle_fileP.name)

            if row is None:
                print('Broken file:', pickle_file)
                continue

        else:
            row['status'] = 'no pickle file'

        csv.append(row)

    df = pd.DataFrame(csv, columns = columns)

    # df['nat'] = df['nat'].apply(np.int64)


    df = df.rename(columns = col_rename)

    df.to_csv('database1.csv')

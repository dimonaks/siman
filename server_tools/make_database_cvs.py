#!/usr/bin/python3
from __future__ import division, unicode_literals, absolute_import, print_function
import sys, glob, os, json, pickle, re
from pathlib import Path
import numpy as np
sys.path.append('/home/aksenov/Simulation_wrapper/siman') #path to siman library
sys.path.append('siman') #path to siman library on sd cluster
from classes import CalculationVasp

"""
NAME
    create .cvs files from .pickle for web visualization

SYNOPSIS
    make_database_cvs.py path_to_database

DESCRIPTION

    Walk through path_to_database, read pickle in Calculation objects and create
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


path2database = sys.argv[1]


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
            try:
                row['name'] = cl.name
                row['energy'] = '{:.3f}'.format(cl.energy_sigma0)
                row['energy_at'] = '{:.3f}'.format(cl.energy_sigma0/cl.end.natom)
                row['status'] = 'OK'
                row['spg'] = cl.end.get_space_group_info()[0]
                row['nat'] = cl.end.natom
                row['vol'] = cl.end.vol
            except:
                row['status'] = 'broken pickle file'

            # print (pickle_fileP.name.split('_')[0])
            c = pickle_fileP.name.split('_')[0]
            print(c) 
            c = re.sub('\D', '', c) # replace any non-digit character with ''
            if c:
                c = int(c)
                if c == 1:
                    c *= 100
            else:
                c = None
            row['conc'] = c


        else:
            row['status'] = 'no pickle file'

        csv.append(row)
    import pandas as pd
    df = pd.DataFrame(csv, columns = ['name', 'spg', 'conc', 'energy_at', 'nat', 'vol',                      'path', 'status'])

    # df['nat'] = df['nat'].apply(np.int64)


    df = df.rename(columns = {'name':'Name', 'spg':'Group', 'conc':'X, %', 'energy_at':'E, eV/at', 'nat':'Nsites', 'vol':'Vol (A3)',   'path':'path', 'status':'status'})
    # print(df)
    df.to_csv('database1.csv')

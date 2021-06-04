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






def plot_graph(graph_folder = '', name = '', axis1 = '', axis2 = '', args=[]):
    import matplotlib.pyplot as plt
    import shutil as S
    import os

    # Make new directory
    try: 
        os.mkdir(os.getcwd()+'/Graphs/'+graph_folder)
    except OSError: pass

    # Removing old graph
    try: 
        os.remove(os.getcwd()+'/Graphs/'+graph_folder+'/'+name)
    except OSError: 
        pass
        
    # Plot new graph
    d_color = dict(zip([0,1,2,3,4,5,6],['r','b','g','c', 'm', 'k', 'y']))
    d_marker = dict(zip([0,1,2,3,4,5,6,7,8,9,10,11],['^','s','o','D','h','H','8','p','*','v', '>','<']))
    fig = plt.figure(1)
    count = 0
    for g in args:
        plt.plot(g[0],  g[1], linewidth=2, linestyle='-', marker = d_marker[count], markersize=3, color=d_color[count], label = g[2]+'_'+g[3])
        count += 1
    plt.xlabel(axis1, fontsize = 22)
    plt.ylabel(axis2, fontsize = 22)
    plt.legend(bbox_to_anchor=(0.9, 0.92), borderaxespad=0., labelspacing=0.3, numpoints=1, frameon=True, markerscale=1., handletextpad=0.3)
    fig.savefig(name, format="pdf")
    S.move(name, os.getcwd()+'/Graphs/'+graph_folder)




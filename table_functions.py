# -*- coding: utf-8 -*- 
from __future__ import division, unicode_literals, absolute_import 
import sys, os

# import header
# from operator import itemgetter
# from classes import res_loop , add_loop
# from pairs import 
# from functions import image_distance, local_surrounding
# from chargeden.functions import chg_at_point, cal_chg_diff
# from dos.functions import plot_dos

# from ase.utils.eos import EquationOfState
import scipy
from scipy import interpolate
from scipy.interpolate import spline 
# print (scipy.__version__)
# print (dir(interpolate))
try:
    from scipy.interpolate import  CubicSpline
except:
    print('scipy.interpolate.CubicSpline is not avail')

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import header
from header import print_and_log, printlog
from header import calc
from inout import write_xyz
from small_functions import makedir, is_list_like

from geo import replic





def latex_table(table, caption, label, header = None, fullpage = '', filename = None, writetype = 'w', header0 = None, size = None,
    replace = None ):
    """
    If header is not provided, table[0] is used as a header

    header0 - additional header0 befor main header for complex tables
    
    path_to_paper should be provided

    replace - list of tuples for replacements

    """
    def myprint(string):
        if filename:
            f.write(string+"\n")
            print( string)
        else:
            print( string)


    if filename:
        # path = path_to_paper+'/tab/'
        path = ''
        f = open(path+filename, writetype)
        print_and_log("Saving table to "+path+filename+'\n')

    for i in range(len(table)):
        if is_list_like(table[i]):
            tab = ' & '.join([str(l) for l in table[i]])
            table[i] = tab


    n = len(table[0].split('&'))-2
    print( 'Number of columns = ', n + 2)
    
    myprint('\\begin{table'+fullpage+'}')
    myprint('\\center')
    if size: myprint('\\'+size)
        


    myprint('\\caption{'+caption+'}')
    myprint('\\label{'+label+'}')

    myprint('\\begin{tabular}{l'+ n*'c'+'r}')
    myprint('\\hline')

    if header0:
        myprint(header0+'\\\\')
        myprint('\\hline')

    if header:
        myprint(header+'\\\\')
        tabbeg = 0
    else:
        myprint(table[0]+' \\\\')
        tabbeg = 1



    myprint('\\hline')
    for r in table[tabbeg:] :
        if '&-' in r:
            r = r.replace('-','--')
        else:
            r = r.replace(' -','--') #to save beautiful columns 
        r+=' '
        if '-- ' in r:
            r = r.replace('-- ',' - ')
        
        for rep in replace:
            # if rep[0] in r:

            r = r.replace(*rep)



        if 'hline' in r: 
            myprint(r)
        else:
            myprint(r + '\\\\')




    myprint('\\hline')
    myprint('\\end{tabular}')
    myprint('\\end{table'+fullpage+'}')

    if filename:
        f.close()
    return

def geo_table_row(cl, name = '', show_alpha = 0):
    #Basic table
    
    st = cl.end
    # if not name:
    spg = st.get_space_group_info()
    spg = spg[0].replace('_1', '$_1$')

    alpha, beta, gamma = st.get_angles()
    elem = np.array(cl.end.get_elements())

    alpha = '& {:5.1f}'.format(alpha)
    if not show_alpha:
        alpha = ''


    return '{:15s} &{:s} & {:5.2f} & {:5.2f} & {:5.2f} '.format(name[0], 'DFT',  cl.vlength[0], 
        cl.vlength[1], cl.vlength[2])+alpha+'& {:5.1f} & {:s}'.format(st.vol, spg)









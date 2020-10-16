# -*- coding: utf-8 -*- 
from __future__ import division, unicode_literals, absolute_import 
import sys, os, re

# import header
# from operator import itemgetter
# from classes import res_loop , add_loop
# from pairs import 
# from functions import image_distance, local_surrounding
# from chargeden.functions import chg_at_point, cal_chg_diff
# from dos.functions import plot_dos

# from ase.utils.eos import EquationOfState

# print (scipy.__version__)
# print (dir(interpolate))
try:
    import scipy
    from scipy import interpolate
    # print (scipy.__version__)
    # print (dir(interpolate))
except:
    print('table_functions.py: scipy is not avail')
try:
    # from scipy.interpolate import spline 

    from scipy.interpolate import  CubicSpline
except:
    print('table_functions.py: scipy.interpolate.CubicSpline is not avail')

import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

from siman import header
from siman.header import print_and_log, printlog, calc, runBash
from siman.inout import write_xyz
from siman.small_functions import makedir, is_list_like, latex_spg, latex_chem, get_common_chemical_base
from siman.geo import replic
from siman.analysis import calc_redox
if header.pymatgen_flag:
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer



def latex_table(table, caption, label, header = None, fullpage = '', filename = None, writetype = 'w', header0 = None, size = None,
    replace = None, float_format = None, tab_type = 'tabular', width = 0.75):
    """
    If header is not provided, table[0] is used as a header

    header0 - additional header0 befor main header for complex tables
    
    path_to_paper should be provided

    replace - list of tuples for replacements

    float_format - list of float numbers

    tab_type (str) - see latex types of tables
        available: 
        'tabular'
        'tabularx'
        'ruled' - phys rev format
    width (float) - in units of textwidth


    """

    table_string = ''

    def myprint(string):
        nonlocal table_string
        if filename:
            f.write(string+"\n")
            print( string)
        else:
            print( string)
        table_string+=string+'\n'


    if filename:
        # path = path_to_paper+'/tab/'
        path = ''
        f = open(path+filename, writetype)
        print_and_log("Saving table to "+path+filename+'\n')

    for i in range(len(table)):
        if float_format:
            formatter = iter(float_format)
        else:
            formatter = (2 for i in range(100))

        if is_list_like(table[i]):
            tab = ''
            for j, l in enumerate(table[i]):
                # print(type(l))
                if type(l) != str:
                    fmt = '{:3.'+str(next(formatter))+'f}'
                    # fmt = 'a'
                    # print(fmt)
                    pos = fmt.format(l)
                    # 
                    # fmt
                    # print(pos)
                    # pos = str(l)
                else:
                    pos = str(l)
                tab+=pos + " & "
            # tab = ' & '.join([str(l) for l in table[i]])
            table[i] = tab[0:-3]


    n = len(table[0].split('&'))-2
    print( 'Number of columns = ', n + 2)
    
    myprint('\\begin{table'+fullpage+'}')
    myprint('\\center')
    if size: myprint('\\'+size)
        


    myprint('\\caption{'+caption+'}')
    myprint('\\label{'+label+'}')

    if tab_type == 'tabular':
        # tabular = 
        myprint('\\begin{tabular}{l'+ n*'c'+'r}')
        myprint('\\hline')
    
    elif tab_type == 'tabularx':
        myprint('\\begin{tabularx}{'+str(width)+'\\textwidth}{X'+ n*'X'+'X}')
        myprint('\\hline')
    elif tab_type == 'ruled':
        myprint('\\begin{ruledtabular}')
        myprint('\\begin{tabular}{l'+ n*'c'+'r}')


    else:
        printlog('Error! Unknown type of tabular env!')


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
        # print(r)
        if '&-' in r:
            r = r.replace('-','--')
        else:
            r = r.replace(' -','--') #to save beautiful columns 
        r+=' '
        if '-- ' in r:
            r = r.replace('-- ',' - ')
        
        if '_' in r:
            r = r.replace('_','\\_')

        if replace:
            for rep in replace:
                # if rep[0] in r:

                r = r.replace(*rep)



        if 'hline' in r: 
            myprint(r)
        else:
            myprint(r + '\\\\')



    if 'ruled' in tab_type:
        myprint('\\end{tabular}')
        myprint('\\end{ruledtabular}')
    else:
        myprint('\\hline')
        myprint('\\end{'+tab_type+'}')
    
    myprint('\\end{table'+fullpage+'}')

    if filename:
        f.close()
    return table_string

def geo_table_row(cl = None, st = None, name = '', show_angle = 0, mnpo4_hack = False, param_order = None):
    #return list of geo data for cl, which can be used for table
    """
    mnpo4_hack (bool) - if true exchange a and c for mnpo4 phase

    param_order - default [0,1,2]
    """
    po = param_order
    if po is None:
        po = [0,1,2]

    if cl:
        st = cl.end
    # if not name:
    if header.pymatgen_flag:
        spg = st.get_space_group_info()
        spg = latex_spg(spg[0])
       
        #transform to standard
        st_mp = st.convert2pymatgen()
        symprec = 0.1
        sf = SpacegroupAnalyzer(st_mp, symprec = symprec)

        st_mp_prim = sf.find_primitive()
        st_mp_conv = sf.get_conventional_standard_structure()
    else:
        spg = '-'


    # st_mp_prim.
    # print('primitive,', st_mp_prim.lattice)
    # print('conventio,', st_mp_conv.lattice)


    # st_mp_prim = sf.get_primitive_standard_structure()
    # st_mp_prim = sf.get_conventional_standard_structure()


    if not name:
       # print(dir(st_mp ))
       # st.printme()
       name = st.get_reduced_formula()



    name = latex_chem(name)




    alpha, beta, gamma = st.get_angles()
    elem = np.array(st.get_elements())

    if 'a' in show_angle:
        angle = '& {:5.2f}'.format(alpha)
    elif 'b' in show_angle:
        angle = '& {:5.2f}'.format(beta)
    elif 'g' in show_angle:
        angle = '& {:5.2f}'.format(gamma)
    else:
        angle = ''

    v = st.vlength
    a, b, c = v
    if mnpo4_hack and 'MnPO' in name:
        c, b, a = v

    ps = [a,b,c]

    p = []
    # print(ps)
    # print(po)
    for o in po:
        p.append(ps[o])
    # print(p)
    # sys.exit()

    return '{:15s} &{:s} & {:5.2f} & {:5.2f} & {:5.2f} '.format(name, 'DFT+U',  p[0], 
        p[1], p[2])+angle+'& {:5.1f} & {:s}'.format(st.vol, spg)




def table_geometry(st_list = None, cl_list = None, show_angle = None, exp = None, param_order = None):
    """
    Produce standart table with lattice constants
    # print(row)
    exp (list) - list of strings with exp data with '& & &' format
    """
    if st_list is None:
        st_list = [cl.end for cl in cl_list]

    rows = []
    for st in st_list:
        # st.printme()
        row = geo_table_row(st = st, show_angle = show_angle, param_order = param_order)
        rows.append(row)
    if exp:
        rows.extend(exp)

    if 'a' in show_angle:
        angle = "$\\alpha$"
    elif 'b' in show_angle:
        angle = "$\\beta$"
    elif 'g' in show_angle:
        angle = "$\\gamma$"
    else:
        angle = ''


    caption = "Lattice parameters (\AA), volume (\AA$^3$), and space group (spg)."
    return latex_table(rows, caption, 'tab:const', 'Structure & src & $a$ & $b$ & $c$ &'+angle+'& Vol. & spg' )

def table_potentials(cl_list):

    cl_b = cl_list[0]
    rows = []
    for cl in cl_list[1:]:
        a = calc_redox(cl_b, cl)
        n_b = cl_b.end.get_name()
        n = cl.end.get_name()

        base = get_common_chemical_base(cl_b.end, cl.end)

        rows.append([n_b.replace(base, 'X')+'/'+n.replace(base, 'X'), a['redox_pot'], a['vol_red']])

    caption = "Redox potential ($U$) and volume change (\\%). X="+base
    return latex_table(rows, caption, 'tab:const', 'Pair & $U$ & $dV$ ', float_format = [1,1] )


def generate_latex_report(text, filename):
    # 

    fn = filename+'.tex'
    makedir(fn)
    head = r"""
\documentclass{article}
% General document formatting
\usepackage[margin=0.7in]{geometry}
\usepackage[parfill]{parskip}
\usepackage[utf8]{inputenc}

% Related to math
\usepackage{amsmath,amssymb,amsfonts,amsthm}

\begin{document}"""
    

    with open(fn, 'w') as f:
        f.write(head+'\n')
        f.write(text+'\n')
        f.write(r'\end{document}')


    runBash('pdflatex '+os.path.basename(fn), cwd = os.path.dirname(fn))
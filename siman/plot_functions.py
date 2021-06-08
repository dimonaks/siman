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



def plot_energy_xy(it,ise,verlist,calc, folder='', type_plot = '3d_contour', labels=(), lab_size=15, tick_size=15, fig_size=(), 
                   fig_title='', xlim = (), ylim = (), zlim = ()):

    acell_list = []
    etotal_list = []    
    for v in verlist:
        id1 = (it,ise,v)
        acell = []
        acell.append(calc[id1].rprimd[0][0]); acell.append(calc[id1].rprimd[1][1]); acell.append(calc[id1].rprimd[2][2])
        acell_list.append( acell )
        etotal_list.append ( calc[id1].energy_sigma0 )


    f0 = open(folder+'/deformation_xy.out', 'w')
  
  
    f0.write('Calc equilibrium acell for\n')

    a=[]; b=[];
    for acell in acell_list:
        if acell[0] not in a: a.append(acell[0])
        if acell[1] not in b: b.append(acell[1])
        
    f0.write('Read number points acell = '+ str(len(a))+'\n')
    f0.write('Read number points ccell = '+ str(len(b))+'\n\n')
    print ('Read number points acell = '+ str(len(a))+'\n')
    print ('Read number points ccell = '+ str(len(b))+'\n\n')

    f0.write('Limits acell read from file = ' + str(min(a))+'  '+ str(max(a)) + '; (max(a)-min(a))/2 = ' + str((max(a)-min(a))/2) +'\n')
    f0.write('Limits ccell read from file = ' + str(min(b))+'  '+ str(max(b)) + '; (max(c)-min(c))/2 = ' + str((max(b)-min(b))/2) +'\n\n')
          
    etot = [None for i in range(len(a)*len(b))]

    for i in range(len(etotal_list)):
        acell = acell_list[i]
        jj=0
        exit = False
        for j in b:
            if exit==True: break
            for k in a:
                if acell[0]==k and acell[1]==j:
                    etot[jj] = etotal_list[i]
                    exit = True
                    break
                else: jj+=1

    # Check etot
    assert None not in etot, 'None in etot'
    f0.close()

    # Chuan bi cac danh sach
    if len(acell_list)!=len(etotal_list): raise RuntimeError ('Не равны длины списков len(acell_list)!=len(etotal_list) для построения 3D рисунка!')
    acell_list1 = [tuple(i) for i in acell_list]
    setca_ab_dict = dict(zip(acell_list1, etotal_list))
    x_setca, y_setca, z_setca = [], [], []
    x1, y1, z1 = [], [], []
    ii = 0
    for key in sorted(setca_ab_dict):
        if ii==0: 
            key_cur = key
            ii=1
        if abs(key[0]-key_cur[0])<10**(-10): 
            x1.append(key[0])
            y1.append(key[1])
            z1.append(setca_ab_dict[key])
        else:
            key_cur = key
            x_setca.append(x1)
            y_setca.append(y1)
            z_setca.append(z1)
            x1 = [key[0]]; y1 = [key[1]]; z1 = [setca_ab_dict[key]];
    else:
        x_setca.append(x1)
        y_setca.append(y1)
        z_setca.append(z1)

    # Tim gia tri it nhat
    min_x_list = [min(i) for i in x_setca]
    min_y_list = [min(i) for i in y_setca]
    min_z_list = [min(i) for i in z_setca]

    min_x = min(min_x_list)
    min_y = min(min_y_list)
    min_z = min(min_z_list)

    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.mplot3d import axes3d
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter, MultipleLocator, AutoLocator

    # Biến thành numpy.array

    x_setca_np = np.array(x_setca)
    y_setca_np = np.array(y_setca)
    z_setca_np = np.array(z_setca)


    # Vẽ bức tranh

    fig = plt.figure()

    three_dim = False
    contour = True

    if type_plot == '3dim':
        # ax = fig.gca(projection='3d')
        ax = fig.add_subplot(111, projection='3d')
        cset = ax.contourf(x_setca_np, y_setca_np, z_setca_np, zdir='z', offset=zlim[0], cmap=cm.rainbow, zorder=0.5)          
        ax.plot_surface(x_setca_np, y_setca_np, z_setca_np, rstride=8, cstride=8, alpha=0.2,  cmap=cm.rainbow, zorder=0.3)
        ax.scatter(x_setca_np, y_setca_np, z_setca_np, marker='o', color='red') 

        cset = ax.contourf(x_setca_np, y_setca_np, z_setca_np, zdir='x', offset=xlim[0], cmap=cm.rainbow )  
        cset = ax.contourf(x_setca_np, y_setca_np, z_setca_np, zdir='y', offset=ylim[1], cmap=cm.rainbow )


        ax.set_xlabel(labels[0], fontsize=lab_size, labelpad=10)
        ax.set_ylabel(labels[1], fontsize=lab_size, labelpad=10)
        ax.set_zlabel(labels[2], fontsize=lab_size, labelpad=10)

        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        ax.set_zlim(zlim[0], zlim[1])

        ax.tick_params(axis = 'both', labelsize=tick_size)
        # ax.xaxis.set_major_locator(LinearLocator(4))
        # ax.yaxis.set_major_locator(LinearLocator(4))
        # ax.zaxis.set_major_locator(LinearLocator(4))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        ax.set_title(fig_title, fontsize=lab_size)

    elif type_plot == 'contourf':
        cs = plt.contourf(x_setca, y_setca, z_setca, locator=LinearLocator(numticks = 50), cmap=cm.rainbow)
        cs.ax.set_xlabel(labels[0], fontsize=lab_size)
        cs.ax.set_ylabel(labels[1], fontsize=lab_size)
        cs.ax.tick_params(axis = 'both', labelsize=tick_size)
        cbar = plt.colorbar(format = FormatStrFormatter('%.2f'))
        cbar.set_label('$'+labels[2]+'$', fontsize = lab_size)
        cbar.ax.tick_params(labelsize=tick_size) 
        cs.set_title(fig_title, fontsize=lab_size)


    # plt.show()
    fig.set_figheight(fig_size[0])
    fig.set_figwidth(fig_size[1])        
    fig.savefig(folder+'/'+it+'_'+ise+'_deform_xy_'+type_plot+'.pdf', format='pdf', dpi=600)
    
    
    
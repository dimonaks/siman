import sys, os
sys.path.append(os.path.dirname(__file__)+'/../alglib_cpython')
# print sys.path
#from read_write_i import ReadWrite as RW
#from plot_3d_i import Plot_3D as Pl3d


import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from siman.analysis_functions import ALGLIB as ALB

def plot_energy_xy(it,ise,verlist,calc, folder='', type_plot = 'contourf', labels=(), lab_size=15, tick_size=15, fig_size=(), 
                   fig_title='', xlim = (), ylim = (), zlim = (), shag_a=0.001,shag_b=0.001,npoint_a=1000,npoint_b=1000):
    """
    This function produces the plot of the energy surface corresponding to the structure deformed along
    the lattice vectors 1 and 2. The result can be presented as a 2D or 3D figure. 
    Also the function has feature finding minimum of the energy surface based on bicubic spline scheme

    INPUT:
        - it (str) - the name of crystal structure
        - ise (str) - name of the parameters' set
        - verslist (list of int) - list of versions of new calculations
        - calc (.gbdm3) - database dictionary; could be provided; if not then are taken from header
        - folder (str) - directory where all the results will be built
        - type_plot (str) - represents kind of plot to build (2D,3D, contour etc.). 
          Possible options are:
            '3D_proj_xyz_fill' - 3D plot with projections on xy, yz and xz planes in form 2D filled contour plots
            '3D_proj_z_contour' - 3D plot with projection on the xy plane only in form of 2D contour plot
            'contourf' - 2D filled contour plot
        - labels (tuple of str) - has form ('label_x', 'label_y', 'label_z') and represents desirable names of axes
        - lab_size (int) - font size of labels
        - tick_size (int) - font size of ticks
        - fig_size (tuple of floats) - has form (height, width) set size of the figure in units of cm
        - fig_title (str) - optional, the title placed above the plot
        - xlim (tuple of floats) - has form (min_x, max_x) and represents limits of the plot along the 'x' axis
        - ylim (tuple of floats) - has form (min_y, max_y) and represents limits of the plot along the 'y' axis
        - zlim (tuple of floats) - has form (min_z, max_z) and represents limits of the plot along the 'z' axis
        - shag_a (float) - step of the fine grid of points along the lattice vector 1 of the structure for the bicubic spline
        - shag_b (float) - step of the fine grid of points along the lattice vector 2 of the structure for the bicubic spline
        - npoint_a (int) - number of points of the fine grid of points along the lattice vector 1 of the structure for the bicubic spline
        - npoint_b (int) - number of points of the fine grid of points along the lattice vector 2 of the structure for the bicubic spline

    RETURN:
        None
    SOURCE:
        None
    TODO:
        - Add more types of plots
        - Add plot based on fine bicubic spline grid
    """
    acell_list = []
    etotal_list = []    
    for v in verlist:
        id1 = (it,ise,v)
        acell = []
        acell.append(calc[id1].rprimd[0][0]); acell.append(calc[id1].rprimd[1][1]); acell.append(calc[id1].rprimd[2][2])
        acell_list.append( acell )
        etotal_list.append ( calc[id1].energy_sigma0 )


    f0 = open(folder+'/deformation_'+it+'_'+ise+'_xy.out', 'w')
  
  
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

    # Find position with minimum energy
    etot_min = min(etot)
    f0.write('Etot_min (without spline) = '+ str(etot_min)+' eV\n\n')
    ii=0; exit = False
    for j in b:
        if exit==True: break
        for k in a:
            if etot[ii]==etot_min:
                amin = k
                bmin = j
                exit = True
                break
            else: ii+=1

    f0.write('xcell_min (without spline) = '+ str(amin)+' Angstrom\n')
    f0.write('ycell_min (without spline) = '+ str(bmin)+' Angstrom\n\n\n')

    e2 = ALB()
    e2.build_2d_bicubic_spline(a, len(a), b, len(b), etot, 1)

    e_min = 10**(8)
    for j in range(-npoint_b, npoint_b):
        b_cur = bmin+i*shag_b
        for k in range(-npoint_a, npoint_a):
            a_cur = amin+k*shag_a
            e_cur = e2.calc(a_cur,b_cur,0)
            if e_cur < e_min:
                e_min = e_cur
                a_min = a_cur
                b_min = b_cur
    f0.write('Etot_min (with spline) = '+ str(e_min)+' eV\n\n')                
    f0.write('xcell_min (with spline) = '+ str(a_min)+' Angstrom\n')
    f0.write('ycell_min (with spline) = '+ str(b_min)+' Angstrom\n\n\n')
    f0.write('Found min etot for limitation:\n')
    f0.write('xcell = '+str(amin)+' +/- '+str(npoint_a*shag_a)+'  Angstrom'+' (step = '+str(shag_a)+')'+'\n')
    f0.write('ycell = '+str(bmin)+' +/- '+str(npoint_b*shag_b)+'  Angstrom'+' (step = '+str(shag_b)+')'+'\n\n\n') 
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



    # Chuan bi tep voi ket qua
    x_list = []
    y_list = []
    z_list = []
    for i in range(len(x_setca)):  
        x_list += x_setca[i]
        y_list += y_setca[i]
        z_list += z_setca[i]

    f = open(folder+'/plot_energy_xy_'+it+'_'+ise+'.out','w')
    f.write('# {0:^15s} {1:^15s} {2:^15s}'.format('x, Angstrom', 'y, Angstrom', 'E, eV')+'\n')
    for i in range(len(x_list)):
        f.write('{0:^15.5f} {1:^15.5f} {2:^15.5f}'.format(x_list[i], y_list[i], z_list[i])+'\n')
    f.close()

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

    if not xlim: xlim = (min(x_setca), max(x_setca))
    if not ylim: ylim = (min(y_setca), max(y_setca))
    if not zlim: zlim = (min(z_setca), max(z_setca))
    # Vẽ bức tranh

    fig = plt.figure()
    font = matplotlib.font_manager.FontProperties()
    font.set_size(20)

    if type_plot == '3D_proj_xyz_fill':
        # ax = fig.gca(projection='3d')
        ax = fig.add_subplot(111, projection='3d')
        cset = ax.contourf(x_setca_np, y_setca_np, z_setca_np, zdir='z', offset=zlim[0], cmap=cm.rainbow, zorder=0.3)          
        ax.plot_surface(x_setca_np, y_setca_np, z_setca_np, rstride=8, cstride=8, alpha=0.2,  cmap=cm.rainbow, zorder=0.5)
        ax.scatter(x_setca_np, y_setca_np, z_setca_np, marker='o', color='red') 

        cset = ax.contourf(x_setca_np, y_setca_np, z_setca_np, zdir='x', offset=xlim[0], cmap=cm.rainbow )  
        cset = ax.contourf(x_setca_np, y_setca_np, z_setca_np, zdir='y', offset=ylim[1], cmap=cm.rainbow )


        ax.set_xlabel(labels[0], fontsize=lab_size, labelpad=10)
        ax.set_ylabel(labels[1], fontsize=lab_size, labelpad=10)
        ax.set_zlabel(labels[2], fontsize=lab_size, labelpad=10,rotation=0)

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

    if type_plot == '3D_proj_z_contour':
        # ax = fig.gca(projection='3d')
        ax = fig.add_subplot(111, projection='3d')

        # norm = plt.Normalize(z_setca_np.min(), z_setca_np.max())
        # colors = cm.viridis(norm(z_setca_np))
        # rcount, ccount, _ = colors.shape


        # cset = ax.contour(x_setca_np, y_setca_np, z_setca_np, zdir='z', offset=zlim[0], cmap=cm.rainbow, zorder=0.3)          
        ax.plot_surface(x_setca_np, y_setca_np, z_setca_np, rstride=1, cstride=1, cmap=cm.rainbow, zorder=0.5)
        # surf = ax.plot_surface(x_setca_np, y_setca_np, z_setca_np, rcount=rcount, ccount=ccount, facecolors=colors, shade=False, zorder=0.5)
        # surf.set_facecolor((0,0,0,0))
        ax.scatter(x_setca_np, y_setca_np, z_setca_np, marker='o', color='red') 
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
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        ax.azim = 225
        ax.zaxis.set_rotate_label(False)  # disable automatic rotation
        ax.set_zlabel(labels[2], fontsize=lab_size, rotation=90)
        ax.text2D(0.17,0.73, fig_title, fontproperties=font, transform=ax.transAxes)

    elif type_plot == 'contourf':
        cs = plt.contourf(x_setca, y_setca, z_setca, locator=LinearLocator(numticks = 50), cmap=cm.rainbow)
        cs.ax.set_xlabel(labels[0], fontsize=lab_size)
        cs.ax.set_ylabel(labels[1], fontsize=lab_size)
        cs.ax.tick_params(axis = 'both', labelsize=tick_size)
        cbar = plt.colorbar(format = FormatStrFormatter('%.2f'))
        cbar.set_label('$'+labels[2]+'$', fontsize = lab_size)
        cbar.ax.tick_params(labelsize=tick_size) 
        # cs.title(fig_title, fontsize=lab_size)

    # plt.show()
    fig.set_figheight(fig_size[0])
    fig.set_figwidth(fig_size[1])        
    fig.savefig(folder+'/'+it+'_'+ise+'_deform_xy_'+type_plot+'.pdf', format='pdf', dpi=600)
    fig.savefig(folder+'/'+it+'_'+ise+'_deform_xy_'+type_plot+'.eps', format='eps', dpi=600)    
    
    

  
def plot_dielectric_functions( dielectric_filepath='', folder=''):
    """
    This function produces 2D plots of real and imaginary parts of the 
    frequency dependent dielectric fuction


    INPUT:
        - dielectric_filepath (str) - path to the OUTCAR of the calculation of the dielectric properties
        - folder (str) - directory where all the results will be built 

    RETURN:
        None
    SOURCE:
        None
    TODO:
        Some improvements
    """
    eV_to_recip_cm = 1.0/(physical_constants['Planck constant in eV s'][0]*speed_of_light*1e2)

 
    # dielectric
    outcar = Outcar(dielectric_filepath)
    outcar.read_freq_dielectric_independent()
    

    freq=outcar.frequencies
    real_dielectric = outcar.real
    imag_dielectric = outcar.imaginary

    # print(dir(outcar))
    # print(outcar.final_energy)
    # Diagonal elements of real dielectric tensor
    real_xx = [real_dielectric[i][0][0] for i in range(len(real_dielectric))]
    real_yy = [real_dielectric[i][1][1] for i in range(len(real_dielectric))]
    real_zz = [real_dielectric[i][2][2] for i in range(len(real_dielectric))]

    # Diagonal elements of imaginary dielectric tensor
    imag_xx = [imag_dielectric[i][0][0] for i in range(len(imag_dielectric))]
    imag_yy = [imag_dielectric[i][1][1] for i in range(len(imag_dielectric))]
    imag_zz = [imag_dielectric[i][2][2] for i in range(len(imag_dielectric))]

    # if ax is None:
    #     fig, ax = plt.subplots(1, 1, figsize=(24.0,12.0))
    # else:
    #     fig = None
    
    # Files
    f = open(folder+'/real_dielectric_'+dielectric_filepath.split('/')[-2]+'.out', 'w')
    for i in range(len(freq)):
        f.write('{0:15.9f} {1:15.9f} {2:15.9f} {3:15.9f}'.format(freq[i],real_xx[i], real_yy[i], real_zz[i])+'\n')
    f.close()

    f = open(folder+'/imaginary_dielectric_'+dielectric_filepath.split('/')[-2]+'.out', 'w')
    for i in range(len(freq)):
        f.write('{0:15.9f} {1:15.9f} {2:15.9f} {3:15.9f}'.format(freq[i],imag_xx[i], imag_yy[i], imag_zz[i])+'\n')
    f.close()

    # Real part of dielectric function
    fig = plt.figure(1)
    plt.plot( freq, real_xx, linewidth=1, linestyle='-', color='red', label='$\\epsilon_1(\\omega)$ xx'  )
    plt.plot( freq, real_yy, linewidth=3, linestyle='--', color='blue', label='$\\epsilon_1(\\omega)$ yy'  )
    plt.plot( freq, real_zz, linewidth=2, linestyle='-', color='green', label='$\\epsilon_1(\\omega)$ zz'  )
    # plt.xlim([0,40])
    # plt.ylim([-20,30])
    plt.xlabel("Energy (eV)", fontsize=18)
    plt.ylabel("$\\epsilon_1$ ($\\omega$)  ", fontsize=25)
    plt.legend(bbox_to_anchor=(0.6, 0.65), borderaxespad=0., labelspacing=0.3, numpoints=1, frameon=True, markerscale=1., handletextpad=0.3, fontsize=17)
    # return fig
    plt.tight_layout()
    plt.savefig(folder+'/'+dielectric_filepath.split('/')[-2]+"_real.png")
    plt.clf()
    plt.cla()


    # Imaginary part of dielectric function
    fig = plt.figure(1)
    plt.plot( freq, imag_xx, linewidth=2, linestyle='-', color='red', label='$\\epsilon_2(\\omega)$ xx'  )
    plt.plot( freq, imag_yy, linewidth=3, linestyle='--', color='blue', label='$\\epsilon_2(\\omega)$ yy'  )
    plt.plot( freq, imag_zz, linewidth=2, linestyle='-', color='green', label='$\\epsilon_2(\\omega)$ zz'  )
    # plt.xlim([0,40])
    # plt.ylim([-20,30])
    plt.xlabel("Energy (eV)", fontsize=18)
    plt.ylabel("$\\epsilon_2$ ($\\omega$)  ", fontsize=25)
    plt.legend(bbox_to_anchor=(0.734, 0.80), borderaxespad=0., labelspacing=0.3, numpoints=1, frameon=True, markerscale=1., handletextpad=0.3, fontsize=17)
    # return fig
    plt.tight_layout()
    plt.savefig(folder+'/'+dielectric_filepath.split('/')[-2]+"_imag.png")
    plt.clf()
    plt.cla()

   


def give_absorption_coeff( dielectric_filepath, folder):
    """
    This function produces 2D plots of frequency dependent optical properties:
    1. Extinction coefficient
    2. Absorption coefficient
    3. Refractive index
    4. Optical reflectivity
    5. Electron energy-loss spectrum


    INPUT:
        - dielectric_filepath (str) - path to the OUTCAR of the calculation of the dielectric properties
        - folder (str) - directory where all the results will be built 

    RETURN:
        None
    SOURCE:
        None
    TODO:
        Some improvements
    """
    eV_to_recip_cm = 1.0/(physical_constants['Planck constant in eV s'][0]*speed_of_light*1e2)
    
    outcar = Outcar(dielectric_filepath)
    outcar.read_freq_dielectric_independent()
   
    
    freq=outcar.frequencies
    real_dielectric = outcar.real
    imag_dielectric = outcar.imaginary
    # Diagonal elements of real dielectric tensor
    real_xx = [real_dielectric[i][0][0] for i in range(len(real_dielectric))]
    real_yy = [real_dielectric[i][1][1] for i in range(len(real_dielectric))]
    real_zz = [real_dielectric[i][2][2] for i in range(len(real_dielectric))]

    # Diagonal elements of imaginary dielectric tensor
    imag_xx = [imag_dielectric[i][0][0] for i in range(len(imag_dielectric))]
    imag_yy = [imag_dielectric[i][1][1] for i in range(len(imag_dielectric))]
    imag_zz = [imag_dielectric[i][2][2] for i in range(len(imag_dielectric))]

    epsilon_1_xx=real_xx
    epsilon_1_yy=real_yy
    epsilon_1_zz=real_zz
    epsilon_2_xx=imag_xx
    epsilon_2_yy=imag_yy
    epsilon_2_zz=imag_zz

    # extinction coefficient k
    k_xx=[]
    k_yy=[]
    k_zz=[]
    for i in range(len(real_dielectric)):
        x=np.sqrt( -epsilon_1_xx[i] + np.sqrt( epsilon_1_xx[i]**2 + epsilon_2_xx[i]**2 ) )*2**(-0.5)
        y=np.sqrt( -epsilon_1_yy[i] + np.sqrt( epsilon_1_yy[i]**2 + epsilon_2_yy[i]**2 ) )*2**(-0.5) 
        z=np.sqrt( -epsilon_1_zz[i] + np.sqrt( epsilon_1_zz[i]**2 + epsilon_2_zz[i]**2 ) )*2**(-0.5) 
        k_xx.append(x)
        k_yy.append(y)
        k_zz.append(z)

    # absorption_coeff in 1e4 (cm-1)
    absor_xx=[]
    absor_yy=[]
    absor_zz=[]
    for i in range(len(real_dielectric)):
        x=2.0 *pi*eV_to_recip_cm*freq[i]*k_xx[i]/1e4
        y=2.0 *pi*eV_to_recip_cm*freq[i]*k_yy[i]/1e4 
        z=2.0 *pi*eV_to_recip_cm*freq[i]*k_zz[i]/1e4 
        absor_xx.append(x)
        absor_yy.append(y)
        absor_zz.append(z)

    # Files
    f = open(folder+'/extinction_coefficient_k_'+dielectric_filepath.split('/')[-2]+'.out', 'w')
    for i in range(len(freq)):
        f.write('{0:15.9f} {1:15.9f} {2:15.9f} {3:15.9f}'.format(freq[i],k_xx[i], k_yy[i], k_zz[i])+'\n')
    f.close()

    f = open(folder+'/absorption_coefficient_'+dielectric_filepath.split('/')[-2]+'.out', 'w')
    for i in range(len(freq)):
        f.write('{0:15.9f} {1:15.9f} {2:15.9f} {3:15.9f}'.format(freq[i],absor_xx[i], absor_yy[i], absor_zz[i])+'\n')
    f.close()

    # Plot absorption_coeff
    fig = plt.figure(1)
    plt.plot( freq, absor_xx, linewidth=2, linestyle='-', color='red', label='$\\alpha(\\omega)$ xx'  )
    plt.plot( freq, absor_yy, linewidth=3, linestyle='--', color='blue', label='$\\alpha(\\omega)$ yy'  )
    plt.plot( freq, absor_zz, linewidth=2, linestyle='-', color='green', label='$\\alpha(\\omega)$ zz'  )
    # plt.xlim([0,40])
    # plt.ylim([-20,30])
    plt.xlabel("Energy (eV)", fontsize=18)
    plt.ylabel("$\\alpha$ ($\\omega$) $10^4cm^{-1}$ ", fontsize=20)
    plt.legend(bbox_to_anchor=(0.55, 0.95), borderaxespad=0., labelspacing=0.3, numpoints=1, frameon=True, markerscale=1., handletextpad=0.3, fontsize=17)
    plt.tight_layout()
    plt.savefig(folder+'/'+dielectric_filepath.split('/')[-2]+"_absorption_coeff.png")
    plt.clf()
    plt.cla()

    # Plot extinction coefficient k
    fig = plt.figure(1)
    plt.plot( freq, k_xx, linewidth=2, linestyle='-', color='red', label='$k(\\omega)$ xx'  )
    plt.plot( freq, k_yy, linewidth=3, linestyle='--', color='blue', label='$k(\\omega)$ yy'  )
    plt.plot( freq, k_zz, linewidth=2, linestyle='-', color='green', label='$k(\\omega)$ zz'  )
    # plt.xlim([0,40])
    # plt.ylim([-20,30])
    plt.xlabel("Energy (eV)", fontsize=18)
    plt.ylabel("$k(\\omega$) ", fontsize=20)
    plt.legend(bbox_to_anchor=(0.55, 0.95), borderaxespad=0., labelspacing=0.3, numpoints=1, frameon=True, markerscale=1., handletextpad=0.3, fontsize=17)
    plt.tight_layout()
    plt.savefig(folder+'/'+dielectric_filepath.split('/')[-2]+"_extinction.png")
    plt.clf()
    plt.cla()

    # Refractive index n
    n_xx=[]
    n_yy=[]
    n_zz=[]
    for i in range(len(real_dielectric)):
        x=np.sqrt( epsilon_1_xx[i] + np.sqrt( epsilon_1_xx[i]**2 + epsilon_2_xx[i]**2 ) )*2**(-0.5)
        y=np.sqrt( epsilon_1_yy[i] + np.sqrt( epsilon_1_yy[i]**2 + epsilon_2_yy[i]**2 ) )*2**(-0.5)
        z=np.sqrt( epsilon_1_zz[i] + np.sqrt( epsilon_1_zz[i]**2 + epsilon_2_zz[i]**2 ) )*2**(-0.5)
        n_xx.append(x)
        n_yy.append(y)
        n_zz.append(z)


    # Files
    f = open(folder+'/refractive_index_n_'+dielectric_filepath.split('/')[-2]+'.out', 'w')
    for i in range(len(freq)):
        f.write('{0:15.9f} {1:15.9f} {2:15.9f} {3:15.9f}'.format(freq[i],n_xx[i], n_yy[i], n_zz[i])+'\n')
    f.close()

    # Plot Refractive index n
    fig = plt.figure(1)
    plt.plot( freq, n_xx, linewidth=2, linestyle='-', color='red', label='$n(\\omega)$ xx'  )
    plt.plot( freq, n_yy, linewidth=3, linestyle='--', color='blue', label='$n(\\omega)$ yy'  )
    plt.plot( freq, n_zz, linewidth=2, linestyle='-', color='green', label='$n(\\omega)$ zz'  )
    # plt.xlim([0,40])
    # plt.ylim([-20,30])
    plt.xlabel("Energy (eV)", fontsize=18)
    plt.ylabel("$n(\\omega$) ", fontsize=20)
    plt.legend(bbox_to_anchor=(0.55, 0.95), borderaxespad=0., labelspacing=0.3, numpoints=1, frameon=True, markerscale=1., handletextpad=0.3, fontsize=17)
    plt.tight_layout()
    plt.savefig(folder+'/'+dielectric_filepath.split('/')[-2]+"_Refractive index.png")
    plt.clf()
    plt.cla()

    # Optical reflectivity R
    R_xx=[]
    R_yy=[]
    R_zz=[]
    for i in range(len(real_dielectric)):
        x=( (n_xx[i]-1)**2 + k_xx[i]**2 ) /( (n_xx[i]+1)**2 + k_xx[i]**2 )
        y=( (n_yy[i]-1)**2 + k_yy[i]**2 ) /( (n_yy[i]+1)**2 + k_yy[i]**2 ) 
        z=( (n_zz[i]-1)**2 + k_zz[i]**2 ) /( (n_zz[i]+1)**2 + k_zz[i]**2 ) 
        R_xx.append(x)
        R_yy.append(y)
        R_zz.append(z)

    # Files
    f = open(folder+'/optical_reflectivity_R_'+dielectric_filepath.split('/')[-2]+'.out', 'w')
    for i in range(len(freq)):
        f.write('{0:15.9f} {1:15.9f} {2:15.9f} {3:15.9f}'.format(freq[i],R_xx[i], R_yy[i], R_zz[i])+'\n')
    f.close()

    # Plot Optical reflectivity R
    fig = plt.figure(1)
    plt.plot( freq, R_xx, linewidth=2, linestyle='-', color='red', label='$R(\\omega)$ xx'  )
    plt.plot( freq, R_yy, linewidth=3, linestyle='--', color='blue', label='$R(\\omega)$ yy'  )
    plt.plot( freq, R_zz, linewidth=2, linestyle='-', color='green', label='$R(\\omega)$ zz'  )
    # plt.xlim([0,40])
    # plt.ylim([-20,30])
    plt.xlabel("Energy (eV)", fontsize=18)
    plt.ylabel("$R(\\omega$) ", fontsize=20)
    plt.legend(bbox_to_anchor=(0.55, 0.95), borderaxespad=0., labelspacing=0.3, numpoints=1, frameon=True, markerscale=1., handletextpad=0.3, fontsize=17)
    plt.tight_layout()
    plt.savefig(folder+'/'+dielectric_filepath.split('/')[-2]+"_optical_reflec.png")
    plt.clf()
    plt.cla()



    # Electron energy-loss spectrum L
    L_xx=[]
    L_yy=[]
    L_zz=[]
    for i in range(len(real_dielectric)):
        x=epsilon_2_xx[i]/( epsilon_1_xx[i]**2 + epsilon_2_xx[i]**2 )
        y=epsilon_2_yy[i]/( epsilon_1_yy[i]**2 + epsilon_2_yy[i]**2 ) 
        z=epsilon_2_zz[i]/( epsilon_1_zz[i]**2 + epsilon_2_zz[i]**2 ) 
        L_xx.append(x)
        L_yy.append(y)
        L_zz.append(z)

    # Files
    f = open(folder+'/eels_L_'+dielectric_filepath.split('/')[-2]+'.out', 'w')
    for i in range(len(freq)):
        f.write('{0:15.9f} {1:15.9f} {2:15.9f} {3:15.9f}'.format(freq[i],L_xx[i], L_yy[i], L_zz[i])+'\n')
    f.close()

    # Plot Electron energy-loss spectrum L
    fig = plt.figure(1)
    plt.plot( freq, L_xx, linewidth=2, linestyle='-', color='red', label='$L(\\omega)$ xx'  )
    plt.plot( freq, L_yy, linewidth=3, linestyle='--', color='blue', label='$L(\\omega)$ yy'  )
    plt.plot( freq, L_zz, linewidth=2, linestyle='-', color='green', label='$L(\\omega)$ zz'  )
    # plt.xlim([0,40])
    # plt.ylim([-20,30])
    plt.xlabel("Energy (eV)", fontsize=18)
    plt.ylabel("$L(\\omega$) ", fontsize=20)
    plt.legend(bbox_to_anchor=(0.55, 0.95), borderaxespad=0., labelspacing=0.3, numpoints=1, frameon=True, markerscale=1., handletextpad=0.3, fontsize=17)
    plt.tight_layout()
    plt.savefig(folder+'/'+dielectric_filepath.split('/')[-2]+"_Elec_energy-loss.png")




    # # print(real_dielectric)
    # print(k_xx)
    # print(absor_xx)
    
   
    

                           
    
    

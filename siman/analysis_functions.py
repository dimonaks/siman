#!/usr/bin/env python3
""" 
Author: Kartamyshev A.I. (Darth Feiwante)
"""

def min_distance(database = None, calculation = ()):
    """
    This function is used to find minimal distance in the lattice presented 
    as a part of the 'Calculation' object

    INPUT:
        - database (.gbdm3) - database dictionary; could be provided; if not then are taken from header
        - calculation (tuple) - tuple describing the Calculation object in form ('structure', 'set', version)
    RETURN:
        None
    SOURCE:
        None
    TODO:
        Some improvements
    """
    c = database[calculation] 
    min_dist = 1000000
    for i in range(c.natom):
        for j in range(c.natom):
            if j == i: continue
            at1 = c.xcart[i]
            at2 = c.xcart[j]
            dist = ((at1[0]-at2[0])**2 + (at1[1]-at2[1])**2 + (at1[2]-at2[2])**2)**0.5
            if min_dist>dist: 
                min_dist = dist
                atom1 = at1
                atom2 = at2
    print('Minimal distance = ', min_dist)
    print('Atom 1 = ', atom1)
    print('Atom 2 = ', atom2)

def formation_energy(database = None, calc_def = (), calc_id = ()):
    """
    This function is used to find minimal distance in the lattice presented 
    as a part of the 'Calculation' object

    INPUT:
        - database (.gbdm3) - database dictionary; could be provided; if not then are taken from header
        - calc_def (tuple) - tuple describing the Calculation object for the 
                             lattice containing a defect in form ('structure', 'set', version)
        - calc_id (tuple) - tuple describing the Calculation object for the 
                             lattice without a defect in form ('structure', 'set', version)
    RETURN:
        None
    SOURCE:
        None
    TODO:
        - Add different type of defects
    """
    defect = database[calc_def]
    ideal = database[calc_id]
    n_at_def = defect.natom
    e_def = defect.energy_free
    e_id_at = ideal.energy_free/ideal.natom
    E_f = e_def - n_at_def*e_id_at
    print('Formation energy for defect '+calc_def[0]+' = '+str(E_f)+' eV')

# Fitting of the E(a,c) dependence for the equilibrium c/a searching

import xalglib

class ALGLIB:

    def build_2d_bicubic_spline(self, x, m, y, n, z, d): self.bicubicv2d = xalglib.spline2dbuildbicubicv(x, m, y, n, z, d)

    def calc(self, x, y, ind): 
        l = xalglib.spline2ddiff(self.bicubicv2d,x,y)
        if ind==0: return l[0]    # z
        elif ind==1: return l[1]  # dz/dx
        elif ind==2: return l[2]  # dz/dy
        elif ind==3: return l[3]  # d2z/dxdy
        else: raise RuntimeError ('Unknown ind = '+str(ind))


class Approximation:

    def aprx_lsq(self, fun_aprx, num_coef, xx, yy):

        from scipy.optimize import leastsq

        coefs = ''
        for i in range(num_coef): coefs += 'a'+str(i)+', '
        coefs+=' = par'
        
        def f_aprx(par, yy1, xx1):
            exec(coefs)
            x1 = []
            for x in xx1: x1.append(eval(fun_aprx))
            return [yy1[i] - x1[i] for i in range(len(x1))]
            
        plsq = leastsq(f_aprx, [1 for i in range(num_coef)], args=(yy, xx))
        self.coefs = list(plsq[0])

def carrier_mobility(calc_init=(), calc_deform_x=(), calc_deform_y=(), 
                     vbm=1, vbm_point=0, cbm=2, cbm_point=0, 
                     deform_x=(), n_points_x=2, deform_y=(), n_points_y=2,
                     temperature_range = (300,900), temperature_ref=300,
                     effective_mass_el={}, effective_mass_xy_el={}, 
                     effective_mass_hole={}, effective_mass_xy_hole={},                      
                     lab_size=15, tick_size=15, leg_size=15, fig_size=(9,17), fig_title='', xlim=(), ylim_elastic=(), ylim_deform=(),
                     expression='Guo2021_JMCC', database=None, folder=''):
    """
    This function is used to calculate carrier mobility for 2D structure

    INPUT:
        - calc_init (tuple) - tuple describing the Calculation object for the 
                              undeformed lattice in form ('structure', 'set', version)
        - calc_deform_x (tuple) - tuple describing the Calculation object for the 
                                  lattice deformed along the 'x' axis in form ('structure', 'set', range(version1, version2))
        - calc_deform_y (tuple) - tuple describing the Calculation object for the 
                                  lattice deformed along the 'y' axis in form ('structure', 'set', range(version1, version2))
        - vbm (int) - number of the last occupied band (valence band) (count starts from '1')
        - vbm_point (int) - number of the k-point in the IBZKPT file, at which the valence band maximum (VBM) is located (count starts from '0')
        - cbm (int) - number of the first unoccupied band (conduction band) (count starts from '1')
        - cbm_point (int) - number of the k-point in the IBZKPT file, at which the conduction band minimum (CBM) is located (count starts from '0')
        - deform_x (tuple of floats) - range of deformations along the 'x' axis in relative units
        - n_points_x (int) - number of different deformations along the 'x' axis
        - deform_y (tuple of floats) - range of deformations along the 'y' axis in relative units
        - n_points_y (int) - number of different deformations along the 'y' axis
        - temperature_range (tuple of floats) - temperature range for the calculations of the temperature-dependent carrier mobility
        - temperature_ref (float) - temperature, for which the detailed information is caculated
        - effective_mass_el (dict) - the dictionary of the effective masses of electrons in different directions of the reciprocal space
                                     has form {'M':0.250,...}, where 'M' - name of the point and 0.250 is the effective mass in 
                                     electron mass units. The point determines the direction 'CBM' -> 'M'.
        - effective_mass_xy_el (dict) - the same as 'effective_mass_el' but for directions 'CBM' -> 'X' and 'CBM' -> 'Y' only                                  
        - effective_mass_hole (dict) - the same as 'effective_mass_el' but for holes
        - effective_mass_xy_el (dict) - the same as 'effective_mass_xy_el' but for holes
        - lab_size (int) - font size of labels
        - tick_size (int) - font size of ticks
        - leg_size (int) - font size of legend
        - fig_size (tuple of floats) - has form (height, width) set size of the figure in units of cm
        - fig_title (str) - optional, the title placed above the plot
        - xlim (tuple of floats) - has form (min_x, max_x) and represents limits of deformation for
                                   both plots of the elastic moduli and deformation potential
        - ylim_elastic (tuple of floats) - has form (min_y, max_y) and represents limits for
                                           the elastic moduli plot, units are N/m
        - ylim_deform (tuple of floats) - has form (min_y, max_y) and represents limits for
                                           the deformation plot, units are eV 
        - expression (str) - key representing the formula to calculate carrier mobility
                             Possible options:
                             'Guo2021_JMCC' - Eq. (4) from J. Mater. Chem. C. 9 (2021) 2464–2473. doi:10.1039/D0TC05649A
        - database (.gbdm3) - database dictionary; could be provided; if not then are taken from header
        - folder (str) - directory where all the results will be built 
    RETURN:
        None
    SOURCE:
        Guo2021_JMCC equation - S.-D. Guo, W.-Q. Mu, Y.-T. Zhu, R.-Y. Han, W.-C. Ren, J. Mater. Chem. C. 9 (2021) 2464–2473. doi:10.1039/D0TC05649A. (Eq. 4)
    TODO:
        - Add calculations for 3D structures
    """
    from math import pi, sqrt, acos, sin
    from pymatgen.io.vasp import Vasprun, Outcar
    from pymatgen.electronic_structure.core import Spin, OrbitalType
    from analysis_functions import Approximation as A
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gspec

    # Transformation coefficients
    ang     = 1e-10          # Angstrom, m
    ev_in_j = 1.60217662e-19 # Electron-volt, J

    # Functions for calculation of the carrier mobility
    def carrier_Guo2021_JMCC(t, c2d, el, m, md):
        # Constants
        e   = 1.60217662e-19 # Electron charge, Coulomb
        m_e = 9.10938356e-31 # Electron mass, kg
        k_B = 1.38064852e-23 # Boltzmann constant, J/K
        h   = 6.62607004e-34 # Plank constant, J*s
        mu = (e * (h/(2*pi))**3 * c2d * 10000)/(k_B*t*m*md*m_e**2*el**2)
        return mu

    # Lists of deformation along x and y axes
    step_x = (deform_x[1] - deform_x[0])/(n_points_x - 1)
    eps_x = [deform_x[0] + step_x*i for i in range(n_points_x)]

    step_y = (deform_y[1] - deform_y[0])/(n_points_y - 1)
    eps_y = [deform_y[0] + step_y*i for i in range(n_points_y)]

    # Reading energies
    energy_x = []
    energy_vbm_x = []
    energy_cbm_x = []
    for i in calc_deform_x[2]:
        calc_cur = database[(calc_deform_x[0], calc_deform_x[1], i)]
        energy_x.append(calc_cur.energy_sigma0)

        path_vasprun = calc_cur.path['output'].replace('OUTCAR','vasprun.xml')
        path_kpoints = calc_cur.path['output'].replace(str(i)+'.OUTCAR','IBZKPT')
        run = Vasprun(path_vasprun, parse_projected_eigen=True)
        dosrun = Vasprun(path_vasprun)
        try:
            bands = run.get_band_structure(path_kpoints,
                                   line_mode=True,
                                   efermi=dosrun.efermi)
        except Exception:
            path_kpoints = calc_cur.path['output'].replace(str(i)+'.OUTCAR','KPOINTS')
            bands = run.get_band_structure(path_kpoints,
                                   line_mode=True,
                                   efermi=dosrun.efermi)

        energy_vbm_x.append(bands.bands[Spin.up][vbm-1][vbm_point])
        energy_cbm_x.append(bands.bands[Spin.up][cbm-1][cbm_point])

    energy_y = []
    energy_vbm_y = []
    energy_cbm_y = []   
    for i in calc_deform_y[2]:  
        calc_cur = database[(calc_deform_y[0], calc_deform_y[1], i)]
        energy_y.append(calc_cur.energy_sigma0) 

        path_vasprun = calc_cur.path['output'].replace('OUTCAR','vasprun.xml')
        path_kpoints = calc_cur.path['output'].replace(str(i)+'.OUTCAR','IBZKPT')
        run = Vasprun(path_vasprun, parse_projected_eigen=True)
        dosrun = Vasprun(path_vasprun)
        try:
            bands = run.get_band_structure(path_kpoints,
                                   line_mode=True,
                                   efermi=dosrun.efermi)
        except Exception:
            path_kpoints = calc_cur.path['output'].replace(str(i)+'.OUTCAR','KPOINTS')
            bands = run.get_band_structure(path_kpoints,
                                   line_mode=True,
                                   efermi=dosrun.efermi)
        
        energy_vbm_y.append(bands.bands[Spin.up][vbm-1][vbm_point])
        energy_cbm_y.append(bands.bands[Spin.up][cbm-1][cbm_point])

    # Calculation of the initial surface area s0
    a = list(database[calc_init].end.rprimd[0])
    b = list(database[calc_init].end.rprimd[1])
    mod_a = sqrt(a[0]**2 + a[1]**2 + a[2]**2)
    mod_b = sqrt(b[0]**2 + b[1]**2 + b[2]**2)
    sc_ab = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    angle = acos(sc_ab/(mod_a * mod_b))
    s0 = mod_a * mod_b * sin(angle)

    # ********************************************* 2D elastic moduli ********************************************
    # 2D Elastic modulus
    # Deformation along x direction
    ek1 = A()
    ek1.aprx_lsq('a0*x**2+a1*x+a2', 3, eps_x, energy_x)
    C_2D_x = ek1.coefs[0]*ev_in_j/(s0*ang**2)
    a0_x = ek1.coefs[0]
    a1_x = ek1.coefs[1]
    a2_x = ek1.coefs[2]

    eps_x_fit = [deform_x[0] + (deform_x[1] - deform_x[0])/999*i for i in range(1000)]
    energy_x_fit = [a0_x*i**2 + a1_x*i + a2_x for i in eps_x_fit]

    # Deformation along y direction    
    ek2 = A()
    ek2.aprx_lsq('a0*x**2+a1*x+a2', 3, eps_y, energy_y)
    C_2D_y = ek2.coefs[0]*ev_in_j/(s0*ang**2)
    a0_y = ek2.coefs[0]
    a1_y = ek2.coefs[1]
    a2_y = ek2.coefs[2]

    eps_y_fit = [deform_y[0] + (deform_y[1] - deform_y[0])/999*i for i in range(1000)]
    energy_y_fit = [a0_y*i**2 + a1_y*i + a2_y for i in eps_y_fit]
    # ************************************************************************************************************

    # ***************************************** Deformation potential ********************************************
    # Deformation potential
    # X direction
    ek3 = A()
    ek3.aprx_lsq('a0*x+a1', 2, eps_x, energy_vbm_x)
    E_l_hole_x = ek3.coefs[0]*ev_in_j
    vbm_x_0 = ek3.coefs[0]
    vbm_x_1 = ek3.coefs[1]

    ek4 = A()
    ek4.aprx_lsq('a0*x+a1', 2, eps_x, energy_cbm_x)
    E_l_el_x = ek4.coefs[0]*ev_in_j
    cbm_x_0 = ek4.coefs[0]
    cbm_x_1 = ek4.coefs[1]

    energy_vbm_x_fit = [vbm_x_0*i + vbm_x_1 for i in eps_x_fit]
    energy_cbm_x_fit = [cbm_x_0*i + cbm_x_1 for i in eps_x_fit]    

    # Y direction
    ek5 = A()
    ek5.aprx_lsq('a0*x+a1', 2, eps_y, energy_vbm_y)
    E_l_hole_y = ek5.coefs[0]*ev_in_j
    vbm_y_0 = ek5.coefs[0]
    vbm_y_1 = ek5.coefs[1]

    ek6 = A()
    ek6.aprx_lsq('a0*x+a1', 2, eps_y, energy_cbm_y)
    E_l_el_y = ek6.coefs[0]*ev_in_j
    cbm_y_0 = ek6.coefs[0]
    cbm_y_1 = ek6.coefs[1]

    energy_vbm_y_fit = [vbm_y_0*i + vbm_y_1 for i in eps_y_fit]
    energy_cbm_y_fit = [cbm_y_0*i + cbm_y_1 for i in eps_y_fit]  

    # Writing files with the energy dependencies on deformation value
    f = open(folder+'/elastic_x.out', 'w')
    for i in range(len(eps_x)):
        f.write('{0:15.5f} {1:15.10f}'.format(eps_x[i], energy_x[i])+'\n')
    f.close()

    f = open(folder+'/elastic_y.out', 'w')
    for i in range(len(eps_y)):
        f.write('{0:15.5f} {1:15.10f}'.format(eps_y[i], energy_y[i])+'\n')
    f.close()

    f = open(folder+'/elastic_x_fit.out', 'w')
    for i in range(len(eps_x_fit)):
        f.write('{0:15.5f} {1:15.10f}'.format(eps_x_fit[i], energy_x_fit[i])+'\n')
    f.close()

    f = open(folder+'/elastic_y_fit.out', 'w')
    for i in range(len(eps_y_fit)):
        f.write('{0:15.5f} {1:15.10f}'.format(eps_y_fit[i], energy_y_fit[i])+'\n')
    f.close()

    f = open(folder+'/deformation_potential_x_hole.out', 'w')
    for i in range(len(eps_x)):
        f.write('{0:15.5f} {1:15.10f}'.format(eps_x[i], energy_vbm_x[i])+'\n')
    f.close()    

    f = open(folder+'/deformation_potential_x_electron.out', 'w')
    for i in range(len(eps_x)):
        f.write('{0:15.5f} {1:15.10f}'.format(eps_x[i], energy_cbm_x[i])+'\n')
    f.close()   

    f = open(folder+'/deformation_potential_x_hole_fit.out', 'w')
    for i in range(len(eps_x_fit)):
        f.write('{0:15.5f} {1:15.10f}'.format(eps_x_fit[i], energy_vbm_x_fit[i])+'\n')
    f.close()    

    f = open(folder+'/deformation_potential_x_electron_fit.out', 'w')
    for i in range(len(eps_x_fit)):
        f.write('{0:15.5f} {1:15.10f}'.format(eps_x_fit[i], energy_cbm_x_fit[i])+'\n')
    f.close()     

    f = open(folder+'/deformation_potential_y_hole.out', 'w')
    for i in range(len(eps_y)):
        f.write('{0:15.5f} {1:15.10f}'.format(eps_y[i], energy_vbm_y[i])+'\n')
    f.close()    

    f = open(folder+'/deformation_potential_y_electron.out', 'w')
    for i in range(len(eps_y)):
        f.write('{0:15.5f} {1:15.10f}'.format(eps_y[i], energy_cbm_y[i])+'\n')
    f.close()   

    f = open(folder+'/deformation_potential_y_hole_fit.out', 'w')
    for i in range(len(eps_y_fit)):
        f.write('{0:15.5f} {1:15.10f}'.format(eps_y_fit[i], energy_vbm_y_fit[i])+'\n')
    f.close()    

    f = open(folder+'/deformation_potential_y_electron_fit.out', 'w')
    for i in range(len(eps_y_fit)):
        f.write('{0:15.5f} {1:15.10f}'.format(eps_y_fit[i], energy_cbm_y_fit[i])+'\n')
    f.close()    
    # ************************************************************************************************************

    # Figure with Elastic moduli and Deformation potential
    interval = 0.05
    npoint = 50
    precision = 0.000001
    shift_left = 100000.15
    shift_right = 100000.20
    # Use GridSpec to set subplots
    gs = gspec.GridSpec(1, 2, width_ratios=[1.0, 1.0],  height_ratios = [1.0]) #  height_ratios = [5, 5, 5]
    gs.update(bottom=0.07, hspace=0.17, wspace=0.25, top=0.93, right=0.97, left=0.1)
    # ========================== Vanadium =================================
    plt1 = plt.subplot(gs[0,0])
    plt1.plot(eps_x,     energy_x,     linewidth=2, linestyle='' , marker='o', markersize = 7, color='red' )
    plt1.plot(eps_x_fit, energy_x_fit, linewidth=2, linestyle='-', color='red' , label='x: $C_{2D} = $'+'{0:7.2f}'.format(C_2D_x)+' N/m')
    plt1.plot(eps_y,     energy_y,     linewidth=2, linestyle='' , marker='o', markersize = 7,  color='blue')
    plt1.plot(eps_y_fit, energy_y_fit, linewidth=2, linestyle='-', color='blue', label='y: $C_{2D} = $'+'{0:7.2f}'.format(C_2D_y)+' N/m')
    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.xlabel('$\epsilon$', fontsize = lab_size)
    plt.ylabel('Energy, eV', fontsize = lab_size)
    if xlim: plt1.set_xlim(xlim[0], xlim[1])
    if ylim_elastic: plt1.set_ylim(ylim_elastic[0], ylim_elastic[1])
    plt1.legend(bbox_to_anchor=(0.5, 0.5), borderaxespad=0., labelspacing=0.3, numpoints=1, frameon=True, fancybox=True, markerscale=1., handletextpad=0.3, fontsize=leg_size)
   
    # ========================== Chromium =================================
    plt2 = plt.subplot(gs[0,1])
    plt2.plot(eps_x,     energy_vbm_x,     linewidth=2, linestyle='' , marker='o', markersize = 7, color='red' )
    plt2.plot(eps_x_fit, energy_vbm_x_fit, linewidth=2, linestyle='-', color='red', label='x: VBM $E_{l} = $'+'{0:7.2f}'.format(E_l_hole_x/ev_in_j)+' eV')
    plt2.plot(eps_x,     energy_cbm_x,     linewidth=2, linestyle='' , marker='o', markersize = 7,  color='blue')
    plt2.plot(eps_x_fit, energy_cbm_x_fit, linewidth=2, linestyle='-', color='blue', label='x: CBM $E_{l} = $'+'{0:7.2f}'.format(E_l_el_x/ev_in_j)+' eV')
    plt2.plot(eps_y,     energy_vbm_y,     linewidth=2, linestyle='' , marker='s', markersize = 7, color='orange' )
    plt2.plot(eps_y_fit, energy_vbm_y_fit, linewidth=2, linestyle='-', color='orange', label='y: VBM $E_{l} = $'+'{0:7.2f}'.format(E_l_hole_y/ev_in_j)+' eV')
    plt2.plot(eps_y,     energy_cbm_y,     linewidth=2, linestyle='' , marker='s', markersize = 7,  color='green')
    plt2.plot(eps_y_fit, energy_cbm_y_fit, linewidth=2, linestyle='-', color='green', label='y: CBM $E_{l} = $'+'{0:7.2f}'.format(E_l_el_y/ev_in_j)+' eV')    
    plt.xticks(fontsize=tick_size)
    plt.yticks(fontsize=tick_size)
    plt.xlabel('$\epsilon$', fontsize = lab_size)
    plt.ylabel('Energy, eV', fontsize = lab_size)
    if xlim: plt2.set_xlim(xlim[0], xlim[1])
    if ylim_deform: plt2.set_ylim(ylim_deform[0], ylim_deform[1])
    plt2.legend(bbox_to_anchor=(0.5, 0.5), borderaxespad=0., labelspacing=0.3, numpoints=1, frameon=True, fancybox=True, markerscale=1., handletextpad=0.3, fontsize=leg_size)
    
    # Make figure with the chosen size
    fig = plt.figure(1, dpi=600)
    fig.set_figheight(fig_size[0])
    fig.set_figwidth(fig_size[1])
    fig.savefig(folder+'/Elastic_moduli_Deformation_potential.pdf', format="pdf", dpi=600)
    plt.clf()
    plt.cla()
    # Calculation of the carrier mobility
    # Check if the temperature range option is used
    if temperature_range:
        step_t = (temperature_range[1] - temperature_range[0])/(1000 - 1)
        t_list = [temperature_range[0] + i*step_t for i in range(1000)]

    # Make calculations
    f = open(folder+'/carrier_mobility_t'+str(temperature_ref)+'.out', 'w')
    if expression == 'Guo2021_JMCC':
        cite = 'The expression for the carrier mobility was taken from \n S.-D. Guo, W.-Q. Mu, Y.-T. Zhu, R.-Y. Han, W.-C. Ren, J. Mater. Chem. C. 9 (2021) 2464–2473. doi:10.1039/D0TC05649A. (Eq. 4)\n'
    f.write(cite)
    f.write('\n')
    f.write('Temperature '+str(temperature_ref)+' K\n')
    f.write('{0:^15s}{1:^10s}{2:^15s}{3:^15s}{4:^15s}{5:^17s}'.format('Carrier type', 'Direction', 'C_2D, N/m', 'm*, m_e', 'E_l, eV', 'mu_2D, cm^2V^(-1)s^(-1)')+'\n')
    
    if expression == 'Guo2021_JMCC':
        function_cm = carrier_Guo2021_JMCC
    else:
        raise RuntimeError('Choose the appropriate expression for carrier mobility!!!')


    # Electrons
    effective_mass_xy_el_md = (effective_mass_xy_el['X']*effective_mass_xy_el['Y'])**0.5

    for j in effective_mass_el.keys():
        f.write('CBM -> '+j+'\n')
        if expression == 'Guo2021_JMCC':
            mu_el_x = function_cm(t=temperature_ref, c2d=C_2D_x, el=E_l_el_x, m=effective_mass_el[j], md=effective_mass_xy_el_md)
            mu_el_y = function_cm(t=temperature_ref, c2d=C_2D_y, el=E_l_el_y, m=effective_mass_el[j], md=effective_mass_xy_el_md)

        f.write('{0:^15s}{1:^10s}{2:^15.2f}{3:^15.2f}{4:^15.2f}{5:^17.2f}'.format('Electron', 'x', C_2D_x, effective_mass_el[j], E_l_el_x/ev_in_j, mu_el_x)+'\n')
        f.write('{0:^15s}{1:^10s}{2:^15.2f}{3:^15.2f}{4:^15.2f}{5:^17.2f}'.format('Electron', 'y', C_2D_y, effective_mass_el[j], E_l_el_y/ev_in_j, mu_el_y)+'\n')
        
        # If the temperature range is set
        if temperature_range:
            mu_el_x_list = [function_cm(t=i, c2d=C_2D_x, el=E_l_el_x, m=effective_mass_el[j], md=effective_mass_xy_el_md) for i in t_list]
            mu_el_y_list = [function_cm(t=i, c2d=C_2D_y, el=E_l_el_y, m=effective_mass_el[j], md=effective_mass_xy_el_md) for i in t_list]
            
            f1 = open(folder+'/carrier_mobility_t'+str(temperature_range[0])+'_'+str(temperature_range[1])+'_el_CBM_'+j+'_X.out', 'w')
            for i in range(len(t_list)):
                f1.write('{0:^15.2f} {1:^15.2f}'.format(t_list[i], mu_el_x_list[i])+'\n')
            f1.close()

            f1 = open(folder+'/carrier_mobility_t'+str(temperature_range[0])+'_'+str(temperature_range[1])+'_el_CBM_'+j+'_Y.out', 'w')
            for i in range(len(t_list)):
                f1.write('{0:^15.2f} {1:^15.2f}'.format(t_list[i], mu_el_y_list[i])+'\n')
            f1.close()

            # Make figure for the temperature dependence of the carrier mobility of electrons
            plt.plot(t_list, mu_el_x_list, linewidth=2, linestyle='-', color='red', label = 'X')
            plt.plot(t_list, mu_el_y_list, linewidth=2, linestyle='-', color='blue', label = 'Y')
            plt.xlabel('Temperature, K', fontsize = lab_size)
            plt.ylabel('$\mu_{2D}, cm^2V^{-1}s^{-1}$', fontsize = lab_size)
            plt.legend(bbox_to_anchor=(0.5, 0.5), borderaxespad=0., labelspacing=0.3, numpoints=1, frameon=True, fancybox=True, markerscale=1., handletextpad=0.3, fontsize=leg_size)
            plt.title('Electron CBM ->'+j)
            fig = plt.figure(1)
            fig.set_figheight(9)
            fig.set_figwidth(9)
            fig.savefig(folder+'/carrier_mobility_t'+str(temperature_range[0])+'_'+str(temperature_range[1])+'_el_CBM_'+j+'_XY.pdf', format="pdf", dpi=600)
            plt.clf()
            plt.cla()

    # Holes
    effective_mass_xy_hole_md = (effective_mass_xy_hole['X']*effective_mass_xy_hole['Y'])**0.5

    for j in effective_mass_hole.keys():
        f.write('VBM -> '+j+'\n')
        if expression == 'Guo2021_JMCC':
            mu_hole_x = function_cm(t=temperature_ref, c2d=C_2D_x, el=E_l_hole_x, m=effective_mass_hole[j], md=effective_mass_xy_hole_md)
            mu_hole_y = function_cm(t=temperature_ref, c2d=C_2D_y, el=E_l_hole_y, m=effective_mass_hole[j], md=effective_mass_xy_hole_md)

        f.write('{0:^15s}{1:^10s}{2:^15.2f}{3:^15.2f}{4:^15.2f}{5:^17.2f}'.format('Hole', 'x', C_2D_x, effective_mass_hole[j], E_l_hole_x/ev_in_j, mu_hole_x)+'\n')
        f.write('{0:^15s}{1:^10s}{2:^15.2f}{3:^15.2f}{4:^15.2f}{5:^17.2f}'.format('Hole', 'y', C_2D_y, effective_mass_hole[j], E_l_hole_y/ev_in_j, mu_hole_y)+'\n')

        # If the temperature range is set
        if temperature_range:
            mu_hole_x_list = [function_cm(t=i, c2d=C_2D_x, el=E_l_hole_x, m=effective_mass_hole[j], md=effective_mass_xy_hole_md) for i in t_list]
            mu_hole_y_list = [function_cm(t=i, c2d=C_2D_y, el=E_l_hole_y, m=effective_mass_hole[j], md=effective_mass_xy_hole_md) for i in t_list]

            f1 = open(folder+'/carrier_mobility_t'+str(temperature_range[0])+'_'+str(temperature_range[1])+'_hole_VBM_'+j+'_X.out', 'w')
            for i in range(len(t_list)):
                f1.write('{0:^15.2f} {1:^15.2f}'.format(t_list[i], mu_hole_x_list[i])+'\n')
            f1.close()

            f1 = open(folder+'/carrier_mobility_t'+str(temperature_range[0])+'_'+str(temperature_range[1])+'_hole_VBM_'+j+'_Y.out', 'w')
            for i in range(len(t_list)):
                f1.write('{0:^15.2f} {1:^15.2f}'.format(t_list[i], mu_hole_y_list[i])+'\n')
            f1.close()
        
            # Make figure for the temperature dependence of the carrier mobility of holes
            plt.plot(t_list, mu_hole_x_list, linewidth=2, linestyle='-', color='red', label = 'X')
            plt.plot(t_list, mu_hole_y_list, linewidth=2, linestyle='-', color='blue', label = 'Y')
            plt.xlabel('Temperature, K', fontsize = lab_size)
            plt.ylabel('$\mu_{2D}, cm^2V^{-1}s^{-1}$', fontsize = lab_size)
            plt.legend(bbox_to_anchor=(0.5, 0.5), borderaxespad=0., labelspacing=0.3, numpoints=1, frameon=True, fancybox=True, markerscale=1., handletextpad=0.3, fontsize=leg_size)
            plt.title('Hole CBM ->'+j)
            fig = plt.figure(1)
            fig.set_figheight(9)
            fig.set_figwidth(9)
            fig.savefig(folder+'/carrier_mobility_t'+str(temperature_range[0])+'_'+str(temperature_range[1])+'_hole_VBM_'+j+'_XY.pdf', format="pdf", dpi=600)
            plt.clf()
            plt.cla()

    f.close()
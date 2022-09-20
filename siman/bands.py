#!/usr/bin/env python
# -*- coding=utf-8 -*-

import sys
import numpy as np


import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec

# import pymatgen.io
# print(dir(pymatgen.io))
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType


from siman.small_functions import makedir

def rgbline(ax, k, e, red, green, blue, alpha=0.2):
    """
    This function is used to provide colour for the line of the band structure plot
    depending on the contribution of different orbitals of the specific element.
    It is used in the 'plot_bands' function. 

    INPUT:
        - ax (matplotlib plot object) - object of the band structure plot
        - k (list of ints) - list of numbers of k-points in the reciprocal space
        - e (list of floats) - list of energies corresponding to the k-points from the parameter 'k'
        - red (<class 'numpy.ndarray'>) - contribution from s orbitals
        - green (<class 'numpy.ndarray'>) - contribution from p orbitals
        - blue (<class 'numpy.ndarray'>) - contribution from d orbitals
        - alpha (float) - 

    RETURN:
        None
    SOURCE:
        http://nbviewer.ipython.org/urls/raw.github.com/dpsanders/matplotlib-examples/master/colorline.ipynb
    TODO:
        Some improvements
    """

    pts = np.array([k, e]).T.reshape(-1, 1, 2)
    seg = np.concatenate([pts[:-1], pts[1:]], axis=1)

    nseg = len(k) - 1
    r = [0.5 * (red[i] + red[i + 1]) for i in range(nseg)]
    g = [0.5 * (green[i] + green[i + 1]) for i in range(nseg)]
    b = [0.5 * (blue[i] + blue[i + 1]) for i in range(nseg)]
    a = np.ones(nseg, np.float) * alpha
    lc = LineCollection(seg, colors=list(zip(r, g, b, a)), linewidth=2)
    ax.add_collection(lc)


def read_kpoint_labels(filename):
    """
    Read commented kpoint labels from VASP KPOINTS file 

    INPUT:
        - filename (str) - path to the KPOINTS file
    RETURN:
        None
    SOURCE:
        None
    TODO:
        Some improvements
    """
    labels = []
    with open(filename, 'r') as f:
        [f.readline() for i in range(4)]
        # next(f)
        lab = ''
        for line in f:
            # print (line)
            if '!' in line:
                lab_next = line.split('!')[1].strip()
                # print (lab_next, 'q')
                if lab_next and lab_next != lab:
                    # print (lab_next)
                    labels.append(lab_next)
                lab = lab_next
    return labels


def plot_bands_old(vasprun_dos, vasprun_bands, kpoints, element, ylim = (None, None)):
    """
    This function is used to build plot of the electronic band structure along with
    the density of states (DOS) plot. It has feature to provide contributions from different
    elements to both band structure and DOS

    INPUT:
        - vasprun_dos (str) - path to the vasprun file of the DOS calculation
        - vasprun_band (str) - path to the vasprun file of the band structure calculation
        - kpoints (str) - path to the KPOINTS file of the band structure calculation
        - element (str) - label of the chemical element, for which the contribution
                          to the band structure and DOS
        - ylim (tuple of floats) - energy range of the band structure and DOS plots, units are eV
    RETURN:
        None
    SOURCE:
        Credit https://github.com/gVallverdu/bandstructureplots
    TODO:
        Some improvements
    """



    labels = read_kpoint_labels(kpoints)


    # density of states
    # dosrun = Vasprun(vasprun_dos)
    dosrun = Vasprun(vasprun_bands)
    spd_dos = dosrun.complete_dos.get_spd_dos()

    # bands
    run = Vasprun(vasprun_bands, parse_projected_eigen=True)
    bands = run.get_band_structure(kpoints,
                                   line_mode=True,
                                   efermi=dosrun.efermi)

    # set up matplotlib plot
    # ----------------------

    # general options for plot
    font = {'family': 'serif', 'size': 24}
    plt.rc('font', **font)

    # set up 2 graph with aspec ration 2/1
    # plot 1: bands diagram
    # plot 2: Density of States
    gs = GridSpec(1, 2, width_ratios=[2, 1])
    fig = plt.figure(figsize=(11.69, 8.27))
    # fig.suptitle("Bands diagram of copper")
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])  # , sharey=ax1)

    # set ylim for the plot
    # ---------------------
    if ylim[0]:
        emin = ylim[0]
    else: 
        emin = -10.
    if ylim[1]:
        emax = ylim[1]
    else: 
        emax =  10.

    ax1.set_ylim(emin, emax)
    ax2.set_ylim(emin, emax)

    # Band Diagram
    # ------------
    name = element 
    pbands = bands.get_projections_on_elements_and_orbitals({name: ["s", "p", "d"]})
    # print(bands)

    # compute s, p, d normalized contributions
    contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    # print(pbands)
    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):
            # print(Spin.up)
            sc = pbands[Spin.up][b][k][name]["s"]**2
            pc = pbands[Spin.up][b][k][name]["p"]**2
            dc = pbands[Spin.up][b][k][name]["d"]**2
            tot = sc + pc + dc
            if tot != 0.0:
                contrib[b, k, 0] = sc / tot
                contrib[b, k, 1] = pc / tot
                contrib[b, k, 2] = dc / tot

    # plot bands using rgb mapping
    for b in range(bands.nb_bands):
        rgbline(ax1,
                range(len(bands.kpoints)),
                [e - bands.efermi for e in bands.bands[Spin.up][b]],
                contrib[b, :, 0],
                contrib[b, :, 1],
                contrib[b, :, 2])

    # style
    ax1.set_xlabel("k-points")
    ax1.set_ylabel(r"$E - E_f$   /   eV")
    ax1.grid()

    # fermi level at 0
    ax1.hlines(y=0, xmin=0, xmax=len(bands.kpoints), color="k", linestyle = '--', lw=1)

    # labels
    nlabs = len(labels)
    step = len(bands.kpoints) / (nlabs - 1)
    for i, lab in enumerate(labels):
        ax1.vlines(i * step, emin, emax, "k")
    ax1.set_xticks([i * step for i in range(nlabs)])
    ax1.set_xticklabels(labels)

    ax1.set_xlim(0, len(bands.kpoints))

    # Density of states
    # ----------------

    ax2.set_yticklabels([])
    ax2.grid()
    ax2.set_xlim(1e-4, 5)
    ax2.set_xticklabels([])
    ax2.hlines(y=0, xmin=0, xmax=5, color="k", lw=2)
    ax2.set_xlabel("Density of States", labelpad=28)

    # spd contribution
    ax2.plot(spd_dos[OrbitalType.s].densities[Spin.up],
             dosrun.tdos.energies - dosrun.efermi,
             "r-", label="3s", lw=2)
    ax2.plot(spd_dos[OrbitalType.p].densities[Spin.up],
             dosrun.tdos.energies - dosrun.efermi,
             "g-", label="3p", lw=2)
    ax2.plot(spd_dos[OrbitalType.d].densities[Spin.up],
             dosrun.tdos.energies - dosrun.efermi,
             "b-", label="3d", lw=2)

    # total dos
    ax2.fill_between(dosrun.tdos.densities[Spin.up],
                     0,
                     dosrun.tdos.energies - dosrun.efermi,
                     color=(0.7, 0.7, 0.7),
                     facecolor=(0.7, 0.7, 0.7))

    ax2.plot(dosrun.tdos.densities[Spin.up],
             dosrun.tdos.energies - dosrun.efermi,
             color=(0.6, 0.6, 0.6),
             label="total DOS")

    # plot format style
    # -----------------
    ax2.legend(fancybox=True, shadow=True, prop={'size': 18})
    plt.subplots_adjust(wspace=0)

    # plt.show()
    makedir("figs/bands.png")
    plt.savefig("figs/bands.png")

def plot_bands(vasprun_dos, vasprun_bands, kpoints, element, ylim = (None, None), folder = '', renew_folder=True, vb_top=0, cb_bottom=1, vbm_pos=0, cbm_pos=0):
    """
    This function is used to build plot of the electronic band structure along with
    the density of states (DOS) plot. It has feature to provide contributions from different
    elements to both band structure and DOS. In addition, the band gap (in eV) is automatically
    calculated.

    INPUT:
        - vasprun_dos (str) - path to the vasprun file of the DOS calculation
        - vasprun_band (str) - path to the vasprun file of the band structure calculation
        - kpoints (str) - path to the KPOINTS file of the band structure calculation
        - element (str) - label of the chemical element, for which the contribution
                          to the band structure and DOS
        - ylim (tuple of floats) - energy range of the band structure and DOS plots, units are eV
        - folder (str) - directory where all the results will be built
        - renew_folder (bool) - if True then the folder will be renewed with removing old one
        - vb_top (int) - number of the last occupied band (valence band) (count starts from '1')
        - vbm_pos (int) - supposed number of the k-point in the IBZKPT file, at which the valence band maximum (VBM) is located (count starts from '0')
        - cb_bottom (int) - number of the first unoccupied band (conduction band) (count starts from '1')
        - cbm_pos (int) - supposed number of the k-point in the IBZKPT file, at which the conduction band minimum (CBM) is located (count starts from '0')
 
    RETURN:
        None
    SOURCE:
        Credit https://github.com/gVallverdu/bandstructureplots
    TODO:
        Some improvements
    """    

    from siman.small_functions import makedir
    # The folder to collect data necessary for graph building
    full_name_folder_data = folder+vasprun_bands.split('/')[-2]
    print('bands.py, string 274, full_name_folder_data ', full_name_folder_data)
    makedir(full_name_folder_data+'/', renew_folder = renew_folder)




    labels = read_kpoint_labels(kpoints)

    # density of states
    # dosrun = Vasprun(vasprun_dos)
    dosrun = Vasprun(vasprun_bands)
    spd_dos = dosrun.complete_dos.get_spd_dos()

    # bands
    run = Vasprun(vasprun_bands, parse_projected_eigen=True)
    bands = run.get_band_structure(kpoints,
                                   line_mode=True,
                                   efermi=dosrun.efermi)

    # set up matplotlib plot
    # ----------------------

    # general options for plot
    font = {'family': 'serif', 'size': 24}
    plt.rc('font', **font)

    # set up 2 graph with aspec ration 2/1
    # plot 1: bands diagram
    # plot 2: Density of States
    gs = GridSpec(1, 2, width_ratios=[2, 1])
    fig = plt.figure(figsize=(11.69, 8.27))
    # fig.suptitle("Bands diagram of copper")
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])  # , sharey=ax1)

    # set ylim for the plot
    # ---------------------
    if ylim[0]:
        emin = ylim[0]
    else: 
        emin = -10.
    if ylim[1]:
        emax = ylim[1]
    else: 
        emax =  10.

    ax1.set_ylim(emin, emax)
    ax2.set_ylim(emin, emax)

    # Band Diagram
    # ------------
    name = element 
    pbands = bands.get_projections_on_elements_and_orbitals({name: ["s", "p", "d"]})
    # print(pbands)

    # compute s, p, d normalized contributions
    contrib = np.zeros((bands.nb_bands, len(bands.kpoints), 3))
    for b in range(bands.nb_bands):
        for k in range(len(bands.kpoints)):
            sc = pbands[Spin.up][b][k][name]["s"]**2
            pc = pbands[Spin.up][b][k][name]["p"]**2
            dc = pbands[Spin.up][b][k][name]["d"]**2
            tot = sc + pc + dc
            if tot != 0.0:
                contrib[b, k, 0] = sc / tot
                contrib[b, k, 1] = pc / tot
                contrib[b, k, 2] = dc / tot

    # plot bands using rgb mapping
    for b in range(bands.nb_bands):
        edif = [e - bands.efermi for e in bands.bands[Spin.up][b]]
        rgbline(ax1,
                range(len(bands.kpoints)),
                [e - bands.efermi for e in bands.bands[Spin.up][b]],
                contrib[b, :, 0],
                contrib[b, :, 1],
                contrib[b, :, 2])
        f = open(full_name_folder_data+'/band_'+str(b), 'w')
        for i in range(len(bands.kpoints)):
            f.write('{0:10d} {1:15.8f}'.format(i, edif[i])+'\n')
        f.close()

    f = open(full_name_folder_data+'/band_structure_full', 'w')
    for j in range(len(bands.kpoints)):
        s = str(j)+3*' '
        for i in range(bands.nb_bands):
            f1 = open(full_name_folder_data+'/band_'+str(i))
            l1 = f1.readlines()
            f1.close()
            s += l1[i].rstrip().split()[1]+3*' '
        s += '\n'
        f.write(s)
    f.close()

    # Calculating the band gap
    print('band_'+str(cb_bottom)+'_'+str(cbm_pos), bands.bands[Spin.up][cb_bottom-1][cbm_pos])
    print('band_'+str(vb_top)+'_'+str(vbm_pos), bands.bands[Spin.up][vb_top-1][vbm_pos])
    print('Eg = ', min(bands.bands[Spin.up][cb_bottom-1])-max(bands.bands[Spin.up][vb_top-1]), 'eV')
    print('Eg_man = ', bands.bands[Spin.up][cb_bottom-1][cbm_pos]-bands.bands[Spin.up][vb_top-1][vbm_pos], 'eV') 
    print('band_'+str(cb_bottom)+'_min', min(bands.bands[Spin.up][cb_bottom-1]))
    print('band_'+str(vb_top)+'_max', max(bands.bands[Spin.up][vb_top-1]))       
    print('band_'+str(cb_bottom)+'_min_point', list(bands.bands[Spin.up][cb_bottom-1]).index(min(bands.bands[Spin.up][cb_bottom-1])))
    print('band_'+str(vb_top)+'_max_point', list(bands.bands[Spin.up][vb_top-1]).index(max(bands.bands[Spin.up][vb_top-1])))


    # style
    ax1.set_xlabel("k-points")
    ax1.set_ylabel(r"$E - E_f$   /   eV")
    ax1.grid()

    # fermi level at 0
    ax1.hlines(y=0, xmin=0, xmax=len(bands.kpoints), color="k", linestyle = '--', lw=1)

    # labels
    nlabs = len(labels)
    step = len(bands.kpoints) / (nlabs - 1)
    for i, lab in enumerate(labels):
        ax1.vlines(i * step, emin, emax, "k")
    ax1.set_xticks([i * step for i in range(nlabs)])
    ax1.set_xticklabels(labels)

    ax1.set_xlim(0, len(bands.kpoints))

    # Density of states
    # ----------------

    ax2.set_yticklabels([])
    ax2.grid()
    ax2.set_xlim(1e-4, 5)
    ax2.set_xticklabels([])
    ax2.hlines(y=0, xmin=0, xmax=5, color="k", lw=2)
    ax2.set_xlabel("Density of States", labelpad=28)

    # spd contribution
    ax2.plot(spd_dos[OrbitalType.s].densities[Spin.up],
             dosrun.tdos.energies - dosrun.efermi,
             "r-", label="s", lw=2)
    ax2.plot(spd_dos[OrbitalType.p].densities[Spin.up],
             dosrun.tdos.energies - dosrun.efermi,
             "g-", label="p", lw=2)
    ax2.plot(spd_dos[OrbitalType.d].densities[Spin.up],
             dosrun.tdos.energies - dosrun.efermi,
             "b-", label="d", lw=2)

    # total dos
    ax2.fill_between(dosrun.tdos.densities[Spin.up],
                     0,
                     dosrun.tdos.energies - dosrun.efermi,
                     color=(0.7, 0.7, 0.7),
                     facecolor=(0.7, 0.7, 0.7))

    ax2.plot(dosrun.tdos.densities[Spin.up],
             dosrun.tdos.energies - dosrun.efermi,
             color=(0.6, 0.6, 0.6),
             label="total DOS")

    edif1 = dosrun.tdos.energies - dosrun.efermi
    f = open(full_name_folder_data+'/dos_s_'+element, 'w')
    for i in range(len(spd_dos[OrbitalType.s].densities[Spin.up])):
        f.write('{0:15.8f} {1:15.8f}'.format(spd_dos[OrbitalType.s].densities[Spin.up][i], edif1[i])+'\n')
    f.close()

    f = open(full_name_folder_data+'/dos_p_'+element, 'w')
    for i in range(len(spd_dos[OrbitalType.p].densities[Spin.up])):
        f.write('{0:15.8f} {1:15.8f}'.format(spd_dos[OrbitalType.p].densities[Spin.up][i], edif1[i])+'\n')
    f.close()

    f = open(full_name_folder_data+'/dos_d_'+element, 'w')
    for i in range(len(spd_dos[OrbitalType.d].densities[Spin.up])):
        f.write('{0:15.8f} {1:15.8f}'.format(spd_dos[OrbitalType.d].densities[Spin.up][i], edif1[i])+'\n')
    f.close()

    f = open(full_name_folder_data+'/dos_tot', 'w')
    for i in range(len(dosrun.tdos.densities[Spin.up])):
        f.write('{0:15.8f} {1:15.8f}'.format(dosrun.tdos.densities[Spin.up][i], edif1[i])+'\n')
    f.close()

    # plot format style
    # -----------------
    ax2.legend(fancybox=True, shadow=True, prop={'size': 18})
    plt.subplots_adjust(wspace=0)
    plt.tight_layout()

    # plt.show()
    plt.savefig(folder+vasprun_bands.split('/')[-2]+".pdf", format='pdf', dpi=600)
    # print('bands.efermi ', bands.efermi)
    # print('dosrun.efermi', dosrun.efermi)
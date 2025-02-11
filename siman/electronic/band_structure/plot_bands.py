"""
Author: Kobernik Tatiana
"""

import numpy as np
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import re
from siman.header import printlog
import os
from collections import defaultdict
from matplotlib.lines import Line2D


def find_vbm_cbm(energies, valent_band_num):
    """
    This function is used to find the global maximum of the valence band
    and the global minimum of the conduction band.

    INPUT:
        - energies (numpy.ndarray) - array of arrays with energy values for all calculated bands
        - valent_band_num (int) - valence band number
    RETURN:
        - vbm (float) - energy value corresponding to the valence band maximum
        - cbm (float) - energy value corresponding to the conduction band minimum
    """
    vbm = energies[valent_band_num-1].max()
    cbm = energies[valent_band_num].min()
    return vbm, cbm


def lattice_from_vasprun(file_path):
    """
    The function parses lattice vectors from vasprun.xml

    INPUT:
        - file_path (str) - path to vasprun.xml file
    RETURN:
        - lattice (np.array) - matrix with lattice vectors
    """
    tree = ET.parse(file_path)
    root = tree.getroot()

    basis = root.find(".//structure/crystal/varray[@name='basis']")
    vectors = basis.findall("v")

    lattice = []
    for vector in vectors:
        vector_values = np.array(list(map(float, vector.text.split())))
        lattice.append(vector_values)

    return np.array(lattice)


def parse_kpoints_for_band(file_path):
    """
    The function parses k-points from vasprun.xml
    Used to plot band structure

    INPUT:
        - file_path (str) - path to vasprun.xml file
    RETURN:
        - kpoints (list of tuples) - list with coordinates of all kpoints at which energies were calculated
    """
    tree = ET.parse(file_path)
    root = tree.getroot()

    kpoints_list = root.find(".//varray[@name='kpointlist']")
    if kpoints_list is not None:
        kpoints = [tuple(map(float, v.text.split())) for v in kpoints_list.findall("v")]

    return kpoints


def parse_efermi(file_path):
    """
    The function parses fermi energy from vasprun.xml

    INPUT:
        - file_path (str) - path to vasprun.xml file
    RETURN:
        - efermi (int) - fermi energy value
    """
    tree = ET.parse(file_path)
    root = tree.getroot()

    efermi_tag = root.find(".//i[@name='efermi']")
    if efermi_tag is not None:
        efermi = float(efermi_tag.text.strip())

    return efermi


def parse_vasprun_total(file_path):
    """
    This function parses total energy values from vasprun.xml output file

    INPUT:
        - file_path (str) - path to vasprun.xml file
    RETURN:
        - kpoints (list of tuples) - list with coordinates of all kpoints at which energies were calculated
        - energies (numpy.ndarray) - array of arrays with energy values for all calculated bands
        - valent_band_num (int) - valence band number
    """
    tree = ET.parse(file_path)
    root = tree.getroot()

    energies = []

    kpoints = parse_kpoints_for_band(file_path)

    spin_set = root.find(".//set[@comment='spin 1']")
    if spin_set is not None:
        kpoint_sets = spin_set.findall("set[@comment]")
        if kpoint_sets:
            for kpoint_set in kpoint_sets:
                if 'kpoint' in kpoint_set.attrib['comment']:
                    values = [float(r.text.split()[0]) for r in kpoint_set.findall("r")]
                    energies.append(values)

    energies = np.array(energies).T if energies else np.array([])

    nelect_tag = root.find(".//i[@name='NELECT']")
    if nelect_tag is not None:
        nelect = float(nelect_tag.text.strip())
        valent_band_num = int(nelect / 2)
        vbm, cbm = find_vbm_cbm(energies, valent_band_num)

    energies = energies - vbm

    return kpoints, energies, valent_band_num


def parse_kpoints(file_path):
    """
    This function is used to perform parsing of KPOINTS file

    INPUT:
        - file_path (str) - path to KPOINTS file
    RETURN:
        - kpoints_list (list of tuples) - list with coordinates of all kpoints in file
        - kpoints_dict (dict) - dictionary that correlates the coordinates of points with their names,
                                (key = (tuple) coord, value = (str) name)
    """
    kpoints_dict = {}
    kpoints_list = []
    pattern = re.compile(r"([-+]?[0-9]*\.?[0-9]+)\s+([-+]?[0-9]*\.?[0-9]+)\s+([-+]?[0-9]*\.?[0-9]+)\s*!\s*(\S+)")

    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            match = pattern.search(line)
            if match:
                coords = tuple(map(float, match.groups()[:3]))
                label = match.group(4).strip()
                kpoints_dict[coords] = label
                kpoints_list.append(coords)

    return kpoints_list, kpoints_dict


def compute_k_distances(kpoints):
    """
    This function calculates the distance between k-points in 3D space and returns an array of x-axis values
    Used to plot band structure

    INPUT:
        - kpoints (list of tuples) - list with coordinates of the kpoints at which energies were calculated
    RETURN:
        - distances (np.array) - list with distances between kpoints
    """
    distances = [0]
    for i in range(1, len(kpoints)):
        dist = np.linalg.norm(np.array(kpoints[i]) - np.array(kpoints[i - 1]))
        distances.append(distances[-1] + dist)
    return np.array(distances)


def total_band(st_name, vasprun_path, kpoints_file_path, ylim=(-10, 10), method=None, vbm_cbm_marker=False,
                  file_name='', debug=False):
    """
    This function is used to build plot of the total electronic band structure
    It saves the plot to .png file

    INPUT:
        - st_name (str) - name of the structure which will be noted in title
        - vasprun_path (str) - path to the vasprun.xml file
        - kpoints_file_path (str) - path to the KPOINTS file for band structure
        - ylim (tuple of floats) - energy range of the band structure and DOS plots, units are eV
        - method (str) - name of the method that was used to find kpath, just for title
        - vbm_cbm_marker (bool) - indicator showing whether vbm and cbm points need to be marked on the plot
        - file_name (str) - global path to the file where the plot should be saved.
                            The default path = os.getcwd() + '/band_structures/' + f'{st_name}_{method}.png'
        - debug (bool) - indicator showing whether the plot should be displayed on the screen
    RETURN:
        None
    """

    kpoints, energies, valent_band_num = parse_vasprun_total(vasprun_path)
    kpoints_list, kpoints_dict = parse_kpoints(kpoints_file_path)

    if energies.size == 0:
        printlog("ERROR: energy array is empty. Check the XML file.")
        return

    kpath = compute_k_distances(kpoints)

    plt.figure(figsize=(8, 6), dpi=100)
    line_color = (0.2, 0.5, 0.8, 0.8)
    font = {'family': 'serif', 'size': 24}
    plt.rc('font', **font)

    for energy_band in energies:
        plt.plot(kpath, energy_band, color=line_color, linewidth=1)

    tick_positions = [kpath[i] for i, k in enumerate(kpoints) if k in kpoints_dict]
    tick_labels = [kpoints_dict[k] for k in kpoints if k in kpoints_dict]

    plt.xticks(tick_positions, tick_labels, fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel("k-points", fontsize=16)
    plt.ylabel("E - E$_f$ (eV)", fontsize=16)
    if method:
        plt.title(f"Total band structure for {st_name}\nMethod of k-path finding {method}\n", fontsize=19)
    else:
        plt.title(f"Total band structure for {st_name}\n", fontsize=19)

    for pos in tick_positions:
        plt.axvline(x=pos, color='black', linestyle='-', linewidth=0.7)

    plt.axhline(y=0, color='black', linestyle='-', linewidth=0.7)
    plt.grid(False)
    plt.xlim(min(kpath), max(kpath))
    plt.ylim(ylim)

    if vbm_cbm_marker and valent_band_num:
        vbm, cbm = find_vbm_cbm(energies, valent_band_num)

        vbm_x_id = np.where(energies[valent_band_num-1] == vbm)[0]
        for id in vbm_x_id:
            plt.scatter(kpath[id], vbm, color='r')

        cbm_x_id = np.where(energies[valent_band_num] == cbm)[0]
        for id in cbm_x_id:
            plt.scatter(kpath[id], cbm, color='g')

    plt.tight_layout()

    if method:
        name_png = f'{st_name}_{method}'
    else:
        name_png = f'{st_name}'
        
    if not file_name:
        file_name = os.getcwd() + '/band_structures/' + f'{name_png}.png'
    plt.savefig(file_name)

    if debug:
        plt.show()


def parse_ion_names(file_path):
    """
    This function analyzes all the ions in the structure for which the energy has been calculated.
    Working with vasprun.xml a file.
    Used in the projected_band function.

    INPUT:
        - file_path (str) - path to vasprun.xml file
    RETURN:
        - ion_names (list of str) - list with names of the ions in the order in which they are calculated in the .xml file
                                    'ion 1' = ion_names[0]
    """
    tree = ET.parse(file_path)
    root = tree.getroot()

    atominfo = root.find(".//atominfo")
    if atominfo is None:
        printlog("ERROR from parse_ion_names: Tag <atominfo> not found!")
        exit()

    atoms_array = atominfo.find(".//array[@name='atoms']")
    if atoms_array is None:
        printlog("ERROR from parse_ion_names: Tag <array name='atoms'> not found!")
        exit()

    ion_names = []
    for rc in atoms_array.findall(".//set/rc"):
        columns = rc.findall("c")
        if len(columns) >= 1:
            ion_names.append(columns[0].text.strip())

    return ion_names


def parse_nband(file_path):
    """
    The function parses number of calculated bands from vasprun.xml

    INPUT:
        - file_path (str) - path to vasprun.xml file
    RETURN:
        - nbands (int) - number of bands
    """
    tree = ET.parse(file_path)
    root = tree.getroot()

    nbands_tag = root.find(".//incar/i[@name='NBANDS']")

    if nbands_tag is not None:
        nbands = int(nbands_tag.text.strip())
    else:
        printlog("ERROR from parse_projected: NBANDS not found")
        exit()

    return nbands


def parse_vasprun_projected(file_path, ion_names):
    """
    The function parses the projected energy of different orbitals of all ions in the structure
    from vasprun.xml output file

    INPUT:
        - file_path (str) - path to vasprun.xml file
        - ion_names (list of str) - from parse_ion_names function - list with names of the ions
                                    in the order in which they are calculated in the .xml file
    RETURN:
        - kpoints (list of tuples) - list with coordinates of all kpoints at which energies were calculated
        - projected_bands (dict) - dictionary with information about the contribution of a certain ion to the energy of each orbital

        # Structure of projected_bands: {'spin 1': [ /for band 1/ :[ /for kpoint 1/ :{element: defaultdict(float, {'s': .., 'p': .., 'd': ..})},
                                                                    /for kpoint 2/ {...}, ... ], /for band 2/: [ ], .. ] }
    """
    tree = ET.parse(file_path)
    root = tree.getroot()

    nbands = parse_nband(file_path)
    kpoints = parse_kpoints_for_band(file_path)

    projected_array = root.find(".//calculation/projected/array")
    if projected_array is None:
        printlog("ERROR from parse_projected: tag projected -> array not found")
        exit()

    projected_bands = dict()

    spin_sets = [child for child in projected_array.findall(".//set")
                 if child.get("comment") and re.match(r"spin\s*\d+", child.get("comment"))]

    for spin in spin_sets:
        projected_bands[spin.get("comment")] = [[] for _ in range(nbands)]
        for kpoint_set in spin.findall("set"):
            for band_id, band_set in enumerate(kpoint_set.findall("set")):
                if band_id >= nbands:
                    break

                kpoint_data = {}
                for ion_idx, r in enumerate(band_set.findall("r")):
                    values = list(map(float, r.text.split()))
                    s_val = values[0]
                    p_val = sum(values[1:4])  # px, py, pz
                    d_val = sum(values[4:])  # dxy, dyz, dz2, dxz, x2-y2

                    ion_name = ion_names[ion_idx]
                    if ion_name not in kpoint_data:
                        kpoint_data[ion_name] = defaultdict(float)

                    kpoint_data[ion_name]['s'] += s_val
                    kpoint_data[ion_name]['p'] += p_val
                    kpoint_data[ion_name]['d'] += d_val

                projected_bands[spin.get("comment")][band_id].append(kpoint_data)

    return kpoints, projected_bands


def rgbline(ax, k, e, red, green, blue, alpha=20):
    """
    This function is used to provide colour for the line of the band structure plot
    depending on the contribution of different orbitals of the specific element.
    It is used in the 'projected_band' function.

    INPUT:
        - ax (matplotlib plot object) - object of the band structure plot
        - k (np.array) - from compute_k_distances func - list with distances between kpoints
        - e (list of floats) - list of energies corresponding to the k-points from the parameter 'k'
        - red (<class 'numpy.ndarray'>) - contribution from s orbitals
        - green (<class 'numpy.ndarray'>) - contribution from p orbitals
        - blue (<class 'numpy.ndarray'>) - contribution from d orbitals
        - alpha (float) - transparancy
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
    a = np.ones(nseg, np.float64) * alpha
    lc = LineCollection(seg, colors=list(zip(r, g, b, a)), linewidth=1.5)
    ax.add_collection(lc)


def projected_band(st_name, vasprun_path, kpoints_file_path, element, ylim=(-10, 10), method=None, file_name='',
                   debug=False):
    """
    This function is used to build plot of the orbital projected electronic band structure for some element
    It saves the plot to .png file

    INPUT:
        - st_name (str) - name of the structure which will be noted in title
        - vasprun_path (str) - path to the vasprun.xml file
        - kpoints_file_path (str) - path to the KPOINTS file for band structure
        - element (str) - name of element for which the projected band structure will be built
        - ylim (tuple of floats) - energy range of the band structure and DOS plots, units are eV
        - method (str) - name of the method that was used to find kpath, just for title
        - file_name (str) - global path to the file where the plot should be saved.
                            The default path = os.getcwd() + '/band_structures/' + f'{st_name}_{method}.png'
        - debug (bool) - indicator showing whether the plot should be displayed on the screen
    RETURN:
        None
    """
    nbands = parse_nband(vasprun_path)
    kpoints_list, kpoints_dict = parse_kpoints(kpoints_file_path)

    ion_names = parse_ion_names(vasprun_path)
    if element not in ion_names:
        print("\n\nERROR: element for projected band structure not in structure's ions list\n\n")
        exit()

    kpoints, pbands = parse_vasprun_projected(vasprun_path, ion_names)
    kpoints, energies, valent_band_num = parse_vasprun_total(vasprun_path)

    if not nbands:
        printlog("ERROR: energy array is empty. Check the XML file.")
        exit()

    kpath = compute_k_distances(kpoints)

    font = {'family': 'serif', 'size': 20}
    plt.rc('font', **font)
    fig = plt.figure(figsize=(8, 6), dpi=100)
    ax1 = plt.subplot()
    ax1.set_ylim(ylim[0], ylim[1])
    ax1.set_xlim(min(kpath), max(kpath))
    if method:
        fig.suptitle(f"Orbital projected band structure for {st_name}\nElement {element}\nMethod of k-path finding {method}", fontsize=19)
    else:
        fig.suptitle(f"Orbital projected band structure for {element}\nElement {element}", fontsize=19)

    contrib = np.zeros((nbands, len(kpoints), 3))
    for b in range(nbands):
        for k in range(len(kpoints)):
            sc = pbands['spin1'][b][k][element]["s"] ** 2
            pc = pbands['spin1'][b][k][element]["p"] ** 2
            dc = pbands['spin1'][b][k][element]["d"] ** 2
            tot = sc + pc + dc
            if tot != 0.0:
                contrib[b, k, 0] = sc / tot
                contrib[b, k, 1] = pc / tot
                contrib[b, k, 2] = dc / tot

    # plot bands using rgb mapping
    for b in range(nbands):
        rgbline(ax1,
                kpath,
                energies[b],
                contrib[b, :, 0],
                contrib[b, :, 1],
                contrib[b, :, 2])


    ax1.set_xlabel("k-points", fontsize=16)
    ax1.set_ylabel("E - E$_f$ (eV)", fontsize=16)
    ax1.grid()

    # fermi level at 0
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.7)

    tick_positions = [kpath[i] for i, k in enumerate(kpoints) if k in kpoints_dict]
    tick_labels = [kpoints_dict[k] for k in kpoints if k in kpoints_dict]

    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(tick_labels, fontsize=16)

    for pos in tick_positions:
        ax1.vlines(pos, ylim[0], ylim[1], color="black", linewidth=0.7)

    legend_elements = [
        Line2D([0], [0], color='r', lw=2, label='s'),
        Line2D([0], [0], color='g', lw=2, label='p'),
        Line2D([0], [0], color='b', lw=2, label='d')
    ]

    # ax1.legend(handles=legend_elements, bbox_to_anchor=(-0.09, 1), loc='upper right', prop={'size': 16})
    ax1.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1), prop={'size': 16})
    # ax1.legend(handles=legend_elements, bbox_to_anchor=(-0.25, 1), loc='upper left', prop={'size': 16})

    plt.subplots_adjust(wspace=0)
    plt.grid(False)
    plt.tight_layout()

    if method:
        name_png = f'{st_name}_{method}'
    else:
        name_png = f'{st_name}'

    if not file_name:
        file_name = os.getcwd() + '/band_structures/' + f'{name_png}.png'
    plt.savefig(file_name)

    path = os.getcwd() + '/band_structures/' + f'{name_png}_proj_{element}.png'
    plt.savefig(path)

    if debug:
        plt.show()


def plot_bands(st_name, vasprun_path, kpoints_file_path, element=None, ylim=(-10, 10), mode='total',
               vbm_cbm_marker=False, file_name='', method="", debug=False):
    """
    The function is a wrapper for total_band and projected_band functions
    It is used to plot the total or projected band structure

    INPUT:
        - st_name (str) - name of the structure which will be noted in title
        - vasprun_path (str) - path to the vasprun.xml file
        - kpoints_file_path (str) - path to the KPOINTS file for band structure
        - element (str) - name of element for which the projected band structure will be built
        - ylim (tuple of floats) - energy range of the band structure and DOS plots, units are eV
        - vbm_cbm_marker (bool) - indicator showing whether vbm and cbm points need to be marked on the plot
        - file_name (str) - global path to the file where the plot should be saved.
                            The default path = os.getcwd() + '/band_structures/' + f'{st_name}_{method}.png'
        - method (str) - name of the method that was used to find kpath (just for title)
        - debug (bool) - indicator showing whether the plot should be displayed on the screen
    RETURN:
        None
    TODO:
        Some improvements
    """
    if mode == 'projected':
        if not element:
            print('Please, enter the element for projected')
        else:
            projected_band(st_name, vasprun_path, kpoints_file_path, element, method=method, ylim=ylim,
                           file_name=file_name, debug=debug)
    elif mode == "total":
        total_band(st_name, vasprun_path, kpoints_file_path, ylim=ylim, method=None, vbm_cbm_marker=vbm_cbm_marker,
                   file_name=file_name, debug=debug)
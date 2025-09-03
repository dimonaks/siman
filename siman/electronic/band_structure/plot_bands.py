"""
Author: Kobernik Tatiana
"""

import numpy as np
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use('Agg')  # for matsolver backend
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import re
from siman.header import printlog
import os
from collections import defaultdict
from matplotlib.lines import Line2D


def _fix_margins(fig=None, *, left=0.15, right=0.82, top=0.90, bottom=0.14):
    """
    The function is designed to align the size of graphs for total and projected band structures. 
    Fixes the same fields around the graph.

    INPUT:
        - fig (matplotlib.figure.Figure) - the object of the matplotlib shape that needs to be indented
        - left (float) - left margin (0–1, fraction of width)
        - right (float) - right margin (0–1, fraction of width); a place for a legend/void (0.80 ≈ 20% of the width)
        - top (float) - top margin (0–1, fraction of height)
        - bottom (float) - bottom margin (0–1, fraction of height) 
    RETURN:
        None
    """
    fig.subplots_adjust(left=left, right=right,
                        bottom=bottom, top=top)
    

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


def parse_vasprun_kpoints(file_path):
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


def parse_total_energies(file_path, valent_band_num):
    """
    The function parses total energy's array from vasprun.xml

    INPUT:
        - file_path (str) - path to vasprun.xml file
        - valent_band_num (int) - valence band number
    RETURN:
        - energies (np.array) - array of arrays with energy values for all calculated bands;
                                the energies are shifted so that the valence band maximum is at 0 eV;
                                the energies are rounded to 4 decimal places;
    """
    tree = ET.parse(file_path)
    root = tree.getroot()

    energies = []

    spin_set = root.find(".//set[@comment='spin 1']")
    if spin_set is not None:
        kpoint_sets = spin_set.findall("set[@comment]")
        if kpoint_sets:
            for kpoint_set in kpoint_sets:
                if 'kpoint' in kpoint_set.attrib['comment']:
                    values = [float(r.text.split()[0]) for r in kpoint_set.findall("r")]
                    energies.append(values)

    energies = np.array(energies).T if energies else np.array([])

    if energies.size == 0:
        printlog("ERROR: energy array is empty. Check the XML file", imp = 'y')
        return

    vbm, _ = find_vbm_cbm(energies, valent_band_num)
    energies = energies - vbm
    energies = np.round(energies, 4)

    return energies


def parse_valense_band_num(file_path):
    """
    The function parses valence band number from vasprun.xml

    INPUT:
        - file_path (str) - path to vasprun.xml file
    RETURN:
        - valent_band_num (int) - valence band number
    """
    tree = ET.parse(file_path)
    root = tree.getroot()

    nelect_tag = root.find(".//i[@name='NELECT']")
    if nelect_tag is not None:
        nelect = float(nelect_tag.text.strip())
        valent_band_num = int(nelect / 2)
    
    return valent_band_num


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


def prepare_kpoints(vasprun_path, kpoints_file_path):
    """
    INPUT:
        - vasprun_path (str) - path to the vasprun.xml file
        - kpoints_file_path (str) - path to the KPOINTS file for band structure
    RETURN:
        kpath, high_symmetric_kpoints, high_symmetric_kpoints_labels, valent_band_num
        - kpath (np.array) - list with linear distances between points
        - high_symmetric_kpoints (list of floats) - list with elements from 'kpath' corresponding to highly symmetrical k-points
        - high_symmetric_kpoints_labels (list of str) - list with names of highly symmetrical k-points
        - valent_band_num (int) - valence band number
    """

    all_kpoints = parse_vasprun_kpoints(vasprun_path) # all_kpoints - list of tuples with coordinates in 3D dimension; all range
    valent_band_num = parse_valense_band_num(vasprun_path)
    _ , kpoints_names = parse_kpoints(kpoints_file_path) # kpoints_names - dictionary
                                                                            
    kpath = compute_k_distances(all_kpoints)

    high_symmetric_kpoints = [kpath[i] for i, k in enumerate(all_kpoints) if k in kpoints_names]
    high_symmetric_kpoints_labels = [kpoints_names[k] for k in all_kpoints if k in kpoints_names]

    return kpath, high_symmetric_kpoints, high_symmetric_kpoints_labels, valent_band_num


def total_band(energies, kpath, high_symmetric_kpoints, high_symmetric_kpoints_labels, valent_band_num, st_name, ylim=(-10, 10), method=None, vbm_cbm_marker=False,
                  file_name='', title=''):
    """
    This function is used to build plot of the total electronic band structure
    It saves the plot to .png file

    INPUT:
        - energies (numpy.ndarray) - array of arrays with energy values for all calculated bands
        - kpath (np.array) - list with linear distances between points
        - high_symmetric_kpoints (list of floats) - list with elements from 'kpath' corresponding to highly symmetrical k-points
        - high_symmetric_kpoints_labels (list of str) - list with names of highly symmetrical k-points
        - valent_band_num (int) - valence band number
        - st_name (str) - name of the structure which will be noted in title
        - ylim (tuple of floats) - energy range of the band structure, units are eV
        - method (str) - name of the method that was used to find kpath, just for title
        - vbm_cbm_marker (bool) - indicator showing whether vbm and cbm points need to be marked on the plot
        - file_name (str) - global path to the file where the plot should be saved.
                            The default path = os.getcwd() + '/band_structures/' + f'{st_name}_{method}.png'
        - title (str) - manual title for figure
    RETURN:
        None
    """
    fig, ax = plt.subplots(figsize=(10, 6), dpi=100)
    line_color = (0.2, 0.5, 0.8, 0.8)
    font = {'family': 'serif', 'size': 24}
    plt.rc('font', **font)

    for energy_band in energies:
        ax.plot(kpath, energy_band, color=line_color, linewidth=1)

    tick_positions = high_symmetric_kpoints
    tick_labels = high_symmetric_kpoints_labels

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=14)
    ax.tick_params(axis='y', labelsize=14)
    ax.set_xlabel("k-points", fontsize=16)
    ax.set_ylabel("E - E$_f$ (eV)", fontsize=16)
    if title:
        ax.set_title(title, fontsize=19)
    else:
        if method:
            ax.set_title(f"Total band structure. Method {method}", fontsize=19)
        else:
            ax.set_title(f"Total band structure", fontsize=19)

    ax.grid(False)
    ax.set_xlim(min(kpath), max(kpath))
    ax.set_ylim(ylim)

    for pos in tick_positions:
        ax.axvline(
            x=pos,
            color="black",
            linestyle="-",
            linewidth=0.7,
            zorder=5
        )

    ax.axhline(
        y=0,
        color="black",
        linestyle="-",
        linewidth=0.7,
        zorder=5
    )

    if vbm_cbm_marker and valent_band_num:
        vbm, cbm = find_vbm_cbm(energies, valent_band_num)

        vbm_x_id = np.where(energies[valent_band_num-1] == vbm)[0]
        for id in vbm_x_id:
            ax.scatter(kpath[id], vbm, color='r')

        cbm_x_id = np.where(energies[valent_band_num] == cbm)[0]
        for id in cbm_x_id:
            ax.scatter(kpath[id], cbm, color='g')

    _fix_margins(fig)

    if method:
        name_png = f'{st_name}_{method}'
    else:
        name_png = f'{st_name}'
        
    if not file_name:
        file_name = os.getcwd() + '/band_structures/' + f'{name_png}.png'
    
    folder = os.path.dirname(file_name)
    if folder != '' and not os.path.exists(folder):
        os.makedirs(folder, exist_ok=True)

    fig.savefig(file_name, dpi=100)


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
        printlog("ERROR from parse_ion_names: Tag <atominfo> not found!", imp = 'y')
        exit()

    atoms_array = atominfo.find(".//array[@name='atoms']")
    if atoms_array is None:
        printlog("ERROR from parse_ion_names: Tag <array name='atoms'> not found!", imp = 'y')
        exit()

    ion_names = []
    for rc in atoms_array.findall(".//set/rc"):
        columns = rc.findall("c")
        if len(columns) >= 1:
            ion_names.append(columns[0].text.strip())

    return ion_names


def parse_vasprun_projected(file_path, ion_names, nbands, selected_k_indices=None):
    """
    The function parses the projected energy of different orbitals of all ions in the structure
    from vasprun.xml output file with optional k-point thinning

    INPUT:
        - file_path (str) - path to vasprun.xml file
        - ion_names (list of str) - from parse_ion_names function - list with names of the ions
                                    in the order in which they are calculated in the .xml file
        - nbands (int) - number of bands
        - selected_k_indices (list of int or np.ndarray) - indices of k-points to keep after thinning
    RETURN:
        - kpoints (list of tuples) - list with coordinates of all kpoints at which energies were calculated
        - projected_bands (dict) - dictionary with information about the contribution of a certain ion to the energy of each orbital

        # Structure of projected_bands: {'spin 1': [ /for band 1/ :[ /for kpoint 1/ :{element: defaultdict(float, {'s': .., 'p': .., 'd': ..})},
                                                                    /for kpoint 2/ {...}, ... ], /for band 2/: [ ], .. ] }
    """
    tree = ET.parse(file_path)
    root = tree.getroot()

    projected_array = root.find(".//calculation/projected/array")
    if projected_array is None:
        printlog("ERROR from parse_projected: tag projected -> array not found", imp='y')
        exit()

    projected_bands = {}

    spin_sets = [
        child for child in projected_array.findall(".//set")
        if child.get("comment") and re.match(r"spin\s*\d+", child.get("comment"))
    ]

    for spin in spin_sets:
        spin_name = spin.get("comment")
        kpoint_sets = list(spin.findall("set"))
        n_kpoints = len(kpoint_sets)

        if selected_k_indices is None:
            selected_k_indices = np.arange(n_kpoints)

        projected_bands[spin_name] = [[] for _ in range(nbands)]

        for k_idx in selected_k_indices:
            if k_idx >= n_kpoints:
                continue  # protection against index out of range

            kpoint_set = kpoint_sets[k_idx]

            for band_id, band_set in enumerate(kpoint_set.findall("set")):
                if band_id >= nbands:
                    break

                kpoint_data = {}
                for ion_idx, r in enumerate(band_set.findall("r")):
                    values = list(map(float, r.text.split()))
                    s_val = round(values[0], 4)
                    p_val = round(sum(values[1:4]), 4)
                    d_val = round(sum(values[4:]), 4)

                    ion_name = ion_names[ion_idx]
                    if ion_name not in kpoint_data:
                        kpoint_data[ion_name] = {'s': 0.0, 'p': 0.0, 'd': 0.0}

                    kpoint_data[ion_name]['s'] += s_val
                    kpoint_data[ion_name]['p'] += p_val
                    kpoint_data[ion_name]['d'] += d_val

                projected_bands[spin_name][band_id].append(kpoint_data)

    return projected_bands


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


def projected_band(energies, project_energies, kpath, high_symmetric_kpoints, high_symmetric_kpoints_labels, st_name, element, 
                           ylim=(-10, 10), valent_band_num=None, vbm_cbm_marker=None, method=None, file_name=None, title=''):
    """
    This function is used to build plot of the orbital projected electronic band structure for some element
    It saves the plot to .png file

    INPUT:
        - energies (numpy.ndarray) - array of arrays with energy values for all calculated bands
        - project_energies (dict) - from parse_vasprun_projected function - dictionary with information about 
                                                the contribution of a certain ion to the energy of each orbital
        - kpath (np.array) - list with linear distances between points
        - high_symmetric_kpoints (list of floats) - list with elements from 'kpath' corresponding to highly symmetrical k-points
        - high_symmetric_kpoints_labels (list of str) - list with names of highly symmetrical k-points
        - st_name (str) - name of the structure which will be noted in title
        - element (str) - name of element for which the projected band structure will be built
        - ylim (tuple of floats) - energy range of the band structure and DOS plots, units are eV
        - valent_band_num (int) - valence band number
        - vbm_cbm_marker (bool) - indicator showing whether vbm and cbm points need to be marked on the plot
        - method (str) - name of the method that was used to find kpath, just for title
        - file_name (str) - global path to the file where the plot should be saved.
                            The default path = os.getcwd() + '/band_structures/' + f'{st_name}_{method}.png'
        - title (str) - manual title for figure
    RETURN:
        None
    """
    
    nbands = len(energies)

    font = {'family': 'serif', 'size': 24}
    plt.rc('font', **font)
    fig = plt.figure(figsize=(10, 6), dpi=100)
    ax1 = plt.subplot()
    ax1.set_ylim(ylim[0], ylim[1])
    ax1.set_xlim(min(kpath), max(kpath))

    if title:
        fig.suptitle(title, fontsize=19)
    else:
        if method:
            fig.suptitle(f"Orbital projected band structure for {element}\nMethod of k-path finding {method}", fontsize=19)
        else:
            fig.suptitle(f"Orbital projected band structure for {element}", fontsize=19)

    contrib = np.zeros((nbands, len(kpath), 3))
    for b in range(nbands):
        for k in range(len(kpath)):
            sc = project_energies['spin1'][b][k][element]["s"] ** 2
            pc = project_energies['spin1'][b][k][element]["p"] ** 2
            dc = project_energies['spin1'][b][k][element]["d"] ** 2
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

    tick_positions = high_symmetric_kpoints
    tick_labels = high_symmetric_kpoints_labels

    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(tick_labels, fontsize=16)

    for pos in tick_positions:
        ax1.vlines(pos, ylim[0], ylim[1], color="black", linewidth=0.7)

    legend_elements = [
        Line2D([0], [0], color='r', lw=2, label='s'),
        Line2D([0], [0], color='g', lw=2, label='p'),
        Line2D([0], [0], color='b', lw=2, label='d')
    ]

    if vbm_cbm_marker and valent_band_num:
        vbm, cbm = find_vbm_cbm(energies, valent_band_num)

        vbm_x_id = np.where(energies[valent_band_num-1] == vbm)[0]
        for id in vbm_x_id:
            ax1.scatter(kpath[id], vbm, color='black', zorder=10)

        cbm_x_id = np.where(energies[valent_band_num] == cbm)[0]
        for id in cbm_x_id:
            ax1.scatter(kpath[id], cbm, color='black', zorder=10)

    ax1.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.15, 1), prop={'size': 16})

    plt.subplots_adjust(wspace=0)
    plt.grid(False)

    _fix_margins(fig)

    if method:
        name_png = f'{st_name}_{method}'
    else:
        name_png = f'{st_name}'

    if not file_name:
        file_name = os.getcwd() + '/band_structures/' + f'{name_png}_proj_{element}.png'
    
    folder = os.path.dirname(file_name)
    if folder != '' and not os.path.exists(folder):
        os.makedirs(folder, exist_ok=True)

    plt.savefig(file_name, dpi=100)


def thinning_total(kpath, high_symmetric_kpoints, energies, thinning_coeff):
    """
    The function is designed to thin out points in the kpoint and total energies arrays. 
    It is necessary to reduce the volume of files while saving band structures data (uses in matsolver).

    INPUT:
        - kpath (np.array) - list with linear distances between points
        - high_symmetric_kpoints (list of floats) - list with elements from 'kpath' corresponding to highly symmetrical k-points
        - energies (numpy.ndarray) - array of arrays with energy values for all calculated bands
        - thinning_coeff (int) - the thinning coefficient, the number of k points between highly 
                                    symmetrical ones, will be reduced by this number of times
    RETURN:
        - new_kpath (np.array) - thinned kpath array
        - new_energies (numpy.ndarray) - thinned energies array
        - selected_indices (np.array) - indices of k-points that were kept after thinning (uses in projected band structure)
    """
    high_symmetric_kpoints_idx = []
    i, j = 0, 0

    while i < len(kpath) and j < len(high_symmetric_kpoints):
        if kpath[i] == high_symmetric_kpoints[j]:
            high_symmetric_kpoints_idx.append(i)
            j += 1
        i += 1

    selected_indices = []

    start = 0
    for finish in high_symmetric_kpoints_idx:
        if finish - start > 1:
            idx = list(range(start, finish, thinning_coeff))
        else:
            idx = list(range(start, finish))

        selected_indices.extend(idx)
        start = finish

    if start < len(kpath):
        idx = list(range(start, len(kpath)))
        selected_indices.extend(idx)

    selected_indices = np.array(sorted(set(selected_indices)), dtype=int)

    new_kpath = np.asarray(kpath)[selected_indices]
    new_energies = energies[:, selected_indices]

    return new_kpath, new_energies, selected_indices


def plot_bands(st_name, vasprun_path, kpoints_file_path, thinning_coeff=1, element=None, ylim=(-15, 10), mode='total',
               vbm_cbm_marker=False, file_name='', method=None, title=''):
    """
    The function is a wrapper for total_band and projected_band functions
    It is used to plot the total or projected band structure

    INPUT:
        - st_name (str) - name of the structure which will be noted in title
        - vasprun_path (str) - path to the vasprun.xml file
        - kpoints_file_path (str) - path to the KPOINTS file for band structure
        - thinning_coeff (int) - the thinning coefficient, the number of k points between highly 
                                    symmetrical ones, will be reduced by this number of times
        - element (str) - name of element for which the projected band structure will be built
        - ylim (tuple of floats) - energy range of the band structure and DOS plots, units are eV
        - mode (str) - shows type of the band structure (total or projected)
        - vbm_cbm_marker (bool) - indicator showing whether vbm and cbm points need to be marked on the plot
        - file_name (str) - global path to the file where the plot should be saved.
                            The default path = os.getcwd() + '/band_structures/' + f'{st_name}_{method}.png'
        - method (str) - name of the method that was used to find kpath (just for title)
        - title (str) - manual title for figure
    RETURN:
        None
    """
    if thinning_coeff < 1:
        printlog('Thinning coefficient must be >= 1')
        return None
    
    kpath, high_symmetric_kpoints, high_symmetric_kpoints_labels, valent_band_num  = prepare_kpoints(vasprun_path, kpoints_file_path)
    energies = parse_total_energies(vasprun_path, valent_band_num)
    full_len_energies = len(energies)

    selected_indices = None
    if thinning_coeff > 1:
        kpath, energies, selected_indices = thinning_total(kpath, high_symmetric_kpoints, energies, thinning_coeff)

    if mode == 'projected':
        if not element:
            printlog('Please, enter the element for projected')
        else:
            ion_names = parse_ion_names(vasprun_path)
            if element not in ion_names:
                printlog("\n\nERROR: element for projected band structure not in structure's ions list\n\n")
                exit()

            projected_energies = parse_vasprun_projected(vasprun_path, ion_names, full_len_energies, selected_k_indices=selected_indices)

            projected_band(energies, projected_energies,
                           kpath, high_symmetric_kpoints, high_symmetric_kpoints_labels, st_name, element, ylim=ylim,
                           vbm_cbm_marker=vbm_cbm_marker, valent_band_num=valent_band_num, file_name=file_name, title=title)
    elif mode == "total":
        total_band(energies, kpath, high_symmetric_kpoints, high_symmetric_kpoints_labels, valent_band_num, st_name, ylim=ylim, method=method, vbm_cbm_marker=vbm_cbm_marker,
                            file_name=file_name, title=title)
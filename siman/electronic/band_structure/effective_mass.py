"""
Author: Kobernik Tatiana
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.io.vasp import BSVasprun
from siman.database import read_database
from siman.header import _update_configuration, db
from siman.electronic.band_structure.plot_bands import lattice_from_vasprun
from siman.header import printlog


def generate_kpoints_for_mass(kp_cart, vasprun_band_path, h, flag='tensor', direct=None, debug=None): #method from emc.py
    """
    The function creates a list of k-points for calculating effective masses.
    Supports tensor calculation (flag='tensor') and unit value calculation mode in a specific direction in the reverse space (flag='direct')

    INPUT:
        - kp_cart (list) - list with cartesian coordinates of k-point (extreme) in which you want to calculate effective mass
        - vasprun_band_path (str) - path to the vasprun.xml file with results of band structure calculations
        - h (float) - the step between points in calculations
        - flag (str) - specify the type of calculation ('tensor'/'direct')
        - direct (tuple) - tuple with vector which specify the directions of calculations (for flag='direct')
        - debug (bool) - if True, shows the obtained k-points
    RETURN:
        k_path (list of tuples) - list with coordinates of k-points in the obtained explicit k-path in the format for siman set
                                    Example: [19 */number of points/*, (1, 0.0, 0.0, 0.0), (1, 0.0, 0.00435, 0.00435), ...]
    """
    v_disp = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 1, 0],
        [1, -1, 0],
        [1, 0, 1],
        [1, 0, -1],
        [0, 1, 1],
        [0, 1, -1],
    ]

    L = lattice_from_vasprun(file_path=vasprun_band_path)

    a1 = np.array(L[0])
    a2 = np.array(L[1])
    a3 = np.array(L[2])

    volume = np.dot(a1, np.cross(a2, a3))

    b1 = 2 * np.pi * np.cross(a2, a3) / volume
    b2 = 2 * np.pi * np.cross(a3, a1) / volume
    b3 = 2 * np.pi * np.cross(a1, a2) / volume

    reciprocal_lattice = np.array([b1, b2, b3])

    kps = [kp_cart]

    if flag=='tensor':
        for v in v_disp:
            kps.append((h * np.array(v) + kp_cart))
            kps.append((- h * np.array(v) + kp_cart))

    if flag == 'direct' and direct:
        kps.append((h * np.array(direct) + kp_cart))
        kps.append((- h * np.array(direct) + kp_cart))

    kps_rec = np.dot(kps, np.linalg.inv(reciprocal_lattice))

    if debug:
        print("\n\nResults k-points_mass in cart:")
        for i in range(len(kps)):
            print(' %5.6f %5.6f %5.6f  1' % (kps[i][0], kps[i][1], kps[i][2]))

        print("\n\nResults k-points_mass in rec:")
        for i in range(len(kps)):
            print(' %5.6f %5.6f %5.6f  1' % (kps_rec[i][0], kps_rec[i][1], kps_rec[i][2]))

    return kpath_mass_to_siman(np.array(kps_rec))


def kpath_mass_to_siman(path_list):
    """
    The function is used to convert the coordinates of k-points from np.array into a format for a set
    Used in effective mass calculations, for explicit KPOINTS file

    INPUT:
        - path_list (np.array) - array with coordinates of k-points
    RETURN:
        k_path (list of tuples) - list with coordinates of k-points in the obtained explicit k-path in the format for siman set
                                    Example: [19 */number of points/*, (1, 0.0, 0.0, 0.0), (1, 0.0, 0.00435, 0.00435), ...]
    """
    k_path = [len(path_list)]
    for point in path_list:
        elem = (1, *point)
        k_path.append(elem)
    return k_path


def calc_num_diriv(bs, flag, band, h):
    """
    The function is used to obtain effective masses by calculating the numerical derivative over 3 points,
    the accuracy O(h^2)

    INPUT:
        - bs <class 'pymatgen.electronic_structure.bandstructure.BandStructureSymmLine'> - class with information about non-self-consit calc;
        can be obtaind from v = BSVasprun(effective_mass_vasprun.xml), bs = v.get_band_structure()
        - flag (str) - specify the type of calculation ('tensor'/'direct')
        - band (int) - number of band for which tou want to calculate effective mass
        - h (float) - the step between points in calculations
    RETURN:
        For 'tensor': - m (np.array) - matrix 3x3 with tensor of effective mass, the values NOT reversed; used to obtain reversed eigenvalues
        For 'direct': - m (np.float64) - value of effective mass expressed in the masses of an electron
    """
    spin = list(bs.bands.keys())[0]  # using spin up
    en = [e for e in bs.bands[spin][band - 1] / 27.21]  # and convert to Hartree units

    h = h * 0.529177  # converted to Bohr-1 units

    if flag == 'tensor':
        m = np.ones((3, 3))
        for j, i in enumerate([1, 3, 5]):
            deriv = (en[i + 1] - 2 * en[0] + en[i]) / h ** 2  # second order derivative with three-point method
            m[j][j] = deriv

        for jk, i in zip([(0, 1), (0, 2), (1, 2)], [7, 11, 15]):
            deriv = (en[i] + en[i + 1] - en[i + 2] - en[
                i + 3]) / 4 / h ** 2  # second order derivative with four-point method
            m[jk[0]][jk[1]] =  deriv
            m[jk[1]][jk[0]] =  deriv

    elif flag == 'direct':
        dir = (en[1] - 2 * en[0] + en[2]) / h ** 2
        m = 1 / dir

    return m


def eff_mass_value(band, h, folder_mass, vasprun_name='1.vasprun.xml', kpt_extr=None, flag='tensor', direct=None, debug=False):
    """
    The function is used to obtain effective masses by calculating the numerical derivative over 3 points,
    the accuracy O(h^2)

    INPUT:
        - band (int) - number of band for which tou want to calculate effective mass
        - h (float) - the step between points in calculations
        - folder_mass (str) - path to the folder with calculations of effective mass
        - vasprun_name (str) - name of vasprun.xml output file in 'folder_mass'
        - kp_extr (list) - list withcoordinates of k-point (extreme) in which you want to calculate effective mass
        - flag (str) - specify the type of calculation ('tensor'/'direct')
        - direct (tuple) - tuple with vector which specify the directions of calculations (for flag='direct')
        - debug (bool) - if True, shows the obtained effective masses

    RETURN:
        For 'tensor': - mass (np.array) - matrix 3x3 with tensor of effective mass, the values NOT reversed; used to obtain reversed eigenvalues
                      - eigenval_mass (list) - list with tensor eigenvalues

        For 'direct': - mass (np.float64) - value of effective mass expressed in the masses of an electron
    TODO:
        Implement without pymatgen
    """
    v = BSVasprun(f'{folder_mass}{vasprun_name}', occu_tol=1e-08)
    bs = v.get_band_structure(kpoints_filename=f'{folder_mass}KPOINTS')

    mass = calc_num_diriv(bs, flag, band, h)
    eigenval_mass = []

    if debug:
        if flag=='direct' and direct:
            print(f'Effective mass in point {kpt_extr} for direction {direct} = {round(mass, 3)} m_e; band = {band}')

        elif flag=='tensor':
            print("\ntensor of masses")
            for row in mass:
                print(" ".join(f"{num:.5f}  " for num in row))
            v, w = np.linalg.eig(mass)
            eigenval_mass = [round(1 / val,3) for val in v]

            print("\neigenvals:\n", eigenval_mass)
            print("\neigenvectors:\n")
            for i in range(3):
                for j in range(3):
                    w[i][j] = round(w[i][j], 3)
            print(w)

        else:
            printlog("ERROR from 'eff_mass_value': flag is wrong")
            return 0

    return mass, eigenval_mass


def eff_mass_from_interpol(band, h, folder_mass, vasprun_name='1.vasprun.xml', kpt_extr=None, flag='tensor', direct=None, debug=False):
    """
    The function is used to obtain effective masses by interpolations of energies values

    INPUT:
        - band (int) - number of band for which tou want to calculate effective mass
        - h (float) - the step between points in calculations
        - folder_mass (str) - path to the folder with calculations of effective mass
        - vasprun_name (str) - name of vasprun.xml output file in 'folder_mass'
        - kp_extr (list) - list withcoordinates of k-point (extreme) in which you want to calculate effective mass
        - flag (str) - specify the type of calculation ('tensor'/'direct')
        - direct (tuple) - tuple with vector which specify the directions of calculations (for flag='direct')
        - debug (bool) - if True, shows the obtained effective masses
    RETURN:
        For 'tensor': - result_tensor (np.array) - matrix 3x3 with tensor of effective mass
        For 'direct': - second_dir (np.float64) - value of effective mass expressed in the masses of an electron
    TODO:
        Some improvements
    """
    from scipy.interpolate import CubicSpline

    try:
        v = BSVasprun(f'{folder_mass}{vasprun_name}', occu_tol=1e-08)
        bs = v.get_band_structure(kpoints_filename=f'{folder_mass}KPOINTS')
        spin = list(bs.bands.keys())[0]  # using spin up
        en = [e for e in bs.bands[spin][band - 1] / 27.21]  # and convert to Hartree units


        h = h * 0.529177

        if flag == 'tensor':
            result_tensor = []

            for i in range(1, 6, 2):
                y = np.array([en[i], en[0], en[i+1]])
                x = np.array([i*h for i in range(len(y))])

                print(f"manual dir {i}: ", 1 / ((y[0] - 2 * y[1] + y[2]) / h ** 2)) #check the manual dir calk

                cs = CubicSpline(x, y)
                second_dir = cs.derivative(2)(x[1])
                result_tensor.append(1 / second_dir)
            print(result_tensor)
            return result_tensor

        if flag == 'direct':

            y = np.array([en[1], en[0], en[2]])
            x = np.array([i * h for i in range(len(y))])
            print(f"manual dir: ", 1 / ((y[0] - 2 * y[1] + y[2]) / h ** 2))  # check the manual dir calk

            cs = CubicSpline(x, y)
            second_dir = 1 / cs.derivative(2)(x[1])

            print('interpol mass:', second_dir)
            return second_dir

            x_fine = np.linspace(0, 2, 100)
            y_fine = cs(x_fine)

            plt.plot(x_fine, y_fine, label=f"Spline between points {y[0]} и {y[2]}")
            plt.scatter(x, y, color='red', zorder=5)


            plt.title("Cubic spline and starting points")
            plt.xlabel("x")
            plt.ylabel("y")
            plt.legend()
            plt.grid(True)
            plt.show()
    except:
        pass
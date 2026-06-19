import numpy as np
from pymatgen.core import Structure


def read_bvse_grid(filename):
    with open(filename, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    a, b, c, alpha, beta, gamma = map(float, lines[1].split())
    nx, ny, nz = map(int, lines[2].split())

    values = np.array([float(x) for x in lines[3:]], dtype=float)
    if len(values) != nx * ny * nz:
        raise ValueError(f"Expected {nx*ny*nz} grid values, got {len(values)}")

    grid = values.reshape((nx, ny, nz))
    return (a, b, c, alpha, beta, gamma), grid


def find_local_minima(grid):
    is_min = np.ones_like(grid, dtype=bool)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if dx == 0 and dy == 0 and dz == 0:
                    continue
                shifted = np.roll(grid, shift=(dx, dy, dz), axis=(0, 1, 2))
                is_min &= (grid <= shifted)
    return is_min


def grid_indices_to_frac(indices, shape):
    nx, ny, nz = shape
    return np.array([[i / nx, j / ny, k / nz] for i, j, k in indices], dtype=float)


def min_image_distance(structure, frac1, frac2):
    d = np.array(frac1) - np.array(frac2)
    d -= np.round(d)
    cart = structure.lattice.get_cartesian_coords(d)
    return np.linalg.norm(cart)


def select_low_energy_sites_with_distance(structure, frac_coords, energies, n_li, li_li_min_dist):
    order = np.argsort(energies)
    frac_coords = frac_coords[order]
    energies = energies[order]

    selected_frac = []
    selected_energies = []

    for fc, en in zip(frac_coords, energies):
        ok = True
        for sfc in selected_frac:
            if min_image_distance(structure, fc, sfc) < li_li_min_dist:
                ok = False
                break
        if ok:
            selected_frac.append(fc)
            selected_energies.append(en)
            if len(selected_frac) == n_li:
                break

    if len(selected_frac) < n_li:
        raise RuntimeError(
            f"Could select only {len(selected_frac)} Li positions out of {n_li}. "
            f"Try decreasing Li_Li_min_dist or using a finer BVSE grid."
        )

    return np.array(selected_frac), np.array(selected_energies)


def choose_li_positions_from_bvse(
    structure,
    bvse_filename,
    n_li=16,
    Li_Li_min_dist=2.2,
    prefer_local_minima=True,
    el = 'Li'
):
    _, grid = read_bvse_grid(bvse_filename)

    if prefer_local_minima:
        minima_mask = find_local_minima(grid)
        candidate_indices = np.argwhere(minima_mask)

        # если минимумов мало, автоматически берем все точки сетки
        if len(candidate_indices) < n_li:
            print(
                f"Warning: only {len(candidate_indices)} local minima found, "
                f"falling back to all grid points."
            )
            candidate_indices = np.argwhere(np.ones_like(grid, dtype=bool))
    else:
        candidate_indices = np.argwhere(np.ones_like(grid, dtype=bool))

    candidate_frac = grid_indices_to_frac(candidate_indices, grid.shape)
    candidate_energies = grid[tuple(candidate_indices.T)]

    selected_frac, selected_energies = select_low_energy_sites_with_distance(
        structure=structure,
        frac_coords=candidate_frac,
        energies=candidate_energies,
        n_li=n_li,
        li_li_min_dist=Li_Li_min_dist,
    )

    new_structures_oneLi = []
    new_structure = structure.copy()
    for fc in selected_frac:
        new_structure.append(el, fc, coords_are_cartesian=False)
        new_structure_oneLi = structure.copy()
        new_structure_oneLi.append(el, fc, coords_are_cartesian=False)

        new_structures_oneLi.append( new_structure_oneLi )
    return selected_frac, selected_energies, new_structure, new_structures_oneLi




def grid_bvlain(file = None, st = None, mobile_ion = 'Li1+' ):
    from bvlain import Lain

    if file is None:
        file = st.write_cif()

    calc = Lain(verbose = False)
    atoms = calc.read_file(file, oxi_check = False)       # alternatively, you can use read_atoms() or read_structure()
    # print(atoms)

    oxi_map = {
        'Li': +1,
        'Na': +1,
        'C': +4,
        'F': -1,
        'La': +3,
        'Y': +3,
        'Ta': +5,
        'Cl': -1,
        'O': -2,
        'K': +1,
        'W': +6,
        'Al': +3,
        'Nb': +5,
    }

    atoms.arrays['oxi_states'] = np.array(
        [oxi_map[symbol] for symbol in atoms.get_chemical_symbols()],
        dtype=int
    )


    calc.atoms_copy = atoms


    params = {'mobile_ion':mobile_ion ,              # mobile specie
              'r_cut': 8.0,  #10                   # cutoff for interaction between the mobile species and framework
              'resolution': 0.4, #0.2                 # distance between the grid points
              'k': 30,               #100           # maximum number of neighbors to be collected for each point
              'use_softbv_covalent_radii': False # default is False, use True to compare results with softBV
    }
    _ = calc.bvse_distribution(**params)
    energies = calc.percolation_barriers(encut = 5.0)
    for key in energies.keys():
        print(f'{key[-2:]} percolation barrier is {round(energies[key], 4)} eV')


    calc.write_grd(file + '_bvse')

    return file + '_bvse.grd'





# =========================
# Example
# =========================
if __name__ == "__main__":
    structure = Structure.from_file("POSCAR")
    bvse_file = "bvse_grid.dat"

    selected_frac, selected_energies, new_structure = choose_li_positions_from_bvse(
        structure=structure,
        bvse_filename=bvse_file,
        n_li=16,
        Li_Li_min_dist=2.2,
        prefer_local_minima=True,
    )

    print("Selected Li sites:")
    for fc, e in zip(selected_frac, selected_energies):
        print(f"{fc}   E = {e:.6f}")

    new_structure.to(fmt="poscar", filename="POSCAR_with_16Li.vasp")
"""
Test to check work of fragment_placer for molecules

Author: Marina Titarenko

To do:
Addability to check reaction sites and place ions/other molecules according to the reaction sites found
"""
from siman.core.molecule import FragmentPlacer


import os
import numpy as np
from ase.io import read
from ase.data import vdw_radii
from scipy.spatial.distance import cdist


def test(output, molecule_file):
    #output - file with final geometry of complex
    #molecule_file - str, file with initial geometry of host molecule

    if os.path.isfile(output):
        print("CORRECT, Output file created")
    else:
        print("Output file missing")
        return

    atoms = read(output)
    molecule_atoms = read(molecule_file)
    host_atoms = len(molecule_atoms)

    host = atoms[:host_atoms]
    frag = atoms[host_atoms:]

    D = cdist(host.positions, frag.positions)

    r_host = np.array([vdw_radii[a.number] for a in host])
    r_frag = np.array([vdw_radii[a.number] for a in frag])

    cutoff = 0.85 * (r_host[:, None] + r_frag[None, :])

    if np.any(D < cutoff):
        print("Atom overlaps detected")
    else:
        print("NO OVERLAPPING")


if __name__ == "__main__":

    molecule = "acn.xyz"
    anion= "pf6.xyz"
    fragment_charge = -1
    output = "acn_pf6_complex.xyz"

    placer = FragmentPlacer()
    acn_comp = placer.place_fragment(host_file = molecule,
                                     fragment_file = anion,
                                     fragment_charge = -1,
                                     output = output)
    
    test(output, molecule)


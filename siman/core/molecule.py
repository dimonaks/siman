# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
import os

from pymatgen.core.structure import Molecule as Molecule_pymatgen
from siman import header
from siman.inout import read_structure
from siman.small_functions import makedir
from siman.header import printlog, runBash

import numpy as np
from ase.io import read
from ase.data import vdw_radii
from ase.optimize import BFGS
from scipy.spatial.distance import cdist
from scipy.spatial.transform import Rotation
from tblite.interface import Calculator
from tblite.ase import TBLite

class Molecule(Molecule_pymatgen):
    """Class for molecule structure representation based on pymatgen Molecule """
    def __init__(self, species=None, coords=None, *args, **kwargs):

        if species is None:
            species = []
        if coords is None:
            coords = []

        super().__init__(species, coords, *args, **kwargs)
    @classmethod
    def cast(cls, some_a: Molecule_pymatgen, filename = None):
        """Cast  Molecule_pymatgen into MyA."""
        assert isinstance(some_a, Molecule_pymatgen)
        some_a.__class__ = cls  # now mymethod() is available
        assert isinstance(some_a, Molecule)
        some_a.filename = filename
        return some_a

    # def __init__(self, filename = None, pymatgen_molecule_object = None):
    #     # super(Molecule, self).__init__(filename)

    #     if filename:
    #         self = read_structure(filename)
    #     elif pymatgen_molecule_object:
    #         self = pymatgen_molecule_object
    #         # print(pymatgen_molecule_object)
    #         # print(self)

    #     # self.len_units = 'Angstrom'

    def write_xyz(self, filename = None):

        if filename:
            filename += '.xyz'
        else:
            if hasattr(self, 'name') and self.name:
                name = self.name
            elif hasattr(self, 'filename') and self.filename:
                name = os.path.basename(self.filename)
            else:
                printlog('Warning! the molecule has no unique name, chemical formula will be used')
                name = self.formula.replace(' ', '')
            filename = 'xyz/'+name+'.xyz'
        makedir(filename)
        self.to(filename=filename, fmt = 'xyz')
        return filename

    def jmol(self, program = 'jmol'):
        ''

        filename = self.write_xyz()
        if 'jmol' in program :
            runBash(header.PATH2JMOL+' "'+filename+'"', detached = True)
        elif 'vesta' in program:
            runBash(header.PATH2VESTA+' "'+filename+'"', detached = True)
        return

    def get_volume(self):
        printlog('Molecule has no volume!', imp = 'n')
        return 0


class FragmentPlacer:
    """
    Places an ionic or molecular fragment near a host molecule

    The placement procedure consists of:
    1. Calculating of atomic charges by xTB.
    2. Selecting candidate binding sites based on atomic charges.
    3. Generating random fragment positions and orientations around the selected sites.
    4. Ranking generated poses using a scoring function.
    5. Performing an optional xTB geometry optimization of the best structure.

    Parameters:
    clash_scale: float, scaling factor applied to vdw radii when detecting steric overlaps
    n_sites: int, number of candidate binding sites considered
    n_trials: int, number of random poses generated for site
    distance_range: tuple(float, float), minimum and maximum placement distances (Å) 
                        between the host site and the fragment center
    """

    def __init__(self, clash_scale=0.9, n_sites=5, n_trials=1000, distance_range=(2.0, 5.0)):

        self.clash_scale = clash_scale
        self.n_sites = n_sites
        self.n_trials = n_trials
        self.distance_range = distance_range

    def place_fragment(self, host_file, fragment_file, fragment_charge, output=None, optimize=True):
        """
        Place a fragment near a host molecule

        Parameters:
        host_file: str, host structure file (all types, that can be read by ASE)
        fragment_file: str, fragment structure file (-//-)
        fragment_charge: int, charge of the fragment
        output: str, file name for output structure
        optimize (optional): If True, perform geometry optimization of final structure

        Returns:
        ase.Atoms - final complex geometry
        """

        host = read(host_file)
        fragment = read(fragment_file)

        host_q = self.get_xtb_charges(host)
        frag_q = self.get_xtb_charges(fragment)

        sites = self.find_sites(host_q, fragment_charge)

        best_score = np.inf
        best_pose = None

        for site in sites:

            for trial in range(self.n_trials):

                pose = self.generate_pose(host, fragment, site)

                score = self.score_pose(host,
                                        host_q,
                                        pose,
                                        frag_q)

                if score < best_score:
                    best_score = score
                    best_pose = pose

        structure = host + best_pose

        if optimize:
            structure = self.optimize(structure)

        structure.write(output)

        return structure

    def get_xtb_charges(self, atoms):

        calc = Calculator(method="GFN2-xTB",
                          numbers=atoms.numbers,
                          positions=atoms.positions)
        res = calc.singlepoint()

        return np.array(res["charges"])

    def find_sites(self, charges, fragment_charge):

        if fragment_charge < 0:
            idx = np.argsort(-charges)

        elif fragment_charge > 0:
            idx = np.argsort(charges)

        else:
            idx = np.argsort(np.abs(charges))[::-1]

        return idx[:self.n_sites]

    def generate_pose(self, host, fragment, site_idx):

        frag = fragment.copy()

        coords = frag.positions - frag.get_center_of_mass()
        coords = Rotation.random().apply(coords)

        direction = np.random.normal(size=3)
        direction /= np.linalg.norm(direction)

        distance = np.random.uniform(*self.distance_range)

        target = host.positions[site_idx] + distance * direction

        frag.positions = coords + target

        return frag

    def score_pose(self, host, host_q, frag, frag_q):

        D = cdist(host.positions, frag.positions)
        D[D < 0.5] = 0.5

        electrostatic = np.sum(host_q[:, None] * frag_q[None, :] / D)

        r_host = np.array([vdw_radii[a.number] for a in host])
        r_frag = np.array([vdw_radii[a.number] for a in frag])

        cutoff = self.clash_scale * (r_host[:, None] + r_frag[None, :])

        overlap = np.clip(cutoff - D, 0.0, None)
        penalty = np.sum(overlap**2)

        dmin = D.min()

        return electrostatic + 1000.0 * penalty + 0.1 * dmin

    def optimize(self, atoms, method="GFN2-xTB", fmax=0.05, steps=300):

        mol = atoms.copy()
        mol.calc = TBLite(method=method)

        opt = BFGS(mol)
        opt.run(fmax=fmax, steps=steps)

        return mol

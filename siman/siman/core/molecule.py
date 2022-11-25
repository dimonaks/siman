# Copyright (c) Siman Development Team.
# Distributed under the terms of the GNU License.
import os

from pymatgen.core.structure import Molecule as Molecule_pymatgen
from siman import header
from siman.inout import read_structure
from siman.small_functions import makedir
from siman.header import printlog, runBash, plt

class Molecule(Molecule_pymatgen):
    """Class for molecule structure representation based on pymatgen Molecule """
    
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
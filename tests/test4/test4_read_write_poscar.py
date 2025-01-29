"""
Test for read_poscar and write_poscar functions

Author: Marina Titarenko

To do:
Add ability read and write POSCAR file for MD simulation 
"""

from siman.inout import read_poscar
from siman.calc_manage import smart_structure_read
from siman.core.structure import Structure


def test(poscar):
	
	#poscar - list of POSCAR's with different structures

	for file in poscar:
		st = smart_structure_read(file)
		st.write_poscar('new_' + file)
		line_init = sum(1 for line in open(file))
		line_after = sum(1 for line in open('new_' + file +'.POSCAR'))
		
		if line_init == line_after:
			print('success')

		else:
			print('failure')


if __name__== "__main__":

	poscar = ['GEN.POSCAR', 'NVT.CONTCAR', 'NPT.CONTCAR']

	test(poscar)

#!/usr/bin/env python
'''
Script for creating migration path using IDPP methode for NEB calculations 
Taken from the Henkelman group's VTST Tools https://theory.cm.utexas.edu/vtsttools/index.html
'''

print(__name__)
if __name__ == "__main__":
	from aselite import read_vasp, write_vasp, NEB
	from sys import argv
	import os
	if '-h' in argv or 4 < len(argv) < 3:
		print('usage: nebmake.py POSCAR1 POSCAR2 num_images [-NOIDPP]')
		print
		exit(1)

	ini_atoms = read_vasp(argv[1])
	fin_atoms = read_vasp(argv[2])
	images = [ini_atoms]
	for i in range(int(argv[3])):
		images.append(ini_atoms.copy())
		images.append(fin_atoms)
		neb = NEB(images)
	if '-NOIDPP' in argv:
		neb.interpolate(mic=True)
	else:
		neb.interpolate('idpp',mic=True)
		dir_names = ['0'+str(i) if i < 10 else str(i) for i in range(len(images))]
		for i, image in zip(dir_names,neb.images):
			if not os.path.isdir(i):
				os.mkdir(i)
			write_vasp(i+'/POSCAR',image)
			print('Ok, all set up here.')
		print('For later analysis, put OUTCARs in folders 00 and ' + dir_names[-1])
else:
	from siman.external_tools.aselite import read_vasp, write_vasp, NEB
	from sys import argv
	import os
	
	def nebmaker(init, fin, n, idpp = True):
		
		ini_atoms = read_vasp(init)
		fin_atoms = read_vasp(fin)
		images = [ini_atoms]
		for i in range(n):
			images.append(ini_atoms.copy())
		images.append(fin_atoms)
		neb = NEB(images)
		

		if idpp == False :
			neb.interpolate(mic=True)
		else:
			neb.interpolate('idpp', mic=True)
			dir_names = ['0'+str(i) if i < 10 else str(i) for i in range(len(images))]
		for i, image in zip(dir_names,neb.images):
			if not os.path.isdir(i):
				os.mkdir(i)
			write_vasp(i+'/POSCAR',image)
			print('Ok, all set up here.')
		print('For later analysis, put OUTCARs in folders 00 and ' + dir_names[-1])


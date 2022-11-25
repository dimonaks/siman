import sys
import numpy as np


sys.path.append('/home/aksenov/Simulation_wrapper/siman') #provide your path to siman package
import header
from header import db # database dictionary
from calc_manage import smart_structure_read, add_loop, res_loop
from geo import supercell
from database      import read_database, write_database


read_database() # read saved results



"""1. Read structure and create supercells with defects"""
if 1:
	st = smart_structure_read(input_geo_file = 'Li2/POSCAR')
	
	sc = supercell(st, [10,10,10]) # create 3x3x3 supercell, sizes close to 10 10 10 A 
	





	sc_vac = sc.del_atom(25) # create vacancy defect in supercell

	# sc_vac.write_poscar('Li_vac/POSCAR') # write POSCAR to check geometry






	sc_oct = sc.add_atom([1/2 / 3, 1/2 / 3, 1 / 3], 'Li')  # create interstitial defect in octahedral position
	sc_tet = sc.add_atom([1/4 / 3, 1/2 / 3, 1 / 3], 'Li')  # create interstitial defect in tetrahedral position

	# sc_oct.write_poscar('Li_oct/POSCAR') 
	# sc_tet.write_poscar('Li_tet/POSCAR')  





	sc_vacum = sc.add_vacuum(0, 5)  # add 5 A of vacuum along first lattice vector - allows to create 100 surface
	# sc_vacum.write_poscar('Li_vacuum/POSCAR')






"""2. Run calculations of pure cell and with point defects"""

if 0:
	it_folder = 'recitationCD/'
	add_loop('Li333',    'opt', 1, up = 'up2', input_st = sc,     it_folder = it_folder, run = 1) 
	add_loop('Li333vac', 'opt', 1, up = 'up2', input_st = sc_vac, it_folder = it_folder, run = 1) 
	add_loop('Li333tet', 'opt', 1, up = 'up2', input_st = sc_tet, it_folder = it_folder, run = 1) 
	add_loop('Li333oct', 'opt', 1, up = 'up2', input_st = sc_oct, it_folder = it_folder, run = 1) 

"""3. Read results"""
if 0:
	res_loop('Li333',    'opt', 1, show = 'fo')
	res_loop('Li333vac', 'opt', 1, show = 'fo')
	res_loop('Li333tet', 'opt', 1, show = 'fo')
	res_loop('Li333oct', 'opt', 1, show = 'fo')









"""4. Create 100 surfaces with different vacuum thickness, run and read"""
thickness_list = [1,2,3,4,5,6,7,8,9,10]
if 0:
	

	sc = supercell(st, [10,3,3]) # create 3x1x1 supercell, sizes close to 10 3 3 A; we need slab!

	# add_loop('Li311',    'opt', 1, up = 'up2', input_st = sc, it_folder = 'recitationCD/',     run = 1) 

	res_loop('Li311',    'opt', 1)

	for thickness in thickness_list:
		''
		# add_loop('Li_suf_'+str(thickness), 'opt', 1, up = 'up2', input_st = sc.add_vacuum(0, thickness), it_folder = 'recitationCD/', run = 0)
		# res_loop('Li_suf_'+str(thickness), 'opt', 1, up = 'up2')









"""5. Calculate defect energies"""
if 0:
	Li311 = db['Li311', 'opt', 1]
	Li333 = db['Li333', 'opt', 1]
	Livac = db['Li333vac', 'opt', 1]
	Litet = db['Li333tet', 'opt', 1]
	Lioct = db['Li333oct', 'opt', 1]


	#1 Point defects formation energy
	natom = Li333.end.natom
	print('Vacancy formation energy = {:3.3f} eV, E(relax) = {:3.3f} eV'.format(  
		Livac.e0 - Li333.e0 * (natom-1)/natom,    Livac.list_e_sigma0[0] - Livac.e0   )    )
	

	print('Tet int formation energy = {:3.3f} eV, E(relax) = {:3.3f} eV'.format(  
		Litet.e0 - Li333.e0 * (natom+1)/natom,   Litet.list_e_sigma0[0] - Litet.e0   )    )
	

	print('Oct int formation energy = {:3.3f} eV, E(relax) = {:3.3f} eV'.format(  Lioct.e0 - Li333.e0 * (natom+1)/natom,
	    Lioct.list_e_sigma0[0] - Lioct.e0   )    )
	
	
	# Evac = 0.34 eV  Experimental value Feder (1970) from [Jacucci1979]
	# Evac = 0.45 eV  Pair potentials from [Jacucci1979]




if 0:
	#2 Surface energy, depending on thickness of vacuum layer

	A = np.linalg.norm( np.cross(Li311.end.rprimd[1] , Li311.end.rprimd[2]) ) # surface area

	gamma_list = []
	for thickness in thickness_list:
		Li_suf = db['Li_suf_'+str(thickness), 'opt', 1]

		gamma = (Li_suf.e0 - Li311.e0)/A * header.eV_A_to_J_m / 2
		gamma_list.append(gamma)
		print('Surface energy = {:3.3f} J/m2; vacuum thickness = {:2d} A'.format(gamma, thickness) )


	#gamma(100) = 0.46 J/m2    PAW PBE from [Jackle2014] 



	if 1: #plot gamma - thickness
		from picture_functions import fit_and_plot
		fit_and_plot(a = (thickness_list, gamma_list, '-o'), filename = 'gamma', ylabel = 'Surface energy, J/m$^2$', xlabel = 'Vacuum thickness, $\AA$')





"""6. Optimize 100 surface"""
if 0:
	sc = supercell(st, [10,10,10]) # create 3x1x1 supercell, sizes close to 10 3 3 A; we need slab!
	thickness = 5

	from set_functions import read_vasp_sets
	read_vasp_sets([('opt_sym0', 'opt',{'ISYM':0})])

	# add_loop('Li333_suf_'+str(thickness), 'opt_sym0', 1, up = 'up2', input_st = sc.add_vacuum(0, thickness), it_folder = 'recitationCD/', run = 1)
	res_loop('Li333_suf_'+str(thickness), 'opt_sym0', 1, up = 'up2')


write_database() # save calculation results in compact format
#!/usr/bin/python
from glob import glob
import os

from shutil import copyfile, rmtree
#print(os.listdir('.'))
#print(glob(os.getcwd()))
for d in os.listdir('.'):
    if os.path.isdir(d):
	#print(d)
	if os.path.exists(d+'/'+'CONTCAR.static') and not os.path.exists(d+'/'+'CONTCAR.relax'):
	    print('no static', d)
            #copyfile(d+'/str_relax.out', d+'/str_hint.out')
	    #open(d+'/wait', 'a').close()
	if not os.path.exists(d+'/CONTCAR.static'):
	   ""
            #rmtree(d)
  	if not os.path.exists(d+'/force.out'):
           print(d, 'is not finished')
	   continue
        if os.path.exists(d+'/vasp.out'):
	   print(d, 'is running')
	   continue
	
        with open(d+'/force.out') as f:
	   maxforce = 0
	   forces = []
           for line in f:
	      force = [abs(float(f)) for f in line.split()]
	      forces.extend(force)
        if not forces:
           print(d, 'force.out is empty')
 	   continue
	maxf = max(forces)*1000
        print('maxforce in {0:3s} is {1:6.1f} meV/A'.format(d, maxf))
	if 1:
            os.chdir(d)
            os.system('gunzip OUTCAR.static.gz.gz')
	    if os.path.exists('energy'):
	        os.remove('energy')
            os.system('extract_vasp')
	    os.chdir('..')
        if maxf > 50:
	   ''
           copyfile(d+'/str_relax.out', d+'/str_hint.out')
        
	   open(d+'/wait', 'a').close()

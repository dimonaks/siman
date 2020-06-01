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
        
        if os.path.exists(d+'/wait'):
            print(d, 'waiting for running')
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

        if 0:
            os.chdir(d)
            if os.path.exists('CONTCAR.static.gz.gz'):
                os.system('gunzip OUTCAR.static.gz.gz')
            if os.path.exists('energy'):
                os.remove('energy')
            os.system('extract_vasp')
            os.chdir('..')
        if d == '1':
            continue #skip 1 as energy is usually provided by hands due to DFT+U problem
	
        if maxf > 700:
            ''
            #get current step
            os.chdir(d)

            step_file_list = glob('step.*')
            if len(step_file_list) > 0:
                step_file = step_file_list[-1]
                step = step_file.split('.')[1]
                os.rename(step_file, 'step.'+str(int(step)+1))
            else:
                step = '1'
                open('step.1', 'a').close()


            if os.path.exists('OUTCAR.static.gz'):
                os.rename('OUTCAR.static.gz', 'OUTCAR-'+step+'.static.gz')
            if os.path.exists('OUTCAR.relax.gz'):
                os.rename('OUTCAR.relax.gz', 'OUTCAR-'+step+'.relax.gz')

            copyfile('str_relax.out', 'str_hint.out')
            open('wait', 'a').close()
            os.chdir('..')

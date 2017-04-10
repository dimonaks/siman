"""
fit_tool.py by Aksyonov Dmitry, Skoltech, Moscow

Reads energies and volumes from provided VASP OUTCAR's files (CONTCAR's should exist and have the same naming pattern), 
fits them using EOS and 
creates file '100.POSCAR' with minimum energy volume.
"""
import sys
# print sys.argv
sys.path.append('/usr/lib64/python2.7/site-packages/numpy')
import numpy as np

class Simple(object):
    def __init__(self):
        self.xcart = []
        self.rprimd = []
        self.xred  = []
        self.nznucl = []
        self.els    = [] #list of elements


sts = []
for filename in sys.argv[1:]:
    st = Simple()
    st.filename = filename
    # print filename
    with  open(filename, 'r') as outcar:
        outcarlines = outcar.readlines()

    for i, line in enumerate(outcarlines):

        if "energy  without entropy=" in line:
            st.energy_sigma0 = float(line.split()[6]) #energy(sigma->0)
        
        if "volume of cell" in line:
            st.vol = float(line.split()[4])


        # if 'POSITION' in line:
        #     xcart = []
        #     j = 0
        #     local_line = outcarlines[i+2]

        #     while '-------' not in local_line:
        #         xcart.append( [ float(x) for x in local_line.split()[0:3]   ]  )
        #         j+=1
        #         local_line = outcarlines[i+j+2]
        #     st.xcart = xcart

    with  open(filename.replace('OUTCAR', 'CONTCAR'), 'r') as f:
        st.name = f.readline().strip()
        mul = float( f.readline() )

        for i in 0, 1, 2:
            vec = f.readline().split()
            st.rprimd.append( np.array([float(vec[0])*mul, float(vec[1])*mul, float(vec[2])*mul])  ) 


        ilist = f.readline().split() #nznucl of elements?
        try:
            int(ilist[0])
            vasp5 = False
        except:
            vasp5 = True

        if vasp5:
            for el in ilist:
                st.els.append(el)
            ilist = f.readline().split() #nznucl of elements?
        
        for z in ilist:
            st.nznucl.append( int(z)  )

        type_of_coordinates = f.readline()

        if "dir" in type_of_coordinates or 'Dir' in type_of_coordinates:
            for nz in st.nznucl:

                for i in range(nz):
                    vec = f.readline().split()
                    st.xred.append( np.array( [float(vec[0]), float(vec[1]), float(vec[2])] ) )


    sts.append(st)

	# print xcart
energies = [st.energy_sigma0  for st in sts]
volumes  = [st.vol            for st in sts]
power    = 3

coeffs   = np.polyfit(volumes, energies, power)        

fit_func = np.poly1d(coeffs)
fine_volumes  = np.linspace(min(volumes), max(volumes))
fine_energie  = fit_func(fine_volumes)

root_vols = fit_func.deriv().r
root_engs = fit_func(root_vols)

i_min = np.argmin(root_engs)
e_min = np.real( root_engs[i_min] )
v_min = np.real( root_vols[i_min] )

# print v_min                              
# print np.real(e_min)
# import matplotlib.pyplot as plt
# plt.plot(fine_volumes, fine_energie)
# plt.show()

#find st with vol most close to v_min
i_tar = 0
d1 = abs(sts[0].vol-v_min)
for i, st in enumerate(sts):
    d2 = abs(st.vol - v_min) 
    if d2 < d1:
        d1 = d2
        i_tar = i

st = sts[i_tar]
# st = sts[7]
#use i_tar
scale = (v_min/st.vol)**(1./3)
# print scale
# print st.vol
for i in (0,1,2):
    st.rprimd[i]*=scale 

print ('diff between new volume and target volume', np.dot( st.rprimd[0], np.cross(st.rprimd[1], st.rprimd[2])  ) - v_min)

with  open('100.POSCAR', 'w') as f:
    f.write(st.name+' '+st.filename+' :fit_tool.py\n')

    f.write("{0:18.15f}\n".format(1.0))
    

    for i in 0, 1, 2:
        f.write('{0:10.6f} {1:10.6f} {2:10.6f}\n'.format(st.rprimd[i][0], st.rprimd[i][1], st.rprimd[i][2] ) )

    for n in st.nznucl:    
        f.write(str(n)+' ')

    f.write('\n')

    f.write('Direct\n')

    for x in st.xred:
        f.write("  {0:18.16f}  {1:18.16f}  {2:18.16f}\n".format(x[0], x[1], x[2]) )
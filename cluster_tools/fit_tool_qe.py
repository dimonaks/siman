"""
tmp version for qe analysis 

fit_tool_qe.py by AG based on fit_tool.py  

Reads energies and volumes from provided QE output files fits them using EOS and creates file '100.POSCAR' (part of the QE input file with atomic positions and unitcell vectors) with minimum energy volume.

"""
import sys
import numpy as np
import re


class Simple(object):
    def __init__(self):
        self.xcart = []
        self.rprimd = []
        self.xred = []
        self.nznucl = []
        self.els = []  # list of elements
        self.atoms = {}


pattern_lat = r'a\(\d+\) = \(\s*([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s*\)'
pattern_celldm = r'celldm\(1\)=(\s*[-+]?\d+\.\d+)'
pattern_atoms = r'(\d+)\s+(\w+)\s+tau\(\s*(\d+)\)\s+=\s+\(\s*([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+([-+]?\d+\.\d+)\s+\)'
sts = []
for filename in sys.argv[1:]:
    st = Simple()
    st.filename = filename
    lat = []
    with open(filename, 'r') as outcar:
        outcarlines = outcar.readlines()
        for i, line in enumerate(outcarlines):
            if "! " in line:
                st.energy_sigma0 = float(line.split()[4])  # energy(sigma->0)
                # print(st.energy_sigma0)
            if "unit-cell volume" in line:
                st.vol = float(line.split()[3])

            if re.findall(pattern_lat, line):
                matches = re.findall(pattern_lat, line)
                for match in matches:
                    a1, a2, a3 = match
                    lat.append([float(a1), float(a2), float(a3)])

            if re.findall(pattern_celldm, line):
                celldm_1_value = re.findall(pattern_celldm, line)
            if 'tau' in line:
                tmp = line.split()
                site_num, atom,  x, y, z = tmp[0], tmp[1], tmp[6], tmp[7], tmp[8]
                st.atoms[site_num] = {'element': atom, 'pos': np.array(
                    [float(x), float(y), float(z)])}
        for i in range(3):
            st.rprimd.append(
                np.array([a*float(celldm_1_value[0]) for a in lat[i]]))

            for site in st.atoms:
                st.xred.append(st.atoms[site]['pos'])
        sts.append(st)

energies = [st.energy_sigma0 for st in sts]
volumes = [st.vol for st in sts]


vectors1 = [np.linalg.norm(st.rprimd[0]) for st in sts]
vectors2 = [np.linalg.norm(st.rprimd[1]) for st in sts]
power = 3
coeffs = np.polyfit(volumes, energies, power)

fit_func = np.poly1d(coeffs)
fine_volumes = np.linspace(min(volumes), max(volumes))
fine_energie = fit_func(fine_volumes)

root_vols = fit_func.deriv().r
root_engs = fit_func(root_vols)

i_min = np.argmin(root_engs)
e_min = np.real(root_engs[i_min])
v_min = np.real(root_vols[i_min])


i_tar = 0
d1 = abs(sts[0].vol-v_min)
for i, st in enumerate(sts):
    d2 = abs(st.vol - v_min)
    if d2 < d1:
        d1 = d2
        i_tar = i

st = sts[i_tar]

if np.std(vectors1)+np.std(vectors2) < 1e-5:  # Vectors 1 and 2 do not change
    scale = (v_min/st.vol)

    i_scale_list = [2]  # only third vector is scaled - c_scale regime
else:
    scale = (v_min/st.vol)**(1./3)

    i_scale_list = [0, 1, 2]  # uniform_scale regime
for i in i_scale_list:
    st.rprimd[i] *= scale

print('diff between new volume and target volume', np.dot(
    st.rprimd[0], np.cross(st.rprimd[1], st.rprimd[2])) - v_min)
# print(st.rprimd)
bohr2angst = 0.529177249
with open('100.POSCAR', 'w') as f:
    f.write('CELL_PARAMETERS angstrom\n')
    for i in 0, 1, 2:
        f.write('{0:10.6f} {1:10.6f} {2:10.6f}\n'.format(
            st.rprimd[i][0]*bohr2angst, st.rprimd[i][1]*bohr2angst, st.rprimd[i][2]*bohr2angst))

    f.write('\n')

    f.write('ATOMIC_POSITIONS crystal\n')

    for x in st.atoms:
        f.write("{}  {}  {}  {}\n".format(
            st.atoms[x]['element'], st.atoms[x]['pos'][0], st.atoms[x]['pos'][1], st.atoms[x]['pos'][2]))

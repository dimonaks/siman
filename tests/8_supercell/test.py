"""
Test supercell

Author: Dmitry Aksyonov
"""

from siman.inout import smart_structure_read
from siman.geo import supercell
import os

files = os.listdir('geo')
# print(files)

try:
    for file in files: 
        st = smart_structure_read('geo/'+file)
        st = st.get_primitive_cell()
        st = supercell(st, [10,10,10], test_natom = 0)
        print(st.get_volume())
        st.write_poscar()
    exit('success')
except:
    exit('failure')

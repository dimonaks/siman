from siman.calc_manage import smart_structure_read
from siman.geo import replace_x_based_on_symmetry
#CuIn2Se3I
st = smart_structure_read('/home/d.aksenov/scientific_projects/boldyreva/CuIn2Se3I_In4Se4.vasp') #считать POSCAR, он получен в VESTA из cif файла, в котором установлены целочисленные заселенности только In и Se

st = st.replic([2,1,1])
# st.printme()
replace_x_based_on_symmetry(st, 'In', 'Cu', x = 0.25, info_mode = 1, silent = 0) # узнать список возможных групп симметрии при замене 25% In на Cu

sts, atrs = replace_x_based_on_symmetry(st, 'In', 'Cu', x = 0.25, sg = 215, silent = 0, mode = 'rep') # сделать замену для группы 215 и вернуть полный список структур

sts[0].write_poscar() # записать POSCAR для первой структуры с группой симметрии 215
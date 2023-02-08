"""

Test for SIMAN plot_dos function

AUTHOR: 

A. Boev

TODO:

- Add ability to check plotted graphs

"""
from siman.calculators.vasp import CalculationVasp
from siman.core.structure import Structure
from siman.dos_functions import plot_dos
from siman import header

header.warnings = ''

cl = CalculationVasp(output = 'LCO_Mn/1.OUTCAR')
cl.read_results(show='')
# cl.end.printme()

def printlog(test_name,local_test_name, status):
    print(f'{test_name:30s} -> {local_test_name:30s} - {status:10}')
    


def test_method_dos(cl, test_name = 'test_method_dos'):
    print('Starting ...')
    
    #1
    local_test_name = 'i_at test'
    try:
        cl.dos(dostype = 'total', i_at = 0,res_show = '')
        cl.dos(dostype = 'total', i_at = 127,res_show = '')
        status = 'ok!'
    except Exception as exception:
        status = f'Error! {type(exception).__name__}\n{exception}'
    printlog(test_name,local_test_name, status)
    
    test_name = ''

    #2
    local_test_name = 'iatoms = [i] test'
    try:
        cl.dos(dostype = 'total', iatoms = [127],res_show = '')
        status = 'ok!'
    except Exception as exception:
        status = f'Error! {type(exception).__name__}\n{exception}'
    printlog(test_name,local_test_name, status)

    #3
    local_test_name = 'iatoms = [i,j] test'
    try:
        cl.dos(dostype = 'total', iatoms = [126,127],res_show = '')
        status = 'ok!'
    except Exception as exception:
        status = f'Error! {type(exception).__name__}\n{exception}'
    printlog(test_name,local_test_name, status)


test_method_dos(cl)

def test_function_plot_dos(cl,test_name = 'test_function_plot_dos'):
    print('Starting ...')
    
    # 1
    local_test_name = 'dostype total'
    try:
        plot_dos(cl, dostype = 'partial', iatom = 128, orbitals = ('d','p6'),show=0)
        status = 'ok!'
    except Exception as exception:
        status = f'Error! {type(exception).__name__}\n{exception}'
    printlog(test_name,local_test_name, status)
    test_name = ''
    
    # 2
    local_test_name = 'dostype partial'
    try:
        plot_dos(cl, dostype = 'partial', iatom = 128, orbitals = ('d','p6'),show=0)
        status = 'ok!'
    except Exception as exception:
        status = f'Error! {type(exception).__name__}\n{exception}'
    printlog(test_name,local_test_name, status)

    # 3
    local_test_name = 'dostype diff_total'
    try:
        plot_dos(cl, cl2 = cl, dostype = 'diff_total', iatom = 128, iatom2 = 128, orbitals = ('d','p6'),show=0)
        status = 'ok!'
    except Exception as exception:
        status = f'Error! {type(exception).__name__}\n{exception}'
    printlog(test_name,local_test_name, status)


    # 4
    local_test_name = 'cl1,cl2 objects'
    try:
        plot_dos(cl, cl2 = cl, dostype = 'partial', iatom = 128, iatom2 = 1, orbitals = ('d','p6'),show=0)
        status = 'ok!'
    except Exception as exception:
        status = f'Error! {type(exception).__name__}\n{exception}'
    printlog(test_name,local_test_name, status)
    

    # 5
    orb_list = ['s', 'p', 'd', 'py', 'pz', 'px', 'dxy', 'dyz', 'dz2', 'dxz', 'dx2', 
                'p6', 'd6', 'p_all', 'd_all']
    for orb in orb_list:
        local_test_name = f'orbitals test ({orb})'
        try:
            if orb in ['s', 'p', 'd', 'py', 'pz', 'px', 'p_all']: iatom = 127
            else: iatom = 128
            plot_dos(cl,  dostype = 'partial', iatom = iatom, orbitals = [orb], show=0, plot_param = {'figsize':(10,12)})
            status = 'ok!'
        except Exception as exception:
            status = f'Error! {type(exception).__name__}\n{exception}'
        printlog(test_name,local_test_name, status)


# plot_dos(cl1, cl2 = None, dostype = None, iatom = None, iatom2= None,
#     orbitals = ('s'), up = None, neighbors = 6, show = 1, labels = None,
#     path = 'dos', xlim = (None, None), ylim = (None,None), savefile = True, plot_param = {}, suf2 = '', nsmooth = 3,
#     lts2 = '--', split_type = 'octa', plot_spin_pol = 1, show_gravity = None, 
#     efermi_origin = True, efermi_shift = 0, invert_spins  = 0, name_suffix = '', image_name = None, color_dict = None):
test_function_plot_dos(cl)
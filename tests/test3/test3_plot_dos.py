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


# printlog might be general function for all the developed tests
def printlog(test_name,local_test_name, status):
    print(f'{test_name:30s} -> {local_test_name:30s} - {status:10}')
    


def test1(cl):
    test_name = 'test_method_dos'
    # testing method self.dos() of Calculation() object
    
    print(f'\n\nStarting {test_name}...')
    
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



def test2(cl):
    test_name = 'test_function_plot_dos'
    # testing the function plot_dos() from siman.dos_functions
    print(f'\n\nStarting {test_name}...')
    
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




if __name__== "__main__":
    
    cl = CalculationVasp(output = 'LCO_Mn/1.OUTCAR')
    cl.read_results(show='')

    test1(cl)
    test2(cl)
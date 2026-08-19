"""
Test for vasp OUTCAR parser, 
1. NPT calculations, read rprimd

Author: Dmitry Aksyonov
"""

from siman.inout import read_poscar
from siman.calc_manage import smart_structure_read
from siman.core.structure import Structure
from siman.calculators.vasp import CalculationVasp
import numpy as np 


def test(outcar):
    
    try:
        cl  = CalculationVasp(output = outcar, )
        cl.res(show = 'md_geo.npt')
        np.set_printoptions(formatter={'float': '{:6.2f}'.format})
        for s, r in zip (cl.list_stress, cl.end.list_rprimd):
            print(np.array(r))
            print('Stress, GPa:', np.array(s)/1000)
        result = 'success'
    except:
        result = 'failure'

    print(result)

if __name__== "__main__":


    test('/ssd2/scientific_projects/siman_tests_datafiles/5/NPT_OUTCAR')

    # return result
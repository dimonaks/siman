"""
Test MAE

Author: Dmitry Aksyonov
"""
from siman.calc_manage import add, res, smart_structure_read
from siman.set_functions import read_vasp_sets
from siman.header import db

def single_volume():

    #Define new sets

    read_vasp_sets([('mag', 'static', 
        {'ISPIN':2, 'LORBIT':11, 'magnetic_moments':{'Fe':5},
        'EDIFF':'1E-7', 'LREAL': '.FALSE.', 'LASPH':'.TRUE.', 'PREC':'Accurate',
        'ISMEAR':0, 'SIGMA':0.01, 'KSPACING':None,  }, 'override')])

    noncoliniar_pack = {'LNONCOLLINEAR':True,'LSORBIT':'.TRUE.', 'LASPH':'.TRUE.', 
    'LORBMOM':'.TRUE.', 'ICHARG':11, 'EDIFF':'1E-8', 'ISYM':-1, 'KGAMMA':False,
    'magnetic_moments':{'Fe':2.186}} # non-selfconsistent

    read_vasp_sets([('ncl', 'mag', noncoliniar_pack, 'override')])
    
    st = smart_structure_read('/ssd2/ydisk/15_Simulation_wrapper/siman/tests/t9_mae/POSCAR_Fe1bcc')
    ngkpt = [20,20,20]
    mk = 2 # dimensionality of shift mesh, e.g. 2x2x2

    run = 0

    if 0:
        #step 1. Calculate magnetic charge density

        if run:
            add('Fe1bcc', 'mag', 1, input_st = st, it_folder = 'Fe/MAE',
            params = {'ngkpt':ngkpt},
            )
        else:
            cl_mag = db['Fe1bcc', 'mag', 1]
            cl_mag.res()
            cl_mag.get_file('1.CHGCAR')


    else:
        #step 2. Calculate MAE using non-selfconsitent calculations
        # be sure that you have copied siman folder to ~/tools; 
        # check that you have ~/tools/siman/cluster_runners/shiftk.py file on cluster
        cl_mag = db['Fe1bcc', 'mag', 1]


        energies = []
        for d in '0 0 1', '1 1 1':
            dns = d.replace(" ", "")
            if run:
                add('Fe1bcc_'+dns, 'ncl', 1, input_st = st, it_folder = 'Fe/MAE',
                    calc_method = 'shiftk_average', 
                    params = {'shiftk_average':{'mk':mk},'ngkpt':ngkpt, 'chgcar':cl_mag.path['charge']},
                    update_set = {'SAXIS':d}
                    )
            else:
                cl_ncl = db['Fe1bcc_'+dns, 'ncl', 1]
                cl_ncl.res(up = 'up2')
                print(cl_ncl.shiftk_energies)
                average_energy = sum(cl_ncl.shiftk_energies) / len(cl_ncl.shiftk_energies)
                energies.append(average_energy)
        if not run:
            print('MAE is {:.1f} ueV'.format( (energies[0]-energies[1])*1e6  ) )

single_volume()

if __name__== "__main__":


    try:
        single_volume()
        exit('9_mae.single_volume: success')
    except:
        exit('9_mae.single_volume: failure')



    try:
        volume_scan()
        exit('9_mae.volume_scan: success')
    except:
        exit('9_mae.volume_scan: failure')
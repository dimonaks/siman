# -*- coding: utf-8 -*- 
#Copyright Aksyonov D.A
from __future__ import division, unicode_literals, absolute_import 
import sys, copy, re, os, shutil
#external
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from siman import header
from siman.header import printlog, calc, db
from siman.picture_functions import fit_and_plot
from siman.small_functions import merge_dics as md, makedir
from siman.calc_manage import add_loop, name_mod_supercell, res_loop, inherit_icalc, push_figure_to_archive, smart_structure_read
from siman.neb import add_neb
from siman.classes import Calculation
from siman.analysis import calc_redox,  matrix_diff
from siman.geo import create_deintercalated_structure, remove_one_atom, remove_half_based_on_symmetry, remove_half, create_replaced_structure, create_antisite_defect3, determine_symmetry_positions
from siman.inout import write_occmatrix
from siman.database import add_to_archive_database

mpl.rcParams.update({'font.size': 22})



def workflow():

    # manually_remove_from_struct_des(struct_des, 'Na111.su')


    #Na
    # add_loop('Na111', '8u', 1, up = 'up1', input_geo_format = 'vasp', it_folder = 'Na/bcc')

    # Li
    # res_loop('Li111', '8Utmb4', 1, up = 'up1',)
    # res_loop('Li111', '2b4', 1, up = 'up1', input_geo_format = 'mat_proj', it_folder = 'Li/bcc')
    # res_loop('Li111', '2b4', 1, up = 'up1')
    # add_loop('K111', '8Utmb4', 1, up = 'up1', input_geo_format = 'mat_proj', it_folder = 'K/bcc')

    # add_loop('Li111', '', 1, ise_new = '1ub4', savefile = 'o', 
    #     calc_method = 'uniform_scale', scale_region = [-4,4], it_folder = 'Li/bcc/scaled/', cluster = 'skol')

    # add_loop('Na111', '', 1, ise_new = '1ub4', savefile = 'o', 
    #     calc_method = 'uniform_scale', scale_region = [-4,4], it_folder = 'Na/bcc/scaled/', cluster = 'skol')

    # add_loop('K111',  '', 1, ise_new = '1ub4', savefile = 'o', 
    #     calc_method = 'uniform_scale', scale_region = [-4,4], it_folder = 'K/bcc/scaled/', cluster = 'skol')

    # res_loop('Li111.su', ['1ub4'], [1, 2, 3, 4, 5, 6, 7, 100], up = 'up2', show = 'fofit', analys_type = 'fit_a'  )     # , on 2016-10-19  
    # res_loop('Na111.su', ['1ub4'], [1, 2, 3, 4, 5, 6, 7, 100], up = 'up2', show = 'fofit', analys_type = 'fit_a'  )     # , on 2016-10-19  
    # res_loop('K111.su',  ['1ub4'], [1, 2, 3, 4, 5, 6, 7, 100], up = 'up2', show = 'fofit' , analys_type = 'fit_a' ) 


    # print (calc)




    """OCC"""
    # add_loop('LiFePO4', '1uOCC', 1, input_geo_format = 'vasp', show = 'fo')



    """Convergence"""
    # add_loop('LiFePO4', '1u', 1, show = 'force')
    # res_loop('LiCoO2', '9uk1', 1, show = 'force')
    # res_loop('LiCoO2', '8ue4', 1, show = 'force')
    # res_loop('LiCoO2', '8ue5', 1, show = 'force')



    """Uniform scale"""
    #pmna from LiFePO4
    # inherit_icalc('replace_atoms', 'NaFePO4.pnma', 1, ('LiFePO4','9u',1), calc, atom_new = 'Na', atom_to_replace = 'Li', it_folder = 'NaFePO4')    #on 2014-02-26
    # inherit_icalc('replace_atoms', 'KFePO4.pnma', 1, ('LiFePO4','9u',1), calc, atom_new = 'K', atom_to_replace = 'Li', it_folder = 'KFePO4')    #on 2014-02-26

    # res_loop('NaFePO4.pnma','1u', 1)
    # res_loop('KFePO4.pnma','1u', 1, show = 'fo')



    # for ise in ['1u', 1ui', '4ui']:
    # for ise in ['4uiWF', '4uiCF', '4uiCWF',]:
    for ise in ['4ui', '4uiWF', ]:
        ''
        # add_loop('KFePO4.pnma','1u', 1, ise_new = '4uiCWF', savefile = 'o', calc_method = 'uniform_scale', scale_region = (4,12), inherit_option = 'inherit_xred' )  
        # add_loop('FePO4.pnma2','1u', 1, ise_new = '4ui', savefile = 'o', calc_method = 'uniform_scale', scale_region = (-4,4), inherit_option = 'inherit_xred'  )  

        # res_loop('KFePO4.pnma.su',ise, [1,2,3,4,5,6,100], up = '1', show = 'fo', analys_type = 'fit_a')  

    #NaFePo
    # add_loop('NaFePO4.pnma','1u', 1, ise_new = '4ui', savefile = 'o', calc_method = 'uniform_scale', scale_region = (-2,6), inherit_option = 'inherit_xred' )  
    # res_loop('NaFePO4.pnma.su','4ui', list(range(1,8) )+[100], analys_type = 'fit_a' )  


    # sys.exit()



    inherit_option = None
    calc_method = 'neb'


    """LiFePO4  """

    """IS  """
    #111
    # add_neb(calc['LiFePO4', '9u', 1], ise_new = '1uD', images = 3, corenum = 15, inherit_option = None, mag_config = None, replicate = None)
    # res_loop('LiFePO4v1.n3Li1v1', '1',     [1],  show = 'mep', old_behaviour = 0 )    


    # for ise in ['1ue6', '1uk1', '1ut' ]:
    for ise in     ['1ue5k15']:
    # for ise in ['1uMM40','1up', '1upA' , '1upAMMD', ]:
    # for ise in ['1up', '1upA']:
    # for ise in ['1uMM40', '1up']:
    # for ise in ['1uMM40', '1uA']:
    # for ise in ['1uMM40','1uB03' , '1uB01' , '1uB005' ,]:
    # for ise in ['1m', '1uA', '1uLF', '1ulm', '1uB05', '1uWC1000', '1ulmA04', '1uMM40']:
    # for ise in [ '1m' ]: #inherit_magmom = False
    # for ise in [ '1uD', '1uMM40' ]: 
        ''
        # add_neb(calc['LiFePO4', '9u', 1], ise_new = ise, images = 3, corenum = 15, inherit_option = None, mag_config = None, inherit_magmom = True)
        # res_loop('LiFePO4v1.n3Li1v1mp', ise,     [5],  show = 'fo', old_behaviour = 0 )    
        # res_loop('LiFePO4v1.n3Li1v1ms', ise,     [1,2,3,4,5],  show = 'me', old_behaviour = 0 )    


    # res_loop('LiFePO4.nebLi1v1','1u',           [5,],  show = 'en', old_behaviour = 0 )    
    # res_loop('LiFePO4v1.n3Li1v1mp','1uMMD',     [5,],  show = 'en', old_behaviour = 0 )    
    # res_loop('LiFePO4.nebLi1v1','1ulong', [1],  show = 'occ', old_behaviour = 0 )   


    if 0: #checking initial magmom - in this case ms and mp are equivalent, 
        ''
        #but still some differences in final for 1uA. the situat. is much better for 1upAMMD, where smaller EDIFF is used
        #for 1m no such problem
        #intersting moment - the elapsed time could different for the same run
        #since the base calc['LiFePO4', '9u', 1] did not have magmom
        # print calc['LiFePO4v1.n3Li1v1mp', '1uA', 1].init.magmom
        # print calc['LiFePO4v1.n3Li1v1ms', '1uA', 1].init.magmom
        # print calc['LiFePO4', '9u', 1].end.magmom
        # res_loop('LiFePO4v1.n3Li1v1mp', '1uA',     [1],  show = '', old_behaviour = 0 )    
        # res_loop('LiFePO4v1.n3Li1v1ms', '1uA',     [1],  show = '', old_behaviour = 0 )    

        # res_loop('LiFePO4v1.n3Li1v1mp', '1upAMMD',     [4],  show = '', old_behaviour = 0 )    
        # res_loop('LiFePO4v1.n3Li1v1ms', '1upAMMD',     [4],  show = '', old_behaviour = 0 )    

        # res_loop('LiFePO4v1.n3Li1v1mp', '1m',     [4],  show = '', old_behaviour = 0 )    
        # res_loop('LiFePO4v1.n3Li1v1ms', '1m',     [4],  show = '', old_behaviour = 0 )  


    if 0: #checking influence of MAXMIX
        ''
        # res_loop('LiFePO4.nebLi1v1','1u',           [5,],  show = 'en', old_behaviour = 0 )    
        # res_loop('LiFePO4v1.n3Li1v1mp','1uMMD',     [5,],  show = 'en', old_behaviour = 0 )    


    if 0: # cheking influence of restart
        ''
        
        #result - doesnot help a lot
        #TODO - make possible restart with  magnetic moments initial for parent calc

        # res_loop('LiFePO4.nebLi1v1',   '1ulong', [1], up = '1', show = 'fo',  )   
        # res_loop('LiFePO4.nebLi1v1',   '1ulong', [1], up = '1', show = 'fo',  )   
        # res_loop('LiFePO4.nebLi1v1',   '1ulong', [4], up = '1', show = 'fo',  )   
        # res_loop('LiFePO4.nebLi1v1',   '1u', [1,2], up = '1', show = 'mep',  )   
        # res_loop('LiFePO4.nebLi1v1m0', '1u', [1,2], up = '1', show = 'mep',  )  
        # res_loop('LiFePO4.nebLi1v1m1', '1u', [1,2], up = '1', show = 'mag',  )   


        # res_loop('LiFePO4.nebLi1v1.ifn',   '1u', [1,2], up = 'up1', show = 'mep', inherit_option = 'full_nomag' )   
        # res_loop('LiFePO4.nebLi1v1m0.ifn', '1u', [1,2], up = 'up1', show = 'mep', inherit_option = 'full_nomag' )  
        # res_loop('LiFePO4.nebLi1v1m1.ifn', '1u', [1,2], up = 'up1', show = 'mag', inherit_option = 'full_nomag' ) 


    if 0: #compare 1urNM and 1urm
        ''    
        # res_loop('LiFePO4.nebLi1v1m0','1urNM', verlist = [2], up = 'up1', show = 'fo', old_behaviour = 0)
        # res_loop('LiFePO4.nebLi1v1m0','1urm',  verlist = [2], up = 'up1', show = 'en', old_behaviour = 1)
        # res_loop('LiFePO4.nebLi1v1m1','1urm', verlist = [2], up = 'up1', show = 'fo', old_behaviour = 0)



    # res_loop('LiFePO4','9uk1', 1, up = 'up1',)    #on 2016-09-03





    #122
    # for ise in [ '1u', '1uA' ]:
    for ise in [ '1m' ]:
        ' '
        # add_neb(calc['LiFePO4', '9u', 1], ise_new = ise, images = 3, corenum = 15, replicate = (1,2,2), it_new_folder = 'LiFePO4/122/neb/')
    # res_loop('LiFePO4.neb122Li1v1',   '1urWF', [1,], up = '1', show = 'me',  )   
    # res_loop('LiFePO4.neb122Li1v1m1',  '1urm', [1,2], up = '1', show = 'mep',  )   
    # res_loop('LiFePO4.neb122Li1v1m0', '1urm', verlist = [1, 2,], up = 'up1', show = 'fo', choose_image = 2, choose_outcar = 3 )    #on 2016-09-09

    # res_loop('LiFePO4v1.n3122Li1v1ms',[u'1u'],  verlist = [1], show ='mep' )
    # res_loop('LiFePO4v1.n3122Li1v1ms',[u'1uA'], verlist = [1], show ='occ' )


    # res_loop('LiFePO4v1.n3122Li1v1ms',[u'1m'], verlist = [1,], up = 'up1', show = 'mep')    #on 2016-09-26




    #NaFePO4, 

    # add_neb(calc['NaFePO4.pnma.su', '4ui', 100], ise_new = '1m', images = 3, corenum = 15, 
    #     inherit_option = None, mag_config = None, replicate = None, inherit_magmom = False)


    # res_loop('NaFePO4.pnma.n3Na1v1mp',[u'1ulong'], verlist = [1,2,4], up = 'up1', show = 'mepfo' )    #on 2016-09-15

    # res_loop('NaFePO4.pnma.suv100.n3Na1v1mp',[u'1u'], verlist = [4,], show = 'me' )    #
    # res_loop('NaFePO4.pnma.suv100.n3Na1v1ms',[u'1ulong'], verlist = [1,2,3,4,5], show = 'mep' )    #on 2016-09-28
    # res_loop('NaFePO4.pnma.suv100.n3Na1v1ms',[u'1uA'], verlist = [1,2,3,4,5], show = 'mep' )    #on 2016-09-28
    # res_loop('NaFePO4.pnma.suv100.n3Na1v1ms',[u'1ur30'], verlist = [1,2,3,4,5], show = 'mep' )    #on 2016-09-28
    # res_loop('NaFePO4.pnma.suv100.n3Na1v1ms', [u'1m'], [1, 2,3,4,5], show = 'fomep'  )     # comment = None, on 2016-09-30  


    # sys.exit()





    #KFePO4
    # add_neb(calc['KFePO4.pnma.su', '4uiWF', 100], ise_new = '1u', 
    #     images = 3, corenum = 15, i_void_start = 1, i_void_final = 2,
    #     inherit_magmom = False)

    # res_loop('KFePO4.pnma.suv100.n3K1v1ms',[u'1u'], verlist = [1],show = 'mep' )    #on 2016-09-23


    # sys.exit()



    """DS  """
    #
    # res_loop('FePO4.pnma','1u', 1, up = 'up1', input_geo_format = 'vasp', it_folder = 'LiFePO4', mat_proj_id = 'mp-777026')    #on 2016-09-03
    # res_loop('FePO4.pnma2','1u', 1, up = 'up1', input_geo_format = 'mat_proj', it_folder = 'LiFePO4', mat_proj_id = 'mp-25001')    #on 2016-09-03
    if 1:
        for i_s in [0,1]:
            for i_f in [0,2]:
                if i_s == 0 and i_f == 2: #stupid path
                    continue
                atom_to_insert = 'Li'
                # add_neb(calc['FePO4.pnma2', '1u', 1], ise_new = '1m', 
                #     images = 3, corenum = 15, atom_to_insert = atom_to_insert, i_void_start = i_s, i_void_final = i_f,
                #     r_impurity = 1.6 , it_new_folder = atom_to_insert+'FePO4/neb/',
                #     inherit_magmom = False)

                # res_loop('FePO4.pnma2.n3i'+str(i_s)+'e'+str(i_f)+atom_to_insert, [u'1'], verlist = [1,2], up = 'up1', show = 'm', old_behaviour = 0 )    #on 2016-09-15
                # res_loop('FePO4.pnma21.n3i'+str(i_s)+'e'+str(i_f)+atom_to_insert, [u'1'], verlist = [1], up = 'up1', show = 'mep', old_behaviour = 0 )    #on 2016-09-15
                # res_loop('FePO4.pnma21.n3i'+str(i_s)+'e'+str(i_f)+atom_to_insert+'ms', [u'1m'], verlist = [1], up = 'up1', show = 'mep', old_behaviour = 0 )    #on 2016-09-15
                # res_loop('FePO4.pnma2.n7i'+str(i_s)+'e'+str(i_f)+atom_to_insert+'ms', [u'1usv'], verlist = [1,2], up = 'up1', show = 'mep' )    #on 2016-09-15
                # res_loop('FePO4.pnma2.n7i'+str(i_s)+'e'+str(i_f)+atom_to_insert+'ms', [u'1ulong'], verlist = [2], up = 'up1', show = 'mep', old_behaviour = 0 )    #on 2016-09-15

    # for ise in ['1ur', '1urWF', '1uA', '1uLF' ]:
    #     res_loop('FePO4.pnma2v1.n3i1e2Lims', ise, verlist = [1], up = 'up1', show = 'mep', old_behaviour = 0 )    #on 2016-09-15
    atom_to_insert = 'Li'

    # add_neb(calc['FePO4.pnma2', '1u', 1], ise_new = '1ur30', 
    #     images = 3, corenum = 15, atom_to_insert = atom_to_insert, i_void_start = 1, i_void_final = 2,
    #     r_impurity = 1.6 , it_new_folder = atom_to_insert+'FePO4/neb/',
    #     inherit_magmom = False)

    # sys.exit()
    #Li

    # res_loop('FePO4.pnma2.n3i1e2Lims', '1u', verlist = [1], up = 'up1', show = 'mep', old_behaviour = 0 )    #on 2016-09-15
    # res_loop('FePO4.pnma2v1.n3i1e2Lims', '1ur', verlist = [1], up = 'up1', show = 'mep', old_behaviour = 0 )    #on 2016-09-15
    # res_loop('FePO4.pnma2v1.n7i1e2Lims',[u'1ur30'], verlist = [1], show = 'mep' )    #on 2016-09-22
    # add_neb(calc['FePO4.pnma2', '1u', 1], ise_new = '1u', images = 3, corenum = 15, 
    #     i_void_start = 1, i_void_final = 2, atom_to_insert = 'Li',
    #     r_impurity = 1.6 ,
    #     inherit_magmom = False, it_new_folder = 'LiFePO4/neb/')

    # res_loop('FePO4.pnma2v1.n3i1e2Lims',[u'1u'], verlist = [1], show ='mep' )    #on 2016-09-26

    # res_loop('FePO4.pnma2.n3i1e2Li',[u'1u'], verlist = [1], up = 'up1', show = 'mep')
    # res_loop('FePO4.pnma2.n7i1e2Lims',[u'1ulong'], verlist = [1,2,3,4,5,6,7,8,9], analys_type = 'neb', show = 'mep')
    # res_loop('FePO4.pnma2v1.n7i1e2Lims',[u'1ur30'], verlist = [1,2,3,4,5,6,7,8,9], analys_type = 'neb', show = 'mep')

    # res_loop('FePO4.pnma2v1.n7i1e2Lims.if',[u'0m'], 6, show = 'fo' )    #???

    # res_loop('FePO4.pnma2v1.n3i1e2Lims', [u'1ur30'], [1, 2,3,4,5], choose_outcar = None, show = 'mag', analys_type = 'neb'  )     # for OMC tutorial!!!









    #122
    # add_neb(calc['FePO4.pnma2', '1u', 1], ise_new = '1ur30s', images = 3, corenum = 15, 
    #     replicate = (1,2,2), i_void_start = 2, i_void_final = 2, atom_to_insert = 'Li',
    #     r_impurity = 1.6 ,
    #     inherit_magmom = False, it_new_folder = 'LiFePO4/122/neb/')

    # res_loop('FePO4.pnma2v1.n3122i2e2Lims',[u'1u'], verlist = [4,], show = 'mep' )    #on 2016-09-26
    # res_loop(u'FePO4.pnma2v1.n7i1e2Lims', u'1ur30', 6, show = 'occ' )    #on 2016-09-26
    # add_loop('FePO4.pnma2v1.n3122i2e2Lims',[u'1u'], verlist = [4,], inherit_option = 'occ', ise_new = '1uos', 
    #     id_from = ('FePO4.pnma2v1.n7i1e2Lims', '1ur30', 6), it_folder = 'LiFePO4/inherited',)

    # add_loop('FePO4.pnma2v1.n3122i2e2Lims.ifo', [u'1uos'], [4],up = 'up1', show = 'fo'  )     # comment = None, on 2016-09-28  
    # sys.exit()
    # res_loop('FePO4.pnma2v1.n3122i2e2Lims.ifo', [u'1uos'], [4], up = 'up2', show = 'occ'  )     # comment = None, on 2016-09-29  



    """hybrid functional!"""

    # add_loop('LiFePO4','1hss', 1,)  

    # res_loop('LiFePO4','1hss', 1,)  
    # res_loop('LiFePO4','1u', 1,)  


    # add_loop('FePO4.pnma2','1hss', 1)

    # add_neb(calc['FePO4.pnma2', '1u', 1], ise_new = '1hse', images = 3, corenum = 15, 
    #     i_void_start = 1, i_void_final = 2, atom_to_insert = 'Li',
    #     r_impurity = 1.6 ,
    #     inherit_magmom = False, it_new_folder = 'LiFePO4/neb/')





    #Na
    #     res_loop('FePO4.pnma2v1.n3i1e2Nams', '1u', verlist = [1], up = 'up1', show = 'mep', old_behaviour = 0 )    #on 2016-09-15

    #K 
    # atom_to_insert = 'K'
    # add_neb(calc['FePO4.pnma2', '1u', 1], ise_new = '1u', 
    #     images = 3, corenum = 15, atom_to_insert = atom_to_insert, i_void_start = 1, i_void_final = 2,
    #     r_impurity = 1.6 , it_new_folder = atom_to_insert+'FePO4/neb/',
    #     inherit_magmom = False)

    # res_loop('FePO4.pnma2v1.n3i1e2Kms',[u'1u'], verlist = [1,2,3,4,5], up = 'up1', show ='fo')



    # sys.exit()



    if 0: #checking influence of compression on neb
        ''
        atom_to_insert = 'Li'

        # add_neb(calc['FePO4.pnma2.su', '1u', 1], ise_new = '1u', 
        #                 images = 3, corenum = 15, atom_to_insert = atom_to_insert, i_void_start = 1, i_void_final = 2,
        #                 r_impurity = 1.3 , it_new_folder = atom_to_insert+'FePO4/neb/',
        #                 inherit_magmom = False)

        # add_neb(calc['FePO4.pnma2.su', '1u', 2], ise_new = '1u', 
        #                 images = 3, corenum = 15, atom_to_insert = atom_to_insert, i_void_start = 2, i_void_final = 3,
        #                 r_impurity = 1.4 , it_new_folder = atom_to_insert+'FePO4/neb/',
        #                 inherit_magmom = False)

        # res_loop('FePO4.pnma2.suv1.n3i1e2Lims', [u'1u'], verlist =  [4], up = 'up1', show = 'me')  
        # res_loop('FePO4.pnma2.n3i1e2Li',        [u'1u'], verlist =  [4], up = 'up1', show = 'me')  
        # res_loop('FePO4.pnma2.suv2.n3i2e3Lims',[u'1u'], verlist = [1,2,3,4,5], up = 'up1', show = 'mep' )    #on 2016-09-23
        
        # sys.exit()



    #LiVPO4F
    # res_loop('LiVPO4F.nebLi1v1','1ur', verlist = [1,2,3,4,5], up = 'up1', show = 'fo', old_behaviour =0)    #on 2016-09-08




    #LiMnPO4
    for it in ['LiMnPO4', 'LiNiO2', ]:
        ''
        # add_neb(calc[it, '9u', 1], ise_new = '1u', images = 3, corenum = 15, inherit_option = None, mag_config = None, inherit_magmom = False)
        # res_loop(it+'v1.n3Li1v1ms',[u'1u'], verlist = [1,2,3,4,5], up = 'up1', show = 'm')    #on 2016-09-23


    # add_neb(calc['LiCoPO4', '8Utmr2-8', 1], ise_new = '1u', images = 3, corenum = 15, inherit_option = None, mag_config = None, inherit_magmom = False)

    #LiMnPO4
    # res_loop('LiCoPO4v1.n3Li1v1ms',[u'1u'], verlist = [1,2,3,4,5], up = 'up1',  )    #on 2016-09-23
    # res_loop('CoPO4v1.n3i3e2Nams', [u'1u'], [1, 2,3,4,5], show = 'mep'  )     # comment = None, on 2016-09-29  



    #DS
    # add_neb(calc['CoPO4', '8Utmr4-2', 1], ise_new = '1u', 
    #     images = 3, corenum = 15,
    #     atom_to_insert = 'Na', i_void_start = 3, i_void_final = 2,
    #     r_impurity = 1.8 , it_new_folder = 'LiCoPO4/neb/'
    #     )

    # sys.exit()
    #NaMnAsO4
    #IS
    # add_neb(calc['NaMnAsO4', '8u', 1], ise_new = '1u', images = 3, corenum = 15, inherit_option = None, mag_config = None, inherit_magmom = False)
    # res_loop('NaMnAsO4v1.n3Na17v1ms', [u'1u'], [1, 2,3,4,5], show = 'mep'  )     # comment = None, on 2016-09-29  


    #DS
    # add_neb(calc['MnAsO4', '8u', 1], ise_new = '1u', 
    #     images = 3, corenum = 15, inherit_option = None, mag_config = None, inherit_magmom = False,
    #     atom_to_insert = 'Na', i_void_start = 0, i_void_final = 0,
    #     r_impurity = 1.6 , it_new_folder = 'NaMnAsO4/neb/'
    #     )

    # res_loop('MnAsO4v1.n3i0e1Nams', [u'1u'], [1, 2,3,4,5], show = 'mep'  )     # comment = None, on 2016-09-29  
    # res_loop('MnAsO4v1.n3i0e0Nams', [u'1u'], [1, 2,3,4,5], show = 'fomep'  )     # comment = None, on 2016-09-29  


    # sys.exit()


    #some additional
    # add_loop('Na2FePO4F', '1u', 1, input_geo_format = 'cee_database', it_folder = 'Na2FePO4F')
    # add_loop('KFeSO4F', '1u', 1, input_geo_format = 'cee_database', it_folder = 'KFeSO4F')
    # add_loop('Na2FeVF7', '1u', 1, input_geo_format = 'cee_database', it_folder = 'Na2FeVF7')











    """Charge den analysis"""
    calc = header.calc
    # print calc['LiFePO4.nebLi1v1',   '1u', 1].get_chg_file('CHG')
    # print calc['LiFePO4v1.n3Li1v1ms',   '1m', 1].get_chg_file('CHG')
    # print calc['LiFePO4.nebLi1v1.if', '0m', 1].get_chg_file('CHG')
    # print calc['FePO4.pnma2v1.n7i1e2Lims', '1ur30', 1].get_chg_file('CHG')

    # print (calc['LiFePO4.nebLi1v1',   '1u', 4].get_chg_file('CHG'))
    # print (calc['LiFePO4.nebLi1v1.if','0m', 4].get_chg_file('CHG', 'asoutcar'))

    # print (calc['LiFePO4.nebLi1v1.ifi',   '0u', 1].get_chg_file('CHG', 'asoutcar'))
    # print (calc['LiFePO4',   '1u', 1].get_chg_file('CHG', 'asoutcar'))


    # add_loop('LiFePO4.nebLi1v1','1u', 1, up = 'up1', ise_new = '0m', inherit_option = 'full' )
    # add_loop('LiFePO4.nebLi1v1','1u', 4, up = 'up1', ise_new = '0m', inherit_option = 'full', savefile = 'ocv', override = 1 )
    # res_loop('LiFePO4.nebLi1v1',   '1u', 1, show = 'mag', )
    # res_loop('LiFePO4.nebLi1v1',   '1u', 1, show = 'mag', )
    # res_loop('LiFePO4','1u', 1, show = 'occ', )


    # add_loop('FePO4.pnma2v1.n7i1e2Lims','1ur30', 6, up = 'up1', ise_new = '0m', inherit_option = 'full' )


    # inherit_icalc('full', 'LiFePO4.nebLi1v1.ifi', 1, ('LiFePO4.nebLi1v1','1u',1), header.calc, id_base_st_type = 'init', it_folder = 'LiFePO4' )    #on 2014-02-26
    # add_loop('LiFePO4.nebLi1v1.ifi', '0u', 1, savefile = 'cov')
    # res_loop('LiFePO4.nebLi1v1.ifi', '0u', 1, show = 'mag')





    # 

    # add_loop('LiFePO4.nebLi1v1','1u', 1, up = 'up1', ise_new = '1urseq', calc_method = 'uniform_scale', it_folder = 'LiFePO4/uscale/', scale_region = [-2,6],  inherit_option = 'inherit_xred' )
    # sys.exit()

    # add_loop('LiFePO4.nebLi1v1','1useq', [1,2], up = 'up1', ise_new = '1urseq', calc_method = 'neb', n_neb_images = 3 )
    # sys.exit()

    """OMC analysis"""
    # res_loop('LiFePO4.nebLi1v1','1u', 4, show = 'occ')

    # res_loop('FePO4.pnma2v1.n3i1e2Lims',[u'1u'],    verlist = [4], show ='occ21' ) 


    # res_loop('FePO4.pnma2.n7i1e2Lims', [u'1ulong'], verlist = [1], show = 'oc')
    # res_loop('FePO4.pnma2v1.n7i1e2Lims',[u'1ur30'], verlist = [6], show = 'occ' ) 



    # add_loop('FePO4.pnma2.n7i1e2Lims','1ulong', 1, up = 'up1',  inherit_option = 'occ', ise_new = '1uos', 
    #     id_from = ('FePO4.pnma2v1.n7i1e2Lims', '1ur30', 1), it_folder = 'LiFePO4/inherited' )
    # res_loop('FePO4.pnma2.n7i1e2Lims.ifo','1uomc', 1, up = 'up2', show = 'occ' )
    # res_loop('FePO4.pnma2.n7i1e2Lims.ifo','1uos', 1, show = 'occ' )



    # res_loop('LiCoO2.s10','1urs',1,   choose_outcar = 0,    show = 'maga')
    # res_loop('LiCoO2.s10.su','4uis',100, choose_outcar = 0, show = 'maga')




    # res_loop('NiO2.id.su.s10.su','4uris',100, choose_outcar = 3, show = 'oc')
    # res_loop('NiO2.id.su.s10.su','4uris',100, choose_outcar = 4, show = 'oc')

    # res_loop('NiO2.id.su','4uis',  4, choose_outcar = 2, show = 'smag')
    # res_loop('NiO2.id.su','4uisWF',4, choose_outcar = 2, show = 'smag')
    # res_loop('NiO2.id.su','4uis',4, choose_outcar = 2, show = 'smag')

    # res_loop('NiO2.id.su','4uisC0W1',4, choose_outcar = 1, show = 'smag')
    # res_loop('NiO2.id.su','4uisWF', list(range(1,7))+[100],  choose_outcar = 2, show = 'fi', analys_type = 'fit_a')
    # res_loop('NiO2.id.su','4uisC0W1', list(range(1,7))+[100], choose_outcar = 1, show = 'fo', analys_type = 'fit_a')
    # ('NiO2.id.su', '4uisC0W1', 7)
    # ('NiO2.id.su', '4uisC0W2', 7)

    # print(calc['NiO2.id.su','4uis',4].dir)



    # res_loop('LiFePO4.s10.suv100.n3Li1v1ms','1urs',4, choose_outcar = 3, show = 'occ102')
    # res_loop('LiFePO4.s10.suv100.n3Li1v1ms','1urs',4, choose_outcar = 4, show = 'occ102')


    # add_loop('Na111', '8u', 1, up = 'up1', input_geo_format = 'vasp', it_folder = 'Na/bcc')




    """Space groups"""

    # res_loop('LiNiO2.r3m','1u', 1, input_geo_format = 'mat_proj', it_folder = 'LiNiO2/r3m', mat_proj_id = 'mp-25592')    #on 2016-09-03
    # add_loop('NiO2.r3m','1u', 1, input_geo_format = 'mat_proj', it_folder = 'LiNiO2/r3m', mat_proj_id = 'mp-35925')  
    # res_loop('NiO2.r3m','1u', 1, input_geo_format = 'mat_proj', it_folder = 'LiNiO2/r3m', mat_proj_id = 'mp-35925', show = 'occ1')  
    # res_loop('NiO2.r3m.id.su','4urisC0W1',100, show = 'occ1')
    # res_loop('NiO2.r3m.id.su','4uiAl',100, show = 'occ1')
    # rprimd = calc['KVPO4F.Pna21','1u', 1].init.rprimd
    # print(rprimd)
    # print(calc['NiO2.r3m.id.su','4urisC0W1', 1].end.rprimd)
    # res_loop('LiTiS2.s10.su','4uis',100, show = 'occ1' )
    # from geo import ortho_vec2, ortho_vec
    # print(ortho_vec(rprimd, [10,10,10])  )
    # print(ortho_vec2(rprimd, [10,10,10])  )

    # add_loop('LiMn2O4.Fd3m.111','1u', 1, input_geo_format = 'mat_proj', it_folder = 'LiMn2O4/Fd3m', mat_proj_id = 'mp-25015', mat_proj_cell = 'conv')    #on 2016-09-03
    # res_loop('LiMn2O4.Fd3m.111','1u', 1, show = 'fomag')    #on 2016-09-03
    # res_loop('LiMn2O4.Fd3m','1urs', 1, choose_outcar  = 3, show = 'mag')    #on 2016-09-03
    # res_loop('LiMn2O4.Fd3m','1u', 1, show = 'fomag')    #on 2016-09-03
    # res_loop('LiMn2O4.Fd3m.s10','1urs', 1, show = 'fomag', choose_outcar = 3)    #on 2016-09-03
    # res_loop('LiMn2O4.Fd3m.s10','1ur10', 1, show = 'fomag', choose_outcar = 10)    #on 2016-09-03
    # print (calc['LiMn2O4.Fd3m','1u',1].end.get_space_group_info() )
    # print (calc['LiMn2O4.Fd3m.s10','1ur10',1].end.get_space_group_info() )
    # print (calc['LiMn2O4.Fd3m.111','1u', 1,].end.get_space_group_info() )
    # res_loop('LiMn2O4.Fd3m.111','1u', 1,)


    # print (calc['Na2FePO4F.s10','1urs', 1,].end.get_space_group_info() )

    # add_loop('KVPO4F.Pna21','1u', 1, input_geo_format = 'cee_database', it_folder = 'KVPO4F/Pna21', cee_file = 'x1aa_exp.cif')    #on 2016-09-03
    # res_loop('KVPO4F.Pna21','1u', 1, input_geo_format = 'cee_database', it_folder = 'KVPO4F/Pna21')    #on 2016-09-03
        # add_loop('Na2FePO4F', '1u', 1, input_geo_format = 'cee_database', it_folder = 'Na2FePO4F')
    # print (calc['KVPO4F.Pna21','1u', 1].init.get_space_group_info() )


    # print (calc['LiVPO4F.s10.su','4uis',100].end.get_space_group_info() )
    # print (calc['LiVP2O7.s10.su','4uis',100].end.get_space_group_info() )
    # print (calc['LiNiO2.s10.su','4uisns',100].end.get_space_group_info() )

    # print (calc['NiO2.id.su','4uisC0W1',1].init.get_space_group_info() )
    # print (calc['NiO2.id.su','4uisC0W1',1].init.get_space_group_info() )

    # print (calc['LiTiO2.s10.su','4uisns',100].end.get_space_group_info() )
    # print (calc['TiS2.id.su','4uis',100].end.rprimd )
    # print (calc['TiS2','8Utm3-5',1].end.rprimd )
    # sys.exit()
    # print (calc['LiTiS2.s10.su','4uis',100].end.get_space_group_info() )
    # print (calc['LiMn2O4.s10.su','4uis',100].end.get_space_group_info() )
    # print (calc['LiTiS2.s10.su','4uis',100].end.get_space_group_info() )
    # print (calc['KVPO4F.ir.su','4uisC0W1',100].end.get_space_group_info() )

    # print (calc['Li2FePO4F.ir.su.s10.su','4uis',100].end.get_space_group_info() )


    # res_loop('LiFePO4.s10.suv100.n3Li1v1ms','1urs',4, choose_outcar = 3   , show = 'occ97')
    # res_loop('LiFePO4.s10.suv100.n3Li1v1ms','1urs',4, choose_outcar = None, show = 'occ97')
     










    """For notebooks"""

    # sys.path.extend(['/home/aksenov/Simulation_wrapper/siman', '/home/aksenov/Simulation_wrapper/SSHTools'])
    # from SSHTools import SSHTools
    # header.ssh_object = SSHTools()
    # header.ssh_object.setup(user="aksenov",host="10.30.16.62",pkey="/home/aksenov/.ssh/id_rsa")
    # header.PATH2PROJECT    = 'topologic'
    # header.varset['static'].potdir = {83:'Bi_pv', 34:'Se'}
    # header.siman_run = 0
    # add_loop('Bi2Se3', 'static', 1, input_geo_file = 'Bi2Se3_mp-541837_computed.cif', it_folder = 'Bi2Se3', run = 0)

    # res_loop('Bi2Se3', 'static', 1, up = 'up')
    # header.warnings = 'yY'
    # add_loop('Bi2Se3', 'static', 1, calc_method = 'uniform_scale', scale_region = (-5, 5), run = 0)



    from siman.dos_functions import plot_dos
    # add_loop('Li2', 'dosb4', 1, run = 1)
    # res_loop('Li2', 'dosb4', 1, )
    # add_loop('Li2', 'fullb4sdos', 1, run = 1)
    # plot_dos(calc['Li2', 'dosb4', 1], dostype = 'partial', orbitals = 'd')
    # plot_dos(calc['Li2', 'fullb4sdos', 1], dostype = 'partial', orbitals = 's', up = 1, show =0 )

    # header.PATH2PROJECT    = 'topologic'

    # add_loop('Bi2Se3', 'dos', 1, input_geo_file = 'Bi2Se3_mp-541837_computed.POSCAR', run = 0)


    # res_loop('Bi2Se3', 'dos', 1)
    # plot_dos(header.calc['Bi2Se3', 'dos', 1], dostype = 'total')




    # add_loop('Li2', 'fullb4', 1, run =0)
    # res_loop('Li2', 'fullb4', 1, show = 'fo', up = 'up')
    # add_loop('Li2', 'fullb4', 1, ise_new = 'ion', inherit_option = 'supercell', mul_matrix = [[3,0,0],[0,3,0],[0,0,3]], run = 1)
    # res_loop('Li2.s333', ['ion'], [1], up = 'up2', show = 'fo' )
    # add_loop('Li2.s333', ['ion'], [1], inherit_option = 'make_vacancy', i_atom_to_remove = 0, run = 1, override = 1)

    # res_loop('Li2.s333.vac', ['ion'], [1], up = 'up2', show = 'fo', analys_type = 'matrix_diff', b_id = ('Li2.s333', 'ion', 1))
    from impurity import add_impurity
    # add_impurity('Li2.s333.oct', 'Li', calc = calc, it_to = 'Li2.s333', put_exactly_to = (0, 1/6, 1/6))

    # add_loop('Li2.s333.tet', ['ion'], [1], run = 1, override = 1)
    # res_loop('Li2.s333.tet', ['ion'], [1], show = 'fo', analys_type = 'matrix_diff', b_id = ('Li2.s333', 'ion', 1))
    # add_loop('Li2.s333.oct', ['ion'], [1], run = 1, override = 1)
    # res_loop('Li2.s333.oct', ['ion'], [1], show = 'fo', analys_type = 'matrix_diff', b_id = ('Li2.s333', 'ion', 1))

    from impurity import find_pores
    # st = calc['Li2.s333', 'ion', 1].end
    # st_pores = find_pores(calc['Li2.s333', 'ion', 1].end, r_matrix = 1.0, r_impurity = 0.5, fine = 1, calctype = 'all_pores')
    # write_xyz(st.add_atoms(st_pores.xcart, 'H'), file_name = st.name+'_possible_positions')




    def test_energy_convergence():
        """Comparison of two Li2FePO4 structures obtained with U-ramping and without it"""


        # res_loop('Li2FePO4F.ir.su.s10.su.ifn','1ur30', 100, show = 'occ35')
        # res_loop('Li2FePO4F.ir.su.s10.su',    '4uis',  100, show = 'occ35')

        # add_loop('Li2FePO4F.ir.su.s10.su.ifn','1ur30', 100, ise_new = '0ur', inherit_option = 'full', savefile = 'ocv', override = True)
        # add_loop('Li2FePO4F.ir.su.s10.su',    '4uis',  100, ise_new = '0ur', inherit_option = 'full', savefile = 'ocv', override = True)
        res_loop('Li2FePO4F.ir.su.s10.su.if',     ['0ur'], [100], show = 'occ', choose_outcar = 3  )
        res_loop('Li2FePO4F.ir.su.s10.su.ifn.if', ['0ur'], [100], show = 'occ', choose_outcar = 3  ) # relaxed in some high energy state
        # print(calc['Li2FePO4F.ir.su.s10.su.if',     '0ur', 100].get_chg_file('100.U00.CHGCAR'))
        # print(calc['Li2FePO4F.ir.su.s10.su.ifn.if',     '0ur', 100].get_chg_file('100.U00.CHGCAR'))

        # res_loop('Li2FePO4F.ir.su',    '4uis',  100)
        # add_loop('Li2FePO4F.ir.su',    '4uis',  100, ise_new = '1urs', inherit_option = 'full', savefile = 'ocv', override = True)

        # add_loop('Li2FePO4F.ir.su.if',    '1urst5',  100)
        # res_loop('Li2FePO4F.ir.su.if',    '1urs',    100, show = 'fo')
        # res_loop('Li2FePO4F.ir.su.if',    '1urst5',  100, show = 'fo', choose_outcar = 3)
        # res_loop('Li2FePO4F.ir.su.if',    '1urst5',  100, show = 'fo', choose_outcar = 4)

        # add_loop('Li2FePO4F.ir.su.if',    '1ursns',  100)


    # test_energy_convergence()











    """Adding to database"""
    # from calc_manage import add_to_database
    # cl = calc['Li2FePO4F.ir.su.s10.su','4uis',100]
    # print(cl.end.natom)
    # print (st.get_space_group_info())
    # add_to_database(cl)

    # calc2['asdf'] = 'test'
    # print(calc2['asdf'])








    """Additional task"""
    if 0:
        from siman.set_functions import InputSet
        folder = '/home/aksenov/scientific_projects/cathode/azh_task/'
        ise = 'azh30'
        s = InputSet('azh30', folder+'POTCAR')
        # varset[ise] =
        s.read_incar(folder+'/INCAR')
        s.kpoints_file = folder+'/KPOINTS'
        s.vasp_params['LDAUL'] = {'Co':2}
        s.vasp_params['LDAUU'] = {'Co':5.7}
        s.vasp_params['LDAUJ'] = {'Co':0}
        s.u_ramping_nstep = 30
        s.save_last_wave = True
        # manually_remove_from_struct_des(struct_des, 'CoPO4F')

        # add_loop('CoPO4F', 'azh', 1, input_geo_file = 'azh_task/POSCAR',  it_folder = 'azh_task')
        # res_loop('CoPO4F', '0ur', 1, input_geo_file = 'azh_task/POSCAR',  it_folder = 'azh_task')

        # res_loop('CoPO4F', 'azh30', 1, input_geo_file = 'azh_task/POSCAR',  it_folder = 'azh_task', show = 'mag')
        # res_loop('CoPO4F', '0u57r', 1, input_geo_file = 'azh_task/POSCAR',  it_folder = 'azh_task', show = 'mag')









        Li2Co_pbcn = calc['Li2CoPO4F.pbcn.su','4uis',100] #'Li2CoPO4F.pbcn.su.s10'
        # v_LiCo = calc['LiCoPO4F.pbcn.id1.su.s10','1u',100] 

        cl = Li2Co_pbcn
        it_folder = header.struct_des[cl.id[0]].sfolder+'/chg'
        for del_pos in 1, 2:
            it_new = cl.id[0]+'.id'+str(del_pos)
            # st = create_deintercalated_structure(cl.end, 'Li', del_pos = del_pos)
            # add_loop(it_new,'0u',1, input_st = st, it_folder = it_folder, override = True) #charge density without relaxation
            res_loop(it_new,'0u',1)        


        v_LiCo_pbcn = calc['Li2CoPO4F.pbcn.su.id1','0u',1] 
        cal_chg_diff(Li2Co_pbcn, v_LiCo_pbcn, wcell = 0)










    from analysis import calc_redox


    # calc_redox(calc['KVPO4F.Pna21.s10.su','4uis',100], calc['VPO4F.Pna21.id.su.s10.su','4uis',100])




    from geo import determine_symmetry_positions

    # st = calc['NaFePO4F.id1.su.s10.su','4uis',100].end
    # st = calc[('KVPO4F.Pna21', '1u', 1)].end
    # st = calc['LiFePO4F.pnma.id2.su.s10','1u',100].end
    # determine_symmetry_positions(st, 'Fe')



    from geo import create_antisite_defect2

    st1 = calc['NaFePO4F.id1.su.s10.su','4uis',100].init
    st2 = calc['NaFePO4F.id2.su.s10.su','4uis',100].end
    # print(st1.convert2pymatgen())
    # print()


    # st = create_antisite_defect2(st_base = st2, st_from = st1, cation = 'Na', trans = 'Fe', mode = 'add_swp')







    from calc_manage import smart_structure_read



    # st  = calc['Li2CoPO4F.pnma','1u', 1].end
    # print(calc['Li2CoPO4F.pnma','1u', 1].end.get_elements())
    # from geo import create_deintercalated_structure
    # create_deintercalated_structure(st, 'Li', del_pos = 1)



    # add_loop('Li2CoPO4F.pnma','1u', 1, input_geo_format = 'mat_proj', it_folder = 'Li2CoPO4F/pnma', mat_proj_id = 'mp-770853')  

    # add_loop('Li2CoPO4F.pbcn','1u', 1, input_geo_format = 'mat_proj', it_folder = 'Li2CoPO4F/pbcn', mat_proj_id = 'mp-770624')  

    # add_loop('Li2CoPO4F.pnma','1ur', 1, input_geo_format = 'mat_proj', it_folder = 'Li2CoPO4F/pnma', mat_proj_id = 'mp-770853', run = 1)  


    # res_loop('Li2CoPO4F.pbcn','1u', 1)
    # res_loop('Li2CoPO4F.pnma','1u', 1)
    # res_loop('Li2CoPO4F.pnma','1ur', 1, show = 'occ')

    # st = smart_structure_read(input_geo_file = 'Li2FePO4F/Li2FePO4F_abakumov.cif')

    # add_loop('Li2FePO4F.abak','1u', 1, input_geo_file = 'Li2FePO4F/Li2FePO4F_abakumov.cif', it_folder = 'Li2FePO4F/pnma' )  
    # res_loop('Li2FePO4F.abak','1u',1, show = 'fo')
    # st = calc['Li2FePO4F.abak','1u',1].end
    # st = calc['Li2FePO4F.pnma','1u', 1].end
    # print(st.convert2pymatgen())
    # print(st.get_space_group_info(0.3))


    # add_loop('Li2FePO4F.pnma','1u', 1, input_geo_format = 'mat_proj', it_folder = 'Li2FePO4F/pnma', mat_proj_id = 'mp-776062')  

    # res_loop('Li2FePO4F.pnma','1u', 1)
    # res_loop('Li2FePO4F.ir.su','4uis',100)


    # calc['Li2FePO4F.abak','1u',1].end.jmol()
    # calc['Li2FePO4F.pnma','1u', 1].end.jmol()






    # add_loop('LiCrPO4.pnma','1u', 1, input_geo_format = 'mat_proj', it_folder = 'LiCrPO4/pnma', mat_proj_id = 'mp-25507')  
    # add_loop('LiCrPO4.pnma','1u', 1, input_geo_format = 'mat_proj', it_folder = 'LiCrPO4/pnma', mat_proj_id = 'mp-25507')  
    if 0:
        param_dic = {'id':('LiCrPO4.pnma','', 1), 'ds':'CrPO4.pnma', 'itfolder':'LiCrPO4/pnma', 'main_set':'1u'  }
        # calc_barriers('normal', up = 0, param_dic = param_dic)
        # calc_barriers('replace',  'Li', 'Na', show_fit = 0, up = 0, upA = 0, param_dic = param_dic)
        calc_barriers('make_ds',  'Li', 'Na', show_fit = 0, up = 0, upA = 0, param_dic = param_dic)



    # from SSHTools import SSHTools
    # header.ssh_object = SSHTools()
    # header.ssh_object.setup(user="aksenov",host="10.30.16.62",pkey="/home/aksenov/.ssh/id_rsa")










    # add_loop('Bi2Se3.hex', '1', 1, input_geo_file = 'Bi2Se3/Bi2Se3_hex.POSCAR', it_folder = 'Bi2Se3', run = 1)
    # add_loop('Bi2Se3.rho', '1', 1, input_geo_file = 'Bi2Se3_mp-541837_computed.POSCAR', it_folder = 'Bi2Se3', run = 0)


    # res_loop('Bi2Se3.hex.su', '1', 100)
    # res_loop('Bi2Se3.rho', '1', 1)

    # add_loop('Bi2Se3.hex', '1', 1, calc_method = 'uniform_scale', inherit_option = 'inherit_xred', scale_region = (-4, 4), run = 1)


    # print(calc['Bi2Se3.hex.su', '1', 100].energy_sigma0/15 - calc['Bi2Se3.rho', '1', 1].energy_sigma0/5)

    # res_loop('Bi2Se3.hex.su', ['1'], [1, 2, 3, 4, 5, 6, 7, 100], show = 'fitfo', analys_type = 'fit_a'  ) 

    # add_loop('Bi2Se3.rho', '1', 1, ise_new = '0bsr', inherit_option = 'full', savefile = 'ocx', override = 1, run = 1)

    # add_loop('Bi2Se3.rho', 'dos', 1, savefile = 'ocx', input_geo_file = 'Bi2Se3_mp-541837_computed.POSCAR', it_folder = 'Bi2Se3', run = 0)

    # res_loop('Bi2Se3.rho.if', '0bsr', 1   ,up = 'x')
    # res_loop('Bi2Se3.rho', 'dos', 1 ,up = 'x')

    # add_loop('Bi2Se3.rho', '1', 1, ise_new = '0bsrsoc', inherit_option = 'full', savefile = 'ocx', override = 1, run = 1)
    # add_loop('Bi2Se3.rho', '1', 1, ise_new = '0bsrsoc100', inherit_option = 'full', savefile = 'ocx', override = 1, run = 1)
    # add_loop('Bi2Se3.rho', '1', 1, ise_new = '0bsrsoc010', inherit_option = 'full', savefile = 'ocx', override = 1, run = 1)
    # add_loop('Bi2Se3.rho', '1', 1, ise_new = '0bsrsoc000', inherit_option = 'full', savefile = 'ocx', override = 1, run = 1)
    # add_loop('Bi2Se3.rho', '1', 1, ise_new = '0bsrsoc0', inherit_option = 'full', savefile = 'ocx', override = 1, run = 1)

    # res_loop('Bi2Se3.rho.if', '0bsrsoc', 1, up = 'x')
    # res_loop('Bi2Se3.rho.if', '0bsrsoc100', 1, up = 'x')
    # res_loop('Bi2Se3.rho.if', '0bsrsoc010', 1, show = 'mag')
    # res_loop('Bi2Se3.rho.if', '0bsrsoc000', 1, up = 'x')
    # res_loop('Bi2Se3.rho.if', '0bsrsoc0', 1, up = 'x')


    # from bands import plot_bands
    # ban = '/home/aksenov/scientific_projects/cathode/Bi2Se3/Bi2Se3.rho.if.0bsrsoc0/'
    # plot_bands('/home/aksenov/scientific_projects/cathode/Bi2Se3/Bi2Se3.rho.dos/1.vasprun.xml', 
    #     vasprun_bands = ban+'1.vasprun.xml', kpoints = ban+'KPOINTS', element = 'Bi', ylim = (-1.5, 2.5))




    # file = 'VO/CONTCAR_mopac'
    # # file = 'VO/Kayunkid10.cif'
    # st = smart_structure_read(input_geo_file = file)

    # st = st.return_atoms_to_cell()
    # st = st.shift_atoms(vector_red  = (0.5,0.5,0.5))
    # st = st.return_atoms_to_cell()

    # st.write_poscar(file+'_conv', vasp5 = True)



    # add_loop('Ti2', '0', 1, input_geo_file = 'Ti/hex111.geo', it_folder = 'Ti/hex/', run = 1, corenum = 1)
    # add_loop('Ti2', '1hssb', 1, input_geo_file = 'Ti/hex111.geo', it_folder = 'Ti/hex/', run = 1, corenum = 4)
    # add_loop('TiFe', '1hssb', 1, input_geo_file = 'TiFe/1.CONTCAR', it_folder = 'TiFe/B1/', run = 1, corenum = 4)
    # add_loop('Fe2', '1hssbm', 1, input_geo_file = 'Fe/1.CONTCAR', it_folder = 'Fe/bcc/', run = 1, corenum = 8)

    # res_loop('Ti2', '0', 1)


    # res_loop('Ti2', '1hssb', 1)
    # res_loop('TiFe', '1hssb', 1)
    # res_loop('Fe2', '1hssbm', 1)


    # print(-20.8478--9.1541--9.2957 )








    def k10():
        V250 = -2235.4046
        V249v = -2224.0141
        V248Tiv = -2223.2459
        V249Ti  = -2233.7968


        print(V249Ti+V249v  - (V248Tiv + V250 )  )

    def k15():
        V250 = -2235.5223
        V249v = -2224.1072
        V248Tiv = -2223.2724
        V249Ti  = -2234.0743
        print(V249Ti+V249v  - (V248Tiv + V250 )  )

    # k15()














    return #end of workflow





















def redox(run = 0, ramp = 0, readfiles = 0, input_geo_format = 'vasp',
    image_name = None, image_format = 'pdf', dpi = 300): 
    """
	Run, analyze, plot redox potentials

    """
    calc = header.calc



    cl_list7 = [
        ['LiVPO4F', '8Utmr3-1', 'VPO4F', '8Utmr3-1', 'LiVPO4F'] ,
        ['LiVP2O7', '8Utmr3-1', 'VP2O7', '8Utmr3-1', 'LiVP2O7'] ,
        ['LiTiO2', '8Utmr4-2', 'TiO2', '8Utmr4-2', 'LiTiO2'] ,
        # ['LiMn2O4', '8Utmr4-9', 'Mn2O4', '8Utmr4-9', 'LiMn2O4'] ,
        ['NaMnAsO4', '8Utmr3-9', 'MnAsO4', '8Utmr3-9', 'NaMnAsO4'] ,
        # ['Na2FePO4F', '8Utmr4', 'FePO4F', '8Utmr4', 'Na2FePO4F'] ,
        # ['KFeSO4F', '8Utmr4', 'FeSO4F', '8Utmr4', 'KFeSO4F'] ,
        # ['KVPO4F', '8Utmr3-1', 'VPO4F', '8Utmr3-1', 'KVPO4F'] ,
    ]




    cl_list8 = []
    for ise in '8u', '8ue4', '8ue5', '8uL':
    # for ise in '9u', '9uk15', '9uk1':
        cl_list8.extend([
        # ['KFeSO4F', ise, 'FeSO4F', ise, 'KFeSO4F'] ,
        # ['KVPO4F', ise, 'VPO4F', ise, 'KVPO4F'] ,
        ['LiCoO2', ise, 'CoO2', ise, 'LiCoO2'] ,
        # ['LiCoPO4', ise, 'CoPO4', ise, 'LiCoPO4'] ,
        ['LiFePO4', ise, 'FePO4', ise, 'LiFePO4'] ,
        # ['LiMn2O4', ise, 'Mn2O4', ise, 'LiMn2O4'] ,
        ['LiMnPO4', ise, 'MnPO4', ise, 'LiMnPO4'] ,
        ['LiNiO2', ise, 'NiO2', ise, 'LiNiO2'] ,
        ['LiTiO2', ise, 'TiO2', ise, 'LiTiO2'] ,
        ['LiTiS2', ise, 'TiS2', ise, 'LiTiS2'] ,
        ['LiVP2O7', ise, 'VP2O7', ise, 'LiVP2O7'] ,
        ['LiVPO4F', ise, 'VPO4F', ise, 'LiVPO4F'] ,
        # ['Na2FePO4F', ise, 'FePO4F', ise, 'Na2FePO4F'] ,
        ['NaFePO4', ise, 'FePO4', ise, 'NaFePO4'] ,
        ['NaMnAsO4', ise, 'MnAsO4', ise, 'NaMnAsO4'] ,
        ]
        )

    # for li in cl_list8:  
    #     print li
    # sys.exit()



    cl_list0 = [ #Shishkin2016
    ['LiCoPO4', '8Utmr2-8', 'LiCoPO4'] ,
    ['LiFePO4', '8Utmr2-1', 'LiFePO4'] ,
    ['FePO4', '8Utmr3-7', 'LiFePO4'] ,
    ['LiTiS2', '8Utmr3-3', 'LiTiS2'] ,
    ['TiS2', '8Utmr3-5', 'LiTiS2'] ,
    ['NiO2', '8Utmr4', 'LiNiO2'] ,
    ['MnPO4', '8Utmr4', 'LiMnPO4'] ,
    ['CoO2', '8Utmr3-9', 'LiCoO2'] ,
    ['LiNiO2', '8Utmr4-6', 'LiNiO2'] ,
    ['LiMnPO4', '8Utmr2-2', 'LiMnPO4'] ,
    ['NaFePO4', '8Utmr2-2', 'NaFePO4'] ,
    ['CoPO4', '8Utmr4-2', 'LiCoPO4'] ,
    ['LiCoO2', '8Utmr3-6', 'LiCoO2'] ,
    ['LiTiO2.s10.su','4uisns', 'LiTiO2'], 
    ['TiO2.id.su','4uis', 'LiTiO2'], 

    ['LiMn2O4.s10.su','4uis',    'LiMn2O4'], 
    ['Mn2O4.id.su.s10','1urs', 'LiMn2O4'], 

    ]


    cl_list1 = [
    # ['LiCoO2', '8Ur3-6', 'LiCoO2'] ,
    # ['CoO2', '8Ur3-9', 'LiCoO2'] ,
    # ['LiTiS2', '8Ur3-3', 'LiTiS2'] ,
    # ['TiS2', '8Ur3-5', 'LiTiS2'] ,
    ['LiFePO4', '8Ur2-1', 'LiFePO4'] ,
    ['FePO4', '8Ur3-7', 'LiFePO4'] ,
    # ['NiO2', '8Ur4', 'LiNiO2'] ,
    # ['MnPO4', '8Ur4', 'LiMnPO4'] ,
    # ['LiCoPO4', '8Ur2-8', 'LiCoPO4'] ,
    # ['LiNiO2', '8Ur4-6', 'LiNiO2'] ,
    # ['LiMnPO4', '8Ur2-2', 'LiMnPO4'] ,
    # ['NaFePO4', '8Ur2-2', 'NaFePO4'] ,
    # ['CoPO4', '8Ur4-2', 'LiCoPO4'] ,
    ]

    # print cl_list8
    # sys.exit()

    DSs = []
    results = []
    final_list = []
    # for ls in  cl_list7+cl_list0, :
    # for ls in  cl_list8+cl_list0+cl_list7, :
    

    # for ls in  cl_list8,:
    for ls in  cl_list0,:
    # for ls in  cl_list7, :
        # print ls
        for cl in ls:
            # print ('hello',cl)
            energy_ref = 0
            if 'Li' in cl[0]: energy_ref = calc['Li111', '8Utmb4', 1].energy_sigma0/2
            if 'Na' in cl[0]: energy_ref = calc['Na111', '8Utmb4', 1].energy_sigma0/2

            if len(cl) == 3:
                if not ramp: 
                    cl[1] = cl[1].replace('r', '')

                if run: 
                    add_loop(cl[0], cl[1], 1, up = 'up1', input_geo_format = 'vasp', it_folder = cl[2])
                
                else:
                    if 'Li' not in cl[0] and 'Na' not in cl[0]: 
                        continue
                    # print 'hello'


                    for clb in ls: #find deintercalated
                        # print clb == cl
                        if clb == cl: 
                            continue
                        if clb[0].split('.')[0] in cl[0]: 
                            bcl = copy.deepcopy(clb)
                            # print (clb[0], cl[0])
                    
                    # print energy_ref
                        
                    if not ramp: 
                        bcl[1] = bcl[1].replace('r', '')
                    
                    v = 1 
                    if 'TiO2' in cl[0] or 'Mn2O4' in cl[0]:
                        v = 100
                        # print (v, cl[0], bcl[0])
                    # else:
                    #     v = 1


                    res_loop(bcl[0], bcl[1],               v, choose_outcar = -1,readfiles = readfiles)
                    final_list, _ = res_loop(cl[0], cl[1], v, choose_outcar = -1,readfiles = readfiles, analys_type = 'redox_pot', b_id = (bcl[0], bcl[1], v), energy_ref = energy_ref)
                        # print final_list


            elif len(cl) == 5:
                if not ramp: cl[1] = cl[1].replace('r', ''); cl[3] = cl[3].replace('r', '')
                
                if run: 
                    add_loop(cl[0], cl[1], 1, up = 'up1', input_geo_format = input_geo_format, it_folder = cl[4])
                    id2 = cl[2]+cl[3]
                    if id2 not in DSs: 
                        add_loop(cl[2], cl[3], 1, up = 'up1', input_geo_format = input_geo_format, it_folder = cl[4])
                    DSs.append(id2)

                else:
                    final_list, _ = res_loop(cl[0], cl[1], 1, readfiles = readfiles, analys_type = 'redox_pot', b_id = (cl[2], cl[3], 1), energy_ref = energy_ref)

            # print 'hello'
            # print (final_list)
            if not run and  final_list:
                results.append(final_list)

    if not run:
        ''
        # print (results)
        df = pd.DataFrame(data = results, columns=['is', 'redox_pot', 'kspacing', 'time', 'mdstep', 'ecut', 'niter', 'set_is'])
        # print df.set_index(['is', 'mdstep', 'set_is']).sortlevel(0).round(3)



        if 1: # compare souces
            df = pd.DataFrame(data = results, columns=['is', 'redox_pot', 'id_is', 'id_ds'])[['is', 'redox_pot']]
            df['source'] = 'GGA+U_ramp'
            df['type'] = 'GGA+U'
            df['is'] = [ d.split('.')[0] for d in df['is']]
            # print df
            dfs = pd.read_csv(r'database/literature.csv')[['is', 'redox_pot', 'source', 'type']]
            dfs.drop(0, inplace=True)
            dfs.drop(1, inplace=True)
            # print (dfs)
            # sys.exit()

            df = pd.concat([df, dfs])
            # df = df.set_index(['is', 'redox_pot']).sortlevel(0)
            df = df.reset_index(drop=True)
            # print (df)

            df['is'] = pd.Series([re.sub("([0-9])", "$_\\1$", s) for s in df['is'] ] ) # #Make low indexes

            
            wc2= 'redox_pot'
            df = df.apply(lambda x: pd.to_numeric(x, errors='ignore'))
            df[[wc2]] = df[[wc2]].apply(lambda x: pd.Series.round(x, 1))


            if 1: #compare with other dft
                df_piv = df.pivot(index='is', columns='source', values='redox_pot')
                # print(df_piv)
                df_piv = df_piv.sort_values(['azh'])
                df_piv[['GGA+U_ramp', 'azh', 'Shishkin2016' ]].plot(kind = 'bar',rot = 70, ylim = (0,6), 
                    colormap = 'RdYlBu')
            else: #compare with experiment
                df = df.dropna()
                # print (df)
                df_piv = df.pivot(index='is', columns='type', values='redox_pot')[['GGA+U', 'exp']]
                df_piv = df_piv.sort_values(['exp'])
                df_piv.plot(kind = 'bar',rot = 70, ylim = (0,6), 
                    colormap = 'RdYlBu')



            

            plt.xlabel('')
            plt.ylabel('Intercalation voltage (V)')


            # plt.legend(loc='best') 
            plt.legend(loc='best', ncol =3, handlelength = 0.25, 
                handleheight = 2, frameon=False, fontsize=20)

        if 1:
            # mpl.rcParams.update({'font.size': 14})
            # mpl.rcParams.update({'font.size': 18})
            # mpl.rc('legend', fontsize= 10) 


            plt.tight_layout()




            if image_name:
                if not os.path.exists(header.path_to_images):
                    os.makedirs(header.path_to_images)
                path2im = header.path_to_images+str(image_name)+'.'+image_format
                printlog( "Saving image ...", path2im)
                # plt.savefig(path2im, dpi = dpi, format=image_format)
                plt.savefig('figs/'+image_name+'.png', dpi = 300)
                plt.savefig('figs/'+image_name+'.eps', dpi = 300)
            else:
                plt.show()






def potential_barriers():
    """
    Creates figures and tables for report
    """
    calc = header.calc
    readfiles = False

    res_loop('LiFePO4.nebLi1v1',[u'1ulong'], verlist = [1], up = '', show = 'm', push2archive = True,
    readfiles = readfiles,
    description_for_archive = 'Diffusion path for Li ')    #mep

    res_loop('LiFePO4.nebLi1v1',   '1u', [1], up = '', show = 'm', push2archive = True,
    readfiles = readfiles,
    description_for_archive = 'Diffusion path for Li '
     )   #unexpectedly slightly lower energy for saddle point

    res_loop('LiFePO4.n7Li1v1',[u'1ulong'], verlist = [1], up = 'up1', show ='m', push2archive = True,
    readfiles = readfiles,
    description_for_archive = 'Diffusion path for Li '
     )


    res_loop('FePO4.pnma2.n3i1e2Lims', '1u', verlist = [1], up = 'up1', show = 'm', push2archive = True,
    readfiles = readfiles,
    description_for_archive = 'Diffusion path for Li without u-ramping '
     )
    res_loop('FePO4.pnma2v1.n3i1e2Lims', '1ur', verlist = [1], up = 'up1', show = 'm',push2archive = True,
    readfiles = readfiles,
    description_for_archive = 'Diffusion path for Li, u-ramping du =  1.3 V'
     )
    res_loop('FePO4.pnma2v1.n7i1e2Lims',[u'1ur30'], verlist = [1], show = 'm', push2archive = True,
    readfiles = readfiles,
    description_for_archive = 'Diffusion path for Li, u-ramping du =  0.13 V '
     )



    res_loop('NaFePO4.pnma.n3Na1v1mp',[u'1ulong'], verlist = [1], up = 'up1', show = 'm', push2archive = True,
    readfiles = readfiles,
    description_for_archive = 'Diffusion path for Na '
     )

    res_loop('FePO4.pnma2.n3i1e2Na',[u'1u'], verlist = [1], up = 'up1', show = 'm', push2archive = True,
    readfiles = readfiles,
    description_for_archive = 'Diffusion path for Na '
     )

    res_loop('KFePO4.pnma.suv100.n3K1v1ms',[u'1u'], verlist = [1],show = '', push2archive = True, 
    readfiles = readfiles,
    description_for_archive = 'Diffusion path for K ')    #mep


    res_loop('LiVPO4F.nebLi1v1','1ur', verlist = [1], up = '', show = '', push2archive = True, 
    readfiles = readfiles,
    description_for_archive = 'Diffusion path for Li ')    #mep

    res_loop('LiMnPO4v1.n3Li1v1ms',[u'1u'], verlist = [1], up = 'up1', show = 'm', push2archive = True, 
    readfiles = readfiles,
    description_for_archive = 'Diffusion path for Li ') 


    mpl.rc('axes', titlesize= 12) 

    mpl.rcParams.update({'font.size': 14})

    res_loop('KFePO4.pnma.su','4uiWF', range(1,8)+[100], up = '1', show = 'fo', analys_type = 'fit_a', 
        push2archive = True,
        readfiles = readfiles,
        description_for_archive = 'Cell shape relaxation and search of minimum')    #mep



def find_matched_calculations(it):
    """
    Look for it in struct_des and find similar calcs.
    Useful to find corresponding intercalated structure
    """
    for key in header.struct_des:
        if it in key:
            print (key)
    sys.exit()




def calc_barriers(mode = '', del_ion = '', new_ion = '', func = 'gga+u', show_fit = False, up = None, upA = None, upB = None, upC = None,
    param_dic = None, run_sc = 1, run_neb = 1, up_add_loop = 'up2', 
    up_res = 'up1', add_loop_dic = None, fitplot_args = None,
    nise = 1,
    choose_outcar_global = None,
    cathodes = None, style_dic = None, readfiles = 1 ):
    """
    Simplifies calculation of intercalation voltages
    and barriers by performing step by step calculations - 
    remove cations or replace cations,
    create supercell, optimize supercell, calculate migration barriers

    mode (str) -
        - 'normal'
            param_dic['start_pos'] (int) - number of starting non-eqiv position
            param_dic['end_pos'] (int)   - number of final non-eqiv position

        - 'replace' - creates additional structures by replacing cations

            - *new_ion* (str) - choose new ion

        - 'make_ds'  - creates deintercalated structures by removing all cations and then inserting in analogy to

            - *del_ion* - ions to be deleted either by replacement or removing, use 'Li Na' if both should be removed

            param_dic['neb_search_voids'] - if exists then in 'make_ds' mode new voids are searched as final positions

        - func - functional either gga or gga+u, if gga is used the energies are taken for u=0, however
        take in mind that optimization was made for 'gga+u'

    nise - if 1 than neb_ise is use, if 0 than main_ise is used. Variable is useful for two-step barrier determination with the second step with fixed occ

    param_dic - dictionary with parameters
        'id' - starting id
        'ds' - name of ds structure
        'itfolder' - folder for calculation
        'main_set' - set used for supercells and neb
        'neb_set' - set for neb used instead of main_set
        'readfiles' - to res_loop
        -del_pos - specific non-eqv position of del_ion to be removed starting from 0
        
        - meps - list of tuples (start_pos, end_pos)
        - meps2 - list of tuples (i_atom_to_move, end_pos)
        - ortho - global target sizes of supercell
        - ortho.mode_id - specific target sizes of supercell
        - old_behaviour (bool) - if True then first optimization of cell is not done for 'normal', second optimization for all modes;
            if False - only first optimization for all cells
        'old.'+mode_id
        ! war

    run_sc - run supercell construction part
    run_neb - run neb construction part 

    style_dic - passed to plot_mep()

    choose_outcar_global - allows to override param_dic values, choose outcar to read from sequence set or u-ramp loop

    cathodes - list of lists - another way to configure the calculation - was initially used

    update




    TODO
    curver - is 1!!! check if it works correctly for other cases


    return


    """
    if not param_dic:
        param_dic = {}

    if not add_loop_dic:
        add_loop_dic = {}


    calc = header.calc
    struct_des = header.struct_des
    varset = header.varset


    def calc_added(cl, attr, state):
        """
        Check whether *cl* has attribute *attr* (tuple) and attr[3] is equal to *state*
        """
        added = False
        hasat = False
        if hasattr(cl, attr):
            hasat = True

            if getattr(cl, attr)[3] == state:
                added = True
        
        return added and hasat

    def calc_added2(cl, attr, state):
        """
        Check whether *cl* has attribute *attr* (tuple) and attr[3] is equal to *state*
        """
        # print(attr, state)
        # print(cl.neb_id)
        # sys.exit()
        added = False
        hasat = False
        if hasattr(cl, attr):
            hasat = True

            if state in getattr(cl, attr):
                added = True
        else:
            setattr(cl, attr, {})
        return added and hasat



    def make_dummy_calc_obj(calc_id):
        """
        if id not in calc saves dummy Calculation() to calc database
        """
        if calc_id:
            if calc_id in calc:
                cl = calc[calc_id]
            else:
                cl = Calculation() #just empty calc to store state iformation
                calc[calc_id] = cl
            return cl


    def cat_name_mod(init_name):
        #'cathode_name'
        if 'replace' in mode:
            name = init_name.replace(del_ion, new_ion)
        elif 'make_ds' in mode:
            name = init_name.replace(del_ion, '')
        elif 'normal' in mode:
            name = init_name
        return name


    def replace_cations(ion_to_replace = None, replace_by =None, datalist = None, ver = None):
        """Step P: Create cells with Na and K from Li cells and optimize them first for small cells"""
        cat = datalist
        LiX    = cat[0]
        curise = cat[1]
        curver = ver

        curit  = LiX.replace(ion_to_replace, replace_by)
        curfol = cat[4].replace(ion_to_replace, replace_by)

        itP = curit+'.ir' #inherit by replacement

        if ion_to_replace in LiX and ion_to_replace != replace_by:

            if itP not in struct_des:
                inherit_icalc('replace_atoms', itP, curver, 
                (LiX, curise, curver), calc, atom_new = replace_by, atom_to_replace = ion_to_replace, it_folder = curfol, use_init = 1) 
        else:
            printlog('Error!', ion_to_replace, 'was not found in name', LiX, '; Use correct names')
            itP = None    
        return  itP


    def remove_cations(ions = None, cat = None, ver = None, update = None): 
        """Step S: remove cations and create DS structures from IS
            ions (list)
        """
        idA0 = (cat[0], cat[1], ver)
        dic = cat[7]


        if 'del_pos' in dic:
            itS = cat[2]+'.id'+str(dic['del_pos']) #inherit delete
            del_pos = dic['del_pos']
        else:
            itS = cat[2]+'.id'
            del_pos = None

        ions_present = False
        for ion in ions: # all ions should be in name
            if ion in cat[0]:
                ions_present = True
            else:
                ions_present = False
        # print(ions_present)

        if ions_present:
            # print(update)

            if update or itS not in struct_des:
                # print(update)
                inherit_icalc('remove_atoms', itS, ver, idA0, calc, 
                    atoms_to_remove = ions, del_pos = del_pos, it_folder = cat[4], use_init = 1) 
        
            return itS



    def optimize_cell(id_base = None, datalist = None, scale_region  = None, update = None, irun = ''):
        """
        Create and read results
        irun (str) - number of run 
        """
        dic = cat[7]
        # print(id_base)
        # sys.exit()

        if id_base and id_base[0]:
            
            if 'scaling_set' in dic:
                ise_new = dic['scaling_set']
            else:
                try:
                    ise_new = dic['scaling_set'+irun+'.'+mode_id] #
                except:
                    ise_new = dic['scaling_set.'+mode_id] #

            if 'save' in dic:
                savefile = dic['save']
            else:
                savefile = 'o'

            if dic['scale_region.'+mode_id]:
                printlog('optimize_cell(): scale_region changed from', scale_region, 'to',  dic['scale_region.'+mode_id], imp ='y')

                scale_region = dic['scale_region.'+mode_id]


            curfol = struct_des[id_base[0]].sfolder

            id_res = (id_base[0]+'.su', ise_new, 100)


            if update or id_res not in calc:
                printlog('Scale region is ', scale_region, imp = 'y')
                add_loop(*id_base, up = up_add_loop, ise_new = ise_new, savefile = savefile, calc_method = 'uniform_scale', scale_region = scale_region, 
                inherit_option = 'inherit_xred', it_folder = curfol+'/scaled/', **add_loop_dic)                    
            
            else:
                # print(dic['scale_outcar.'+mode_id], up_res)
                # print(up_res)
                # sys.exit()
                res_loop(*id_res, up = up_res, readfiles= readfiles, choose_outcar = dic['scale_outcar.'+mode_id], show = 'e', check_job = 0)
                if show_fit:
                    res_loop(*id_res[0:2],list(range(0+1,0+8))+[100], up = up_res, readfiles= readfiles, choose_outcar = dic['scale_outcar.'+mode_id], analys_type = 'fit_a', show = 'fitfo', check_job = 0)
                # sys.exit()
                if '2'  in calc[id_res].state or '5' in calc[id_res].state:
                    ''
                    # del calc[id_res]

                elif '4' in calc[id_res].state:
                    # print(np.real(calc[id_res].B))
                    try:
                        curres.update({'B':float(calc[id_res].B)})
                    except:
                        printlog('No Bulk modulus')
                    
                        
                    curres['vol'] = calc[id_res].end.vol

                    curres['id'] = id_res

                    if '2' in irun:
                        ''
                        if add_to_database:
                            add_to_database(calc[id_res])



                    return id_res




    def make_supercell(base_id, cat, update):
        """Step A: Create supercells and calculate their total energy"""
        dic = cat[7]
        pd = dic
        ise_new = cat[6]

        if 'ortho' in pd:
            ortho = pd['ortho']
            printlog('ortho found in pd:', ortho)

        elif 'ortho.'+mode_id in pd:
            ortho = pd['ortho.'+mode_id]
            printlog('ortho mode_id found in pd:', ortho)
        else:
            ortho = [10,10,10]
            printlog('Default ortho is used:', ortho)
        nonlocal support_dict_key
        
        support_dict_key = 'support_'+name_mod_supercell(ortho)
        

        sc_unique_id = (name_mod_supercell(ortho), ise_new)


        if support_dict_key not in calc:
            calc[support_dict_key] = {} #create supportive dictionary, which depends on ortho parameter
        



        if base_id:
            # print(base_id)
            clA0 = db[base_id]

            mul_matrix = None
            ngkpt      = None




            # print(clA0.inh_id, sc_unique_id)
            # sys.exit()
            # if update or not calc_added(clA0, 'inh_id', 1):
            if update or not calc_added2(clA0, 'inh_id', sc_unique_id):
                

                if 'make_ds' in mode: #for DS structures use the same parameters as for IS structures
                    if 1:
                        printlog ('Cell obtained by removing atoms', imp = 'y')
                        printlog ('I use mul_matrix and ngkpt from intercalated structure', imp = 'y')
                        try:
                            mul_matrix = calc[support_dict_key][cat[0], 'mul_matrix']
                            ngkpt_dict      = calc[support_dict_key][cat[0], 'ngkpt_dict']

                        except KeyError:
                            printlog('Error! You are in IS->DS inheritance mode, please first make IS calculation to fill in supportive dictionary with mul_matrix and ngkpt')
                        ks = varset[ise_new].kspacing
                        ngkpt = ngkpt_dict[ks]
                        # printlog ('I use mul_matrix and ngkpt from intercalated structure', imp = 'y')


                sfolder = struct_des[base_id[0]].sfolder
                # print (sfolder)
                # sys.exit()
                itA = add_loop(*base_id, ise_new = ise_new, up = up_add_loop, 
                it_folder = sfolder+'/super/', inherit_option = 'supercell', ortho = ortho,
                mul_matrix = mul_matrix, ngkpt = ngkpt , **add_loop_dic)
                
                # print([itA, ise_new, base_id[2], support_dict_key, 'exist'])
                # clA0.inh_id = [itA, ise_new, base_id[2], support_dict_key, 'exist']
                if not hasattr(clA0, 'inh_id') or type(clA0.inh_id) == list:# temporary for compat
                    clA0.inh_id = {} 

                clA0.inh_id[sc_unique_id] = (itA, ise_new, base_id[2])  # additional unique id and state, TODO: maybe move to add_neb 

                # print(clA0.inh_id)
                # print(db[('Cu.su', '1lo', 100)].inh_id)
                # header.db = db
                if 'normal' in mode:
                    calc[support_dict_key][cat[0], 'mul_matrix'] = struct_des[itA].mul_matrix               #save for intercalated structure
                    calc[support_dict_key][cat[0], 'ngkpt_dict'] = struct_des[itA].ngkpt_dict_for_kspacings #save for intercalated structure



            else:
                idA = clA0.inh_id[sc_unique_id]
                
                if func == 'gga':
                    choose_outcar = 1
                else: #gga+u with additional control 
                    choose_outcar = dic['neb_outcar.'+mode_id]
                # print(dic['neb_outcar.'+mode_id])
                # print('suprecell: chosen outcar:', choose_outcar)
                # print(up_res)
                res_loop(*idA, up = up_res, readfiles = readfiles, choose_outcar = choose_outcar, show = 'm', check_job = 0)
                # res_loop(*idA, choose_outcar = 3, show = 'fo')
                # res_loop(*idA, choose_outcar = 4, show = 'fo')
                curres['id_sc'] = idA
                if '4' in calc[idA].state:
                    


                    return idA


    def make_neb(base_id, atom_to_insert, cat, update):
        """Step C: make neb
        base_id - calc id - where to study migration
        cat - special list of parameters, at cat[7] contains dict with paramerters
            dict
                'start_pos' - i_void_start
                'end_pos' - i_void_final
                'meps' - list of tuples ('start_pos', 'end_pos')

        TODO:
        1. Transform cat to pure dicts?
        2. Make names consistent, like 'start_pos' and i_void_start is the same!

        """
        pd = cat[7]#param_dic

        # print(atom_to_insert)
        # sys.exit()


        if 'neb_set' in pd and pd['neb_set']:
            ise_new = pd['neb_set']
        else:
            ise_new = cat[6]

        if base_id:
            idB = base_id
            clB = calc[idB]
            # print(pd)
            default = {'i_atom_to_move':'', 'start_pos':1, 'end_pos':1}
            for key in default:
                if key not in pd:
                    pd[key] = default[key]


            printlog('Start:', pd['start_pos'], '; End:', pd['end_pos'])

            sup_key = (cat[0], pd['start_pos'], pd['end_pos'], pd['i_atom_to_move']) # support key - should be unique for current path, used to transfer xr from IS to DS 

            # if 'replace' in mode:
            #     sup_key = (cat[0], pd['start_pos'], pd['end_pos'], pd['i_atom_to_move'])


            neb_unique_id = (pd['start_pos'], pd['end_pos'], pd['i_atom_to_move'], atom_to_insert, ise_new)
            if 'rep_moving_atom' in pd:
                # sys.exit()
                neb_unique_id = (pd['start_pos'], pd['end_pos'], pd['i_atom_to_move'], atom_to_insert, ise_new, pd['rep_moving_atom'])



            printlog('neb_unique_id = ', neb_unique_id)

            if 'images' in pd:
                images = pd['images']
            else:
                images = 3
                printlog('Number of images set to', images)

            ok = 'occmatrix_id' # general

            if ok not in pd:
                ok = 'occmatrix_id.'+mode_id





            if ok in pd:
                cl = db[pd[ok]]
                occfile = write_occmatrix(cl.occ_matrices, cl.dir)

                add_loop_dic['params'] = {'occmatrix':occfile}


                # up_res = 'up1'

            # sys.exit()
            if update or not calc_added2(clB, 'neb_id', neb_unique_id):

                other_param = {}
                if 'r_impurity' in pd:
                    other_param['r_impurity'] = pd['r_impurity']

                if 'i_atom_to_move' in pd:
                    other_param['i_atom_to_move'] = pd['i_atom_to_move']

                if 'atom_to_move' in pd:
                    other_param['atom_to_move'] = pd['atom_to_move']
                    # print(other_param['atom_to_move'])
                if 'end_pos_types_z' in pd:
                    other_param['end_pos_types_z'] = pd['end_pos_types_z']

                if 'rep_moving_atom' in pd:
                    other_param['rep_moving_atom'] =  pd['rep_moving_atom']
                       


                # print('sdfsaf')

                if 'normal' in mode or 'replace' in mode:
                
                    if 'neb_search_voids' in pd and pd['neb_search_voids'] == 1:
                        search_type = 'existing_voids'
                    else:
                        search_type = 'vacancy_creation'

                    # print('make_neb:', add_loop_dic)
                    it = add_neb(clB, up = up_add_loop, ise_new = ise_new, images = images, 
                        i_void_start = pd['start_pos'], i_void_final = pd['end_pos'], 
                        atom_to_insert = atom_to_insert,
                        search_type = search_type, add_loop_dic = add_loop_dic, old_behaviour = old_behaviour, **other_param)                
                
                elif 'make_ds' in mode:

                    # printlog('made_ds mode')
                    if 'neb_search_voids' in pd and pd['neb_search_voids'] == 1:
                        it = add_neb(clB, up = up_add_loop, ise_new = ise_new, images = images, 
                            i_void_start = pd['start_pos'], i_void_final = pd['end_pos'],
                            atom_to_insert = atom_to_insert, search_type = 'existing_voids', 
                            add_loop_dic = add_loop_dic, old_behaviour = old_behaviour, **other_param)                          
                    else:
                        # use the same positions as was used in normal
                        xr_start = calc[support_dict_key][sup_key, 'xr_m_ion_start']
                        xr_final = calc[support_dict_key][sup_key, 'xr_m_ion_final']
                        printlog('Using start and final positions from', support_dict_key, sup_key, xr_start, xr_final, imp = 'Y')
                        # sys.exit()
                        it = add_neb(clB, up = up_add_loop, ise_new = ise_new, images = images,  
                            xr_start = xr_start,
                            xr_final = xr_final,
                            i_void_start = pd['start_pos'], i_void_final = pd['end_pos'], #just for name
                            atom_to_insert = atom_to_insert, add_loop_dic = add_loop_dic,
                            old_behaviour = old_behaviour)                

                id_n = (it, ise_new, 1)
                
                if 'map' not in calc: # attempt to create map of calculations of graph of calculations, not working yet
                    calc['map'] = {}
                    mp = calc['map']
                else:
                    mp = calc['map']
                # print(mp)
                # sys.exit()
                if idB not in mp:
                    mp[idB] = []
                
                try:
                    id_nold = clB.neb_id[0:3]
                    
                    if id_nold not in mp[idB]:
                        mp[idB].append(id_nold)
                except:
                    pass


                if id_n not in mp[idB]:
                    mp[idB].append(id_n)


                if not hasattr(clB, 'neb_id') or type(clB.neb_id) == list:# temporary for compat
                    clB.neb_id = {} 
                # print(neb_unique_id)
                clB.neb_id[neb_unique_id] = (it, ise_new, 1)  # additional unique id and state, TODO: maybe move to add_neb 
                # print(clB.id, clB.neb_id)

                if 'normal' in mode:
                    calc[support_dict_key][sup_key, 'xr_m_ion_start'] = struct_des[it].xr_m_ion_start # xred coordinate of migrating ion in starting position
                    calc[support_dict_key][sup_key, 'xr_m_ion_final'] = struct_des[it].xr_m_ion_final # xred coordinate of migrating ion in final    position
                    # print (struct_des[it].x_m_ion_start, struct_des[it].x_m_ion_final)
                    # sys.exit()




            else:
                idC =  clB.neb_id[neb_unique_id]



                # res_loop(*idC[0:2],[1,2], show = 'me', readfiles = 1, choose_outcar = 3)
                # res_loop(*idC[0:2],[4], up = 'up1',show = 'mag', readfiles = 1, choose_outcar = 3, analys_type = 'ne')
                # res_loop(*idC[0:2],[4], up = 'up1',show = 'mag', readfiles = 1, choose_outcar = None, analys_type = 'ne')
                # print (dic)
                # print (dic['neb_outcar.'+mode_id])
                if choose_outcar_global:
                    choose_outcar = choose_outcar_global

                elif func == 'gga':
                    choose_outcar = 1
                else: #gga+u with additional control 
                    okey = 'neb_outcar.'+mode_id
                    printlog('choosing outcar from dic ', okey)
                    choose_outcar = dic[okey]
                printlog('chosen outcar = ', choose_outcar)


                if neb_fit:
                    ''
                    printlog('up key of res_loop = ', up_res)
                    # printlog('choose_outcar = ', choose_outcar)

                    # sys.exit()

                    res, _ = res_loop(idC[0],ise_new,range(1, images+3), up = up_res, show = 'fomep', readfiles = readfiles, 
                        choose_outcar = choose_outcar, 
                        analys_type = 'neb', fitplot_args = fitplot_args, 
                        check_job = add_loop_dic['check_job'], style_dic = style_dic, params = pd)
                else:
                    res_loop(*idC[0],ise_new, [1,2], show = 'me', up = up_res, readfiles = readfiles, choose_outcar = choose_outcar, check_job = add_loop_dic['check_job'])



                # print (res)

                # if '5' in calc[idC].state: #some error
                #     clB.neb_id[3] = None
                
                if '4' in calc[idC].state:
                    # try:
                        # curres['dAO_change'] = res['dAO_change']
                    curres['sts'] = res['sts']
                    # except:
                    #     pass
                    

                    try:
                        curres.update({'barrier':calc[idC].barrier})
                    
                    except:
                        curres.update({'barrier':0})
                    
                    curres['dEm1'] = res['dEm1'] #difference of energies between middle and first image of NEB

                    # print(curres['barrier'], curres['dEm1'])



                    curres['id_start'] = (idC[0], idC[1],1)
                    curres['id_end']   = (idC[0], idC[1],2)


                    curres['atom_pos']     = res['atom_pos']
                    curres['mep_energies'] = res['mep_energies']

                    # dn = calc[idC].end.natom - calc[idB].end.natom
                    # v_int1 = calc[idC[0], idC[1],1].energy_sigma0 - calc[idB].energy_sigma0 - dn * curres['Eref']
                    # v_int2 = calc[idC[0], idC[1],2].energy_sigma0 - calc[idB].energy_sigma0 - dn * curres['Eref']
                    printlog('Calculating atom-wise intercalation potentials', imp = 'Y', end = '\n')
                    v_int1 = calc_redox(calc[idC[0], idC[1],1], calc[idB])['redox_pot']
                    v_int2 = calc_redox(calc[idC[0], idC[1],2], calc[idB])['redox_pot']

                    
                    # print('Redox pot:', dn, v_int1, v_int2, curres['Eref'])


                    curres.update({'V2':v_int1}) # the batterry is fully discharged
                    return




    if 1: # prepare parameters
        main_set    = '1urs'   
        results = [] # result list, for each structure we have dict
        scale_regions = {}

        neb_fit = 1 # if to make neb fitting; the barries is saved to cl1.barrier
        run_scale = 1
        add_to_database = 0

        curver = 1
        local_state = None
        
        support_dict_key = ''

        pd = param_dic





        if 'check_job' not in add_loop_dic:
            add_loop_dic['check_job'] = 1


        # print(pd)
        if 'scaling_set' in pd:
            scaling_set = pd['scaling_set']
            # print(scaling_set)
        else:
            scaling_set = '4uis'
        # print(scaling_set)
        # sys.exit()

        mode_id = mode+'.'+del_ion+'.'+new_ion
        dic = {
        'neb_outcar.' + mode_id:None, #by default the last one is used
        'scale_outcar.' + mode_id:None, #by default the last one is used
        'scaling_set.'+ mode_id:scaling_set, #
        'skip.'+mode_id:False,
        'scale_region.'+mode_id:None,
        }# additional specific control  
        # 'neb_outcar' - (int) - choose_outcar in res_loop of make_neb

        if 'scale_region' in pd:
            dic['scale_region.'+mode_id] = pd['scale_region']

        scale_regions['optimal'] = (-4,4)
        scale_regions['Li'] = (-5,3)
        scale_regions['Na'] = (-1,7)
        scale_regions['K']  = (4,16)
        scale_regions['Rb']  = (-4,4)
        scale_regions['Li_removal']  = (-5,3)
        scale_regions['Na_removal']  = (-5,3)
        scale_regions['K_removal']   = (-8,0)
        scale_regions['Rb_removal']   = (-6,2)

        scale_regions[new_ion] = (-4,4)
        
        # print(new_ion)
        # sys.exit()

        
        if up:
            update  = up # initial update for removing and replacing atoms
        else:
            update  = 0
        
        if upA:
            updateA = upA
        else:
            updateA = 0 # make supecell 
        
        if upB:
            updateB = upB
        else:    
            updateB = 0 # volume optimization, su
        
        if upC:
            updateC = upC
        else:
            updateC = 0 # NEB
        
        if not up_res:
            up_res = 'up1' # update files

        if 'readfiles' in pd:
            readfiles = pd['readfiles']
        else:
            readfiles = readfiles
        
        if not show_fit:
            show_fit  = 0 #control if to show fitting of supercell sizes
        

    ise = '8u' # old
    if pd:
        if 'ds' not in pd:
            pd['ds'] = 'None'
        cathodes = [ [pd['id'][0], pd['id'][1], pd['ds'], '', pd['itfolder'], '', pd['main_set'], md(dic, pd) ], ]
        curver = pd['id'][2]

    else:
        if nise:
            neb_set = nise #'1uos'
        else:
            neb_set = None
        
        if cathodes:
            for c in cathodes:
                c[-1] = md(dic, c[-1])
        else:
        # cathodes == None:

            cathodes = [  #cat[3] and cat[5] eventually was not used                                  
               # ['LiCoO2', ise, 'CoO2', ise, 'LiCoO2'       ,'', main_set, merge_dics(dic,{'scaling_set.make_ds.Li.Li':'4uris',   }  ) ], old
               # ['LiCoO2', ise, 'CoO2', ise, 'LiCoO2'       ,'', '1u', merge_dics(dic,{'scaling_set.normal.Li.Li':'1u',   }  ) ],
               # ['LiCoO2', ise, 'CoO2', ise, 'LiCoO2'       ,'', '1uafm', merge_dics(dic,{ })  ], 
               # ['LiTiO2', ise, 'TiO2', ise, 'LiTiO2',         '',  main_set, merge_dics(dic,{'scaling_set.normal.Li.Li':'4uisns', 'scaling_set.replace.Li.Na':'4uisns', 'scaling_set2.make_ds.Li.Li':'4uisns'}) ], # only for IS Li Na
               # ['LiNiO2', ise, 'NiO2', ise, 'LiNiO2'       ,'', main_set, merge_dics(dic,{'scale_outcar.make_ds.Li.Li':None, 'neb_outcar.replace.Li.Na':3,'neb_outcar.replace.Li.K':3, 'scaling_set.normal.Li.Li':'4uisns', 'scaling_set2.make_ds.Li.Li':'4uris', 'scaling_set1.make_ds.Li.Li':'4uisC0W1', 'save':'oc'}) ] , # for Li use larger symprec, or isym = 0
               # ['LiTiS2', ise, 'TiS2', ise, 'LiTiS2'       ,'', main_set, md(dic, {'neb_outcar.make_ds.Li.Li':None}) ] ,
               # ['LiMn2O4',   '', 'Mn2O4',   '', 'LiMn2O4' , '', main_set, md(dic, {'neb_outcar.normal.Li.Li':3} )  ] ,
               # ['LiVP2O7', ise, 'VP2O7', ise, 'LiVP2O7'    ,'', main_set, md(dic, {'neb_outcar.normal.Li.Li':None, 'scale_region.replace.Li.K':(0,8) } )   ] ,

               # ['Li2CoPO4F.pnma',   '1u', 'CoPO4F',   '', 'Li2CoPO4F/pnma' , '', main_set, md(dic, {} )  ],
               # ['NaMnAsO4',  ise, 'MnAsO4', ise, 'NaMnAsO4' ,'', main_set, md(dic, {'neb_outcar.normal.Li.Li':None} ) ] ,
               # ['Na2FeVF7',   '', 'FeVF7', '', 'Na2FeVF7' , '', main_set, dic ] ,
               # ['KFeSO4F',   '', 'FeSO4F',  '', 'KFeSO4F' , '', '1urS', md(dic,{'scaling_set.normal.Li.Li':'4uSi', 'scale_outcar.normal.Li.Li':1, 'scaling_set.make_ds.K.K':'4uSi'  }) ] ,
               # ['NaLiCoPO4F' ,'', 'CoPO4F','','NaLiCoPO4F' , '', main_set, md(dic,{'scaling_set.make_ds.Li Na.Li':'4uris',   })  ] ,

               # ['LiNiO2.r3m',   '1u', 'NiO2.r3m',  '', 'LiNiO2/r3m' , '', '1urAl', md(dic,{'neb_outcar.normal.Li.Li':1, 'scaling_set.normal.Li.Li':'4uisC0W1', 'scaling_set.make_ds.Li.Li':'4uiAl',    }) ] , # for IS main_set is '1urns'


               # for paper exchange, from 06.09.2017 moved to file paper7
               # ['LiMn2O4.Fd3m.111',   '1u', 'Mn2O4.Fd3m.111', '', 'LiMn2O4/Fd3m' , '', main_set, md(dic,{'neb_set':neb_set, 'old.normal.Li.Li':1, 'neb_outcar.normal.Li.Li':None,'neb_outcar.replace.Li.Na':3, 'occmatrix_id.replace.Li.Na':'NaMn2O4.Fd3m.111.ir.su.s10v100.n3Na1v1ms.1urs.4'}) ] ,
               # ['LiFePO4', ise, 'FePO4', ise, 'LiFePO4'    ,'', main_set, md(dic,{'neb_set':neb_set, 'neb_outcar.normal.Li.Li':3, 'neb_outcar.replace.Li.K':3, 'old.normal.Li.Li':1,'old.make_ds.Li.Li':0,  'occmatrix_id.normal.Li.Li':'LiFePO4.s10.suv100.n3Li1v1ms.1urs.3', 'occmatrix_id.make_ds.Li.Li':'FePO4.id.su.s10v100.n3Li1v1ms.1urs.4' }) ] ,
               # ['LiMnPO4', ise, 'MnPO4', ise, 'LiMnPO4'    ,'', main_set, md(dic, {'neb_outcar.make_ds.Li.Li':None, 'old.normal.Li.Li':1, 'old.make_ds.Li.Li':1})  ] ,
               # ['LiVPO4F', ise, 'VPO4F', ise, 'LiVPO4F'    ,'',  main_set, md(dic,{'ortho.replace.Li.Na':[10,10,13],'ortho.replace.Li.K':[10,10,13], 'old.normal.Li.Li':1, 'scale_outcar.replace.Li.Na':None,  'neb_outcar.normal.Li.Li':None, 'neb_outcar.make_ds.Li.Li':None,'scaling_set.normal.Li.Li':'4uris', 'scaling_set.replace.Li.K':'4uisC0W1', 'scaling_set.replace.Li.Na':'4uris', 'scale_region.replace.Li.K':(-1,7),   'skip.replace.Li.K':0, 'skip.replace.Li.Na':0}) ]  ,# only for IS
               # ['KVPO4F.Pna21',   '1u', 'VPO4F.Pna21', '', 'KVPO4F/Pna21',  '', '1urs', md(dic,{'neb_outcar.normal.Li.Li':None, 'old.normal.K.K':1, 'old.normal..':1, 'images':3}) ],
               # ['Li2FePO4F.pnma',   '1u', 'LiFePO4F.pnma',   '', 'Li2FePO4F/pnma' , '', '1u', md(dic, {'ortho.replace.Li.Na':[13,13,13],'ortho.replace.Li.K':[13,13,13],'meps':[(1,1), (2,1), (1,3), (2,3), (3,2), (3,3)], 'del_pos':1} )  ], #(1,1), 
               # ['Li2FePO4F.pnma',   '1u', 'LiFePO4F.pnma',   '', 'Li2FePO4F/pnma' , '', '1u', md(dic, {'neb_set':'1urs', 'ortho.replace.Li.Na':[13,13,13],'ortho.replace.Li.K':[13,13,13],'meps':[(3,2), (3,5)], 'del_pos':1} )  ], #(1,1), 
               # ['Na2FePO4F', '', 'NaFePO4F', '', 'Na2FePO4F' , '', '1u', md(dic, {'neb_set':nise, 'meps':[(1,1), (1,2), (1,3), (2,1), (2,2), (2,3)], 'scaling_set':'4uis', 'del_pos':2})  ] , #neb_set: 1urs, 1u
               # ['Na2FePO4F', '', 'NaFePO4F', '', 'Na2FePO4F' , '', '1u', md(dic, {'neb_set':nise, 'meps':[(2,2),  ], 'scaling_set':'4uis', 'del_pos':2})  ] , #neb_set: 1urs, 1u, DS

               # for checks
               # ['Li2FePO4F.pnma',   '1u', 'LiFePO4F.pnma',   '', 'Li2FePO4F/pnma' , '', '1u', md(dic, {'ortho.replace.Li.Na':[13,13,13],'ortho.replace.Li.K':[13,13,13],'meps':[(1,3)], 'del_pos':1} )  ], #(1,1), 




               # other try
               # ['LiFePO4', ise, 'FePO4', ise, 'LiFePO4'    ,'', '1urs10', merge_dics(dic,{'neb_outcar.normal.Li.Li':None, 'neb_outcar.replace.Li.K':4}) ] ,
               # ['LiFePO4', ise, 'FePO4', ise, 'LiFePO4'    ,'', '1m', merge_dics(dic,{'scaling_set.normal.Li.Li':'4mi'}) ] ,
               # ['NaFePO4', ise, 'FePO4', ise, 'NaFePO4'    ,scaling_set, main_set ] , # maricite structure, barrier is high

               # ['LiNiO2.r3m',   '1u', 'NiO2.r3m',  '', 'LiNiO2/r3m' , '', '1urns', md(dic,{'neb_outcar.normal.Li.Li':3, 'scaling_set.normal.Li.Li':'4uisC0W1', 'scaling_set.make_ds.Li.Li':'4uiAl',    }) ] , # for IS main_set is '1urns'
               # ['LiMn2O4.Fd3m',   '1u', 'Mn2O4.Fd3m', '', 'LiMn2O4/Fd3m' , '', '1ur10', md(dic,{'neb_outcar.normal.Li.Li':None, 'scaling_set.'+mode_id:'1ur10',}) ] ,
               # ['LiMn2O4.Fd3m',   '1u', 'Mn2O4.Fd3m', '', 'LiMn2O4/Fd3m' , '', '1ur10', md(dic,{'neb_outcar.normal.Li.Li':None, 'scaling_set.'+mode_id:'4uis',}) ] ,
               # ['LiVPO4F', ise, 'VPO4F', ise, 'LiVPO4F'    ,'',  main_set, md(dic,{'neb_outcar.normal.Li.Li':None, 'scaling_set.normal.Li.Li':'4uris', 'scaling_set.replace.Li.K':'4uris',  'scale_region.replace.Li.K':(-5,3),   'skip.replace.Li.K':0, 'skip.replace.Li.Na':1}) ]  ,# old
         
         ]


    for cat in cathodes:
        curres = {} #current result dict
        dic = cat[7] # dictionary with parameters



        name = cat_name_mod(cat[0])
        curres['name'] = name
        curres['DS']   = cat[2]
        curres['proto']= 'X'+cat[2]
        curres['ion']  = name.replace(cat[2], '')
        if not curres['ion']:
            # print('no')
            curres['ion'] = new_ion

        if 'make_ds' in mode:
            curres['x'] = '0%'
        else:
            curres['x'] = '100%'


        curres['vol'] = 0.001
        curres['id'] = ('','',1)


        if 'old.'+mode_id in dic: 
            old_behaviour = dic['old.'+mode_id]
        elif 'old_behaviour' in dic:
            old_behaviour = dic['old_behaviour']
        else:
            old_behaviour = 0
        # print ('old_behaviour', old_behaviour)

        if not dic['skip.'+ mode_id]:
            # continue

            if   'replace' in mode:
                
                itP0 = replace_cations(del_ion, new_ion, cat, curver)  
                
                idA0 = optimize_cell( (itP0, dic['scaling_set.'+mode_id], curver), cat, scale_regions[new_ion], update, irun = '1')

            elif 'make_ds'  in mode:
                
                itP0 = remove_cations(del_ion.split(), cat, curver, update)

                idA0 = optimize_cell( (itP0, dic['scaling_set.'+mode_id], curver), cat, scale_regions[del_ion.split()[-1]+'_removal'], update, irun = '1')

            elif 'normal'  in mode: #nothing done, use structures from provided list
                
                if old_behaviour:
                    idA0 = (cat[0], cat[1], curver)
                else:
                    #9.04.2017# idA0 = optimize_cell( (cat[0], cat[7]['scaling_set.'+mode_id], curver), cat, scale_regions['optimal'], update, irun = '1')
                    idA0 = optimize_cell( (cat[0], cat[1], curver), cat, scale_regions['optimal'], update, irun = '1')



            if run_sc and idA0:# and name != 'NiO2':
                
                make_dummy_calc_obj(idA0)


                idA = make_supercell(idA0, cat, updateA)
                
                # continue
                if run_scale:
                    if old_behaviour:
                        idB = optimize_cell(idA, cat, scale_regions['optimal'], updateB, irun = '2')
                    else:
                        idB = idA


                    if run_neb and idB:


                        if 'meps' in dic: # start_pos, end_pos
                            for mep in dic['meps']:
                                printlog('MEP ', mep, ':')
                                dic['start_pos'] = mep[0]
                                dic['end_pos']   = mep[1]
                                make_neb(idB, new_ion, cat, updateC)
                        elif 'meps2' in dic: # (i_atom_to_move, end_pos)
                            for mep in dic['meps2']:

                                dic['i_atom_to_move'] = mep[0]
                                dic['end_pos']   = mep[1]
                                make_neb(idB, new_ion, cat, updateC)

                        else:
                            make_neb(idB, new_ion, cat, updateC)


        results.append(curres)


    return results







def calc_antisite_defects(dpi = 300, image_format = 'eps', update = 0):
    struct_des = header.struct_des
    calc = header.calc
    readfiles = 1
    # update = 0
    main_set = '1urs'
    # ise_new = '1ur30'
    gga = 'gga'

    if gga == 'gga+u':
        choose_outcar = 0
    elif gga == 'gga':
        choose_outcar = 1


    cathodes = [


    # {'id':('LiFePO4.s10.su','4uis',100), 'cluster':'cee', 'main_set':'1ur30'},
    # {'id':('NaFePO4.ir.su.s10.su','4uis',100), 'cluster':'skol'},
    # {'id':('LiMnPO4.s10.su','4uis',100), 'cluster':'skol'},
    # {'id':('LiVPO4F.s10.su','4uis',100), 'cluster':'skol','main_set':'1ur30'},
    # {'id':('LiVP2O7.s10.su','4uis',100), 'cluster':'skol'},
    # {'id':('NaMnAsO4.s10.su','4uis', 100)  , 'cluster':'cee'    },
    # {'id':('LiTiS2.s10.su','4uis', 100)    , 'cluster':'cee'    },



    # {'id':('Li2FePO4F.ir.su.s10.su','4uis',100), 'main_set':'1u', 'cluster':'cee'},  #14 meV!!!
    # {'id':('Li2FePO4F.ir.su.s10.su','4uis',100), 'main_set':'1uAlf', 'cluster':'cee'}, #50 meV!!

    # {'id':('LiFePO4.s10.su','4uis',100), 'cluster':'cee', 'main_set':'1urZhang'},
    
    # {'id':('LiFePO4.s10.su','4mi',100), 'cluster':'cee', 'main_set':'1m'},
    # {'id':('LiFePO4.s10.su','4mi',100), 'cluster':'cee', 'main_set':'1msv'},
    # {'id':('LiFePO4.s10.su','4mi',100), 'cluster':'cee', 'main_set':'1mk15'},
    # {'id':('LiFePO4.s10.su','4mi',100), 'cluster':'cee', 'main_set':'1uAlf'},
    # {'id':('LiFePO4.s10.su','4mi',100), 'cluster':'cee', 'main_set':'1mAlf'},


    # {'id':('KVPO4F.Pna21.s10.su', '4uis', 100) , 'cluster':'cee'}, # energy is too high

    # {'id':('Na2FePO4F.s10.su','4uis', 100) , 'cluster':'cee'    }, #also incorrect groundstate
    # {'id':('Li2FePO4F.ir.su.s10.su','4uis',100), 'main_set':'1ur30', 'cluster':'skol'},  #incorrect groundstate
    # {'id':('KFeSO4F.s10.su','4uSi', 100), 'main_set':'1urS', 'cluster':'cee'}, #too much forces
    # {'id':('LiMn2O4.s10.su','4uis',100), 'cluster':'cee'},
    # {'id':('LiCoO2.s10.su','4uis',100), 'cluster':'cee'},
    # {'id':('Li2FePO4F.ir.su','4uis',100), 'cluster':'cee'},

    #pbcn
    # {'id':('Na2FePO4F.s10.su','4uis', 100) , 'main_set':'1uAlf', 'cluster':'cee'    }, 
    # {'id':('Na2FePO4F.s10.su','4uis', 100) , 'main_set':'1mAlf', 'cluster':'cee'    }, 
    # {'id':('NaFePO4F.id2.su.s10.su','4uis',100) , 'main_set':'1mAlf', 'cluster':'cee'    }, 
    # {'id':('NaFePO4F.id2.su.s10.su','4uis',100) , 'main_set':'1uAlf', 'cluster':'cee'    }, 

    #pbcn
    # {'id':('NaFePO4F.id2.su.s10.su','4uis',100),'st_from':calc['NaFePO4F.id1.su.s10.su','4uis',100].init, 'cation':'Na', 'trans':'Fe', 'main_set':'1mAlf', 'cluster':'cee'    }, 
    # {'id':('NaFePO4F.id2.su.s10.su','4uis',100),'st_from':calc['NaFePO4F.id1.su.s10.su','4uis',100].init, 'cation':'Na', 'trans':'Fe', 'main_set':'1uAlf', 'cluster':'cee'    }, 
    # {'id':('NaFePO4F.id2.su.s10.su','4uis',100),'st_from':calc['NaFePO4F.id1.su.s10.su','4uis',100].init, 'cation':'Na', 'trans':'Fe', 'main_set':'1urs', 'cluster':'cee'    }, 

    # {'id':('LiCoPO4F.pbcn.id2.su.s10','1u',100),'st_from':calc['LiCoPO4F.pbcn.id1.su.s10','1u',100].init, 'cation':'Li', 'trans':'Co', 'main_set':'1urs', 'cluster':'cee'    }, 
    # {'id':('LiFePO4F.id2.su.s10','1u',100),'st_from':calc['LiFePO4F.id1.su.s10','1u',100].init, 'cation':'Li', 'trans':'Fe', 'main_set':'1urs', 'cluster':'cee'    }, 
    # {'id':('LiFePO4F.id2.su.s10','1u',100),'st_from':calc['LiFePO4F.id1.su.s10','1u',100].init, 'cation':'Li', 'trans':'Fe', 'main_set':'1m', 'cluster':'cee'    }, 
    # {'id':('LiFePO4F.id2.su.s10','1u',100),'st_from':calc['LiFePO4F.id1.su.s10','1u',100].init, 'cation':'Li', 'trans':'Fe', 'main_set':'1mAlf', 'cluster':'cee'    }, 
    # {'id':('LiFePO4F.id2.su.s10','1u',100),'st_from':calc['LiFePO4F.id1.su.s10','1u',100].init, 'cation':'Li', 'trans':'Fe', 'main_set':'1uAlf', 'cluster':'cee'    }, 
    # {'id':('LiFePO4F.id2.su.s10','1u',100),'st_from':calc['LiFePO4F.id1.su.s10','1u',100].init, 'cation':'Li', 'trans':'Fe', 'main_set':'1urs', 'cluster':'cee'    }, 




    #pnma Co!

    # {'id':('LiCoPO4F.pnma.id2.su.s10','1urs',100),'st_from':calc['LiCoPO4F.pnma.id1.su.s10','1urs',100].init, 'cation':'Li', 'trans':'Co', 'main_set':'1urs', 'cluster':'cee'    }, 
    
    # {'id':('LiFePO4F.pnma.id2.su.s10','1u',100),'st_from':calc['LiFePO4F.pnma.id1.su.s10','1u',100].init, 'cation':'Li', 'trans':'Fe', 'main_set':'1urs', 'cluster':'cee'    }, 
    
    # {'id':('LiFePO4F.pnma.id2.su.s10','1u',100),'st_from':calc['LiFePO4F.pnma.id1.su.s10','1u',100].init, 'cation':'Li', 'trans':'Fe', 'main_set':'1m', 'cluster':'cee'    }, 





    ]

    results = {} # result dictionary, for each structure we have dict
    for cat in cathodes:
        # name = cat_name_mod(l[0])
        name = cat['id'][0].split('.')[0]
        results[name] = {}

    for cat in cathodes:
        # print(cat['id'])
        name = cat['id'][0].split('.')[0]

        if 'main_set' in cat:
            ise_new = cat['main_set']
        else:
            ise_new = main_set


        itA = cat['id'][0]+'.as'
        itB = cat['id'][0]+'.ifn'
        idA = (itA, ise_new, cat['id'][2])
        idB = (itB, ise_new, cat['id'][2])

        if 'st_from' in cat:
            #for creating in partly deintercalated structures
            ids = []
            # sufs = ['.a1', '.a2', '.a3']
            sufs = ['.a1', '.a3']
            for suf in sufs:
                it = cat['id'][0]+suf
                idd = (it, ise_new, cat['id'][2])
                ids.append( idd )
                if 'a1' in suf:
                    ida1 = idd
                elif 'a3' in suf:
                    ida3 = idd


            idsuf = ids[0]
        
        else:
            sufs = ['.as']
            idsuf = idA

        # print(idA not in calc, idsuf not in calc )
        # sys.exit()

        # if update or idA not in calc or idsuf not in calc:
        if update or idsuf not in calc:
            ''
            add_loop(*cat['id'], ise_new = ise_new, inherit_option = 'full_nomag', 
                it_folder = struct_des[cat['id'][0]].sfolder, 
                cluster = cat['cluster'], override = 1)

            # cat['st_from'] = calc[cat['id_frm']].init
            for suf in sufs:
                cat['mode'] = suf
                # print(suf)
                # sys.exit()
                add_loop(*cat['id'], ise_new = ise_new, 
                    inherit_option = 'antisite'+suf, confdic = cat,
                    it_folder = struct_des[cat['id'][0]].sfolder+'/anti/',
                    cluster = cat['cluster'], override = 1)
      



        else:
            res = None
            try:
                for idx in ids:
                    res_loop(*idx, up = 'up1',choose_outcar = choose_outcar,  show = 'fo', readfiles = readfiles)
                
                _, res = res_loop(*ida3, up = 'up1', analys_type = 'diff', show = 'fo', b_id = ida1,  choose_outcar = choose_outcar, readfiles = readfiles) #a1 - a3 - antisite energy
                
                # _, res = res_loop(*ids[1], up = 'up1', analys_type = 'diff', show = 'fo', b_id = idB,  choose_outcar = choose_outcar, readfiles = readfiles) #A - a2 move iron - energy of moving iron

                calc_redox(calc[ida1],  calc[idB] ) #voltage for alk ion intercalation
           
            except:
                _, res = res_loop(*idA, up = 'up1', analys_type = 'diff', show = 'fo', b_id = idB,  choose_outcar = choose_outcar, readfiles = readfiles)

            if res:
                results[name]['e_as'] = res[0]










    df = pd.DataFrame(results)
    df = df.transpose()
    df['source'] = gga
    print(df)
    if 1:
        dfs = pd.read_csv(r'database/literature.csv')[['is', 'e_as', 'source', 'functional']]
        dfs.drop(0, inplace=True)
        dfs.drop(1, inplace=True)
        dfs.dropna(inplace=True)
        dfs.set_index('is', inplace = True)
        dfs = dfs.apply(lambda x: pd.to_numeric(x, errors='ignore'))
        
        dfs['source'] = dfs['source']+' '+dfs['functional']
        del dfs['functional']
        # print(dfs)
        df = df.append(dfs)
        # print(df)
        

        df_piv = df.pivot(columns='source', values='e_as')
        print(df_piv)
        
        df_piv = df_piv.reindex(index = ['Li2FePO4F', 'LiFePO4', 'LiMnPO4', 'NaFePO4', 'NaMnAsO4', 'LiTiS2', 'LiVPO4F', 'LiVP2O7' ])

        df_piv.index = pd.Series([re.sub("([0-9])", "$_\\1$", s) for s in df_piv.index ] ) # #Make low indexes


        df_piv.plot(kind = 'bar', rot = 50, colormap = 'RdYlBu_r')

        plt.ylabel('Antisite formation energy (eV)')


        # plt.legend(loc='best') 
        plt.legend(loc=2, ncol =2, handlelength = 0.25, 
            handleheight = 2, frameon=False, fontsize=14)


        plt.tight_layout()

        # plt.show()
        figname = 'figs/antisite_'+gga#+'.'+image_format
        plt.savefig(figname+'.'+image_format, dpi = dpi, format=image_format)
        plt.savefig(figname+'.png', dpi = 300)

        # push_figure_to_archive(figname, 
        #     caption = gga+r""" Formation energies of antisite defects Li$_M^`M_{ \rm Li}^\bullet$ ($M$ = Fe, Mn, V) for several cathode materials""", 
        #     figlabel = 'as_'+gga, autocompl = False )

    return








def calc_antisite_defects3(update = 0, cathodes = None, param_dic = None, add_loop_dic = None, only = None):
    """

    'add' - element to be inserted in void!
    if provided
    only - only these configurations are considered
    """
    from impurity import insert_atom
    from analysis import find_polaron
    from geo import image_distance

    from current_structures import Na2X

    struct_des = header.struct_des
    calc = header.calc

    if not cathodes:
        cathodes = [
        # {'cl':calc[('Li2FePO4F.pnma','1u', 1)], 'el1':'Li', 'el2':'Fe', 'max_sep':4, 'set':'1uAlf', 'cluster':'cee'    }, 
        # {'cl':calc['Li2FePO4F.pnma.su.s10','1u', 100], 'el1':'Li', 'el2':'Fe', 'max_sep':4, 'set':'1uAlf', 'cluster':'cee'    }, 

        {'cl':Na2X, 'el1':'Na', 'el2':'Fe', 'max_sep':4, 'set':'1uAlf', 'cluster':'cee'    }, 
        # {'cl':calc['Na2FePO4F.s10.su','4uis', 100], 'el1':'Na', 'el2':'Fe', 'max_sep':4, 'set':'1uAlf', 'cluster':'cee'    }, 
        # {'cl':Na_X, 'el1':'Na', 'el2':'Fe', 'max_sep':4, 'set':'1uAlf', 'cluster':'cee'    }, 
        # {'cl':Na_X, 'el1':'Na', 'el2':'Fe', 'max_sep':4, 'set':'1uAlf', 'cluster':'cee', 'add':'Na' ,'i_void':3  }, #i_void: 0, 1
        ]

    if update:
        up = 'up2'
    else:
        up = 'up1'


    if not add_loop_dic:
        add_loop_dic = {}

    if param_dic:
        cathodes = [param_dic]

    # print(cathodes)
    for c in cathodes:

        st = c['cl'].end
        it = c['cl'].id[0]
        it_base = it

        if 'add' in c:
            el_add = c['add']
            st, i_add = insert_atom(st, el_add, i_void = c['i_void'],  r_imp = 1.7)
            it+=el_add+str(c['i_void'])
            add_loop(it, c['set'], 1, input_st = st, it_folder = struct_des[it_base].sfolder+'/as', up = up, **add_loop_dic)
        else:
            if 'iatom' in c:
                i_add = c['iatom']
            else:
                i_add = None

        sts = create_antisite_defect3(st, c['el1'], c['el2'], max_sep = c['max_sep'], iatom = i_add  )
        


        cl_base = c['cl'].run(c['set'], iopt = 'full_nomag', up = up, add = update, **add_loop_dic)

        header.show = 'fo'

        for i, st_as in enumerate(sts):
            if only is not None:
                if i not in only:
                    continue
            suf = 'as'+str(i)
            st_as.name+=suf
            st_as.write_poscar()
            add_loop(it+'.'+suf, c['set'], 1, input_st = st_as, it_folder = struct_des[it_base].sfolder+'/as', up = up, **add_loop_dic)
            cl_as = calc[it+'.'+suf, c['set'], 1]
            # print(cl_as.path['output'])
            try:
                print('dE(as) = {:.0f} meV'.format( (cl_as.energy_sigma0-cl_base.energy_sigma0)*1000))
                st = cl_as.end
                find_polaron(st, st_as.i_el1, out_prec = 2)
                print('Separation after relax = {:.2f} A'.format(   image_distance(st.xcart[st_as.i_el1], st.xcart[st_as.i_el2], st.rprimd )[0] )  )

            except:
                pass

            # if 'as0' in suf:
            #     break

    # sys.exit()
    return


def alkali_bar(m, param, option = ''):
    # bulk mod, volume change, average voltage

    if option in [ '1', '3']:
        m.set_index('name', inplace = True)
        try:
            m.drop('KNiO2', inplace =1)
        except:
            pass

    # print (m)

    if option != '3':
        m.set_index('proto', inplace = True)
    
    if option == '2': 

        m_piv = m.pivot(columns='x', values=param['val'])
    elif option == '3':
        # m_piv = m.pivot(columns='ion', values=param['val'])
        m_piv = m[['barrier']]
        # m_piv = m.pivot(columns='ion', values=param['val'])[['Li', 'Na', 'Na2', 'NaLi']]

    else:
        m_piv = m.pivot(columns='ion', values=param['val'])[['Li', 'Na', 'K',]]

    # if option == '2': 
    #     m.set_index('name', inplace = True)
    if option == '1':
        m_piv = m_piv.reindex(index = ['XFePO4', 'XMnPO4', 'XTiS2', 'XNiO2', 'XVP2O7', 'XTiO2' ])
    elif option == '2':
        m_piv = m_piv.reindex(index = ['XFePO4', 'XMnPO4', 'XVPO4F',  'XTiS2', 'XVP2O7', 'XTiO2', 'XMn2O4', ])
    
    elif option == '3':
        m_piv = m_piv.sort('barrier')
        pass

    else:
        print(m_piv)
        m_piv = m_piv.reindex(index = ['XFeSO4F', 'XVPO4F.Pna21','XVPO4F', 'XVP2O7',  'XTiS2', 'XMnPO4', 'XFePO4', 'XNiO2', 'XMn2O4',  'XTiO2' ])
    
    # m_piv.index['XFeSO4F'] = 
    m_piv.rename(index = {'XFeSO4F':'KFeSO4F', 'XVPO4F.Pna21':'KVPO4F.Pna21',  'XVPO4F':'LiVPO4F'}, inplace = 1)
    m_piv.index = pd.Series([re.sub("([0-9])", "$_\\1$", s.split('.')[0]) for s in m_piv.index ] ) # #Make low indexes
    if option == '2':
        m_piv.index = pd.Series([s.replace('X', 'Li') for s in m_piv.index ] ) # #Make low indexes



    m_piv.plot(kind = 'bar', rot = 60, colormap = 'RdYlBu')

    plt.legend(loc='best', ncol =3, handlelength = 0.25, 
        handleheight = 2, frameon=1, fontsize=20)

    # mpl.rcParams.update({'font.size': 22})
    plt.axhline(color='black')
    plt.axvline(x=1.5, color='black')

    if 0:
        plt.axvline(x=3.5, color='black')
        plt.axvline(x=5.5, color='black')
        s = 0.1; h = param['h']
        plt.text(0+s,  h,'Anatase\n         P1', ha='left', va='top', fontsize = 14)
        plt.text(1.5+s,h,'Olivine', ha='left', va='top', fontsize = 14)
        plt.text(3.5+s,h,r'P$\bar{1}$      Layered', ha='left', va='top', fontsize = 14)
        plt.text(5.5+s,h,'Pyrophospate\n          Pristine', ha='left', va='top', fontsize = 14)

    plt.xlabel('')
    plt.ylabel(param['yl'])
    plt.tight_layout()

    # plt.show()
    plt.savefig('figs/'+param['fig']+'.png', dpi = 300)
    plt.savefig('figs/'+param['fig']+'.eps', dpi = 300)
    # plt.close('all')
    # plt.clf()
    plt.show()

def cathode_screening():
    """
    Analyze results

    Notes:
    Eref is removed
    replace res_loop with calc_redox !!!

    """
    calc = header.calc

    def sort_func(s):
        s = s.replace('Li', 'A')
        s = s.replace('Na', 'B')
        s = s.replace('K',  'C')
        s = s[::-1]
        # print (s)
        return s


    if 'main' not in calc:
        calc['main'] = {}
    df = calc['main']

    if 1:
        if 1:
            ''
            # df['Li_IS'] = calc_barriers('normal') #Basic intercalated structures
            # df['IS'] = calc_barriers('normal') # all, except LiCoO2, KFeSO4F

            # df['Li_IS_gga'] = calc_barriers('normal', func = 'gga') # U = 0, LiFePO4 LiMnPO4 LiNiO2 LiVPO4F NaMnAsO4 Na2FePO4F 
            # df['Li_IS2'] =calc_barriers('normal') #NaMnAsO4, Na2FePO4F, Na2FeVF7, NaLiCoPO4F, LiNiO2, LiMn2O4
            
            # df['K_IS2'] =  calc_barriers('normal') #KVPO4, KFeSO4

            calc_barriers('normal') 

            # df['Na_IS'] = calc_barriers('replace', 'Li', 'Na') #Na IS obtained from Li=base
            # df['K_IS'] = calc_barriers('replace', 'Li', 'K')  #K  IS obtained from Li=base
            # calc_barriers('replace', 'Li', 'Na')  #K  IS obtained from Li=base
        # calc_barriers('replace', 'Na', 'K')  #K  IS obtained from Li=base
        # calc_barriers('replace', 'Na', 'Li')  #K  IS obtained from Li=base

        if 1:
            ''
            # df['Li_DS'] =  calc_barriers('make_ds',  'Li', 'Li') #Basic deintercalated structures
            # calc_barriers('make_ds',  'Li', 'Li') #Basic deintercalated structures
            # df['DS_gga'] =  calc_barriers('make_ds',  'Li', 'Li', func = 'gga') #Basic deintercalated structures
            # sys.exit()
            # calc_barriers('make_ds',  'Na', 'Na') #Basic deintercalated structures
            # df['K_DS'] = calc_barriers('make_ds',  'K', 'K') #Basic deintercalated structures
            # calc_barriers('make_ds',  'Li Na', 'Li') #Basic deintercalated structures
        # save = 1
        # print(df['Li_DS'][0]['vol'])
    else:
        #calculate change of volume due to deintercalation
        # print(df)
        al = df['IS']
        a = df['Li_IS']
        a2=  df['Li_IS2']
        b = df['Na_IS']
        c = df['K_IS']+df['K_IS2']
        
        d = df['Li_DS']+df['K_DS']
        e = df['Li_IS_gga']
        f = df['DS_gga']



        for x in a,b,c:
            for ix, jd  in zip(x,d):
                
                if jd['name'] in ix['name']:
                    # ix['vol']
                    # ix['vol_red'] = (ix['vol']/jd['vol']-1)*100 #reduction of volume in % due to deintercalation
                    dic, _ = res_loop(*ix['id'], b_id = jd['id'], energy_ref = ix['Eref'], analys_type = 'redox_pot', readfiles = 0)
                    if 'vol_red' in dic:
                        ix['vol_red'] = dic['vol_red']
                        ix['Vav'] = dic['redox_pot']
                        print('Volumes are', calc[ix['id']].end.vol, calc[jd['id']].end.vol,'; ids are',ix['id'], jd['id'] )

                    else:
                        ix['vol_red'] = None
                    ix['vol_ds']  = jd['vol']

                    # print (ix['name'])
                    if ix['name'] == 'KVP2O7':
                        ix['vol_red'] = None
                        ix['Vav']     = None


                else:
                    printlog('DS structure  ', jd['name'], 'is not in intercalated structure', ix['name'])
                    ix['vol_red'] = None
                ''




        al = pd.DataFrame(al)
        a = pd.DataFrame(a)
        a2 = pd.DataFrame(a2)
        # a.dropna('LiMn2O4')
        b = pd.DataFrame(b)
        c = pd.DataFrame(c)
        d = pd.DataFrame(d)
        e = pd.DataFrame(e)
        f = pd.DataFrame(f)


        # print(b)
        # sys.exit()


        # print(a.round(1))
        m = pd.concat([a, b, c])
        md = pd.concat([a, b, c, d])
        # print(m)
        # m = pd.DataFrame(m, index =  sorted(m.index, key = sort_func  )  ) #sort
        # print(m)
        
        # print(m)
        # sys.exit()



        # alkali_bar(m, {'val':'B',       'h':135, 'yl':'Bulk modulus (GPa)', 'fig':'bulk' } ) #plot Bulk modulus
        # alkali_bar(m, {'val':'vol_red', 'h':440,  'yl':'Volume change (%)' , 'fig':'volume' } )
        # alkali_bar(m, {'val':'Vav', 'h':40,  'yl':'Intercalation voltage (V)' , 'fig':'av_pot' })
        # alkali_bar(m, {'val':'barrier',       'h':135, 'yl':'Migration barrier (eV)', 'fig':'barrier_Li_Na_K' }, option = '1' ) 
        # alkali_bar(pd.concat([a, d]), {'val':'barrier',       'h':135, 'yl':'Migration barrier (eV)', 'fig':'barrier_Li_IS_DS' }, option = '2' ) 
        # alkali_bar(a2, {'val':'barrier',       'h':135, 'yl':'Migration barrier (eV)', 'fig':'barrier_Li_IS2' }, option = '3' ) 

        #comparison of gga and gga+
        if 1:
            if 1: #IS
                al['func'] = 'GGA+U'
                e['func']  = 'GGA'
                mc = pd.concat([al, e])
                filename = 'barriers_gga_ggau'
                # print(mc)
            else: #DS 
                d['func'] = 'GGA+U'
                f['func'] = 'GGA'
                mc = pd.concat([d, f])
                filename = 'barriers_gga_ggau_DS'

            mc.set_index('name', inplace = True)
            m_piv = mc.pivot(columns='func', values='barrier')
            m_piv.dropna(inplace =1)

            # print(m_piv)
            # sys.exit()
            m_piv = m_piv.reindex(index = ['Na2FePO4F', 'NaMnAsO4', 'LiVPO4F', 'LiFePO4', 'LiMnPO4',  'LiNiO2', 'LiVP2O7' ])


            m_piv.index = pd.Series([re.sub("([0-9])", "$_\\1$", s) for s in m_piv.index ] ) # #Make low indexes

            m_piv.plot(kind = 'bar', rot = 60, colormap = 'RdYlBu_r')

            plt.legend(loc='best', ncol =3, handlelength = 0.25, 
                handleheight = 2, frameon=1, fontsize=20)
            plt.xlabel('')
            plt.ylabel('Migration barrier (eV)')
            plt.ylim(0, 1.2)

            plt.tight_layout()

            plt.savefig('figs/'+filename+'.png', dpi = 300)
            plt.savefig('figs/'+filename+'.eps', dpi = 300)
            plt.show()        


        if 0:
            dfs = pd.read_csv(r'database/literature.csv')[['is', 'barrier', 'source', 'type']]
            dfs.drop(0, inplace=True)
            dfs.dropna(inplace=True)
            dfs.set_index('is', inplace = True)
            dfs = dfs.apply(lambda x: pd.to_numeric(x, errors='ignore'))
            # dfs['owner'] = 'world'    
            # print(dfs.sort_index())
            # print(dfs)
            # sys.exit()

            # m = 

            m['source'] = 'our'
            # m['owner'] = 'Skoltech'    
            m['type'] = 'our'   
            # print(m)

            # m.dropna(inplace =1)
            m = m[['name', 'barrier', 'source', 'type']]
            m.set_index('name', inplace = 1)
            if 1:
                m = m.append(dfs)

            m.dropna(inplace =1)


            # m.sort_index(inplace = 1)

            # print (m)
            m_piv = m.pivot(columns='type', values='barrier')[['our', 'theory']]
            # print( m_piv)
            m_piv.dropna(inplace =1)


            # df_piv = pd.DataFrame(df_piv, index =  sorted(df_piv.index, key = sort_func  )  )

            if 1:
                m_piv.drop('LiMn2O4', inplace =1)
                m_piv = m_piv.reindex(index = ['LiFePO4', 'NaFePO4', 'LiMnPO4', 'LiVPO4F', 'LiTiS2', 'LiTiO2' ])

                m_piv.rename(columns = {'theory':'Literature DFT'}, inplace =1)

                m_piv.index = pd.Series([re.sub("([0-9])", "$_\\1$", s) for s in m_piv.index ] ) # #Make low indexes

                m_piv.plot(kind = 'bar', rot = 60, colormap = 'RdYlBu')

                plt.legend(loc='best', ncol =3, handlelength = 0.25, 
                    handleheight = 2, frameon=1, fontsize=20)

                plt.ylim(0, 1.6)

                plt.xlabel('')
                plt.ylabel('Migration barrier (eV)')
                plt.tight_layout()

                plt.savefig('figs/barriers_theory.png', dpi = 300)
                plt.savefig('figs/barriers_theory.eps', dpi = 300)
                plt.show()


def alkali_bar2(m, param = None, option = '', plot_type = 'bar', suf = '', reindex = None):
    """
     bulk mod, volume change, average voltage
    plot bar figures for paper7

    m - dataframe object with data build from output of calc_barriers()
    param
        - val - y value to plot
            - barrier
        - yl - y legend
        - h - high of figure
        - fig - figure name

    option - 
        1 - compare Li Na K
        2 - compare IS and DS on one plot

    reindex - change the order of index - for each project is added additionaly


    plot_type 
        bar
        dAO_barrier - specific plot
    
    suf (str) - arbirtary suffix

    param['val'] - ['redox', 'vol_red']

    """
    m = copy.deepcopy(m)
    if param == None:
        param = {'val':None}

    # if plot_type == 'dAO_barrier':
    #     param = None

    if option in [ '1', '3']:
        # print(m['name'])
        # sys.exit()
        m.set_index('name', inplace = True)
        # sys.exit()
        try:
            m.drop('KNiO2', inplace =1)
        except:
            pass



    # print(m)

    if 'barrier' in param['val']:
        printlog('Attention! Barriers were converted to absolute values', imp = 'y')
        # convert negative barriers to positive!
        m.loc[:, param['val']] = m[param['val']].abs()
        # print(m[param['val']])

        # sys.exit()


    m_x = m.set_index('x')
    mIS = m_x.ix['100%']
    mDS = m_x.ix['0%']


    # print(mIS)
    # print(mDS)



    ms = []
    if option in ['1', '2']:
        m.set_index('proto', inplace = True)
        mIS.set_index('proto', inplace = True)
        mDS.set_index('proto', inplace = True)

        


    #determine redox pot and volume change
    if plot_type == 'bar':

        redox, volred = [], []
        if param['val'] in ['redox', 'vol_red']:
            
            # mIS.ix[0]['new'] = 'go'
            # print(mIS)
            for rowIS, rowDS in zip(mIS.iterrows(), mDS.iterrows()):
                idIS = rowIS[1]['id']
                idDS = rowDS[1]['id']
                # print(idIS, idDS)
                out = calc_redox(db[idIS], db[idDS], silent = 1)
                redox.append(out['redox_pot'])
                volred.append(out['vol_red'])


            mIS = mIS.assign(redox=pd.Series(redox,    index=mIS.index)) # add column
            mIS = mIS.assign(vol_red=pd.Series(volred, index=mIS.index)) # add column
        # print(mIS)
        # sys.exit()



        if option == '1':
            m_pivIS = mIS.pivot(columns='ion', values=param['val'])[['Li', 'Na', 'K',]]
            ms.append(m_pivIS)

            if param['val'] in ['barrier', ]: #use negative axis for DS state


                mDS.loc[:, param['val']] = mDS[param['val']]*-1

                m_pivDS = mDS.pivot(columns='ion', values=param['val'])[['Li', 'Na', 'K',]]

                ms.append(m_pivDS)





        if option == '2': 
            m_piv = m.pivot(columns='x', values=param['val'])
            ms.append(m_piv)
        

        elif option == '3':
            # m_piv = m.pivot(columns='ion', values=param['val'])
            m_piv = m[['barrier']]
            # m_piv = m.pivot(columns='ion', values=param['val'])[['Li', 'Na', 'Na2', 'NaLi']]
            ms.append(m_piv)

        else:
            ''
            # m_piv = m.pivot(columns='ion', values=param['val'])[['Li', 'Na', 'K',]]

            # ms.append(m_piv)
        # if option == '2': 
        #     m.set_index('name', inplace = True)
        
        # print(m_pivDS)


        if option == '1':
            ''
            # m_piv = m_piv.reindex(index = ['XFePO4', 'XMnPO4', 'XTiS2', 'XNiO2', 'XVP2O7', 'XTiO2' ])
        elif option == '2':
            m_piv = m_piv.reindex(index = ['XFePO4', 'XMnPO4', 'XVPO4F',  'XTiS2', 'XVP2O7', 'XTiO2', 'XMn2O4', ])
        
        elif option == '3':
            m_piv = m_piv.sort('barrier')
            pass

        else:
            # print(m_piv)
            m_piv = m_piv.reindex(index = ['XFeSO4F', 'XVPO4F.Pna21','XVPO4F', 'XVP2O7',  'XTiS2', 'XMnPO4', 'XFePO4', 'XNiO2', 'XMn2O4',  'XTiO2' ])
        

        for i in range(len(ms)):
            mp = ms[i]
            mp.rename(index = {'XMn2O4.Fd3m.111':'XMn2O4','XFeSO4F':'KFeSO4F', 'XVPO4F.Pna21':'KVPO4F.Pna21', }, inplace = 1)
            mp.rename(index = {'XMO':'XMn2O4','XFPO':'XFePO4', 'XMPO':'XMnPO4', 'XVPOF':'XVPO4F', }, inplace = 1)
     
            ms[i] = mp.reindex(index = ['XMn2O4', 'XFePO4', 'XMnPO4', 'XVPO4F'])
            mp = ms[i]

            mp.index = pd.Series([re.sub("([0-9])", "$_\\1$", s.split('.')[0]) for s in mp.index ] ) # #Make low indexes

            # print(mp.index)

            mp.index = pd.Series([s.replace('X', 'A') for s in mp.index ] ) #
            # print(mp.index)


        if option == '2':
            m_piv.index = pd.Series([s.replace('X', 'Li') for s in m_piv.index ] ) #


    # print(m_piv)
    # mpl.style.use('ggplot')

    # fig, ax = plt.subplots()
    

    if 'ylim' in param:
        ylim = param['ylim']
    else:
        ylim = None

    if plot_type == 'bar':
        ax = None


        for mp in ms:
            if ax != None:
                ax.legend('')
            # print(mp)
            ax = mp.plot(kind = 'bar', rot = 60, colormap = 'RdYlBu', ax = ax, ylim = ylim)
        
        # plt.legend(['Li', 'Na', 'K'])

        ax.legend(['Li', 'Na', 'K'],loc='best', ncol =3, handlelength = 0.25, 
            handleheight = 2, frameon=1, fontsize=20)

        # mpl.rcParams.update({'font.size': 22})
        plt.axhline(color='black')
        plt.axvline(x=1.5, color='black')

        if 0:
            plt.axvline(x=3.5, color='black')
            plt.axvline(x=5.5, color='black')
            s = 0.1; h = param['h']
            plt.text(0+s,  h,'Anatase\n         P1', ha='left', va='top', fontsize = 14)
            plt.text(1.5+s,h,'Olivine', ha='left', va='top', fontsize = 14)
            plt.text(3.5+s,h,r'P$\bar{1}$      Layered', ha='left', va='top', fontsize = 14)
            plt.text(5.5+s,h,'Pyrophospate\n          Pristine', ha='left', va='top', fontsize = 14)

        plt.xlabel('')
        plt.ylabel(param['yl'])

        filename = param['fig']

    elif plot_type == 'dAO_barrier':
        # print(m_x)
        m_dAO_IS = m_x.set_index('dAO_change')
        # dAO_bar = mIS[['dAO_change', 'barrier']]
        # print(dAO_bar)  

        # dAO_bar.plot(x= 'dAO_change',y = 'barrier', xlim = (0, 0.3), style = ['o'])
        m_dAO_IS_piv = m_dAO_IS.pivot(columns='ion', values='barrier')[['Li', 'Na', 'K',]]


        # m_dAO_IS_piv.columns.name = 'A-(O,F) distance, $\AA$'
        m_dAO_IS_piv.index.name =  'Change of A-(O,F) distance, $\AA$'
        # print(m_dAO_IS_piv)  



        header.mpl.rc('legend', fontsize= 16) 

        m_dAO_IS_piv.plot(xlim = (0, 0.3), style = ['ro', 'gs', 'bv'], ms = 10, )
        
        plt.ylabel('Migration barrier, eV')
        # plt.show()
        filename = param['fig']

        # filename = plot_type+suf

    plt.tight_layout()


    # plt.show()
    plt.savefig('figs/'+filename+'.png', dpi = 300)
    plt.savefig('figs/'+filename+'.pdf', dpi = 300)
    # plt.close('all')
    # plt.clf()
    # plt.show()
    return

def e_bind(cl_bulk, cl_vac, cl_sol, cl_compl):
    #binding energy of vacancy - solute complex
    #cl_compl - complex
    dE = cl_sol.e0 + cl_vac.e0 - (cl_bulk.e0 + cl_compl.e0)
    
    if cl_compl.end.nznucl[1] > 1:
        n_sol = cl_compl.end.nznucl[1]
        print('Number of solute atoms = ', n_sol)
        dE = n_sol*cl_sol.e0 + cl_vac.e0 - (n_sol*cl_bulk.e0 + cl_compl.e0)



    return dE




def plot_UvsD(it1, it2, it_b, ise, sc_reg, n_imag, ise_b = None, v_b = None, suf = '', 
    up = 0, legend = None, invert_x = False, linetypes = None,  
    st_start = None, ver_lines = None, ylim = None):
    """
    Intercalation potential vs deformation
    sc_reg - scale region
    n_imag - n_scale_images
    suf - suffix
    v_b (int) - if provided than version of base fixed to this number
    """
    plot_dist = False

    if plot_dist:
        show = 'sur'
    else:
        show = 'fo'


    scales = np.linspace(*sc_reg, n_imag)
    U1, U2 = [], []

    if not ise_b:
        ise_b = ise

    for it in [it1, it2, it_b]:
        # print(calc[it, ise, 2].state)
        if it == it_b:
            iset = ise_b
        else:
            iset = ise
        if '4' not in calc[it, iset, 1].state or up :
            res_loop(it, iset, range(1, n_imag+1), show = show)
            # print(calc[it, ise, 2].state)
    # sums = {'Li-O':[], 'Na-O':[], 'Li-O(b)':[], 'Na-O(b)':[]}
    sums = {}
    for v in range(1,n_imag+1):
        cl1 = calc[it1, ise, v]
        cl2 = calc[it2, ise, v]
        if v_b:
            cl_b = calc[it_b, ise_b, v_b]
        else:
            cl_b = calc[it_b, ise_b, v]
        U1.append(calc_redox(cl1, cl_b, value = -scales[v-1])['redox_pot'] )
        U2.append(calc_redox(cl2, cl_b, value = -scales[v-1])['redox_pot'] )

        if plot_dist:
            for k, cl in ('1', cl1), ('2', cl2), ('b', cl_b):
                for bond in ('Li-O', 'Na-O', 'Fe-O', 'O-O', 'Li-Fe', 'Na-Fe'):
                    if bond in cl.sumAO:
                        key = bond+'('+str(k)+')'
                        if key not in sums:
                            sums[key]= []
                        sums[key].append(cl.sumAO[bond])

        # cl1.end.write_xyz()
        # cl2.end.write_xyz()
        # cl_b.end.write_xyz()

    if not legend:
        legend = ('pos1', 'pos2', 'lower left')
    elif len(legend) == 2:
        legend.append('lower left')

    if invert_x:
        scales = [-s for s in scales ] 
        # U1 = U1[::-1]
        # U2 = U2[::-1]


    if not linetypes:
        linetypes = ['b-o', 'g-o']

    if 1:
        fit_and_plot(U1 = (scales, U1, linetypes[0], legend[0]), U2 = (scales, U2, linetypes[1], legend[1]), 
            image_name = 'figs/strain_'+suf, legend = legend[2], ylabel = 'Redox potential, V', xlabel = 'Compression, $\delta$, %',
            ylim = ylim,
            ver_lines = ver_lines)









    if plot_dist:
        header.mpl.rc('legend', fontsize= 10) 

        fit_and_plot(
                    # l1 = (scales, sums['Li-O(1)'], 'g--', 'Li-O'), 
                    # l2 = (scales, sums['Li-O(b)'], 'g-', 'Li-O(b)'), 
                    # l3 = (scales, sums['Na-O(2)'], 'r--', 'Na-O'), 
                    # l4 = (scales, sums['Na-O(b)'], 'r-', 'Na-O(b)'), 
                    l5 = (scales, sums['Fe-O(1)'], 'g--', 'Fe-O(1)'), 
                    l6 = (scales, sums['Fe-O(2)'], 'r--', 'Fe-O(2)'), 
                    l7 = (scales, sums['Fe-O(b)'], 'r-', 'Fe-O(b)'), 
                    # l5 = (scales, sums['O-O(1)'], 'g--', 'O-O(1)'), 
                    # l6 = (scales, sums['O-O(2)'], 'r--', 'O-O(2)'), 
                    # l7 = (scales, sums['O-O(b)'], 'k-',  'O-O(b)'), 
                    
                    # l5 = (scales, sums['Li-Fe(1)'], 'g--', 'Li-Fe(1)'), 
                    # l6 = (scales, sums['Na-Fe(2)'], 'r--', 'Na-Fe(2)'), 
                    # l7 = (scales, sums['Li-Fe(b)'], 'g-',  'Li-Fe(b)'), 
                    # l8 = (scales, sums['Na-Fe(b)'], 'r-',  'Na-Fe(b)'), 


                    # ld = (scales, np.array(sums['Na-O'])  - np.array(sums['Na-O(b)']), 'r--', 'd Na-O'), 
                    # ld2 = (scales, np.array(sums['Li-O']) - np.array(sums['Li-O(b)']), 'g-', 'd Li-O'), 


            image_name = 'figs/AO_dist_'+suf, legend = 'best', ncol = 1, ylabel = 'average A-O, A', xlabel = 'Compression, $\delta$, %', ylim = None)



    return scales, U1, U2, legend


def calc_strain_influence(cl_b, replace = False, rep_pos = None, 
    half = False, del_el_list = None, del_pos_list = None, 
    invert_x = None, vac = False, it_folder = 'Na2FePO4F/scaled', 
    plot = 0, ise = '', mul_matrix = None, sreg = None, 
    suf = None, up = 0, n_scale_images = None,
    ver_lines = None, ylim = None):
    """
    Helper function
    Creates 3 calculations
    replace - if replacement of one Na by Li is needed
    rep_pos - replacement pos of Na with Li
    del_pos_list - del pos
    vac - if False all atoms are removed, if True only one atom is removed
    half - (bool) remove half of atoms (half charged)
    sreg - scale_region
    """
    scale_region = sreg
    if not mul_matrix:
        mul_matrix = [[0.99, 0, 0], [0, 1.0025, 0], [0, 0, 0.98]]
    if not scale_region:
        scale_region = (-3, 8)

    if not ylim:
        ylim = (2.6, 3.5)


    if not n_scale_images:
        n_scale_images = 12
    
    cl_list = []
    if not ise:
        ise = '1uis'

    if not del_el_list:
        del_el_list = ['Na', 'Li']

    if not del_pos_list:
        del_pos_list = [1,]


    if replace:
        st_b = create_replaced_structure(cl_b.end, el1 = 'Na', el2 = 'Li', rep_pos = rep_pos)
        it_b = cl_b.id[0]+'.ir'+str(rep_pos)+'Li'


    else:
        st_b = copy.deepcopy(cl_b.end)
        it_b = cl_b.id[0]
    

    st_b.magmom = [None]
    cl_list.append([it_b, st_b])
    # names = []
    for el in del_el_list:
        for del_pos in del_pos_list:
            lsuf = ''
            if vac:
                # del_pos = 1
                st_del = remove_one_atom(st_b, el, del_pos = del_pos)
                it_del = it_b+'.vac'+el+str(del_pos)

            elif half:
                lsuf = 'half'
                if len(del_el_list) == 2:
                    if el == 'Na':
                        at_to_remove = [61, 60, 57,64]
                        conf = 'c6'
                        st_del = st_b.remove_atoms(at_to_remove, from_one = 1)
                    elif el == 'Li':
                        at_to_remove = [65, 72, 69, 68]
                        conf = 'c5'
                        st_del = st_b.remove_atoms(at_to_remove, from_one = 1)
                    it_del = it_b+'.id'+el+lsuf+'.'+conf # 


                elif len(del_el_list) == 1:
                    if del_pos == 1:
                        at_to_remove = [61,62,57,59]
                        conf = 'c13'
                        st_del = st_b.remove_atoms(at_to_remove, from_one = 1)
                    elif del_pos == 2:
                        at_to_remove = [65,72, 69,68]
                        conf = 'c5'
                        st_del = st_b.remove_atoms(at_to_remove, from_one = 1)                        

                    it_del = it_b+'.id'+el+str(del_pos)+lsuf+'.'+conf #



            else:
                st_del = create_deintercalated_structure(st_b, el, del_pos = del_pos)
                
                if len(del_pos_list) ==1:
                    it_del = it_b+'.id'+el # for compatibility 
                else:
                    it_del = it_b+'.id'+el+str(del_pos)
        



            cl_list.append([it_del, st_del])






    if not suf:
        suf = '.sm'
    scales, U1, U2, legend = None, None, None, None

    if not plot:
        for clt in cl_list:
            idd = (clt[0], ise, 1)
            iddsm = (clt[0]+suf, ise, 1)
            if iddsm not in calc or up:
                add_loop(*idd, input_st = clt[1], it_suffix = suf.replace('.', ''), calc_method = 'scale', mul_matrix = mul_matrix,
                 scale_region = scale_region, n_scale_images = n_scale_images,  inherit_option = 'inherit_xred', it_folder = it_folder)
    else:





        if rep_pos == 1:
            pos2 = 2
            # lin1 = '-b'
            # lin2 = '-g'
        elif rep_pos == 2:
            pos2 = 1
            # l1 = lin1
            # lin1 = lin2
            # lin2 = l1




        if len(del_el_list) == 2:
            leg2 = 'Li'+str(rep_pos) +' removal'
            leg1 = 'Na'+str(pos2)    +' removal'
        
        elif len(del_el_list) == 1:
            el = del_el_list[0]
            leg1 = el+str(del_pos_list[0]) +' removal'
            leg2 = el+str(del_pos_list[1]) +' removal'

        lin1 = '-ob'
        lin2 = '-og'
        
        # if 'Li' in leg2:
        #     lin2 = '-og'
        # elif 'Na' in leg2:
        #     lin2 = '-ok'


        scales, U1, U2, legend = plot_UvsD(it1 = cl_list[1][0]+suf , it2 = cl_list[2][0]+suf, it_b = cl_list[0][0]+suf, ise = ise, sc_reg = scale_region, 
            n_imag = n_scale_images, 
            suf = cl_list[2][0]+suf+ise, legend = (leg1, leg2, 'lower left'), linetypes = (lin1, lin2), up = 1 , 
            invert_x = invert_x, st_start = cl_b.end, ver_lines = ver_lines, ylim = ylim)


# ise, sc_reg, n_imag
    return scales, U1, U2, legend



def replace(update, cl, el1, el2, reptype = None, show_fit = 0):
    """
    cl - calculation to work with
    el1 - element to be replaced
    el2 - replace by this
    reptype - type of replacement
        'full' - all atoms of specific symmetry type
        'one' - only one atom of specific symmetry type

    update - update 
    """

    st = cl.end
    # print(determine_symmetry_positions(st, el1))
    n_noneq_pos = len(determine_symmetry_positions(st, el1))
    printlog('Number of non-equiv pos', n_noneq_pos)


    # sys.exit()
    if not reptype:
        printlog('Error! replace(): Please provide reptype = "full" or "one"')

    for pos in range(n_noneq_pos):
        ''
        pos+=1
        if 'one' in reptype:
            suf = 'one'
            ise = '1uis'
            scale_region = (-4, 4)
            one = 1
        else:
            suf = ''
            ise = '4uis'
            scale_region = (-6, 2)
            one = 0

        it_new = cl.id[0]+'.ir'+str(pos)+suf+el2
        id_new = (it_new, ise, 1)
        it_folder = header.struct_des[cl.id[0]].sfolder + '/' + 'replaced'

        # print((it_new+'.su', ise, 1))
        # print((it_new+'.su', ise, 1) not in header.calc)
        # sys.exit()
        if one:
            id_new2 = id_new
        else:
            id_new2 = (it_new+'.su', ise, 1)

        if update or id_new2 not in header.calc:
            st_rep = create_replaced_structure(st, el1 = el1, el2 = el2, rep_pos = pos, only_one = one)
            st_rep.write_xyz()
            
            if one:
                add_loop(*id_new, input_st = st_rep, it_folder = it_folder, override = True) #charge density without relaxation

            else:
                add_loop(*id_new, calc_method = 'uniform_scale', inherit_option = 'inherit_xred', scale_region = scale_region, input_st = st_rep, it_folder = it_folder)
        



        else:
            ''
            if 'full' in reptype:
                if show_fit:
                    res_loop(it_new+'.su', ise, list(range(1,8))+[100], analys_type = 'fit_a', show = 'fitfo')
                else:
                    res_loop(it_new+'.su', ise, 100, show = 'fo')

            else:
                res_loop( *id_new, show = 'fo')


    return


calc_material = calc_barriers



def neb_wrapper( param_dic = None, paths = None, run_neb = 0, read = 0, plot = 0, mode = 'IS', DS = None, ylim = None, 
    special_case = None, add_loop_args = None, first =1, last = 1, style_dic = None, substitute = None, up_res = 'up2'):
    """
    Allows to run, read, and plot composition diffusion path;
    was used for RbVPO4F and KVPO4F 


    paths - list of tuples; each tuple contains (A, B, C, D, E),
        where A - #atom number from 0 - start position for migration
              B - end position number according add_neb
              C - name of created it
              D - (not nessesary) either path to occupation matrix or
                               id, from which the occupation matrix should be taken for IS
              E - (not nessesary) either path to occupation matrix or
                               id, from which the occupation matrix should be taken for DS

    special_case: # for KVP additional static run
    mode - 
        'IS', 'DS'

    DS - calc used as initial ds





    substitute - manually substitute specific points of MEP

    """
    from picture_functions import plot_mep



    struct_des = header.struct_des
    calc = header.calc
    
    if not add_loop_args:
        add_loop_args = {}

    if 'upC' not in param_dic:
        param_dic['upC'] = None

    n_set = param_dic['neb_set']
    imag = param_dic['images']

    nameadd = '_'+param_dic['id'][0]+'_'+param_dic['main_set']+'_'+n_set
    tkeys = {'1':'t_mep'+nameadd, '2':'t_it'+nameadd,'3':'t_mepDS'+nameadd}
    for key in tkeys:
        if tkeys[key] not in calc:
            # print(key)
            calc[tkeys[key]] = {}


    pos = []
    mep = []

    el = param_dic['el']
    # if 'old_behaviour' in param_dic:
    #     old_behaviour = param_dic['old_behaviour']
    # else:
    #     old_behaviour = 0

    if run_neb or read:
        plot = 0

    for p in paths:
       
        if len(p)>3:
            if 'OCCMATRIX' in p[3]: # the file is provided
                add_loop_args['params'] = {'occmatrix':p[3]}
            else: # assuming that id of calculation is provided
                param_dic['occmatrix_id'] = p[3]
        pk = p[0:3]

        if mode == 'IS':
            """intercalated"""

            if len(p)>3:
                if 'OCCMATRIX' in p[3]: # the file is provided
                    add_loop_args['params'] = {'occmatrix':p[3]}
                else: # assuming that id of calculation is provided
                    param_dic['occmatrix_id'] = p[3]
            pk = p[0:3]

            t_key = tkeys['1']
            # print(t_key)
            param_dic['i_atom_to_move'] = p[0]
            param_dic['start_pos'] = None
            param_dic['end_pos']   = p[1]  
            if run_neb or read:
                calc_barriers('normal',  el, el, up_res = up_res, run_neb = 1, show_fit = 0, up = 0, upA = 0, upC = param_dic['upC'], param_dic = param_dic,  add_loop_dic = add_loop_args) # here 1 is 4.47 A; 3 is 6.64 A
                # print(t_key)
                calc[t_key][pk] = calc['_mep']
                # print(calc[t_key][p])
            
            # print (t_key) 
            # print (calc[t_key]) 
            # print(calc[t_key][p][1])
            print('t_key', t_key)
            pos.extend(list((calc[t_key][pk][0])))

            if special_case: # for KVP additional static run
                e_stat = [calc[pk[2]+'.ifn', '0m', i].e0 for i in [1, 3,4,5,6,7, 2]]
                mep.extend(list(reversed(e_stat)))
                suf2= '_spec'

            else:
                mep.extend(list(reversed(calc[t_key][pk][1])))
                suf2 = ''



        elif mode == 'DS':


            if len(p)>4:
                if 'OCCMATRIX' in p[4]: # the file is provided
                    add_loop_args['params'] = {'occmatrix':p[4]}
                else: # assuming that id of calculation is provided
                    param_dic['occmatrix_id'] = p[4]

                cl = db[p[4]]
                # print(p[4])

                occfile = write_occmatrix(cl.occ_matrices, cl.dir)
                # continue
                add_loop_args['params'] = {'occmatrix':occfile}


            pk = p[0:3]


            suf2 = ''
            """deintercalated"""
            for cl in (
            DS, 
            # Rb05V12, 
            # Rb05V22
            ):
                ''
                t_key_it = tkeys['2']
                t_key   = tkeys['3']

                if run_neb:
                    it = add_neb(cl, up = 'up2', ise_new = n_set, images = imag,  
                            xr_start = struct_des[p[2]].xr_m_ion_start,
                            xr_final = struct_des[p[2]].xr_m_ion_final, 
                            i_atom_to_move = p[0], i_void_final = p[1], #just needed for name
                            atom_to_insert = el, add_loop_dic = add_loop_args) 
                    calc[t_key_it][pk] = it




                if read:
                    res, _ = res_loop(calc[t_key_it][pk],n_set,range(1, imag+3), up = 'up1', check_job = 0, show = 'fomep', readfiles = 1, analys_type = 'neb', fitplot_args = None)
                    calc[t_key][pk] = calc['_mep']
                elif not run_neb:
                    pos.extend(list((calc[t_key][pk][0])))
                    # mep.extend(list(reversed(calc[t_key][p][1])))
                    mep.extend(list(calc[t_key][pk][1]))




    if substitute:
        for k in substitute:
            mep[k] = substitute[k]

    # print(pos)
    # print(mep)

    if plot:
        name = t_key.replace("'","").replace('(','').replace(')','').replace(' ', '').replace(',','_')
        # print(pos, mep)
        plot_mep(pos, mep, filename = 'figs/path_'+name+suf2, fitplot_args = {'figsize':(8,6), 'hor':1, 'ylim':ylim, 'legend':1, 'first':first, 'last':last}, style_dic = style_dic)

    return




def calc_charged(cl, del_dic, name = None, run = 0, ise = '4uis', it_folder = None ):

    """
    create charged by removing specific sets of atoms provided in del_dic manually starting from 1

    """

    # if not del_dic:
        
    if not it_folder:
        it_folder = name+'/'

    for key in del_dic:
        del_pos = del_dic[key]
        it_new = name+'.c'+str(key)

        if run: 
            st = cl.end.remove_atoms(del_pos, from_one = 1)
            st.name+='c'+str(key)
            # st_halfLi.write_xyz()
            add_loop(it_new, ise, 1, calc_method = 'uniform_scale', inherit_option = 'inherit_xred', scale_region = (-5, 3), input_st = st, it_folder = it_folder)
            # print(st_halfLi.get_space_group_info())
            # add_loop(it_new,'0u',1, input_st = st_halfLi, it_folder = 'Na2FePO4F/chg', override = True) #charge density without relaxation

        else:
            idd = (it_new+'.su', ise, 100)
            res_loop(*idd, show = 'maga', up = 'up1')#list(range(1,8))+[100], analys_type = 'fit_a', show = 'fitfo', up = '1')
            # st = calc[idd].end
            # print(st.get_space_group_info())
            # alpha, beta, gamma = st.get_angles()

            # print(alpha, beta, gamma)  
    return


def optimize(st, name = None, run = 0, ise = '4uis', it_folder = None, fit = 0 ):

    """

    """

    # if not del_dic:
        
    if not it_folder:
        it_folder = 'optimization/'+name.split('.')[0]+'/'

    it_new = name
    if run: 
        add_loop(it_new, ise, 1, calc_method = 'uniform_scale', inherit_option = 'inherit_xred', scale_region = (-5, 3), input_st = st, it_folder = it_folder)

    else:
        idd = (it_new+'.su', ise, 100)
        if fit:
            res_loop(*idd[0:2], list(range(1,8))+[100], analys_type = 'fit_a', show = 'fitfo', up = '1')

        else:
            res_loop(*idd, show = 'maga', up = 'up1')#list(range(1,8))+[100], analys_type = 'fit_a', show = 'fitfo', up = '1')
        # st = calc[idd].end
        # print(st.get_space_group_info())
        # alpha, beta, gamma = st.get_angles()

        # print(alpha, beta, gamma)  
    return


def create_project_from_geofile(filename):
    """
    empty project is added to database
    get name from geofile and create required folder and put file into it  
    Various rules to create project name

    """
    db = header.db
    if 1:
        #simple - just name of file
        basename = os.path.basename(filename)
        projectname = basename.split('.')[0]
        makedir(projectname+'/temp')
        startgeofile = projectname+'/'+ basename


    if projectname not in db:

        shutil.copyfile(filename, startgeofile)
        db[projectname]  = {}
        db[projectname]['startgeofile'] = startgeofile
        db[projectname]['steps'] = [] 
        printlog('Project ', projectname, 'was created', imp = 'y')
    
    else:
        printlog('Error! project already exist')

    return projectname

def get_alkali_ion(st):

    for_diffusion = []
    for el in st.get_elements():
        if el in ['Li', 'Na', 'K', 'Rb', 'Mg']:
            if el not in for_diffusion:
                for_diffusion.append(el)

    el = for_diffusion[0]
    if len(for_diffusion) > 1:
        printlog('Warning! More than one candidate for NEB was found, I use first', el)

    return el




def process_cathode_material(projectname, step = 1, target_x = 0, update = 0, params = None ):
    """
    AI module to process cif file and automatic calculation of standard properties of cathode material 
    
    step 1 - read geo and run simple relaxation

    step 2 - calc barriers, IS

    step 3 - calc barriers, DS

    step 4 - make table with lattice constants for IS and DS

    step 5 - make table with intercalation potential 

    INPUT:
    target_x (float) - required concentration of Na in DS state
    update - allows to rewrite service table
    params (dic)
        show_fit
        run_neb
        up_SC
        up_res
        atom_to_move


    """
    from siman.geo import determine_symmetry_positions, primitive, remove_x
    pn = projectname
    
    m_set = '1u'; n_set  = '1u'; clust = 'cee' # 
    
    p = params
    show_fit  = p.get('show_fit')
    up_scale  = p.get('up_scale')
    up_SC  = p.get('up_SC')
    up_res    = p.get('up_res') or 'up1'
    sc_set    = p.get('sc_set') or '4uis'
    run_neb   = p.get('run_neb')
    # atom_to_move 





    if update or 'res' not in db[pn]:
        db[pn]['res'] = []

    if 'latex' not in db[pn]:
        db[pn]['latex'] = {}


    service_list = db[pn]['res']

    if step == 1:
        ''
        if 1 not in db[pn]['steps']:
            startgeofile = db[pn]['startgeofile'] 
            print('geo file is ', startgeofile)
            st = smart_structure_read(startgeofile)
            st = primitive(st)
            add_loop(pn, m_set, 1, input_st = st, it_folder = pn)
            startgeofile = db[pn]['steps'].append(1) 
        else:
            res_loop(pn, m_set, 1)

    if step in [2, 3]:

        # if 2 not in db[pn]['steps']:
        it = pn
        cl = db[it, m_set, 1]

        el  = get_alkali_ion(cl.end)

        it_ds = it.replace(el, '')
        printlog('Name for DS is', it_ds)

        pd = {'id':cl.id, 'el':el, 'ds':it_ds, 'itfolder':cl.sfolder, 
        'images':5, 'neb_set':n_set, 'main_set':m_set, 'scaling_set':sc_set, 'scale_region':(-3, 5), 'readfiles':1, 'ortho':[10,10,10]}


        pd['atom_to_move'] = p.get('atom_to_move')
        p = p.get('path')


        pd['start_pos'] = p[0] 
        pd['end_pos']   = p[1] 

        add_loop_dic = { 'check_job':1, 'cluster':clust}
        fitplot_args = {'ylim':(-0.02, 1.8)}

        if step == 2:
            style_dic  = {'p':'bo', 'l':'-b', 'label':'IS'}
            a = calc_barriers('normal', up_res = up_res, show_fit = show_fit, up = up_scale, upA = up_SC, upC = 0, param_dic = pd, add_loop_dic = add_loop_dic,
            fitplot_args = fitplot_args, style_dic = style_dic, run_neb = run_neb) 
            
            # print(a)
            if a[0] not in service_list:
                service_list.append(a[0])
            # db[pn]['B'] = [ a[0]['B'] ]

        if step == 3:
            # cl.res()
            style_dic  = {'p':'bo', 'l':'-b', 'label':'DS'}

            st = cl.end

            pos = determine_symmetry_positions(st, el)
            # cl.me()

            if target_x == 0:
                a = calc_barriers('make_ds', el, el, up_res = 'up1', show_fit = show_fit, up = 0, upA = 0, upC = 0, param_dic = pd, add_loop_dic = add_loop_dic,
                fitplot_args = fitplot_args, style_dic = style_dic, run_neb = run_neb) 
                if a[0] not in service_list:
                
                    service_list.append(a[0])

            else:
                x_str = str(target_x).replace('.', '')
                x_vac = 1 - target_x # concentration of vacancies
                name = el+x_str+it_ds

                syms =  remove_x(st, el, info_mode = 1, x = x_vac)

                printlog('The following syms are found', syms, 'I check all of them', imp = 'y')


                # sys.exit()
                for sg in syms:
                    st_rem  =  remove_x(st, el, sg = sg, x = x_vac)

                    id_new = (name+'sg'+str(sg), m_set, 1)
                    add_loop(*id_new, input_st = st_rem, it_folder = cl.sfolder+'/ds', up = 'up1')
                    
                    pd['id'] = id_new
                    a = calc_barriers('normal', el, el, up_res = 'up1', show_fit = show_fit, up = 0, upA = 0, upC = 0, param_dic = pd, add_loop_dic = add_loop_dic,
                    fitplot_args = fitplot_args, style_dic = style_dic, run_neb = run_neb) 
                    info = a[0]
                    info['x'] = target_x
                    if info not in service_list:
                        service_list.append(info)

                    # print (service_list)


    if step == 4:
        ''
        #Lattice constants, and intercalation potentials 
        #the first calculation now is considered as intercalated 
        # IS = 
        from table_functions import table_geometry, table_potentials

        # print('service list is ', service_list)
        """Lattice constants"""

        sts = []
        cll = []
        for a in service_list:
            # print(a)
            cl = db[a['id']]
            try:
                cl.end
            except:
                service_list.remove(a)
                continue
            st = cl.end
            st.x = a['x']
            sts.append(st)
            cll.append(cl)
            if 1:
                """Plot figures"""
                if 'id_sc' in a:
                    db[a['id_sc']].end.write_xyz(jmol = 1)

        table = table_geometry(sts)
        db[pn]['latex']['t1'] = table


        """Intercalation potentials"""
        table = table_potentials(cll)
        db[pn]['latex']['t2'] = table





    if step == 4:
        #create report
        from table_functions import generate_latex_report

        latex_text = ''

        for key in ['t1', 't2']:
            latex_text+=db[pn]['latex'][key] +'\n'

        generate_latex_report(latex_text, filename = 'tex/'+pn+'/'+pn)
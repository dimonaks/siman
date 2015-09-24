# -*- coding: utf-8 -*- 
from header import *
import header
from classes import res_loop, inherit_icalc, add_loop
from functions import write_xyz, replic, write_jmol, image_distance
# from matplotlib import rc
# rc('font',**{'family':'serif'})
# rc('text', usetex=True)
# rc('text.latex',unicode=True)
# rc('text.latex',preamble='\usepackage[utf8]{inputenc}')
# rc('text.latex',preamble='\usepackage[russian]{babel}')

import matplotlib.pyplot as plt
import matplotlib
# calc = header.calc
# conv = header.conv
# varset = header.varset


def checking_mechanical_contribution_Aksyonov2013():
    e0 = -3.1672524186E+03 # puregb
    #Cells without impurity: set 500
    e606 = -3.1672134863E+03 #max force 0.7 eV/A; in VASP for csl71sgCi6Ov.m.83.2 is 0.58 eV/A
    e601 = -3.1672486876E+03 
    # print (e606 - e601)*27.2
    e706 = -3.1672376411E+03
    e701 = -3.1672505836E+03
    # print (e706- e701)*27.2

    print (e601- e0)*27.2
    print (e606- e0)*27.2
    # print (e706- e0)*27.2


    #noimp from DOS folder; set5003
    e606 = -3.1672052644E+03
    e601 = -3.1672404559E+03
    print (e606- e601)*27.2 #0.96 eV         


    #Cells with impurity set 502
    e606 = -3.1730275975E+03 # impurity at gb
    e601 = -3.1730519928E+03 # impurity in grain bulk
    # print (e606 - e601)*27.2
    # e706 = 
    # e701 = 
    return

def calculate_formation_energies():
    """
    Simple snippet to calculate formation energy for TiC(225) from alpha-Ti and graphite or carbon in octahedral void.
    """
    Etic = -18.705 # TiC.83
    Eoct = -10.826 #from hs443C.f - hs443
    Egra = -9.208 #gr221.83
    Eti  = -7.831 #hs221; the same from hs443
    EoctO= -10.520 #from hs443O.f - hs443



    dE1 = Etic - Eti - Eoct # -0.048 eV
    dE2 = (Etic - Eti - Egra)/2. # -0.833 eV
    print dE1, 'eV; relative to octahedral'
    print dE2, 'eV; enthalpy from graphite'
    #Ti2C227
    E = -107.438# Ti2C227.93.1; Ti8C4
    dE1 = (E - 8*Eti - 4*Eoct)/4
    # dE2 = (E - 8*Eti - 4*Egra)/12
    # print dE1

    #Ti2C164
    E = -26.701# Ti2C164.93.1; 
    dE1 = (E - 2*Eti - Eoct)/1
    print dE1, 'Ti2C164'




    #Ti2O227
    E = -102.811# Ti2O227.93.1; Ti8O4
    dE1 = (E - 8*Eti - 4*EoctO)/4
    # dE2 = (E - 8*Eti - 4*Egra)/12
    print dE1, "Ti2O227"

    #Ti2O164
    E = -26.071# Ti2O164.23.1; 
    dE1 = (E - 2*Eti - Eoct)/1
    print dE1, 'Ti2O164'


    return
# calculate_formation_energies()


def seg_energy_vs_voronoi(calc, conv, varset, plot = 0):
    """
    Was used to calculate all segregation energies vs voronoi volume to find some dependence
    plot - if 1: than plot dependence
    """


    aggregate_list = []




    if 0:
        """1. Start processing of calculation specific information"""

        #base = 't111gCOv2'
        base = 't111gCvOvms'
        coseg_list = ['t111gCOi1.4-3', 't111gOCi1.4-3', 't111gCOi2.2-1', 't111gOCi2.2-1', 't111gCOi3.4-1', 't111gOCi3.4-1', \
        't111gCOi4.4-4is','t111gCOi5.2-2is', 't111gCOi6.1-1is', 't111gCOi7.3-3is','t111gCOi8.2-3','t111gOCi8.2-3',  \
        't111gCOi9.4-4ms','t111gCOi10.2-2ms','t111gCOi11.4-3', 't111gOCi11.4-3', \
        't111gCOi12.1-3', 't111gOCi12.1-3', 't111gCOi13.4-2', 't111gOCi13.4-2', 't111gCOi14.2-1', 't111gOCi14.2-1', 't111gCOi15.4-4is','t111gCOi16.2-2is']

        seg_list  = ['t111gCi1Ov', 't111gCi2Ov','t111gCi3Ov','t111gCi4Ov','t111gOi1Cv','t111gOi2Cv','t111gOi3Cv','t111gOi4Cv',] 

        bulk_list = ['t111gCOv2', 't111gCOv5' ]

        list = res_loop(base,'93kp9',range(2,3),calc,conv,varset, 'coseg', (base,) , voronoi = True )
        aggregate_list.append(list)
        # os._exit(0)
        
        for key in seg_list:
            # print key
            # if 'Oi' in key: continue
            list = res_loop(key,'93kp9',range(2,3),calc,conv,varset, 'coseg', (base,) , voronoi = True )
            aggregate_list.append(list)
            
            pass

        # res_loop('t112gCOi6.1-1is','93kp9',range(1,5),calc,conv,varset, 'coseg', ('t112gCvOvms',)  )


        base = 'c1gCvOvms'
        coseg_list = [
        'c1gCOi1.2-1',
        'c1gCOi2.2-1',
        'c1gCOi3.1-2',
        'c1gCOi4.2-1',
        'c1gCOi5.2-1',
        'c1gCOi6.1-2',
        'c1gCOi7.1-1',
        'c1gCOi8.2-2',
        'c1gOCi1.2-1',
        'c1gOCi2.2-1',
        'c1gOCi3.1-2',
        'c1gOCi4.2-1',
        'c1gOCi5.2-1',
        'c1gOCi6.1-2',
        'c1gOCi7.1-1',
        'c1gOCi8.2-2',
        'c1gCOi10.1'
        ]
        bulk_list = ['c1gCOv2',
        'c1gCOv4',
        'c1gCOv7',]

        seg_list = ['c1gCi1Ov', 'c1gCi2Ov', 'c1gOi1Cv', 'c1gOi2Cv', ]

        list = res_loop(base,'93kp7',range(2,3),calc,conv,varset, 'coseg', (base,), voronoi = True  )
        aggregate_list.append(list)

        for key in seg_list:
            # print key
            # if 'Oi' in key: continue
            list = res_loop(key,'93kp7',range(2,3),calc,conv,varset, 'coseg', (base,) , voronoi = True )
            aggregate_list.append(list)

            pass

        base = 't111sgCvOvms'
        seg_list = ['t111sgCi1Ov', 't111sgCi2Ov', 't111sgCi3Ov', 't111sgCi4Ov', 't111sgCi5Ov', 't111sgCi6Ov', 't111sgCi7Ov', 't111sgOi1Cv', 't111sgOi2Cv', 't111sgOi3Cv', 't111sgOi4Cv', 't111sgOi5Cv', 't111sgOi6Cv', 't111sgOi7Cv', ]
        list = res_loop(base,'93kp9',range(2,3),calc,conv,varset, 'coseg', (base,) , voronoi = True )
        aggregate_list.append(list)

        for key in seg_list:
            # print key
            # if 'Oi' in key: continue
            list = res_loop(key,'93kp9',range(2,3),calc,conv,varset, 'coseg', (base,) , voronoi = True )
            aggregate_list.append(list)

            pass

        base = 't21gCvOvms'
        seg_list = ['t21gCi1Ov', 't21gCi2Ov', 't21gCi3Ov', 't21gCi4Ov', 't21gOi1Cv', 't21gOi2Cv', 't21gOi3Cv', 't21gOi4Cv', ]
        list = res_loop(base,'93',range(2,3),calc,conv,varset, 'coseg', (base,), voronoi = True  )
        aggregate_list.append(list)

        for key in seg_list:
            # print key
            # if 'Oi' in key: continue
            list = res_loop(key,'93',range(2,3),calc,conv,varset, 'coseg', (base,) , voronoi = True )
            aggregate_list.append(list)
            write_xyz( replic(calc[(key, '93', 1)].end, (1,2,2)), )
            pass
        base = 'csl71sgCvOvms'
        seg_list = ['csl71sgCi1Ov', 'csl71sgCi2Ov', 'csl71sgCi3Ov', 'csl71sgCi4Ov', 'csl71sgCi5Ov', 'csl71sgCi6Ov', 'csl71sgOi1Cv', 'csl71sgOi2Cv', 'csl71sgOi3Cv', 'csl71sgOi4Cv', 'csl71sgOi5Cv', 'csl71sgOi6Cv', ]
        list = res_loop(base,'93',range(2,3),calc,conv,varset, 'coseg', (base,) , voronoi = True )
        aggregate_list.append(list)

        for key in seg_list:
            # print key
            # if 'Oi' in key: continue
            list = res_loop(key,'93',range(2,3),calc,conv,varset, 'coseg', (base,) , voronoi = True )
            
            aggregate_list.append(list)

            pass
        calc['voronoi'] = aggregate_list

    else:
        aggregate_list = calc['voronoi']
        pass 




    if 1:

        # print res_loop('hs221C','83',101,calc,conv,varset, voronoi = True )
        # # print res_loop('hs221.f','83',1,calc,conv,varset, voronoi = True )
        # # print res_loop('csl71b_r','83kp8',1,calc,conv,varset, voronoi = True )
        # print res_loop('hs221C.f','93',1,calc,conv,varset, voronoi = True ) #9.25 But the sum of Voronoi volume for hexagonal cell is wrong!!!
        # print res_loop('hs332C.f','93',1,calc,conv,varset, voronoi = True ) #9.18 But volume for internal atoms seems to be OK, 
        # print res_loop('hs443C.f','93',1,calc,conv,varset, voronoi = True ) #9.13 But for edge atoms it is even zero.
        # print res_loop('hs554C.f','93',1,calc,conv,varset, voronoi = True ) #9.14 It seems that Voronoi module does not recognize oblique cell geometry.


        """2. Analysis"""

        esegC    = [float(x[1]) for x in aggregate_list if 'Oi' not in x[0]] # C segregation including CvOvms
        vorovolC = [float(x[2]) for x in aggregate_list if 'Oi' not in x[0]]
        # print esegC, 
        esegO    = [float(x[1]) for x in aggregate_list if 'Ci' not in x[0]] # O segregation including CvOvms
        vorovolO = [float(x[4]) for x in aggregate_list if 'Ci' not in x[0]]

        """Plot dependence"""

        name = 'Voronoi volume - segregation energy'
        



        
        set = 'Eng'
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
        
        if set == 'Rus':
            matplotlib.rcParams.update({'font.size': 16})

            xlabel = u'Энергия сегрегации (мэВ)'
            # print xlabel
            ylabel = u'Объем Вороного  ($\AA^3$)'
            label1 = u'Рассчет'
            label2 = u'Аппроксимация'
            title1 = u'Углерод'
            title2 = u'Кислород'
        else:
            xlabel = 'Segregation energy (meV)'
            ylabel = 'Voronoi volume  ($\AA^3$)'
            label1 = 'Calculated'
            label2 = 'Fitted'
            title1 = 'Carbon'
            title2 = 'Oxygen'
        # plt.figure()
        #plt.title(name)


        ax1.plot(esegC, vorovolC, 'o',label = label1)
        
        k, b = np.polyfit(esegC, vorovolC, 1)        

        ax1.plot(esegC, np.asarray(esegC) * k + b, '-r',label = label2)
        ax1.set_title(title1)

        # ax1.legend(loc =4)
        # image_name = 'eseg_vorovol_C'

        # plt.savefig('images/'+str(image_name)+'.png',format='png', dpi = 300)            
        
        # plt.cla()
        
        #plt.title(name)
        # plt.ylabel(ylabel)
        # plt.xlabel(xlabel)
        ax2.plot(esegO, vorovolO, 'o',label = label1)

        k, b = np.polyfit(esegO, vorovolO, 1)        

        ax2.plot(esegO, np.asarray(esegO) * k + b, '-r',label = label2)
        ax2.set_title(title2)
        ax2.set_xlabel(xlabel)
        # ax2.legend(loc =4)

        # image_name = 'eseg_vorovol_O'
        plt.figtext(0.03, 0.5, ylabel, ha='center', va='center', rotation='vertical')

        # fig.set_ylabel(ylabel, va ='center')
        # plt.yaxis.set_label_coords(0, 0.5)
        # plt.xlabel(xlabel)
        plt.tight_layout()
        plt.subplots_adjust(left=0.15, bottom = 0.15)

        image_name = 'eseg_vorovol_CandO'

        plt.savefig('images/'+str(image_name)+'.png',format='png', dpi = 300) 

    return





def make_pictures_for_disser():
    """
    Allows to create semi-automatically beatifull pictures of structures for my dissertation using Jmol.
    """
    path = '/home/aksenov/dim/visualisation/ti-c/disser_2014/'
    path_to_TiC_xyz = path + 'Ti-C_xyz/'

    if 1:
        """Converter to angstroms and cut along z"""
        xyzfilelist = glob.glob(path_to_TiC_xyz+'initial_in_bohr/*.xyz*')
        for xyz in xyzfilelist:
            mul = 0.529177
            with open(xyz,'r') as f:
                lines = f.readlines()
                # lines = l.splitlines()
            name = lines[1]
            natom = lines[0]
            xcart = []
            typ = []
            v = np.zeros((3))
            for l in lines[2:]:
                at = l.split()
                # print at
                v[0] = float(at[1])*mul
                v[1] = float(at[2])*mul
                v[2] = float(at[3])*mul

                if v[2]>12:continue #remove not needed
                typ.append(at[0])
                xcart.append(v.copy())
            name = xyz.split('/')[-1].split('.')[0]+'.xyz'

            out = path_to_TiC_xyz+name
            with open(out,'w') as f:
                f.write(str(len(xcart))+'\n')
                f.write(name+'\n')
                for i, x in enumerate(xcart):
                    f.write(typ[i])
                    f.write(' {0:f} {1:f} {2:f}\n'.format(x[0],x[1],x[2]) )


        # print xcart



    xyzfilelist = glob.glob(path_to_TiC_xyz+'*.xyz*')
    # print geofilelist
    scriptfile = '/home/aksenov/dim/visualisation/ti-c/disser_2014/temp_jmol_sc.jmol'

    # orientation = """moveto /* time, axisAngle */ 1.0 { 50 726 -686 174.13} /* zoom, translation */  100.0 0.0 0.0  /* center, rotationRadius */ {7.9565997 4.6085496 5.5436106} 16.738642 /* navigation center, translation, depth */ {0 0 0} 0 0 0 /* cameraDepth, cameraX, cameraY */  3.0 0.0 0.0;"""
    orientation = "moveto /* time, axisAngle */ 1.0 { -24 -687 -726 177.36} /* zoom, translation */  132.25 0.0 0.0  /* center, rotationRadius */ {7.9565997 4.6085496 5.5436106} 16.74 /* navigation center, translation, depth */ {0 0 0} 0 0 0 /* cameraDepth, cameraX, cameraY */  3.0 0.0 0.0;"

    for xyz in xyzfilelist:#[0:1]:
        name = xyz.split('/')[-1].split('.')[0]+'.png'
        write_jmol(xyz, path+name, scriptfile, atomselection = None, topview = False, orientation = orientation, axis = True, bonds = False)

    return










def remove_imp_segreg(calc, conv, varset):
    """Creates new calculation with removed impurities for determination of """

    calc_list = []

    ver = 2
    set_T1   = '93kp9'
    set_T1_m = '83kp9'
    seg_list_T1  = ['t111gCvOvms', 't111gCi1Ov', 't111gCi2Ov','t111gCi3Ov','t111gCi4Ov','t111gOi1Cv','t111gOi2Cv','t111gOi3Cv','t111gOi4Cv',] 
    seg_list_T1sg = ['t111sgCvOvms', 't111sgCi1Ov', 't111sgCi2Ov', 't111sgCi3Ov', 't111sgCi4Ov', 't111sgCi5Ov', 't111sgCi6Ov', 't111sgCi7Ov', 't111sgOi1Cv', 't111sgOi2Cv', 't111sgOi3Cv', 't111sgOi4Cv', 't111sgOi5Cv', 't111sgOi6Cv', 't111sgOi7Cv', ]

    for el in seg_list_T1:
        calc_list.append( (el, set_T1, ver, 'T1/CO/t111g_segreg_m')   )

    for el in seg_list_T1sg:
        calc_list.append( (el, set_T1, ver, 'T1/CO/t111sg_segreg_m')   )

    set_C1   = '93kp7'
    set_C1_m = '83kp7'
    seg_list_C1 = ['c1gCvOvms', 'c1gCi1Ov', 'c1gCi2Ov', 'c1gOi1Cv', 'c1gOi2Cv', ]

    for el in seg_list_C1:
        calc_list.append( (el, set_C1, ver, 'C1/CO/c1g_segreg_m')   )

    set_T2_CSL7   = '93'
    set_T2_CSL7_m = '83'
    seg_list_T2 = ['t21gCvOvms', 't21gCi1Ov', 't21gCi2Ov', 't21gCi3Ov', 't21gCi4Ov', 't21gOi1Cv', 't21gOi2Cv', 't21gOi3Cv', 't21gOi4Cv', ]
    seg_list_CSL7 = ['csl71sgCvOvms', 'csl71sgCi1Ov', 'csl71sgCi2Ov', 'csl71sgCi3Ov', 'csl71sgCi4Ov', 'csl71sgCi5Ov', 'csl71sgCi6Ov', 'csl71sgOi1Cv', 'csl71sgOi2Cv', 'csl71sgOi3Cv', 'csl71sgOi4Cv', 'csl71sgOi5Cv', 'csl71sgOi6Cv', ]
    
    for el in seg_list_T2:
        calc_list.append( (el, set_T2_CSL7, ver, 'T2/CO/t21g_segreg_m')   )
    for el in seg_list_CSL7:
        calc_list.append( (el, set_T2_CSL7, ver, 'CSL7/CO/csl71sg_segreg_m')   )

    # print len(calc_list)
    for el in calc_list:
        # print "struct_des['{0:s}.m'] = des('{3:s}', 'without imp; made from {0:s}.{1:s}; ver are the same'   )".format(el[0],el[1], el[2], el[3])
        # print el[0:3]

        # inherit_icalc('remove_imp', el[0]+'.m', ver, el[0:3], calc)
        s8 = '8'+el[1][1:]# rename all sets from 93* to 83*
        # add_loop(el[0]+'.m', s8, ver, calc, conv, varset, 'up1' )
        res_loop(el[0]+'.m', s8, ver, calc, conv, varset,  )



    return






def run_dos_segregations(calc, conv, varset):
    list = [
    # ('c1gCi1Ov'    , '93kp7', 2),
    # ('c1gCi2Ov'    , '93kp7', 2),
    # ('c1gOi1Cv'    , '93kp7', 2),
    # ('c1gOi2Cv'    , '93kp7', 2),
    # ('t111gCi2Ov'  , '93kp9', 2),
    # ('t111gCi3Ov'  , '93kp9', 2),
    # ('t111gOi2Cv'  , '93kp9', 2),
    # ('t111gOi3Cv'  , '93kp9', 2),
    # ('t111gCi1Ov', '93kp9', 2),
    # ('t111gCi4Ov', '93ksp9', 2),
    # ('t111gOi1Cv', '93kp9', 2),
    # ('t111gOi4Cv', '93kp9', 2),



    # ( 't111sgCi6Ov', '93kp9', 2),
    # ( 't111sgOi6Cv', '93kp9', 2),
    # ('t111sgCi1Ov', '93kp9', 2),
    # ('t111sgCi2Ov', '93kp9', 2),
    # ('t111sgCi3Ov', '93kp9', 2),
    
    # ('t111sgCi4Ov', '93kp9', 2),
    # ('t111sgCi5Ov', '93kp9', 2),
    # ('t111sgCi7Ov', '93kp9', 2),
    # ('t111sgOi1Cv', '93kp9', 2),
    # ('t111sgOi2Cv', '93kp9', 2),
    # ('t111sgOi3Cv', '93kp9', 2),
    # ('t111sgOi4Cv', '93kp9', 2),
    # ('t111sgOi5Cv', '93kp9', 2),
    # ('t111sgOi7Cv', '93kp9', 2),




    # ('csl71sgCi4Ov', '93'   , 2),
    # ('csl71sgOi4Cv', '93'   , 2),
    # ('t21gCi1Ov','93', 2 ),
    # ('t21gCi4Ov','93', 2 ),
    # ('t21gOi1Cv','93', 2 ),
    # ('t21gOi4Cv','93', 2 ),
    
    # ('t21gCi2Ov','93', 2 ),
    ('t21gCi3Ov','93', 2 ),
    # ('t21gOi2Cv','93', 2 ),
    # ('t21gOi3Cv','93', 2 ),

    # ('c1gCvOvms'    , '93kp7', 2),
    # ( 't111sgCvOvms', '93kp9', 2),
    # ( 't111gCvOvms', '93kp9', 2),
    # ('csl71sgCvOvms', '93'   , 2),
    # ('t21gCvOvms','93', 2 ),
    ]

    for key in list[:]:
        # print "struct_des['{}'] = des('{}', 'relaxed'   )".format(key[0]+'.r', struct_des[key[0]].sfolder)
        # inherit_icalc('full', key[0]+'.r', key[2], (key[0], key[1], key[2]), calc, '', '', 0 )
        # res_loop(key[0]+'.r'    ,'dos',2,calc,conv,varset, 'up2', )
        # res_loop(key[0]+'.m'    ,'dos',2,calc,conv,varset, 'up2', )
        pass




def make_scaled_versions_of_segen(calc, conv, varset):


    list = [
    ('c1gCi1Ov'    , '93kp7', 2  , 'c1gCvOvms' ),
    ('c1gCi2Ov'    , '93kp7', 2  , 'c1gCvOvms' ),
    ('c1gOi2Cv'    , '93kp7', 2  , 'c1gCvOvms' ),
    ('c1gCvOvms'    , '93kp7', 2 , 'c1gCvOvms' ),
    ( 't111sgCi6Ov', '93kp9', 2  , 't111sgCvOvms' ),
    ( 't111sgOi6Cv', '93kp9', 2  , 't111sgCvOvms' ),
    ( 't111sgCvOvms', '93kp9', 2 , 't111sgCvOvms' ),
    ('csl71sgCi4Ov', '93'   , 2  , 'csl71sgCvOvms' ),
    ('csl71sgOi4Cv', '93'   , 2  , 'csl71sgCvOvms' ),
    ('csl71sgCvOvms', '93'   , 2 , 'csl71sgCvOvms' ),
    ('t21gCi4Ov','93', 2         , 't21gCvOvms' ),
    ('t21gOi4Cv','93', 2         , 't21gCvOvms' ),
    ('t21gCvOvms','93', 2        , 't21gCvOvms' ),
    ]

    suffix = ".r2d"; sr2 = 2; sr3 = 1
    # suffix = ".r3d"; sr2 = 1; sr3 = 2
    verlist = [1,2,3,4]
    for id in list[:]:
        for v in verlist[:]:
            cl = calc[(id[0], id[1], v)]
            # cl_doubled = copy.deepcopy(cl)

            # cl_doubled.init = replic(cl.end, (1,sr2,sr3), only_atoms= "only_matrix")
            new_name = id[0]+suffix
            # cl_doubled.des = "The relaxed "+cl.name+" was doubled along "+suffix+" only for matrix; the sizes are the same"
            # cl_doubled.path["input_geo"] = 'geo/scaled/'+new_name+"/"+new_name+"."+str(v)+".geo"
            # cl_doubled.write_geometry(geotype = "init")

        # print "struct_des['{}'] = des('{}', 'scaled only matrix; the sizes are doubled without any manipulations'   )".format(new_name, "scaled")
        

        # add_loop(new_name,'93',verlist ,calc,conv,varset, 'up1', savefile = 'only_outcar', inherit_option = 'inherit_xred' )
        # print calc[]
        res_loop(new_name,'93',verlist ,calc,conv,varset, 'e_seg', (id[3]+suffix,'93',) )

    return



def thesis_presentation_pictures():
    #Pictures for presentation
    # write_xyz(calc[('hs332C.f','93', 1)].end, repeat = 2, analysis = 'imp_surrounding',nnumber = 6,topview = 0, boundbox = None,
    # specialcommand = 'select all\nwireframe 0.1\nselect Ti*\ncpk 200\nselect C*\ncpk 230\ncolor red\n',
    # orientation = 'moveto /* time, axisAngle */ 1.0 { 233 510 828 162.65} /* zoom, translation */  75.61 0.0 0.0  /* center, rotationRadius */ {3.396805 2.94223 3.50884} 3.7937615 /* navigation center, translation, depth */ {0 0 0} 0 0 0 /* cameraDepth, cameraX, cameraY */  3.0 0.0 0.0;' )
    write_xyz(calc[('hs443C.f','93', 1)].init, repeat = 2, topview = 0 )
    return
def testing_average_deviation():
    l = [2.08, 2.09, 2.08, 2.06, 2.09, 2.11, ] #for C1(2) end from jmol

    av = sum(l)/6
    print av
    print sum([ abs(av - d) for d in l ])/6 * 1000
    print max([ abs(av - d) for d in l ]) * 1000
    return

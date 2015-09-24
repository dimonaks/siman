# -*- coding: utf-8 -*-
from header import *
import header
from operator import itemgetter
from classes import res_loop , add_loop
# from pairs import 
from functions import image_distance, local_surrounding, write_xyz
from chargeden.functions import chg_at_point, cal_chg_diff
from dos.functions import plot_dos
import matplotlib.pyplot as plt
plt.rcParams['mathtext.fontset'] = "stix"
from picture_functions import fit_and_plot, plot_bar
from scipy.interpolate import interp1d
from scipy import interpolate
from ase.utils.eos import EquationOfState
import inspect
def latex_table(table, caption, label, header = None, fullpage = '', filename = None, writetype = 'w', header0 = None, size = None ):
    """
    If header is not provided, table[0] is used as a header

    header0 - additional header0 befor main header for complex tables
    """
    def myprint(string):
        if filename:
            f.write(string+"\n")
            print string
        else:
            print string


    if filename:
        path = "/home/dim/gdrive/Наука/paper5/tab/"
        f = open(path+filename, writetype)



    n = len(table[0].split('&'))-2
    
    myprint('\\begin{table'+fullpage+'}')
    myprint('\\center')
    if size: myprint('\\'+size)
        


    myprint('\\caption{'+caption+'}')
    myprint('\\label{'+label+'}')

    myprint('\\begin{tabular}{l'+ n*'c'+'r}')
    myprint('\\hline')

    if header0:
        myprint(header0+'\\\\')
        myprint('\\hline')

    if header:
        myprint(header+'\\\\')
        tabbeg = 0
    else:
        myprint(table[0]+'\\\\')
        tabbeg = 1



    myprint('\\hline')
    for r in table[tabbeg:] :
        if '&-' in r:
            r = r.replace('-','--')
        else:
            r = r.replace(' -','--') #to save beautiful columns 

        if 'hline' in r: 
            myprint(r)
        else:
            myprint(r + '\\\\')

    myprint('\\hline')
    myprint('\\end{tabular}')
    myprint('\\end{table'+fullpage+'}')

    if filename:
        f.close()
    return


def make_table_puregb(calc, conv, varset):


    if 0: #can be updated if needed
        list_of_boundaries = [
        ('c1g',     '929',    range(1,11),   ('c1b_r','830',1)      , 'C1',    ),
        ('t111g',   '9292',   range(1,11),   ('t111b_r', '8302', 1) , 'T1 m',  ),
        ('t111sg',  '93kp9',  range(1,5),    ('t111b_r', '8302', 1) , 'T1 g',  ),
        ('t21g',    '93',     range(1,6),    ('t21b_r','83',1)      , 'T2',    ),
        ('csl71g',  '93',     range(3,7),    ('csl71b_r','83kp8',1) , 'S7',     ),
        ('csl71sg5','93',     range(1,6),    ('csl71b_r','83kp8',1) , 'S7.5',   ),
        ('csl71sg10','93',    range(1,6),    ('csl71b_r','83kp8',1) , 'S7.10',   ),
        ('csl71sg15','93',    range(1,6),    ('csl71b_r','83kp8',1) , 'S7.15',   ),
        ]

        outlist = []
        for c in list_of_boundaries:
            outlist.append( [c[4]]+res_loop(c[0], c[1], c[2], calc, conv, varset, 'gbep', c[3], readfiles = 1 ) )
            # print out
            pass
        calc['puregb_table'] = outlist

    else:
        outlist = calc['puregb_table']


    print outlist #contains gb energies and volumes at -500 MPa
    s1 = 8
    s2 = 11
    print '\n'+'GB'.ljust(s2),
    for l in outlist:
        print "&"+l[0].center(s1), 
    print '\\\\\n'+'$\gamma$'.ljust(s2),
    for l in outlist:
        print "&"+l[2].center(s1), 
    print '\\\\\n'+'$\delta$'.ljust(s2),
    for l in outlist:
        print "&"+l[3].center(s1), 
    print '\\\\'
    #res_loop('c1g','929',range(1,10), calc, conv, varset, 'gbep',('c1b_r','830',1) )
        
        # res_loop('t111g','9292',range(1,11), calc, conv, varset, 'gbep',('t111b_r','8302',1) )
    #res_loop('t21g','93',range(1,6), calc, conv, varset,'gbep',('t21b_r','83',1))
    #res_loop('t21g','935',range(11,15), calc, conv, varset,'gbep',('t21b_r','835',1))
    #res_loop('csl71g','93',range(3,7), calc, conv, varset,'gbep',('csl71b_r','83kp8',1))
    return

def table_imp(calc, conv, varset):


    if 1: #can be updated if needed
        E_O2_in_box = -4.926 #per one oxygen atom in dimer; from res_loop('OO','9dim_spin',1, calc, conv, varset, 'dimer')
        list_of_structures = [
        ('hs221C.f',   '93',   1,    ('hs221.f', '83', 1)   , 'HCP' , ('gr221', '83', 1)  ),
        ('hs332C.f',   '93',   1,    ('hs332.f', '83', 1)   , 'HCP' , ('gr221', '83', 1)  ),
        ('hs443C.f',  '93',   1,     ('hs443.f', '83', 1)   , 'HCP' , ('gr221', '83', 1)  ),
        # ('hs443C.f',  '93f',   1,     ('hs443.f', '83', 1)   , 'H.f' , ('gr221', '83', 1)  ), #fine relaxation
        # ('hs443C',  '93',   100,     ('hs443.f', '83', 1)   , 'H.ac(Ti)' , ('gr221', '83', 1)  ), #lattice constants as for pure Ti
        ('hs554C.f',   '93',   1,    ('hs554.f', '83', 1)     , 'HCP' , ('gr221', '83', 1)  ),
        ('t111bC.f',  '93kp9',   1,  ('t111b_r', '8302', 1) , 'T1', ('gr221', '83', 1)  ),
        ('t21bC.f',   '93',   1,     ('t21b_r', '83', 1)    , 'T2', ('gr221', '83', 1)  ),
        ('c1bC.f',  '93kp7',   1,    ('c1b_r', '830', 1)    , 'C1', ('gr221', '83', 1)  ),
        ('hs221O.f',   '93',   1,    ('hs221.f', '83', 1)   , 'HCP' , E_O2_in_box  ), 
        ('hs332O.f',   '93',   1,    ('hs332.f', '83', 1)   , 'HCP' , E_O2_in_box  ),
        ('hs443O.f',  '93',   1,     ('hs443.f', '83', 1)   , 'HCP' , E_O2_in_box   ),
        ('hs554O.f',  '93',   1,     ('hs554.f', '83', 1)   , 'HCP' , E_O2_in_box   ),       
        ('t111bO.f',  '93kp9',   1,  ('t111b_r', '8302', 1) , 'T1', E_O2_in_box  ),
        ('t21bO.f',   '93',   1,     ('t21b_r', '83', 1)      , 'T2', E_O2_in_box  ),   # These last two was fitted; Initially was not fitted, the lattice constants was used from t111bO.f; The difference after fitting is really small.
        ('c1bO.f',  '93kp7',   1,    ('c1b_r', '830', 1)      , 'C1', E_O2_in_box  ),   # v1 and v20 are the same. apparently It was made  for convenience of use.
        #
        # ('t111gCvOvms',  '93kp9', 2, ('t111g', '9292', 1)   , 'T1mCO',     -10.768 -10.513 ), #-21.346),#   relative to T1bC and T1bO # or H96C and H96O     
        # ('t111sgCvOvms',  '93kp9',2, ('t111sg', '93kp9', 2) , 'T1gCO',     -10.768 -10.513 ), #-21.346),#
        # ('t21gCvOvms',   '93',    2, ('t21g', '93', 2)      , 'T2CO' ,     -10.804 -10.556 ), #-21.346),#   
        # ('c1gCvOvms',  '93kp7',   2, ('c1g', '929', 2)      , 'C1CO' ,     -10.786 -10.590 ), #-21.346),#  

        ]

        outlist = []
        for c in list_of_structures:
            st = calc[(c[0], c[1], c[2])]
            nTi = st.nznucl[0]
            nX  = st.nznucl[1]
            if st.znucl[1] == 6: X = 'C'
            elif st.znucl[1] == 8: X = 'O'
            name = 'Ti$_{'+str(nTi)+'}$'+X+' ('+c[4]+')'

            outlist.append( [name]+res_loop(c[0], c[1], c[2], calc, conv, varset, 'e_imp', c[3], r_id = c[5], readfiles = True ) )
            # print out
            pass
        calc['imp_table'] = outlist

    else:
        outlist = calc['imp_table']


    # print outlist
    s1 = 8
    s2 = 22


    table = []
    caption = """Solution energy $E_{sol}$ (eV), absolute  change of volume $\Delta V$ and  relative volume $\delta V$ ({\AA $^3$ } and \%) of impurity 
    and changes of lattice constants $\delta a$, $\delta c$, $\delta (c/a)$ (\%)
depending on the sizes $L_1$, $L_2$ and $L_3$ (\AA) of the simulation cell."""

    header =  'Cell'.ljust(22)+'&'+'$E_{sol}$'.center(8)+'&'+'$\Delta V$'.center(8)+'&'+'$\delta V$'.center(8)+'&'+'$L_1$ & $L_2$ & $L_3$'.center(22)+\
                "&"+'$\delta a$'.center(8)+"&"+'$\delta c$'.center(8)+"&"+'$\delta c/a$'.center(8)+'\\\\'


    for l in outlist:
        # table.append(l[0].ljust(22)+"&"+l[2].center(8)+"&"+l[21].center(8)+"&"+l[3].center(8)+"&"+l[4].center(8)+"&"+l[8].replace(';',' & ').center(22)+\
        # "&"+l[18].center(8)+"&"+l[19].center(8)+"&"+l[20].center(8))

        table.append(l[0].ljust(22)+"&"+l[21].center(8)+"&"+l[3].center(8)+"&"+l[4].center(8)+"&"+l[8].replace(';',' & ').center(22)+\
        "&"+l[18].center(8)+"&"+l[19].center(8)+"&"+l[20].center(8))

    latex_table(table, caption, label = 'tab:e_imp', header = header, filename = 'tab_e_imp.tex', fullpage = '*', size = 'small')





    return





def segen_analysis():
    """
    To study the physical reasons of negative segregation energies



    Смысл методики заключается в следующем.

    1. Наиболее ярким фактором, влияющем на полные энергии, являются значения ПЭС d-орбиталей титана на уровне Ферми.
    2. При сегрегации наблюдается два эффекта. 
        1. Изменение энергии при образовании связей примеси с окружающими атомами -  I  = x(with imp) - x(without imp), где x - энергия, ПЭС на уровне Ферми и т.д. - 
        определеяется при замороженных положениях атомов - соответствует химическому вкладу.
        2. Изменение энергии при релаксации положений окружающих атомов  - R = x(relaxed) - x(ideal).
        Оба эффекта могут происходить как в объеме, так и на границе I(v), R(v),  I(gb), R(gb) , 
        C этими эффектами связано изменение ПЭС на уровне Ферми. Тогда можно рассчитать
        два значения: I(gb) - I(v)  и  R(gb) - R(v) 






    """ 
    calc = header.calc
    # it_list = [
    # ('c1gCi1Ov'    ,'c1gCvOvms'),
    # ('c1gCi2Ov'    ,'c1gCvOvms'),
    # ('c1gOi1Cv'    ,'c1gCvOvms'),
    # ('c1gOi2Cv'    ,'c1gCvOvms'),
    # ('t111gCi2Ov'  ,'t111gCvOvms'),
    # ('t111gCi3Ov'  ,'t111gCvOvms'),
    # ('t111gOi2Cv'  ,'t111gCvOvms'),
    # ('t111gOi3Cv'  ,'t111gCvOvms'),
    # ( 't111sgCi6Ov','t111sgCvOvms'),
    # ( 't111sgOi6Cv','t111sgCvOvms'),
    # ('t21gCi1Ov', 't21gCvOvms'),
    # ('t21gCi4Ov', 't21gCvOvms'),
    # ('t21gOi1Cv', 't21gCvOvms'),
    # ('t21gOi4Cv', 't21gCvOvms'),
    # # ('csl71sgCi4Ov','csl71sgCvOvms'),
    # # ('csl71sgOi4Cv','csl71sgCvOvms'),
    # ]

    it_list = [
    ('c1gCi1Ov'    ,'c1gCvOvms'),
    ('c1gCi2Ov'    ,'c1gCvOvms'),
    ('c1gOi1Cv'    ,'c1gCvOvms'),
    ('c1gOi2Cv'    ,'c1gCvOvms'),

    ('t111gCi2Ov'  ,'t111gCvOvms'),
    ('t111gCi3Ov'  ,'t111gCvOvms'),
    ('t111gOi2Cv'  ,'t111gCvOvms'),
    ('t111gOi3Cv'  ,'t111gCvOvms'),

    ('t111gCi1Ov',  't111gCvOvms'),

    ('t111gCi4Ov',  't111gCvOvms'),
    ('t111gOi1Cv',  't111gCvOvms'),
    ('t111gOi4Cv',  't111gCvOvms'),


    ('t111sgCi6Ov', 't111sgCvOvms'),
    ('t111sgOi6Cv', 't111sgCvOvms'),
    
    ('t111sgCi1Ov', 't111sgCvOvms'),
    ('t111sgCi2Ov', 't111sgCvOvms'),
    ('t111sgCi3Ov', 't111sgCvOvms'),
    
    ('t111sgCi4Ov', 't111sgCvOvms'),
    ('t111sgCi5Ov', 't111sgCvOvms'),
    ('t111sgCi7Ov', 't111sgCvOvms'),
    ('t111sgOi1Cv', 't111sgCvOvms'),
    ('t111sgOi2Cv', 't111sgCvOvms'),
    ('t111sgOi3Cv', 't111sgCvOvms'),
    ('t111sgOi4Cv', 't111sgCvOvms'),
    ('t111sgOi5Cv', 't111sgCvOvms'),
    ('t111sgOi7Cv', 't111sgCvOvms'),
    ('t21gCi1Ov',   't21gCvOvms' ),
    ('t21gCi4Ov',   't21gCvOvms' ),
    ('t21gOi1Cv',   't21gCvOvms' ),
    ('t21gOi4Cv',   't21gCvOvms' ),
    ('t21gCi2Ov',   't21gCvOvms' ),
    
    # ('t21gCi3Ov',   't21gCvOvms' ),
    
    ('t21gOi2Cv',   't21gCvOvms' ),
    ('t21gOi3Cv',   't21gCvOvms' ),
    ('csl71sgCi4Ov','csl71sgCvOvms'),
    ('csl71sgOi4Cv','csl71sgCvOvms'),
    ]






    ise = "dos"; v = 2
    results = {}
    for key in 'Cgb', 'Cbk', 'Ogb', 'Obk':
        results[key] = []
    # imp = 'O'
    # if imp == 'O': impatom = -1 # around oxygen
    # if imp == 'C': impatom = -2 # around carbon
    # for cl_r, cl_m in [  (calc[(it+'.r', ise, v)], calc[(it+'.m', ise, v)]) for it in it_list if 't111sg' in it   ]:
    # for cl_r, cl_m in [  (calc[(it+'.r', ise, v)], calc[(it+'.m', ise, v)]) for it in it_list if 't111g' in it   ]:
    # for cl_r, cl_m in [  (calc[(it+'.r', ise, v)], calc[(it+'.m', ise, v)]) for it in it_list if 'c1g' in it   ]:
    if 0:
        added_bulkC,added_bulkO  = [], []
        for cl_r, cl_m, clv_r, clv_m in [  (calc[(iti+'.r', ise, v)], calc[(iti+'.m', ise, v)], calc[(itv+'.r', ise, v)], calc[(itv+'.m', ise, v)]  ) for iti, itv in it_list]:# if 'c1gCi1Ov' in iti]:
            neighbors = 10

            # print cl_r.name

            """0. Geometry of void"""
            # write_xyz(cl_r, analysis = 'imp_surrounding', nnumber = 8   )

            """1. Charge differences"""
            # print cal_chg_diff( cl_r.dir,  cl_r.get_chg_file(),  cl_m.get_chg_file() )
            # if 'c1gCi1' in cl_r.name or 'c1gOi1' in cl_r.name: neighbors = 8
            # print clv_r.name; cl_r.read_results()
            # continue
            """2. DOS """
            Esegch = (cl_r.energy_sigma0 - clv_r.energy_sigma0) - (cl_m.energy_sigma0 - clv_m.energy_sigma0) #chemical contribution of segregation energy
            Eseg = (cl_r.energy_sigma0 - clv_r.energy_sigma0)
            listsegen = [Eseg, Esegch]

            plot = True



            if 1: #with impurity vs without  I

                if 'Ci' in cl_r.name:
                    results['Cgb'].append( plot_dos(cl_r, cl_m,   'partial', impatom = -2, neighbors = neighbors, plot = plot, 
                        # xlim = (-1,1), ylim = (None, 20),
                         )+listsegen )
                    try: #produce corresponding  list with calculation for impurity in the bulk; to save time if plot_dos was already called for this bulk cell it is copied
                        i = added_bulkC.index(clv_r.name)
                        results['Cbk'].append(results['Cbk'][i])
                        added_bulkC.append("copy_of_"+clv_r.name)
                    except ValueError:
                        results['Cbk'].append( plot_dos(clv_r, clv_m, 'partial', impatom = -2, neighbors = neighbors, plot = plot, 
                            # xlim = (-1,1), ylim = (None, 20)  
                            ) )
                        added_bulkC.append(clv_r.name)
                
                if 'Oi' in cl_r.name:
                    results['Ogb'].append( plot_dos(cl_r, cl_m,   'partial', impatom = -1, neighbors = neighbors, plot = plot, 
                        # xlim = (-1,1), ylim = (None, 20)  
                        )+listsegen )
                    try:
                        i = added_bulkO.index(clv_r.name)
                        results['Obk'].append(results['Obk'][i])
                        added_bulkO.append("copy_of_"+clv_r.name)
                    except ValueError:
                        results['Obk'].append( plot_dos(clv_r, clv_m, 'partial', impatom = -1, neighbors = neighbors, plot = plot, 
                            # xlim = (-1,1), ylim = (None, 20) 
                            ) )
                        added_bulkO.append(clv_r.name)
                
                calc['segen_analysis'] = results

            else: #without impurity: after relax vs before relax,  R

                if 'Ci' in cl_r.name:
                    results['Cgb'].append( plot_dos(cl_m, clv_m,   'partial_both_empty', impatom = cl_r.end.xcart[-2], neighbors = neighbors,     plot = plot,   ) + listsegen )
                    try: #produce corresponding  list with calculation for impurity in the bulk
                        i = added_bulkC.index(clv_r.name)
                        results['Cbk'].append(results['Cbk'][i])
                        added_bulkC.append("copy_of_"+clv_r.name)
                    except ValueError:
                        results['Cbk'].append( plot_dos(clv_m, cl_m, 'partial_both_empty', impatom = clv_r.end.xcart[-2], neighbors = neighbors, plot = plot,   )   )
                        added_bulkC.append(clv_r.name)
                
                if 'Oi' in cl_r.name:
                    results['Ogb'].append( plot_dos(cl_m, clv_m,   'partial_both_empty', impatom = cl_r.end.xcart[-1], neighbors = neighbors,     plot = plot,   ) + listsegen )
                    try:
                        i = added_bulkO.index(clv_r.name)
                        results['Obk'].append(results['Obk'][i])
                        added_bulkO.append("copy_of_"+clv_r.name)
                    except ValueError:
                        results['Obk'].append( plot_dos(clv_m, cl_m, 'partial_both_empty', impatom = clv_r.end.xcart[-1], neighbors = neighbors, plot = plot,   )   )
                        added_bulkO.append(clv_r.name)
                
                calc['segen_analysis2'] = results





    else:

        if 1:
            print "With impurity vs without impurity comparison"
            results  = calc['segen_analysis']  # with impurity vs without
        else:
            print "After relax vs before relax comparison"

            results  = calc['segen_analysis2'] # after relax vs before relax







    """DOS analysis"""
    #e_at_Ef_shift - shift of dos at Fermy energy due to adding of impurity
    #d_shift - Shift of Ti d center of gravity after adding impurity
    #pd_drms - C p- Ti d hybridization,  the higher the number the higher hybridization

    if 0:
        for imp in 'C', 'O':
        #     gb, bulk = [], []
            print "Cell: around ", imp, '| E_seg(meV) | E_seg_ch(meV), Delta (gb-bk) for e_at_Ef_shift, d_shift, pd_drms      |:'
            for rgb, rbk in zip(results[imp+'gb'], results[imp+'bk']):
                print "{:15s} & {:5.0f} & {:5.0f} & {:5.2f} & {:5.2f} & {:5.2f}".format(rgb[0].split('.')[0], rgb[-2]*1000, rgb[-1]*1000, rgb[1] - rbk[1], rgb[2] - rbk[2], rbk[3] - rgb[3]  )

            # for rgb, rbk, rgb2, rbk2 in zip(results[imp+'gb'], results[imp+'bk'], results2[imp+'gb'], results2[imp+'bk']):
            #     print "{:15s} & {:5.0f} & {:5.0f} & {:5.2f} & {:5.2f} & {:5.2f}".format(rgb[0].split('.')[0], rgb[-2]*1000, rgb[-1]*1000, rgb[1] - rbk[1] + rgb2[1] - rbk2[1], rgb[2] - rbk[2] + rgb2[2] - rbk2[2], rbk[3] - rgb[3]  )





        #         if imp+'i' in r[0]:
        #             gb.append(  "{:25s} {:5.2f} & {:5.2f}".format(r[0].split('.')[0], r[1], r[2])  )
        #         if imp+'v' in r[0]: 
        #             bulk.append("{:5.2f} & {:5.2f} {:25s} ".format(r[1], r[2], r[0].split('.')[0])   )
        #     # print 
        #     for g,b in zip(gb, bulk):
        #         print g, b


        imp = 'C'; C = [ [rgb[0].split('.')[0], rgb[-2]*1000, rgb[-1]*1000, rgb[1] - rbk[1], rgb[2] - rbk[2], rbk[3] - rgb[3] ] for  rgb, rbk in zip(results[imp+'gb'], results[imp+'bk']) ]
        imp = 'O'; O = [ [rgb[0].split('.')[0], rgb[-2]*1000, rgb[-1]*1000, rgb[1] - rbk[1], rgb[2] - rbk[2], rbk[3] - rgb[3] ] for  rgb, rbk in zip(results[imp+'gb'], results[imp+'bk']) ]
        C = zip(*C)
        O = zip(*O)
        print C[1]
        
        Ofun = [red + 200*shift for red, shift in zip(O[2], O[3])]
        
        fit_and_plot(image_name = "eseg_delta_dos_C",
            xlabel = "Change at Fermi (DOS)", ylabel = "Segregation energy (meV)",
            # title = "Segregation energy vs change",
            legend = 2,
            # xlim = (-5,1.5), ylim = (None, 2),#0.4
            # xlim = (-0.5,+0.5), 
            ylim = (-350, 220),#0.4
            carbon   = (C[2], C[1], 'ok'),
            # oxygen   = (O[2], O[1], 'or'),
            # oxygen   = (Ofun, O[1], 'ob'),
            )

        fit_and_plot(image_name = "eseg_shift_C",
            xlabel = "Shift of DOS", ylabel = "Segregation energy (meV)",
            legend = 2,
            # xlim = (-5,1.5), ylim = (None, 2),#0.4
            # xlim = (-0.5,+0.5), 
            # ylim = (-350, 220),#0.4
            carbon   = (C[3], C[1], 'ok'),
            ) 
        """ Зависимость полной энергии сегрегации от сдвига центра тяжести ПЭС d орбиталей титана. Видно два области, которые можно
        подогнать прямыми. Для каждой из областей наблюдается логичная тенденция, согласно которой, чем сильнее плотность сдвигается
        в область более низких энергий при добавлении примеси, тем меньше энергия сегрегации. Остается не ясным, почему первая область
        соответствует более сильным отрицательным сдвигам и при этом данные положения имеют положительные энергии сегрегации. Эта область
        состоит из сильнодеформированных междоузлий, поэтому можно предположить, что несмотря на попытки примеси стабилизировать эти
        междоузлия, путем сдвига плотности и уменьшения её значений на уровне Ферми этого оказывается недостаточным, а с другой стороны,
        такая оптимизация приводит к деоптимизации в других местах. Чтобы увидеть эту деоптимизацию возможно нужно провести разделение
        на связывающие-антисвязывающие орбитали, а также более внимательно изучить изменение ПЭС для более удаленных атомов. К примеру
        для T2(1) сдвиг орбиталей компенсируется сильным увеличением ПЭС на уровне ферми, что приводит к положительной энергии.
        Также и для T1(2). 

        Разительное выбивается из тренда лишь T1g(1) - при большой положительной энергии сегрегации для него наблюдается
        как сдвиг центра тяжести ПЭС в область низких энергий, так и существенное снижение ПЭС на уровне Ферми.
        Однако построение ПЭС на уровне Ферми для пустых ячеек с замороженными атомами и вычисление e_at_Ef_shift
        показывает, что для этого положения оптимизация границы после добавления примеси существенно увеличивает
        значения ПЭС на уровне Ферми, что нивелирует все выигрышы и делает сегрегацию не выгодной. 

        Остается
        не понятным, почему два использованных метода определения e_at_Ef_shift, которые должны коррелировать с химическим
        и механическими вкладами энергию сегрегации иногда выбиваются из тренда. Возможно к этому приводит то, что
        подобное разделение является условным и в химическом вкладе все равно присутствует упругий, а в упругом вкладе химический.
        Нужно подумать над этим. Ответ: Все таки судить лишь по одному значению ПЭС на уровне Ферми не совсем корректно.
        Так для t21gCi1Ov разница между e_at_Ef_shift для R эффекта является отрицательной, что указывает на отрицательный
        механический вклад в энергию сегрегации (Но он сильно положительный). Однако, одного взгляда на ПЭС достаточно, чтобы
        увидеть существенной увеличение всей псевдощели ниже уровня Ферми, что сильно не выгодно. 

        """


        # plot_dos(calc[('t111sgCi1Ov.m', 'dos', 2)], calc[('t111sgCvOvms.m', 'dos', 2)],   'partial_both_empty', impatom = calc[('t111sgCi1Ov.r', 'dos', 2)].end.xcart[-2], neighbors = 10, plot = 1,
        # # xlim = (-1,1), ylim = (None, 20),
        #   )
        # plot_dos(calc[('t111sgCvOvms.m', 'dos', 2)], calc[('t111sgCi1Ov.m', 'dos', 2)],   'partial_both_empty', impatom = calc[('t111sgCvOvms.r', 'dos', 2)].end.xcart[-2], neighbors = 10, plot = 1,
        # #xlim = (-1,1), ylim = (None, 20),
        #   )
        # plot_dos(calc[('t111sgCi1Ov.m', 'dos', 2)], calc[('t111sgCvOvms.m', 'dos', 2)],   'partial_both_empty', impatom = calc[('t111sgCi1Ov.r', 'dos', 2)].end.xcart[-2], neighbors = 10, plot = 1,
        # # xlim = (-1,1), ylim = (None, 20),
        #   )
    return








def table_pressure_energy(calc, conv, varset):
    """
    To study the influence of pressure on energy of impurity in octahedral site

    it_p - pure
    it_i - with impurity
    it_e - empty, removed impurity, frozen
    """
    table = []
    implist = {}
    it_p = 'hs443';    ise_p = '83';
    it_i = 'hs443C';   ise_i = '93'; 
    it_e = 'hs443C.m'; ise_e = '83'
    
    implist['C'] = [
    (it_i, ise_i, 1  ,  it_p,     ise_p, 1,     it_e, ise_e, 1  ),
    (it_i, ise_i, 6  ,  it_p,     ise_p, 6,     it_e, ise_e, 6  ),
    ('hs443C.f','93',1, 'hs443.f','83',1,       'hs443C.fm','83',1  ),
    (it_i, ise_i, 11 ,  it_p,     ise_p, 11,    it_e, ise_e, 11 ),
    (it_i, ise_i, 16 ,  it_p,     ise_p, 16,    it_e, ise_e, 16 ),
    (it_i, ise_i, 38 ,  it_p,     '83', 38,  it_e, ise_e, 38 ),
    (it_i, ise_i, 39 ,  it_p,     '83', 39,  it_e, ise_e, 39 ),
    (it_i, ise_i, 40 ,  it_p,     '83', 40,  it_e, ise_e, 40 ),
    (it_i, ise_i, 41 ,  it_p,     '83', 41,  it_e, ise_e, 41 ),
    (it_i, '95', 50 ,  it_p,     '85', 50,  it_e, '85', 50 ),
    (it_i, '95', 511 ,  it_p,    '85', 511, it_e, '85', 511 ),
    (it_i, '95', 52 ,  it_p,     '85', 52,  it_e, '85', 52 ),
    (it_i, '95', 53 ,  it_p,     '85', 53,  it_e, '85', 53 ),
    ]
    


    it_p = 'hs443';    ise_p = '83';
    it_i = 'hs443O';   ise_i = '93'; 
    it_e = 'hs443O.m'; ise_e = '83'
    
    implist['O'] = [
    (it_i, ise_i, 1  ,  it_p,     ise_p, 1,     it_e, ise_e, 1  ),
    (it_i, ise_i, 6  ,  it_p,     ise_p, 6,     it_e, ise_e, 6  ),
    ('hs443O.f','93',1, 'hs443.f','83',1,       'hs443C.fm','83',1  ),
    (it_i, ise_i, 11 ,  it_p,     ise_p, 11,    it_e, ise_e, 11 ),
    (it_i, ise_i, 16 ,  it_p,     ise_p, 16,    it_e, ise_e, 16 ),
    (it_i, ise_i, 38 ,  it_p,     '83', 38,  it_e, ise_e, 38 ),
    (it_i, ise_i, 39 ,  it_p,     '83', 39,  it_e, ise_e, 39 ),
    (it_i, ise_i, 40 ,  it_p,     '83', 40,  it_e, ise_e, 40 ),
    (it_i, ise_i, 41 ,  it_p,     '83', 41,  it_e, ise_e, 41 ),
    (it_i, '95', 50 ,  it_p,     '85', 50,  it_e, '85', 50 ),
    (it_i, '95', 51 ,  it_p,     '85', 51, it_e,  '85', 51 ),
    (it_i, '95', 52 ,  it_p,     '85', 52,  it_e, '85', 52 ),
    (it_i, '95', 53 ,  it_p,     '85', 53,  it_e, '85', 53 ),
    ]




    
    plt.figure()

    for imp in 'CO': # make table

        if imp == 'C':
            e_at = -9.208 #graphite
            Vsol = 7 * 1e-30 #m^3
            color = 'ko-'
        elif imp == 'O':
            e_at = -4.926 #O2
            Vsol = 4 * 1e-30 #m^3
            color = 'ro-'



        a, c = [], []
        volumes_i = []
        volumes_p = []
        etot_i, etot_e = [], []
        etot_p = []
        dElist = []
        pres_i = []
        pres_p = []
        voro_i, daver_i = [], []
        names = []
        esol_m = []




        temp =[0,]; i = 0  

        for it_i, ise_i, v_i, it_p, ise_p, v_p, it_e, ise_e, v_e   in implist[imp]:



            try:
                id_i = (it_i,ise_i,v_i) # with impurity
                id_p = (it_p,ise_p,v_p) # pure Ti
                cl_i = calc[id_i]
                cl_p = calc[id_p]
            except:
                print 'Something from ', id_i, id_p, ' not found, continue'
                continue

            print id_i, id_p
            try:
                pres_i.append(calc[id_i].extpress)
                pres_p.append(calc[id_p].extpress)
            except:
                continue


            try:
                vvol = calc[id_i].vorovol[0]
            except:
                pass
                vvol = 0
            else:
                vvol = 0

            voro_i.append(vvol)

            volumes_i.append(calc[id_i].vol)
            volumes_p.append(calc[id_p].vol)
            #calculate average distanced between impurity and 6 most close atoms
            average_d = local_surrounding(calc[id_i].end.xcart[-1], calc[id_i].end, 6, control = 'av')
            daver_i.append(average_d)
            names.append(calc[id_i].name)

            # print average_d, 4./3 * math.pi * (average_d - 1.44)**3

            a.append(calc[id_i].hex_a)
            c.append(calc[id_i].hex_c)
            # print (a[i]**2 + (a[i]*2.88/2.93)**2)**0.5 / 2
            daver_p = ( (4./3 * a[i]**2 + c[i]**2 / 4)**0.5) / 2
            



            # cl_un = calc[(it_un,ise_un,v)] #unrelaxed; only for oxygen
            print calc[id_i].hex_a, calc[id_i].hex_c, calc[id_i].hex_c/calc[id_i].hex_a

            etot_i.append(calc[id_i].energy_sigma0)
            etot_p.append(calc[id_p].energy_sigma0)
            # etot_e.append(calc[id_e].energy_sigma0)
            # print 'I am here'
            # dE = etot_i[i]-etot_p[i] - e_at
            # dE = etot_i[i]#--751.78701 - e_at
            dE = calc[id_i].energy_sigma0 - calc[id_p].energy_sigma0 - e_at
            dElist.append(dE)

            id_e = (it_e,ise_e,v_e) # empty cell with frozen Ti
            try:            
                cl_e = calc[id_e]
                dEm = calc[id_e].energy_sigma0 - calc[id_p].energy_sigma0
                dEch = dE - dEm
                esol_m.append(dEm)
            except:
                print id_e, 'not found'
                dEm = 0
                dEch = 0



            # print 'I am here2'

            # dEch = calc[id_i].energy_sigma0 - calc[id_e].energy_sigma0 - e_at

            # dEr = cl_un.energy_sigma0 - cl_i.energy_sigma0 # almost the same as dEm!!!
            # print 'dEr', dEr, "dEm", dEm
            #Calculate charge density at impurity position for id_e
            density = 0
            if 0:
                calc[id_i].get_chg_file()
                xred1 = calc[id_i].end.xred[-1] #position of impurity
                chgfile = calc[id_e].get_chg_file()
                density = chg_at_point(chgfile, xred1, )
                print density

            #charge differences
            # print cal_chg_diff( calc[id_i].dir,  calc[id_i].get_chg_file(),  calc[id_e].get_chg_file() )

            """Bader"""
            # calc[id_i].bader_analysis()
            
            """DOS"""
            # print plot_dos(calc[id_i], calc[id_e], 'partial', path = 'dos/hs')
            # print "{:15s} & {:5.2f} & {:5.2f} & {:5.2f}".format( *plot_dos(calc[id_i], calc[id_e], 'partial')  )
            # plot_dos(calc[id_p], 'total')
            
            """Energy contributions"""
            if 0:
                cl_i.read_results(); cl_p.read_results();
                # print  cl_p.energy.ewald
                sum = 0
                sumdE = 0
                dEdict = {}
                for key in cl_i.energy.__dict__:
                    # print key, cl_i.energy.__dict__[key]
                    sum+=cl_i.energy.__dict__[key]
                    dEdict[key] = cl_i.energy.__dict__[key] - cl_p.energy.__dict__[key]
                    sumdE+=dEdict[key]
                dEdict['pawdc'] = dEdict['pawdc1']+dEdict['pawdc2']
                print 'ewald+alpha', 'hartree', 'xc', 'bands', 'pawdc', 'alpha'
                print "{:5.1f} {:5.1f} {:5.1f} {:5.1f} {:5.1f} {:5.1f}".format(dEdict['ewald']+dEdict['alpha'],\
                    dEdict['hartree'], dEdict['xc'], dEdict['bands'], dEdict['pawdc'],  dEdict['alpha'])#+dEdict['hartree']+dEdict['bands']+dEdict['xc']
                
                # A = dEdict['ewald']
                # B = dEdict['hartree']+dEdict['xc']+dEdict['bands']+dEdict['pawdc']+dEdict['alpha']+dEdict['atomic']
                # temp.append(A+B) 
                # print "{:5.1f} {:5.1f} {:5.1f} {:5.1f}".format(A, B, A+B, temp[i+1]-temp[i])
                
                print sumdE
                # print 'sum and sigma0', sum, cl_i.energy_sigma0

                # print cl_p.energy.ewald, cl_i.vol
                # print cl_p.energy.hartree
                


            table.append( "{:20s} & {:5.2f} & {:5.2f}& {:5.2f}      & {:5.2f} & {:5.2f} & {:5.2f} & {:8.0f} & {:8.0f} & {:5.4f} & {:5.4f}".format(
                names[i], dE, dEch, dEm, voro_i[i], daver_i[i],  daver_p, pres_i[i], pres_p[i],  a[i], c[i], density ) )
            i+=1

        #find minimum of mech energy:
        redu = 9
        power = 3
        an = a[:redu]
        coeffs1 = np.polyfit(an, esol_m[:redu], power)
        fit_func = np.poly1d(coeffs1)
        # fit_and_plot(a = (np.linspace(min(an), max(an)) , fit_func( np.linspace(min(an), max(an))  ), '-' ),
        #  b = (an, esol_m[:redu], 'o' ) )
        
        # print fit_func.deriv()
        # print fit_func.deriv().r
        minimum_m =  fit_func.deriv().r
        print minimum_m
        print fit_func(minimum_m)
        # os.exit()

        caption = """ Solution energy depending on pressure for carbon and oxygen. $\Delta E$ = E(Ti$_{96}$X) - E(Ti$_{96}$) - E(C$_g$) ({\AA}),
        where E(Ti$_{96}$) is the energy of cell with 
        Ti atoms wit the same lattice constants as in Ti$_{96}$X, E(C$_g$) is energy of one atom in graphite
        or one oxygen atom in O$_2$ molecule..
        $\Delta E_{ch}$ and $\Delta E_{m}$ are chemical and mechanical contributions of solution energy.
        V$_{Vor}$ is Voronoi volume {\AA$^3$},
        $d_{ai}$ is average distance between impurity and six neighbouring Ti atoms ({\AA}) and $d_{ab}$ is corresponding distance in bulk cell,  P1 is external pressure in cell with impurity
        P2 is pressure in cell without impurity (MPa), a and c are hexagonal lattice constants ({\AA}).  """
        #CHG is charge density in the center of octahedral site, from where the impurity was removed ( el \AA$^-3$).
        


        header = "Name & $\Delta E$ & $\Delta E_{ch}$ & $\Delta E_{m}$ & V$_{Vor}$ & $d_{ai}$ & $d_{ab}$ & P1 & P2 & a & c"  


        table.append('\\hline')


        #plot
        if 1:
            power = 7
            name = ''
            xlabel = 'Pressure, MPa'
            ylabel = 'Solution energy, eV'
            x1 = pres_i
            y1 = etot_i
            x2 = pres_p
            y2 = etot_p
            # print y1
        
            x_range1 = np.linspace(min(x1), max(x1))
            x_range2 = np.linspace(min(x2), max(x2))
            v_range = np.linspace(min(volumes_i), max(volumes_i))

            PdV = np.asarray([P*Vsol*1e6 / 1.6 / 1e-19 for P in x_range1])

            plt.ylabel(ylabel)
            plt.xlabel(xlabel)

            if 1:
                plt.ylabel(ylabel)
                plt.xlabel('hcp lattice constant $a$, $\AA$')
                # fit_func = np.poly1d( np.polyfit(volumes_i, dElist, 6)  )
                plt.plot(a, dElist, color, label = imp)
                # plt.plot(v_range,   fit_func(v_range), 'r-', label = 'fit')
                # tck = interpolate.splrep(volumes_i, dElist, s = 0)
                # fit_spline = interpolate.splev(v_range, tck, der=0)
                # plt.plot(v_range,   fit_spline, 'b-', label = 'fit_spline')
                # fit_func1 = interp1d(volumes_i, dElist, kind='zero')
                # plt.plot(v_range,   fit_func1(v_range), 'b-', label = 'fit_spline')
                plt.axvline(2.9371, color='k')

            elif 0:
                coeffs1 = np.polyfit(x1, y1, power)        
                coeffs2 = np.polyfit(x2, y2, power)
                fit_func1 = np.poly1d(coeffs1)
                fit_func2 = np.poly1d(coeffs2)
                f = fit_func1.deriv()
                print "The minimum energy  ", fit_func1(f.r[power-2]).real

            elif 0:
                x1.reverse()
                y1.reverse()
                x2.reverse()
                y2.reverse()
                # fit_func1 = interp1d(x1, y1, kind='cubic')
                # fit_func2 = interp1d(x2, y2, kind='cubic')
                # print x1
                # print y1
                tck1 = interpolate.splrep(x1, y1,)
                tck2 = interpolate.splrep(x2, y2,)
                fit_y1 = interpolate.splev(x_range1, tck1, der=0)
                fit_y2 = interpolate.splev(x_range1, tck2, der=0)

                # print fit_y1
                df = fit_y1 - fit_y2 - e_at

                print "The minimum energy  ", interpolate.sproot(tck1)
                plt.plot(x_range1 , df,             'g-', label = 'diff1_imp-pure') #take Pulay

            else:
                # eos1 = EquationOfState(volumes_i, y1, 'birchmurnaghan')
                eos1 = EquationOfState(np.asarray(x1)+25000, y1, 'murnaghan')
                v0, e0, B = eos1.fit()
                f, xr1, yr1 = eos1.plot('images/eos1.png', xr = x_range1+25000)
                # print yr1

                eos2 = EquationOfState(np.asarray(x2)+25000, y2, 'murnaghan')
                v0, e0, B = eos2.fit()
                f, xr2, yr2 = eos2.plot('images/eos2.png', xr = x_range1+25000)
                # print eos1.e
                # #print "a = ", alist[2]
                # print '''
                # v0 = {0} A^3
                # E0 = {1} eV
                # B  = {2} eV/A^3'''.format(v0, e0, B)
                vnew = 0
                # for v in xr1:
                #     print vnew 
                #     vnew = v

                dif = 8
                # df = yr1[dif:] - yr2[:-dif] - e_at
                df = yr1 - yr2 - e_at
                # xr1 = xr1[dif:]
                # print "Volume diff is ", xr1[dif] - xr1[0] 
                # xlabel = 'Volume, A^3'
                # print yr1[1:]
                # print yr2[:-1]
                # print df
                # print xr1
                # print xr2
                plt.figure()
                plt.title(name)
                plt.ylabel(ylabel)
                plt.xlabel(xlabel)

                plt.plot(xr1 - 25000, df,             'g-', label = 'diff1_imp-pure') #take Pulay

            print "Energy of hs443C.f", calc[('hs443C.f', '93', 1)].energy_sigma0


            # fit_y1 = fit_func1(x_range1[50:-50]); 
            # fit_y2 = fit_func2(x_range1[50:-50]); 

            # dff = fit_func1 - fit_func2# - e_at



            # plt.plot(x1, y1, 'ro', label = 'init')
            # plt.plot(x_range1, fit_func1(x_range1), 'r-', label = 'fit_imp')
            # plt.plot(x_range1, fit_func2(x_range1), 'b-', label = 'fit_pure')

            # plt.plot(x_range1, fit_y1 , 'k-', label = 'fit_imp')
            # plt.plot(x_range1, fit_y2 , 'b-', label = 'fit_pure')
            # plt.plot(x1, y1, 'ko', label = 'diff_imp-pure')
            # plt.plot(x2, y2, 'bo', label = 'diff_imp-pure')


            # plt.plot(x_range1 + 0, df+PdV,            'b-', label = 'PdV') #take Pulay
            # plt.plot(x_range1 + 0, dff(x_range1), 'r-', label = 'diff2_imp-pure') #take Pulay
    plt.legend(loc = 2)
    # plt.show()
    # plt.ion()
    plt.tight_layout()

    plt.savefig(path_to_images+"esol_a"+'.png', dpi = 300, format='png')



    latex_table(table, caption, label = 'tab:solepres', header = header, filename = 'tab_solepres.tex', fullpage = '*')


    return










def table_and_picture_segen(calc, conv = {}, varset = []):
    """
    Table with segregation energies, volumes;
    By default, voronoi volume, and chemical and mechanical decomposition were calculated for version 2. 
    You may need to run seg_energy_vs_voronoi() function first to fill self.vorovol list.

    """


    def cycle(name = 'empty'):
        """Read results and aggregate them
        Need seg_list, set, ver, calc,conv,varset and  base,  variables defined before running the function

        readfiles = False in res_loop allows significantly improve speed by skipping unneeded reading of OUTCAR and running of plot_conv.

        """
        for key in seg_list:#[0:1]:
            # print key
            # if 'Oi' in key: continue
            # print "id is", key, set, ver
            list = res_loop(key,set,ver,calc,conv,varset, 'e_seg', (base,) , voronoi = False , readfiles = 0)# matrix_id = , matrix_id_b =  )
            # list = ['','']
            #Make name
            if 'i' in list[0]:
                n1 = list[0].index('i');
                list[0]= name + '('+list[0][n1+1]+')';


            id2 = (key, set, 2)
            id_m2 = (key+'.m', '8'+set[1:], 2)




            #Add Voronoi volume to the list 
            ivor = -1           
            if 'Ci' in key:
                if 'C' not in list[ivor]: #additional check -  use last word segimp
                    raise RuntimeError 
                listkey = 'C'
                iimp_v = 0; iimp_x = -2
                print 'coordinates of impurities', calc[id2].end.xcart[-2] #print coordinates of impurities


            elif 'Oi' in key:
                if 'O' not in list[ivor]:
                    raise RuntimeError  
                listkey = 'O'
                iimp_v = 1; iimp_x = -1

            else:
                #For volume cases; Voronoi volume and distances are shown for carbon
                listkey = 'C' #adding to C list
                iimp_v = 0; iimp_x = -2 #information for carbon

            if hasattr(calc[id2].init, 'vorovol'):
                list[ivor] = calc[id2].init.vorovol[iimp_v] # replace it with init vorovol
                list.append(calc[id2].vorovol[iimp_v])#end vorovol
            else:
                print 'Warning!!! Calculation', id2, 'does not have .init.vorovol'
                list[ivor] = 0
                list.append(0)

            #calculate distortions
            av_dev = local_surrounding(calc[id2].init.xcart[iimp_x], calc[id2].init, 6, control = 'av_dev', periodic = True) #The init deviation of segregation interstitial geometry from ideal octahedron
            list+=[av_dev[0],]
            av_dev = local_surrounding(calc[id2].end.xcart[iimp_x], calc[id2].end, 6, control = 'av_dev', periodic = True) #The end deviation of segregation interstitial geometry from ideal octahedron
            list+=[av_dev[0],]


            density = 0; den_diff = 0
            if 0:
                # Calculate charge density at impurity position for id_e
                chgfile = calc[id_m2].get_chg_file("CHG")

                if 'vms' in key:
                    xredC = calc[id2].end.xred[-2] #position of impurity
                    xredO = calc[id2].end.xred[-1] #position of impurity
                    vmsdensityC = chg_at_point(chgfile, xredC, )
                    vmsdensityO = chg_at_point(chgfile, xredO, )
                    print 'vms density C O', vmsdensityC, vmsdensityO



                xred1 = calc[id2].end.xred[iimp_x] #position of impurity
                density = chg_at_point(chgfile, xred1, )
                if os.path.exists(chgfile): 
                    print 'Density is ',  density
                else: 
                    print "!No", chgfile
                

                if iimp_x is -2:
                    den_diff = density - vmsdensityC # minus the density in vms
                elif iimp_x is -1:
                    den_diff = density - vmsdensityO # minus the density in vms

            list.append(density)

            list.append(den_diff)







            """Energy contributions"""
            if 0:
                cl_i =  calc[(key, set, 2)]; #impurity at interface
                cl_p = calc[(base, set, 2)]; #impurity in the bulk
                cl_i.read_results(); cl_p.read_results();
                # print  cl_p.energy.ewald
                sum = 0
                sumdE = 0
                dEdict = {}
                for dictkey in cl_i.energy.__dict__:
                    print dictkey, cl_i.energy.__dict__[dictkey]
                    sum+=cl_i.energy.__dict__[dictkey]
                    dEdict[dictkey] = cl_i.energy.__dict__[dictkey] - cl_p.energy.__dict__[dictkey]
                    sumdE+=dEdict[dictkey]
                dEdict['pawdc'] = dEdict['pawdc1']+dEdict['pawdc2']
                dEdict['sumdE'] = sumdE 
                dEdict['id_i']  = (key, set, 2)
                print sum
                print 'ewald+alpha', 'hartree', 'xc', 'bands', 'pawdc', 'alpha'
                print "{:5.1f} {:5.1f} {:5.1f} {:5.1f} {:5.1f} {:5.1f}".format(dEdict['ewald']+dEdict['alpha'],\
                    dEdict['hartree'], dEdict['xc'], dEdict['bands'], dEdict['pawdc'],  dEdict['alpha'])#+dEdict['hartree']+dEdict['bands']+dEdict['xc']
                
                print sumdE

                list.append(dEdict)








            # print list
            aggregate_list[listkey].append(list)


        return



    """2. Read results"""
    if 1:
        aggregate_list = {'C':[],'O':[]}; 

        # name = 'T1m'
        # base = 't111gCvOvms'
        # seg_list  = ['t111gCvOvms','t111gCi1Ov', 't111gCi2Ov','t111gCi3Ov','t111gCi4Ov','t111gOi1Cv','t111gOi2Cv','t111gOi3Cv','t111gOi4Cv',] 
        # set = '93kp9';  ver = range(1,5)
        # cycle(name)
        # # os._exit(1)

        # name = 'T1g'
        # base = 't111sgCvOvms'
        # seg_list = ['t111sgCvOvms','t111sgCi1Ov', 't111sgCi2Ov', 't111sgCi3Ov', 't111sgCi4Ov', 't111sgCi5Ov', 't111sgCi6Ov', 't111sgCi7Ov', 't111sgOi1Cv', 't111sgOi2Cv', 't111sgOi3Cv', 't111sgOi4Cv', 't111sgOi5Cv', 't111sgOi6Cv', 't111sgOi7Cv', ]
        # set = '93kp9';  ver = range(1,6)
        # cycle(name)

        # name = 'C1'
        # base = 'c1gCvOvms'
        # seg_list = ['c1gCvOvms','c1gCi1Ov', 'c1gCi2Ov', 'c1gOi1Cv', 'c1gOi2Cv', ] 
        # set ='93kp7';  ver = range(1,6)
        # cycle(name)

        # name = 'T2'
        # base = 't21gCvOvms'
        # seg_list = ['t21gCvOvms','t21gCi1Ov', 't21gCi2Ov', 't21gCi3Ov', 't21gCi4Ov', 't21gOi1Cv', 't21gOi2Cv', 't21gOi3Cv', 't21gOi4Cv', ]
        # set = '93';  ver = range(1,5)
        # cycle(name)

        name = 'S7'
        base = 'csl71sgCvOvms'
        seg_list = ['csl71sgCvOvms','csl71sgCi1Ov', 'csl71sgCi2Ov', 'csl71sgCi3Ov', 'csl71sgCi4Ov', 'csl71sgCi5Ov', 'csl71sgCi6Ov', 'csl71sgOi1Cv', 'csl71sgOi2Cv', 'csl71sgOi3Cv', 'csl71sgOi4Cv', 'csl71sgOi5Cv', 'csl71sgOi6Cv', ]
        set = '93';  ver = range(1,6)
        cycle(name)

    
        calc['segreg_table'] = aggregate_list
    else:
        aggregate_list = copy.deepcopy(calc['segreg_table'])


    """2. Produce table"""
    if 0:
        s1 = 8
        s2 = 22
        header0 =  'Void'.ljust(s1)+'&'+'\multicolumn{2}{c}{$E_{seg}$}'.center(s1)+'&'+'\multicolumn{2}{c}{$V_{seg}$}'.center(s1)+'&'+\
        '\multicolumn{2}{c}{$V_{vor}$}'.center(s1)+'&'+'Distance'.ljust(s1)+'&'+'NN'.ljust(s1)+'&'+\
        '\multicolumn{2}{c}{$E_{seg}(2)$}'.ljust(s1)+'&'+'\multicolumn{2}{c}{$E_{ch}(2)$}'.ljust(s1)+'&'+'\multicolumn{2}{c}{$E_{m}(2)$}'.ljust(s1)

        header = ' '.ljust(s1)+'&'+'C'.center(s1)+'&'+'O'.center(s1)+'&'+'C'.center(s1)+'&'+'O'.center(s1)+'&'+'C'.center(s1)+'&'+'O'.center(s1)+'&'+\
        'C,O '.center(s1)+'&'+'C,O '.center(s1)+'&'+\
        'C'.center(s1)+'&'+'O'.center(s1)+'&'+'C'.center(s1)+'&'+'O'.center(s1)+'&'+'C'.center(s1)+'&'+'O'.center(s1)
        

        caption = """
        $dev_{av}$  is average deviation from average value of six Ti-X distances ({m\AA})
        $dev_{max}$  is maximum deviation from average value of six Ti-X distances ({m\AA})
        """
        header_new = \
        "Void & $E_{seg}$ & $V_{seg}$  & D & NN & $E_{seg}(2)$  & $E_{ch}(2)$  & $E_{m}(2)$ & $V_{vor, init}$ & $V_{vor, end}$ & $dev_{av, init}$ & $dev_{av, end}$" #+" & chg & chgdiff "

        table = []


        if 0: #sort by carbon only
            aggregate_listCO = [C[:-1]+O[:-1] for C, O in zip(aggregate_list['C'], aggregate_list['O'])]


            for el in aggregate_listCO:
                print el

            print 'After sort:'

            aggregate_listCO.sort(key = itemgetter(1)) #sort carbon and oxygen by carbon energies
            
            for el in aggregate_listCO:
                print el


            aggregate_list['C'] = []
            aggregate_list['O'] = []
            for CO in aggregate_listCO:
                # print CO[0], CO[13:15]
                aggregate_list['C'].append(CO[0:len(CO)/2])
                aggregate_list['O'].append(CO[len(CO)/2:])
            # print len(aggregate_list['C']), len(aggregate_list['O'])
        else:
            #sort separatly
            aggregate_list['C'].sort(key = itemgetter(1))
            aggregate_list['O'].sort(key = itemgetter(1))




        caption ="""Optimized segregation energies $E_{seg}$ (meV) and excess volumes $V_{seg}$ ({\AA}) 
                    obtained by fitting of segregation energies of five different expansions of boundary.
                    Other characteristics are provided for one
                    specific cell, the expansion of which is almost concide with optimized one. 
                    D is the distance {\AA} from the impurity to the grain boundary plane.
                    NN is the number of neighbouring titanium atoms at distances smaller than 3 {\AA}.
                    The $E_{seg}(2)$ is segregation enrergy and its decomposition
                    into chemical $E_{ch}(2)$ and mechanical $E_{m}(2)$  contributions. 
                    The $V_{vor, init}$ and  $V_{vor, end}$ ({\AA$^3$}) are Voronoi volumes of the pore before
                    relaxation and after relaxation. The $dev_{av, init}$ and $dev_{av, end}$ are deviations
                    of pores from the ideal octahedron before and after relaxation."""
        
        """Заметки: S7(4) как буд-то бы выбивается из тренда. Отклонения от идеального октаэдра малы в начале и в конце, а энергии сегрегации
        все равно положительные. На самом деле это положение находится несколько в стороне от границы и несмотря на близкое к идеальному 
        междоузлие, все равно оказывает возмущения на границу зерна, что делает сегрегацию не выгодной - однако это не подкрепляется изменениями
        электронной плотности. В реальности трудно понять не выгодность этого положения. Возможно здесь (накладываются) играют роль периодические граничные условия"""



        for C,O in zip(aggregate_list['C'], aggregate_list['O']):
            # if float(C[7]) < 46: continue



            # table.append( 
            # C[0].ljust(s1)+"&"+C[1].center(s1)+"&"+O[1].center(s1)+"&"+C[2].center(s1)+"&"+\
            # O[2].center(s1)+"&"+'{0:.1f}'.format(C[-1]).center(s1)+"&"+'{0:.1f}'.format(O[-1]).center(s1)+"&"+\
            # (C[3]+','+O[3]).center(s1)+"&"+(C[4]+','+O[4]).center(s1)+"&"+\
            # (C[5]+'&'+O[5]).center(s1)+"&"+(C[6]+'&'+O[6]).center(s1)+"&"+(C[7]+'&'+O[7]).center(s1)  
            # )
            # print C

            #skip vms
            if 'vms' in C[0]: continue
            if 'vms' in O[0]: continue



            table.append( "{:15s}C & {:5.0f} & {:5.1f} & {:5.2f} & {:5d} & {:5.0f} & {:5.0f} & {:5.0f} & {:5.1f} & {:5.1f} & {:5.0f} & {:5.0f} ".format(*C)  ) # without chg and chgdiff
            table.append( "{:15s}O & {:5.0f} & {:5.1f} & {:5.2f} & {:5d} & {:5.0f} & {:5.0f} & {:5.0f} & {:5.1f} & {:5.1f} & {:5.0f} & {:5.0f} ".format(*O)  )

     


        # latex_table(table, caption = '', label = 'tab:e_seg', header = header, header0 = header0  )
        latex_table(table, caption = caption,  label = 'tab:e_seg', header = header_new, fullpage = '*', filename = 'tab_e_seg.tex', writetype = 'w' )






        if 0: #horisontal table - should be repaired, because aggregate_listC was changed
            namel = ['Void              ']
            esegC = ['$\Delta E_{seg}$ C']
            esegO = ['$\Delta E_{seg}$ O']
            evolC = ['$\Delta V_{seg}$ C']
            evolO = ['$\Delta V_{seg}$ O']
            evvC  = ['$V_{vor}$ C       ']
            evvO  = ['$V_{vor}$ O       ']

            print '\\begin{table}'
            print '\\label{tab:seg_twins}'
            print '\\caption{}'
            print '\\begin{tabular}{l'+ 15*'c'+'r}'
            print '\\hline'
            # print header
            # print '\\hline'
            # print header2
            # print '\\hline'
            for C,O in zip(aggregate_listC, aggregate_listO):
                # print len(l)
                if 't111g' in C[0]: 
                    name = 'T1m('
                elif 't111sg' in C[0]:
                    name = 'T1g('
                elif 'c1g' in C[0]:
                    name = 'C1('
                elif 't21g' in C[0]:
                    name = 'T2('
                elif 'csl71sg' in C[0]:
                    name = 'S7('

                n1 = C[0].index('i');
                n2 = O[0].index('i'); 
                assert int(n1), int(n2)
                name += C[0][n1+1]+')'; 
                i =  int(name.split('(')[1].split(')')[0])
                if i > 4: continue
                if 'S' in name: continue
                # print name
                namel.append(name)
                esegC.append(C[1])
                esegO.append(O[1])
                evolC.append(C[2])
                evolO.append(O[2])
                evvC.append('{0:.1f}'.format(C[-1]))
                evvO.append('{0:.1f}'.format(O[-1]))

            for e in namel:
                print e.ljust(s1)+"&", 
            print '\\\\'
            print '\\hline'

            for e in esegC:
                print e.center(s1)+"&", 
            print '\\\\'

            for e in esegO:
                print e.center(s1)+"&", 
            print '\\\\'

            for e in evolC:
                print e.center(s1)+"&", 
            print '\\\\'

            for e in evolO:
                print e.center(s1)+"&", 
            print '\\\\'

            for e in evvC:
                print e.center(s1)+"&", 
            print '\\\\'

            for e in evvO:
                print e.center(s1)+"&", 
            print '\\\\'
                # print name.ljust(s1)+"&"+C[1].center(s1)+"&"+O[1].center(s1)+"&"+C[2].center(s1)+"&"+O[2].center(s1)+"&"+'{0:.1f}'.format(C[4]).center(s1)+"&"+'{0:.1f}'.format(O[4]).center(s1)+   '\\\\'

            print '\\hline'
            print '\\end{tabular}'
            print '\\end{table}'
        # print aggregate_listC
        # print aggregate_listO



    """2. Make some transformations for convenience"""
    if 0:
        rev_list = {}
        # for el in 'CO':
        #     aggregate_list[el][:] = [x+[i] for i, x in enumerate(aggregate_list[el]) if not 'vms' in x[0]  ] #remove unneeded entries
        #     # aggregate_list[el][:] = [x for x in aggregate_list[el] if not 'C1(1)' in x[0]]
        #     aggregate_list[el].sort(key = itemgetter(5))

        #     i = [i for i, x in enumerate(aggregate_list[el]) if 'C1(1)' in x[0]][0]
        #     # print aggregate_list[el]
        #     print i
        #     aggregate_list[el].insert(-3, aggregate_list[el].pop(i))



        aggregate_listCO = [C+O for C, O in zip(aggregate_list['C'], aggregate_list['O'])]
        aggregate_listCO.sort(key = itemgetter(5)) #sort carbon and oxygen by carbon energies

        aggregate_list['C'] = []
        aggregate_list['O'] = []
        for CO in aggregate_listCO:
            # print CO[0], CO[13:15]
            aggregate_list['C'].append(CO[0:len(CO)/2])
            aggregate_list['O'].append(CO[len(CO)/2:])
        # print len(aggregate_list['C']), len(aggregate_list['O'])


        for el in 'CO': 
            i = [i for i, x in enumerate(aggregate_list[el]) if 'C1(1)' in x[0]][0] #move C1(1)
            aggregate_list[el].insert(-3, aggregate_list[el].pop(i)) #move C1(1)
            
            rev_list[el] =  zip(*aggregate_list[el])


    """2. Produce plot: segregation energies - chemical and mechanical parts"""
    el1= 'C'; el2 = 'O' 
    if 0:
        # name =  zip(*aggregate_list['C'])[0]
        # print aggregate_list['C']
        # print rev_list['C'][5]
        # print rev_list['C'][6]
        # print rev_list['C'][7]

        # plot_bar(xlabel = 'Segregation site', ylabel = 'Segregation energy (meV)',
        #     C = (rev_list['C'][0], rev_list['C'][1], "#4078D3" ),
        #     O = (rev_list['O'][0], rev_list['O'][1], "#F46529")
        #     )
        to = -4 #divide into two blocks
        # for el in 'C':
        #     plot_bar(xlabel = 'Segregation site', ylabel = 'Segregation energy (meV)',
        #         # C_seg_eq     = (rev_list[el][0], rev_list[el][1], 'k' ),
        #         b_seg_2_ch   = (rev_list[el][0][:to], rev_list[el][6][:to], "#F46529", el+' chemical'    ),
        #         c_seg_2_mech = (rev_list[el][0][:to], rev_list[el][7][:to], "#4078D3",   el+' elastic'  ),
        #         a_seg_2      = (rev_list[el][0][:to], rev_list[el][5][:to],   "#7C144D",    el+' total'  ),

        #         )

        if 1:
            xlabel = 'Segregation site'; ylabel = 'Segregation energy (meV)'
            chemical = ' Chemical'; elastic = ' Elastic'
        else:
            xlabel = u""; ylabel = u"Энергия сегрегации (мэВ)"
            chemical = u' Хим. вклад'; elastic = u' Мех. вклад'


        plot_bar (image_name = 'seg_energy_full',
            xlabel = xlabel, ylabel = ylabel,
            data1 = [
                (rev_list[el1][0][:to], rev_list[el1][5][:to], "#F46529",   el1    ),
                (rev_list[el2][0][:to], rev_list[el2][5][:to], "#7FB005",   el2     ),
                ],
            data2 = [
                (rev_list[el1][0][to:], rev_list[el1][5][to:], "#F46529",   el1    ),
                (rev_list[el2][0][to:], rev_list[el2][5][to:], "#7FB005",   el2    ),
                ],
            )

           
        plot_bar (image_name = 'seg_energy_elastic_chemical',
            xlabel = xlabel, ylabel = ylabel,
            data1 = [
                (rev_list[el1][0][:to], rev_list[el1][6][:to], "#F46529",   el1+ chemical    ),
                (rev_list[el1][0][:to], rev_list[el1][7][:to], "#4078D3",   el1+ elastic     ),
                ],
            data2 = [
                (rev_list[el1][0][to:], rev_list[el1][6][to:], "#F46529",   el1+ chemical    ),
                (rev_list[el1][0][to:], rev_list[el1][7][to:], "#4078D3",   el1+ elastic     ),
                ],
            data3 = [
                (rev_list[el2][0][:to], rev_list[el2][6][:to], "#7FB005",   el2+ chemical    ),
                (rev_list[el2][0][:to], rev_list[el2][7][:to], "#7C144D",   el2+ elastic     ),
                ],
            data4 = [
                (rev_list[el2][0][to:], rev_list[el2][6][to:], "#7FB005",   el2+ chemical    ),
                (rev_list[el2][0][to:], rev_list[el2][7][to:], "#7C144D",   el2+ elastic     ),
                ],
            )






    """3. Produce plot: segregation energies - distortion"""
    if 0:

        if 1:
            decontib = {};
            for el in el1,el2:
                decontib[el] = {}
                for key in 'ewald', 'alpha', 'hartree', 'xc', 'bands', 'pawdc', 'alpha', 'sumdE', 'id_i':
                    decontib[el][key] = []
                    for dic in rev_list[el][14]:
                        decontib[el][key].append( dic[key] )

            for eb, eseg, id_i  in zip(decontib['C']['bands'], decontib['C']['sumdE'], decontib['C']['id_i']) :
                if eb > 2 and eseg < -0.1:
                    print eb, eseg, id_i

            fit_and_plot(image_name = 'eseg_deband',
            xlabel = "Contribution of band energy in segregation energy (eV)", ylabel = "Segregation energy (eV)",
            xlog = 0, ylog =0, legend = 2,
            power = 1,
            # title = "Segregation energy vs change",
            # xlim = (-5,1.5), ylim = (None, 2),#0.4
            # xlim = (4,250), 
            # ylim = (None, 12),#0.4
            carbon   = (decontib[el1]['bands'], decontib[el1]["sumdE"], 'ko'),
            oxygen   = (decontib[el2]['bands'], decontib[el2]["sumdE"], 'ro'),
            )




        if 0:
            fit_and_plot(image_name = 'eseg_chgdiff',
            xlabel = "Difference of charge densities in interstital sites", ylabel = "Chemical segregation energy (meV)",
            xlog = 0, ylog =0, legend = 1,
            # title = "Segregation energy vs change",
            # xlim = (-5,1.5), ylim = (None, 2),#0.4
            # xlim = (4,250), 
            # ylim = (None, 12),#0.4
            carbon   = (rev_list[el1][13], rev_list[el1][6], 'ko'),
            oxygen   = (rev_list[el2][13], rev_list[el2][6], 'ro'),
            )



        if 0:
            vor_vol = {}
            for el in 'CO':
                vor_vol[el] = [(v - 8)**3 * 50 for v in rev_list[el][8]] #operation to make differences between volumes on picture more pronounced

            fit_and_plot(image_name = 'eseg_dist',
            xlabel = "Average deviation from ideal octahedron ($m\AA$)", ylabel = "Segregation energy (meV)",
            xlog = 1, ylog =1, scatter = 1, legend = 7,
            # title = "Segregation energy vs change",
            # xlim = (-5,1.5), ylim = (None, 2),#0.4
            xlim = (4,250), 
            # ylim = (None, 12),#0.4
            carbon   = (rev_list[el1][10], rev_list[el1][5], vor_vol[el1], 'k'),
            oxygen   = (rev_list[el2][10], rev_list[el2][5], vor_vol[el2], 'r'),
            )

        if 0:
            fit_and_plot(image_name = 'eseg_vorvol',
            xlabel = "Voronoi volume ($\AA^3$)", ylabel = "Segregation energy (meV)",
            xlog = 0, ylog =1, legend = 7,
            # title = "Segregation energy vs change",
            # xlim = (-5,1.5), ylim = (None, 2),#0.4
            xlim = (8.8, 10.6), 
            # ylim = (None, 12),#0.4
            carbon   = (rev_list[el1][8], rev_list[el1][5], 'ok'),
            oxygen   = (rev_list[el2][8], rev_list[el2][5], 'or'),
            )

            fit_and_plot(image_name = 'vorvol_dist',
            xlabel = "Average deviation from ideal octahedron ($m\AA$)", ylabel = "Voronoi volume ($\AA^3$)",
            xlog = 0,ylog =0, legend = 7,
            # title = "Segregation energy vs change",
            # xlim = (-5,1.5), ylim = (None, 2),#0.4
            # xlim = (8.8, 10.5), 
            # ylim = (None, 12),#0.4
            carbon   = (rev_list[el1][10], rev_list[el1][8], 'ok'),
            oxygen   = (rev_list[el2][10], rev_list[el2][8], 'or'),
            )

    return















def table_pairs(calc, conv, varset):
    """
    Interaction energies of atomic pairs in bulk Ti.
    """


    eCC, eOO, eCO = [], [], []
    names  = [], []
    distsCC, distsOO, distsCO = [], [] ,[]
    types = ['$0.5c$', '$a$', '$a+0.5c$', '$c$', '$\sqrt{3}a$', '$a+c$','$\sqrt{3}a +0.5c$', '$2a$','$2a+0.5c$','$\sqrt{3}a+c$','$1.5c$','$2a+c$','$a+1.5c$','$\sqrt{3}a+1.5c$','$2a+1.5c$']
    for v in range(1,16):
        CC = calc[('hs443CC', '93ns', v)]
        OO = calc[('hs443OO', '93ns', v)]
        CO = calc[('hs443CO', '93ns', v)]
        eCC.append(CC.energy_sigma0 )
        eOO.append(OO.energy_sigma0 )
        eCO.append(CO.energy_sigma0 )

        distsCC.append( image_distance(CC.end.xcart[-2], CC.end.xcart[-1], CC.end.rprimd, 2)  )
        distsOO.append( image_distance(OO.end.xcart[-2], OO.end.xcart[-1], OO.end.rprimd, 2)  )
        distsCO.append( image_distance(CO.end.xcart[-2], CO.end.xcart[-1], CO.end.rprimd, 2)  )

        # names.
    # print len(eCC)
    table = []
    # table.append( "N & type  &  d1  &  d2  &  $\Delta E(C-C)$  & $\Delta E(O-O)$  & $\Delta E(C-O)$  &  $\Delta_{av}$ & $E_{for}$"   )
    header = "\\textnumero & Type  &  $d_1$(C-C) &  $d_1$(O-O) &  $d_1$(C-O) &  $d_2$(C-O)  &  $\Delta E$(C-C)  & $\Delta E$(O-O)  & $\Delta E$(C-O)  &  $\Delta E_{av}$ & $E_{for}$(C-O)"
    caption = """
    Interaction energies for C-C, O-O and C-O pairs relative to configuration 15 (meV). 
    The first and second minimum distances $d_1$ and $d_2$ between pair atoms are provided in {\AA}
    for C-C. $\Delta E_{av}$ = ($\Delta E$(C-C) + $\Delta E$(O-O))/2. 
    $E_{for}$ is formation energy of C-O pair from C and O in Ti$_{96}$ cells (meV/atom).
    The type denotes the direction in hcp lattice along which the pair is located.
    """
    eC = -10.826 #h443C
    eO = -10.520 #h443O 
    eTi  = -7.831 #h443
    dedCO = []
    defCC_OO_CO = []
    for v in range(15):
        deCC = (eCC[v]-eCC[14])*1000
        deOO = (eOO[v]-eOO[14])*1000
        deCO = (eCO[v]-eCO[14])*1000
        efCC = (eCC[v]- 2*eC - eTi*96)*1000/2
        efOO = (eOO[v]- 2*eO - eTi*96)*1000/2
        efCO = (eCO[v]- eC - eO - eTi*96)*1000/2
        # table.append("{0:2d} & {1:20s} & {2:4.2f} & {3:5.2f} & {4:3.0f} & {5:3.0f} & {6:3.0f} & {7:3.0f} & {8:3.0f} ".format(
        #                v+1, types[v], distsCC[v][0], distsCC[v][1], deCC ,deOO ,deCO, 0.5*(deCC+deOO), efCO ) )  #only one distance 

        table.append("{:2d} & {:20s} & {:4.2f} & {:4.2f}& {:4.2f} & {:5.2f} & {:3.0f} & {:3.0f} & {:3.0f} & {:3.0f} & {:3.0f} ".format(
                       v+1, types[v], distsCC[v][0],distsOO[v][0],distsCO[v][0], distsCC[v][1], deCC ,deOO ,deCO, 0.5*(deCC+deOO), efCO ) )

        # dedCO.append( (round(distsCO[v][0],2), round(deCO, 0) ) )
        defCC_OO_CO.append( (round(distsCO[v][0],2), round(efCC, 0), round(efOO, 0), round(efCO, 0) ) )
    latex_table(table, caption, header = header, label = 'tab:epairs', fullpage = '*', filename = 'tab_epairs.tex')
    print defCC_OO_CO 
    return


def table_coseg(calc, conv, varset):
    """
    The main table with co-segregation information

    """


    def plot_coseg(data, image_name, split = -4):
        data.sort(key = itemgetter(1))

        rev_data =  zip(*data)



        if 1:
            to = split #divide into two blocks

            if 1:
                xlabel = 'Segregation site'; ylabel = 'Energy (meV)'
            else:
                xlabel = u""; ylabel = u"Энергия сегрегации (мэВ)"


            plot_bar (image_name = image_name,
                xlabel = xlabel, ylabel = ylabel, bottom = 0.3, barwidth = 0.2,

                data1 = [
                    (rev_data[0][:to], rev_data[1][:to], "#F46529",   '$\Delta E_2  $'           ),
                    (rev_data[0][:to], rev_data[2][:to], "#7FB005",   '$E_{pair, GB}$'     ),
                    (rev_data[0][:to], rev_data[3][:to], "#7C144D",   '$E_{coseg}   $'        ),
                    ],
                data2 = [
                    (rev_data[0][to:], rev_data[1][to:], "#F46529",   '$\Delta E_2  $'          ),
                    (rev_data[0][to:], rev_data[2][to:], "#7FB005",   '$E_{pair, GB}$'    ),
                    (rev_data[0][to:], rev_data[3][to:], "#7C144D",   '$E_{coseg}   $'       ),
                    ],
                )
        return




    coseg_list_T1m = [      
    't111gCOi1.4-3'   ,
    't111gCOi2.2-1'   ,
    't111gCOi3.4-1'   ,
    't111gCOi4.4-4is' ,
    't111gCOi5.2-2is' ,
    't111gCOi6.1-1is' ,
    't111gCOi7.3-3is' ,
    't111gCOi8.2-3'   ,
    't111gCOi9.4-4ms' ,
    't111gCOi10.2-2ms',
    't111gCOi11.4-3'  ,
    't111gCOi12.1-3'  ,
    't111gCOi13.4-2'  ,
    't111gCOi14.2-1'  ,
    't111gCOi15.4-4is',
    't111gCOi16.2-2is',
    't111gOCi1.4-3'   ,
    't111gOCi2.2-1'   ,
    't111gOCi3.4-1'   ,
    't111gOCi8.2-3'   ,
    't111gOCi11.4-3'  ,
    't111gOCi12.1-3'  ,
    't111gOCi13.4-2'  ,
    't111gOCi14.2-1'  ,
    ]
    coseg_list_T1g = [     
    't111sgCOi17.7-4',
    't111sgCOi28.4-5', 
    't111sgCOi29.4-6',
    't111sgCOi30.7-5',
    't111sgOCi8.4-7' ,
    't111sgOCi9.6-2' ,
    't111sgOCi15.6-6',
    't111sgOCi22.6-3',
    't111sgOCi25.6-7',
    't111sgOCi27.6-2',
     ]


    coseg_list_C1 = [
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
    'c1gCOi10.1' ,
    ]

    coseg_list_T2 = [ 
    't21gCOi1.4-3' ,
    't21gCOi2.3-2' ,
    't21gCOi3.2-2' ,
    't21gCOi4.2-1' ,
    't21gCOi5.3-1' ,
    't21gCOi6.4-2' ,
    't21gCOi7.2-1' ,
    't21gCOi8.4-1' ,
    't21gCOi9.4-2' ,
    't21gCOi10.3-3',
    't21gCOi11.4-3',
    't21gCOi12.3-2',
    't21gCOi13.4-1',
    't21gCOi14.1-4',
    't21gCOi15.3-2',
    't21gCOi16.4-2',
    't21gOCi1.4-3' ,
    't21gOCi2.3-2' ,
    't21gOCi3.2-2' ,
    't21gOCi4.2-1' ,
    't21gOCi5.3-1' ,
    't21gOCi6.4-2' ,
    't21gOCi7.2-1' ,
    't21gOCi8.4-1' ,
    't21gOCi9.4-2' ,
    't21gOCi10.3-3',
    't21gOCi11.4-3',
    't21gOCi12.3-2',
    't21gOCi13.4-1',
    't21gOCi14.1-4',
    't21gOCi15.3-2',
    't21gOCi16.4-2',
    ]

    coseg_list_S7 = [ 
    'csl71sgCOi1.6-4' ,
    'csl71sgCOi11.4-4',
    'csl71sgCOi22.1-2',
    'csl71sgOCi6.5-4' ,
    'csl71sgOCi7.1-6' ,
    ]




    def process_list(it_coseglist, it_base, ise, verlist, calc, conv, varset, run_res = False, bname = None, eCObulk = None):
        #bname - boundary name
        # eCObulk - distances and interaction energies for CO in corresponding cell; !mix from cells with boundaries and hs443CO
        table, e_seg = [], []
        tablelist = []
        for it in it_coseglist:
            if run_res:
                res_loop(it,ise,verlist, calc, conv, varset,  'coseg', (it_base,ise,1))

            id =(it,ise, 2)

            st = calc[id].end
            gbpos2 = calc[id].gbpos 
            gbpos1 = gbpos2 - st.rprimd[0][0]/2.
            d1 = st.xcart[-2][0] - gbpos2
            d2 = st.xcart[-1][0] - gbpos2

            dd1, dd2 = image_distance(st.xcart[-2], st.xcart[-1], st.rprimd, order = 3)
            # print st.xcart[-2][0], st.xcart[-1][0]
            e = round(calc[id].e_seg)
            if (e in e_seg) or (e+1 in e_seg) or (e-1 in e_seg) : 
                # print id, "was not included"
                continue 

            e_seg.append(e)
            # print e_seg

            #Name parser
            naml = it.split('.')
            
            num = naml[0].split('i')[1] # number of configuration
            
            voids = naml[1]
            if 'is' in voids: voids = voids.split('i')[0]
            if 'ms' in voids: voids = voids.split('m')[0]

            imp = naml[0].split('i')[0].split('g')[1]
            twin = naml[0].split('g')[0] + 'g'
            # print it
            # print bname+' '+num+' '+voids+' '+imp
            
            if len(voids) == 3:
                name = bname+' ('+imp[0]+voids[0:2]+imp[1]+voids[2]+')'

                seg1 = twin+imp[0]+'i'+voids[0] + imp[1]+'v'
                seg2 = twin+imp[1]+'i'+voids[2] + imp[0]+'v'
                e_seg_sum = calc[(seg1, ise, 1)].e_seg +  calc[(seg2, ise, 1)].e_seg #sum of independent segregation energies
                # de3 = '{:5.0f}'.format(calc[id].e_seg - e_seg_sum)
                e_pair_gb = calc[id].e_seg - e_seg_sum
                #Find $\Delta E$ (C-O) for bulk hs443CO cell  with approximately the same dd1 distance between impurities
                mindr = 100
                for t in eCObulk:
                    dr = abs(t[0]-dd1)
                    if dr < mindr: 
                        mindr = dr
                        e_pair_bulk = t[1]
                # print dd1, de
                # ecoseg = '{:5.0f}'.format(de3 - de)
                ecoseg = e_pair_gb - e_pair_bulk

            else: #For C1 (CO1) configuration; cosegregation energies was not calculated
                name = bname+' ('+imp+voids+')'
                # print naml
                e_pair_gb = 0.001
                ecoseg = 0.001

            tablelist.append(  [name , calc[id].e_seg, e_pair_gb, ecoseg,  calc[id].v_seg, d1, d2, dd1, dd2]  )

        tablelist.sort(key = itemgetter(1))
        for s in tablelist:
            table.append( 
            "{:25s} & {:5.0f} & {:5.0f} & {:5.0f} & {:5.1f} & {:5.2f} & {:5.2f} & {:5.2f} & {:5.2f}".format( *s  )   )

        return table, tablelist






    table1, table2 = [], []
    header = "Name & $\Delta E_2$ & $E_{pair, GB}$ & $E_{coseg}$ & $\Delta V_2$ & $d_{gb}$(C) & $d_{gb}(O)$ & $d_1$(C-O) & $d_2$(C-O)"
    
    caption = ['', '']
    boundaries = ["T1m, T1g  and S7", "C1 and T2"]
    for i in 0, 1:
        caption[i] = """ Co-segregation of C and O atoms at """ + boundaries[i]  +""" boundaries.  
        The types of segregation sites and impurities are shown in parenthesis.
        The $\Delta E_2 = E({\\rm CO}, gb) - E( {\\rm C}, bulk; {\\rm O}, bulk)$.
        $E_{pair, GB}$ = $\Delta E_2$ -- ($E_{seg}$(C) + $E_{seg}$(O)).
        $E_{coseg}$ = $E_{pair, GB}$ -- $E_{pair, bulk}$, 
        where $E_{pair, bulk} = \Delta E = E( {\\rm CO}, bulk - E( {\\rm C}, bulk; {\\rm O}, bulk)$
        is interaction energy between C and O atoms in bulk for corresponding inter-atomic distances 
        (obtained from cells with twin boundaries and Ti$_{96}$CO cell). All energies are in meV.
        $\Delta V_2$ is corresponding change of volume ({\AA$^3$}). 
        The $d_{gb}$(C) and $d_{gb}$(O) are distances from carbon and oxygen to gb plane ({\AA}); the negative and positive values shows different sides of gb plane.
        The $d_1$(C-O) and  $d_2$(C-O) are the first and second closest distances between C and O  ({\AA}). 
        """


    #distances and energies in hs443CO calculated using table_pairs(calc, conv, varset)
    dedCO    = [(2.57, 404.0), (2.99, -24.0), (3.8, -28.0), (4.67, -52.0), (5.09, -62.0), (5.53, -47.0)]#, (5.61, -15.0), (5.88, -47.0), (6.33, -1.0), (6.92, 6.0), (7.0, -23.0), (7.51, -21.0), (7.6, -3.0), (8.66, -5.0), (9.14, 0.0)]
    #The values of interaction for specific cell are given only for 2.99 and 3.8 distances; In other cases, the values from hs443 are used
    dedCO_T1 = [(2.57, 404.0), (2.99, -9),   (3.8, -43),  (4.67, -52.0), (5.09, -62.0), (5.53, -47.0)]
    dedCO_C1 = [(2.57, 404.0), (2.99, 17),   (3.8, 57),   (4.67, -52.0), (5.09, -62.0), (5.53, -47.0)]
    dedCO_T2 = [(2.57, 404.0), (2.99, -108), (3.8, -101), (4.67, -52.0), (5.09, -62.0), (5.53, -47.0)]
    dedCO_S7 = [(2.57, 404.0), (2.99, 48),   (3.8, -29),  (4.67, -52.0), (5.09, -62.0), (5.53, -47.0)]





    """ Produce plot and tables: co-segregation energies """
    data = []

    str, list = process_list(coseg_list_T1m, 't111gCvOvms', '93kp9', range(1,5), calc, conv, varset, run_res = False, bname = 'T1m', eCObulk =  dedCO_T1) 
    table1.extend(str     )
    table1.append('\hline')
    data.extend(list)

    str, list =  process_list(coseg_list_T1g, 't111sgCvOvms', '93kp9', range(1,5), calc, conv, varset, run_res = False, bname = 'T1g', eCObulk =  dedCO_T1 )
    table1.extend(  str   )
    table1.append('\hline')
    data.extend(list)
    

    str, list = process_list(coseg_list_S7, 'csl71sgCvOvms', '93', range(1,5), calc, conv, varset , run_res = False, bname = 'S7', eCObulk =  dedCO_S7) 
    table1.extend( str )
    data.extend(list)

    plot_coseg(data, 'coseg_T1', split = -5 )

    data = []

    str, list = process_list(coseg_list_C1, 'c1gCvOvms', '93kp7', range(1,6), calc, conv, varset, run_res = False, bname = 'C1', eCObulk =  dedCO_C1 )
    table2.extend(  str    )
    table2.append('\hline')
    data.extend(list)

    str, list = process_list(coseg_list_T2, 't21gCvOvms', '93', range(1,5), calc, conv, varset, run_res = False, bname = 'T2', eCObulk =  dedCO_T2 )
    table2.extend(   str   )
    # table2.append('\hline')
    data.extend(list)

    plot_coseg(data, 'coseg_C1_T2', split = -6)




    """ Produce tables """
    if 1:
        latex_table(table1, caption[0], label = 'tab:coseg1', header = header, fullpage = '*', filename = 'tab_coseg.tex', writetype = 'w')
        latex_table(table2, caption[1], label = 'tab:coseg2', header = header, fullpage = '*', filename = 'tab_coseg.tex', writetype = 'a')

    # print data
    return








def table_triple_and_tetrad(calc, conv, varset):



    def process(it,  ise,  verlist, nC, nO):
        """
        nC, nO - numbers of C and O atoms 
        """
        # struct_des['hs443CO2']   = des("H/CO2",   "obtained from template by adding carbon and two oxygen; v1-15 - different topo configurations; 99 atoms; ")
        # struct_des['hs443C2O2']   = des("H/C2O2",   "obtained from template by adding two carbon and two oxygen; v1-42 - different topo configurations; 100 atoms; ")
        #Formation of triple from atoms in ideal solid solution
        EoctC = -10.826305 #from hs443C.f - hs443
        EoctO = -10.51961 #from hs443O.f - hs443
        Ehs443_83 = -751.78701
        table = []
        #the sequence of atoms in typat is C O O C , but in reality C C O O, because of the error while converting Abinit structure to POSCAR
        for v in verlist:
            # res_loop(it,  '93', v, calc, conv, varset)
            cl = calc[(it, '93', v)]

            e3 =  (cl.energy_sigma0 - Ehs443_83 - nC*EoctC - nO*EoctO) * 1000 /(nC+nO)# Solution energy

            #Find distances d1, d2, d3
            st = cl.end
            assert st.nznucl[1] == nC 
            assert st.nznucl[2] == nO 

            for i, t in enumerate(st.typat):
                if t > 1:
                    i_imp1 = i #number of first impurity atom
                    break
            
            #Distances: C1-C2, C1-O1, C1-O2; based on distances each configuration can be described as a1+c, a2-c, where a1 and a2 are lattice constants a but distinguished by 120 degrees
            out = ''
            # x1 = st.xcart[i_imp1        ]
            # for j in range(nC+nO-1):
            #     x2 = st.xcart[i_imp1 + j + 1]
            #     dist, d2 = image_distance(x1, x2, st.rprimd, order = 3)
            #     out  += ' & {:5.2f}'.format(dist)


            # xcartC = st.xcart[i_imp1:i_imp1+nC]
            # xcartO = st.xcart[i_imp1+nC:i_imp1+nC+nO]
            xcartCO = st.xcart[i_imp1:i_imp1+nC+nO]
            if nC == 1: 
                typCO = ['C', 'O', 'O']
            elif nC == 2: 
                typCO = ['C', 'C', 'O', 'O']
            
            # print xcartC
            # print xcartO
            keylist = []
            bonds = []
            distlist = []
            for i, x1 in enumerate(xcartCO ):
                for j, x2 in enumerate(xcartCO):
                    if all(x1 == x2): continue 
                    dist, d2 = image_distance(x1, x2, st.rprimd, order = 3)
                    k = str(i)+str(j); ki = str(j)+str(i)
                    if k in keylist or ki in keylist: continue
                    keylist.append(k)
                    keylist.append(ki)
                    out  += ' & {:5.2f}'.format(dist)

                    distlist.append(dist)
                        
                    if   typCO[i]+typCO[j] == 'CC': ib = 0
                    elif typCO[i]+typCO[j] == 'OO': ib = 1
                    elif typCO[i]+typCO[j] == 'CO': ib = 2 
                    bonds.append(ib)
            # print bonds
            # print distlist


            #Take energies of independent bonds and sum them
            defCC_OO_CO = [(2.57, 136.0, 267.0, 208.0), (2.99, 3.0, 44.0, -5.0), (3.8, 15.0, 6.0, -8.0), 
            (4.67, 10.0, -34.0, -19.0), (5.09, -14.0, -27.0, -25.0), (5.53, -0.0, -20.0, -17.0), (5.61, 2.0, -2.0, -1.0), 
            (5.88, -9.0, -20.0, -17.0), (6.33, 9.0, 5.0, 6.0), (6.92, 17.0, 2.0, 9.0), (7.0, -1.0, -6.0, -5.0), 
            (7.51, 1.0, -4.0, -4.0), (7.6, 10.0, 4.0, 5.0), (8.66, 6.0, 4.0, 4.0), (9.14, 12.0, 4.0, 6.0)] #formation energies of independent CC, OO and CO pairs;
            ef_sum = 0
            for ibond, d in zip(bonds, distlist):
                mindr = 100
                for t in defCC_OO_CO:
                    dr = abs(t[0]-d)
                    if dr < mindr: 
                        mindr = dr
                        ef_pair_bulk = t[1+ibond]
                # print d, ibond, ef_pair_bulk
                ef_sum += ef_pair_bulk
            # print ef_sum



            geoname = cl.path["input_geo"]
            w = geoname.split('/')[-1].split('.')[0]
            conf = w.split('COv')[1]

            if len(conf) == 3:
                name = 'CO'+conf[0]+'-'+conf[1:3]
            elif len(conf) == 5:
                name = 'CC'+conf[0]+'-'+conf[1:3]+'-O'+conf[4] #Manually rearranged to take into account error of constructing POSCAR files

            name =  str(v)+" ("+name+')'

            table.append("{:25s} & {:5.0f} & {:5.0f}".format(name, e3, ef_sum/3)+out) # ef_sum/3 - formation energy of  separate pairs with the same types of bonds as in cluster


        return table

    caption = """ Formation energies $E_{form}$ (meV) of clusters. 
    Distances between first carbon atom and other impurity atoms define configuration of cluster.
    The first part of name in parenthesis, for example CO1, denotes C-O pair
    configuration given in table \\ref{tab:epairs}
    The equally numbered atoms in the second part of the name usually belongs to different positions in the bulk cell.
    """


    header = "Name & $E_{form}$ & d(C1-O1) & d(C1-O2) & d(O1-O2)"
    table1 = process('hs443CO2', '93', range(1,16), 1, 2, )
    latex_table(table1, caption, label = 'tab:CO2', header = header, filename = 'tab_clusters.tex')
    

    header = "Name            & $E_{form}$ & d(C1-C2) & d(C1-O1) & d(C1-O2) & d(C2-O1) & d(C2-O2) & d(O1-O2)"
    table2 = process('hs443C2O2', '93',range(1,43), 2, 2, )
    latex_table(table2, caption, label = 'tab:C2O2', header = header,  filename = 'tab_clusters.tex', writetype = 'a')






def table_convergence(calc, conv, varset):
    """Convergence according to the size of cells for gb energy, segregation energy and co-segregation"""


    pure_gb = [('t111g', '9292', 'T1m [2x1x2] '),  ('t112g','93kp9', 'T1m [3x1x2]'),]# ('t114g', '93kp9') ] - gives reduced gbe (336)- the problem may be due to the differet kpoint numbers for the bulk cell
    segreg1  = [
    ('t111gCi1',      '93kp9', 't111gCv'    ,  'T1m(1) C [2x1x2]'),     
    ('t113gCi1',      '93kp9', 't113gCv'    ,  'T1m(1) C [2x1x3]'),     
    ('t114gCi1',      '93kp9', 't114gCv'    ,  'T1m(1) C [2x1x4]'), 
    ('t115gCi1',      '93kp9', 't115gCv'    ,  'T1m(1) C [2x2x2]')
]

    segreg2 = [
    ('c1gCi1Ov'    , '93kp7',  'c1gCvOvms'      , 'C1(1) C'),
    ('c1gCi2Ov'    , '93kp7',  'c1gCvOvms'      , 'C1(2) C'),
    ('c1gOi2Cv'    , '93kp7',  'c1gCvOvms'      , 'C1(2) O'),
    ( 't111sgCi6Ov', '93kp9',  't111sgCvOvms'   , 'T1g(6) C'),
    ( 't111sgOi6Cv', '93kp9',  't111sgCvOvms'   , 'T1g(6) O'),
    ('csl71sgCi4Ov', '93'   ,  'csl71sgCvOvms'  , 'S7(4) C'),
    ('csl71sgOi4Cv', '93'   ,  'csl71sgCvOvms'  , 'S7(4) O'),
    ('t21gCi4Ov','93',         't21gCvOvms'   , 'T2(4) C'),
    ('t21gOi4Cv','93',         't21gCvOvms'   ,'T2(4) O' ),
    ]


    coseg   = [('t111g',         'T1m (1-1) C-O [2x1x2]'), ('t112g', 'T1m (1-1) C-O [3x1x2]'), ('t113g', 'T1m (1-1) C-O [2x1x3]')]
    ise = '93kp9'

    table = []
    for t in pure_gb:
        # res_loop(t[0],t[1],range(1,5),calc,conv,varset, 'gbep', ('t111b_r','8302',1))
        cl1 = calc[(t[0], t[1], 1)] #Fitted egb and vgb were stored here on previous step
        cl2 = calc[(t[0], t[1], 2)] #The closest to equillibrium; use r1 r2 and r3
        # print st2.rprimd[0][0]
        table1 = "{8:20s} & {0:5} & {1:5.2f} & {2:5.2f} & {3:5.2f} & {4:5.0f} & {5:5.0f} & {6:5s} & {7:5s}".format(
            cl2.natom, cl2.vlength[0], cl2.vlength[1], cl2.vlength[2], cl1.egb, cl1.vgb, '-', '-',  t[2], )
        table.append(table1)
    table.append('\hline')

    for it in segreg1:

        itC =  it[0]
        itCv = it[2]

        # res_loop(itC,ise,range(1,5),calc,conv,varset, 'e_seg', (itCv,))
        
        cl1 = calc[(itC, it[1], 1)] #Fitted egb and vgb were stored here on previous step
        cl2 = calc[(itC, it[1], 2)] #The closest to equillibrium; use r1 r2 and r3
        # print st2.rprimd[0][0]
        table1 =  "{8:20s} & {0:5} & {1:5.2f} & {2:5.2f} & {3:5.2f} & {4:5s} & {5:5s} & {6:5.0f} & {7:5.1f}".format(
             cl2.natom, cl2.vlength[0], cl2.vlength[1], cl2.vlength[2], '-', '-', cl1.e_seg, cl1.v_seg,  it[-1] )
        table.append(table1)
    
    table.append('\hline')


    for it in segreg2:
        for suffix in  [("", ' [1x1x1]'), (".r2d", ' [1x2x1]'), (".r3d", ' [1x1x2]')]:

            itC =  it[0]+suffix[0]
            itCv = it[2]+suffix[0]
            if suffix[0]:
                ise = '93'
            else:
                ise = it[1]

            # res_loop(itC,ise,range(1,5),calc,conv,varset, 'e_seg', (itCv,), plot = 1)
            
            cl1 = calc[(itC, ise, 1)] #Fitted egb and vgb were stored here on previous step
            cl2 = calc[(itC, ise, 2)] #The closest to equillibrium; use r1 r2 and r3
            # print st2.rprimd[0][0]
            try:
                table1 =  "{8:20s} & {0:5} & {1:5.2f} & {2:5.2f} & {3:5.2f} & {4:5s} & {5:5s} & {6:5.0f} & {7:5.1f}".format(
                 cl2.natom, cl2.vlength[0], cl2.vlength[1], cl2.vlength[2], '-', '-', cl1.e_seg, cl1.v_seg,  it[-1]+suffix[1] )
            except:
                continue

            table.append(table1)
        table.append('\hline')






    # table.append('\hline')

    ise = '93kp9'

    for it in coseg:
        itC = it[0]+'COi6.1-1is'
        itCv = it[0]+'CvOvms'
        # res_loop(itC,ise,range(1,5),calc,conv,varset, 'coseg', (itCv,))
        cl1 = calc[(itC, ise, 1)] #Fitted egb and vgb were stored here on previous step
        cl2 = calc[(itC, ise, 2)] #The closest to equillibrium; use r1 r2 and r3
        # print st2.rprimd[0][0]
        table1 =  "{8:20s} & {0:5} & {1:5.2f} & {2:5.2f} & {3:5.2f} & {4:5s} & {5:5s} & {6:5.0f} & {7:5.1f}".format(
             cl2.natom, cl2.vlength[0], cl2.vlength[1], cl2.vlength[2], '-', '-', cl1.e_seg, cl1.v_seg,  it[1] )
        table.append(table1)


    caption = """ Influence of cell sizes $L_1$, $L_2$, and $L_3$ ({\AA}) on GB energy $\gamma$ (mJ/m$^2$) and excess volume $\delta $  (m{\AA}),
    segregation energy $E_{seg}$ (or interaction energy $\Delta E_2$ in the case of co-segregation) (meV) and segregation volume $V_{seg}$ ($\Delta V_{2}$ )({\AA$^3$}) due to the segregation and co-segregation of carbon and oxygen.
    $N_{at}$ is number of atoms. Multipliers in brackets characterize relative cell sizes.
    The sizes are provided for the expansion \\textnumero 2, while the segregation energies and volumes are obtained from the approximation.
    """
    header = "Name &  $N_{at}$  &   $L_1$  &   $L_2$   &  $L_3$  &  $\gamma$  & $\delta $  &  $E_{seg}$  & $\Delta V$ "


    latex_table(table, caption, label = 'tab:conv', header = header, fullpage = '*', filename = 'tab_conv.tex')


def table_pairs_in_grains(calc, conv, varset):
    """
    C-O pairs in grains of cells with boundaries

    """
    coseg_list_T1m = [      
    # 't111gCOv1'  ,
    't111gCOv2'  ,
    # 't111gCOv3is',
    # 't111gCOv4'  ,
    't111gCOv5'  ,
    # 't111gCOv6'  ,
    ]
    coseg_list_T1g = [     
    't111sgCOv6',  
    't111sgCOv2'  
     ]


    coseg_list_C1 = [
    # 'c1gCOv1',
    'c1gCOv2',
    # 'c1gCOv3',
    # 'c1gCOv4', #almost the same as 'c1gCOv2' 
    # 'c1gCOv5',
    # 'c1gCOv6',
    'c1gCOv7',
    ]

    coseg_list_T2 = [ 
    't21gCOv2',
    't21gCOv6',
    ]

    coseg_list_S7 = [ 
    'csl71sgCOv5',
    'csl71sgCOv9',
    ]



    def process_list(it_coseglist, it_base, ise, verlist, calc, conv, varset, run_res = False, bname = None):
        #bname - boundary name
        #distances and energies in hs443CO calculated using table_pairs(calc, conv, varset)
        dedCO = [(2.57, 404.0), (2.99, -24.0), (3.8, -28.0), (4.67, -52.0), (5.09, -62.0), (5.53, -47.0)]#, (5.61, -15.0), (5.88, -47.0), (6.33, -1.0), (6.92, 6.0), (7.0, -23.0), (7.51, -21.0), (7.6, -3.0), (8.66, -5.0), (9.14, 0.0)]

        # run_res = True
        table, e_seg = [], []
        for it, dt in zip(it_coseglist, ['$a$', '$a+0.5c$'] ):
            if run_res:
                # add_loop(it,ise,verlist, calc, conv, varset,  'up1')

                res_loop(it,ise,verlist, calc, conv, varset,  'coseg', (it_base,ise,1))

            id =(it,ise, 2)

            st = calc[id].end
            gbpos2 = calc[id].gbpos 
            gbpos1 = gbpos2 - st.rprimd[0][0]/2.
            d1 = -(st.xcart[-2][0] - gbpos2)  #Here inverted in comparison to coseg table 
            d2 = -(st.xcart[-1][0] - gbpos2)

            dd1, dd2 = image_distance(st.xcart[-2], st.xcart[-1], st.rprimd, order = 3)

            name = bname +' ('+it.split('v')[1]+')'
            name = bname + ' '+dt
            table.append( 
                "{:25s} & {:5.0f} & {:5.1f} & {:5.2f} & {:5.2f} & {:5.2f} & {:5.2f}".format(
                name , calc[id].e_seg,  calc[id].v_seg, d1, d2, dd1, dd2) )

        return table

    table1, table2 = [], []
    header = "Name & $\Delta E$ & $\Delta V$ & $d_{gb}$(C) & $d_{gb}(O)$ & $d_1$(C-O) & $d_2$(C-O)"

    caption = """ Interaction of C and O atoms in bulk region of cells with T1, T2, C1, and S7 boundaries.  The $\Delta E = E( {\\rm CO }, bulk) - E( {\\rm C}-bulk, {\\rm O}-bulk)$.
    The $d_{gb}$(C) and $d_{gb}$(O) are distances from carbon and oxygen to gb plane;
    The $d_1$(C-O) and  $d_2$(C-O) are the first and second closest distances between C and O. 
    """

    table1.extend( process_list(coseg_list_T1m, 't111gCvOvms', '93kp9', range(1,5), calc, conv, varset, run_res = False, bname = 'T1m' )     )
    table1.append('\hline')
    # table1.extend( process_list(coseg_list_T1g, 't111sgCvOvms', '93kp9', range(1,5), calc, conv, varset, run_res = False, bname = 'T1g' )     )
    # table1.append('\hline')
    table1.extend( process_list(coseg_list_C1, 'c1gCvOvms', '93kp7', range(1,6), calc, conv, varset, run_res = False, bname = 'C1' )     )
    table1.append('\hline')
    table1.extend( process_list(coseg_list_T2, 't21gCvOvms', '93', range(1,5), calc, conv, varset, run_res = False, bname = 'T2' )     )
    table1.append('\hline')
    table1.extend( process_list(coseg_list_S7, 'csl71sgCvOvms', '93', range(1,5), calc, conv, varset , run_res = False, bname = 'S7')  )

    latex_table(table1, caption, label = 'tab:epairs_grains', header = header, fullpage = '*', filename = 'tab_epairs_grains.tex')

    return
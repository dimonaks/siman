
        if imp == 'C':
            it_p = 'hs443';  ise_p = '83'
            it_i = 'hs443C'; ise_i = '93'; it_e = 'hs443C.m'; ise_e = '83'
            # it_i = 'hs443C.r'; ise_i = '83dos'; it_e = 'hs443C.m'; ise_e = '83dos'
            e_at = -9.208 #graphite
            Vsol = 7 * 1e-30 #m^3

        elif imp == 'C_bigger':
            it_p = 'hs443';  ise_p = '83dos'
            it_i = 'hs443C'; ise_i = '83dos'; it_e = 'hs443C.m'; ise_e = '83dos'
            e_at = -9.208 #graphite
        

        elif imp == 'O':
            it_p = 'hs443';  ise_p = '83'
            # it_i = 'hs443O'; ise_i = '93'; it_e = 'hs443O.m'; ise_e = '83'
            it_i = 'hs443O.r'; ise_i = '83dos'; 
            it_e = 'hs443O.m'; ise_e = '83dos'; it_un = 'hs443O.ur'; ise_un = '83dos'
            e_at = -4.926 #O2
            Vsol = 4 * 1e-30 #m^3


    #Additional
    # calc[('hs443C.r','dos', 1)].bader_analysis()
    # calc[('hs443O.r','dos', 1)].bader_analysis()
    # calc[('hs221C.f','83dos', 1)].bader_analysis()
    # calc[('hs221C.f','83acdos', 1)].bader_analysis()
    
    # id1 = ('hs443O.r','dos', 1)
    # # id1 = ('hs443C','93', 100)
    # id2 = ('hs443O.m','83',  1)

    #dos analysis


    # plot_dos(calc[('hs221.f','dos',   1)])
    # plot_dos(calc[('hs221.f','83dos', 1)])

    # plot_dos(calc[('hs221O.fr','83dos', 1)])
    # plot_dos(calc[('hs221O.fr','dos', 1)])
    # plot_dos(calc[('hs221O.fr','dosf', 1)])
    # plot_dos(calc[('hs221O.fr','dosm', 1)])

    # plot_dos(calc[('hs221C.f','83acdos',   1)])
    # plot_dos(calc[('hs443O.r','dos',   1)])
    # plot_dos(calc[('hs443C.r','83dos',   1)], calc[('hs443C.r','83dos',   16)], "small_and_big")



# res_loop('t111gCvOvms',  '93kp9', range(1,5), calc, conv, varset)
# st = 
# st_r  = replic(st,   mul = (1,2,2), inv =  1 )
# st = replic(st_r, mul = (1,2,2), inv = -1 )
# write_xyz(calc[('c1gCi2Ov','93kp7', 2)].init, analysis = 'imp_surrounding'   )
# write_xyz(calc[('t111sgCi6Ov','93kp9', 2)].end, analysis = 'imp_surrounding', nnumber = 14   )


# write_xyz(calc[('hs443C.f','93', 1)].end, analysis = 'imp_surrounding',nnumber = 18   )
# write_xyz(calc[('hs221O.f','93', 1)].end,   replications = [2,2,2]  )

# thesis_presentation_pictures()

# )

# res_loop('t111g',  '9292', range(1,11), calc, conv, varset)
# res_loop('t111sg',  '93kp9', range(1,5), calc, conv, varset)
# res_loop('t21g',  '93', range(1,5), calc, conv, varset)
# res_loop('c1g',  '929', range(1,5), calc, conv, varset)




"""1) Section for inheritance""" 
for v in range(1,2):
    ver_new = v + 0
    #inherit_icalc('r1r2r3', it_new = 't113gC_template', ver_new = ver_new, id_base = ('t113gC_template','93kp9',2), calc = calc, id_from = geo_folder+'/T1/C/t113gC_template/target/t113gC_templatetarget.'+str(v)+'.in.geo')
    #inherit_icalc('full', 'hs221O.fr', ver_new, ('hs221O.f','93',v), calc)# 'hs221.f','83',1)
    #inherit_icalc('replace_atoms', 't111bO', ver_new, ('t111bC','93kp9',v), calc, atom_new = 'O', atom_to_replace = 'C')
    #inherit_icalc('replace_atoms', 'hs443CC', ver_new, ('hs443CO','93',v), calc, atom_new = 'C', atom_to_replace = 'O')
    #inherit_icalc('replace_atoms', 'hs443OO', ver_new, ('hs443CO','93',v), calc, atom_new = 'O', atom_to_replace = 'C')
    #inherit_icalc('remove_imp', 'hs443CO.m', ver_new, ('hs443CO','93',v), calc)
    # inherit_icalc('remove_imp', 'hs443C.m', ver_new, ('hs443C.f','93',v), calc)
    # inherit_icalc('remove_imp', 'hs443O.m', ver_new, ('hs443O.f','93',v), calc)


# inherit_icalc('remove_imp', 'csl77gCi6'+'.m', 1, ('csl77gCi6', '93kp9', 1) , calc)
# inherit_icalc('remove_imp', 'csl77gCv'+'.m', 1, ('csl77gCv', '93kp9', 1) , calc)
    # inherit_icalc('replace_atoms', 't21bO', v, ('t21bC', '93', v) , calc, atom_new = 'O', atom_to_replace = 'C')
    # inherit_icalc('replace_atoms', 'c1gBi1', v, ('c1gCi1', '93kp7', v) , calc, atom_new = 'B', atom_to_replace = 'C')
    # inherit_icalc('replace_atoms', 'c1gBv', v,  ('c1gCv', '93kp7', v) , calc, atom_new = 'B', atom_to_replace = 'C')
    # inherit_icalc('replace_atoms', 'hs554O', v, ('hs554C', '93', v) , calc, atom_new = 'O', atom_to_replace = 'C')
    # inherit_icalc('full', 'hs443C.r', v, ('hs443C', '93', v), calc, '', '', 0 )   #at 2014-01-21
    # inherit_icalc('full',       'hs443C.r', v, ('hs443C', '93', v), calc)
    # inherit_icalc('full',       'hs443O.r', v, ('hs443O', '93', v), calc)
    # inherit_icalc('remove_imp', 'hs443C.m', v, ('hs443C', '93', v), calc)
    # inherit_icalc('remove_imp', 'hs443O.m', v, ('hs443O', '93', v), calc)

# inherit_icalc('full', 'hs443C.fr', 1, ('hs443C.f', '93', 1), calc)  
# inherit_icalc('full', 'hs443O.fr', 1, ('hs443O.f', '93', 1), calc)  

# inherit_icalc('full', 'hs221O.fr', v, ('hs221O.f', '93', v), calc, '', '', 0 )   #at 2014-01-21

# inherit_icalc('r1r2r3', it_new = 'hs443C', ver_new = 100, id_base = ('hs443C.f','93',1), calc = calc, id_from = ('hs443.f', '83', 1) )



"""2) Section for creation of the structures with impurities"""

# add_impurity('t114gCi1', 'C', 'gb', calc, 0.56, find_close_to = (0.77, 0.6, 0.25))    #on 2014-01-24    #on 2014-02-08
# add_impurity('t114gCv', 'C', 'grain_vol', calc, 0.56,find_close_to = (0.51, 0.65, 0.25))
# add_impurity('c1gCOi10.1', 'O', 'gb', calc, 0.25,  find_close_to = (0.740608513355, 0.480618625879, 0.278868079185))  

#add_impurity('t112gCOi6.1-1is', 'C', 'gb', calc, 0.45, find_close_to = (0.74, 0.585767209530, -0.000012566336))  
# add_impurity('t112gCOi6.1-1is', 'O', 'gb', calc, 0.45, find_close_to = (0.74, 0.585767209530,  0.5))  

# add_impurity('t112gCvOvms', 'C', 'grain_vol', calc, 0.50, find_close_to = (0.481233566999, 0.617050468922, 0.249992430210) )    
# add_impurity('t112gCvOvms', 'O', 'grain_vol', calc, 0.50, find_close_to = (0.481233566999 + 0.5, 0.617050468922, 0.249992430210 ))    
# add_impurity('csl71bC', 'C', 'central', calc, 0.50)



#st = calc[('hs443','83',20)].init


#find_pairs('hs443CO_test', gb4_geo_folder+'hps443CO_template/grainA_s/hps443CO_templategrainA_s.1.in.geo', central_atoms = 57  )

#calc_ac( calc[('t111gCOv2','93kp9',1)].hex_a, calc[('t111gCOv2','93kp9',1)].hex_c, calc[('t111b_r','8302',1)].hex_a, calc[('t111b_r','8302',1)].hex_c , type = "double_cell")

# calc_ac( calc[('t111bC.f','93kp9',1)].hex_a, calc[('t111bC.f','93kp9',1)].hex_c, calc[('t111b_r','8302',1)].hex_a, calc[('t111b_r','8302',1)].hex_c , type = "double_cell")

# calc_ac( calc[('hs443CO','93',1)].hex_a, calc[('hs443CO','93',1)].hex_c, calc[('hs443O.f','93',1)].hex_a, calc[('hs443O.f','93',1)].hex_c, calc[('hs443','83',20)].hex_a, calc[('hs443','83',20)].hex_c, type = "two_atoms")
# calc_ac( calc[('hs443CO2','93',1)].hex_a, calc[('hs443CO2','93',1)].hex_c, calc[('hs443C.f','93',1)].hex_a, calc[('hs443C.f','93',1)].hex_c, calc[('hs443','83',20)].hex_a, calc[('hs443','83',20)].hex_c, type = "two_atoms")

#calc_ac( calc[('t111bC.f','93kp9',1)].hex_a, calc[('t111bC.f','93kp9',1)].hex_c, calc[('t111bO.f','93kp9',1)].hex_a, calc[('t111bO.f','93kp9',1)].hex_c, calc[('t111b_r','8302',1)].hex_a, calc[('t111b_r','8302',1)].hex_c , type = "two_atoms")

# calc_ac( calc[('t21bC.f','93',1)].hex_a, calc[('t21bC.f','93',1)].hex_c, calc[('t21bO','93',1)].hex_a, calc[('t21bO','93',1)].hex_c, calc[('t21b_r','83',1)].hex_a, calc[('t21b_r','83',1)].hex_c , type = "two_atoms")

# print calc[('t21gCi1','93',1)].end.vol
# print calc[('t21gOi1','93',1)].end.vol
# print calc[('t21g','93',1)].end.vol

# print calc[('t21bC.f','93',1)].end.vol - calc[('t21b_r','83',1)].end.vol
# print calc[('t21bO','93',20)].end.vol - calc[('t21b_r','83',1)].end.vol
# # print calc[('hs443CO','93',1)].end.vol
# # print calc[('hs443O.f','93',1)].end.vol - calc[('hs443','83',20)].end.vol
# # print calc[('hs443C.f','93',1)].end.vol - calc[('hs443','83',20)].end.vol

# #print calc[('t111bC.f','93kp9',1)].end.vol - calc[('t111b_r','8302',1)].end.vol
# #print calc[('t111bO.f','93kp9',1)].end.vol - calc[('t111b_r','8302',1)].end.vol
# print calc[('t111gCi1','93kp9',1)].end.vol
# print calc[('t111gOi1','93kp9',1)].end.vol
# print calc[('t111g','9292',6)].end.vol


# add_loop('csl71sgCOv5','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgCOv9','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# add_loop('csl71sgOCi6.5-4','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgOCi7.1-6','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgCOi1.6-4','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgCOi11.4-4','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgCOi22.1-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# add_loop('t111sgCOi17.7-4','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgCOi28.4-5','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgCOi29.4-6','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgCOi30.7-5','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgOCi8.4-7','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgOCi9.6-2','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgOCi15.6-6','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgOCi22.6-3','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgOCi25.6-7','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgOCi27.6-2','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')


# add_loop('t111sgCOv2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgCOv6','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')






# add_loop('t21gCvOvms','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# add_loop('t21gCi1Ov','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCi2Ov','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCi3Ov','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCi4Ov','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOi1Cv','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOi2Cv','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOi3Cv','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOi4Cv','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')



# add_loop('hs443CO.m','83',[4,11],calc,conv,varset, 'up1')
# add_loop('hs443C.m','83',1,calc,conv,varset, 'up1')
# add_loop('hs443O.m','83',1,calc,conv,varset, 'up1')

# add_loop('csl71sgCvOvms','93',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# add_loop('csl71sgCi1Ov','93',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgCi2Ov','93',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgCi3Ov','93',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgCi4Ov','93',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgCi5Ov','93',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgCi6Ov','93',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgOi1Cv','93',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgOi2Cv','93',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgOi3Cv','93',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgOi4Cv','93',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgOi5Cv','93',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sgOi6Cv','93',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')



# add_loop('t111sgCvOvms','93kp9',range(5,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# add_loop('t111sgCi1Ov','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgCi2Ov','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgCi3Ov','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgCi4Ov','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgCi5Ov','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgCi6Ov','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgCi7Ov','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgOi1Cv','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgOi2Cv','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgOi3Cv','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgOi4Cv','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgOi5Cv','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgOi6Cv','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111sgOi7Cv','93kp9',range(4,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')



# add_loop('hs554C','93',range(12,17),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# add_loop('t114gCi1','93kp9',range(2,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t114gCv','93kp9',range(2,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# add_loop('hs554','83',range(1,17),calc,conv,varset, 'up1')

# add_loop('hs554','83',1,calc,conv,varset, 'up1')
# add_loop('hs554C','93kp9',1,calc,conv,varset, 'up1')
# add_loop('hs554O','93kp9',1,calc,conv,varset, 'up1')

# add_loop('hs443CO.m','83',range(1,16),calc,conv,varset, 'up1')
# add_loop('t111sg','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')




# add_loop('t21gCOv2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOv6','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCvOvms','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')


# add_loop('t21gCOi1.4-3','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi2.3-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi3.2-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi4.2-1','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi5.3-1','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi6.4-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi7.2-1','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi8.4-1','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi9.4-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi10.3-3','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi11.4-3','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi12.3-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi13.4-1','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi14.1-4','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi15.3-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCOi16.4-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi1.4-3','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi2.3-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi3.2-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi4.2-1','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi5.3-1','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi6.4-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi7.2-1','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi8.4-1','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi9.4-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi10.3-3','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi11.4-3','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi12.3-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi13.4-1','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi14.1-4','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi15.3-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOCi16.4-2','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')


# add_loop('t21gCi1Ov','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCi2Ov','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCi3Ov','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gCi4Ov','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOi1Cv','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOi2Cv','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOi3Cv','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t21gOi4Cv','93',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')





# add_loop('csl71sg5','93',range(1,6)[::-1],calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sg10','93',range(1,6)[::-1],calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sg15','93',range(1,6)[::-1],calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('csl71sg10z2y','93',range(1,6)[::-1],calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# write_xyz(calc[('csl71sg15','93', 5)].init)

# add_loop('t114g','93kp9',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# add_loop('hs443C2O2','93',range(1,43),calc,conv,varset, 'up1')

# add_loop('c1gCOi10.1','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

#add_loop('c1gCi1Ov','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gCi2Ov','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gOi1Cv','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gOi2Cv','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# add_loop('c1gCOi1.2-1','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gCOi2.2-1','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gCOi3.1-2','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gCOi4.2-1','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gCOi5.2-1','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gCOi6.1-2','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gCOi7.1-1','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gCOi8.2-2','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gOCi1.2-1','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gOCi2.2-1','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gOCi3.1-2','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gOCi4.2-1','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gOCi5.2-1','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gOCi6.1-2','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gOCi7.1-1','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gOCi8.2-2','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')


# add_loop('c1gCOv2','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gCOv4','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gCOv7','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1gCvOvms','93kp7',range(1,6),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')


#add_loop('hs443CO','93',2, calc, conv, varset, 'up1')
#add_loop('hs443OO','93ns',range(1,16), calc, conv, varset, 'up1')
#add_loop('hs443CO','93hd',[1, 2, 5, 15], calc, conv, varset, 'up1')

# add_loop('t111gCOi1.4-3','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCOi2.2-1','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCOi3.4-1','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCOi4.4-4is','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCOi5.2-2is','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# add_loop('t112gCOi6.1-1is','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# add_loop('t111gCOi7.3-3is','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCOi8.2-3','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCOi9.4-4ms','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCOi10.2-2ms','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCOi11.4-3','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCOi12.1-3','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCOi13.4-2','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCOi14.2-1','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCOi15.4-4is','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCOi16.2-2is','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gOCi1.4-3','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gOCi2.2-1','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gOCi3.4-1','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gOCi8.2-3','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gOCi11.4-3','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gOCi12.1-3','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gOCi13.4-2','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gOCi14.2-1','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# add_loop('t111gCi1Ov','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCi2Ov','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCi3Ov','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gCi4Ov','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gOi1Cv','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gOi2Cv','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gOi3Cv','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('t111gOi4Cv','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')

# add_loop('t112gCvOvms','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')


# res_loop('c1gBi1', '93kp7', range(1,6), calc, conv, varset, 'up1')    
# res_loop('c1gBv', '93kp7', range(1,6), calc, conv, varset, 'up1')    

# res_loop('c1gBi1', '93kp7', range(1,6), calc, conv, varset, 'e_seg', ('c1gBv', '93kp7', 1))    

# res_loop('c1gCv', '93kp7', range(1,6), calc, conv, varset)    
# res_loop('c1gCi1', '93kp7', range(1,6), calc, conv, varset, 'e_seg', ('c1gCv', '93kp7', 1))    
# res_loop('c1gCi1Ov', '93kp7', range(1,5), calc, conv, varset, 'e_seg', ('c1gCvOvms', '93kp7', 1))    

# res_loop('Ti2C227'  ,  '93',1 , calc, conv, varset, 'up1',)
# res_loop('Ti2O227'  ,  '93',1 , calc, conv, varset, 'up1',)
# res_loop('Ti2O227'  ,  '23',1 , calc, conv, varset, 'up1',)
# res_loop('Ti2O164'  ,  '23',1 , calc, conv, varset, 'up1',)
# res_loop('Ti2C164'  ,  '93',1 , calc, conv, varset, 'up1',)

# add_loop('t111sgCi6Ov.r'  ,  'lobster',2 , calc, conv, varset, 'up1', savefile = 'allw')
# add_loop('t111sgOi6Cv.r'  ,  'lobster',2 , calc, conv, varset, 'up1',savefile = 'allw')

# add_loop('t111sgCi6Ov.m'  ,  'lobster',2 , calc, conv, varset, 'up1',savefile = 'allw')

# # add_loop('t111gCi2Ov.r'  ,  'lobster',2 , calc, conv, varset, 'up1', savefile = 'allw')
# add_loop('t111gCi2Ov.m'  ,  'lobster',2 , calc, conv, varset, 'up1',savefile = 'allw')
# res_loop('hs443C'  ,  '93',[38,39,40,41] , calc, conv, varset, 'up1',)
# res_loop('hs443O'  ,  '93',[38,39,40,41] , calc, conv, varset, 'up1',)

# add_loop('hs554O',  '93', range(12,17) , calc, conv, varset, 'up1', inherit_option = 'inherit_xred', savefile = "min")
# res_loop('hs443.f',  '83dos', 1, calc, conv, varset, 'up1',)
# res_loop('hs443C.fr','83dos', 1, calc, conv, varset, 'up1',)
# res_loop('hs443O.fr','83dos', 1, calc, conv, varset, 'up1',)



# res_loop('hs443'  ,  '83dos', [1,6,11,16], calc, conv, varset, 'up1',)

# res_loop('hs443C.r'  ,  '83dos', [38,39,40,41], calc, conv, varset, 'up1',)
# res_loop('hs443O.r'  ,  '83dos', [38,39,40,41], calc, conv, varset, 'up1',)
# res_loop('hs443C.m'  ,  '83dos', [38,39,40,41], calc, conv, varset, 'up1',)
# res_loop('hs443O.m'  ,  '83dos', [38,39,40,41], calc, conv, varset, 'up1',)




# run_dos_segregations()

# res_loop('hs443O.ur'  ,  '83dos',[1,6,11,16] , calc, conv, varset, 'up1', )

# add_loop('hs221O.fr',    'dosm',  1 , calc, conv, varset,  'up1')
# add_loop('hs221O.fr',    'dos',  1 , calc, conv, varset,  'up1')
# add_loop('hs221O.fr',    'dosf',  1 , calc, conv, varset,  'up1')


# add_loop('hs221.f',    'dos',  1 , calc, conv, varset,  'up1')
# add_loop('hs221.f',    '83dos',  1 , calc, conv, varset,  'up1')

# add_loop('hs221O.f',    'dosh',  1 , calc, conv, varset,  'up1')
# # add_loop('hs221C.f',    '83acdos',  1 , calc, conv, varset,  'up1')

# res_loop('hs443C.r',    '83dos',[1,6,11,16], calc, conv, varset,  'up1', voronoi = True)
# res_loop('hs443C.m'  ,  '83dos',[1,6,11,16] , calc, conv, varset, 'up1',)
# res_loop('hs443O.m'  ,  '83dos',[1,6,11,16] , calc, conv, varset, 'up1')
# add_loop('hs443C.f',  '93f', 1, calc, conv, varset, 'up1')
# add_loop('hs443C',  '93', 100, calc, conv, varset, 'up1')

# add_loop('t111g', '93kp9p4', range(1,5), calc, conv, varset,"up1" , inherit_option = 'inherit_xred')    
# add_loop('t111b_r', '93kp9p4', 1, calc, conv, varset, "up1" )  

# add_loop('t21bO.f','93',1, calc, conv, varset, 'up1' )
# add_loop('c1bO.f','93kp7',1, calc, conv, varset,'up1')

# add_loop('t21bO','93',range(1,17), calc, conv, varset, 'up1', inherit_option = 'inherit_xred')
# add_loop('c1bO','93kp7',range(1,17), calc, conv, varset,'up1' , inherit_option = 'inherit_xred')

# add_loop('csl77gCv.m','83kp9',1,calc,conv,varset,'up1' )
# add_loop('csl77gCi6.m','83kp9',1,calc,conv,varset,'up1' )

# add_loop('hs554.f', '83', 1, calc, conv, varset, 'up1' )
# add_loop('gr221','83is1',1, calc, conv, varset,'up1' )
# add_loop('TiC','83is1lf',1, calc, conv, varset,'up1' )
# add_loop('TiC','83lf',1, calc, conv, varset,'up1' )
# add_loop('O','8at',1, calc, conv, varset,'up1', coord = 'car' )
# add_loop('OO','9dim',1, calc, conv, varset,'up1', coord = 'car' )

# add_loop('hs443O.m','83',range(1,17), calc, conv, varset, 'up1')







"""4) Section for reading results"""





# res_loop('t111g', '93kp9p4', range(1,5), calc, conv, varset, 'gbep', ('t111b_r', '93kp9p4', 1))    
# res_loop('t111b_r', '93kp9p4', 1, calc, conv, varset, )  


# res_loop('t21bO','93',range(1,17), calc, conv, varset,  )
# res_loop('c1bO','93kp7',range(1,17), calc, conv, varset, 'fit_ac')
# res_loop('c1bO','93kp7',range(1,17), calc, conv, varset,)
# res_loop('hs443O.m','83',range(1,17), calc, conv, varset, )
# res_loop('hs443O','93',range(1,17), calc, conv, varset, voronoi = True)

# res_loop('csl77gCv.m','83kp9',1,calc,conv,varset )
# res_loop('csl77gCv','93kp9',1,calc,conv,varset )
# res_loop('csl77gCi6.m','83kp9',1,calc,conv,varset )

# res_loop('t21bC','93',range(1,17), calc, conv, varset, )
# res_loop('c1bC','93kp7',range(1,17), calc, conv, varset, )



# calc_ac( calc[('hs443O.f','93',1)].hex_a, calc[('hs443O.f','93',1)].hex_c, calc[('hs443','83',20)].hex_a, calc[('hs443','83',20)].hex_c , type = "double_cell")    #on 2014-04-01

# print calc[('t114g','93kp9',1)].end.vol

# for i in range(1,6):
    # write_xyz( calc[('t21g','93',i)].end )
# print calc[('t112g','93kp9',1)].gbpos



# # add_loop('csl71g','93',    range(3,6),calc,conv,    varset,'up1')
# # res_loop('csl71g','93',    range(3,6),calc,conv,    varset,)

# temp_conc_simple()
# add_loop('c1g',     '929',    range(1,11), calc, conv, varset, 'up1')
# add_loop('csl71g',  '93',     range(3,7), calc, conv, varset, 'up1')
# add_loop('hs221C.f',  '93', 1, calc, conv, varset, 'up1')
# add_loop('hs221.f',  '83', 1, calc, conv, varset, 'up1')
# res_loop('hs221.f',  '83', 1, calc, conv, varset)
# add_loop('hs332C.f',  '93', 1, calc, conv, varset, 'up1')
# add_loop('hs332.f',  '83', 1, calc, conv, varset, 'up1')
# res_loop('hs443.f',  '83', 1, calc, conv, varset)
# res_loop('c1bO',  '93kp7', 1, calc, conv, varset)
# res_loop('c1bO',  '93kp7', 20, calc, conv, varset)


# res_loop('t111sgCi6Ov',  '93kp9', range(1,5), calc, conv, varset, 'coseg', ('t111sgCvOvms',) )
# res_loop('t111gCi3Ov',  '93kp9', range(1,5), calc, conv, varset, 'coseg', ('t111gCvOvms',) )
# res_loop('t111gCi1Ov',  '93kp9', range(1,5), calc, conv, varset, 'coseg', ('t111gCvOvms',) )
# res_loop('t21gCi4Ov',  '93', range(1,5), calc, conv, varset, 'coseg', ('t21gCvOvms',) )
# res_loop('c1gCi2Ov',  '93kp7', range(1,5), calc, conv, varset, 'coseg', ('c1gCvOvms',) )
# res_loop('t111sgCvOvms','93kp9',range(1,6), calc,conv,varset, 'e_seg', ('t111sgCvOvms',) , voronoi = False , readfiles = True)


# ('c1bO',  '93kp7',   1,    ('c1b_r', '830', 1)   )
# add_loop('csl71bC', '83kp8', 1, calc,conv,varset, 'up1')
# print calc[('hs554','83',1)].end.vol
# res_loop('hs332C.f',  '93', 1, calc, conv, varset, 'e_imp', ('hs332.f',  '83', 1) ,r_id = ('gr221', '83', 1))

# res_loop('hs554','83',range(1,17),calc,conv,varset, 'fit_ac')
# res_loop('hs554',  '83', 20, calc, conv, varset)
# print calc[('hs554','83',20)].energy_sigma0 #-  calc[('hs554.f','83',1)].energy_sigma0
# print calc[('hs554.f','83',1)].energy_sigma0

# print calc[('hs443C.f','93',1)].energy_sigma0 -  calc[('hs443.f','83',1)].energy_sigma0
# print calc[('hs443O.f','93',1)].energy_sigma0 -  calc[('hs443.f','83',1)].energy_sigma0
# print calc[('hs443.f','83',1)].energy_sigma0

# print calc[('hs554C.f','93',1)].energy_sigma0  -  calc[('hs554','83',20)].energy_sigma0
# print calc[('hs554C.f','93',1)].energy_sigma0 -  calc[('hs554.f','83',1)].energy_sigma0
# print calc[('csl77gCi6.m','83kp9',1)].energy_sigma0 -  calc[('csl77gCv.m','83kp9',1)].energy_sigma0

# print calc[('csl71sgCi6Ov.m','83',2)].energy_sigma0 -  calc[('csl71sgCvOvms.m','83',2)].energy_sigma0
# print calc[('csl71sgCi6Ov.m','83',2)].energy_sigma0 -  calc[('csl71sg10','93',2)].energy_sigma0
# print calc[('csl71sgCvOvms.m','83',2)].energy_sigma0 -  calc[('csl71sg10','93',2)].energy_sigma0
# res_loop('csl71sgCi6Ov.m','83',2,calc,conv,varset, )
# res_loop('csl71sgCi6Ov','93',2,calc,conv,varset, )
# res_loop('csl77gCv','93kp9',1,calc,conv,varset, )
# res_loop('csl77gCi6','93kp9',1,calc,conv,varset, )




# res_loop('hs443CO','93ns',2,calc,conv,varset,)
# add_loop('csl71b_r','83kp8',1,calc,conv,varset, 'up1' )
# res_loop('csl71b_r','83kp8',1,calc,conv,varset, 'up1' )

# res_loop('O','8at',1, calc, conv, varset)
# res_loop('gr221','83',1, calc, conv, varset )
# res_loop('gr221','83is1',1, calc, conv, varset )
# res_loop('OO','9dim',1, calc, conv, varset, 'dimer')

# res_loop('t111sgCvOvms','93kp9',    range(1,6),calc,conv,    varset,)


# res_loop('TiC','83',1, calc, conv, varset)
# res_loop('TiC','83lf',1, calc, conv, varset,)
# res_loop('TiC','83is1',1, calc, conv, varset,)
# res_loop('TiC','83is1lf',1, calc, conv, varset )










# res_loop('t111gCi1','93kp9',range(1,2),calc,conv,varset, 'e_seg',('t111gCv','93kp9',1))
# res_loop('t111gCi1Ov','93kp9',1,calc,conv,varset)
# res_loop('c1gCi1Ov','93kp7',1,calc,conv,varset)
# add_loop('csl71b_r','83kp8',1,calc,conv,varset, 'up1')
# res_loop('csl71b_r','83kp8',1,calc,conv,varset)


# res_loop('hs443CO.m','83',range(1,16),calc,conv,varset, 'e_2imp', ('hs443CO.m','83',15))

# res_loop('hs443','83',range(1,17),calc,conv,varset, 'fit_ac')

# res_loop('t114gCv','93kp9',range(1,5),calc,conv,varset)
# res_loop('t114gCi1','93kp9',range(1,5),calc,conv,varset, 'e_seg',('t114gCv','93kp9',1))
# res_loop('t115gCi1','93kp9',range(1,5),calc,conv,varset)
# res_loop('hs554','83kp9',1,calc,conv,varset)
# res_loop('hs554','83',1,calc,conv,varset,)

# print calc[('hs443','83',20)].energy_sigma0 
# print calc[('hs443.f','83',1)].end.vol/96/1000

# # print calc[('hs554','83kp9',1)].energy_sigma0 
# res_loop('hs443','83',20,calc,conv,varset)
# res_loop('hs443.f','83',1,calc,conv,varset,)

# res_loop('hs554C.f','93',1,calc,conv,varset, 'e_imp', ('hs554','83',1))
# res_loop('hs554O','93kp9',1,calc,conv,varset, 'e_imp', ('hs554','83kp9',1))

# res_loop('t111sg','93kp9',range(1,5),calc,conv,varset, 'gbep',('t111b_r','8302',1))

# res_loop('c1gCOi10.1','93kp7',range(1,6),calc,conv,varset)


# res_loop('csl71sg5','93',    range(1,6),calc,conv,    varset, 'gbep',('csl71b_r','83kp8',1))
# res_loop('csl71sg10','93',   range(1,6),calc,conv,   varset, 'gbep',('csl71b_r','83kp8',1))
# res_loop('csl71sg15','93',   range(1,6),calc,conv,   varset, 'gbep',('csl71b_r','83kp8',1))
# res_loop('csl71sg10z2y','93',range(1,6),calc,conv,varset, 'gbep',('csl71b_r','83kp8',1))



# for i in range(858155,858165):
#     print i, ', ',








# res_loop('t112g','93kp9',range(1,5),calc,conv,varset, 'gbep', ('t111b_r','8302',1))
#plot_conv( conv['hs443OO'], calc,  'e_2imp', [conv['hs443CO'],conv['hs443CC']])

#add_loop('hs221','83',1, calc, conv, varset, 'up1')
#res_loop('hs221','83',1, calc, conv, varset)
# write_xyz( replic( calc[('hs221C', '93', 1)].end, (3,3,3)) )
# write_xyz( replic( calc[('hs443CC', '93ns', 8)].end, (2,2,2)) )


#res_loop('t111gCi2Ov','93kp9',range(1,5),calc,conv,varset)
#res_loop('t111gCi3Ov','93kp9',range(1,5),calc,conv,varset)
#res_loop('t111gCi4Ov','93kp9',range(1,5),calc,conv,varset)

#res_loop('t111gOi1Cv','93kp9',range(1,5),calc,conv,varset)


#res_loop('t111gCvOvms','93kp9',range(1,5),calc,conv,varset)

#res_loop('t111gCOv2','93kp9',range(1,5),calc,conv,varset)
#res_loop('t111gCOv5','93kp9',range(1,5),calc,conv,varset)

#res_loop('t21gOi1', '94',    2, calc, conv, varset, "up1")   
# res_loop('t21gOi1', '93fe1D',    2, calc, conv, varset, "up")    
# res_loop('t21gOi1', '93fe1MD',     2, calc, conv, varset, "up") 
# res_loop('t21gOi1', '93',       1, calc, conv, varset) 
# res_loop('t21gOi1', '93fe1',    2, calc, conv, varset) 
# res_loop('t21gOi1', '93acc',    2, calc, conv, varset) 
# res_loop('t21gOi1', '93PT1',    2, calc, conv, varset) 
# res_loop('t21gOi1', '93CG',     2, calc, conv, varset) 


# res_loop('hs443CO','93ns',range(1,16), calc, conv, varset, 'e_2imp', ('hs443CO','93ns',15))
# res_loop('hs443CC','93ns',range(1,16), calc, conv, varset, 'e_2imp', ('hs443CC','93ns',15))
# res_loop('hs443OO','93ns',range(1,16), calc, conv, varset, 'e_2imp', ('hs443OO','93ns',15))




#add_loop('t111gCOv2','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
#add_loop('t111gCOv5','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')


#add_loop('t111gCvOvms','93kp9',range(1,5),calc,conv,varset, 'up1', inherit_option = 'inherit_xred')
#add_loop('csl77gCv','93kp9',1,calc,conv,varset, 'up1') #at 2014-01-24
#add_loop('csl77gCi6','93kp9',1,calc,conv,varset, 'up1') #at 2014-01-24
#write_xyz(calc[('csl77gCi6','93kp9',1)].init )
#write_xyz(calc[('csl77gCv','93kp9',1)].init )

#add_loop('t115gCi1','93kp9',[1,2,3,4],calc,conv,varset, 'up1') #at 2014-01-24
#res_loop('t113gC_template','93kp9',[2],calc,conv,varset )

#add_loop('hs443CO',  '93',    range(1,16), calc, conv, varset, "up1") 
# add_loop('t21gOv',  '93',    range(1,6), calc, conv, varset, "up1") 




# add_loop('c1gOv',   '93kp7', range(1,6), calc, conv, varset, "up1") 
# add_loop('c1gOi1',  '93kp7', range(1,6), calc, conv, varset, "up1") 
# add_loop('t111gOv', '93kp9', range(1,6), calc, conv, varset, "up1") 
# add_loop('t111gOi1','93kp9', range(1,6), calc, conv, varset, "up1") 

#del header.history[130:]
#add_loop('t111b_ml','8301',1,calc,conv,varset,'up1')
#calc[('t111b_ml','8301',1)].path["input_geo"] = "geo/test2"
#calc[('t111b_ml','8301',1)].write_geometry('init',override = True)
#add_loop('t111bO','93kp9',[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],calc,conv,varset,"up1",)



#add_loop('t111g','9292',1, calc, conv, varset, 'up')#, 'kpoint_conv')

#add_loop('t113gCi1','93kp9',2, calc, conv, varset, 'up1')
#add_loop('t111gCv','93kp9',12, calc, conv, varset, 'up1')

#add_loop('c1gCi1','93kp7',4, calc, conv, varset, 'up1')

#add_loop('hs443O','93',range(11,17), calc, conv, varset,'up1' )




#res_loop('csl77gCv','93kp9',1,calc,conv,varset) #at 2014-01-24
#res_loop('csl77gCi6','93kp9',1,calc,conv,varset) #at 2014-01-24
#res_loop('csl77gCi6','93kp9',1, calc, conv, varset,'e_seg', ('csl77gCv',) )
#write_xyz(replic(calc[('csl77gCv','93kp9',1 )].end, (2,1,2)) )

#res_loop('t111bO','93kp9',range(1,17), calc, conv, varset,'fit_ac') 
#res_loop('t111b_ml','8301',1,calc,conv,varset) 
#res_loop('hs221O.f','93',1, calc, conv, varset, 'e_imp_ts', ('hs221.f','83',1))
#for t in 'is1',:#'kp1','kp2':#,'','ec1','ec2', :#
#    res_loop('hs221O.fr','83'+t,1, calc, conv, varset, 'e_imp_ts', ('hs221.f','83'+t,1))

#res_loop('c1gOi1','93kp7',range(1,2), calc, conv, varset) 
#res_loop('t21gOi1','93',1, calc, conv, varset) 
#res_loop('t115gCv','93kp9',range(1,5), calc, conv, varset) 
#res_loop('t111gOv','93kp9',range(1,6), calc, conv, varset) 


#res_loop('c1gCi1','93kp7',range(1,6), calc, conv, varset,'e_seg', ('c1gCv',) )
#res_loop('t21gCi1','93',range(1,6), calc, conv, varset,'e_seg', ('t21gCv',) )
#res_loop('t111gOi1','93kp9',range(1,6), calc, conv, varset,'e_seg', ('t111gOv',) )
#res_loop('t111gCi1','93kp9acc',4, calc, conv, varset,'e_seg', ('t111gCv',) )
#res_loop('t111gCi1','93kp9',12, calc, conv, varset,'e_seg', ('t111gCv',) )
#res_loop('t115gCi1','93kp9',range(1,5), calc, conv, varset,'e_seg', ('t115gCv',) )
#res_loop('t111g','9292',2, calc, conv, varset)

#write_xyz(  replic( calc[('t111gCi1','93kp9',2)].end, (1,2,2)  ), repeat = 1  )

#a = replic( calc[('c1bC.f','93kp7',1)].end, 2)
#add_loop('t21g','93',range(1,6), calc, conv, varset,'up1') 
#res_loop('t21g','93',range(1,6), calc, conv, varset) 
#a = replic( calc[('t21g','93',3)].init, (1,1,2) )

#write_xyz(a, repeat = 2  )

    #
#res_loop('hs443C.f','93',1, calc, conv, varset)    
   
   
# conv['e_imp'] = []
# res_loop('hs221C.f','93',1, calc, conv, varset,'e_imp',('hs221.f','83',1) )
# res_loop('hs332C.f','93',1, calc, conv, varset,'e_imp',('hs332.f','83',1) )
# res_loop('hs443C.f','93',1, calc, conv, varset,'e_imp',('hs443','83',20) )
# res_loop('t21bC.f','93',1, calc, conv, varset,'e_imp',('t21b_r','83',1) )
# res_loop('c1bC.f','93kp7',1, calc, conv, varset,'e_imp',('c1b_r','830',1) )
# res_loop('t111bC.f','93kp9',1, calc, conv, varset,'e_imp',('t111b_r','8302',1) )

# res_loop('hs221O.f','93',1, calc, conv, varset,'e_imp',('hs221.f','83',1) )
# res_loop('hs332O.f','93',1, calc, conv, varset,'e_imp',('hs332.f','83',1) )
# res_loop('hs443O.f','93',1, calc, conv, varset,'e_imp',('hs443','83',20) )
# res_loop('t21bO','93',20, calc, conv, varset,'e_imp',('t21b_r','83',1) )
# res_loop('c1bO','93kp7',20, calc, conv, varset,'e_imp',('c1b_r','830',1) )
# res_loop('t111bO.f','93kp9',1, calc, conv, varset,'e_imp',('t111b_r','8302',1) )
# print conv['e_imp']
# plot_conv(conv['e_imp'],calc,'e_imp')


#res_loop('hs221O.f','93',1, calc, conv, varset,'e_imp',('hs221.f','83',1) )
#res_loop('hs332O.f','93',1, calc, conv, varset,'e_imp',('hs332.f','83',1) )




#res_loop('c1g','929',range(1,10), calc, conv, varset, 'gbep',('c1b_r','830',1) )
#res_loop('t111g','9292',range(1,11), calc, conv, varset, 'gbep',('t111b_r','8302',1) )
#res_loop('t21g','93',range(1,6), calc, conv, varset,'gbep',('t21b_r','83',1))
#res_loop('t21g','935',range(11,15), calc, conv, varset,'gbep',('t21b_r','835',1))
#res_loop('csl71g','93',range(3,7), calc, conv, varset,'gbep',('csl71b_r','83kp8',1))

#res_loop('t111gCv','93kp9',range(1,6), calc, conv, varset, 'gbep',('t111bC.f','93kp9',1) )

#res_loop('hs443O','93',range(1,17), calc, conv, varset,'fit_ac' )
#res_loop('csl71b_r','83kp8',1, calc, conv, varset)

# res_loop('t111b_r','8301',1, calc, conv, varset)
# res_loop('t111b_r','8302',1, calc, conv, varset)

# res_loop('t111g','9291',12, calc, conv, varset)
#res_loop('t111g','9292',1, calc, conv, varset)
#
#res_loop('hs221.f','83',1, calc, conv, varset)
#res_loop('hs443C.f','93',1, calc, conv, varset)


"""6) Section for insterting reduced coordinates to gb cells with new lattice constants"""
# res_loop('hs554C','93kp9',1, calc, conv, varset)
# res_loop('t114g','93kp9',range(1,6), calc, conv, varset)

# insert('t114g', '93kp9', gb4_geo_folder+'T1/t114gC_template/', 't114gCi1', calc, "xred" )
# insert('t114g', '93kp9', gb4_geo_folder+'T1/t114gC_template/', 't114gCv', calc, "xred" )

# insert('hs554C', '93', gb4_geo_folder+'H/hs554C.f.93_template/grainA_s', 'hs554C.f', calc, "xred" )

# insert('t112g', '93kp9', gb4_geo_folder+'T1/t112gCO_template/target', 't112gCOi6.1-1is', calc, "xred" )
# insert('t112g', '93kp9', gb4_geo_folder+'T1/t112gCO_template/target', 't112gCvOvms', calc, "xred" )

# t111gC_r = copy.deepcopy(calc[('t111g','9292',2)])
# t111gC_r.end = replic( t111gC_r.end, (1,2,1))
# t111gC_r.name +=".replic_by_r2"
# for i in range(1,6):
#     calc[('t111g_replic_r2_from_v2','9292',i)] = t111gC_r

#write_xyz(  t111gC_r.end, repeat = 1  )
#insert('t111g_replic_r2_from_v2', '9292', gb4_geo_folder+'T1/t115gC_template/target', 't115gCi1', calc, "xred" )
#insert('t111g_replic_r2_from_v2', '9292', gb4_geo_folder+'T1/t115gC_template/target', 't115gCv', calc, "xred" )


#insert('hs332C', '93', 'geo/hs443C/from', 'hs443C')

#insert('t21g', '93', 'geo/t21gC/', 't21gCi', calc, "xred" )
#insert('c1g', '929', 'geo/c1gC_template/target', 'c1gCi1', calc, "xred" )
#insert('c1g', '929', 'geo/c1gC_template/target', 'c1gCv', calc, "xred" )

#insert('t111g', '9292', gb4_geo_folder+'T1/t111gO_template/target', 't111gOi1', calc, "xred" )
# insert('t111g', '9292', gb4_geo_folder+'T1/t111gO_template/target', 't111gOv', calc, "xred" )
#insert('c1g', '929', gb4_geo_folder+'C1/c1gO_template/target', 'c1gOi1', calc, "xred" )
#insert('c1g', '929', gb4_geo_folder+'C1/c1gO_template/target', 'c1gOv', calc, "xred" )

# insert('t21g', '93', gb4_geo_folder+'T2/t21gO_template/target', 't21gOi1', calc, "xred" )
# insert('t21g', '93', gb4_geo_folder+'T2/t21gO_template/target', 't21gOv', calc, "xred" )







"""Test total drift"""

#add_loop('bcc','83hsp',      1,calc,conv,varset,'up1') #ideal
#add_loop('bcc','83hsp_nosym',2,calc,conv,varset,'up1') #energy higher due to the shift
#add_loop('bcc','83hsp',      2,calc,conv,varset,'up1') #error



#for ise in varset:
#    try:
#        varset[ise].add_nbands
#    except AttributeError:
#        print "added for",ise
#        varset[ise].add_nbands = None
#    t = varset[ise].potdir[0]
#    varset[ise].potdir = {}
#    varset[ise].potdir[22] = t
#    varset[ise].ise = ise




#s = calc[('o1b','216',1)]
#s = calc[('c1b','230',1)]
#s = calc[('c1g','929',1)]
#s = calc[('t111g','9292',1)]
#s = calc[('t21g','93',1)]
#import cProfile

#cProfile.run('find_pores(s, 1.4, 0.56, 0.02)')
#st = find_pores(s, 1.4, 0.6, 0.05, 0.1)
#pores = find_pores(s, 1.4, 0.55, 0.1, 0.4,'central')
#write_xyz(st)


#add_loop('hs221C','93',range(1,17), calc, conv, varset, 'up1', "kpoint_conv")
#add_loop('hs221C.f','93',1, calc, conv, varset, 'up1', "kpoint_conv",'93') #from 93
#add_loop('hs221C.f','83is1',1, calc, conv, varset, 'up1', "",) #from 93



#add_loop('hs221C','93',range(1,17), calc, conv, varset, 'up1', "ecut_conv")
#add_loop('hs221C','93',range(1,17), calc, conv, varset, 'up1')#, "kpoint_conv")


#add_loop('hs221','83',range(1,17), calc, conv, varset, 'up1', "ecut_conv")




#add_loop('hs332','83',range(1,17), calc, conv, varset, 'up1')


#add_loop('t111bC.f','93kp9',1, calc, conv, varset, 'up1')

#add_loop('hs332C','93',range(1,17), calc, conv, varset, 'up1')
#add_loop('csl71b_r','83kp8',range(1,2), calc, conv, varset, 'up1')


#add_loop('hs443','83',20, calc, conv, varset, 'up1')

#add_loop('hs443C','93',range(1,17), calc, conv, varset, 'up1')

#last addings
#add_loop('c1bC','93kp7',range(1,17), calc, conv, varset, 'up1')

#add_loop('c1bC.f','93kp7',1, calc, conv, varset, 'up1')

#add_loop('c1gCv','93kp7',range(1,6), calc, conv, varset, 'up1')







#results_loop('t21b_r','835',1, calc)
#results_loop('o1b','216',1, calc)
#results_loop('t21g','935',range(11,15), calc,conv,varset,('t21b_r','835',1))

#add_loop('t21b_r','835',1, calc,conv, varset, 'up1')
#add_loop('t21g','93',range(1,6), calc, conv, varset,'up1')
#add_loop('t21g','935',range(11,15), calc, conv, varset,'up1')



#it = 'c1b_ml'; ise = '830'; verlist = range(1,17); setlist = (ise,) 
#it = 'c1b'; setlist = ("230",); verlist = range(1,2);
#it = 'c1g'; ise = '929'; verlist = range(1,11); setlist = (ise,)

#results_loop('c1b', '230', 1,  calc, conv)
#results_loop('c1b_r', '830', 1, calc, conv)
#results_loop('t111b', '2302', 1, calc, conv)
#results_loop('t111b_r', '8302', 1, calc, conv)
#results_loop('t111b_r', '8301', 1, calc, conv)

#results_loop('t111g', '9292', range(1,11), calc, conv, ('t111b_r','8302',1))
#results_loop('t111b_r', '9292', range(11,16), calc,conv, ('t111b_r','8302',1))
#results_loop('t111g', '9291', range(11,15), calc, conv, ('t111b_r','8301',1))

#results_loop('t111g',['9282',],[1,] ,  calc,conv,('t111b_r','8302',1))
#results_loop('t111g',['9292',],[1,] , conv, calc,('t111b_r','8302',1))
#results_loop('t111g',['9292',],[7,8] , conv, calc,('t111b_r','8302',1))
#results_loop('t111g',['9292',],range(21,25)+[26] ,  calc,conv,('t111b','2302',1))
#results_loop('t111g',['9292',],range(1,5) , conv, calc,('t111b_r','8302',1))
#results_loop('t111b_r',['8301',], [1,] , conv, calc,('t111b_r','8302',1))

#results_loop('c1g',['929',],range(1,11),calc,conv,('c1b_r','830',1))











#varset["900"].blockfolder = ""
#ise = "802"
#varset[ise].blockfolder = "npar_conv"
#varset[ise].blockfolder = "nband_conv"
#varset[ise].blockfolder = "kpoint_conv"
#varset[ise].blockfolder = "ecut_conv"
#varset[ise].blockfolder = "other"
#typconv = "kpoint_conv"
#varset[ise].conv = {}
#varset[ise].conv[typconv] = varset[ise].conv_kpoint

#it = 'o1b'; typconv = "ecut_conv"; verlist = range(1,2);
#it = 'o1b'; ise = '217'; typconv = "kpoint_conv"; verlist = range(1,2);
#it = 'o1b_ml'; ise = '801'; typconv = "ecut_conv"; verlist = range(26,51);
#it = 'o1b_ml'; ise = '801'; typconv = ""; verlist = range(1,26);  
#it = 'o1b'; ise = '218'; typconv = "kpoint_conv"; verlist = range(1,2);
#it = 'o1b'; ise = '203'; typconv = ""; verlist = range(1,2);





#analys section
#typconv = "kpoint_conv"
#typconv = "npar_conv"
#typconv = "ecut_conv"
#typconv = "nband_conv"
#typconv = "tsmear_conv"

#it = 'o1b'; ise = '216';
#setlist = ("203k2",); verlist = range(1,2);
#setlist = varset["215"].conv[typconv] #+varset["204"].conv[typconv]+varset["206"].conv[typconv]+["207",]
#setlist = varset["213"].conv[typconv]+varset["214"].conv[typconv]+varset["210"].conv[typconv]
#setlist = varset["216"].conv[typconv] #+varset["217"].conv[typconv]+varset["218"].conv[typconv]
#it = 'o1b'; setlist = varset["216"].conv[typconv]+varset["217"].conv[typconv]
#setlist = varset["802"].conv["ecut_conv"]+varset["801"].conv["kpoint_conv"]
#it = 'o1b_ml'; ise = '802'; verlist = range(26,51); setlist = (ise,) 
#it = 'o1b_ml'; ise = '801'; verlist = range(26,51); setlist = varset[ise].conv[typconv] 
#it = 'o1b_ml'; ise = '801'; verlist = range(1,26); setlist = (ise,)

#it = 'o1b_ml';
#verlist = range(26,51);  
   
#for v in verlist:
    #id1 = (it,"801",v)
    #id2 = (it,"801ec1",v)
    #print "Diff in pressure for version",v,calc[id1].extpress-calc[id2].extpress


#plot_kconv(conv[n],calc,"contour")


#plot_kconv( conv['o1b.204k'][0:1]+conv['o1b.204k'][4:7]+conv['o1b.204k'][1:3] )
#plot_kconv( conv['o1b.204k'][:5])
#plot_kconv( conv['o1b.203k'][7:] )
#plot_kconv( conv['o1b.203k'][4:7] )
#plot_kconv( conv['o1b.206k'] )


#for ise in '801','801kp1','801kp2','801ec1','802','802ec1':
    #print ("%s\_fit & %.4f & %.3f & %.3f"  % (fit_hex(0.0001,0.0001,150,230, 'o1b_ml', ise, range(26,51), calc) )  )
#print ("%s\_fit & %.4f & %.3f & %.3f"  % (fit_hex(0.0005,0.0005,150,230, 'o1b_ml', '800', range(1,26), calc) )  )
#print ("%s\_fit & %.4f & %.3f & %.3f"  % (fit_hex(0.0005,0.0005,150,230, 'o1b_ml', '801', range(1,26), calc) )  )
#print fit_hex(0.0001,0.0001,150,230, 'o1b_ml', '801', range(36,51), calc)
#print fit_hex(0.00015,0.00023,100,100, 'o1b_ml', '801', range(26,51), calc)
#print fit_hex(0.0005,0.0005,150,230, 'o1b_ml', '800', range(1,26), calc)
#print fit_hex(0.0005,0.0005,150,230, 'o1b_ml', '800', range(6,26), calc)


#try:
#    set_des[ise]
#except NameError:
#    print "Error, no description for this set. First write description"

#log.write( "There is no " +it+'.'+ise+'k'+" convergence list, I create new\n")

#s =CalculationVasp(1)
#s.read_geometry("geo/O1bulk.in.geo")
#log.write( str(s.rprimd[1][1]) )
#s.calc_kspacings()

#del calc[id]

#Section for tsmear convergences
#it = 'o1b'; ise = '203k10'; #conv[it+'.'+ise+'t'] = []
#for inputset in varset[ise].conv_tsmear:
#    id = (it,inputset,1)
    #add_calculation(it,inputset,1,"geo/O1bulk.in.geo","pure/kconv",calc,varset)
    #calc[id].read_results(1)
    #conv[it+'.'+ise+'t'].append(id)

#Section for kpoints convergences
#it = 'o1b'; ise = '203'; #conv[it+'.'+ise+'k'] = []
#for inputset in varset[ise].conv_kpoint:
#    id = (it,inputset,1)
    #add_calculation(it,inputset,1,"geo/O1bulk.in.geo","pure/kconv",calc,varset)
    #calc[id].calc_kspacings()
    #calc[id].read_results()
    #conv[it+'.'+ise+'k'].append(id)


    #Section for separate runs
#it = 'o1b';
#for inputset in '209',:
#    id = (it,inputset,1)
    #add_calculation(it,inputset,1,"geo/O1bulk.in.geo","pure/kconv",calc,varset,"up")
    #calc[id].read_results(1)













    # """2. Produce plot: co-segregation energies """
    # def plot_coseg(data):
    #     data.sort(key = itemgetter(1))

    #     rev_data =  zip(*data)



    #     if 1:
    #         to = -4 #divide into two blocks

    #         if 1:
    #             xlabel = 'Segregation site'; ylabel = 'Energy (meV)'
    #         else:
    #             xlabel = u""; ylabel = u"Энергия сегрегации (мэВ)"


    #         plot_bar (image_name = 'coseg',
    #             xlabel = xlabel, ylabel = ylabel, bottom = 0.3,
    #         #     CO1 = (rev_data[0][:to], rev_data[1][:to], "#F46529",   'dE2'           ),
    #         #     CO2 = (rev_data[0][:to], rev_data[2][:to], "#7FB005",   'Epair, gb'     ),
    #         #     CO3 = (rev_data[0][:to], rev_data[3][:to], "#7C144D",   'Ecoseg'        ),

    #             data1 = [
    #                 (rev_data[0][:to], rev_data[1][:to], "#F46529",   'dE2'           ),
    #                 (rev_data[0][:to], rev_data[2][:to], "#7FB005",   'Epair, gb'     ),
    #                 (rev_data[0][:to], rev_data[3][:to], "#7C144D",   'Ecoseg'        ),
    #                 ],
    #             data2 = [
    #                 (rev_data[0][to:], rev_data[1][to:], "#F46529",   'dE2'          ),
    #                 (rev_data[0][to:], rev_data[2][to:], "#7FB005",   'Epair, gb'    ),
    #                 (rev_data[0][to:], rev_data[3][to:], "#7C144D",   'Ecoseg'       ),
    #                 ],
    #             )






"""For find pores """
    #stepsi = np.linspace(0, 1, nstepsi)
    #stepsj = np.linspace(0, 1, nstepsj)
    #stepsk = np.linspace(0, 1, nstepsk)
    # stepsi = list(np.linspace(0, 1, nstepsi) )
    # stepsj = list(np.linspace(0, 1, nstepsj) )
    # stepsk = list(np.linspace(0, 1, nstepsk) )

    # stepsi = [ctypes.c_float(k) for k in stepsi]
    # stepsj = [ctypes.c_float(k) for k in stepsj]
    # stepsk = [ctypes.c_float(k) for k in stepsk]
    # #print stepsi[0]
    # print type(stepsi)
    # #print stepsj
    # scansi = 1/nstepsi # step of scaning of reduced coordinates
    # scansj = 1/nstepsj
    # scansk = 1/nstepsk
    # scansif = scansi / 1. #fine steps
    # scansjf = scansj / 1.
    # scanskf = scansk / 1.
    # v = np.zeros((3))
    # xredpores = []
    # npores = 0
    # natomlist = xrange(natom)
    #ii = []
    # = np.array
        #nstepsi = math.ceil(np.linalg.norm(rprimd[0])/step_dec)
    #nstepsj = math.ceil(np.linalg.norm(rprimd[1])/step_dec)
    #nstepsk = math.ceil(np.linalg.norm(rprimd[2])/step_dec)
    #print "Number of points", nstepsi*nstepsj*nstepsk
    #for i in np.nditer(stepsi):
        #for j in np.nditer(stepsj):
            #for k in np.nditer(stepsk):

    #a = ctypes.c_float(0.5)





    # for i in stepsi:
    #     for j in stepsj:
    #         for k in stepsk:
                #c_k = ctypes.c_float(k)

            
    """dsqrmin = 100
                #x = rprimd[0][0] * i + rprimd[1][0] * j + rprimd[2][0] * k
                #y = rprimd[0][1] * i + rprimd[1][1] * j + rprimd[2][1] * k
                #z = rprimd[0][2] * i + rprimd[1][2] * j + rprimd[2][2] * k
                #ii = np.subtract(i, xred1)
                #ii = [i - x for x in xred1]
                #print type(i)
                #print type(c_xred1[0])
                #print type(natom)
                
                #
                for n in natomlist:
                    #xx = x - xcart[n][0]
                    #ii = i - xred[n][0]
                    ii = i.value - xred1[n]
                    jj = j.value - xred2[n]
                    kk = k.value - xred3[n]
                    #
                    
                    if ii >  0.5: ii = ii - 1
                    elif ii < -0.5: ii = ii + 1
                    elif jj >  0.5: jj = jj - 1
                    elif jj < -0.5: jj = jj + 1
                    elif kk >  0.5: kk = kk - 1
                    elif kk < -0.5: kk = kk + 1 
                    xx = rprimd[0][0] * ii + rprimd[1][0] * jj + rprimd[2][0] * kk
                    yy = rprimd[0][1] * ii + rprimd[1][1] * jj + rprimd[2][1] * kk
                    zz = rprimd[0][2] * ii + rprimd[1][2] * jj + rprimd[2][2] * kk
                    #print xx,xx2
                    dsqr = xx * xx + yy * yy + zz * zz
                    if dsqr < dsqrmin: dsqrmin = dsqr #find distance between current point and nearest atom
                    #print ii,"py"
                #print xred1
                    #print dsqr, " ",
                print "dsqrmin python", dsqrmin#

            
                #return
                if False:#dsqrmin > refsqr: #this point situated in the spherical pore with radius larger than ri, considering spherical matrix atoms 
                    print "Pore was found"
                    npores+=1
                    v[0] = i; v[1] = j; v[2] = k;
                    #xredpores.append( v.copy() ) #add coordinates of pore
                    vc = v.copy() #coordinates of approximate center of pore
                    #Find all nearest points up to the surfaces of nearest spherical atoms.
                    #To prevent accounting of points we control rate of growth
                    #We assume that if the rate of growth starts to increase after decreasing, this means coalescence with neihgbour pore
                    iterat = 1 #
                    prate = 0 #previous rate
                    rate = 0
                    reincreasing = False
                    decreasing = False
                    #print "center of pore is",i,j,k
                    #print "scansif",scansif
                    #print "scansjf",scansjf
                    #print "scanskf",scanskf
                    while True:
                        localstepsi = np.linspace( -iterat * scansif, iterat * scansif, iterat*2)
                        localstepsj = np.linspace( -iterat * scansjf, iterat * scansjf, iterat*2)
                        localstepsk = np.linspace( -iterat * scanskf, iterat * scanskf, iterat*2)
                        #print localstepsi
                        #print localstepsj
                        #print localstepsk
                        n_localpoints = 8*iterat**3 #Total number of local points
                        n_porepoints = 0 #number of points inside the current pore
                        print "Number of local points", n_localpoints
                        xredporeslocal = []
                        for i in localstepsi:
                            for j in localstepsj:
                                for k in localstepsk:
                                    dsqrmin = cal_min_d(vc[0]+i, vc[1]+j, vc[2]+k, rprimd, xred, xcart, natom)
                                    if dsqrmin > rmsqr:
                                        #print "neihgbour pore added"
                                        v[0] = vc[0]+i; v[1] = vc[1]+j; v[2] = vc[2]+k;
                                        xredporeslocal.append( v.copy() )
                                        n_porepoints += 1
                        prate = rate
                        rate = n_porepoints * 1. / n_localpoints # part from local points inside the pore for this iteration
                        iterat += 1
                        print rate
                        if rate < prate: decreasing = True
                        if decreasing and rate > prate: break
                        if abs(rate - prate)<0.2 :break
                        xredporeslocal_prev = copy.copy(xredporeslocal)
                    #xredpores.append ( xredporeslocal ) #list of coordinate lists
                    print 'Number of points inside pore', n_porepoints

                    st.end.xcart.extend( xred2xcart(xredporeslocal_prev,rprimd) )"""
    #print "Number of found pores is",npores
    #print "Number of atoms is",natom

def cal_min_d(i,j,k,rprimd, xred,xcart, natom):
    dsqrmin = 100
    #x = rprimd[0][0] * i + rprimd[1][0] * j + rprimd[2][0] * k
    #y = rprimd[0][1] * i + rprimd[1][1] * j + rprimd[2][1] * k
    #z = rprimd[0][2] * i + rprimd[1][2] * j + rprimd[2][2] * k
    for n in range(natom):
        #xx = x - xcart[n][0]
        ii = i - xred[n][0]
        jj = j - xred[n][1]
        kk = k - xred[n][2]
        
        if ii >  0.5: ii = ii - 1
        if ii < -0.5: ii = ii + 1
        if jj >  0.5: jj = jj - 1
        if jj < -0.5: jj = jj + 1
        if kk >  0.5: kk = kk - 1
        if kk < -0.5: kk = kk + 1
        xx = rprimd[0][0] * ii + rprimd[1][0] * jj + rprimd[2][0] * kk
        yy = rprimd[0][1] * ii + rprimd[1][1] * jj + rprimd[2][1] * kk
        zz = rprimd[0][2] * ii + rprimd[1][2] * jj + rprimd[2][2] * kk
        #print xx,xx2
        dsqr = xx * xx + yy * yy + zz * zz
        if dsqr < dsqrmin: dsqrmin = dsqr #find distance between current point and nearest atom
    return dsqrmin
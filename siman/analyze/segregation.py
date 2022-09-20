import sys
from siman.analysis import fit_a
from siman.header import db
def inloop_segreg_analysis(outst, id, b_id, cl, analys_type, conv, n, base, readfiles, loadflag, choose_outcar):
    """
    analysis inside the res_loop; Used for grain boundary segregation project 

    """
    calc = db
    outst2 = outst
    if analys_type in ('e_seg', 'coseg'):
        try:
            b_id[1] 
            b_id = (b_id[0], b_id[1], id[2] + b_ver_shift)
        except:
            b_id = (b_id[0], id[1], id[2] + b_ver_shift)
        printlog('b_id', b_id)


    conv[n].append(id)
    # print base
    conv[base].append(b_id)


    if   b_id :
        n_m = cl.end.nznucl[0] # number of matrix atoms

        # if "4" not in calc[b_id].state:
        if readfiles:    
            # print(b_id)
            calc[b_id].read_results(loadflag, choose_outcar = choose_outcar)

        # print(b_id)
        if "4" in calc[b_id].state:    
            # print(id, b_id)
            # sys.exit()
            if calc[id].set.ngkpt != calc[b_id].set.ngkpt:
                printlog("Warning! you are trying to compare calcs with "+str(calc[id].set.ngkpt)+" and "+str(calc[b_id].set.ngkpt)+"\n")
                pass

            if calc[id].NKPTS != calc[b_id].NKPTS:
                printlog("Warning! you are trying to compare calcs with "+str(calc[id].NKPTS)+" and "+str(calc[b_id].NKPTS)+" kpoints \n")

            if 'gbe' in analys_type:      
                outst2 = gb_energy_volume(calc[id], calc[b_id])
        
            elif 'e_imp' in analys_type:
                calc[id].e_imp = e - e_b
                calc[id].v_imp = v - v_b

                #calc[id].v_imp = e - e_b
                outst2 += ("%.3f & %.2f  & %.2f &"% (e - e_b, v - v_b, (v - v_b)/v_b*100 ) )
                conv['e_imp'].append(id)
                a    = calc[id].hex_a;                     a_b  = calc[b_id].hex_a
                c    = calc[id].hex_c;                     c_b  = calc[b_id].hex_c
                ca = c/a;                                   ca_b = c_b/a_b
                outst_end = " & {0:.1f} & {1:.1f} & {2:.1f} & {3:.3f}".format((a - a_b)/a_b*100,  (c - c_b)/c_b*100, (ca - ca_b)/ca_b*100, (e - e_b - e1_r) )

            elif analys_type == 'e_2imp': #"""For calculation of energies of two impurities in big cell"""
                calc[id].e_imp = e - e_b            
                outst2 += ("%.0f "% ( (e - e_b)*1000 ) )  
                conv[it].append(id)


            elif analys_type in ('e_seg', 'coseg'): #"""For calculation of segregation and cosegregation energies"""
                e_b = calc[b_id].energy_sigma0 * bulk_mul
                n_m_b = calc[b_id].end.nznucl[0]
                v_b = calc[b_id].end.vol * bulk_mul
                # diffE = e - e_b/n_m_b*n_m
                diffE = e - e_b

                # outst2 += ("%.0f & %.2f "% ( (e - e_b)*1000, v - v_b ) )
                
                outst2 += " {:.3f} & {:.2f} ".format( (diffE - energy_ref), (v - v_b) ).center(6)
                outst2 +='&'
                # write_xyz(calc[id].end)
                # write_xyz(calc[b_id].end)
                result_list = [diffE - energy_ref, v - v_b]


            elif analys_type == 'matrix_diff': #
                printlog( 'Calculating matrix_diff...')
                

                diffE, diffV = matrix_diff(calc[id], calc[b_id], energy_ref)
                outst2 += " {:.3f} & {:.2f} &".format( diffE, diffV ).center(6)
                result_list = [diffE, diffV]


            elif analys_type == 'diff': #
                printlog( 'Calculating diff...')
                e_b = calc[b_id].energy_sigma0
                v_b = calc[b_id].end.vol
                diffE = e - e_b
                
                outst2 += " {:.3f} & {:.2f} &".format( (diffE - energy_ref), (v - v_b) ).center(6)
                result_list = [diffE - energy_ref, v - v_b]


    if analys_type == 'clusters':
        e1  = calc[imp1].energy_sigma0
        e2  = calc[imp2].energy_sigma0
        e_m = calc[matr].energy_sigma0
        n1 = calc[id].init.nznucl[1]
        if len(calc[id].init.nznucl) == 3:
            n2 = calc[id].init.nznucl[2]
        else:
            n2 = 0
        # print n1,n2
        outst2 += ("%.0f "% ( (e - n1*e1 - n2*e2 + (n1+n2-1)*e_m )*1000 / (n1+n2) ) )

    return outst2, conv



def outloop_segreg_analysis(b_id, analys_type, conv, n, description_for_archive, show, push2archive):


    if b_id: 
        bcl = calc[b_id]

    if analys_type == 'gbe':
        printlog("\nGrain boundary energy and excess volume fit:")
        plot_conv( conv[n], calc, "fit_gb_volume")

    elif analys_type == 'gbep':
        printlog("\nGrain boundary energy and excess volume fit:")
        # plot_conv( conv[n], calc, "fit_gb_volume")
        final_outstring = plot_conv( conv[n], calc, "fit_gb_volume_pressure")

    elif analys_type in ('e_seg', 'coseg') and len(verlist) > 3:

        #Test lateral sizes
        A   = calc[id].end.yzarea
        A_b = calc[b_id].end.yzarea
        
        if A != A_b: 
            printlog("Warning! you are trying to compare calcs with different lateral sizes: "+str(A)+" "+str(A_b))
            printlog( "Areas are ", A, A_b," A^3")
        
        #Show results 
        id1 = (it,inputset,verlist[0]) #choosen to save calculated values at first version of version set
        
        if readfiles and plot:           
            #print " \n\nImpurity at the interface :"
            e, v, emin, vmin       = plot_conv( conv[n], calc,  "fit_gb_volume2")
            #print " \n\nImpurity in the volume    :"
            e_b, v_b, e_bmin, v_bmin = plot_conv( conv[base], calc, "fit_gb_volume2")
            e_segmin = (emin - e_bmin) * 1000
            v_segmin =  vmin - v_bmin

            
            e_seg = (e - e_b * bulk_mul) * 1000
            v_seg =  v - v_b * bulk_mul


            # print(e_seg, e_segmin)

            calc[id1].e_seg = e_seg
            calc[id1].v_seg = v_seg
        
        if not hasattr(calc[id1], 'e_seg'): 
            printlog( "Warning! Calculation ", id1, 'does not have e_seg and v_seg. Try to run with readfiles = True to calculate it.')
            calc[id1].e_seg = 0; calc[id1].v_seg = 0
        


        natom = calc[id1].natom
        calc[id1].X = 1./natom
        v1 = v / natom
        calc[id1].Xgb = v1 / A # for grain boundary with 1 A width. For other boundaries should be divided by width. 
        #print ("__________________________________________________________________________")
        
        # print (" At zero pressure: segregation energy is %.0f  meV; Seg. volume is %.1f A^3; excess seg. vol. is %.2f A" %(e_seg, v_seg, v_seg/A ) )
        print ("At min: segregation energy is %.0f  meV; Seg. volume is %.1f A^3; excess seg. vol. is %.2f A" %(e_segmin, v_segmin, v_segmin/A ) )
        # print ("%s.fit.pe & %.0f & %.1f & %.2f & %.3f & %.1f" %(id[0]+'.'+id[1], e_seg, v_seg, v_seg/A, 1./A, 1./calc[id].natom * 100  ) )
        
        #Calculate distance from impurity to boundary and number of neighbours for version 2!
        id2 =(it,inputset, 2)
        st = calc[id2].end
        gbpos2 = calc[id2].gbpos 
        printlog( id2, 'is id2')
        # print gbpos2
        # print st.rprimd[0][0]/2.
        if gbpos2 == None:
            gbpos2 = 100
        gbpos1 = gbpos2 - st.rprimd[0][0]/2.
        d1 = abs(st.xcart[-2][0] - gbpos2)
        d2 = abs(st.xcart[-1][0] - gbpos2)
        dgb = d1; 
        iimp = -2
        if d2 < d1: 
            dgb = d2
            iimp = -1
        
        t = st.typat[iimp]
        z = st.znucl[t-1]
        segimp = element_name_inv(z) #Type of impurity closest to gb
        # print segimp, d

        id_m2   = (it+'.m',      '8'+inputset[1:], 2)
        
        if analys_type == 'e_seg':

            #calc e_seg2 and decomposition to mechanical and chemical contributions

            if id_m2 in calc: #additional analysis
                b_id2 = (b_id[0],inputset, 2)
                b_id_m2 = (b_id[0]+'.m', '8'+inputset[1:], 2)

                e_seg2 = (calc[id2].energy_sigma0 - calc[b_id2].energy_sigma0) * 1000
                e_m2   = (calc[id_m2].energy_sigma0 - calc[b_id_m2].energy_sigma0) * 1000
                e_ch2  = e_seg2 - e_m2
            else:
                e_seg2 = 0
                e_m2   =0
                e_ch2  =0


            #calculate number of close neibours around closest to gb imp
            x_central = st.xcart[iimp]
            st_r  = replic(st,   mul = (1,2,2), inv =  1 )
            st_rr = replic(st_r, mul = (1,2,2), inv = -1 ) # to be sure that impurity is surrounded by atoms

            dmax = 3
            nlist = [ x  for x, t  in zip(st_rr.xcart, st_rr.typat) if np.linalg.norm(x_central - x) < dmax and t == 1]
            nneigbours =  len(nlist)


            final_outstring = ("%s.fit.pe & %.0f & %.1f & %.2f & %.d & %4.0f & %4.0f & %4.0f & %s " %(
                id2[0]+'.'+id2[1], calc[id1].e_seg, calc[id1].v_seg, dgb, nneigbours, e_seg2, e_ch2, e_m2, segimp  ))

            # final_outstring = ("%s.fit.pe & %.0f & %.0f & %.1f & %.1f & %.2f & %.d & %4.0f & %4.0f & %4.0f & %s " %(
            #     id2[0]+'.'+id2[1], calc[id1].e_seg, e_segmin, calc[id1].v_seg, v_segmin , dgb, nneigbours, e_seg2, e_ch2, e_m2, segimp  )) #e_segmin and v_segmin are minimum energy (but at some pressure) and corresponing volume



            results_dic = [id2[0]+'.'+id2[1], calc[id1].e_seg, calc[id1].v_seg, dgb, nneigbours, e_seg2, e_ch2, e_m2, segimp]
        
        elif analys_type == 'coseg' :
            calc[id2].e_seg = calc[id1].e_seg #save in version 2
            calc[id2].v_seg = calc[id1].v_seg
            final_outstring = ("%s.fit.pe & %.0f & %.1f & %.1f & %.1f" %(id[0]+'.'+id[1], calc[id2].e_seg, calc[id2].v_seg, d1, d2 ))




        printlog(  final_outstring)
        printlog( '\\hline')


    elif analys_type == 'e_2imp':
        # plot_conv( conv[it], calc,  analys_type, conv_ext) #instead use plot_conv( conv['hs443OO'], calc,  'e_2imp', [conv['hs443CO'],conv['hs443CC']]) 
        pass



    elif analys_type == 'fit_ac':

        printlog ("name %s_template          acell  %.5f  %.5f  %.5f # fit parameters are &%.5f &%.5f &%i &%i"  % (fit_hex(0.00002,0.00003,4000,6000, it, inputset, verlist, calc) )  )    

    elif 'fit_a' in analys_type:
        # print(conv)
        fit_a(conv, n, description_for_archive, analys_type, show, push2archive)


    elif analys_type == 'dimer':
        """Fit of md steps obtained from constant speed; see vasp description for dimer"""
        # print calc[id].list_e_sigma0
        # calc[id].list_dE = []
        # for E in calc[id].list_e_without_entr:

        #     calc[id].list_dE.append(E - 2*calc[b_id].e_without_entr)

        # print calc[id].list_e_without_entr
        # print calc[id].list_dE
        calc[id].e_ref = calc[b_id].e_without_entr

        plot_conv( [id], calc,  "dimer")

    return
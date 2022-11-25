import sys, os
sys.path.append(os.path.dirname(__file__)+'../../savelyev')
sys.path.append(os.path.dirname(__file__)+'../../alglib_cpython')
# print sys.path
#from read_write_i import ReadWrite as RW
try:
    from SPL_i import ALGLIB as ALB
except:
    print('Cant import SPL_I, check paths to savelyev')
#from plot_3d_i import Plot_3D as Pl3d

def fit_hex(shag_a,shag_c,npoint_a,npoint_c,it,ise,verlist,calc,gb_volume = False, type_plot = 'contourf' ):




    acell_list = []
    etotal_list = []
    #i = 1     
    for v in verlist:
        #if i in (1,6,11,16,21): i+=1; continue #Cheking for reduced meshes
        #if i in (2,7,12,17,22): i+=1; continue
        #if i in (3,8,13,18,23): i+=1; continue     
        #i+=1
        #if v-25 in (11,12,13,14,15): continue
        #if v-25 in (3,8,13,18,23): continue
        id1 = (it,ise,v)
        acell = []
        if gb_volume:
            acell.append(calc[id1].v_gb); acell.append(calc[id1].v_gb); acell.append(0.1)
            acell_list.append( acell )
            etotal_list.append ( calc[id1].energy_sigma0 )
            acell = []
            acell.append(calc[id1].v_gb); acell.append(calc[id1].v_gb); acell.append(0.2)
        else:            
            acell.append(calc[id1].hex_a); acell.append(calc[id1].hex_a); acell.append(calc[id1].hex_c)
        acell_list.append( acell )
        etotal_list.append ( calc[id1].energy_sigma0 )

    # print acell_list
    # print etotal_list



    if not os.path.exists('a_c_convergence'):
        os.makedirs('a_c_convergence')

    f0 = open('a_c_convergence/fit_hex.out', 'w')
  
  
    f0.write('Calc equilibrium acell for\n')

    
    
    # 
    a=[]; c=[];
    for acell in acell_list:
        assert acell[0]-acell[1]==0, 'acell[0]-acell[1]!=0'
        if acell[0] not in a: a.append(acell[0])
        if acell[2] not in c: c.append(acell[2])
        
    f0.write('Read number points acell = '+ str(len(a))+'\n')
    f0.write('Read number points ccell = '+ str(len(c))+'\n\n')
    print ('Read number points acell = '+ str(len(a))+'\n')
    print ('Read number points ccell = '+ str(len(c))+'\n\n')

    f0.write('Limits acell read from file = ' + str(min(a))+'  '+ str(max(a)) + '; (max(a)-min(a))/2 = ' + str((max(a)-min(a))/2) +'\n')
    f0.write('Limits ccell read from file = ' + str(min(c))+'  '+ str(max(c)) + '; (max(c)-min(c))/2 = ' + str((max(c)-min(c))/2) +'\n\n')
          
    
    etot = [None for i in range(len(a)*len(c))]

    #      

    for i in range(len(etotal_list)):
        acell = acell_list[i]
        jj=0
        exit = False
        for j in c:
            if exit==True: break
            for k in a:
                if acell[0]==k and acell[2]==j:
                    etot[jj] = etotal_list[i]
                    exit = True
                    break
                else: jj+=1

    # print etot
    assert None not in etot, 'None in etot'
    
    # 
    e2 = ALB()
    e2.build_2d_bicubic_spline(a, len(a), c, len(c), etot, 1)
    
    #
    etot_min = min(etot)
    f0.write('Etot_min (without spline) = '+ str(etot_min)+' eV\n\n')
    ii=0; exit = False
    for i in c:
        if exit==True: break
        for j in a:
            if etot[ii]==etot_min:
                amin = j
                cmin = i
                exit = True
                break
            else: ii+=1
    f0.write('acell_min (without spline) = '+ str(amin)+' Angstrom\n')
    f0.write('ccell_min (without spline) = '+ str(cmin)+' Angstrom\n\n\n')
    f0.write('Found min etot for limitation:\n')
    f0.write('Acell = '+str(amin)+' +/- '+str(npoint_a*shag_a)+'  Angstrom'+' (shag = '+str(shag_a)+')'+'\n')
    f0.write('Ccell = '+str(cmin)+' +/- '+str(npoint_c*shag_c)+'  Angstrom'+' (shag = '+str(shag_c)+')'+'\n\n\n')    
    
    # 
    assert abs(e2.calc(amin,cmin,0)-etot_min)<10**(-10), 'e2.calc(amin,cmin,0)!=etot_min'
 
    e_min = 10**(8)
    for i in range(-npoint_c, npoint_c):
        c_cur = cmin+i*shag_c
        for j in range(-npoint_a, npoint_a):
            a_cur = amin+j*shag_a
            e_cur = e2.calc(a_cur,c_cur,0)
            if e_cur < e_min:
                e_min = e_cur
                a_min = a_cur
                c_min = c_cur
    f0.write('Etot_min (with spline) = '+ str(e_min)+' eV\n\n')                
    f0.write('acell_min (with spline) = '+ str(a_min)+' Angstrom\n')
    f0.write('ccell_min (with spline) = '+ str(c_min)+' Angstrom\n\n\n')
    
    
    # 
    al = amin-npoint_a*shag_a; ar = amin+npoint_a*shag_a
    cl = cmin-npoint_c*shag_c; cr = cmin+npoint_c*shag_c
    aspl_l = a_min - shag_a; aspl_r = a_min + shag_a
    cspl_l = c_min - shag_c; cspl_r = c_min + shag_c
    
    if (aspl_l-10**(-3))<=al or (aspl_r+10**(-3))>=ar:
        f0.write('Предупреждение!!!\n Найденный из сплайна минимум acell лежит на краю исследованной области!!!\n Увеличте исследуемую область по данному параметру!!!\n\n')
        print ('Предупреждение!!!\n Найденный из сплайна минимум acell лежит на краю исследованной области!!!\n Увеличте исследуемую область по данному параметру!!!\n\n')
    if (cspl_l-10**(-3))<=cl or (cspl_r+10**(-3))>=cr: 
        f0.write('Предупреждение!!!\n Найденный из сплайна минимум ccell лежит на краю исследованной области!!!\n Увеличте исследуемую область по данному параметру!!!\n\n')    
        print ('Предупреждение!!!\n Найденный из сплайна минимум ccell лежит на краю исследованной области!!!\n Увеличте исследуемую область по данному параметру!!!\n\n')    
    if a_min<min(a) or a_min>max(a): f0.write('Предупреждение!!!\n Найденный минимум энергии с использованием сплайна по параметру acell лежит вне рассчитанной области, считанной из файла!!!!\n Поэтому полученный из сплайн интерполяции результат не может считаться надежным.\n\n\n')
    if c_min<min(c) or c_min>max(c): f0.write('Предупреждение!!!\n Найденный минимум энергии с использованием сплайна по параметру ccell лежит вне рассчитанной области, считанной из файла!!!!\n Поэтому полученный из сплайн интерполяции результат не может считаться надежным.\n\n\n')
    
    f0.close()
    
    
    # 
        
    # 
    if len(acell_list)!=len(etotal_list): raise RuntimeError ('Не равны длины списков len(acell_list)!=len(etotal_list) для построения 3D рисунка!')
    acell_list1 = [tuple(i) for i in acell_list]
    setca_ac_dict = dict(zip(acell_list1, etotal_list))
    x_setca, y_setca, z_setca = [], [], []
    x1, y1, z1 = [], [], []
    ii = 0
    for key in sorted(setca_ac_dict):
        if ii==0: 
            key_cur = key
            ii=1
        if abs(key[0]-key_cur[0])<10**(-10): 
            x1.append(key[0])
            y1.append(key[2])
            z1.append(setca_ac_dict[key])
        else:
            key_cur = key
            x_setca.append(x1)
            y_setca.append(y1)
            z_setca.append(z1)
            x1 = [key[0]]; y1 = [key[2]]; z1 = [setca_ac_dict[key]];
    else:
        x_setca.append(x1)
        y_setca.append(y1)
        z_setca.append(z1)
    
    # 
    for i in range(len(z_setca)):
        for j in range(len(z_setca[i])): z_setca[i][j] -= e_min
       
    # 
    #e5 = Pl3d()
    #if show_plot_3d: show = True
    #else: show = False
    #e5.plot_surface('plot_ac_cell_energ', x_setca, y_setca, z_setca, show = show, write_axes=['a_cell, A','c_cell, A','e, eV'], title_graph = 'Файл данных '+file_data)
    import matplotlib.pyplot as plt
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter, MultipleLocator, AutoLocator
    
    fig = plt.figure()
    if type_plot == '3dim':
        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(x_setca, y_setca, z_setca, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)       

        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        fig.colorbar(surf, shrink=0.5, aspect=5)
    elif type_plot == 'contourf':
        cs = plt.contourf(x_setca, y_setca, z_setca, locator=LinearLocator(numticks = 50), cmap=cm.rainbow)
        cs.ax.set_xlabel('$a, \AA$')
        cs.ax.set_ylabel('$c, \AA$')
        cbar = plt.colorbar(format = FormatStrFormatter('%.0f'))
        cbar.set_label('$E, eV$', fontsize = 21)

    plt.show()           
        
    #return it+'.f.'+ise, e_min/calc[id1].natom, a_min, c_min, shag_a, shag_c, npoint_a, npoint_c
    return "{:s}_fit {:.4f} {:.4f} {:.4f}".format(it+'.f.'+ise, a_min, a_min, c_min, shag_a, shag_c, npoint_a, npoint_c)
    
    
    
    
    
    
    

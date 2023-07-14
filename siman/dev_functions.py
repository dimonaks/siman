#from geo.py
def stoichiometry_criteria(st1,st2):

    natom1 = st1.get_natom()
    natom2 = st2.get_natom()

    tra1 = st1.get_transition_elements()
    tra2 = st2.get_transition_elements()
    ntra1 = len(tra1)
    if ntra1 == 0: 
        ntra1 = natom1
    ntra2 = len(tra2)
    if ntra2 == 0: 
        ntra2 = natom2
    rat1 = natom1/ntra1
    rat2 = natom2/ntra2
    mul = ntra1/ntra2

    if rat1 == rat2:
        return 1
    else:
        return 0

def stoichiometry_criteria2(st1,st2, silent = 1):
    atoms1 = st1.get_elements()
    atoms2 = st2.get_elements()

    from collections import Counter
    el_dict1 = Counter(atoms1)
    el_dict2 = Counter(atoms2)
    el1 = list(el_dict1.keys())[0]
    el2 = list(el_dict1.keys())[1]
    # print(el_dict1)
    # print(el_dict2)
    ratio1 = el_dict1[el1]/el_dict1[el2]
    ratio2 = el_dict2[el1]/el_dict2[el2]

    if ratio1 == ratio2:
        if not silent:
            print('Stoichiometric')
        return 1
    else:
        if not silent:
            print('Non-stoichiometric')
            print(round(ratio1,2), round(ratio2,2))
        return 0

def symmetry_criteria(st):
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    st = st.convert2pymatgen()
    sym_criteria = SpacegroupAnalyzer(st).is_laue()
    if sym_criteria == True:
        print('Symmetric')
        return 1
    else:
        print('Non-symmetric')
        return 0

def symmetry_criteria_at(st):
    from collections import Counter
    
    els = Counter(st.get_elements())
    sym_criteria = 0
    for el in els:
        suf_at1 = st.get_surface_atoms(el, surface = 0, surface_width=1.5)
        suf_at2 = st.get_surface_atoms(el, surface = 1, surface_width=1.5)
        print(el, suf_at1, suf_at2)
        if len(suf_at1) == len(suf_at2):
            sym_criteria += 0
        else:
            sym_criteria += 1

    if sym_criteria == 0:
        print('Symmetric')
        return 1
    else:
        print('Non-symmetric')
        return 0



def sl_misfit(st1,st2, silent = 0):
    size1 = st1.rprimd_len()
    size2 = st2.rprimd_len()
    misfit = [(j-i)*100/j for i,j in zip(size1,size2)]
    # print('\n\nSize 1: {},{},{} A'.format(round(size1[0],2),round(size1[1],2),round(size1[2],2)))
    # print('Size 2: {},{},{} A'.format(round(size2[0],2),round(size2[1],2),round(size2[2],2)))
    if silent == 0:
        print('Misfit: {},{} % \n\n'.format(round(misfit[0],2),round(misfit[1],2)))
    return misfit

def fit2host(st_host, st_oxide):
    replic = [1,1,1]
    misf = sl_misfit(st_host,st_oxide, silent = 1)
    for m in (0,1):
        if 60 < abs(misf[m]) < 150:
            replic[m] +=1
        elif 150 < abs(misf[m]) < 250:
            replic[m] +=2
        elif 250 < abs(misf[m]) < 350:
            replic[m] +=3
        elif 350 < abs(misf[m]) < 450:
            replic[m] +=3
    st_oxide = st_oxide.replic(replic)

    return st_oxide


def hkl_slab(st, st_host, hkl, i_suf = None):

    os.remove('/home/anton/media/vasp/log_best1')
    f = open('/home/anton/media/vasp/log_best1', 'a')
    if slabs2:
        if not i_suf:
            for sl_i in range(0,len(slabs2)):
                # print(hkl)
                st2_new = st.update_from_pymatgen(slabs2[sl_i])
                misf = sl_misfit(st_host,st2_new, silent = 0)

                replic = [1,1,1]
                for m in (0,1):
                    if 80 < abs(misf[m]) < 110:
                        replic[m] +=1
                st2_new = st2_new.replic(replic)
                # print(replic)
                misf = sl_misfit(st_host,st2_new, silent = 1)
                print(hkl, sl_i)
                string = str(hkl) + ' ' + str(sl_i) + '  Misfit: {},{} % \n\n'.format(round(misf[0],2),round(misf[1],2))
                if abs(misf[0]) < 20 and abs(misf[1])<20:
                    f.write(string)
    # else:
    f.close()
    return misf, slabs2


def create_interface_solid(st_host, st_oxide, suf_host, i_suf_host = 0, 
    seek_mode = 0, seek_range = [0,2], check_shift = None, 
    hkl_lio = None, i_suf_lio = None, size = [5,5], ads_pos = None, z_shift = 1.5, lio_thick = 8):


    st1_init = st_host.copy()
    st2_init = st_oxide.copy()
    
    if suf_host:
        sc1 = st1_init.get_conventional_cell()#.replic([2,2,1])
        slabs1 = create_surface2(sc1, suf_host, shift = None, min_slab_size = 10, min_vacuum_size = 15, 
                 surface_i = 0, oxidation = None, suf = '', 
                 symmetrize = 1, cut_thickness = 0, return_one = 0, lll_reduce = 1, primitive = 1)
        slab1 = sc1.update_from_pymatgen(slabs1[i_suf_host])
        # slab1 = slabs1
        mul_matrix = ortho_vec(slab1.rprimd, [5,5,15], silent = 1) # matrix which allows to obtain supercell close to 10x10x10 A cube 
        st1 = create_supercell(slab1, mul_matrix, silent = 1)
    else:
        st1 = st_host


    
    if seek_mode:
        for h in range(seek_range[0],seek_range[1]):
            for k in range(seek_range[0],seek_range[1]):
                for l in range(seek_range[0],seek_range[1]):
                    hkl = [h,k,l]

                    if hkl != [0,0,0]:
                        slabs2 = create_surface2(st2_init, hkl, shift = None, min_slab_size = 10, min_vacuum_size = 25, surface_i = 0, oxidation = None, suf = '', 
                        primitive = 1, symmetrize = 0, cut_thickness = None, return_one = 0, write_poscar = 0, lll_reduce  = 1)
                        for sl in range(0,len(slabs2)):
                            print(hkl, sl)
                            st2 = st2_init.update_from_pymatgen(slabs2[sl])
                            slab_2 = fit2host(st1, st2)
                            sl_misfit(st1, slab_2, silent = 0)
    else:

        slabs2 = create_surface2(st2_init, hkl_lio, shift = 0, min_slab_size = lio_thick, min_vacuum_size = 15, surface_i = i_suf_lio, oxidation = None, suf = '', 
                primitive = 0, symmetrize = 0, cut_thickness = 0, return_one = 1, write_poscar = 0, lll_reduce  = 1)
        # st2 = st2_init.update_from_pymatgen(slabs2[i_suf_lio])
        st2 = slabs2
        # st2.jmol()
        slab_2 = fit2host(st1, st2)
        st2_init = slab_2

        from siman.impurity import make_interface
        z_max1 = 0
        z_max2 = 50
        st1 = st1.add_vacuum(2,40)
        for r in st1.xcart:
            if r[2] > z_max1:
                z_max1 = r[2]
        for r in slab_2.xcart:
            if r[2] < z_max2:
                z_max2 = r[2]
        

        if 0:
            elements1 = list(set(st1.get_elements()))
            elements2 = list(set(slab_2.get_elements()))
            # st1.jmol()
            suf_ats1 = []
            suf_ats2 = []
            for el in elements1:
                try:
                    suf_ats1.extend(st1.get_surface_atoms(el, surface = 1, surface_width=1.5))
                    print(el, suf_ats1)
                except TypeError:
                    None
            for el in elements2:
                try:
                    suf_ats2.extend(slab_2.get_surface_atoms(el, surface = 0, surface_width=1.5))
                    print(el, suf_ats2)
                except TypeError:
                    None


            for el1 in suf_ats1[0:1]:
                for el2 in suf_ats2[1:3]:
                    print('\n\nStarting configuration with {} and {} atoms\n'.format(el1, el2))
                    # xc1 = st1.xcart[23]
                    # xc2 = slab_2.xcart[8]
                    xc1 = st1.xcart[el1]
                    xc2 = slab_2.xcart[el2]
                    xc2[2] -= 1.8
                    # print(xc1, xc2)
                    mat0, pas_scaled2 = make_interface(st1, xc1, slab_2, xc2)
                    # mat0.replic([2,2,1]).jmol(r=2)
                    
                    interface_z_position = min([m[2] for m in mat0.xcart])/2 + max([m[2] for m in mat0.xcart])/2
                    print(interface_z_position)
                    av_pack = []
                    
                    for st in [st1, slab_2, mat0]:
                        av_dist = 0
                        n_atoms_i = 0
                        for i in range(0,st.natom):
                            if st == mat0: 
                                if (interface_z_position - 4) < st.xcart[i][2]< (interface_z_position+4):
                                    n_atoms_i += 1
                                    # print('ok')
                                    xxx = st.nn(i, n = 6, ndict = None, silent = 1, 
                                    from_one = 0, more_info = 0, oxi_state = 0, print_average = 0)
                                    d = 0
                                    for n in xxx['dist']:
                                        d+=n
                                    d = d/6
                                    av_dist += d
                                    # print(d, n_atoms_i)
                            else:
                                n_atoms_i = st.natom
                                xxx = st.nn(i, n = 6, ndict = None, silent = 1, 
                                from_one = 0, more_info = 0, oxi_state = 0, print_average = 0)
                                d = 0
                                for n in xxx['dist']:
                                    d+=n
                                d = d/6
                                av_dist += d

                        try:
                            av = round(av_dist/n_atoms_i,2)
                        except ZeroDivisionError:
                            av = 999
                        av_pack.append(av)
                        print('Average bond lengths is {} A'.format(av))

                    # mat0.replic([2,2,1]).jmol()
                    print(av_pack)
                    if abs(av_pack[2] - av_pack[0]) < 0.03:
                        print('\nGood interface!\n\n')


        print('\n\n\n', st1.rprimd, slab_2.rprimd, '\n\n\n',st1.get_angles(), slab_2.get_angles(),)
        if ads_pos:
            interface_list = []
            suf_ats2 = (slab_2.get_surface_atoms('OLiNa', surface = 0, surface_width=0.5))
            # print(suf_ats2)
            xc1 = ads_pos
            for at in suf_ats2[0:]:
                xc2 = slab_2.xcart[at]
                mat0, pas_scaled2 = make_interface(st1, xc1, slab_2, xc2)
                # mat0.jmol(rep=[3,3,1])
                interface_list.append(mat0)
                interface_list.append(pas_scaled2)
            return interface_list







        mat2, pas_scaled2 = make_interface(st1, [0, 4, z_max1+z_shift], slab_2, [0, 0.626, z_max2],)
       
        # st_wide = st1.add_vacuum(2,15)
        # mat3, pas_scaled3 = make_interface(st_wide, [0, 4, z_max1+10.5], slab_2, [0, 0.626, z_max2],)
        return [mat2, pas_scaled2]



from __future__ import division, unicode_literals, absolute_import 
import os, copy, shutil, sys
import numpy as np

try:
    # sys.path.append('/home/aksenov/Simulation_wrapper/ase') #path to ase library
    from ase.eos import EquationOfState
    ase_flag = True
except:
    print('ase is not avail; run   pip install ase')
    ase_flag = False


from siman import header
from siman.header import printlog, mpl, db
from siman.functions import element_name_inv, invert, get_from_server
from siman.picture_functions import plot_mep
from siman.geo import determine_symmetry_positions, local_surrounding, find_moving_atom, image_distance
from siman.database import push_figure_to_archive
from siman.small_functions import is_list_like, makedir
from siman.inout import write_xyz, read_xyz, write_occmatrix






def determine_barrier(positions = None, energies = None):

    """
    The sign of barrier determined by the curvuture at saddle point. Minimum at saddle point corresponds to negative barrier
    The saddle point is determined as maximum deviation from energy in initial position

    """

    import scipy

    if positions is None:
        positions = range(len(energies))

    if energies is None:
        printlog('Error! Please provide at least energies')


    spl = scipy.interpolate.PchipInterpolator(positions, energies)

    spl_der = spl.derivative()
    spl_der2 = spl_der.derivative()
    mi = min(positions)
    ma = max(positions)
    r = spl_der.roots()


    # print(r)
    r = r[ np.logical_and(mi<r, r<ma) ] # only roots inside the interval are interesting


    e_at_roots = spl(r)
    if len(e_at_roots) > 0:
        # diff_barrier = max( e_at_roots ) # the maximum value 
        printlog('roots are at ', r, e_at_roots)

        #find r for saddle point. the energy at saddle point is most far away from the energy at initial position by definition
        de_s = np.abs(e_at_roots-energies[0])
        i_r_de_max = np.argmax(de_s)
        # print(de_s)
        # print(i_r_de_max)


        r_de_max = r[i_r_de_max]
        e = spl(r_de_max)
        curvuture_at_saddle = spl_der2(r_de_max)
        
        sign =  - np.sign(curvuture_at_saddle)
        if curvuture_at_saddle < 0:
            critical_point_type = 'maximum'
        elif curvuture_at_saddle > 0:
            critical_point_type = 'minimum'
        else:
            critical_point_type = 'undefined'

        # print(type(r_de_max), type(e), critical_point_type)
        print('Saddle point at {:.2f} {:.2f} is a local {:}'.format(r_de_max, float(e), critical_point_type)  )


    else:
        print_and_log('Warning! no roots')
        # diff_barrier = 0
        sign = 1


    mine = min(energies)
    maxe = max(energies)
    de = abs(mine - maxe)
    # if de > diff_barrier:
    diff_barrier = de * sign

    print('Migration barrier is {:.2f}'.format( diff_barrier))
    # plt.plot(spl(np.linspace(0, ma, 1000)))
    # plt.show()

    return diff_barrier







def calc_redox(cl1, cl2, energy_ref = None, value = 0, temp = None, silent = 0):
    """
    Calculated average redox potential and change of volume
    cl1 (Calculation) - structure with higher concentration
    cl2 (Calculation) - structure with lower concentration
    energy_ref (float) - energy in eV per one alkali ion in anode; default value is for Li; -1.31 eV for Na, -1.02 eV for K
    
    temp(float) - potential at temperature, self.F is expected from phonopy calculations
    """
    if cl1 is None or cl2 is None:
        printlog('Warning! cl1 or cl2 is none; return')
        return
    if not hasattr(cl1.end, 'znucl') or not hasattr(cl2.end, 'znucl') :
        printlog('Warning! cl1 or cl2 is bad')
        return

    energy_ref_dict = {3:-1.9,  11:-1.31,  19:-1.02, 37:-0.93}
    z_alk_ions = [3, 11, 19, 37]


    #normalize numbers of atoms by some element except Li, Na, K
    alk1l = [] 
    alk2l = []
    # print cl1.end.znucl
    for i, z in enumerate(cl1.end.znucl):
        # print i, z
        if z in z_alk_ions: 
            alk1l.append(i)
            # print 'i_alk is found'
            continue
        # print i, z

        for j, zb in enumerate(cl2.end.znucl):
            if zb in z_alk_ions: 
                # j_alk = j
                alk2l.append(j)
                continue

            if z == zb:
                # print "I use ", z, " to normalize"
                i_n1 = i
                i_n2 = j

    n1  = cl1.end.nznucl[i_n1]
    n2  = cl2.end.nznucl[i_n2]

    # print(n1,n2)

    nz1_dict = {}
    nz2_dict = {}
    n_alk1 = 0
    n_alk2 = 0
    for z in z_alk_ions:
        nz1_dict[z] = 0 
        nz2_dict[z] = 0 

    for i in alk1l:
        nz1_dict[ cl1.end.znucl[i] ] = cl1.end.nznucl[i]
    for i in alk2l:
        nz2_dict[ cl2.end.znucl[i] ] = cl2.end.nznucl[i]

    for z in z_alk_ions:
        mul = (nz1_dict[z] / n1 - nz2_dict[z] / n2)
        # print(mul)
        if abs(mul) > 0: #only change of concentration of one ion type is allowed; the first found is used
            printlog('Change of concentration detected for ', element_name_inv(z))
            if not energy_ref: #take energy ref from dict
                energy_ref = energy_ref_dict[ z ]
            break

    # print(energy_ref)
    # print(cl1.energy_sigma0, cl2.energy_sigma0, mul)
    
    e1 = cl1.energy_sigma0 
    e2 = cl2.energy_sigma0
    if temp != None:
        #temperature corrections
        e1 += cl1.F(temp)
        e2 += cl2.F(temp)
        # print(cl1.F(temp), cl2.F(temp))
        # print(e1, cl1.energy_sigma0)
        # print(e2, cl2.energy_sigma0)


    if abs(mul) > 0:
        redox = -(  ( e1 / n1 -  e2 / n2 ) / mul  -  energy_ref  )
    else:
        redox = 0

    # print(n1, n2)

    dV = cl1.end.vol / n1 - cl2.end.vol / n2 

    vol_red = dV / (cl1.end.vol/n1) * 100 # %

    # final_outstring = ("{:} | {:.2f} eV \n1".format(cl1.id[0]+'.'+cl1.id[1], redox  ))
    final_outstring = ("{:45} | {:30} | {:10.2f} V | {:10.1f} % | {:6.2f}| {:6.2f}| {:6.0f}| {:6.0f} | {:3.0f}".format(cl1.name,cl2.name, redox, vol_red, cl1.energy_sigma0, cl2.energy_sigma0, cl1.maxforce, cl2.maxforce, value ))
    
    if not silent:
        printlog( final_outstring, end = '\n', imp = 'y' )
    
    try:
        cl1.set.update()

        results_dic = {'is':cl1.id[0], 'redox_pot':redox, 'id_is':cl1.id, 'id_ds':cl2.id, 
        'kspacing':cl1.set.kspacing, 'time':cl1.time/3600.,
        'mdstep':cl1.mdstep, 'ecut':cl1.set.ecut, 'niter':cl1.iterat/cl1.mdstep,
        'set_is':cl1.id[1], 'vol_red':vol_red }
    except:
        results_dic = {'redox_pot':redox, 'vol_red':vol_red}


    return results_dic



def matrix_diff(cl1, cl2, energy_ref = 0):
    #energy of substitutional impurity
    e = cl1.energy_sigma0
    v = cl1.end.vol
    n_m = cl1.end.nznucl[0]
    e_b = cl2.energy_sigma0
    n_m_b = cl2.end.nznucl[0]
    v_b = cl2.end.vol
    print(n_m_b, n_m)
    diffE = e - e_b/n_m_b*n_m - energy_ref
    
    return diffE, v - v_b



def form_en(sources, products, norm_el = None):
    """
    Calculate formation energy of reaction.

    sources, products - list of tuples (x, cl), where x is multiplier and cl is calculation
    norm_el  - which element to use for normalization

        'all' - normalize by total number of atoms
        'el'
        int - divide by this number
    """

    El = []
    Nzl = []


    for ls in [sources, products]:
        E = 0
        Nz = {}
        for x, cl in ls:
            E += x*cl.e0 
            for i, z in enumerate(cl.end.znucl):
                if z not in Nz:
                    Nz[z] = 0
                Nz[z] += x*cl.end.nznucl[i]
        El.append(E)
        Nzl.append(Nz)

    for z in Nzl[0]:
        if abs(Nzl[0][z] - Nzl[1][z]) > 1e-5:
            printlog('Error! Number of', invert(z), 'atoms in source and product are different!')

    # norm = 1
    if 'all' == norm_el:
        norm = sum(Nzl[0].values())
    elif type(norm_el) == str:
        norm = Nzl[0][invert(norm_el)]
    elif norm_el != None:
        norm = norm_el
    else:
        norm = 1
    # print('Normalizing by ', norm_el, norm, 'atoms')

    dE = (El[1]-El[0])/norm
    print('dE = {:4.2f} eV'.format(dE))


    return dE

def chgsum(cll, el, site, silent = 1):
    """
    calculate sum of Bader charges for particular atoms
    """


    for cl in cll:
        # print(cl.id, end = '  ')

        try:
            cl.chgsum[(el, site)] = 0
        except:
            pass
        if not hasattr(cl, 'charges') or len(cl.charges) == 0:
            cl.get_bader_ACF()
        # determine_symmetry_positions(cl.end, el, silent = 0)

    # print('')
    try:
        pos = determine_symmetry_positions(cll[0].end, el, silent = 1)
    except:
        printlog('chgsum() Warning!', cll[0].id, 'is broken!')
        return 0

    for p in pos[site]:
        ''
        for cl in cll:
            if not hasattr(cl, 'chgsum'):
                cl.chgsum = {}
                cl.chgsum[(el, site)] = 0

            cl.chgsum[(el, site)] += cl.charges[p]
        
            # print('{:5.3f}'.format(cl.charges[p]), end = '  ')
        # print('')
    if not silent:
        print('Sum of charges for ', el+str(site+1), ':')
    

    el_ind = cl.init.znucl.index(invert(el)) # index of element in znucl and zval and nznucl
    zval = cl.init.zval[el_ind] # number of electrons in chosen potential

    for cl in cll:
        cl.chgsum[(el, site)]/=len(pos[site])
        
        chgsum = zval - cl.chgsum[(el, site)]



        if cl == cll[0]:
            chgsum_ref = chgsum
        if not silent:

            print('{:5.2f}({:4.2f})'.format(chgsum, chgsum_ref-chgsum), )
    if not silent:
        print('\n')

    # print(cl.charges)
    return chgsum






def fit_a(conv, n, description_for_archive, analysis_type, show, push2archive):

    """Fit equation of state for bulk systems.

    The following equation is used::

       sjeos (default)
           A third order inverse polynomial fit 10.1103/PhysRevB.67.026103

                           2      3        -1/3
       E(V) = c + c t + c t  + c t ,  t = V
               0   1     2      3

       taylor
           A third order Taylor series expansion about the minimum volume

       murnaghan
           PRB 28, 5480 (1983)

       birch
           Intermetallic compounds: Principles and Practice,
           Vol I: Principles. pages 195-210

       birchmurnaghan
           PRB 70, 224107

       pouriertarantola
           PRB 70, 224107

       vinet
           PRB 70, 224107

       antonschmidt
           Intermetallics 11, 23-32 (2003)

       p3
           A third order polynomial fit

        Use::

           eos = EquationOfState(volumes, energies, eos='sjeos')
           v0, e0, B = eos.fit()
           eos.plot()

    """
    # e, v, emin, vmin       = plot_conv( conv[n], calc,  "fit_gb_volume2")


    from picture_functions import fit_and_plot

    alist = []
    vlist = []
    etotlist  = []
    magn1 = []
    magn2 = []
    alphas= []
    for id in conv[n]:
        cl = db[id]
        st = cl.end
        alist.append(cl.end.rprimd[0][0])
        etotlist.append(cl.energy_sigma0)
        vlist.append(cl.end.vol)
        magn1.append(cl.magn1)
        magn2.append(cl.magn2)
        alpha, beta, gamma = st.get_angles()
        alphas.append(alpha)
        print('alpha, energy: {:4.2f}, {:6.3f}'.format(alpha, cl.energy_sigma0))

    fit_and_plot(U1 = (alphas, etotlist, 'o-r'), 
        image_name = 'figs/angle', ylabel = 'Total energy, eV', xlabel = 'Angle, deg', xlim = (89, 92.6))

    if ase_flag:
        if 'angle' in analysis_type:
            eos = EquationOfState(alphas, etotlist, eos = 'sjeos')
        else:
            eos = EquationOfState(vlist, etotlist, eos = 'sjeos')
        # import inspect

        # print (inspect.getfile(EquationOfState))

        v0, e0, B = eos.fit()
        #print "c = ", clist[2]
        printlog( '''
        v0 = {0} A^3
        a0 = {1} A
        E0 = {2} eV
        B  = {3} eV/A^3'''.format(v0, v0**(1./3), e0, B), imp = 'Y'  )

        savedpath = 'figs/'+cl.name+'.png'
        makedir(savedpath)


        cl.B = B*160.218
        # plt.close()
        # plt.clf()
        # plt.close('all')
        if 'fit' in show:
            mpl.rcParams.update({'font.size': 14})

            eos.plot(savedpath, show = True)
            printlog('fit results are saved in ',savedpath, imp = 'y')
        else:
            printlog('To use fitting install ase: pip install ase')
    # plt.clf()

    if push2archive:
        push_figure_to_archive(local_figure_path = savedpath, caption = description_for_archive)

    return






def around_alkali(st, nn, alkali_ion_number):
    #return numbers and distances to 

    n_neighbours = nn
    alkali_ions = []
    dist = []
    ifmaglist = st.get_maglist()

    for i, typ, x in zip(range(st.natom), st.typat, st.xcart):
        z = st.znucl[typ-1]
        if z in header.ALKALI_ION_ELEMENTS:
            alkali_ions.append([i, z, x])

    if len(alkali_ions) > 0:
        if alkali_ion_number:
            kk = alkali_ion_number-1

            chosen_ion = (kk, st.znucl[st.typat[kk]-1], st.xcart[kk])
        else:
            chosen_ion = alkali_ions[0] #just the first one is used
                # alkali_ions[min(alkali_ions)]

        sur   = local_surrounding(chosen_ion[2], st, n_neighbours = n_neighbours, control = 'atoms', 
        periodic  = True, only_elements = header.TRANSITION_ELEMENTS)

        # print (sur)
        dist = np.array(sur[3]).round(2)
        numb = np.array(sur[2])

    else:
        numb = ifmaglist # if no alk ions show for all mag atoms
        chosen_ion = None

    return numb, dist, chosen_ion



def find_polaron(st, i_alk_ion, out_prec = 1):
    """
    #using magmom, find the transition atoms that have different magnetic moments
    #i_alk_ion - number of ion from 0 to calculate distances to transition metals
    out_prec (int) - precision of magmom output


    # maglist = cli.end.get_maglist()
    # magm = np.array(cli.end.magmom)

    # n_tm = len(magm[maglist])
    # # print(len(maglist))
    # numb, dist, chosen_ion = around_alkali(cli.end, n_tm, atom_num)
    # # print(magm[numb][1:]) 
    # mtm = magm[numb][1:] # the first is alkali

    # m_av = sum(mtm)/len(mtm)
    # print(mtm-m_av)
    """


    def zscore(s):
        # print(np.std(s))
        return (s - np.mean(s)) / np.std(s)



    magmom = np.array(st.magmom)
    if len(magmom) == 0 :
        printlog('Warning! magmom is empty')

    _, mag_numbers = st.get_maglist()

    pol = {}
    # for z in mag_numbers:
    #     pos = determine_symmetry_positions(st, invert(z))


    # sys.exit()
    magmom_tm = None
    for key in mag_numbers:
        printlog('Looking at polarons on transition atoms: ',invert(key) )
        numbs = np.array(mag_numbers[key])
        # print(numbs)
        # print(magmom)
        magmom_tm = magmom[numbs]
        dev = np.absolute(  zscore(magmom_tm) )
        # print(magmom_tm)
        # print(list(zip(magmom_tm, dev.round(1))))
        # p = np.where(dev>2)[0] # 2 standard deviations
        # print(dev>2)
        # print (type(numbs))
        nstd = 1.5
        # nstd = 4
        i_pols = numbs[dev>nstd]

        if len(i_pols) > 0:
            x1 = st.xcart[i_alk_ion]
            d_to_pols = []
            for j in i_pols:
                x2 = st.xcart[j]
                d, _ = st.image_distance(x1, x2, st.rprimd)
                d_to_pols.append(d)
            print('polarons are detected on atoms', [i+1 for i in i_pols], 'with magnetic moments:', magmom[i_pols], 'and distances: '+', '.join('{:2.2f}'.format(d) for d in d_to_pols), 'A'  )
            print('mag moments on trans. atoms:', magmom_tm.round(out_prec))
            
            pol[key] = i_pols
        else:
            print('no polarons is detected with nstd', nstd)
            print('mag moments on trans. atoms:', magmom_tm.round(out_prec))
            # print(' deviations                :', dev.round(1))
            pol[key] = None
    return pol, magmom_tm



def calc_oxidation_states(cl = None, st = None, silent = 1):

    if cl:
        st = cl.end
        ch = cl.charges
    if st:
        ch  = st.charges
    

    z_vals = []
    for j, z_val, el in zip(range(st.natom), st.get_elements_zval(), st.get_elements()):
        ox = z_val - ch[j]

        z_vals.append(ox)
        if not silent:
            ''
            print(el, '{:3.1f}'.format(ox))
    # print(list(zip(z_vals, self.end.get_elements())))
    # print(z_vals)
    return z_vals





def neb_analysis(cl, show, up = None, push2archive = None, old_behaviour = None, results_dic = None, fitplot_args = None, style_dic = None, params = None):
    """
    Analyse traectories and polarons

    params
        mep_shift_vector
    """
    if results_dic is None:
        results_dic = {}

    calc = header.calc
    path2mep_s = cl.project_path_cluster+'/'+cl.dir+'/mep.eps'
    itise = cl.id[0]+'.'+cl.id[1]
    # print(cl.ldauu)
    # sys.exit()
    name_without_ext = 'mep.'+itise+'.U'+str(max(cl.ldauu))
    path2mep_l = cl.dir+name_without_ext+'.eps'
    # print(path2mep_l)
    if not os.path.exists(path2mep_l) or '2' in up:
        ''
        get_from_server(files = path2mep_s, to_file = path2mep_l, addr = cl.cluster_address, )
        movie_to = cl.dir+'/movie.xyz'
        get_from_server(files = cl.project_path_cluster+'/'+cl.dir+'/movie.xyz', to_file = movie_to, addr = cl.cluster_address, )
        
        if os.path.exists(movie_to):
            makedir('figs/'+name_without_ext+'.xyz')
            shutil.copyfile(movie_to, 'figs/'+name_without_ext+'.xyz')



    # trying to get one image closest to the saddle point
    if old_behaviour and cl.version == 2: #old behaviour, now created automatically in add callc
        im = cl.set.vasp_params['IMAGES']
        # if im % 2 > 0: #odd
        #     i = im//2 + 1
        # else:
        #     i = im/2
        # if choose_image:
        #     i = choose_image

        for i in range(im):
            i+=1
            cl_i = copy.deepcopy(cl)
            cl_i.version+=i
            cl_i.id = (cl.id[0], cl.id[1], cl_i.version)
            cl_i.name = str(cl_i.id[0])+'.'+str(cl_i.id[1])+'.'+str(cl_i.id[2])
            # print cl_i.name
            cl_i.path["output"] = cl_i.dir+'0'+str(i)+"/OUTCAR"
            # for i in range():

            cl_i.associated_outcars = [ aso[2:] for aso in cl_i.associated_outcars  ]

            # print cl_i.path["output"] 
            cl_i.state = '2. Ready to read outcar'
            # if not os.path.exists(cl_i.path["output"]):
            #     load = 'o'
            outst2 = ("%s"%cl_i.name).ljust(name_field_length)
            if readfiles:
                print(outst2+'|'+cl_i.read_results(loadflag, show = show, choose_outcar = choose_outcar) )
            else:
                print_and_log(outst2+' | File was not read')
            

            if cl_i.id in calc: #move creation of calcs with images to add_neb
                ''
                # print_and_log('Please test code below this message to save prev calcs')
                # if cl_i != calc[cl_i.id]
                #     if hasattr(calc[cl_i.id], 'prev') and calc[cl_i.id].prev:
                #         prevlist = calc[cl_i.id].prev
                #     else:
                #         prevlist = [calc[cl_i.id]]
                #     cl_i.prev = prevlist
                #     calc[cl_i.id] = cl_i
            else:
                calc[cl_i.id] = cl_i






    # print path2mep_l
    if 0:
        if os.path.exists(path2mep_l):
            # get_from_server(file = path2mep_s, to = path2mep_l, addr = cluster_address)

            runBash('evince '+path2mep_l)
        else:
            a =  glob.glob(cl.dir+'*mep*')
            if a:
                runBash('evince '+a[0])


    


    cl1 = calc[cl.id[0], cl.id[1], 1]
    cl2 = calc[cl.id[0], cl.id[1], 2]
    

    atom_num = find_moving_atom(cl1.end, cl2.end)
    # cl1.poscar()
    # cl2.poscar()

    # print('atom_num',atom_num)
    # sys.exit()

    #prepare lists
    ni = cl.set.vasp_params['IMAGES']
    vlist = [1]+list(range(3, ni+3) )+[2]
    # print( vlist)
    mep_energies = []
    atom_pos     = []



    pols = []
    sts = []
    sts_loc = []
    dAO = [] # A-(O,F) distance for each image

    for v in vlist:
        cli = calc[cl.id[0], cl.id[1], v]
        # print(cl.id[0], cl.id[1], v, cli.state)
        if '4' not in cli.state and 'un' not in up:
            printlog('Attention! res_loop(): analys_type == neb, Calc',cli.id,'is not finished; return')
            return {}, []
        # print cli.id
        # cli.end = return_to_cell(cli.end)
        # mep_energies.append(  min(cli.list_e_sigma0)   ) #use minimum energy - not very good, sometimes unconverged energy could be lower! 
        mep_energies.append(  cli.energy_sigma0   ) #use last energy 
        atom_pos.append( cli.end.xcart[atom_num] )

        # Find polaron positions
        if 1 or 'polaron' in show:
            pol, mag = find_polaron(cli.end, atom_num)
            if pol:
                for key in pol:
                    if np.any(pol[key]):
                        for n in pol[key]:
                            if n not in pols:
                                pols.append(n)
            else:
                ''
                # print('Mag_moments on trans,', mag.round(1))
        
        if 1 or 'neb_geo' in show:
            #visualization of path
            # print(atom_num)
            st = copy.deepcopy(cli.end)
            # print('moving_atom', st.xcart[atom_num])
            info = st.nn(atom_num, 15, from_one = False, silent = 1)
            

            st.moving_atom_i = atom_num
            st_loc = info['st']


            # print(st_loc.xcart)
            # st_loc = st_loc.shift
            
            if v == vlist[0]:

                st1 = copy.deepcopy(st)


                vec = st.center_on(atom_num)
                # vec = np.asarray([0.,0.,0.])

            
                if params is not None and 'mep_shift_vector' in params:
                    # vec += np.array([0.11,0.11,0]) # path4
                    # print(params['mep_shift_vector'])
                    vec += np.array(params['mep_shift_vector']) # path4

            # print(vec)
            st_loc = st_loc.shift_atoms(vec)
            if 0:
                st_loc.write_xyz()
            # st.write_cif('xyz/'+st.name)
            st.shift_atoms(vec).write_xyz()
            sts_loc.append(st_loc)

            st1 = st1.add_atom(st.xred[atom_num], 'Rb')

            sts.append(st.shift_atoms(vec))


            if 0:
                info1 = st.nn(atom_num, 2, from_one = False, silent = 1)
                print('Average_distance A-2(O,F)', info1['av(A-O,F)'], 'A')
                dAO.append (info1['av(A-O,F)'])


            if 0 or 'neb_geo' in show:
                av = st.nn(atom_num, 2, from_one = False, silent = 1, more_info = 1)['avsq(A-O,F)']
                print('Average squared distance A-2(O,F)', av, 'A')

                info2 = st.nn(atom_num, 4, from_one = False, silent = 1)
                print('Average_distance A-4(O,F)', info2['av(A-O,F)'], 'A')
                print('Elements are ', info2['el'])

                info3 = st.nn(atom_num, 6, from_one = False, silent = 1, more_info = 1)
                print('Average_distance A-6(O,F)', info3['av(A-O,F)'], 'A')
                print('Average_deviation A-6(O,F)', info3['avdev(A-O,F)'], 'mA')
                print('Elements are ', info3['el'])



    write_xyz(sts = sts) # write traectory
    write_xyz(sts = sts_loc) # write traectory

    st1 = st1.shift_atoms(vec)
    st1.name +='_all'
    # st1.write_cif('xyz/'+st1.name)
    st1.write_xyz()


    if dAO: # find maximum change of distance during migration
        dAO_change = abs(min(dAO) - max(dAO))
        results_dic['dAO_change'] = dAO_change

    results_dic['sts_loc'] = sts_loc # list of local structures, each structure contains dlist - distances from central cation to anions, and ellist - types of elements
    results_dic['sts'] = sts # list of mep structures, each structure contains moving_atom_i - number of moving atom


    if len(pols) > 0:
        print('During migration of alkali ions polarons are detected on atoms:', pols)
    elif len(pols) > 1:
        printlog('Attention! polaron is moving during migration! Obtained barrier is ambiguous')
    else:
        printlog('Compare magnetic moments above! In principle should be the same!')


    # print np.array(atom_pos)

    #test if the distances between points are not spoiled by PBC 
    nbc = range(-1, 2)
    jj=0
    for x in atom_pos:

        x2 = atom_pos[jj+1]
        r = cl.end.rprimd
        d1, _ = image_distance(x, x2, r, order = 1) #minimal distance
        x2_gen = (x2 + (r[0] * i  +  r[1] * j  +  r[2] * k) for i in nbc for j in nbc for k in nbc) #generator over PBC images
        x2c = copy.deepcopy(x2)
        ii = 0
        while  np.linalg.norm(x - x2c) > d1: #find the closest PBC image position
            if ii > 100:
                break
            ii+=1
            x2c = next(x2_gen)
        atom_pos[jj+1] = x2c
        jj+=1
        if jj == len(atom_pos)-1: # the last point is not needed, we could not use slice since we need to use changed atom_pos in place
            break
        # print np.linalg.norm(x - x2c), d1



    _, diff_barrier = plot_mep(atom_pos, mep_energies, plot = 0, show = 0, fitplot_args = fitplot_args, style_dic = style_dic)

    results_dic['barrier'] = diff_barrier
    
    middle_image = len(vlist) // 2
    results_dic['dEm1']    = mep_energies[middle_image] - mep_energies[0]
    


    cl1.barrier = diff_barrier
    cl2.barrier = diff_barrier


    results_dic['atom_pos'] = atom_pos
    results_dic['mep_energies'] = mep_energies


    if 'mep' in show:
        if 'mepp' in show:
            show_flag = True
        else:
            show_flag = False
        # sys.exit()

        plot_mep(atom_pos, mep_energies, image_name = 'figs/'+name_without_ext+'_my.eps', show = show_flag, fitplot_args = fitplot_args,  style_dic = style_dic)






    if push2archive:
        path2saved, _ = plot_mep(atom_pos, mep_energies, image_name = 'figs/'+name_without_ext+'_my', fitplot_args = fitplot_args, style_dic = style_dic)
        push_figure_to_archive(local_figure_path = path2saved, caption = description_for_archive)








        if 0: #copy files according to chosen outcar to run nebresults locally 
            wd = cl_i.dir
            out_i = cl_i.associated_outcars[choose_outcar-1]
            out_1 = calc[cl.id[0],cl.id[1], 1].associated_outcars[choose_outcar-1]
            out_2 = calc[cl.id[0],cl.id[1], 2].associated_outcars[choose_outcar-1]
            # print out_1
            # print out_2 
            shutil.copyfile(wd+out_1, wd+'00/OUTCAR')
            shutil.copyfile(wd+out_2, wd+'04/OUTCAR')
            for d in ['01/','02/','03/' ]:
                shutil.copyfile(wd+d+out_i, wd+d+'OUTCAR')

                # print wd+d+out_i

    return results_dic



def set_oxidation_states(st):
    pm = st.convert2pymatgen()
    pm.add_oxidation_state_by_guess()
    st = st.update_from_pymatgen(pm)
    # print(pm)
    return st



def suf_en(cl1, cl2, silent = 0):
    """Calculate surface energy
    cl1 - supercell with surface
    cl2 - comensurate bulk supercell
    the area is determined from r[0] and r[1];- i.e they lie in surface

    """
    st1 = cl1.end
    st2 = cl2.end
    # pm = st1.convert2pymatgen(oxidation = {'Y':'Y3+', 'Ba':'Ba2+', 'Co':'Co2.25+', 'O':'O2-'})

    A = np.linalg.norm( np.cross(st1.rprimd[0] , st1.rprimd[1]) )
    # print(A)

    if st1.natom%st2.natom > 0:
        printlog('Warning! non-stoichiometric slab, check system sizes: natom1 = ', st1.natom, 'natom2 = ', st2.natom, st1.natom/st2.natom)

    mul = st1.natom/st2.natom
    gamma = (cl1.e0 - cl2.e0*mul)/2/A* header.eV_A_to_J_m
    if not silent:
        print('Surface energy = {:3.2f} J/m2   | {:} | {:} '.format(gamma, cl1.id, cl2.id))
    
    return gamma
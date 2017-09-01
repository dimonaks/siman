from __future__ import division, unicode_literals, absolute_import 
import numpy as np

import header
from header import printlog, mpl
from functions import element_name_inv, invert
from geo import determine_symmetry_positions, local_surrounding
from database import push_figure_to_archive
from picture_functions import fit_and_plot
from header import db
from small_functions import is_list_like, makedir

try:
    # sys.path.append('/home/aksenov/Simulation_wrapper/ase') #path to ase library
    from ase.utils.eos import EquationOfState
    ase_flag = True
except:
    print('ase is not avail; run   pip install ase')
    ase_flag = False


def calc_redox(cl1, cl2, energy_ref = None, value = 0):
    """
    Calculated average redox potential and change of volume
    cl1 (Calculation) - structure with higher concentration
    cl2 (Calculation) - structure with lower concentration
    energy_ref (float) - energy in eV per one alkali ion in anode; default value is for Li; -1.31 eV for Na, -1.02 eV for K
    """
    if cl1 is None or cl2 is None:
        printlog('cl1 or cl2 is none; return')
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
        if abs(mul) > 0: #only change of concentration of one ion type is allowed; the first found is used
            printlog('Change of concentration detected for ', element_name_inv(z))
            if not energy_ref: #take energy ref from dict
                energy_ref = energy_ref_dict[ z ]
            break

    # print(energy_ref)
    # print(cl1.energy_sigma0, cl2.energy_sigma0, mul)
    if abs(mul) > 0:
        redox = -(  ( cl1.energy_sigma0 / n1 - cl2.energy_sigma0 / n2 ) / mul  -  energy_ref  )
    else:
        redox = 0

    # print(n1, n2)

    dV = cl1.end.vol / n1 - cl2.end.vol / n2 

    vol_red = dV / (cl1.end.vol/n1) * 100 # %

    # final_outstring = ("{:} | {:.2f} eV \n1".format(cl1.id[0]+'.'+cl1.id[1], redox  ))
    final_outstring = ("{:45} | {:30} | {:10.2f} V | {:10.1f} % | {:6.2f}| {:6.2f}| {:6.0f}| {:6.0f} | {:3.0f}".format(cl1.name,cl2.name, redox, vol_red, cl1.energy_sigma0, cl2.energy_sigma0, cl1.maxforce, cl2.maxforce, value ))
    
    printlog( final_outstring, end = '\n', imp = 'y' )
    try:
        cl1.set.update()

        results_dic = {'is':cl1.id[0], 'redox_pot':redox, 'id_is':cl1.id, 'id_ds':cl2.id, 
        'kspacing':cl1.set.kspacing, 'time':cl1.time/3600.,
        'mdstep':cl1.mdstep, 'ecut':cl1.set.ecut, 'niter':cl1.iterat/cl1.mdstep,
        'set_is':cl1.id[1], 'vol_red':vol_red }
    except:
        results_dic = {}


    return results_dic



def matrix_diff(cl1, cl2, energy_ref = 0):
    
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
    Calculate formation energy of reaction
    sources, products - list of tuples (x, cl), where x is multiplier and cl is calculation
    norm_el  - which element to use for normalization
        'all' - normalize by total number of atoms

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


    print('dE = {:4.2f} eV'.format((El[1]-El[0])/norm))



def chgsum(cll, el, site):
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

    print('')
    try:
        pos = determine_symmetry_positions(cll[0].end, el, silent = 0)
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
    print('Sum of charges for ', el+str(site+1), ':')
    

    el_ind = cl.init.znucl.index(invert(el)) # index of element in znucl and zval and nznucl
    zval = cl.init.zval[el_ind] # number of electrons in chosen potential

    for cl in cll:
        cl.chgsum[(el, site)]/=len(pos[site])
        
        chgsum = zval - cl.chgsum[(el, site)]



        if cl == cll[0]:
            chgsum_ref = chgsum

        print('{:5.2f}({:4.2f})'.format(chgsum, chgsum_ref-chgsum), end = '  ')
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



def find_polaron(st, i_alk_ion):
    #using magmom, find the transition atoms that have different magnetic moments
    #i_alk_ion - number of ion from 0 to calculate distances to transition metals



    # maglist = cli.end.get_maglist()
    # magm = np.array(cli.end.magmom)

    # n_tm = len(magm[maglist])
    # # print(len(maglist))
    # numb, dist, chosen_ion = around_alkali(cli.end, n_tm, atom_num)
    # # print(magm[numb][1:]) 
    # mtm = magm[numb][1:] # the first is alkali

    # m_av = sum(mtm)/len(mtm)
    # print(mtm-m_av)

    def zscore(s):
        # print(np.std(s))
        return (s - np.mean(s)) / np.std(s)

    magmom = np.array(st.magmom)
    _, mag_numbers = st.get_maglist()

    pol = {}
    # for z in mag_numbers:
    #     pos = determine_symmetry_positions(st, invert(z))


    # sys.exit()

    for key in mag_numbers:
        printlog('Looking at polarons on transition atoms: ',invert(key) )
        numbs = np.array(mag_numbers[key])
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
                d, _ = st.image_distance(x1, x2)
                d_to_pols.append(d)
            print('polarons are detected on atoms', [i+1 for i in i_pols], 'with magnetic moments:', magmom[i_pols], 'and distances: '+', '.join('{:2.2f}'.format(d) for d in d_to_pols), 'A'  )
            print('mag moments on trans. atoms:', magmom_tm.round(1))
            
            pol[key] = i_pols
        else:
            print('no polarons is detected with nstd', nstd)
            print('mag moments on trans. atoms:', magmom_tm.round(1))
            # print(' deviations                :', dev.round(1))
            pol[key] = None
    return pol, magmom_tm

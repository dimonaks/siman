# -*- coding: utf-8 -*-
import sys, os


import numpy as np
import numpy


try:
    from scipy.integrate import simps, trapz
    from scipy import signal
    from scipy.interpolate import interp1d
    from scipy import interpolate
    from scipy.stats.stats import pearsonr 
    from scipy.spatial.distance import sqeuclidean
except:
    print('Warning scipy is not available, some dos functions will not work')

try:
    from ase.calculators.vasp import VaspDos
except:
    print('Warning ase is not installed. I need ase to parse to DOSCAR. install ase with    pip install ase')


from siman import header
from siman.header import printlog
from siman.picture_functions import fit_and_plot
from siman.functions import element_name_inv, smoother
from siman.geo import local_surrounding, determine_symmetry_positions
from siman.small_functions import latex_chem




def det_gravity2(energy, dos, Erange = (-100, 0)):
    """Determine center of gravity for DOS and return values of energy for dos
    INPUT:
    energy - list of energies
    dos    - list of corresponding dos
    Erange - window of energy to determine center of gravity
    """

    sum_dos   = 0
    sum_dos_E = 0

    for i, E in enumerate(energy):
        
        if E < Erange[0]: 
            continue

        if E > Erange[1]: 
            break
        
        sum_dos   += dos[i] 
        sum_dos_E += dos[i]*E 

    gc = sum_dos_E/sum_dos

    return gc

def det_gravity(dos, Erange = (-100, 0), key = None):
    """Determine center of gravity for DOS and return values of energy for d6 orbitals in list
    INPUT:
    dos - ase dos type with added d6 - sum of d orbitals over neighbors to impurity atoms
    Erange - window of energy to determine center of gravity
    """

    if key is None:
        key = 'd6'

    sum_dos_E = {}
    sum_dos   = {}
    sum_dos[key] = 0
    sum_dos_E[key] = 0

    for i, E in enumerate(dos.energy):
        
        if E < Erange[0]: 
            continue

        if E > Erange[1]: 
            break
        if key == 'd6':
            sum_dos['d6']   += dos.d6[i] 
            sum_dos_E['d6'] += dos.d6[i]*E 

        else:
            pass

    d6_gc = sum_dos_E[key]/sum_dos[key]





    if 0: #old
        nn =13
        sum_dos_E = [0 for i in range(nn)]
        sum_dos = [0 for i in range(nn)]

        i=-1
        for E in st.dos[0]:

            i+=1
            if E < -6: continue
            if E >0: break
            l = 1 #s
            sum_dos_E[l] += st.dos[l][i] * E; sum_dos[l] += st.dos[l][i]
            l = 2 #p
            sum_dos_E[l] += st.dos[l][i] * E; sum_dos[l] += st.dos[l][i]
            l = 3 #d
            sum_dos_E[l] += st.dos[l][i] * E; sum_dos[l] += st.dos[l][i]
            l = 0 #p+d
            sum_dos_E[l] += (st.dos[2][i] + st.dos[3][i]) * E; sum_dos[l] += st.dos[2][i] + st.dos[3][i]
            l = 11 #full
            sum_dos_E[l] += st.dos[l][i] * E; sum_dos[l] += st.dos[l][i]
            l = 12 #s+p+d
            sum_dos_E[l] += (st.dos[1][i] + st.dos[2][i] + st.dos[3][i]) * E; sum_dos[l] += st.dos[1][i] + st.dos[2][i] + st.dos[3][i]


        #Determine center of DOS
        st.Ec = [0 for i in range(nn)]
        #    if debug: print "Gravity Centers for s,p,d are (eV):"

        for l in 1,2,3,0,11,12:
            if sum_dos[l] <=0 : 
                printlog('err 123'); 
                sum_dos[l] = 0.00001
            st.Ec[l] = sum_dos_E[l] / sum_dos[l]
            # if debug: print st.Ec[l]


    return d6_gc







def plot_dos(cl1, cl2 = None, dostype = None, iatom = None, iatom2= None,
    orbitals = ('s'), up = None, neighbors = 6, show = 1, labels = None,
    path = 'dos', xlim = (None, None), ylim = (None,None), savefile = True, plot_param = {}, suf2 = '', nsmooth = 3,
    lts2 = '--', split_type = 'octa', plot_spin_pol = 1, show_gravity = None, 
    efermi_origin = True, efermi_shift = 0, invert_spins  = 0, name_suffix = '', image_name = None, color_dict = None):
    """
    cl1 (CalculationVasp) - object created by add_loop()
    dostype (str) - control which dos to plot:
        'total'   - plot total dos
        'diff_total' - difference of total dos, use cl2 for second calculation
        'partial' - partial dos

    orbitals (list of str) - 
        any from 's, p, d, py, pz, px, dxy, dyz, dz2, dxz, dx2' where 'p' and 'd' are sums of projections
        also to sum around neigbours use p6 and d6 and neighbors parameter
        p_all - sum over all atoms, p states 
        d_all - sum over all atoms, d states

    up - 'up2' allows to download the file once again
    labels - two manual labels for cl1 and cl2 instead of auto

    iatom (int) - number of atom starting from 1 to plot DOS;
    iatom ([float]*3) - cartesian coordinates of point around which atoms will be found
    show (bool) - whether to show the dos 
    path (str)  - path to folder with images

    neighbors - number of neighbours around iatom to plot dos on them using p6 or d6; only p6 is implemented to the moment in plot section

    xlim, ylim (tuple)- limits for plot

    color_dict (dict) - custom dict of colors for orbitals. eg: {'s':'g', 'p':}

    plot_param - dict of parameters to fit_and_plot
        dashes - control of dahsed lines

    suf2 - additional suffix for label
    name_suffix - modify name

    image_name - user image name

    # nsmooth = 15 # smooth of dos
    lts2 - style of lines for cl2

    split_type - 
        octa  - the names are t2g and eg
        tetra - the names are t2 and e


    plot_spin_pol -
        0 - spin-polarized components are summed up

    show_gravity (list) - print gravity centers (i, type, range, ); i - 1 or 2 cl
        type (str)
            'p6' - for p orbitals of neighbors
            'p'


    efermi_origin 
        True - e-fermi is zero energy
        False - e-fermi is left, its value is shown
    efermi_shift (float) - additional shift of fermi energy in case if smearing is too large


    invert_spins
        invert spin up and spin down, now only for partial d and p


    #0 s     1 py     2 pz     3 px    4 dxy    5 dyz    6 dz2    7 dxz    8 dx2 
    #In all cases, the units of the l- and site projected DOS are states/atom/energy.

    """

    if dostype == 'partial'  :
        eld1, eld2 = {}, {}
        for i, el in enumerate(cl1.end.get_elements()):
            eld1[i+1] = el
        
        if cl2:
            for i, el in enumerate(cl2.end.get_elements()):
                eld2[i+1] = el

        if not iatom:
            printlog('Warning! Please choose atom number *iatom* from the following list:\n')
            print(eld1)
            sys.exit()
        else:
            printlog('cl1: Atom', iatom, 'of type', eld1[iatom], 'is choosen', imp = 'y')
            printlog('cl1: Atom numbers:', eld1, imp = 'y')
            printlog('cl1:', determine_symmetry_positions(cl1.end, eld1[iatom]), imp = 'y')

            # print(cl2)
            if cl2:
                if not iatom2:
                    printlog('Error! provide iatom2!')
                printlog('cl2: Atom', iatom2, 'of type', eld2[iatom2], 'is choosen', imp = 'y')


                printlog('cl2:', determine_symmetry_positions(cl2.end, eld2[iatom2]), imp = 'y')

    if iatom:
        iatom-=1
    
    if cl2:
        if not iatom2:
            printlog('Error!, provide *iatom2*!')
        iatom2-=1



    if 'figsize' not in plot_param:
        plot_param['figsize'] = (4,6)
    if 'legend' not in plot_param:
        ''
        plot_param['legend'] = 'best'

    pm = plot_param
    lw = pm.get('linewidth') or 0.8


    """1. Read dos"""
    printlog("------Start plot_dos()-----", imp = 'Y')
    

    dos = [] # main list for cl1 and cl2


    for cl in cl1, cl2:
        if cl == None: 
            continue

        if not hasattr(cl, "efermi"):
            cl.read_results('o')

        printlog(cl.name, 'e_fermi', cl.efermi, imp = 'Y')
     
        DOSCAR = cl.get_file('DOSCAR', nametype = 'asoutcar'); 
        printlog('DOSCAR file is ', DOSCAR)
        if efermi_origin:
            dos.append( VaspDos(DOSCAR, cl.efermi) )
        else:
            dos.append( VaspDos(DOSCAR, 0) )

    

    #determine number of zero energy    
    i_efermi = int(len(dos[0].energy)  *  -dos[0].energy[0] / (dos[0].energy[-1] - dos[0].energy[0])) # number of point with zero fermi energy
    if cl2:
        i_efermi_e = int(len(dos[1].energy)  *  -dos[1].energy[0] / (dos[1].energy[-1] - dos[1].energy[0])) # number of point with zero fermi energy

    


    if len(dos[0].dos) == 2:
        spin_pol = True
    else:
        spin_pol = False

    gc = None

    """2. Plot dos for different cases"""
    if dostype == 'total':
        # print(dos[0].dos)
        ylabel = "DOS (states/eV)"
        del plot_param['legend']
        if spin_pol:
            dosplot = {'Tot up':{'x':dos[0].energy,    'y':smoother(dos[0].dos[0], nsmooth), 'c':'b', 'ls':'-'}, 
                        'Tot down':{'x':dos[0].energy, 'y':-smoother(dos[0].dos[1], nsmooth),'c':'r', 'ls':'-'}}
        else:
            dosplot = {'Total':{'x':dos[0].energy, 'y':smoother(dos[0].dos, nsmooth), 'c':'b', 'ls':'-'}}


        # args[nam_down] = {'x':d.energy, 'y':-smoother(d.site_dos(iat, i_orb_down[orb]), nsmooth), 'c':color[orb], 'ls':l, 'label':None}


# xlabel = "Energy (eV)", ylabel = "DOS (states/eV)"
        # print(plot_param)
        image_name = os.path.join(path, cl1.name+'.dosTotal')
        fit_and_plot(show = show, image_name = image_name, hor = True,
            **plot_param,
            **dosplot)



    elif dostype == 'diff_total': #no spin-polarized!!!!
        ylabel = "DOS (states/eV)"

        if len(dos) > 1:    
            #calculate dos diff 
            dosd = [(d0 - d1)*e for d0, d1, e in zip(dos[0].dos, dos[1].dos, dos[0].energy)] #calculate difference
            area = trapz(dosd[:i_efermi], dx=1)
            printlog("area under dos difference = ", -area, imp = 'Y')

            fit_and_plot(show = show, image_name = cl1.name+'--'+cl2.name+'.dosTotal_Diff', xlabel = "Energy (eV)", ylabel = "DOS (states/eV)", hor = True,
                **plot_param,
                Diff_Total = (dos[0].energy, smoother(dosd, nsmooth), 'b-'))
        else:
            printlog('You provided only one calculation; could not use diff_total')





    elif 'partial' in dostype:
        #Partial dos
        #1  p carbon,  d Ti
        #0 s     1 py     2 pz     3 px    4 dxy    5 dyz    6 dz2    7 dxz    8 dx2 
       
        ylabel = "PDOS (states/atom/eV)"


        try:
            dos[0].site_dos(0, 4)
        except:
            printlog('Error! No information about partial dxy dos in DOSCAR; use LORBIT 12 to calculate it')



        """Determine neighbouring atoms """
        # printlog("Number of considered neighbors is ", neighbors, imp = 'y')

        if type(iatom) == int: #for the cases when we need to build surrounding around specific atom in this calculation - just use number of atom
            t = cl1.end.typat[iatom]
            z = cl1.end.znucl[t-1]
            el = element_name_inv(z)
            printlog('Typat of chosen imp atom in cl1 is ', el, imp = 'y')
            surround_center = cl1.end.xcart[iatom]
        else: #for the case when coordinates of arbitary point are provided.
            surround_center = iatom
            el = 'undef'

        local_atoms = local_surrounding(surround_center, cl1.end, neighbors, control = 'atoms', periodic = True)

        numbers = local_atoms[2]
        els = cl1.end.get_elements()
        els_sur = [els[i] for i in numbers]
        el_sur = ', '.join(list(set(els_sur[1:]))) # element of surrounding type
        
        printlog("Numbers of local atoms (from one):", [n+1 for n in numbers], imp = 'Y' )
        
        printlog("Elements of local atoms:", [els[i] for i in numbers], imp = 'Y' )

        printlog("List of distances", [round(d,2) for d in local_atoms[3]], imp = 'Y' )


        iX = numbers[0]# first atom is impurity if exist
        # printlog
        numbers_list = [numbers] # numbers_list is list of lists; exclude first 
        calcs = [cl1]
        

        if cl2:
            numbers_list.append([iatom2, iatom2]) # for cl2 only one atom is supported
            printlog('Warning! for cl2 p6 and d6 doesnot work and will show DOS for central atom')
            calcs.append(cl2)


        for cl, d, numbers in zip(calcs, dos, numbers_list):
            
            d.p = [] #central and surrounding
            d.d = []
            d.p_up = []
            d.p_down = []
            d.p_down = [] #central and and surrounding
            d.d_up = [] #central atom and surrounding atoms
            d.d_down = [] #central atom and surrounding atoms
            d.t2g_up = []
            d.t2g_down = []
            d.eg_up = []
            d.eg_down = []



            d.p_all = [] #sum over all atoms
            d.d_all = [] #sum over all atoms

            d.p_all_up = [] #sum over all atoms
            d.d_all_up = [] #sum over all atoms
            d.p_all_down = [] #sum over all atoms
            d.d_all_down = [] #sum over all atoms



            if 'p_all' in orbitals or 'd_all' in orbitals:
                #sum over all atoms
                p = []
                p_up = []
                p_down = []
                dd = []
                d_up = []
                d_down = []
                els = cl.end.get_elements()
                for i in range(cl.end.natom):
                    # if 'O' not in els[i]:
                    #     continue
                        
                    if spin_pol:
                        plist_up   = [d.site_dos(i, l)  for l in (2,4,6) ]
                        plist_down = [d.site_dos(i, l)  for l in (3,5,7) ]
                        plist = plist_up + plist_down

                        p_up.append(   [ sum(x)  for x in zip(*plist_up)   ] )
                        p_down.append( [ sum(x)  for x in zip(*plist_down) ] )
                        p.append(  [ sum(x) for x in zip(*plist) ] )


                    else:
                        plist = [d.site_dos(i, l)  for l in (1,2,3) ]
                        p.append(  [ sum(x) for x in zip(*plist) ] )


                    if spin_pol:
                        dlist_up   = [d.site_dos(i, l)  for l in (8,10,12,14,16) ] #
                        dlist_down = [d.site_dos(i, l)  for l in (9,11,13,15,17) ] #
                        dlist = dlist_up + dlist_down

                        d_up.append(  [ sum(x) for x in zip(*dlist_up) ]   )
                        d_down.append(  [ sum(x) for x in zip(*dlist_down) ]   )
                        dd.append(  [ sum(x) for x in zip(*dlist) ] )

                    else:
                        dlist = [d.site_dos(i, l)  for l in (4,5,6,7,8) ] #
                        dd.append(  [ sum(x) for x in zip(*dlist) ] )


                d.p_all = [ sum(pi) for pi in zip(*p) ] #sum over all atoms
                d.d_all = [ sum(di) for di in zip(*dd) ] 
                
                if spin_pol:
                    d.p_all_up   = [ sum(pi) for pi in zip(*p_up)   ] 
                    d.p_all_down = [ sum(pi) for pi in zip(*p_down) ] 
                
                    d.d_all_up   = [ sum(pi) for pi in zip(*d_up)   ] 
                    d.d_all_down = [ sum(pi) for pi in zip(*d_down) ] 
                











            #sum by surrounding atoms atoms
            n_sur = len(numbers)-1 # number of surrounding atoms
            printlog('Number of surrounding atoms:', n_sur, imp = 'Y')
            # sys.exit()
            for i in numbers: #Now for central and surrounding atoms in numbers list:

                if spin_pol:
                    plist_up   = [d.site_dos(i, l)  for l in (2,4,6) ]
                    plist_down = [d.site_dos(i, l)  for l in (3,5,7) ]
                    d.p_up.append(   [ sum(x)  for x in zip(*plist_up)   ] )
                    d.p_down.append( [ sum(x)  for x in zip(*plist_down) ] )
                    plist = plist_up + plist_down

                    d.p.append(  [ sum(x) for x in zip(*plist) ] )


                else:
                    plist = [d.site_dos(i, l)  for l in (1,2,3) ]
                    d.p.append(  [ sum(x) for x in zip(*plist) ] )



                if spin_pol:
                    dlist_up   = [d.site_dos(i, l)  for l in (8,10,12,14,16) ] #
                    dlist_down = [d.site_dos(i, l)  for l in (9,11,13,15,17) ] #
                    
                    dlist = dlist_up + dlist_down

                    d.d.append(  [ sum(x) for x in zip(*dlist) ] )


                    d.d_up.append(  [ sum(x) for x in zip(*dlist_up) ]   )
                    d.d_down.append(  [ sum(x) for x in zip(*dlist_down) ]   )

                    t2g_down = [d.site_dos(i, l)  for l in (9, 11, 15) ]
                    eg_down  = [d.site_dos(i, l)  for l in (13, 17) ]
                    
                    t2g_up   = [d.site_dos(i, l)  for l in (8, 10, 14) ]
                    eg_up    = [d.site_dos(i, l)  for l in (12, 16) ]
                    
                    d.t2g_down.append(  [ sum(x) for x in zip(*t2g_down) ]   )
                    d.eg_down.append(  [ sum(x) for x in zip(*eg_down) ]   )
                    d.t2g_up.append(  [ sum(x) for x in zip(*t2g_up) ]   )
                    d.eg_up.append(  [ sum(x) for x in zip(*eg_up) ]   )




                else:
                    dlist = [d.site_dos(i, l)  for l in (4,5,6,7,8) ] #
                    d.d.append(  [ sum(x) for x in zip(*dlist) ] )



            d.p6 = [ sum(pi[1:])/n_sur for pi in zip(*d.p) ] #sum over neighbouring atoms now only for spin up
            
            if spin_pol:
                d.p6_up   = [ sum(pi[1:])/n_sur for pi in zip(*d.p_up)   ] #sum over neighbouring atoms now only for spin up
                d.p6_down = [ sum(pi[1:])/n_sur for pi in zip(*d.p_down) ] #sum over neighbouring atoms now only for spin up
            

            d.d6 = [ sum(di[1:])/n_sur for di in zip(*d.d) ]#sum over neighbouring atoms













        """Plotting"""
        # nsmooth = 15 # smooth of dos
        d1 = dos[0]
        ds = [d1]
        names = []



        names = [cl1.id[0]+'_at_'+eld1[iatom+1]+str(iatom+1)]
        
        atoms = [iatom]
        els   = [eld1[iatom+1]]
        lts = ['-',] #linetypes
        if cl2:
            ds.append(dos[1])
            d2 = dos[1]
            
            # if labels:
            #     names.append(labels[1])
            # else:
            names.append(cl2.id[0]+'_at_'+eld2[iatom2+1]+str(iatom2+1))
            

            lts.append(lts2)
            atoms.append(iatom2)
            els.append(eld2[iatom2+1])

 

        if not spin_pol:
            plot_spin_pol = 0 # could not plot spin polarization for non-spin polarization plot

        if 'dashes' in plot_param:
            dashes = plot_param['dashes']
            del plot_param['dashes']
        else:
            dashes=(5, 1)

        dds = [(None,None), dashes, (None,None), dashes] # loop over orbitals and atoms
        # print(dds)
        # sys.exit()
        if invert_spins:
            mul = -1
        else:
            mul = 1



        energy1 = dos[0].energy
        args = {}
        if spin_pol:
            i_orb =      {'s':0, 'py':2, 'pz':4, 'px':6, 'dxy':8, 'dyz':10, 'dz2':12, 'dxz':14, 'dx2':16}
            i_orb_down = {'s':1, 'py':3, 'pz':5, 'px':7, 'dxy':9, 'dyz':11, 'dz2':13, 'dxz':15, 'dx2':17}

        else:
            i_orb = {'s':0, 'py':1, 'pz':2, 'px':3, 'dxy':4, 'dyz':5, 'dz2':6, 'dxz':7, 'dx2':8}
        # color = {'s':'k', 'p':'#F14343', 'd':'#289191', 'py':'g', 'pz':'b', 'px':'c', 'dxy':'m', 'dyz':'c', 'dz2':'k', 'dxz':'r', 'dx2':'g', 't2g':'b', 'eg':'g', 'p6':'k'}
        
        if color_dict:
            color = color_dict
        else:
            #default dict
            color = {'s':'k', 'p':'#FF0018', 'd':'#138BFF', 'py':'g', 'pz':'b', 'px':'c', 'dxy':'m', 'dyz':'c', 'dz2':'k', 'dxz':'r', 'dx2':'g', 't2g':'#138BFF', 'eg':'#8E12FF', 'p6':'#FF0018', 'p_all':'r', 'd_all':'b'} #http://paletton.com/#uid=54-100kwi++bu++hX++++rd++kX
        

        # color = {'s':'k', 'p':'r', 'd':'g', 'py':'g', 'pz':'b', 'px':'c', 'dxy':'m', 'dyz':'c', 'dz2':'m', 'dxz':'r', 'dx2':'g'}
        j = 0
        # print(orbitals)
        for orb in orbitals:
            printlog('Orbital is ', orb, imp = 'y')
            # sys.exit()
            i = 0
            for n, l, iat, el, d in zip(names, lts, atoms,els, ds):
                if el in ['Ti','Fe', 'Co', 'V', 'Mn', 'Ni'] and orb in ['p', 's', 'p_all']:
                    continue
                if el == 'O' and orb in ('d', 't2g', 'eg', 'dxy', 'dyz', 'dxz', 'dz2', 'dx2', 'd_all'):
                    continue
                nam = orb
                nam_down = nam+'_down'
                # print('name', n)
                # print('lts', l)
                if labels:
                    formula = labels[i]
                else:
                    formula = latex_chem(n.split('.')[0])

                dashes = dds[j]
                # print('dashes ',dashes,j,'\n\n\n\n\n\n\n\n')
                i+=1
                j+=1
                if spin_pol:
                    nam+=''
                suf = '; '+n
                nam+=suf
                nam_down+=suf

                printlog('Plotting for',orb, el, imp = 'y')

                # print('debug: Label is',formula, el, suf2)
                # sys.exit()
                if orb == 'p':

                    # print('debug: p is chosen')
                    if plot_spin_pol:
                        args[nam] = {'x':d.energy, 'y':mul*smoother(d.p_up[0], nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+el+suf2+' '+orb, 'dashes':dashes}

                        args[nam_down] = {'x':d.energy, 'y':mul*-smoother(d.p_down[0], nsmooth), 'c':color[orb], 'ls':l, 'label':None, 'dashes':dashes}
                        color[orb] = 'c'

                    else:
                        args[nam] = {'x':d.energy, 'y':smoother(d.p[0], nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+el+suf2+' '+orb, 'dashes':dashes}

                elif orb == 'p6':

                    # now spin-polarized components could not be shown
                    if plot_spin_pol:
                        args[nam]      = {'x':d.energy, 'y':smoother(d.p6_up, nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+el_sur+suf2+' '+orb, 'dashes':dashes}
                        args[nam_down] = {'x':d.energy, 'y':-smoother(d.p6_down, nsmooth), 'c':color[orb], 'ls':l, 'label':None, 'dashes':dashes}


                    else:
                        args[nam] = {'x':d.energy, 'y':smoother(d.p6, nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+el_sur+suf2+' p', 'dashes':dashes}




                elif orb == 'd':
                    
                    if plot_spin_pol:
                        args[nam] = {'x':d.energy, 'y':mul*smoother(d.d_up[0], nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+el+suf2+' '+orb, 'dashes':dashes}
                        args[nam_down] = {'x':d.energy, 'y':mul*-smoother(d.d_down[0], nsmooth), 'c':color[orb], 'ls':l, 'label':None, 'dashes':dashes}
                        color[orb] = 'm'

                    else:
                        args[nam] = {'x':d.energy, 'y':smoother(d.d[0], nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+el+suf2+' '+orb, 'dashes':dashes}



                elif orb == 't2g':
                    if split_type == 'octa':
                        orb_name = orb
                    elif split_type == 'tetra':
                        orb_name = 't2'

                    args[nam] = {'x':d.energy, 'y':smoother(d.t2g_up[0], nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+el+suf2+' '+orb_name, 'dashes':dashes}
                    if spin_pol:
                        args[nam_down] = {'x':d.energy, 'y':-smoother(d.t2g_down[0], nsmooth), 'c':color[orb], 'ls':l, 'label':None, 'dashes':dashes}
                
                elif orb == 'eg':
                    if split_type == 'octa':
                        orb_name = orb
                    elif split_type == 'tetra':
                        orb_name = 'e'


                    args[nam] = {'x':d.energy, 'y':smoother(d.eg_up[0], nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+el+suf2+' '+orb_name, 'dashes':dashes}
                    if spin_pol:
                        args[nam_down] = {'x':d.energy, 'y':-smoother(d.eg_down[0], nsmooth), 'c':color[orb], 'ls':l, 'label':None, 'dashes':dashes}


                elif orb == 'p_all':
                    
                    if plot_spin_pol:
                        args[nam] = {'x':d.energy, 'y':smoother(d.p_all_up, nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+suf2+' '+orb, 'dashes':dashes}
                        args[nam_down] = {'x':d.energy, 'y':-smoother(d.p_all_down, nsmooth), 'c':color[orb], 'ls':l, 'label':None, 'dashes':dashes}
                        # color[orb] = 'm'

                    else:
                        args[nam] = {'x':d.energy, 'y':smoother(d.p_all, nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+suf2+' '+orb, 'dashes':dashes}



                elif orb == 'd_all':
                    
                    if plot_spin_pol:
                        args[nam] = {'x':d.energy, 'y':smoother(d.d_all_up, nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+suf2+' '+orb, 'dashes':dashes}
                        args[nam_down] = {'x':d.energy, 'y':-smoother(d.d_all_down, nsmooth), 'c':color[orb], 'ls':l, 'label':None, 'dashes':dashes}
                        # color[orb] = 'm'

                    else:
                        args[nam] = {'x':d.energy, 'y':smoother(d.d_all, nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+suf2+' '+orb, 'dashes':dashes}






                else:
                    # args[nam] = (d.energy, smoother(d.site_dos(iat, i_orb[orb]), nsmooth), color[orb]+l)
                    # print(i_orb.keys(), color.keys())
                    # print('debug: else is chosen')

                    # sys.exit()
                    args[nam] = {'x':d.energy, 'y':smoother(d.site_dos(iat, i_orb[orb]), nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+el+suf2+' '+orb, 'dashes':dashes}
                    
                    if spin_pol:
                        args[nam_down] = {'x':d.energy, 'y':-smoother(d.site_dos(iat, i_orb_down[orb]), nsmooth), 'c':color[orb], 'ls':l, 'label':None, 'dashes':dashes}

                        # args[nam_down] = (d.energy, -smoother(d.site_dos(iat, i_orb_down[orb]), nsmooth), color[orb]+l)



        """Additional dos analysis; to be refined"""
        gc = None
        # print(plot_param['ver_lines'])
        if 'ver_lines' not in plot_param or plot_param['ver_lines'] is None:
            plot_param['ver_lines'] = []

        # print(plot_param['ver_lines'])

        if show_gravity:
            if show_gravity[0] == 1:
                d = d1
            elif show_gravity[0] == 2:
                d = d2

            if show_gravity[2]:
                erange = show_gravity[2]
            else:
                erange = (-100, 0)

            mod = show_gravity[1]

            if mod == 'p6':
                
                gc = det_gravity2(d.energy, d.p6, erange)
                # gc = det_gravity2(d.energy, d.d[0], erange)
                # printlog('Gravity center for cl1 for d for {:} is {:5.2f}'.format(erange, gc), imp = 'Y')
            

            elif show_gravity[1] == 'p':
                gc = det_gravity2(d.energy, d.p[0], erange) # for first atom for cl2

            elif show_gravity[1] == 'p_all':
                gc = det_gravity2(d.energy, d.p_all, erange) # for first atom for cl2

            elif mod == 'd':
                gc = det_gravity2(d.energy, d.d[0], erange) # for first atom for cl2

            printlog('Gravity center for cl1 for {:} for {:} is {:5.2f}'.format(mod , erange, gc), imp = 'Y')


            
            plot_param['ver_lines'].append({'x':gc, 'c':'k', 'ls':'--'})

        if efermi_origin:
            if plot_param.get('ver'):
                plot_param['ver_lines'].append({'x':efermi_shift, 'c':'k', 'ls':'-', 'lw':lw})

        else:
            #fermi levels
            plot_param['ver_lines'].append({'x':cl1.efermi + efermi_shift, 'c':'k', 'ls':'-', 'lw':lw})
            if cl2:
                plot_param['ver_lines'].append({'x':cl2.efermi + efermi_shift, 'c':'k', 'ls':'-', 'lw':lw})
        plot_param['ver'] = False
        """Plot everything"""
        if image_name is None:
            image_name = os.path.join(path, '_'.join(names)+'.'+''.join(orbitals)+'.'+el+str(iat+1))+name_suffix


        if 'xlabel' not in plot_param:
            plot_param['xlabel'] = "Energy (eV)"

        if 'ylabel' not in plot_param:
            plot_param['ylabel'] = ylabel


        fit_and_plot(show = show, image_name = image_name, hor = True,
        # title = cl1.name.split('.')[0]+'; V='+str(round(cl1.vol) )+' $\AA^3$; Impurity: '+el,
        **plot_param, 
        **args
        )
        # printlog("Writing file", image_name, imp = 'Y')




        if 0:
            """Calculate d DOS at Fermi level"""
            nn = 50 #number to integrate in both directions 
            x1 = dos[0].energy[i_efermi-nn:i_efermi+nn]
            y1 = smoother(dos[0].d6, nsmooth)[i_efermi-nn:i_efermi+nn]
            x2 = dos[1].energy[i_efermi_e-nn:i_efermi_e+nn]
            y2 = smoother(dos[1].d6, nsmooth)[i_efermi_e-nn:i_efermi_e+nn]
            f1 = interp1d(x1, y1, kind='cubic')
            f2 = interp1d(x2, y2, kind='cubic')
            # if debug: print '\n'
            # if debug: print dos[0].d6[i_efermi] - dos[1].d6[i_efermi_e], " - by points; change of d Ti DOS at the Fermi level due to the carbon"
            # if debug: print f2(0), f1(0)
            e_at_Ef_shift = f1(0) - f2(0)

            printlog("{:5.2f} reduction of d dos at Fermi level; smoothed and interpolated".format( e_at_Ef_shift ), imp = 'Y' )
        
        if 0:
            """Calculate second derivative of d at the Fermi level"""
            # tck1 = interpolate.splrep(x1, y1, s=0)
            # tck2 = interpolate.splrep(x2, y2, s=0)
            # e_at_Ef_shift_spline = interpolate.splev(0, tck1, der=0) - interpolate.splev(0, tck2, der=0)
            # if debug: print "{:5.2f} smoothed and interpolated from spline".format( e_at_Ef_shift_spline )
            # # if debug: print type(interpolate.splev(0, tck1, der=2))
            # if debug: print "{:5.2f} {:5.2f} gb and bulk second derivative at Fermi from spline".format( float(interpolate.splev(0, tck1, der=2)), float(interpolate.splev(0, tck2, der=2)) )




        if 0:
            """Calculate shift of d orbitals after adding impurity"""
            d_shift = det_gravity(dos[0],  Erange = (-2.8, 0)) - det_gravity(dos[1],  Erange = (-2.8, 0) ) #negative means impurity shifts states to negative energies, which is favourable
            printlog( "{:5.2f} Shift of Ti d center of gravity".format( d_shift ), imp = 'Y' )
            # if debug: print det_gravity(dos[0],  Erange = (-2.8, 0)), det_gravity(dos[1],  Erange = (-2.8, 0))
            """Calculate correlation between imp p and matrix d"""
            def rmsdiff(a, b):
                # rms difference of vectors a and b:
                #Root-mean-square deviation
                rmsdiff = 0
                for (x, y) in zip(a, b):
                    rmsdiff += (x - y) ** 2  # NOTE: overflow danger if the vectors are long!
                
                return math.sqrt(rmsdiff / min(len(a), len(b)))

            pd_drms = 1/rmsdiff(dos[0].p, dos[0].d6) # the higher the number the higher hybridization
            printlog("{:5.2f} p-d hybridization estimate".format( pd_drms ) , imp = 'Y')

            # if debug: print "sqeuclidean", sqeuclidean(dos[0].p, dos[0].d6)/len(dos[0].d6)

            # if debug: print "pearsonr", pearsonr(dos[0].p, dos[0].d6) #Pearson correlation coefficient; only shape; the larger number means more similarity in shape

            # def autocorr(x):
            #     result = np.correlate(x, x, mode='full')
            #     return result[result.size/2:]


        


    printlog("------End plot_dos()-----\n\n")


    return {'name':cl1.name, 'filename':image_name, 'gc':gc}
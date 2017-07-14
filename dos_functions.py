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


import header
from header import printlog
from picture_functions import fit_and_plot
from functions import element_name_inv, smoother
from geo import local_surrounding, determine_symmetry_positions








def det_gravity(dos, Erange = (-100, 0)):
    """Determine center of gravity for DOS and return values of energy for d6 orbitals in list
    INPUT:
    dos - ase dos type with added d6 - sum of d orbitals over neighbors to impurity atoms
    Erange - window of energy to determine center of gravity
    """
    sum_dos_E = {}
    sum_dos   = {}
    sum_dos['d6'] = 0
    sum_dos_E['d6'] = 0

    for i, E in enumerate(dos.energy):
        if E < Erange[0]: continue
        if E > Erange[1]: break
        sum_dos['d6']   += dos.d6[i] 
        sum_dos_E['d6'] += dos.d6[i]*E 

    d6_gc = sum_dos_E['d6']/sum_dos['d6']

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
    path = 'dos', xlim = (None, None), ylim = (None,None), savefile = True, plot_param = {}, suf2 = '', fontsize = 6 ):
    """
    cl1 (CalculationVasp) - object created by add_loop()
    dostype (str) - control which dos to plot:
        'total'   - plot total dos
        'diff_total' - difference of total dos, use cl2 for second calculation
        'partial' - partial dos

    orbitals (list of str) - any from 's, p, d, py, pz, px, dxy, dyz, dz2, dxz, dx2' where 'p' and 'd' are sums of projections
    up - 'up2' allows to download the file once again
    labels - two manual labels for cl1 and cl2 instead of auto

    iatom (int) - number of atom starting from 1 to plot DOS;
    iatom ([float]*3) - cartesian coordinates of point around which atoms will be found
    show (bool) - whether to show the dos 
    path (str)  - path to folder with images

    neighbors - number of neighbours around iatom to plot dos on them

    xlim, ylim (tuple)- limits for plot

    plot_param - dict of parameters to fit_and_plot
    suf2 - additional suffix


    #0 s     1 py     2 pz     3 px    4 dxy    5 dyz    6 dz2    7 dxz    8 dx2 
    #In all cases, the units of the l- and site projected DOS are states/atom/energy.

    """
    if fontsize:
        header.mpl.rcParams.update({'font.size': fontsize+4})
        header.mpl.rc('legend', fontsize= fontsize) 



    if dostype == 'partial'  :
        eld1, eld2 = {}, {}
        for i, el in enumerate(cl1.end.get_elements()):
            eld1[i+1] = el
        
        if cl2:
            for i, el in enumerate(cl2.end.get_elements()):
                eld2[i+1] = el

        if not iatom:
            printlog('Warning! Please choose atom number *iatom* from the following list:\n')
            printlog(eld)
            sys.exit()
        else:
            printlog('cl1: Atom', iatom, 'of type', eld1[iatom], 'is choosen', imp = 'y')
            printlog('cl1: Atom numbers:', eld1, imp = 'y')
            printlog('cl1:', determine_symmetry_positions(cl1.end, eld1[iatom]), imp = 'y')

            if cl2:
                printlog('cl2: Atom', iatom2, 'of type', eld2[iatom2], 'is choosen', imp = 'y')


                printlog('cl2:', determine_symmetry_positions(cl2.end, eld2[iatom2]), imp = 'y')

    iatom-=1
    if cl2:
        if not iatom2:
            printlog('Error!, provide *iatom2*!')
        iatom2-=1



    if 'figsize' not in plot_param:
        plot_param['figsize'] = (4,6)
    if 'legend' not in plot_param:
        plot_param['legend'] = 'best'


    """1. Read dos"""
    printlog("------Start plot_dos()-----", imp = 'Y')
    dos = []
    for cl in cl1, cl2:
        if cl == None: 
            continue

        if not hasattr(cl, "efermi"):
            cl.read_results('o')

        printlog(cl.name, 'e_fermi', cl.efermi, imp = 'Y')
     
        DOSCAR = cl.get_file('DOSCAR', update = up); 
        printlog('DOSCAR file is ', DOSCAR)

        dos.append( VaspDos(DOSCAR, cl.efermi) )
    

    #determine number of zero energy    
    i_efermi = int(len(dos[0].energy)  *  -dos[0].energy[0] / (dos[0].energy[-1] - dos[0].energy[0])) # number of point with zero fermi energy
    if cl2:
        i_efermi_e = int(len(dos[1].energy)  *  -dos[1].energy[0] / (dos[1].energy[-1] - dos[1].energy[0])) # number of point with zero fermi energy

    


    if len(dos[0].dos) == 2:
        spin_pol = True
    else:
        spin_pol = False

    """2. Plot dos for different cases"""
    if dostype == 'total':
        # print(dos[0].dos)

        if spin_pol:
            dosplot = {'Tot up':(dos[0].energy, smoother(dos[0].dos[0], 10), 'b-'), 'Tot down':(dos[0].energy, -smoother(dos[0].dos[1], 10), 'r-')}
        else:
            dosplot = {'Total':(dos[0].energy, smoother(dos[0].dos, 10), 'b-')}

        fit_and_plot(show = show, image_name = os.path.join(path, cl1.name+'.dosTotal'), xlabel = "Energy (eV)", ylabel = "DOS (states/eV)",hor = True,
            **plot_param,
            **dosplot)

    elif dostype == 'diff_total': #no spin-polarized!!!!

        if len(dos) > 1:    
            #calculate dos diff 
            dosd = [(d0 - d1)*e for d0, d1, e in zip(dos[0].dos, dos[1].dos, dos[0].energy)] #calculate difference
            area = trapz(dosd[:i_efermi], dx=1)
            printlog("area under dos difference = ", -area, imp = 'Y')

            fit_and_plot(show = show, image_name = cl1.name+'--'+cl2.name+'.dosTotal_Diff', xlabel = "Energy (eV)", ylabel = "DOS (states/eV)", hor = True,
                **plot_param,
                Diff_Total = (dos[0].energy, smoother(dosd, 15), 'b-'))
        else:
            printlog('You provided only one calculation; could not use diff_total')





    elif 'partial' in dostype:
        #Partial dos
        #1  p carbon,  d Ti
        #0 s     1 py     2 pz     3 px    4 dxy    5 dyz    6 dz2    7 dxz    8 dx2 
       
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

        numbers = local_atoms[2] # first atom is impurity if exist
        printlog("Numbers of local atoms:", [n+1 for n in numbers], imp = 'Y' )
        printlog("List of distances", [round(d,2) for d in local_atoms[3]], imp = 'Y' )


        iX = numbers[0]
        # printlog
        numbers_list = [numbers]
        if cl2:
            numbers_list.append([iatom2])



        for d, numbers in zip(dos, numbers_list):
            
            d.p = [] #central and and surrounding
            d.d = []
            d.p_down = [] #central and and surrounding
            d.d_down = [] #central atom and surrounding atoms
            d.t2g_up = []
            d.t2g_down = []
            d.eg_up = []
            d.eg_down = []
            d.d6 = 0 #sum by six atoms

            for i in numbers: #Now for surrounding atoms in numbers list:

                if spin_pol:
                    plist = [d.site_dos(i, l)  for l in (2,4,6) ]
                    plist_down = [d.site_dos(i, l)  for l in (3,5,7) ]
                    d.p_down.append( [ sum(x) for x in zip(*plist_down) ] )

                else:
                    plist = [d.site_dos(i, l)  for l in (1,2,3) ]
                d.p.append( [ sum(x) for x in zip(*plist) ] )




                if spin_pol:
                    dlist      = [d.site_dos(i, l)  for l in (8,10,12,14,16) ] #
                    dlist_down = [d.site_dos(i, l)  for l in (9,11,13,15,17) ] #
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
                
                d.d.append(  [ sum(x) for x in zip(*dlist) ]   )
            

            d.p6 = [ sum(pi) for pi in zip(*d.p) ] #sum over neighbouring atoms now only for spin up
            d.d6 = [ sum(di) for di in zip(*d.d) ] #sum over neighbouring atoms

            # t2g = [dos[0].site_dos(iTi, l)  for l in 4,5,7] #  Now only for first Ti atom
            # dos[0].t2g =  [ sum(x) for x in zip(*t2g) ]  

            # eg = [dos[0].site_dos(iTi, l)  for l in 8, 6] #  Now only for first Ti atom
            # dos[0].eg =  [ sum(x) for x in zip(*eg) ] 
           


        """Plotting"""
        nsmooth = 15 # smooth of dos
        d1 = dos[0]
        ds = [d1]
        names = []
        # if labels:
        #     names.append(labels[0])
        # else:
        names = [cl1.id[0]+'_at_'+eld1[iatom+1]+str(iatom+1)]
        
        atoms = [iatom]
        els   = [eld1[iatom+1]]
        lts = ['-',] #linetypes
        if cl2:
            ds.append(dos[1])
            
            # if labels:
            #     names.append(labels[1])
            # else:
            names.append(cl2.id[0]+'_at_'+eld2[iatom2+1]+str(iatom2+1))
            

            lts.append('-')
            atoms.append(iatom2)
            els.append(eld2[iatom2+1])

        energy1 = dos[0].energy
        args = {}
        if spin_pol:
            i_orb =      {'s':0, 'py':2, 'pz':4, 'px':6, 'dxy':8, 'dyz':10, 'dz2':12, 'dxz':14, 'dx2':16}
            i_orb_down = {'s':1, 'py':3, 'pz':5, 'px':7, 'dxy':9, 'dyz':11, 'dz2':13, 'dxz':15, 'dx2':17}

        else:
            i_orb = {'s':0, 'py':1, 'pz':2, 'px':3, 'dxy':4, 'dyz':5, 'dz2':6, 'dxz':7, 'dx2':8}
        color = {'s':'k', 'p':'#F14343', 'd':'#289191', 'py':'g', 'pz':'b', 'px':'c', 'dxy':'m', 'dyz':'c', 'dz2':'m', 'dxz':'r', 'dx2':'g', 't2g':'b', 'eg':'g'}
        # color = {'s':'k', 'p':'r', 'd':'g', 'py':'g', 'pz':'b', 'px':'c', 'dxy':'m', 'dyz':'c', 'dz2':'m', 'dxz':'r', 'dx2':'g'}

        for orb in orbitals:
            i = 0
            for n, l, iat, el, d in zip(names, lts, atoms,els, ds):
                if el in ['Fe', 'Co', 'V'] and orb == 'p':
                    continue
                if el == 'O' and orb in ('d', 't2g', 'eg', 'dxy', 'dyz', 'dxz', 'dz2', 'dx2'):
                    continue
                nam = orb
                nam_down = nam+'_down'
                print('name', n)
                print('lts', l)
                if labels:
                    formula = labels[i]
                else:
                    formula = n.split('.')[0]
                i+=1
                if spin_pol:
                    nam+=''
                suf = '; '+n
                nam+=suf
                nam_down+=suf

                if orb == 'p':
                    dashes=(5, 1)
                    # dashes=None
                    args[nam] = {'x':d.energy, 'y':smoother(d.p[0], nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+el+suf2+' '+orb}#, 'dashes':dashes}
                    if spin_pol:
                        args[nam_down] = {'x':d.energy, 'y':-smoother(d.p_down[0], nsmooth), 'c':color[orb], 'ls':l, 'label':None,}# 'dashes':dashes}
                        color[orb] = 'c'
                

                elif orb == 'd':
                    args[nam] = {'x':d.energy, 'y':smoother(d.d[0], nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+el+suf2+' '+orb}
                    if spin_pol:
                        args[nam_down] = {'x':d.energy, 'y':-smoother(d.d_down[0], nsmooth), 'c':color[orb], 'ls':l, 'label':None}
                        color[orb] = 'm'
                
                elif orb == 't2g':
                    args[nam] = {'x':d.energy, 'y':smoother(d.t2g_up[0], nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+el+suf2+' '+orb}
                    if spin_pol:
                        args[nam_down] = {'x':d.energy, 'y':-smoother(d.t2g_down[0], nsmooth), 'c':color[orb], 'ls':l, 'label':None}
                
                elif orb == 'eg':
                    args[nam] = {'x':d.energy, 'y':smoother(d.eg_up[0], nsmooth), 'c':color[orb], 'ls':l, 'label':formula+' '+el+suf2+' '+orb}
                    if spin_pol:
                        args[nam_down] = {'x':d.energy, 'y':-smoother(d.eg_down[0], nsmooth), 'c':color[orb], 'ls':l, 'label':None}

                else:
                    args[nam] = (d.energy, smoother(d.site_dos(iat, i_orb[orb]), nsmooth), color[orb]+l)
                    if spin_pol:
                        args[nam_down] = {'x':d.energy, 'y':-smoother(d.d_down[0], nsmooth), 'c':color[orb], 'ls':l, 'label':None}

                        # args[nam_down] = (d.energy, -smoother(d.site_dos(iat, i_orb_down[orb]), nsmooth), color[orb]+l)

        image_name = os.path.join(path, '_'.join(names)+'.'+''.join(orbitals)+'.'+el+str(iat+1))

        fit_and_plot(show = show, image_name = image_name, xlabel = "Energy (eV)", ylabel = "DOS (states/eV)", hor = True,
        # title = cl1.name.split('.')[0]+'; V='+str(round(cl1.vol) )+' $\AA^3$; Impurity: '+el,
        **plot_param, 
        **args
        )
        # printlog("Writing file", image_name, imp = 'Y')

        """Additional dos analysis; to be refined"""

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


    return {'name':cl1.name}
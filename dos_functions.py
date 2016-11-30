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



from header import printlog
from picture_functions import fit_and_plot
from functions import local_surrounding, element_name_inv





def smoother(x, n, mul = 1):
    #mul - additionally multiplies values
     x_smooth = []
     L = len(x)
     store = numpy.zeros((n,1),float)
     for u in range(L-n):
          for v in range(n):
               store[v] = x[u+v]
          av = float(sum(store)) / n
          x_smooth.append(av*mul)
     
     for u in range(L-n,L):
          for v in range(L-u-1):
               store[v] = x[u+v]
          av = float(sum(store)) / n
          x_smooth.append(av*mul)
     return x_smooth




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







def plot_dos(cl1, cl2 = None, dostype = None, iatom = 0, orbitals = ('s'), up = None, neighbors = 6, show = 1, path = 'dos', xlim = (None, None), ylim = (None,None) ):
    """
    cl1 (CalculationVasp) - object created by add_loop()
    dostype (str) - control which dos to plot:
        'total'   - plot total dos
        'diff_total' - difference of total dos, use cl2 for second calculation
        'partial' - partial dos

    orbitals (list of str) - any from 's, p, d, py, pz, px, dxy, dyz, dz2, dxz, dx2' where 'p' and 'd' are sums of projections
    up - 'up2' allows to download the file once again


    iatom (int) - number of atom starting from 1 to plot DOS by default it is assumed that the last atom is used
    iatom ([float]*3) - cartesian coordinates of point around which atoms will be found
    show (bool) - whether to show the dos 
    path (str)  - path to folder with images

    neighbors - number of neighbours around iatom to plot dos on them

    xlim, ylim (tuple)- limits for plot

    #0 s     1 py     2 pz     3 px    4 dxy    5 dyz    6 dz2    7 dxz    8 dx2 
    #In all cases, the units of the l- and site projected DOS are states/atom/energy.

    """
    iatom -= 1





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

    


    

    """2. Plot dos for different cases"""
    if dostype == 'total':

        fit_and_plot(show = show, image_name = os.path.join(path, cl1.name+'.dosTotal'), xlabel = "Energy (eV)", ylabel = "DOS (states/eV)",
            xlim = xlim, ylim = ylim,
            Total = (dos[0].energy, smoother(dos[0].dos, 10), 'b-'))

    elif dostype == 'diff_total':

        if len(dos) > 1:    
            #calculate dos diff 
            dosd = [(d0 - d1)*e for d0, d1, e in zip(dos[0].dos, dos[1].dos, dos[0].energy)] #calculate difference
            area = trapz(dosd[:i_efermi], dx=1)
            printlog("area under dos difference = ", -area, imp = 'Y')

            fit_and_plot(show = show,image_name = cl1.name+'--'+cl2.name+'.dosTotal_Diff', xlabel = "Energy (eV)", ylabel = "DOS (states/eV)",
                xlim = xlim, ylim = ylim,
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
        printlog("Number of considered neighbors is ", neighbors)

        if type(iatom) == int: #for the cases when we need to build surrounding around specific atom in this calculation - just use number of atom
            t = cl1.end.typat[iatom]
            z = cl1.end.znucl[t-1]
            el = element_name_inv(z)
            printlog('Typat of chosen imp atom in cl1 is ', el)
            surround_center = cl1.end.xcart[iatom]
        else: #for the case when coordinates of arbitary point are provided.
            surround_center = iatom
            el = 'undef'

        local_atoms = local_surrounding(surround_center, cl1.end, neighbors, control = 'atoms', periodic = True)

        numbers = local_atoms[2] # first atom is impurity if exist
        printlog("Numbers of local atoms:", [n+1 for n in numbers] )
        printlog("List of distances", [round(d,2) for d in local_atoms[3]] )


        iX = numbers[0]

        for j in range(len(dos)):
            
            dos[j].p = [] #central and and surrounding
            dos[j].d = [] #central atom and surrounding atoms
            dos[j].d6 = 0 #sum by six atoms

            for i in numbers: #Now for surrounding atoms in numbers list:

                plist = [dos[j].site_dos(i, l)  for l in (1,2,3) ]
                dos[j].p.append( [ sum(x) for x in zip(*plist) ] )


                dlist = [dos[j].site_dos(i, l)  for l in (4,5,6,7,8) ] # For total d T96C
                dsum = [ sum(x) for x in zip(*dlist) ] 
                dos[j].d.append(  dsum   )
            

            dos[j].p6 = [ sum(pi) for pi in zip(*dos[j].p) ] #sum over neighbouring atoms
            dos[j].d6 = [ sum(di) for di in zip(*dos[j].d) ] #sum over neighbouring atoms

            # t2g = [dos[0].site_dos(iTi, l)  for l in 4,5,7] #  Now only for first Ti atom
            # dos[0].t2g =  [ sum(x) for x in zip(*t2g) ]  

            # eg = [dos[0].site_dos(iTi, l)  for l in 8, 6] #  Now only for first Ti atom
            # dos[0].eg =  [ sum(x) for x in zip(*eg) ] 
           


        """Plotting"""
        nsmooth = 15 # smooth of dos
        d1 = dos[0]
        energy1 = dos[0].energy
        args = {}
        i_orb = {'s':0, 'py':1, 'pz':2, 'px':3, 'dxy':4, 'dyz':5, 'dz2':6, 'dxz':7, 'dx2':8}
        color = {'s':'k-', 'p':'g-', 'd':'b-', 'py':'r-', 'pz':'b-', 'px':'c-', 'dxy':'m-', 'dyz':'c-', 'dz2':'m-', 'dxz':'r-', 'dx2':'g-'}

        for orb in orbitals:
            if orb == 'p':
                args[orb] = (d1.energy, smoother(d1.p[0], nsmooth), color[orb])
            elif orb == 'd':
                args[orb] = (d1.energy, smoother(d1.d[0], nsmooth), color[orb])
            else:
                args[orb] = (d1.energy, smoother(d1.site_dos(iX, i_orb[orb]), nsmooth), color[orb])

        image_name = os.path.join(path, cl1.name+'.'+''.join(orbitals)+'.'+el+str(iX))

        fit_and_plot(show = show, image_name = image_name, figsize = (4,6), xlabel = "Energy (eV)", ylabel = "DOS (states/eV)", 
        # title = cl1.name.split('.')[0]+'; V='+str(round(cl1.vol) )+' $\AA^3$; Impurity: '+el,
        xlim = xlim, ylim = ylim, legend = 2,
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
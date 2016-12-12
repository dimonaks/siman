from __future__ import division, unicode_literals, absolute_import 
from header import printlog


def calc_redox(cl1, cl2, energy_ref = None):
    """
    Calculated average redox potential and change of volume
    cl1 (Calculation) - structure with higher concentration
    cl2 (Calculation) - structure with lower concentration
    energy_ref (float) - energy in eV per one alkali ion in anode; default value is for Li; -1.31 eV for Na, -1.02 eV for K
    """
    if cl1 is None or cl2 is None:
        printlog('cl1 or cl2 is none; return')
        return

    energy_ref_dict = {3:-1.9,  11:-1.31,  19:-1.02}
    z_alk_ions = [3, 11, 19]


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
            if not energy_ref: #take energy ref from dict
                energy_ref = energy_ref_dict[ z ]
            break

    # print(energy_ref)
    # print(cl1.energy_sigma0, cl2.energy_sigma0, mul)
    if abs(mul) > 0:
        redox = -(  ( cl1.energy_sigma0 / n1 - cl2.energy_sigma0 / n2 ) / mul  -  energy_ref  )
    else:
        redox = 0


    dV = cl1.end.vol / n1 - cl2.end.vol / n2 

    vol_red = dV / (cl1.end.vol/n1) * 100 # %

    # final_outstring = ("{:} | {:.2f} eV \n1".format(cl1.id[0]+'.'+cl1.id[1], redox  ))
    final_outstring = ("{:30} | {:10.2f} eV | {:10.1f} %".format(cl1.name, redox, vol_red  ))
    
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


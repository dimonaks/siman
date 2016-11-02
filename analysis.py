from __future__ import division, unicode_literals, absolute_import 
from header import printlog


def calc_redox(cl1, cl2, energy_ref = -1.9):
    """
    Calculated average redox potential and change of volume
    cl1 (Calculation) - structure with higher concentration
    cl2 (Calculation) - structure with lower concentration
    energy_ref (float) - energy in eV per one alkali ion in anode; default value is for Li; -1.31 eV for Na, 1.02 eV for K
    """
    cl = cl1
    bcl = cl2

    #normalize numbers of atoms by some element except Li, Na, K
    iLi = None; jLi = None
    # print cl.end.znucl
    for i, z in enumerate(cl.end.znucl):
        # print i, z
        if z in [3, 11, 19]: 
            iLi = i
            # print 'iLi is found'
            continue
        # print i, z

        for j, zb in enumerate(bcl.end.znucl):
            if zb in [3, 11, 19]: 
                jLi = j
                continue

            if z == zb:
                # print "I use ", z, " to normalize"
                i_n = i
                j_n = j


    # print "i, j",i, j
    # print 'nznucl cl',  cl.end.nznucl
    # print 'znucl cl',  cl.end.znucl
    n  = cl.end.nznucl[i_n]
    bn = bcl.end.nznucl[j_n]
    if iLi != None:
        nLi  = cl.end.nznucl[iLi]
    else:
        raise RuntimeError

    if jLi != None:
        bnLi  = bcl.end.nznucl[jLi]
    else:
        bnLi  = 0

    # print n, bn, nLi

    # print nLi/n

    mul = 1. / (float(nLi) / n)             

    # print mul


    redox = -(  ( cl.energy_sigma0 / n - bcl.energy_sigma0 / bn ) * mul  -  energy_ref  )

    dV = cl.end.vol / n - bcl.end.vol / bn 

    vol_red = dV / (cl.end.vol/n) * 100 # %

    # final_outstring = ("{:} | {:.2f} eV \n".format(cl.id[0]+'.'+cl.id[1], redox  ))
    final_outstring = ("{:30} | {:10.2f} eV | {:10.1f} %".format(cl.name, redox, vol_red  ))
    
    printlog( final_outstring )

    try:
        results_dic = {'is':cl.id[0], 'redox_pot':redox, 'id_is':cl.id, 'id_ds':bcl.id, 
        'kspacing':cl.set.kspacing, 'time':cl.time/3600.,
        'mdstep':cl.mdstep, 'ecut':cl.set.ecut, 'niter':cl.iterat/cl.mdstep,
        'set_is':cl.id[1], 'vol_red':vol_red }
    except:
        results_dic = {}


    return results_dic


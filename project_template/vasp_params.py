# the title of potential should correspond to the subfolder name
pot_pack = {'set_potential':{1:'H', 3:"Li",  5:'B', 6:'C',  8:"O", 9:"F", 11:'Na', 15:"P", 16:'S', 19:'K_pv',              22:"Ti",        23:"V", 24:'Cr',   25:"Mn",       26:"Fe",        27:"Co_new", 28:"Ni_new", 33:'As', 37:'Rb_pv', 39:'Y_sv', 56:'Ba_sv',   83:'Bi_pv', 34:'Se',    }  }
sv_pot_pack = {'set_potential':{3:"Li_sv2",    8:"O", 9:"F", 11:'Na_sv', 37:'Rb_sv', 15:"P", 16:'S', 19:'K_sv', 22:"Ti_sv_new", 23:"V_sv_new", 25:"Mn_sv",    26:"Fe_sv",     27:"Co_sv" , 28:"Ni_pv", 33:'As_d'  }} #except O_sv, which requires 1000 eV ecut at least

ion_relax_packet = {'NSW':25, 'EDIFFG':-0.025, 'EDIFF':0.0001, 'ISIF':2}
single_point_packet = {'NSW':0, 'EDIFF'     : 6e-06, 'NELM':50}
accurate_packet  = {'PREC':'Accurate', 'ADDGRID':'.TRUE.', 'EDIFF':6e-6, 'EDIFFG':-0.010, 'NSW':50, 'ISTART':1, 'ICHARG':0 }

mag_packet = {
    'GGA_COMPAT': '.FALSE.',
    'ISPIN':2,
    'LORBIT':11, #more info
    'magnetic_moments':{'Ti':0.6, 'V':5, 'Fe':5, 'Co':5, 'Mn':5, 'Ni':5, 'Cr':5 }

}

dftu_packet = {'ISTART'   :1,   'ICHARG':1,  'LDAUTYPE':2, 'LASPH':'.TRUE.', 
                'LDAUPRINT':2, 'LMAXMIX' :4, 'LDAU' :'.TRUE.',
                'LDAUL':{'Ti':2,   'Co':2  , 'Fe':2  , 'Ni':2  , 'Mn':2  , 'V':2   , 'Cr':2 },
                'LDAUU':{'Ti':0,   'Co':3.4, 'Fe':4.0, 'Ni':6.2, 'Mn':3.9, 'V':3.1 , 'Cr':3.5, 'Fe/S':1.9 },
                'LDAUJ':{'Ti':0.0, 'Co':0.0, 'Fe':0.0, 'Ni':0.0, 'Mn':0.0, 'V':0.0 , 'Cr':0.0, 'Fe/S':0   } } # universal set, Jain2011 azh values, Ni from genome

dftu_packet_off = {'LDAU' :None, 'LASPH':None, 'LDAUPRINT':None, 'LDAUTYPE':None,  'LDAUL':None, 'LDAUU':None, 'LDAUJ':None, }



dos_pack = {'NSW':0, 'LORBIT':12, 'ISMEAR':-5, 'SIGMA':None, 'LAECHG':'.TRUE.', 'EMIN':-10, 'EMAX':14, 'NEDOS':2000, 'KSPACING':0.15, 'savefile':'dox'}
bader_pack = {'PREC':'Accurate', 'ADDGRID':'.TRUE.', 'EDIFF':1e-08, 'LAECHG':'.TRUE.', 'NELM':100, 'savefile' : 'acox'}
surface_pack = {'AMIN':0.01, 'AMIX':0.2, 'BMIX':0.001, 'NELMIN':8, 'IDIPOL':3, 'LDIPOL':'.TRUE.'} # from pymatgen


#hybrid packets
hse6_pack = {'ISTART':1, 'LHFCALC':'.TRUE.', 'HFSCREEN':0.2, 'add_nbands':1.1, 'ALGO':'All', 'TIME':0.4}
hse6_pack_low = hse6_pack.copy()
hse6_pack_low.update({'PRECFOCK':'Fast', 'NKRED':2})

dftu_packet_coord = {'ISTART'   :1,   'ICHARG':1,  'LDAUTYPE':2, 'LASPH':'.TRUE.', 
                'LDAUPRINT':2, 'LMAXMIX' :4, 'LDAU' :'.TRUE.',
                'LDAUL':{'Cu':2, 'Ti':2,   'Nb':2,   'Co':2  ,   'Fe/C':2  , 'Fe/N':2, 'Ni':2  , 'Mn':2  , 'V':2   , 'Cr':2, 'Mo':2 },
                'LDAUU':{'Cu':4.0, 'Ti':0,   'Nb':1.5, 'Co':3.4, 'Fe/C':4.0, 'Fe/N':5.0, 'Ni':6.2, 'Mn':3.9, 'V':3.1 , 'Cr':3.5, 'Fe':1.9, 'Mo':3 },
                'LDAUJ':{'Cu':0.0, 'Ti':0.0, 'Nb':0.0, 'Co':0.0, 'Fe/C':0.0, 'Fe/N':0.0, 'Ni':0.0, 'Mn':0.0, 'V':0.0 , 'Cr':0.0, 'Fe':0 , 'Mo':0  } } # universal set, Jain2011 azh values, Ni from genome

mag_packet_coord = {
    'GGA_COMPAT': '.FALSE.',
    'ISPIN':2,
    'LORBIT':11, #more info
    'magnetic_moments':{'Ti':0.6, 'Nb':0.6, 'V':5, 'Fe/C':0.23, 'Fe/N':5, 'Co':5, 'Mn':5, 'Ni':5, 'Cr':5, }
}

magu_coord_pack = dftu_packet_coord.copy()
magu_coord_pack.update(mag_packet_coord)

user_vasp_sets = [
('1u_crd' ,'1u', magu_coord_pack, 'over'), #coordination_dependent
]
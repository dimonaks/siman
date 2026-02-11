# -*- coding: utf-8 -*- 
#Copyright Aksyonov D.A
#Thermochemical data

import numpy as np

from siman import header

def H2O(T, c2ev = 0, ref0K = 0):
    # enthalpy, entropy for ideal gas water!
    #https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Units=SI&Mask=1#Thermo-Gas
    # though the interval is given as 500-1700 it works quite well for for T = 300, 400 compared with original tables 
    #of Chase1998  
    #c2ev - convert to eV
    #T in K

    # ref0K - if one, then from OK, if 0 than from 298.15K

    #return in kJ/mol entalpy and Gibss and J/mol/K for entropy or in eV and eV/K
    #return G, H, S

    t  = T / 1000

    # for T in 500 - 1700 interval
    A = 30.09200   
    B = 6.832514  
    C = 6.793435   
    D = -2.534480 
    E = 0.082139    
    F = -250.8810   
    G = 223.3967    
    H = -241.8264



    dH = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H # H° − H°298.15
    S  = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
    if ref0K:
        dH += 9.904 # Chase1998 - adding enthalpy due to 298.15 K contribution

    G = dH - (T * S)/1000

    if c2ev:
        S  = header.J_mol_T2eV_T * S
        dH = header.kJ_mol2eV * dH
        G = header.kJ_mol2eV * G
        # print(G, dH - T*S )
    return G, dH, S


def O2(T, c2ev = 0, P = 1):
    """ 
    Enthalpy, entropy for gas oxygen relative to H°298.15 and at 1 bar pressure
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI&Mask=1#Thermo-Gas
    T - temperature in K
    c2ev - convert to eV
    P - pressure in bar

    return in kJ/mol entalpy and Gibss and J/mol/K for entropy or in eV and eV/K
    return G, H, S
    """

    t  = T / 1000

    if 100 < T <= 700: # 100 - 700 
        A =    31.32234    
        B =    -20.23531  
        C =    57.86644    
        D =    -36.50624 
        E =    -0.007374    
        F =    -8.903471   
        G =    246.7945     
        H =    0.0      
    elif 700 < T < 2000:
        A  =   30.03235  
        B  =   8.772972  
        C  =   -3.988133 
        D  =   0.788313  
        E  =   -0.741599 
        F  =   -11.32468 
        G  =   236.1663  
        H  =   0.0       


    dH = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H # H° − H°298.15
    S  = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G


    #influence of pressure
    S_P = -header.R * np.log(P) # P is assumed in bar
    # print(P, np.log(P))
    # print('S_P is {:.3f} J/mol/K compared to S at 1 Bar = {:.3f} J/mol/K'.format(S_P, S) )
    S +=S_P


    G = dH - (T * S)/1000 # convert to kj/mol from j/mol

    if c2ev:
        #convert to eV
        S  = header.J_mol_T2eV_T * S
        dH = header.kJ_mol2eV * dH
        G = header.kJ_mol2eV * G
        # print(G, dH - T*S )
    return G, dH, S


def H2(T, c2ev=0, P=1):
    """
    Shomate (NIST) for H2(g) 
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Table=on&Type=JANAFG
    Returns G, dH, S where:
      dH = H°(T) - H°(298.15) in kJ/mol
      S  in J/mol/K (includes pressure correction)
      G  = dH - T*S/1000 in kJ/mol
    P in bar, referenced to 1 bar.
    """
    t = T / 1000

    if 298.0 <= T <= 1000.0:
        A, B, C, D, E, F, Gc, Hc = 33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0.0
    elif 1000.0 < T <= 2500.0:
        A, B, C, D, E, F, Gc, Hc = 18.563083, 12.257357, -2.859786, 0.268238, 1.977990, -1.147438, 156.288133, 0.0
    elif 2500.0 < T <= 6000.0:
        A, B, C, D, E, F, Gc, Hc = 43.413560, -4.293079, 1.272428, -0.096876, -20.533862, -38.515158, 162.081354, 0.0
    else:
        raise ValueError("T out of supported range for H2 Shomate (298–6000 K).")

    dH = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - Hc
    S  = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + Gc

    # pressure correction (ideal gas), reference 1 bar
    S_P = -header.R * np.log(P)
    S = S + S_P

    G = dH - (T * S)/1000 # convert to kj/mol from j/mol

    if c2ev:
        S  = header.J_mol_T2eV_T * S
        dH = header.kJ_mol2eV * dH
        G  = header.kJ_mol2eV * G

    return G, dH, S


def N2(T, c2ev=0, P=1):
    """
    Shomate (NIST) for N2(g)
    https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Table=on&Type=JANAFG
    Returns G, dH, S where:
      dH = H°(T) - H°(298.15) in kJ/mol
      S  in J/mol/K (includes pressure correction)
      G  = dH - T*S/1000 in kJ/mol
    P in bar, referenced to 1 bar.
    """

    t = T / 1000.0

    if 100.0 <= T <= 500.0:
        A, B, C, D, E, F, Gc, Hc = 28.98641, 1.853978, -9.647459, 16.63537, 0.000117, -8.671914, 226.4168, 0.0
    elif 500.0 < T <= 2000.0:
        A, B, C, D, E, F, Gc, Hc = 19.50583, 19.88705, -8.598535, 1.369784, 0.527601, -4.935202, 212.3900, 0.0
    elif 2000.0 < T <= 6000.0:
        A, B, C, D, E, F, Gc, Hc = 35.51872, 1.128728, -0.196103, 0.014662, -4.553760, -18.97091, 224.9810, 0.0
    else:
        raise ValueError("T out of supported range for N2 Shomate (100–6000 K).")

    dH = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - Hc
    S  = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + Gc

    # pressure correction (ideal gas), reference 1 bar
    S_P = -header.R * np.log(P)
    S = S + S_P

    G = dH - (T * S)/1000 # convert to kj/mol from j/mol

    if c2ev:
        S  = header.J_mol_T2eV_T * S
        dH = header.kJ_mol2eV * dH
        G  = header.kJ_mol2eV * G

    return G, dH, S


def ZPE_for_molecule(molecule):
    """
    Returns the zero-point energy correction for N2, H2 molecule in eV
    Data from https://cccbdb.nist.gov
    """

    h = 4.135667696e-15  # Planck constant in eV·s
    c = 2.99792458e10    # Speed of light in cm/s

    if molecule == 'N2':
        wave_vector = 1165.0 # in cm^-1
    if molecule == 'H2':
        wave_vector = 2080.6 # in cm^-1

    zpe = h * c * wave_vector
    return zpe




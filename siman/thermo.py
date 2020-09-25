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




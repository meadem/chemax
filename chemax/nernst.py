# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 19:22:17 2022

@author: maxtm

"""

import numpy as np
import matplotlib.pyplot as plt

def Potential(E_0=0.0, T=298.15, n=1, RE=0.0, C_O1=1.0, C_R1=1.0, C_O2=1.0, C_R2=1.0, coeff_O1=1,  coeff_R1=1,  coeff_O2=1,  coeff_R2=1):
    '''
    Calculates the Nernst electrode potential for the provided arguments. Maximum of two oxidant and two reductant species allowed.
    
    Simplifies activities as concentrations.
    
    Parameters
    ----------
    E_0 : FLOAT, optional
        Standard Reduction potential (V). The default is 0.
    T : FLOAT, optional
        Temperature (K). The default is 298.15.
    n : INT, optional
        Number of electrons transferred in reaction. The default is 1.
    RE : FLOAT, optional
        Reference electrode potential v. SHE. The default is 0.
    C_O1 : FLOAT, optional
        Bulk concentration of oxidant 1, O1 (M). The default is 1.
    C_R1 : FLOAT, optional
        Bulk concentration of reductant 1, R1 (M). The default is 1.
    C_O2 : FLOAT, optional
        Bulk concentration of oxidant 2, O2. The default is 1.
    C_R2 : FLOAT, optional
        Bulk concentration of reductant 2, R2. The default is 1.
    coeff_O1 : FLOAT, optional
        Stoichiometric coefficient of oxidant 1, O1 (M). The default is 1.
    coeff_R1 : FLOAT, optional
        Stoichiometric coefficient of reductant 1, R1 (M). The default is 1.
    coeff_O2 : FLOAT, optional
        Stoichiometric coefficient of oxidant 2, O2. The default is 1.
    coeff_R2 : FLOAT, optional
        Stoichiometric coefficient of reductant 2, R2. The default is 1.
    
    Returns
    -------
    E : FLOAT
        Nernstian reduction potential relative to provided reference potential; default is v. SHE (0 V).

    '''
    F = 96485
    R = 8.3145
    E = E_0 + (R*T)/(n*F) * np.log( ( C_O1**coeff_O1 * C_O2**coeff_O2 ) / ( C_R1**coeff_R1 * C_R2**coeff_R2 ) ) - RE
    return E




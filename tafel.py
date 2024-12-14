# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 20:46:15 2022

@author: maxtm
"""

# =============================================================================
# Import packages
# =============================================================================

import numpy as np
from scipy import stats

import regression


# =============================================================================
# Define constants
# =============================================================================

F = 96485           # Faraday's constant (C/mol)
R = 8.3145          # Gas Constant (8.3145 L atm /mol /K)


# =============================================================================
# Define functions
# =============================================================================

'''
How to perform Tafel analysis using the functions defined herein:
    
# =============================================================================
#   Values used in Function calls
# =============================================================================
    E = np.asarray(potential_data)
    i = np.asarray(current_data)
    R_u = uncompensated series resistance (FLOAT)
    RE_potential = RE potential v. SHE
    E_eq = equilibrium potential v. SHE
    
# =============================================================================
#   Correcting Potential data prior to analysis (voltage unit conversion, e.g., V to mV, not included)
# =============================================================================
    E_corrected = Tafel.Correct_E_for_Ru(E, i, R_u)
    E_corrected = Tafel.Correct_E_for_RE(E_corrected, RE_potential)
    
# =============================================================================
#   Plot Tafel plot and Tafel regression
# =============================================================================
    Tafel1 = Tafel.Plot(E_corrected, I, E_eq=E_eq)
    overpotential = Tafel1[0]
    log_i = Tafel1[1]
    
    Tafel1_fit = Tafel.Regression(overpotential, log_i)
    Tafel_slope = Tafel1_fit[0]
    Tafel_yint = Tafel1_fit[1]

# =============================================================================
#   Extract Exchange Current and Exchange Rate Constant
# =============================================================================
    i_0 = 10**Tafel_yint
    OR
    i_0 = Tafel.ExchangeCurrent(E_corrected, I, E_eq=E_eq)
    
    k_0 = Tafel.ExchangeRateConstant(E_corrected, i, alpha=0.5, A=1.0, C_O=1.0, C_R=1.0, E_eq=E_eq, T=298.15)
    
    
'''


def Correct_E_for_Ru(E, i, Ru):
    E_corrected = E - i*Ru
    return E_corrected


def Correct_E_for_RE(E, RE_potential):
    E_corrected = E + RE_potential
    return E_corrected


def Region(df, E_1, E_2):
    if not E_1 < E_2:
        # check to ensure that E_1 is the smaller E value; if not, switch them
        smaller_potential = E_2
        larger_potential = E_1
        E_1 = smaller_potential
        E_2 = larger_potential
    
    # Filter data for reductive sweep of first CV cycle
    # NOTE: EC-lab appears to label oxidative sweep as 1 and reductive sweep as 0.
    # Also filter by defined potential region where current != 0
    filtered_df = df
    filtered_df = filtered_df[df['ox/red'] == 0]
    filtered_df = filtered_df[df['cycle number'] == 1.0]
    filtered_df = filtered_df[df['Ewe/V'] >= E_1]
    filtered_df = filtered_df[df['Ewe/V'] <= E_2]
    filtered_df = filtered_df[df['<I>/mA'] != 0]
    return filtered_df


def Plot(E, i, E_eq=0.0):
    '''
    Accepts potential (corrected for R_u) and current for CV data.
    Returns overpotential and log(i) arrays.

    Parameters
    ----------
    i : FLOAT
        current.
    E : FLOAT
        potential.
    E_eq : FLOAT (optional)
        equilibrium potential. Default is 0.0.

    Returns
    -------
    E : FLOAT
        Tafel potential data.
    i : FLOAT
        Tafel current data.

    '''
    overpotential = E - E_eq
    # Test functions:
    # print('-------potential-------')
    # print(E)
    # print('-------overpotential-------')
    # print(overpotential)
    log_i = np.log10(abs(i))
    return overpotential, log_i


def Regression(overpotential, log_i, step=False):
    if step == False:
        SLOPE, INTERCEPT, RVALUE, PVALUE, STDERR = stats.linregress(overpotential, log_i)
    else:
        SLOPE, INTERCEPT, RVALUE, PVALUE, STDERR = regression.ReverseStepRegression(overpotential, log_i)
    return SLOPE, INTERCEPT


def ExchangeCurrent(E_corrected, i, E_eq=0.0, step=False):
    '''
    Returns the exchange current in the same units in which the current array is provided.
    E_corrected represents the potential corrected for uncompensated series resistance, R_u.
    If E_eq is not provided it is assumed that E_corrected is v. E_eq (i.e., is overpotential)

    Parameters
    ----------
    i : FLOAT
        current, can be a current density, j.
    E : FLOAT
        potential (V).
    E_eq : FLOAT (optional)
        equilibrium potential (V).

    Returns
    -------
    i_0 : FLOAT
        exchange current, can be an exchange current density, j_0.

    '''
    overpotential, log_i = Plot(E_corrected, i, E_eq=E_eq)
    if step == False:
        slope, intercept, rvalue, pvalue, stderr = stats.linregress(overpotential, log_i)
        b, a = slope, intercept
    else:
        b, a = regression.Regression(overpotential, log_i)
    i_0 = 10**a
    return i_0


def ExchangeRateConstant(E_corrected, i, alpha=0.5, A=1.0, C_O=1.0, C_R=1.0, E_eq=0.0, T=298.15):
    '''
    Returns the exchange rate constant (cm/s). 
    E_corrected represents the potential corrected for uncompensated series resistance, R_u.
    If the current parameter is a current density, then an area, A, must also be 
    provided to return an accurate exchange rate constant.

    Parameters
    ----------
    i : FLOAT
        current, can be a current density, j.
    E : FLOAT
        potential (V).
    alpha : FLOAT, optional
        symmetry factor (unitless). The default is 0.5.
    A : FLOAT, optional
        electrode area (cm^2). The default is 1.
    C_O : FLOAT, optional
        Bulk concentration of oxidant, O (M). The default is 1.
    C_R : FLOAT, optional
        Bulk concentration of reductant, R (M). The default is 1.
    T : FLOAT, optional
        Temperature (K). The default is 298.15.

    Returns
    -------
    k_0 : FLOAT
        exchange rate constant (cm/s).

    '''
    i_0 = ExchangeCurrent(E_corrected, i, E_eq=E_eq)
    k_0 = ( i_0 ) / ( F * A * C_O**(1-alpha) * C_R**(alpha) )
    return k_0


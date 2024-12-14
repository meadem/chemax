# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 15:43:37 2022

@author: maxtm
"""

# =============================================================================
# Import packages
# =============================================================================

import numpy as np


# =============================================================================
# Define functions
# =============================================================================


def DiffusionCoefficient(slope, n=1, A=1, C=1):
    '''
    Accepts Cottrell slope. Keyword arguments for # of electrons transferred, electrode area, and bulk species concentration.
    Returns diffusion coefficient in cm^2 /s.
    
    NOTE: Cottrell slope MUST be derived from i-t data in milliamperes (mA) and seconds (s) in order to provide accurate results. (idea for expansion: add keyword argument for unit type and use case logic to calculate units appropriately.)
    
    Parameters
    ----------
    slope : FLOAT
        Cottrell slope. Must be derived from i-t data in milliamperes (mA) and seconds (s).
    n : INT, optional
        Number of electrons transferred (unitless). The default is 1.
    A : FLOAT, optional
        Area of electrode surface (cm^2). The default is 1.
    C : FLOAT, optional
        Bulk concentration of species for which the diffusion coefficient is being calculated (mol/cm^3, i.e., mol/mL). The default is 1.
    
    Returns
    -------
    D : FLOAT
        Empirically derived diffusion coefficient for species in system (cm^2 /s).
    
    '''    
    F = 96485*1000                                   # Faraday's constant (mC/mol); units of mC/mol are used due to the expectation of the Cottrell slope having been derived using mA current units.
    D = np.pi * ( slope / ( n*F*A*C ) )**2
    return D


def Plot(t, i):
    '''
    Accepts time and current data from chronoamperometry experiments.
    Returns Cottrell plot x,y data.

    Parameters
    ----------
    t : FLOAT
        time elapsed (s).
    i : FLOAT
        current response (mA).

    Returns
    -------
    t_inverse_root : FLOAT
        inverse square root of time.
    i : FLOAT
        current.

    '''
    t_inverse_root = 1/np.sqrt(t)
    return t_inverse_root, i


def Regression(data, x_limit):
    '''
    Accepts 2d-array of Cottrell x,y data and an x-limit for plotting.
    Return x,y data for regression of linear domain, defined [0, x_limit].
    
    Parameters
    ----------
    data : 2D ARRAY
        Cottrell x,y data (where x is inverse root of time and y is current).
    x_limit : INT
        Regression domain will be returned as [0, x_limit].
    
    Returns
    -------
    x_domain : 1D ARRAY
        domain of regression to plot.
    regression : 1D ARRAY
        regression values as solved for using the specified x_domain.
    '''
    
    x = data[0]
    y = data[1]
    
    # Find point at which current is current.max / e
    regression_domain_max = x[y <= (y.max() / np.e)].max()
    
    x_linear = x[x <= regression_domain_max]
    y_linear = y[len(y) - len(x_linear):]
    
    # Fit curve to linear domain of Cottrell plot
    slope, intercept = np.polyfit(x_linear, y_linear, 1)
    
    # Define domain over which to plot regression values
    x_domain = np.linspace(0, x_limit, x_limit*2)
    
    # Plot curve over specified data domain (not just linear region)
    regression = slope*x_domain + intercept
    
    return x_domain, regression


def Slope(data):
    '''
    Accepts 2d-array of Cottrell x,y data.
    Return Cottrell slope, defined by the provided domain max.
    (Function is essentially the first half of the Regression function, above. Seems redundant, but it's what I need. Could probably be more elegant?)
    
    Parameters
    ----------
    data : 2D ARRAY
        Cottrell x,y data (where x is inverse root of time and y is current).
    regression_domain_max : INT
        The highest value of x (inverse root t) to which the regression should be fit.
    
    Returns
    -------
    slope : FLOAT
        Cottrell slope.
    '''
    
    x = data[0]
    y = data[1]
    
    regression_domain_max = x[y <= (y.max() / np.e)].max()
    
    x_linear = x[x <= regression_domain_max]
    y_linear = y[len(y) - len(x_linear):]
    
    # Fit curve to linear domain of Cottrell plot
    slope, intercept = np.polyfit(x_linear, y_linear, 1)
    
    return slope


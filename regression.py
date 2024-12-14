# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 2022

@author: maxtm
"""

# =============================================================================
# Import packages
# =============================================================================

import numpy as np
from scipy import stats


# =============================================================================
# Define functions
# =============================================================================

def IsLinear(R_squared, alpha):
    '''
    Provides criteria for defining the linear region. Compares the provided R^2 value to the tolerated error (alpha).
    
    Parameters
    ----------
    R_squared : FLOAT
        The coefficient of determination of the linear regression (aka r-squared value).
    alpha : FLOAT
        Error tolerance of regression.

    Returns
    -------
    bool
        Returns True (i.e., regression is linear) if the r-squared value is greater than or equal to alpha.
        Return False (i.e., regression is not linear) if the r-squared value is less than alpha.
    '''
    if abs(R_squared) < alpha:
        return False
    else:
        return True


def StepRegression(x, y, certainty=0.98):
    '''
    Step through a 2d data set, starting with the first 2 data points only and incorporating one additional datum per step.
    Check regression linearity for each step using IsLinear function.
    When the coefficient of determination of the regression falls below the certainty argument, return the most recent regression stats.

    Parameters
    ----------
    x : ARRAY
        Pandas series or numpy array of the independent variable data.
    y : ARRAY
        Pandas series or numpy array of the dependent variable data.
    certainty : FLOAT, optional
        The alpha value to which the coefficient of determination is compared as a linearity criterion. The default is 0.98.

    Returns
    -------
    SLOPE : FLOAT
        Slope of linear region.
    INTERCEPT : FLOAT
        Y-intercept of linear region.
    RVALUE : FLOAT
        Coefficient of determination of linear region.
    PVALUE : FLOAT
        P-value of linear region.
    STDERR : FLOAT
        Standard error of linear region.
    '''
    SLOPE = 0.0         #default
    INTERCEPT = 0.0     #default
    RVALUE = 1.0        #default
    PVALUE = 0.1        #default
    STDERR = 0.1        #default
    for i in np.arange(2, (len(x)-1)):
        # Define x-data indices to iterate through using np.arange and the length of the x-data array.
        # Iterate through the x-data, calculating regression stats for each successively larger domain.
        #print(i) # test function
        slope, intercept, rvalue, pvalue, stderr = stats.linregress(x[:i], y[:i])
        if IsLinear(rvalue, certainty):
            # As long as the linearity test passes, update the stats and continue iterating
            SLOPE = slope
            INTERCEPT = intercept
            RVALUE = rvalue
            PVALUE = pvalue
            STDERR = stderr
        else:
            # When linearity test fails, return the regression values, which have not changed since the last successful linearity test.
            return SLOPE, INTERCEPT, RVALUE, PVALUE, STDERR


def ReverseStepRegression(x, y, certainty=0.98):
    '''
    Step through a 2d data set, starting with the LAST 2 data points only and incorporating one additional datum per step.
    Works backwards through data set relative to StepRegression function behavior.
    Check regression linearity for each step using IsLinear function.
    When the coefficient of determination of the regression falls below the certainty argument, return the most recent regression stats.

    Parameters
    ----------
    x : ARRAY
        Pandas series or numpy array of the independent variable data.
    y : ARRAY
        Pandas series or numpy array of the dependent variable data.
    certainty : FLOAT, optional
        The alpha value to which the coefficient of determination is compared as a linearity criterion. The default is 0.98.

    Returns
    -------
    SLOPE : FLOAT
        Slope of linear region.
    INTERCEPT : FLOAT
        Y-intercept of linear region.
    RVALUE : FLOAT
        Coefficient of determination of linear region.
    PVALUE : FLOAT
        P-value of linear region.
    STDERR : FLOAT
        Standard error of linear region.
    '''
    SLOPE = 0.0         #default
    INTERCEPT = 0.0     #default
    RVALUE = 1.0        #default
    PVALUE = 0.1        #default
    STDERR = 0.1        #default
    for i in np.flip(np.arange(0, (len(x)-2))):
        # Define x-data indices to iterate through backwards using np.flip on np.arange and the length of the x-data array.
        # Iterate through the x-data, calculating regression stats for each successively larger domain.
        #print(i) # test function
        slope, intercept, rvalue, pvalue, stderr = stats.linregress(x[i:], y[i:])
        if IsLinear(rvalue, certainty):
            SLOPE = slope
            INTERCEPT = intercept
            RVALUE = rvalue
            PVALUE = pvalue
            STDERR = stderr
        else:
            # When linearity test fails, return the regression values, which have not changed since the last successful linearity test.
            return SLOPE, INTERCEPT, RVALUE, PVALUE, STDERR
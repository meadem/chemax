# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 08:29:29 2022

@author: maxtm
"""

###############################################
# Import packages
###############################################

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# The display function is great for displaying dataframes. Use it to redo table returns in RDE project.
# from IPython.display import display, HTML

# =============================================================================
# Define functions (package == levitch)
# =============================================================================

def LimitingCurrent(dataframe_dictionary, redox='ox'):
    '''
    Accepts a dictionary formatted as {key:dataframe} (from load() function).
    Returns a dictionary with the same keys and each df replaced with its respective limiting current.

    Parameters
    ----------
    dataframe_dictionary : DICT
        Dictionary of {key:dataframe} pairs. 
        Dataframes are expected to be cyclic voltammogram data.
    redox : STRING, optional
        Either 'red' or 'ox' to specify the oxidative (max) or reductive (min) limiting current. 
        The default is 'ox.'

    Raises
    ------
    ValueError
        redox parameter is limited to values 'ox' or 'red.'

    Returns
    -------
    current_dictionary : DICT
        A dictionary with the same keys as the provided dict and each df replaced with its respective limiting current.
    '''
    VALID_REDOX = {'ox', 'red'}
    if redox not in VALID_REDOX:
        raise ValueError("redox must be one of %r." % VALID_REDOX)
    
    current_dictionary = {}
    for key in dataframe_dictionary:
        if redox == 'ox':
            # Get limiting oxidative (+) current value for each experiment
            current_dictionary[key] = dataframe_dictionary[key].max()
        elif redox == 'red':
            # Get limiting reductive (-) current value for each experiment
            current_dictionary[key] = dataframe_dictionary[key].min()

    return current_dictionary


def Potentials(dataframe_dictionary, potentials):
    # Accepts a dictionary of CV dataframes {experiment#:df} and an array or list of potentials. 
    # potentials example: potentials = np.asarray([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    # Cyclic voltammetry dataframes are intended to be at different disc rotation speeds.
    # Returns a df organized as df[potential][CV_data_at_different_rotation_rates]. (for levitch plots)
    
    import copy.deepcopy
    dataframe_dictionary_copy = copy.deepcopy(dataframe_dictionary)
    # The two lines above ensure operations on the returned object do not effect the original object.
    # I am not sure if otherwise the filtered_df_dictionary would contain references to dataframe_dictionary, or independent copies of that data.
    
    sampled_data = {}
    for df in dataframe_dictionary_copy:
        CV = dataframe_dictionary_copy[df]
        for E in potentials:
            CV_slice = CV[CV['Ewe/V'] <= E][CV['Ewe/V'] >= E-0.05]
            CV_datum = CV_slice[CV_slice['Ewe/V'] == CV_slice['Ewe/V'].max()]
            sampled_data[E].append(CV_datum)
    filtered_df_dictionary = pd.DataFrame(sampled_data)
    
    return filtered_df_dictionary


def Regression(levitch_lines, potential, angular_frequency):
    # Experiment order needs to line up with angular frequency order.
    # NOTE: Some regressions will need to be performed on a sub-domain that is linear.
    # How to use for plotting:
    #
    # levitch_lines = levitch.Potentials(dfs, potentials)
    # for potential in potentials:
    #     current = levitch_lines[potential]['<I>/mA']
    #     levitch_x, levitch_y = levitch.Regression(levitch_lines, potential, angular_frequency)   
    #     ax1.scatter(1/np.sqrt(angular_frequency), 1/current, label=str(potential) + ' V')
    #     ax2.plot(levitch_x, levitch_y)
    #
    levitch_current = levitch_lines[potential]['<I>/mA']
    regression_y = 1/levitch_current
    regression_x = 1/np.sqrt(angular_frequency)
    SLOPE, YINT = np.polyfit(regression_x, regression_y, 1)
    domain = np.insert(regression_x, 0, 0)
    y = SLOPE*domain + YINT
    return domain, y


###############################################
# Extract diffusion coefficient values from regressions and format for display (cell 1)
###############################################

def DiffusionCoefficient():
    pass

n = 1                # Number of electrons transferred
F = 96485            # Faraday's constant (C/mol)
d = 0.5              # Disk diameter (cm)
A = np.pi*(d/2)**2   # Disk area (cm^2)
C1= 0.010            # Bulk concentration of O in cell 1 (M)
C2= 0.0              # Bulk concentration of O in cell 2 (M)
v = 1.01e-6          # Kinematic viscosity (m^2 s^-1) (provided in instructions)

cell1_D_o_E7 = ( (1/E7_slope) / (0.62*n*F*A*(v**(-1/6))*C1))**(3/2)
cell1_D_o_E8 = ( (1/E8_slope) / (0.62*n*F*A*(v**(-1/6))*C1))**(3/2)
cell1_D_o = ( cell1_D_o_E7 + cell1_D_o_E8 ) / 2


###############################################
# Extract kinetic current values from regressions (cell 1)
###############################################

def KineticCurrent():
    pass

cell1_i_k_E4 = 1 / E4_yint
cell1_i_k_E5 = 1 / E5_yint
cell1_i_k_E6 = 1 / E6_yint
cell1_i_k_E7 = 1 / E7_yint
cell1_i_k_E8 = 1 / E8_yint
cell1_i_k = np.asarray([cell1_i_k_E4, cell1_i_k_E5, cell1_i_k_E6, cell1_i_k_E7, cell1_i_k_E8])


###############################################
# Tafel analysis (Cell 1)
###############################################

# This should be done using the Tafel package
Tafel_slope1, Tafel_yint1 = np.polyfit(E_list[3:], np.log10(abs(cell1_i_k)), 1)
def TafelRegression_cell1(x):
    y = Tafel_slope1*x + Tafel_yint1
    return y

cell1_i0 = 10**Tafel_yint1


###############################################
# Plot CV data and Levitch plots (Cell 1)
###############################################

# Define x range for Levitch plot(s); 0 included here using np.insert() but not in regression
Levitch_X = np.insert(1/(np.sqrt(angular_frequency)),0,0)

# Define figure 1 object
fig1 = plt.figure(1, figsize=(16, 13))

# Top-Left plot
ax1 = fig1.add_subplot(2,2,1)
ax1.plot(experiment5['Ewe/V'], experiment5['<I>/mA'], label='2000 Hz')
ax1.plot(experiment4['Ewe/V'], experiment4['<I>/mA'], label='1600 Hz')
ax1.plot(experiment3['Ewe/V'], experiment3['<I>/mA'], label='900 Hz')
ax1.plot(experiment2['Ewe/V'], experiment2['<I>/mA'], label='400 Hz')
ax1.plot(experiment1['Ewe/V'], experiment1['<I>/mA'], label='100 Hz')
ax1.set_xlabel(r'${E}_{WE} \: (V)$', fontsize = 16)
ax1.set_ylabel(r'$I \: (mA)$', fontsize = 16)
ax1.axhline(y=0, color='k')
ax1.axvline(x=0, color='k')
ax1.legend()

# Top-Right plot
ax2 = fig1.add_subplot(2,2,2)
ax2.scatter(1/np.sqrt(angular_frequency), 1/E1, label='0.1 V')
ax2.scatter(1/np.sqrt(angular_frequency), 1/E2, label='0.2 V')
ax2.scatter(1/np.sqrt(angular_frequency), 1/E3, label='0.3 V')
ax2.scatter(1/np.sqrt(angular_frequency), 1/E4, label='0.4 V')
ax2.scatter(1/np.sqrt(angular_frequency), 1/E5, label='0.5 V')
ax2.scatter(1/np.sqrt(angular_frequency), 1/E6, label='0.6 V')
ax2.scatter(1/np.sqrt(angular_frequency), 1/E7, label='0.7 V')
ax2.scatter(1/np.sqrt(angular_frequency), 1/E8, label='0.8 V')
ax2.set_xlabel(r'${\omega}^{-\frac{1}{2}}$', fontsize = 16)
ax2.set_ylabel(r'$\frac{1}{i} \: ({mA}^{-1})$', fontsize = 16)
ax2.axhline(y=0, color='k')
ax2.axvline(x=0, color='k')
ax2.legend()

# Bottom-Left plot
ax3 = fig1.add_subplot(2,2,3)
ax3.scatter(1/np.sqrt(angular_frequency), 1/E4, label='0.4 V')
ax3.scatter(1/np.sqrt(angular_frequency), 1/E5, label='0.5 V')
ax3.scatter(1/np.sqrt(angular_frequency), 1/E6, label='0.6 V')
ax3.scatter(1/np.sqrt(angular_frequency), 1/E7, label='0.7 V')
ax3.scatter(1/np.sqrt(angular_frequency), 1/E8, label='0.8 V')
ax3.plot(Levitch_X, E4_regression(Levitch_X))
ax3.plot(Levitch_X, E5_regression(Levitch_X))
ax3.plot(Levitch_X, E6_regression(Levitch_X))
ax3.plot(Levitch_X, E7_regression(Levitch_X))
ax3.plot(Levitch_X, E8_regression(Levitch_X))
ax3.set_xlabel(r'${\omega}^{-\frac{1}{2}}$', fontsize = 16)
ax3.set_ylabel(r'$\frac{1}{i} \: ({mA}^{-1})$', fontsize = 16)
ax3.axhline(y=0, color='k')
ax3.axvline(x=0, color='k')
ax3.legend()

# Bottom-Right plot
ax4 = fig1.add_subplot(2,2,4)
ax4.scatter(E_list[3:], np.log10(abs(cell1_i_k)))
ax4.plot(E_list[3:], TafelRegression_cell1(E_list[3:]))
ax4.set_xlabel(r'${E}_{WE} \: (V)$', fontsize = 16)
ax4.set_ylabel(r'$log|{i}_{k}| \: (i \: in \: mA)$', fontsize = 16)
ax4.axhline(y=0, color='k')
ax4.axvline(x=0, color='k')

# Caption and Figure formatting
caption1 = "Figure 1. Top-Left: Cyclic Voltammograms of ferri-/ferro-cyanide cell collected between -0.2 V and +0.8 V v. Ag/AgCl reference \n electrode with a glassy carbon RDE and Pt CE in 18 mL 200 mM KCl(aq). Redox pair both present at 10 mM in situ and cell was \n purged w/ N2 gas. Top-Right: Inverse current v. inverse root of angular frequency for ferri-/ferro-cyanide cell at selected potentials. \n Bottom-Left: Selected lines from Levitch plot from Top-Right graph with their respective regressions. Bottom-Right: Tafel \n plot of kinetic current values at several potentials. \n\n"
fig1.text(.5, -0.07, caption1, fontsize=14, ha='center');
fig1.suptitle(r'${{Fe(CN)}_{6}}^{3-} + {e}^{-} \rightleftharpoons {{Fe(CN)}_{6}}^{4-}$', fontsize=24);

print('\n\n')



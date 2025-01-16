# CHEMAX
Working on expanding, refactoring, and combining several pieces of code I wrote as a student for electrochemical data analysis. 
I've since taught myself python Classes, and my goal is to build a class object that can do the following with minimal typing on the user end:
- Import and hold an arbitrarily large number of experiment data files
- Manipulate data, e.g., iR-correct voltage, or calculate current density automatically when current and area are known
- Extract useful values, like series resistances, diffusion coefficients, and Tafel slopes
- plot well-formatted curves with minimal user input (CV, IV, nyquist, bode, Tafel, Cottrell)

All of the files here are planned for eventual merging and inclusion into the chemax file.

Below is an example of how to use chemax. In the future I'd like to expand this readme to explain more thoroughly, but for now these demos will have to suffice.

--------------------------------------

```python

###############################################
# Import packages
###############################################

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import chemax


###############################################
# Import data
###############################################

FeCN = chemax.Experiment(name='FeCN', 
                         description='CV and CA data for FeCN chemistry')

FeCN.load(20, technique='CA', file='CA_data_', folder='data/wk1', step=2, filetype='.txt', silent=True)
FeCN.load(7, technique='CV', file='CV_data_', folder='data/wk1', filetype='.txt', silent=True)
FeCN.load(3, technique='CV sim', file='FeCN_', folder='data/CVsim', filetype='.txt', silent=True)


FeSO4 = chemax.Experiment(name='FeSO4', 
                          description='CV and CA data for FeSO4 chemistry')

FeSO4.load(8, technique='CA', file='CA_FeSO4_', folder='data/wk2', filetype='.txt', silent=True)
FeSO4.load(6, technique='CV', file='CV_FeSO4_', folder='data/wk2', filetype='.txt', silent=True)
FeSO4.load(3, technique='CV sim', file='FeSO4_', folder='data/CVsim', filetype='.txt', silent=True)


FePhen = chemax.Experiment(name='FePhen', 
                           description='CV and CA data for FePhen chemistry')

FePhen.load(8, technique='CA', file='CA_FePhen_', folder='data/wk2', filetype='.txt', silent=True)
FePhen.load(6, technique='CV', file='CV_FePhen_', folder='data/wk2', filetype='.txt', silent=True)
FePhen.load(3, technique='CV sim', file='FePhen_', folder='data/CVsim', filetype='.txt', silent=True)


###############################################
# Cottrell analysis
###############################################

FeCN.cottrell()
FeSO4.cottrell()
FePhen.cottrell()


###############################################
# Plot Fe(CN)
###############################################

fig1 = plt.figure(1, figsize=(16,9))
FeCN.plot("CA", technique="CA", position=221, current_units="mA")
FeCN.plot("Cottrell", technique="CA", position=222, current_units="mA")
FeCN.plot("IV", technique="CV", position=223, current_units="mA")
FeCN.plot("IV", technique="CV sim", position=224, current_units="mA")

caption1 = "Figure 1. Top-Left: Chronoamperometry (potential step) experiment data; Top-Right: Cottrell plots of CA data; \n Bottom-Left: experimental cyclic voltammetry data; Bottom-Right: simulated cyclic voltammetry data."
fig1.text(0.15, -0.03, caption1, fontsize=14, va='bottom', ha='left', wrap=True);
fig1.suptitle(r'${{Fe(CN)}_{6}}^{3-/4-} \ Redox \ Couple$', fontsize=24);
#print('\n\n')


###############################################
# Plot Fe2+
###############################################

fig2 = plt.figure(2, figsize=(16,9))
FeSO4.plot("CA", technique="CA", position=221, current_units="mA")
FeSO4.plot("Cottrell", technique="CA", position=222, current_units="mA")
FeSO4.plot("IV", technique="CV", position=223, current_units="mA")
FeSO4.plot("IV", technique="CV sim", position=224, current_units="mA")

caption2 = "Figure 2. Top-Left: Chronoamperometry (potential step) experiment data; Top-Right: Cottrell plots of CA data; \n Bottom-Left: experimental cyclic voltammetry data; Bottom-Right: simulated cyclic voltammograms."
fig2.text(0.15, -0.03, caption2, fontsize=14, va='bottom', ha='left', wrap=True);
fig2.suptitle(r'${{Fe}_{(aq)}}^{2+/3+} \ Redox \ Couple$', fontsize=24);
#print('\n\n')


###############################################
# Plot Fe(phen)
###############################################

# Define figure object
fig3 = plt.figure(3, figsize=(16,9))
FePhen.plot("CA", trials=[3,4,5], position=221, current_units="mA")
FePhen.plot("Cottrell", technique="CA", position=222, current_units="mA")
FePhen.plot("IV", technique="CV", position=223, current_units="mA")
FePhen.plot("IV", technique="CV sim", position=224, current_units="mA")

caption3 = "Figure 3. Top-Left: Chronoamperometry (potential step) experiment data; Top-Right: Cottrell plots of CA data; \n Bottom-Left: experimental cyclic voltammetry data; Bottom-Left: simulated cyclic voltammograms."
fig3.text(0.15, -0.03, caption3, fontsize=14, va='bottom', ha='left', wrap=True);
fig3.suptitle(r'${{Fe(phen)}_{3}}^{2+/3+} \ Redox \ Couple$', fontsize=24);
#print('\n\n')
```
![Untitled](https://github.com/user-attachments/assets/6ff634f5-50d5-42a5-aa4a-c6a7f6754c2e)

![Untitled-1](https://github.com/user-attachments/assets/662690a3-3033-4d3f-bd57-8c874c58d95f)

![Untitled](https://github.com/user-attachments/assets/c3782611-cacc-481e-a5e1-e040e0c45ea9)



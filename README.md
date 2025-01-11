# CHEMAX
Working on expanding, refactoring, and combining several pieces of code I wrote as a student for electrochemical data analysis. 
I've since taught myself python Classes, and my goal is to build a class object that can do the following with minimal typing on the user end:
- Import and hold an arbitrarily large number of experiment data files
- Manipulate data, e.g., iR-correct voltage, or calculate current density automatically when current and area are known
- Extract useful values, like series resistances, diffusion coefficients, and Tafel slopes
- plot well-formatted curves with minimal user input (CV, IV, nyquist, bode, Tafel, Cottrell)

All of the files here are planned for eventual merging and inclusion into the chemax file.

Below are some demonstrations of how to use chemax, taken from projects I did as a student where I've replaced the code with updated code. 
In the future I'd like to expand this readme to explain more thoroughly, but for now these demos will have to suffice.

--------------------------------------

```python
import numpy as np
import matplotlib.pyplot as plt
import chemax

fig1 = plt.figure()

CV = chemax.Experiment()
CV.load(num=2, filetype='.csv', file='101422_MM_CV_', folder='data')

CV.IV(label="250 mV/s scan rate CV", position=111, trials=[1])
CV.IV(label="500 mV/s scan rate CV", position=111, trials=[2])
```
![Untitled](https://github.com/user-attachments/assets/f9f4d216-ac7e-4aeb-a632-9595fa10d556)

(this CV data was acquired using a potentiostat I made on a breadboard, which explains the resolution)

--------------------------------------

```python
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import chemax


exp = chemax.Experiment()
exp.load(18, file='101822_', filetype='.txt')

fig2 = plt.figure(2, figsize=(14,7))
exp.IV(label='500 mV/s', trials=[1], position=121, trim=[250,None])
exp.IV(label='100 mV/s', trials=[2], position=121, trim=[100,None])
exp.IV(label='10 mV/s', trials=[3], position=121, trim=[25,None])
exp.IV(label='250 mV/s', trials=[6], position=122, trim=[50,None])
caption2 = "Figure 2. CVs collected in 0.2 M aqueous KCl using a GC electrode swept from an initial voltage of 0 to +1 V to -1 V \n v. SCE (left) and a Pt electrode swept from an initial voltage of 0 to +2 V to -2 V v. SCE (right). Data corrected to \n exclude aliasing effects."
fig2.text(0.1, -0.08, caption2, fontsize=14, va='bottom', ha='left', wrap=True)

fig3 = plt.figure(3, figsize=(15, 10))
exp.IV(label='10 mV/s', trials=[12])
exp.IV(label='50 mV/s', trials=[15])
exp.IV(label='500 mV/s', trials=[7])
caption3 = "Figure 3. Cyclic voltammograms of GC WE in moderately stirred 0.2 M KCl and 10 mM ferrocyanide at \n varied voltage scan rates. Experiments started at 0 V, swept to -2 V then swept positive to 1 V."
fig3.text(.5, -0.01, caption3, fontsize=14, ha='center')

fig5 = plt.figure(5, figsize=(15, 10))
exp.IV(label='10 mM R; no O', trials=[11])
exp.IV(label='10 mM R; 5 mM O', trials=[16])
exp.IV(label='10 mM R; 10 mM O', trials=[17])
caption5 = "Figure 5. Cyclic voltammogram of glassy carbon working electrode in vigorously stirred 0.2 M KCl and different relative \n proportions of R and O."
fig5.text(.5, -0.01, caption5, fontsize=14, ha='center')

fig6 = plt.figure(6, figsize=(15, 10))
exp.IV(label='Stirred', trials=[17])
exp.IV(label='Unstirred', trials=[18])
caption6 = "Figure 6. Cyclic voltammogram of glassy carbon working electrode in moderately stirred 0.2 M KCl and 10 mM \n ferrocyanide/10 mM ferricyanide, both stirred and unstirred."
fig6.text(.5, -0.01, caption6, fontsize=14, ha='center')
```
![Untitled](https://github.com/user-attachments/assets/5b5c029d-2982-4092-9662-5eb0a7e134ac)

![Untitled](https://github.com/user-attachments/assets/295a0599-0824-4461-b35b-4815f5f8a78b)

![Untitled](https://github.com/user-attachments/assets/75b0d1f3-60c3-4e1b-9cfe-3d05ea735f25)

![Untitled](https://github.com/user-attachments/assets/1ee2efe1-12d0-4443-b15d-ef1764e9a699)


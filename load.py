# -*- coding: utf-8 -*-
"""
Created on Fri Dec 30 15:50:42 2022

@author: maxtm
"""


# =============================================================================
# Import packages
# =============================================================================

import os
import numpy as np
import pandas as pd


# =============================================================================
# Define functions
# =============================================================================

def files(num, file='data', folder='data', step=1):
    '''
    Unpacks data in txt files and creates a dictionary to hold their respective pandas dataframe objects.
    Assumes tab-delimited data (should be updated in the future with another keyword argument).
    
    Example use:
        from chempy import load
        CA_FeCN = load.files(20, file='CA_data_', folder='data/wk1', step=2)
        print(CA_FeCN[1])

    Parameters
    ----------
    num : INT
        Number of data files to load (not all files need be loaded; see the [step] keyword argument.
    file : STR, optional
        File name. 
        This function assumes files are named similarly with an incremented counter at the end. 
        The default is 'data'.
    folder : STR, optional
        Path to data folder. 
        The default is 'data' and assumes folder is in same directory as notebook.
    step : INT, optional
        Step size to use for np.arange (i.e., when iterating through files). 
        If not set =1 then some files will be skipped (i.e., if set =2 then even numbered files are skipped). 
        If every other scan is a relaxation this should be set to 2. 
        The default is 1.

    Returns
    -------
    One dictionary object with key:value pairs {file_number:dataframe}.
    '''
    file_numbers = np.arange(1, num+1, step)
    dictionary = {}
    for file_number in file_numbers:
        file_name = file + str(file_number) + '.txt'
        file_path = os.path.join(folder, file_name)
        dictionary[file_number] = pd.read_csv(file_path, delimiter='\t', encoding='unicode_escape')
    return dictionary


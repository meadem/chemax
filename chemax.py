import os
from datetime import datetime
import numbers
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt




# REFERENCE ELECTRODE DICTIONARY
REFERENCE_ELECTRODES = {
    # "Name": potential v. SHE,
      "SHE" : 0.0,
      "Ag/AgCl" : 0.2225,
      "SCE" : 0.242,
    }




class Experiment():
    '''
    # in plots that have voltage on an axis, it would be nice to label the axis with the reference electrode (e.g., V v. SHE...)
    
    # should I default to current units of mA or A?
    
    # Regression() is a duplicate with linear_function(), which is newer; need to compare and downselect
    
    # BUG: the various input() queries do not print out in a consistent order relative to other print-outs (maybe import datetime?)
    
    # would be nice if tafel() method printed something out when it ran, like the tafel plot at least and maybe also the corrected I-V curve
    
    # add method for resetting the data (and/or undoing the most recent load operation?..)
    
    # add try catch in multiple places where it's missing
    
    '''
    
    def __init__(self, name=None, description=None, RE=None, display_RE=None, area=None, Ru=None):
        '''
        name = experiment name (could be same as variable used for class instance. Mostly used for error and print msgs.)
        description = experiment description
        RE = Reference Electrode used to acquire the data
        display_RE = Reference Electrode used to correct the voltage for display
        area = working electrode area
        Ru = uncompensated series resistance
        '''
        self.history = {}
        self.update_history("Experiment object initialized.")
        self.name = name
        self.description = description
        self.RE = RE
        self.REFERENCE_ELECTRODE_POTENTIAL = None
        self.DISPLAY_REFERENCE_ELECTRODE_POTENTIAL = None
        self.area = area
        self.Ru = Ru
        self.technique = {}
        
        
        # Set the reference electrode potential, for both the RE used to acquire the data and the one used to display it
        if RE == None:
            print(f'Reference electrode not specified for {self.name}.')
        elif RE in REFERENCE_ELECTRODES:
            for entry in REFERENCE_ELECTRODES:
                if RE == entry:
                    self.REFERENCE_ELECTRODE_POTENTIAL = REFERENCE_ELECTRODES[entry]
        else:
            print(f'Specified reference electrode {self.RE} not found for {self.name}.'
                  'Unexpected behavior may occur when using an unrecognized reference electrode.')
        
        
        # Set the display RE as the data acquisition RE, if display RE has not been specified
        if display_RE is None:
            self.display_RE = RE
        else:
            self.display_RE = display_RE
            
        # Print message if area is not specified.
        if area == None:
            print(f'Working electrode area not specified for {self.name}.')
    
    
    # Setting these attributes below as 'properties' ensures the values are updated if needed when other attributes are changed.
    # e.g. it allows the object history to update when an attribute is changed.
    
    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, arg):
        self._name = arg
        self.update_history(f"Name set to {arg}.")
    
    
    @property
    def description(self):
        return self._description
    @description.setter
    def description(self, arg):
        self._description = arg
        self.update_history(f"Description set to {arg}.")
    
    
    @property
    def RE(self):
        return self._RE
    @RE.setter
    def RE(self, arg):
        self._RE = arg
        self.update_history(f"RE set to {arg}.")
    
    
    @property
    def display_RE(self):
        return self._display_RE
    @display_RE.setter
    def display_RE(self, arg):
        self._display_RE = arg
        if hasattr(self, 'data'):
            self.correct_voltage("RE", RE=arg)
        self.update_history(f"Display RE set to {arg}.")
    
    
    @property
    def area(self):
        return self._area
    @area.setter
    def area(self, arg):
        self._area = arg
        self.update_history(f"Area set to {arg}.")
    
    
    @property
    def Ru(self):
        return self._Ru
    @Ru.setter
    def Ru(self, arg):
        self._Ru = arg
        self.update_history(f"Ru set to {arg}.")
    
    
    @property
    def tafel_slope_V(self):
        return self._tafel_slope_V
    @tafel_slope_V.setter
    def tafel_slope_V(self, arg):
        self._tafel_slope_V = arg
        self.update_history(f"Tafel slope (in Volts) set to {arg}.")
    
    
    @property
    def tafel_slope_mV(self):
        return self._tafel_slope_mV
    @tafel_slope_mV.setter
    def tafel_slope_mV(self, arg):
        self._tafel_slope_mV = arg
        self.update_history(f"Tafel slope (in mV) set to {arg}.")
    
    
    @property
    def tafel_yint(self):
        return self._tafel_yint
    @tafel_yint.setter
    def tafel_yint(self, arg):
        self._tafel_yint = arg
        self.update_history(f"Tafel y-intercept set to {arg}.")
    
    
    @property
    def i_0(self):
        return self._i_0
    @i_0.setter
    def i_0(self, arg):
        self._i_0 = arg
        self.update_history(f"Exchange current (i_0) set to {arg}.")
        if self.area is not None:
            self.j_0 = self.i_0 / self.area
    
    
    @property
    def j_0(self):
        return self._j_0
    @j_0.setter
    def j_0(self, arg):
        self._j_0 = arg
        self.update_history(f"Exchange current density (j_0) set to {arg}.")
    
    
    
    
    def cottrell(self, trials=[], n=1, C=1):
        '''
        perform cottrell analysis. current converted to mA for analysis
        '''
        F = 96485
        for trial in trials:
            # Check if Corrected Time for trial data has already been calculated; do so if not (is this check even necessary?)
            # NOTE: it may be more accurate to correct the start time w/o zeroing, similar to how I did this originally (can do easily??)
            if "Corrected Time" not in self.data[trial]:
                self.correct_time(trials=[trial])

            # Check if t^-0.5 for trial data has already been calculated; do so if not (is this check even necessary?)
            if "Time-Inverse-Root" not in self.data[trial]:
                self.data[trial]["Time-Inverse-Root"] = 1 / np.sqrt(self.data[trial]["Corrected Time"])
                self.update_history(f"Time-Inverse-Root calculated from Corrected Time for trial {trial} in {self.name}.")

            # filter out data where current <= peak current / e
            TIME_INVERSE_ROOT = self.data[trial]["Time-Inverse-Root"]
            CURRENT = self.data[trial]["Current"]*1000
            FILTERED_CURRENT = CURRENT[CURRENT <= (CURRENT.max() / np.e)]
            COTTRELL_INDEX_MIN = FILTERED_CURRENT.idxmax()
            COTTRELL_DOMAIN = TIME_INVERSE_ROOT[COTTRELL_INDEX_MIN:]
            COTTRELL_RANGE = CURRENT[COTTRELL_INDEX_MIN:]

            # linear regression
            COTTRELL_SLOPE, COTTRELL_INTERCEPT = np.polyfit(COTTRELL_DOMAIN, 
                                                            COTTRELL_RANGE, 
                                                            1)

            # save as attributes of trial's dataframe object
            self.data[trial].COTTRELL_INDEX_MIN = COTTRELL_INDEX_MIN
            self.data[trial].COTTRELL_SLOPE = COTTRELL_SLOPE
            self.data[trial].COTTRELL_INTERCEPT = COTTRELL_INTERCEPT
            
            # (cm^2 /s) (area expected in cm^2)
            # Bulk conc of species for which the diffusion coefficient is being calculated has units of (mol/cm^3, i.e., mol/mL, or MM)
            DIFFUSION_COEFFICIENT = np.pi * ( COTTRELL_SLOPE / ( n*F*self.area*concentration*1000 ) )**2
            print(f"A diffusion coefficient of {DIFFUSION_COEFFICIENT} was extracted using Cottrell analysis from {self.name} trial {trial}")
    
    
    
    
    def correct_time(self, trials=[], silent=True):
        '''
        Zero-correct time (i.e., shift all time values equally and such that the first value is zero)
        '''
        try:
            for trial in trials:
                self.data[trial]["Corrected Time"] = self.data[trial]["Time"] - self.data[trial]["Time"][0]
                self.update_history(f"Corrected Time calculated from raw Time data for trial {trial} in {self.name}.")
        except:
            raise Exception(f"Error while zero-correcting time for {self.name}.")
    
    
    
    
    def correct_voltage(self, correction=None, trials=[], silent=True, RE=None):
        '''
        Adds a new column to the specified trials' dataframe(s) containing the voltage corrected for either:
        a) uncompensated series resistance (new column referenced as self.data[trial]["iR-Corrected Voltage"])
        b) display reference electrode (new column referenced as self.data[trial]["RE-Corrected Voltage"])
        
        Adds another new column referenced as self.data[trial]["Corrected Voltage"]
        This column is one of the following:
        1) the iR-Corrected Voltage (if this is the only correction performed on the voltage so far)
        2) the RE-Corrected Voltage (if this is the only correction performed on the voltage so far)
        3) Voltage corrected for both iR and display RE (if both corrections have been performed)
        '''
        if not trials:
            trials = self.data
        
        if correction == None:
            raise Exception(f"Please specify the voltage correction to be made for {self.name}: 'iR' or 'RE'")
        
        elif correction == "iR":
            try:
                for trial in trials:
                    if ("Voltage" in self.data[trial]) and ("Current" in self.data[trial]):
                        self.data[trial]["iR-Corrected Voltage"] = self.data[trial]["Voltage"] - ( self.data[trial]["Current"] * self.Ru )
                        self.update_history(f"'iR-Corrected Voltage' calculated (Ru = {self.Ru:.3} ohm) for trial {trial} .")
                        if not silent:
                            print(f"Voltage successfully iR-corrected (Ru = {self.Ru:.3} ohm) for {self.name} trial {trial}.")
                        # Check if data was previousy RE-corrected, and re-perform that correction if needed
                        if "RE-Corrected Voltage" in self.data[trial].columns:
                            ADJUSTMENT = self.REFERENCE_ELECTRODE_POTENTIAL - self.DISPLAY_REFERENCE_ELECTRODE_POTENTIAL
                            self.data[trial]["Corrected Voltage"] = self.data[trial]["iR-Corrected Voltage"] + ADJUSTMENT
                            self.update_history(f"'Corrected Voltage' has been adjusted for Ru = {self.Ru:.3} ohm for trial {trial}.")
                        else:
                            self.data[trial]["Corrected Voltage"] = self.data[trial]["iR-Corrected Voltage"]
                            self.update_history(f"'Corrected Voltage' set equal to 'iR-Corrected Voltage' for trial {trial}.")
            except:
                raise ValueError(f"Error while attempting to iR correct voltage for {self.name} trial {trial}.")
        
        elif correction == "RE":
            try:
                # Check first if RE used to acquire data is the one being requested for display
                if RE == self.RE:
                    if not silent:
                        print(f"{RE} is already the reference electrode implicitly referenced by the acquired voltage data.")
                else:
                    # Fetch the RE potential from its name using the dictionary, or print an error
                    if RE in REFERENCE_ELECTRODES:
                        for entry in REFERENCE_ELECTRODES:
                            if RE == entry:
                                DISPLAY_REFERENCE_ELECTRODE_POTENTIAL = REFERENCE_ELECTRODES[entry]
                    else:
                        print(f"Specified reference electrode {RE} not found.")
                    
                    # Calculate the adjustment to be made
                    ADJUSTMENT = self.REFERENCE_ELECTRODE_POTENTIAL - DISPLAY_REFERENCE_ELECTRODE_POTENTIAL
                    
                    # apply the adjustment to all specified trials
                    for trial in trials:
                        if "Voltage" in self.data[trial].columns:
                            self.data[trial]["RE-Corrected Voltage"] = self.data[trial]["Voltage"] + ADJUSTMENT
                            self.update_history(f"Voltage corrected by {ADJUSTMENT:.4} V for RE '{RE}' for trial {trial}.")
                            if not silent:
                                # BUG: this message prints even for trials where no correction was performed (need to fix)
                                print(f"Voltage successfully RE-corrected by {ADJUSTMENT:.4} V for RE '{RE}' for {self.name} trial {trial}.")
                            if "iR-Corrected Voltage" in self.data[trial].columns:
                                self.data[trial]["Corrected Voltage"] = self.data[trial]["iR-Corrected Voltage"] + ADJUSTMENT
                                self.update_history(f"'Corrected Voltage' has been adjusted for RE '{RE}' for trial {trial}.")
                            else:
                                self.data[trial]["Corrected Voltage"] = self.data[trial]["RE-Corrected Voltage"]
                                self.update_history(f"'Corrected Voltage' set equal to 'RE-Corrected Voltage' for trial {trial}.")
                    
                    self.DISPLAY_REFERENCE_ELECTRODE_POTENTIAL = REFERENCE_ELECTRODES[RE]
            
            except:
                raise Exception(f"Error while attempting to RE correct voltage for {self.name} trial {trial}.")
        
    
    
    
    
    def update_history(self, entry):
        '''
        Append a datetime:entry entry to the object's .history dict object
        '''
        self.history[datetime.today().strftime('%Y-%m-%d %H:%M:%S.%f')] = entry
    
    
    
    
    def load(self, num=1, step=1, filetype='.txt', file='data', folder='data', file_numbers=None, technique=None, silent=False):
        '''
        Upload data files into an Experiment object as a pandas dataframe.
        
        NEXT: allow for this to be called more than once, and to append data files (instead of overwriting)
        '''
        # tag if this operation is an instantiation (alternative is an append operation)
        _IS_INSTANTIATION = not hasattr(self, "data")
        
        # set up and execute file loading
        try:
            # create temporary dicts that will be incorporated into self.source and self.data at the end of the method call 
            source = {}
            data = {}

            # create a temporary dict that holds key:value pairs of INTERNAL:EXTERNAL data indices during import
            # used to force consecutive integer values for internal data reference
            # e.g., files "data_1" and "data_3" would be indexed as 1 and 2, as would files "CV_data_1" and "RDE_data_1"
            # ----------------------------------------------------------------------------------------------------------
            EXTERNAL_DATA_INDICES = None
            INTERNAL_DATA_INDICES = None
            
            # EXTERNAL DATA INDICES
            # if provided, use file_numbers argument in place of the num argument for defining the external data indices
            if file_numbers is None:
                EXTERNAL_DATA_INDICES = np.arange(1, num+1, step)
            else:
                EXTERNAL_DATA_INDICES = file_numbers
            
            # NUMBER OF NEW TRIALS
            NUMBER_OF_NEW_TRIALS = len(EXTERNAL_DATA_INDICES)
            
            # NUMBER OF EXISTING TRIALS
            NUMBER_OF_EXISTING_TRIALS = None
            if _IS_INSTANTIATION:
                NUMBER_OF_EXISTING_TRIALS = 0
            else:
                NUMBER_OF_EXISTING_TRIALS = len(self.data)
            
            # INTERNAL DATA INDICES
            if _IS_INSTANTIATION:
                INTERNAL_DATA_INDICES = np.arange(1, 
                                                  NUMBER_OF_NEW_TRIALS+1)
            else:
                INTERNAL_DATA_INDICES = np.arange(NUMBER_OF_EXISTING_TRIALS + 1, 
                                                  NUMBER_OF_EXISTING_TRIALS + NUMBER_OF_NEW_TRIALS + 1)
            
            # DATA INDICES
            DATA_INDICES = dict(zip(INTERNAL_DATA_INDICES,
                                     EXTERNAL_DATA_INDICES))
            # ----------------------------------------------------------------------------------------------------------
            
            # iteratively load files
            for i in DATA_INDICES:
                INTERNAL_DATA_INDEX = i                     
                EXTERNAL_DATA_INDEX = DATA_INDICES[i]      
                FILE_NAME = file + str(EXTERNAL_DATA_INDEX) + filetype
                FILE_PATH = os.path.join(folder, FILE_NAME)
                source[INTERNAL_DATA_INDEX] = FILE_PATH
                data[INTERNAL_DATA_INDEX] = pd.read_csv(FILE_PATH, sep=None, engine='python', encoding='unicode_escape')
                if not silent:
                    print(f"Data file {FILE_PATH} successfully imported (INDEX={INTERNAL_DATA_INDEX}) for {self.name}.")
        except:
            raise Exception(f"Error while importing data files from {folder} into {self.name}.")
        
        # tag data with the acquisition technique, if it was provided
        try:
            if technique is not None:
                if technique not in self.technique.keys():
                    self.technique[technique] = list(DATA_INDICES.keys())
                else:
                    self.technique[technique] = sorted(set(self.technique[technique] + list(DATA_INDICES.keys())))
        except:
            raise Exception(f"Error while tagging data files imported from {folder} with technique {technique} for {self.name}.")
        
        # standardize imported data headers
        try:
            # List of standardized column names
            VOLTAGE_HEADER = "Voltage"
            CURRENT_HEADER = "Current"
            TIME_HEADER = "Time"
            CYCLE_HEADER = "Cycle"
            POWER_HEADER = "Power"
            FREQUENCY_HEADER = "Frequency"
            ZREAL_HEADER = "Z_real"
            ZIMAG_HEADER = "Z_imag"
            Z_HEADER = "Z"
            PHASE_ANGLE_HEADER = "Phase Angle"

            # A dictionary of headers I've observed previously and what to replace them with when found.
            HEADER_DICT = {"VOLTAGE":VOLTAGE_HEADER,
                           "EWE/V":VOLTAGE_HEADER,
                           "<EWE>/V":VOLTAGE_HEADER,
                           "E (V)":VOLTAGE_HEADER,
                           "EWE/MV":VOLTAGE_HEADER,
                           "E (MV)":VOLTAGE_HEADER,
                           "CURRENT":CURRENT_HEADER,
                           "<I>/A":CURRENT_HEADER,
                           "<I>/MA":CURRENT_HEADER,
                           "I/MA":CURRENT_HEADER,
                           "I/A":CURRENT_HEADER,
                           "TIME (S)":TIME_HEADER,
                           "TIME":TIME_HEADER,
                           "TIME/S":TIME_HEADER,
                           "CYCLE":CYCLE_HEADER,
                           "CYCLE NUMBER":CYCLE_HEADER,
                           "POWER":POWER_HEADER,
                           "POWER (W)":POWER_HEADER,
                           "P (W)":POWER_HEADER,
                           "P/W":POWER_HEADER,
                           "FREQ/HZ":FREQUENCY_HEADER,
                           "RE(Z)/OHM":ZREAL_HEADER,
                           "-IM(Z)/OHM":ZIMAG_HEADER,
                           "|Z|/OHM":Z_HEADER,
                           "PHASE(Z)/DEG":PHASE_ANGLE_HEADER,
                          }

            # Apply header standardization (also: convert current data from mA to A)
            for sheet in data:
                for column_name in data[sheet]:
                    for entry in HEADER_DICT:
                        if column_name.upper() == entry:
                            if entry == "<I>/MA":
                                data[sheet][column_name] /= 1000
                                if not silent:
                                    print(f"Converted current from mA to A for experiment {sheet}.")
                            data[sheet].rename(columns={column_name: HEADER_DICT[entry]}, inplace=True)
        
        except:
            raise Exception(f"Error during header standardization while importing data from {folder} for {self.name}.")
        
        # save data and data sources in namespace
        # These needs to become append statements
        try:
            if _IS_INSTANTIATION:
                self.data = data
                self.source = source
            else:
                self.data.update(data)
                self.source.update(source)
        except:
            raise Exception(f"Error while saving data imported from {folder} for {self.name}.")
        
        # If a display reference electrode has been specified previously, prompt whether it should be applied during loading
        try:
            if self.display_RE is not None:
                check = input(f"Correct voltages for {self.name} using display RE '{self.display_RE}'? "
                              f"Values will be saved in a new column. (Y/N)")
                if check.upper() == "Y":
                    self.correct_voltage("RE", RE=self.display_RE, silent=silent)
                elif check.upper() == "N":
                    pass
                else:
                    print("Input not recognized. "
                          "Change the display_RE attribute to RE correct the voltage for display.")
        except:
            raise Exception(f"Error while prompting/correcting for RE voltage adjustment for {self.name}.")
        
        # Print a success message
        if not silent:
            print(f"Data successfully imported from '{folder}'")
    
    
    
    
    def extract_Ru_from_HFR(self, trial=None, silent=False):
        '''
        Extract uncompensated series resistance, R_u, from EIS data by extracting the minimum real impedance value.
        '''
        try:
            # Extract Ru value (and its index) from the real impedance data
            result = self.data[trial]["Z_real"].min()
            result_index = np.argmin(self.data[trial]["Z_real"])
            
            # Nyquist plot zoomed in around the extracted value
            print("-------------------------------------------------------------------------------------")
            plt.figure(figsize=([10,4]))
            self.plot("Nyquist", position=121, trials=[trial], legend=False)
            self.plot("Nyquist", position=122, trials=[trial], legend=False, zoom=[15, result_index])
            plt.gcf().supxlabel(f"The Ru value extracted from the impedance data for {self.name} trial {trial} is {result} ohms.", 
                                fontsize=12, y=-0.1)
            plt.show(block=False)
            
            # Check extracted value with user
            check = input("Accept this value? (Y/N)")
            print("-------------------------------------------------------------------------------------")
            if check.upper() == "N":
                print("REJECTED")
                self.update_history(f"Failed to extract uncompensated series resistance from {self.name} trial {trial} impedance data.")
            elif check.upper() == "Y":
                print("ACCEPTED")
                self.Ru = result
                self.update_history(f"Extracted uncompensated series resistance "
                                    f"of {result} ohms from trial {trial} impedance data.")
            else:
                print(r"Please respond with 'Y' or 'N' to ACCEPT or REJECT the extracted value of ${R}_{u}$.")
                print("-------------------------------------------------------------------------------------")
                self.extract_Ru_from_HFR(trial=trial)
        except:
            raise KeyError("Trial not found, EIS data not available, or other error occurred while attempting to extract R_u from EIS data.")
        
        # query user if they want to apply extracted Ru value to iR-correct the voltage data, then do so if requested
        try:
            unresolved = True
            while unresolved:
                check = input("iR-correct the voltage data? "
                              "It will be saved into a new column 'iR-Corrected Voltage' without overwriting 'Voltage'. (Y/N)")
                print("-------------------------------------------------------------------------------------")
                if check.upper() == "N":
                    unresolved = False
                elif check.upper() == "Y":
                    self.correct_voltage("iR", silent=silent)
                    unresolved = False
                else:
                    print(r"Please respond with 'Y' to iR-correct the voltage data, or 'N' to not correct it.")
                    print("-------------------------------------------------------------------------------------")
        except:
            raise Exception("Error while iR-correcting voltage data for {self.name}.")
    
    
    
    
    def plot(self, type_, title="", label="", trials=[], position=111, trim=[None,None], zoom=[None,None],
             xlabel=None, ylabel=None, legend=True, voltage_units="V", current_units="A"):
        '''
        Plot the specified curve type. 
        
        Can be used to plot multiple traces at once, but this is probably only useful when labeling is not necessary.
        '''
        
        if not trials:
            trials = self.data
        
        FILTER_START = trim[0]
        FILTER_STOP = trim[1]
        
        plt.gcf()
        plt.subplot(position)
        plt.title(title)
        plt.axhline(color='k', linewidth=0.25)
        plt.axvline(color='k', linewidth=0.25)
        
        if type_ == "IV":
            
            if xlabel == None:
                xlabel = "Potential (V)"
            if ylabel == None:
                ylabel="Current (A)"
            
            plt.xlabel(xlabel, fontsize=16)
            plt.ylabel(ylabel, fontsize=16)
            
            for trial in trials:
                # Use RE-Corrected Voltage if available
                VOLTAGE = None
                CURRENT = self.data[trial]["Current"]
                if "RE-Corrected Voltage" in self.data[trial].columns:
                    VOLTAGE = self.data[trial]["RE-Corrected Voltage"]
                else:
                    VOLTAGE = self.data[trial]["Voltage"]
                # Plot using appropriate units and labeling
                if voltage_units.upper() == "V" and current_units.upper() == "A":
                    plt.plot(VOLTAGE[FILTER_START:FILTER_STOP], 
                             CURRENT[FILTER_START:FILTER_STOP], 
                             label=label)
                elif voltage_units.upper() == "V" and current_units.upper() == "MA":
                    plt.ylabel("Current (mA)")
                    plt.plot(VOLTAGE[FILTER_START:FILTER_STOP], 
                             CURRENT[FILTER_START:FILTER_STOP]*1000, 
                             label=label)
                elif voltage_units.upper() == "MV" and current_units.upper() == "A":
                    plt.xlabel("Voltage (mV)")
                    plt.plot(VOLTAGE[FILTER_START:FILTER_STOP]*1000, 
                             CURRENT[FILTER_START:FILTER_STOP], 
                             label=label)
                elif voltage_units.upper() == "MV" and current_units.upper() == "MA":
                    plt.xlabel("Voltage (mV)")
                    plt.ylabel("Current (mA)")
                    plt.plot(VOLTAGE[FILTER_START:FILTER_STOP]*1000, 
                             CURRENT[FILTER_START:FILTER_STOP]*1000, 
                             label=label)
                else:
                    raise NameError("Units not found or other error while plotting I-V curve.")
        
        if type_ == "IV-corrected":
            
            if xlabel == None:
                xlabel = "Corrected Potential (V)"
            if ylabel == None:
                ylabel="Current (A)"
            
            plt.xlabel(xlabel, fontsize=16)
            plt.ylabel(ylabel, fontsize=16)
            
            for trial in trials:
                try:
                    # Check if iV-corrected voltage has already been calculated; do it if it hasn't been
                    if "iR-Corrected Voltage" not in self.data[trial]:
                        self.correct_voltage("iR", trials=[trial])
                    
                    # Plot corrected voltage using the specified current and voltage units
                    # Note that voltage plots from "Corrected Voltage" and not "iR-Corrected Voltage"; 
                    # this incorporates the RE adjustment if one has been made, and otherwise is just the iR-corrected voltage anyways
                    if voltage_units.upper() == "V" and current_units.upper() == "A":
                        plt.plot(self.data[trial]["Corrected Voltage"][FILTER_START:FILTER_STOP], 
                                 self.data[trial]["Current"][FILTER_START:FILTER_STOP], 
                                 label=label)
                    elif voltage_units.upper() == "V" and current_units.upper() == "MA":
                        plt.ylabel("Current (mA)", fontsize=16)
                        plt.plot(self.data[trial]["Corrected Voltage"][FILTER_START:FILTER_STOP], 
                                 self.data[trial]["Current"][FILTER_START:FILTER_STOP]*1000, 
                                 label=label)
                    elif voltage_units.upper() == "MV" and current_units.upper() == "A":
                        plt.xlabel("Corrected Voltage (mV)", fontsize=16)
                        plt.plot(self.data[trial]["Corrected Voltage"][FILTER_START:FILTER_STOP]*1000, 
                                 self.data[trial]["Current"][FILTER_START:FILTER_STOP], 
                                 label=label)
                    elif voltage_units.upper() == "MV" and current_units.upper() == "MA":
                        plt.xlabel("Voltage (mV)", fontsize=16)
                        plt.ylabel("Current (mA)", fontsize=16)
                        plt.plot(self.data[trial]["Corrected Voltage"][FILTER_START:FILTER_STOP]*1000, 
                                 self.data[trial]["Current"][FILTER_START:FILTER_STOP]*1000, 
                                 label=label)
                    else:
                        raise NameError("Units not found or other error while plotting I-V curve.")
                except:
                    # it would be nice if this and other error readouts were more thorough, maybe including the object name
                    raise ValueError(f"Error while making corrected IV plot for trial {trial}.")
        
        if type_ == "Tafel":
            
            if xlabel == None:
                xlabel = xlabel = r'$\eta \: (V)$'
            if ylabel == None:
                ylabel = r'$log|i| \: (i \ in \ A)$'
            
            plt.xlabel(xlabel, fontsize=16)
            plt.ylabel(ylabel, fontsize=16)
            
            xlim = [None,None]
            ylim = [None,None]
            
            for trial in trials:
                try:
                    X = self.data[trial]["Tafel overpotential"]
                    Y = self.data[trial]["Tafel current"]
                    
                    # Plot log(i) v. iR-corrected overpotential
                    plt.scatter(X, 
                                Y,
                                label=label)
                    
                    # Plot linear regression of above data
                    plt.plot(X, 
                             self.linear_function(X, 
                                                  self.tafel_slope_V, 
                                                  self.tafel_yint))
                    
                    # Check and readjust plotting window for each trace
                    if xlim[0] == None or xlim[0] > X.min():
                        xlim[0] = X.min()
                    if xlim[1] == None or xlim[1] < X.max():
                        xlim[1] = X.max()
                    if ylim[0] == None or ylim[0] > Y.min():
                        ylim[0] = Y.min()
                    if ylim[1] == None or ylim[1] < Y.max():
                        ylim[1] = Y.max()
                    
                except:
                    raise ValueError(f"Error while making Tafel plot for trial {trial}.")
            
            # Apply 10% margins to axis limits
            MARGIN = 0.10
            X_RANGE = abs(xlim[1] - xlim[0])
            Y_RANGE = abs(ylim[1] - ylim[0])
            xlim = [xlim[0] - X_RANGE*MARGIN,
                    xlim[1] + X_RANGE*MARGIN]
            ylim = [ylim[0] - Y_RANGE*MARGIN,
                    ylim[1] + Y_RANGE*MARGIN]
            
            # Set axis limits
            plt.xlim(xlim)
            plt.ylim(ylim)
            
        
        if type_ == "Cottrell":
            
            if xlabel == None:
                xlabel = r"${t}^{-\frac{1}{2}} \: ({s}^{-\frac{1}{2}})$"
            if ylabel == None:
                ylabel = r"$Current \: (mA)$"
            
            plt.xlabel(xlabel, fontsize=16)
            plt.ylabel(ylabel, fontsize=16)
            
            try:
                # list of plotting colors (for matching scatter color with regression color)
                COLOR = plt.rcParams['axes.prop_cycle'].by_key()['color']
                
                for trial in trials:
                    COTTRELL_INDEX_MIN = self.data[trial].COTTRELL_INDEX_MIN
                    X = self.data[trial]["Time-Inverse-Root"][COTTRELL_INDEX_MIN:]
                    Y = self.data[trial]["Current"][COTTRELL_INDEX_MIN:]*1000
                    
                    # Plot time-inverse-root v. current (mA)
                    NUMBER_OF_EXISTING_SCATTER_TRACES = len(plt.gca().collections)
                    COLOR_INDEX = NUMBER_OF_EXISTING_SCATTER_TRACES
                    COLOR = str(COLOR[COLOR_INDEX])
                    plt.scatter(X, Y, label=label, color=COLOR)
                    
                    # Plot linear regression of same
                    plt.plot(X, 
                             self.linear_function(X, 
                                                  self.data[trial].COTTRELL_SLOPE, 
                                                  self.data[trial].COTTRELL_INTERCEPT)
                            , '--', color=COLOR)
                    
            except:
                raise ValueError(f"Error while making Cottrell plot for trial {trial}.")
        
        
        if type_ == "CA":
            
            if xlabel == None:
                xlabel = r"$Time \: (s)$"
            if ylabel == None:
                ylabel = r"$Current \: (A)$"
            
            plt.xlabel(xlabel, fontsize=16)
            plt.ylabel(ylabel, fontsize=16)
            
            for trial in trials:
                TIME = None
                CURRENT = self.data[trial]["Current"]
                if "Corrected Time" in self.data[trial].columns:
                    TIME = self.data[trial]["Corrected Time"]
                else:
                    self.correct_time(trials=[trial])
                    TIME = self.data[trial]["Corrected Time"]
                # Plot using appropriate units and labeling
                if current_units.upper() == "A":
                    plt.plot(TIME[FILTER_START:FILTER_STOP], 
                             CURRENT[FILTER_START:FILTER_STOP], 
                             label=label)
                elif current_units.upper() == "MA":
                    plt.ylabel("Current (mA)")
                    plt.plot(TIME[FILTER_START:FILTER_STOP], 
                             CURRENT[FILTER_START:FILTER_STOP]*1000, 
                             label=label)
                else:
                    raise NameError("Units not found or other error while plotting CA data.")
        
        
        if type_ == "Nyquist":
            
            if xlabel == None:
                xlabel = r"${Z}_{re} \: (\Omega)$"
            if ylabel == None:
                ylabel = r"${-Z}_{im} \: (\Omega)$"
            
            plt.xlabel(xlabel, fontsize=16)
            plt.ylabel(ylabel, fontsize=16)
            
            for trial in trials:
                plt.scatter(self.data[trial]["Z_real"][FILTER_START:FILTER_STOP],
                            self.data[trial]["Z_imag"][FILTER_START:FILTER_STOP],
                            label=label)
            
            if (zoom[0] is not None) and (zoom[1] is not None):
                ZOOM_RATIO = zoom[0]
                CENTER_INDEX = zoom[1]
                
                X_CENTER = self.data[trial]["Z_real"][CENTER_INDEX]
                Y_CENTER = self.data[trial]["Z_imag"][CENTER_INDEX]
                
                CURRENT_LOWER_BOUND_X, CURRENT_UPPER_BOUND_X = plt.xlim()
                DOMAIN = abs(CURRENT_UPPER_BOUND_X - CURRENT_LOWER_BOUND_X)
                NEW_DOMAIN = DOMAIN / ZOOM_RATIO
                
                CURRENT_LOWER_BOUND_Y, CURRENT_UPPER_BOUND_Y = plt.ylim()
                RANGE = abs(CURRENT_UPPER_BOUND_Y - CURRENT_LOWER_BOUND_Y)
                NEW_RANGE = RANGE / ZOOM_RATIO
                
                if abs(NEW_DOMAIN) > abs(NEW_RANGE):
                    NEW_RANGE = NEW_DOMAIN
                elif abs(NEW_RANGE) > abs(NEW_DOMAIN):
                    NEW_DOMAIN = NEW_RANGE
                
                NEW_LOWER_BOUND_X = X_CENTER - (NEW_DOMAIN / 2)
                NEW_UPPER_BOUND_X = X_CENTER + (NEW_DOMAIN / 2)
                NEW_LOWER_BOUND_Y = Y_CENTER - (NEW_RANGE / 2)
                NEW_UPPER_BOUND_Y = Y_CENTER + (NEW_RANGE / 2)
                
                plt.xlim(NEW_LOWER_BOUND_X, NEW_UPPER_BOUND_X)
                plt.ylim(NEW_LOWER_BOUND_Y, NEW_UPPER_BOUND_Y)
                
            else:
                # Set the x and y axis lengths equal by increasing the limits of the shorter axis to equal the longer
                if abs(self.data[trial]["Z_real"]).max() > abs(self.data[trial]["Z_imag"]).max():
                    # set the y axis limits the same as the x axis limits
                    plt.ylim(plt.xlim())
                else:
                    # set the x axis limits the same as the y axis limits
                    plt.xlim(plt.ylim())
            
            # Set a square aspect ratio
            plt.gca().set_aspect('equal', 'box')
        
        if type_ == "Bode":
            
            if xlabel == None:
                xlabel = r"$log(freq/Hz)$"
            if ylabel == None:
                ylabel = r"$log(|Z|/\Omega)$"
            
            ylabel2 = r"$Phase(Z)/deg$"
            
            plt.xlabel(xlabel, fontsize=16)
            plt.ylabel(ylabel, fontsize=16)

            plt.scatter(np.log10(self.frequency), np.log10(abs(self.Z.real)), color='blue')
            
            ax.spines['left'].set_color('blue')
            ax.yaxis.label.set_color('blue')
            ax.tick_params(axis='y', colors='blue')

            plt.twinx(gca())
            plt.ylabel(ylabel2, fontsize=16)
            
            plt.scatter(np.log10(self.frequency), np.angle(self.Z, deg=True), color='red')
            
            gca().spines['right'].set_color('red')
            plt.yaxis.label.set_color('red')
            gca().tick_params(axis='y', colors='red')
        
        
        if legend:
            plt.legend()
    
    
    
    
    def tafel(self, trial=None, voltage_bounds=[None,None], E_eq=0.0):
        '''
        perform tafel analysis
        '''
        # Define some constants and other values
        F = 96485                   # Faraday's constant (C/mol)
        R = 8.3145                  # Gas Constant (8.3145 L atm /mol /K)
        E_1 = voltage_bounds[0]     # Lower voltage boundary for analysis
        E_2 = voltage_bounds[1]     # Upper voltage boundary for analysis
        
        # check to ensure that E_1 is the smaller E value; if not, switch them
        if not E_1 < E_2:
            smaller_potential = E_2
            larger_potential = E_1
            E_1 = smaller_potential
            E_2 = larger_potential
        
        # Check if trial data has already been iR-corrected, and correct it if not
        if "iR-Corrected Voltage" not in self.data[trial]:
            self.correct_voltage("iR", trial=[trial])
        
        # Define dataset for Tafel analysis
        tafel_data = self.data[trial]
        tafel_data = tafel_data[tafel_data['ox/red'] == 0]            # Filter for cycle's reductive sweep (EC-lab convention: ox.=1, red.=0)
        tafel_data = tafel_data[tafel_data['Cycle'] == 1.0]           # First cycle only
        tafel_data = tafel_data[tafel_data['Voltage'] >= E_1]         # Apply lower voltage bound
        tafel_data = tafel_data[tafel_data['Voltage'] <= E_2]         # Apply upper voltage bound
        tafel_data = tafel_data[tafel_data['Current'] != 0]           # Filter for non-zero current only
        
        overpotential = self.data[trial]["Tafel overpotential"] = tafel_data['Corrected Voltage'] - E_eq
        self.update_history(f"Tafel overpotential calculated and saved for trial {trial}.")
        log_i = self.data[trial]["Tafel current"] = np.log10(abs(tafel_data['Current']))
        self.update_history(f"Tafel current calculated and saved for trial {trial}.")
        
        # Define highest index x value on which to perform Tafel regression (this fudge factor will ideally be eliminated in the future)
        index_limit = round(len(overpotential)/4)
        
        # Extract regression parameters
        tafel_slope, tafel_yint = np.polyfit(overpotential[:index_limit], log_i[:index_limit], 1)
        
        # Save extracted values
        # Tafel slope (V/dec)
        self.tafel_slope_V = tafel_slope
        self.update_history(f"Tafel slope (V) calculated and saved from Tafel analysis of trial {trial}.")
        # Tafel y-intercept
        self.tafel_yint = tafel_yint
        self.update_history(f"Tafel y-intercept calculated and saved from Tafel analysis of trial {trial}.")
        # Tafel slope (mv/dec)
        self.tafel_slope_mV = tafel_slope * 1000
        self.update_history(f"Tafel slope (mV) calculated and saved from Tafel analysis of trial {trial}.")
        # exchange current (mA)
        self.i_0 = 10**tafel_yint * 1000
        self.update_history(f"Exchange current (mA) calculated and saved from Tafel analysis of trial {trial}.")
    
    
    
    
    def linear_function(self, x, slope=None, y_int=None):
        y = slope*x + y_int
        return y
    
    
    
    
    def Correct_E_for_Ru(E, i, Ru):
        E_corrected = E - i*Ru
        return E_corrected
    
    
    
    
    # This is a duplicate with linear_function above, which is newer; need to compare and downselect
    def Regression(overpotential, log_i, step=False):
        if step == False:
            SLOPE, INTERCEPT, RVALUE, PVALUE, STDERR = stats.linregress(overpotential, log_i)
        else:
            SLOPE, INTERCEPT, RVALUE, PVALUE, STDERR = ReverseStepRegression(overpotential, log_i)
        return SLOPE, INTERCEPT
    
    
    
    
    # =============================================================================
    #   Extract Exchange Current and Exchange Rate Constant (from original pre-chemax Tafel package)
    # =============================================================================
    '''
    i_0 = 10**Tafel_yint
    OR
    i_0 = Tafel.ExchangeCurrent(E_corrected, I, E_eq=E_eq)
    
    k_0 = Tafel.ExchangeRateConstant(E_corrected, i, alpha=0.5, A=1.0, C_O=1.0, C_R=1.0, E_eq=E_eq, T=298.15)
    '''
    
    def ExchangeCurrent(E_corrected, i, E_eq=0.0):
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
        b, a = Regression(overpotential, log_i)
        # Remember, np.polyfit returns highest power coefficients first
        i_0 = 10**a
        return i_0


    def ExchangeRateConstant(E_corrected, i, alpha=0.5, A=1.0, C_O=1.0, C_R=1.0, E_eq=0.0, T=298.15):
        '''
        Returns the exchange rate constant (cm/s). 
        E_corrected represents the potential corrected for uncompensated series resistance, R_u.
        If the current argument is not a current density, then an area, A, must be provided to return an accurate exchange rate constant.
        If the current argument is a current density, no area needs to be provided. 

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
    
    
    
    # =============================================================================
    # Auto-regression functions
    # =============================================================================

    def IsLinear(R_squared, alpha):
        #provides criteria for defining the linear region
        #tests R^2 value
        if abs(R_squared) < alpha:
            return False
        else:
            return True


    def StepRegression(x, y, certainty=0.98):
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


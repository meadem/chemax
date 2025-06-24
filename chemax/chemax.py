# Standard library imports
from datetime import datetime
import numbers
import time
import os

# PyPI dependencies
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Local imports
from chemax.reference_electrodes import REFERENCE_ELECTRODES


class Experiment():
    '''
    Class object for importing, processing, and plotting echem data.
    '''
    
    def __init__(self, 
                 name=None, 
                 description=None, 
                 RE=None, 
                 display_RE=None, 
                 area=None, 
                 Ru=None
                 ):
        '''
        name = experiment name (could be same variable used for class instance.
                                Mostly used for error and print msgs.)
        description = experiment description
        RE = Reference Electrode used to acquire the data
        display_RE = Reference Electrode used to correct the voltage for display
        area = working electrode area (cm^2)
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
        

        try:
            # Set the reference electrode potential, for both the RE used to 
            # acquire the data and the one used to display it
            if RE == None:
                print(f'Reference electrode not specified for {self.name}.')
            elif RE in REFERENCE_ELECTRODES:
                for entry in REFERENCE_ELECTRODES:
                    if RE == entry:
                        self.REFERENCE_ELECTRODE_POTENTIAL = REFERENCE_ELECTRODES[entry]
            else:
                print(f'Specified reference electrode {self.RE} not found for {self.name}.'
                    'Unexpected behavior may occur when using an unrecognized reference electrode.')
        except:
            raise Exception(f"Error while setting the reference electrode as {self.RE} during initialization for {self.name}.")
        
        
        try:
            # Set the display RE as the data acquisition RE, if display RE has 
            # not been specified
            if display_RE is None:
                self.display_RE = RE
            else:
                self.display_RE = display_RE
        except:
            raise Exception(f"Error while setting the display reference electrode as {self.display_RE} during initialization for {self.name}.")
        

        # Print message if area is not specified.
        if area == None:
            print(f'Working electrode area not specified for {self.name}.')
    
    
    def calculate_voltage_decay(self, trial=None, interval=[None,None]):
        '''
        

        Parameters
        ----------
        trial : TYPE, optional
            DESCRIPTION. The default is None.
        interval : TYPE, optional
            Specified in hours. The default is [None,None].

        Returns
        -------
        result : TYPE
            DESCRIPTION.

        '''
        TIME = self.data[trial]["Time"]
        VOLTAGE = self.data[trial]["Voltage"]
        START_TIME = interval[0]
        FINAL_TIME = interval[1]
        TIME_SECONDS_FILTERED = TIME[START_TIME*3600:FINAL_TIME*3600]
        TIME_HOURS_FILTERED = TIME_SECONDS_FILTERED/3600
        VOLTAGE_V_FILTERED = VOLTAGE[TIME == TIME_SECONDS_FILTERED]
        VOLTAGE_mV_FILTERED = VOLTAGE_V_FILTERED*1000
        START_VOLTAGE = VOLTAGE_V_FILTERED[0]
        FINAL_VOLTAGE = VOLTAGE_V_FILTERED[len(VOLTAGE_V_FILTERED)-1]
        
        coef = np.polyfit(TIME_HOURS_FILTERED, VOLTAGE_mV_FILTERED, 1)
        DECAY_RATE = coef[0]     # mV/hour
        
        _fig = plt.figure(0, figsize=(15,10))
        self.plot("CP", trials=[trial], x_zoom=interval, legend=False)
        plt.gca().plot(TIME_HOURS_FILTERED, np.poly1d(coef))
        
        result = print(f"Voltage decay rate from {START_TIME} hours ({START_VOLTAGE:.4} V) to {FINAL_TIME} hours ({FINAL_VOLTAGE:.4} V) = {DECAY_RATE:.4} mV/hr")
        
        return result
    
    
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
        Adds a new column to the specified trials' dataframe(s) containing the 
        voltage corrected for either:
        a) uncompensated series resistance (new column referenced as 
                                            self.data[trial]["iR-Corrected Voltage"])
        b) display reference electrode (new column referenced as 
                                        self.data[trial]["RE-Corrected Voltage"])
        
        Adds another new column referenced as 
        self.data[trial]["Corrected Voltage"]
        This column is one of the following:
        1) the iR-Corrected Voltage (if this is the only correction performed 
                                     on the voltage so far)
        2) the RE-Corrected Voltage (if this is the only correction performed 
                                     on the voltage so far)
        3) Voltage corrected for both iR and display RE (if both corrections 
                                                         have been performed)
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
                # Check first if RE used to acquire data is the one being 
                # requested for display
                if RE == self.RE:
                    if not silent:
                        print(f"{RE} is already the reference electrode implicitly referenced by the acquired voltage data.")
                else:
                    # Fetch the RE potential from its name using the 
                    # dictionary, or print an error
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
                                # BUG: this message prints even for trials 
                                # where no correction was performed (need to fix)
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
    
    
    def cottrell(self, trials=[], n=1, C=1):
        '''
        Perform cottrell analysis. (current values converted to mA for analysis)
        '''

        F = 96485    # Faraday constant (C/mol)
        
        try:
            # if a list of trials is not passed, default is to perform 
            # analysis on all CA data
            if not trials:
                for trial in self.metadata:
                    if self.metadata[trial]["technique"] == "CA":
                        trials.append(trial)
        except:
            raise Exception(f"Error while selecting trials on which to perform Cottrell analysis for {self.name}.")
        
        for trial in trials:
            try:
                # Check if Corrected Time for trial data has already been 
                # calculated; do so if not (is this check even necessary?)
                # NOTE: it may be more accurate to correct the start time w/o 
                # zeroing, similar to how I did this originally 
                # (can do easily??)
                if "Corrected Time" not in self.data[trial]:
                    self.correct_time(trials=[trial])
            except:
                raise Exception(f"Error while zero-correcting time data during Cottrell analysis for {self.name} trial {trial}.")

            try:
                # Check if t^-0.5 for trial data has already been calculated; 
                # do so if not (is this check even necessary?)
                if "Time-Inverse-Root" not in self.data[trial]:
                    self.data[trial]["Time-Inverse-Root"] = 1 / np.sqrt(self.data[trial]["Corrected Time"])
                    self.update_history(f"Time-Inverse-Root calculated from Corrected Time for trial {trial} in {self.name}.")
            except:
                raise Exception(f"Error while calculating time-inverse-root data during Cottrell analysis for {self.name} trial {trial}.")

            try:
                # filter out data where current <= peak current / e
                TIME_INVERSE_ROOT = self.data[trial]["Time-Inverse-Root"]
                CURRENT = self.data[trial]["Current"]*1000
                FILTERED_CURRENT = CURRENT[CURRENT <= (CURRENT.max() / np.e)]
                COTTRELL_INDEX_MIN = FILTERED_CURRENT.idxmax()
                COTTRELL_DOMAIN = TIME_INVERSE_ROOT[COTTRELL_INDEX_MIN:]
                COTTRELL_RANGE = CURRENT[COTTRELL_INDEX_MIN:]
            except:
                raise Exception(f"Error while filtering/pre-processing data during Cottrell analysis for {self.name} trial {trial}.")

            try:
                # linear regression
                COTTRELL_SLOPE, COTTRELL_INTERCEPT = np.polyfit(COTTRELL_DOMAIN, 
                                                                COTTRELL_RANGE, 
                                                                1)
            except:
                raise Exception(f"Error performing linear regression during Cottrell analysis for {self.name} trial {trial}.")

            try:
                # save as attributes of trial's dataframe object
                self.data[trial].COTTRELL_INDEX_MIN = COTTRELL_INDEX_MIN
                self.data[trial].COTTRELL_SLOPE = COTTRELL_SLOPE
                self.data[trial].COTTRELL_INTERCEPT = COTTRELL_INTERCEPT
            except:
                raise Exception(f"Error while saving calculated values as attributes during Cottrell analysis for {self.name} trial {trial}.")
            
            try:
                # (cm^2 /s) (area expected in cm^2)
                # Bulk conc of species for which the diffusion coefficient is 
                # being calculated has units of (mol/cm^3, i.e., mol/mL, or MM)
                DIFFUSION_COEFFICIENT = np.pi * ( COTTRELL_SLOPE / ( n*F*self.area*C*1000 ) )**2
                print(f"A diffusion coefficient of {DIFFUSION_COEFFICIENT} was extracted using Cottrell analysis from {self.name} trial {trial}")
            except:
                raise Exception(f"Error while calculating/returning the diffusion coefficient during Cottrell analysis for {self.name} trial {trial}.")
    

    def extract_Ru_from_HFR(self, trial=None, silent=False):
        '''
        Extract uncompensated series resistance, R_u, from EIS data by 
        extracting the minimum real impedance value.
        '''
        try:
            # Extract Ru value (and its index) from the real impedance data
            result = self.data[trial]["Z_real"].min()
            result_index = np.argmin(self.data[trial]["Z_real"])
            
            # Nyquist plot zoomed in around the extracted value
            print("-------------------------------------------------------------------------------------")
            plt.figure(figsize=([10,4]))
            self.plot("Nyquist", position=121, trials=[trial], legend=False)
            self.plot("Nyquist", position=122, trials=[trial], legend=False, _zoom=[15, result_index])
            plt.gcf().supxlabel(f"The Ru value extracted from the impedance data for {self.name} trial {trial} is {result} ohms.", 
                                fontsize=12, y=-0.1)
            plt.show(block=False)
            
            # Check extracted value with user
            time.sleep(0.001)
            check = input("Accept this value? (Y/N)")
            time.sleep(0.001)
            print("-------------------------------------------------------------------------------------")
            if check.upper() == "N":
                print("REJECTED")
                self.update_history(f"Failed to extract uncompensated series resistance from {self.name} trial {trial} impedance data.")
                return
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
        
        # query user if they want to apply extracted Ru value to iR-correct 
        # the voltage data, then do so if requested
        try:
            unresolved = True
            while unresolved:
                time.sleep(0.001)
                check = input("iR-correct the voltage data? "
                              "It will be saved into a new column 'iR-Corrected Voltage' without overwriting 'Voltage'. (Y/N)")
                time.sleep(0.001)
                print("-------------------------------------------------------------------------------------")
                if check.upper() == "N":
                    unresolved = False
                elif check.upper() == "Y":
                    self.correct_voltage("iR", silent=silent)
                    unresolved = False
                else:
                    print("Please respond with 'Y' to iR-correct the voltage data, or 'N' to not correct it.")
                    print("-------------------------------------------------------------------------------------")
        except:
            raise Exception("Error while iR-correcting voltage data for {self.name}.")
    

    def extract_E0_from_CV(self, trial=None, silent=False):
        '''
        Extract the standard reduction potential, E0, from the CV's anodic and cathodic peak current.
        '''
        x = self.data[trial]['Voltage']
        x = self.data[trial]['Current']
        Eo = x[y == y.max()].max()
        Er = x[y == y.min()].max()
        E0 = (Eo + Er) / 2
        self.data[trial].E0 = E0

    
    def extract_IV_values(self):
        '''
        Extract voltage and power from last minute of each CP trial. 
        
        (Current is already extracted and saved during load)
        
        TO DO: 
            - Need conditional for getting to power v. power density
            - Save Power to self.data
        '''
        for trial in self.data:
            if self.metadata[trial]['technique'].upper() == "CP":
                # Time
                TIME = self.data[trial]["Time"]
                EXPERIMENT_DURATION = TIME.max()                                                    # (in seconds)
                # Voltage
                VOLTAGE_SAMPLE = self.data[trial]["Voltage"][TIME > (EXPERIMENT_DURATION - 60)]     # Samples voltage from final 60 seconds
                VOLTAGE_MEAN = round(VOLTAGE_SAMPLE.mean(), 3)                                      # Averages sampled voltage
                self.metadata[trial]['I-V Voltage (V)'] = VOLTAGE_MEAN
                # Current & Power
                CURRENT = self.metadata[trial]['I-V Current (mA)']
                POWER = CURRENT * VOLTAGE_MEAN / 1000
                self.metadata[trial]['I-V Power (W)'] = POWER
                if self.area != None:
                    CURRENT_DENSITY = self.metadata[trial]['I-V Current Density (mA/cm2)']
                    POWER_DENSITY = CURRENT_DENSITY * VOLTAGE_MEAN / 1000
                    self.metadata[trial]['I-V Power Density (W/cm2)'] = POWER_DENSITY
    
    
    def find_header(self, file_path, search_string):
        '''
        Find line number for header during file load operations.

        Parameters
        ----------
        file_path : STR
            relative path to file, including file name and suffix
        search_string : STR
            Bit of text to locate in file

        Returns
        -------
        Integer result representing line number for header in specified file

        '''
        
        result = None
        with open(file_path, 'r') as working_file:
            for line_number, line in enumerate(working_file, start=1):
                if search_string in line:
                    result = line_number
        return result
    
    
    def linear_function(self, x, slope=None, y_int=None):
        y = slope*x + y_int
        return y
    

    def load(self, num=1, step=1, filetype='.txt', file='data', folder='data', 
             file_numbers=None, technique=None, potentiostat=None, 
             silent=False, tag=None):
        '''
        Upload data files into an Experiment object as a pandas dataframe.
        
        Built for use with sequentially-numbered data files with a shared base
        file name, e.g., 'CV_thursday_1', 'CV_thursday_2', etc.

        Parameters
        ----------
        num : INT, optional
            DESCRIPTION. The default is 1.
        step : INT, optional
            DESCRIPTION. The default is 1.
        filetype : STR, optional
            DESCRIPTION. The default is '.txt'.
        file : STR, optional
            DESCRIPTION. The default is 'data'.
        folder : STR, optional
            DESCRIPTION. The default is 'data'.
        file_numbers : LIST, optional
            A list that can be passed to indicated specific file numbers. 
            The default is None.
        technique : STR, optional
            Used to label the loaded file(s) with metadata label indicating if
            it is data for a specific electrochemical technique, e.g., EIS, CV,
            etc. Useful for bulk plotting later using the plot() method. The 
            technique can be any arbitrary label.
            The default is None.
        potentiostat : STR, optional
            Tailors the data loading method for the potentiostat producing the
            data file. Current options are 'Biologic' and 'Gamry'. 
            The default is None.
        silent : BOOL, optional
            Will print status messages and experiment/file mapping list if
            False. 
            The default is False.
        tag : STR, optional
            Used to tag the loaded file(s) with a metadata label that can be 
            used later-- e.g., during plotting-- to process the group together.
            The default is None.

        Returns
        -------
        None.

        '''
        
        # tag if this operation is an instantiation (alternative is an 
        # append operation)
        _IS_INSTANTIATION = not hasattr(self, "data")
        
        try:
            # create temporary dicts that will be incorporated into self.data 
            # and self.metadata at the end of the method call
            data = {}
            metadata = {}

            # create a temporary dict that holds key:value pairs of 
            # INTERNAL:EXTERNAL data indices during import
            # used to force consecutive integer values for internal data reference
            # e.g., files "data_1" and "data_3" would be indexed as 1 and 2, 
            # as would files "CV_data_1" and "RDE_data_1"
            # ----------------------------------------------------------------------------------------------------------
            EXTERNAL_DATA_INDICES = None
            INTERNAL_DATA_INDICES = None
            
            # EXTERNAL DATA INDICES
            # if provided, use file_numbers argument in place of the num 
            # argument for defining the external data indices
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
            
            # iteratively load file data
            for INTERNAL_DATA_INDEX in DATA_INDICES:
                EXTERNAL_DATA_INDEX = DATA_INDICES[INTERNAL_DATA_INDEX]
                FILE_NAME = file + str(EXTERNAL_DATA_INDEX) + filetype
                FILE_PATH = os.path.join(folder, FILE_NAME)
                
                # Read file into memory; check if potentiostat info is given 
                # and modify call to pd.read_csv accordingly
                if potentiostat.upper() == "GAMRY":
                    HEADER = self.find_header(FILE_PATH, "TABLE")
                    data[INTERNAL_DATA_INDEX] = pd.read_csv(FILE_PATH, header=HEADER, sep='\t', engine='python', encoding='unicode_escape', skiprows=[HEADER+1])
                
                elif potentiostat.upper() == "BIOLOGIC":
                    data[INTERNAL_DATA_INDEX] = pd.read_csv(FILE_PATH, sep=None, engine='python', encoding='unicode_escape')
                
                else:
                    data[INTERNAL_DATA_INDEX] = pd.read_csv(FILE_PATH, sep=None, engine='python', encoding='unicode_escape')
                
                metadata[INTERNAL_DATA_INDEX] = dict(source=FILE_PATH, technique=technique, tag=tag)
                
                if not silent:
                    print(f"Data file {FILE_PATH} successfully imported (INDEX={INTERNAL_DATA_INDEX}) for {self.name}.")
        
        except:
            raise Exception(f"Error while importing data files from {folder} into {self.name}.")
        
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
            GEIS_DC_CURRENT_HEADER = "GEIS DC Current (A)"

            # A dictionary of headers I've observed previously and what to 
            # replace them with when found.
            HEADER_DICT = {"VOLTAGE":VOLTAGE_HEADER,
                           "EWE/V":VOLTAGE_HEADER,
                           "<EWE>/V":VOLTAGE_HEADER,
                           "E (V)":VOLTAGE_HEADER,
                           "EWE/MV":VOLTAGE_HEADER,
                           "E (MV)":VOLTAGE_HEADER,
                           "VF":VOLTAGE_HEADER,
                           "CURRENT":CURRENT_HEADER,
                           "<I>/A":CURRENT_HEADER,
                           "<I>/MA":CURRENT_HEADER,
                           "I/MA":CURRENT_HEADER,
                           "I/A":CURRENT_HEADER,
                           "IM":CURRENT_HEADER,
                           "TIME (S)":TIME_HEADER,
                           "TIME":TIME_HEADER,
                           "TIME/S":TIME_HEADER,
                           "T":TIME_HEADER,
                           "CYCLE":CYCLE_HEADER,
                           "CYCLE NUMBER":CYCLE_HEADER,
                           "POWER":POWER_HEADER,
                           "POWER (W)":POWER_HEADER,
                           "P (W)":POWER_HEADER,
                           "P/W":POWER_HEADER,
                           "FREQ/HZ":FREQUENCY_HEADER,
                           "FREQ":FREQUENCY_HEADER,
                           "RE(Z)/OHM":ZREAL_HEADER,
                           "ZREAL":ZREAL_HEADER,
                           "-IM(Z)/OHM":ZIMAG_HEADER,
                           "ZIMAG":ZIMAG_HEADER,
                           "|Z|/OHM":Z_HEADER,
                           "ZMOD":Z_HEADER,
                           "PHASE(Z)/DEG":PHASE_ANGLE_HEADER,
                           "ZPHZ":PHASE_ANGLE_HEADER,
                           "IDC":GEIS_DC_CURRENT_HEADER
                          }
            
            # Apply header standardization 
            # (also standardize electrical current data as amperes)
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
        
        # Standaridize units and assign more metadata 
        # (This is easier to do after header standardization, above)
        try:
            for trial in data:
                
                # CA
                # calculate applied voltage from data for CA experiments; set as default label
                # written for Biologic potentiostat data; may need modification for other files
                if metadata[trial]['technique'].upper() == "CA":
                    VOLTAGE = data[trial]["Voltage"]*1000
                    applied_voltage = str(int(VOLTAGE.mean())) + " mV"
                    metadata[trial]['applied voltage'] = applied_voltage
                    metadata[trial]['default label'] = applied_voltage
                
                # CV, CV SIM, LSC
                # calculate scan rate from data for CV, etc., experiments; set as default label
                # written for Biologic potentiostat data; may need modification for other files
                elif metadata[trial]['technique'].upper() in ["CV", "CV SIM", "LSV"]:
                    TIME = data[trial]["Time"]
                    VOLTAGE = data[trial]["Voltage"]*1000
                    VOLTAGE = VOLTAGE[data[trial]["Cycle"] == 1]
                    if "ox/red" not in data[trial]:
                        VOLTAGE = VOLTAGE[VOLTAGE.idxmin():VOLTAGE.idxmin()+100]
                    else:
                        VOLTAGE = VOLTAGE[data[trial]["ox/red"] == 1]
                    INDEX1 = VOLTAGE.idxmax()
                    INDEX2 = VOLTAGE.idxmin()
                    VOLTAGE_RANGE = abs(VOLTAGE[INDEX1] - VOLTAGE[INDEX2])
                    TIME_RANGE = abs(TIME[INDEX1] - TIME[INDEX2])
                    SCAN_RATE = str(int(VOLTAGE_RANGE / TIME_RANGE)) + " mV/s"
                    metadata[trial]['scan rate'] = SCAN_RATE
                    metadata[trial]['default label'] = SCAN_RATE
                
                # CP
                # Determine current and current density (if electrode area provided) for CP experiments; save for I-V curve analysis and set as default label
                # Written for Gamry potentiostat data; may need modification for other files
                elif metadata[trial]['technique'].upper() == "CP":
                    TIME = data[trial]["Time"]
                    CURRENT = data[trial]["Current"][TIME > 0].mean()*1000      # average current (mA)
                    CURRENT = int(round(CURRENT))
                    CURRENT_LABEL = str(CURRENT) + " mA"
                    metadata[trial]['I-V Current (mA)'] = CURRENT
                    metadata[trial]['applied current'] = CURRENT_LABEL
                    
                    # Use current density when electrode area is provided
                    if self.area != None:
                        AREA = self.area
                        CURRENT_DENSITY = int(round(CURRENT / AREA))
                        CURRENT_DENSITY_LABEL = str(CURRENT_DENSITY) + r"$\, mA \, {cm}^{-2}$"
                        metadata[trial]['I-V Current Density (mA/cm2)'] = CURRENT_DENSITY 
                        metadata[trial]['applied current density'] = metadata[trial]['default label'] = CURRENT_DENSITY_LABEL
                    
                    # Use current when electrode area is not provided
                    else:
                        metadata[trial]['default label'] = CURRENT_LABEL
                
                # GEIS
                # Determine DC current used for GEIS experiment and save as default label
                # Written for Gamry potentiostat data; may need modification for other files
                elif metadata[trial]['technique'].upper() == "GEIS":
                    GEIS_DC_CURRENT = data[trial]["GEIS DC Current (A)"] * 1000 # convert to mA
                    CURRENT_MEAN = GEIS_DC_CURRENT.mean()
                    CURRENT_LABEL = None
                    if self.area != None:
                        CURRENT_DENSITY_MEAN = CURRENT_MEAN / self.area
                        CURRENT_LABEL = str(int(round(CURRENT_DENSITY_MEAN))) + r"$\, mA \, {cm}^{-2}$"
                    else:
                        CURRENT_LABEL = str(int(round(CURRENT_MEAN))) + " mA"
                    metadata[trial]['default label'] = CURRENT_LABEL
                    
                
                else:
                    metadata[trial]['applied voltage'] = None
                    metadata[trial]['scan rate'] = None
                    metadata[trial]['default label'] = None
                    metadata[trial]['applied current density'] = None
                    metadata[trial]['applied current'] = None
                    
        except:
            raise Exception(f"Error while assigning metadata during loading for {self.name} trial {trial}.")
        
        
        try:
            # save data and metadata in namespace
            if _IS_INSTANTIATION:
                self.data = data
                self.metadata = metadata
            else:
                self.data.update(data)
                self.metadata.update(metadata)
        except:
            raise Exception(f"Error while saving data and/or metadata imported from {folder} for {self.name}.")
        
        
        try:
            # If a display reference electrode has been specified previously, prompt whether it should be applied during loading
            if self.display_RE is not None:
                time.sleep(0.001)
                check = input(f"Correct voltages for {self.name} using display RE '{self.display_RE}'? "
                              f"Values will be saved in a new column. (Y/N)")
                time.sleep(0.001)
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
            print(f"Data successfully imported from '{folder}' for {self.name}.")
    
    
    def plot(self, type_, technique=None, title="", label=None, trials=None, 
             position=111, trim=[None,None], xlabel=None, ylabel=None,
             legend=True, voltage_units="V", current_units="A", density=True, 
             series=False, group=False, fontsize=24, tag=None, 
             x_zoom=[None,None], y_zoom=[None,None], _zoom=[None,None]):
        '''
        Plot the specified curve type in the pyplot figure that currently has 
        focus. Method is designed to make it easy to plot many well-formatted 
        trials together with few commands. Labels and formatting are 
        pre-specified for each plot type.
        
        Parameters
        ----------
        type_ : STR
            Type of plot. Currently available type options are:
                - "CV"
                    - Cyclic Voltammetry (should also work for LSV)
                    - X-axis: Voltage        Y-axis: Current
                - "CV-corrected"
                    - CV corrected for time and/or iR
                    - X-axis: Voltage        Y-axis: Current
                - "CP"
                    - Chronopotentiometry, aka galvanostatic
                    - X-axis: Time           Y-axis: Voltage
                - "CA"
                    - Chronoamperometry, aka potential-step
                    - X-axis: Time           Y-axis: Current
                - "IV"
                    - I-V curve, aka Polarization curve
                    - X-axis: Current        Y-axis: Voltage
                                             Y2-axis: Wattage
                - "Tafel"
                    - Tafel analysis of CV or LSV data
                    - X-axis: Overpotential  Y-axis: log(Current)
                - "Cottrell"
                    - Cottrell analysis of CA data
                    - X-axis: time^(-1/2)    Y-axis: Current
                - "Nyquist"
                    - Nyquist analysis of EIS data
                    - X-axis: Z_real         Y-axis: -Z_imag
                - "Bode"
                    - Bode analysis of EIS data
                    - X-axis: log(Hz)        Y-axis: log(Z)
                                             Y2-axis: Phase Angle
        technique : STR, optional
            Can be used in place of trials param to fetch and plot several 
            traces with one call. (Technique param must be used when loading 
                                   data using load().)
            The default is None.
        title : STR, optional
            Text used for title of figure/plot.
            The default is "".
        label : STR, optional
            Text used for legend label. Trial labels are usually set by default
            from information extracted during load(), but can be overridden 
            here, if desired.
            The default is None.
        trials : LIST, optional
            List of integers. Can be used to fetch and plot several traces with
            one call. Refer to trials' internal indices in the self.data object
            (also displayed during load() operation if silent=False.)
            The default is None.
        position : INT, optional
            Three digit integer used to specify plotting position within the 
            pyplot figure. Used the same way as the three-digit integers passed
            to plt.subplot() in the example at the link:
                https://matplotlib.org/stable/tutorials/pyplot.html#logarithmic-and-other-nonlinear-axes
            The default is 111.
        trim : TUPLE, optional
            Can be used to slice the trial data using known indices. 
            EXAMPLE: data[ trim[0] : trim[1] ]
            The default is [None,None].
        xlabel : STR, optional
            Default x-axis labels used for each type, but can be overridden 
            using xlabel and ylabel params. 
            The default is None. (defaults to the pre-specified default label)
        ylabel : STR, optional
            Default x-axis labels used for each type, but can be overridden 
            using xlabel and ylabel params. 
            The default is None. (defaults to the pre-specified default label)
        legend : BOOL, optional
            If True, the figure legend is displayed.
            The default is True.
        voltage_units : STR, optional
            Used to specify voltage is plotted in units of "V" (volts) or "mV" 
            (millivolts).
            The default is "V".
        current_units : STR, optional
            Used to specify current is plotted in units of "A" (amperes) or 
            "mA" (milliamperes).
            The default is "A".
        density : BOOL, optional
            If True, current and power data will be displayed as a 'current 
            density' or 'power density' if possible. Will fail if electrode 
            area is not provided.
            The default is True.
        series : BOOL, optional
            If True, plot types with Time as the independent variable will be 
            plotted in chronological sequence instead of graphically 'stacked.'
            The default is False.
        group : BOOL, optional
            If True, all trials passed to plot() will be concatenated together 
            as one data array, and can be labeled and presented as though a 
            single trace. Useful for bundling and labelling sequential tests, 
            e.g., constant-current step-up for conditioning.
            The default is False.
        fontsize : INT, optional
            fontsize is used to specify the axis labels font size; fontsize-6 
            is used for the axis ticks and legend labels; fontsize+6 is used 
            for the figure title.
            The default is 24.
        tag : STR, optional
            Can be used to specify a 'tag' that was used to label the data
            during loading with the load() method. If used, all data with the
            specified tag will be plotted, similar to using 'technique.'
            The default is None.
        x_zoom : TUPLE, optional
            Can be used to zoom the figure using specified x-axis value 
            (not index) boundaries (x_zoom[0]: x_min; x_zoom[1]: x_max).
            The default is [None,None].
        y_zoom : TUPLE, optional
            Can be used to zoom the figure using specified y-axis value 
            (not index) boundaries (y_zoom[0]: y_min; y_zoom[1]: y_max).
            The default is [None,None].
        _zoom : TUPLE, optional
            Can be used to re-scale the figure (i.e., zoom). Usually only used 
            by self.extract_Ru_from_HFR() method. Note that the last plotted 
            trial has focus when zooming.
            _zoom[0]: zoom ratio; passing a value of, e.g., 10, will zoom in 
                      10x around the specified index
            _zoom[1]: zoom center index
            The default is [None,None].

        Returns
        -------
        None.

        '''
        
        # setup/pre-processing for plotting
        try:
            
            # Label voltage as being relative to the reference electrode (implement later)
            display_RE = ''
            if self.display_RE is None:
                display_RE = '???'
            else:
                display_RE = self.display_RE
            
            # apply specified trim, if applicable
            FILTER_START = trim[0]
            FILTER_STOP = trim[1]
            
            # retrieve current figure and prep for plotting
            plt.gcf()
            plt.subplot(position)
            plt.title(title, fontsize=fontsize+6)
            plt.axhline(color='k', linewidth=0.25)
            plt.axvline(color='k', linewidth=0.25)
            plt.xticks(fontsize=fontsize-6)
            plt.yticks(fontsize=fontsize-6)
            
            # Determine trials to plot, using params in this priority order:
            # 1) 'trials' is passed and this entire loop is skipped
            # 2) 'technique'
            # 3) 'tag' (can be used in combo with technique to sub-filter)
            # 4) if nothing is specified, plot all data in experiment object
            if trials is None:
                
                # Plot all data in experiment object
                if technique is None:
                    if tag is None:
                        trials = self.data
                    
                # Use Tag filter only
                    else:
                        trials = []
                        for trial in self.metadata:
                            if self.metadata[trial]["tag"] == tag:
                                trials.append(trial)
                
                # Use Technique + Tag filter in combo
                elif technique != None and tag != None:
                    trials = []
                    for trial in self.metadata:
                        if self.metadata[trial]["tag"] == tag and self.metadata[trial]["technique"] == technique:
                            trials.append(trial)
                
                # Use Technique filter only
                else:
                    trials = []
                    for trial in self.metadata:
                        if self.metadata[trial]["technique"] == technique:
                            trials.append(trial)
        
        except:
            raise Exception(f"Error during setup/pre-processing while attempting to plot a {type_} figure for {self.name}.")
        
        if type_ == "CV":
            
            if xlabel == None:
                xlabel = f"Potential (V v. {display_RE})"
            if ylabel == None:
                ylabel = "Current (A)"
            
            plt.xlabel(xlabel, fontsize=fontsize)
            plt.ylabel(ylabel, fontsize=fontsize)

            for trial in trials:
                try:
                    # Use RE-Corrected Voltage if available
                    VOLTAGE = None
                    CURRENT = self.data[trial]["Current"]
                    if "RE-Corrected Voltage" in self.data[trial].columns:
                        VOLTAGE = self.data[trial]["RE-Corrected Voltage"]
                    else:
                        VOLTAGE = self.data[trial]["Voltage"]
                    
                    # if not label is provided, use default label, if available
                    _label = label
                    if label is None:
                        _label = self.metadata[trial]["default label"]
                    
                    # Plot using appropriate units and labeling
                    if voltage_units.upper() == "V" and current_units.upper() == "A":
                        plt.plot(VOLTAGE[FILTER_START:FILTER_STOP], 
                                CURRENT[FILTER_START:FILTER_STOP], 
                                label=_label)
                    elif voltage_units.upper() == "V" and current_units.upper() == "MA":
                        plt.ylabel("Current (mA)")
                        plt.plot(VOLTAGE[FILTER_START:FILTER_STOP], 
                                CURRENT[FILTER_START:FILTER_STOP]*1000, 
                                label=_label)
                    elif voltage_units.upper() == "MV" and current_units.upper() == "A":
                        plt.xlabel("Voltage (mV)")
                        plt.plot(VOLTAGE[FILTER_START:FILTER_STOP]*1000, 
                                CURRENT[FILTER_START:FILTER_STOP], 
                                label=_label)
                    elif voltage_units.upper() == "MV" and current_units.upper() == "MA":
                        plt.xlabel("Voltage (mV)")
                        plt.ylabel("Current (mA)")
                        plt.plot(VOLTAGE[FILTER_START:FILTER_STOP]*1000, 
                                CURRENT[FILTER_START:FILTER_STOP]*1000, 
                                label=_label)
                    else:
                        raise NameError("Units not found or other error while plotting CV curve.")
                
                except:
                    raise Exception(f"Error while plotting {type_} curve for trial {trial} in {self.name}")
        
        if type_ == "CP":
            
            # If labels are not provided in method call, use default axis labels for chronopotentiometric data
            if xlabel == None:
                xlabel = "Time (hr)"
            if ylabel == None:
                ylabel = "Voltage (V)"
            plt.xlabel(xlabel, fontsize=fontsize)
            plt.ylabel(ylabel, fontsize=fontsize)
            
            # counter used for 'series' form CP plotting
            runtime_to_this_point = 0
            
            # Group data together as one set for plotting if group == True
            # Current method of copying arrays in place to force append is very inefficient and slow; need to revise later
            if group == True:
                group_array = {1:{"Voltage":np.array([]), "Time":np.array([])}}
                for trial in trials:
                    TIME = self.data[trial]["Time"]
                    DURATION = TIME.max()
                    UPDATED_TRIAL_TIME = [point + runtime_to_this_point for point in TIME]
                    group_array[1]["Time"] = np.append(group_array[1]["Time"], UPDATED_TRIAL_TIME)
                    VOLTAGE = self.data[trial]["Voltage"]
                    group_array[1]["Voltage"] = np.append(group_array[1]["Voltage"], VOLTAGE)
                    if series == True:
                        runtime_to_this_point += DURATION
                trials = group_array
            
            # Iterative treatment of desired trials
            for trial in trials:
                try:
                    
                    # Fetch data series to plot
                    TIME = None
                    VOLTAGE = None
                    DURATION = None
                    if group == False:
                        TIME = self.data[trial]["Time"] / 3600
                        DURATION = TIME.max()
                        VOLTAGE = self.data[trial]["Voltage"]
                    else:
                        TIME = trials[trial]["Time"] / 3600
                        VOLTAGE = trials[trial]["Voltage"]
                        runtime_to_this_point = 0
                    
                    # if no label is provided, use default label, if available
                    _label = label
                    if label is None:
                        _label = self.metadata[trial]["default label"]
                    
                    # Plot using appropriate units and labeling
                    TIME_PLOTTED = [point + runtime_to_this_point for point in TIME[FILTER_START:FILTER_STOP]]
                    if voltage_units.upper() == "V":
                        plt.plot(TIME_PLOTTED,
                                 VOLTAGE[FILTER_START:FILTER_STOP],
                                 label=_label)
                    elif voltage_units.upper() == "MV":
                        plt.xlabel("Voltage (mV)")
                        plt.plot(TIME_PLOTTED,
                                 VOLTAGE[FILTER_START:FILTER_STOP]*1000,
                                 label=_label)
                    
                    # Update the runtime variable used to array CP data in series instead of stacking
                    if series == True and group == False:
                        runtime_to_this_point += DURATION
                    
                    # Set x-axis time ticks to every 60 minutes
                    length_of_longest_trial = 0
                    for line in plt.gca().get_lines():
                        time_line = line.get_xdata()
                        if max(time_line) > length_of_longest_trial:
                            length_of_longest_trial = max(time_line)
                    
                    TIME_TICK_INTERVAL = None
                    
                    if length_of_longest_trial+1 > 100:
                        TIME_TICK_INTERVAL = 50.0
                    elif length_of_longest_trial+1 > 48:
                        TIME_TICK_INTERVAL = 24.0
                    elif length_of_longest_trial+1 > 12:
                        TIME_TICK_INTERVAL = 6.0
                    else:
                        TIME_TICK_INTERVAL = 1.0
                    
                    plt.xticks(np.arange(0, length_of_longest_trial+1, TIME_TICK_INTERVAL))
                
                except:
                    raise Exception(f"Error while plotting {type_} curve for trial {trial} in {self.name}")
        
        if type_ == "CV-corrected":
            
            if xlabel == None:
                xlabel = f"Corrected Potential (V v. {display_RE})"
            if ylabel == None:
                ylabel="Current (A)"
            
            plt.xlabel(xlabel, fontsize=fontsize)
            plt.ylabel(ylabel, fontsize=fontsize)
            
            for trial in trials:
                try:
                    # Check if iV-corrected voltage has already been calculated; do it if it hasn't been
                    if "iR-Corrected Voltage" not in self.data[trial]:
                        self.correct_voltage("iR", trials=[trial])
                    
                    # if not label is provided, use default label, if available
                    _label = label
                    if label is None:
                        _label = self.metadata[trial]["default label"]
                    
                    # Plot corrected voltage using the specified current and voltage units
                    # Note that voltage plots from "Corrected Voltage" and not "iR-Corrected Voltage"; 
                    # this incorporates the RE adjustment if one has been made, and otherwise is just the iR-corrected voltage anyways
                    if voltage_units.upper() == "V" and current_units.upper() == "A":
                        plt.plot(self.data[trial]["Corrected Voltage"][FILTER_START:FILTER_STOP], 
                                 self.data[trial]["Current"][FILTER_START:FILTER_STOP], 
                                 label=_label)
                    elif voltage_units.upper() == "V" and current_units.upper() == "MA":
                        plt.ylabel("Current (mA)", fontsize=16)
                        plt.plot(self.data[trial]["Corrected Voltage"][FILTER_START:FILTER_STOP], 
                                 self.data[trial]["Current"][FILTER_START:FILTER_STOP]*1000, 
                                 label=_label)
                    elif voltage_units.upper() == "MV" and current_units.upper() == "A":
                        plt.xlabel("Corrected Voltage (mV)", fontsize=16)
                        plt.plot(self.data[trial]["Corrected Voltage"][FILTER_START:FILTER_STOP]*1000, 
                                 self.data[trial]["Current"][FILTER_START:FILTER_STOP], 
                                 label=_label)
                    elif voltage_units.upper() == "MV" and current_units.upper() == "MA":
                        plt.xlabel("Voltage (mV)", fontsize=16)
                        plt.ylabel("Current (mA)", fontsize=16)
                        plt.plot(self.data[trial]["Corrected Voltage"][FILTER_START:FILTER_STOP]*1000, 
                                 self.data[trial]["Current"][FILTER_START:FILTER_STOP]*1000, 
                                 label=_label)
                    else:
                        raise NameError("Units not found or other error while plotting corrected CV curve.")
                
                except:
                    raise Exception(f"Error while plotting {type_} curve for trial {trial} in {self.name}")
        
        if type_ == "Tafel":
            
            if xlabel == None:
                xlabel = xlabel = r'$\eta \: (V v. {E}^{0})$'
            if ylabel == None:
                ylabel = r'$log|i| \: (i \ in \ A)$'
            
            plt.xlabel(xlabel, fontsize=fontsize)
            plt.ylabel(ylabel, fontsize=fontsize)
            
            xlim = [None,None]
            ylim = [None,None]
            
            for trial in trials:
                try:
                    X = self.data[trial]["Tafel overpotential"]
                    Y = self.data[trial]["Tafel current"]

                    # if not label is provided, use default label, if available
                    _label = label
                    if label is None:
                        _label = self.metadata[trial]["default label"]
                    
                    # Plot log(i) v. iR-corrected overpotential
                    plt.scatter(X, 
                                Y,
                                label=_label)
                    
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
                    raise Exception(f"Error while plotting {type_} curve for trial {trial} in {self.name}")
            
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
            
            plt.xlabel(xlabel, fontsize=fontsize)
            plt.ylabel(ylabel, fontsize=fontsize)
            
            # list of plotting colors (for matching scatter color with regression color)
            # needs to be expanded beyond 10 values to prevent overlap with >10 traces
            COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']
            
            for trial in trials:
                try:
                    COTTRELL_INDEX_MIN = self.data[trial].COTTRELL_INDEX_MIN
                    X = self.data[trial]["Time-Inverse-Root"][COTTRELL_INDEX_MIN:]
                    Y = self.data[trial]["Current"][COTTRELL_INDEX_MIN:]*1000

                    # if no label is provided, use default label, if available
                    _label = label
                    if label is None:
                        _label = self.metadata[trial]["default label"]
                    
                    # determine appropriate color for trace
                    NUMBER_OF_EXISTING_SCATTER_TRACES = len(plt.gca().collections)
                    COLOR_INDEX = NUMBER_OF_EXISTING_SCATTER_TRACES % len(COLORS)
                    COLOR = str(COLORS[COLOR_INDEX])

                    # Plot time-inverse-root v. current (mA)
                    plt.scatter(X, Y, label=_label, color=COLOR)
                    
                    # Plot linear regression of same
                    plt.plot(X, 
                             self.linear_function(X, 
                                                  self.data[trial].COTTRELL_SLOPE, 
                                                  self.data[trial].COTTRELL_INTERCEPT)
                            , '--', color=COLOR)
                
                except:
                    raise Exception(f"Error while plotting {type_} curve for trial {trial} in {self.name}")
        
        if type_ == "CA":
            
            if xlabel == None:
                xlabel = r"$Time \: (s)$"
            if ylabel == None:
                ylabel = r"$Current \: (A)$"
            
            plt.xlabel(xlabel, fontsize=fontsize)
            plt.ylabel(ylabel, fontsize=fontsize)
            
            for trial in trials:
                try:
                    TIME = None
                    CURRENT = self.data[trial]["Current"]
                    
                    if "Corrected Time" in self.data[trial].columns:
                        TIME = self.data[trial]["Corrected Time"]
                    else:
                        self.correct_time(trials=[trial])
                        TIME = self.data[trial]["Corrected Time"]
                    
                    # if not label is provided, use default label, if available
                    _label = label
                    if label is None:
                        _label = self.metadata[trial]["default label"]
                    
                    # Plot using appropriate units and labeling
                    if current_units.upper() == "A":
                        plt.plot(TIME[FILTER_START:FILTER_STOP], 
                                CURRENT[FILTER_START:FILTER_STOP], 
                                label=_label)
                    elif current_units.upper() == "MA":
                        plt.ylabel("Current (mA)")
                        plt.plot(TIME[FILTER_START:FILTER_STOP], 
                                CURRENT[FILTER_START:FILTER_STOP]*1000, 
                                label=_label)
                    else:
                        raise NameError("Units not found or other error while plotting CA data.")
                
                except:
                    raise Exception(f"Error while plotting {type_} curve for trial {trial} in {self.name}")
        
        if type_ == "Nyquist":
            
            if xlabel == None:
                xlabel = r"${Z}_{re} \: (\Omega)$"
            if ylabel == None:
                ylabel = r"${-Z}_{im} \: (\Omega)$"
            
            plt.xlabel(xlabel, fontsize=fontsize)
            plt.ylabel(ylabel, fontsize=fontsize)
            
            try:
                for trial in trials:
                    # if no label is provided, use default label, if available
                    _label = label
                    if label is None:
                        _label = self.metadata[trial]["default label"]
                    
                    plt.scatter(self.data[trial]["Z_real"][FILTER_START:FILTER_STOP],
                                abs(self.data[trial]["Z_imag"][FILTER_START:FILTER_STOP]),
                                label=_label)
            except:
                raise Exception(f"Error while plotting {type_} curve for trial {trial} in {self.name}")
            
            try:
                
                # auto-formatting of Nyquist plot using _zoom param (should normally only be called by self.extract_Ru_from_HFR())
                if (_zoom[0] is not None) and (_zoom[1] is not None):
                    ZOOM_RATIO = _zoom[0]
                    CENTER_INDEX = _zoom[1]
                    
                    X_CENTER = self.data[trial]["Z_real"][CENTER_INDEX]
                    Y_CENTER = abs(self.data[trial]["Z_imag"][CENTER_INDEX])
                    
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
                    '''
                # auto-formatting of Nyquist plot using x_zoom and y_zoom params
                elif (x_zoom != None,None) or (y_zoom != None,None):
                    if (x_zoom != None,None):
                        plt.xlim(x_zoom)
                    if (y_zoom != None,None):
                        plt.ylim(y_zoom)
                    '''
                # In absence of any zoom being specified, Set the x and y axis lengths equal by increasing the limits of the shorter axis to equal the longer
                else:
                    if abs(self.data[trial]["Z_real"]).max() > abs(self.data[trial]["Z_imag"]).max():
                        # set the y axis limits the same as the x axis limits
                        plt.ylim(plt.xlim())
                    else:
                        # set the x axis limits the same as the y axis limits
                        plt.xlim(plt.ylim())
                
                # Set a square aspect ratio
                plt.gca().set_aspect('equal', 'box')
                plt.gca().set_adjustable("box")
            
            except:
                raise Exception(f"Error while auto-zooming {type_} curve for trial {trial} in {self.name}")
        
        if type_ == "Bode":

            if xlabel == None:
                xlabel = "log(freq /Hz)"
            if ylabel == None:
                ylabel = r"log(|Z| /$\Omega)$"
            
            ylabel2 = r"$Phase(Z)/deg$"
            
            plt.xlabel(xlabel, fontsize=fontsize)
            plt.ylabel(ylabel, fontsize=fontsize)
            
            for trial in trials:
                plt.scatter(np.log10(self.data[trial]["Frequency"]), 
                            np.log10(abs(self.data[trial]["Z_real"])), 
                            color='blue')
            
            ax = plt.gca()
            ax.spines['left'].set_color('blue')
            ax.yaxis.label.set_color('blue')
            ax.tick_params(axis='y', colors='blue')

            ax2 = plt.twinx(ax)
            plt.ylabel(ylabel2, fontsize=fontsize)
            
            for trial in trials:
                plt.scatter(np.log10(self.data[trial]["Frequency"]), 
                            self.data[trial]["Phase Angle"], 
                            color='red')
            
            ax2.spines['right'].set_color('red')
            ax2.yaxis.label.set_color('red')
            ax2.tick_params(axis='y', colors='red')
            plt.yticks(fontsize=fontsize-6)
        
        if type_ == "IV":
            
            self.extract_IV_values()
            
            if xlabel == None:
                if density:
                    xlabel = r"Current Density ($mA \, {cm}^{-2}$)"
                else:
                    xlabel = "Current (mA)"
            if ylabel == None:
                ylabel = "Voltage (V)"
            ylabel2 = None
            if density:
                ylabel2 = r"Power Density ($W \, {cm}^{-2}$)"
            else:
                ylabel2 = "Power (W)"
            
            plt.xlabel(xlabel, fontsize=fontsize)
            plt.ylabel(ylabel, fontsize=fontsize)
            
            try:
                CURRENT = []
                VOLTAGE = []
                POWER = []
                
                for trial in trials:
                    if self.metadata[trial]['technique'].upper() == "CP":
                        if density:
                            CURRENT.append(self.metadata[trial]['I-V Current Density (mA/cm2)'])
                            VOLTAGE.append(self.metadata[trial]['I-V Voltage (V)'])
                            POWER.append(self.metadata[trial]['I-V Power Density (W/cm2)'])
                        else:
                            CURRENT.append(self.metadata[trial]['I-V Current (mA)'])
                            VOLTAGE.append(self.metadata[trial]['I-V Voltage (V)'])
                            POWER.append(self.metadata[trial]['I-V Power (W)'])
                
                plt.scatter(CURRENT[FILTER_START:FILTER_STOP],
                            VOLTAGE[FILTER_START:FILTER_STOP], 
                            label=label,
                            color="blue",
                            alpha=0.6,
                            edgecolors='black'
                            )
                plt.plot(CURRENT[FILTER_START:FILTER_STOP],
                         VOLTAGE[FILTER_START:FILTER_STOP], 
                         label=None,
                         color="blue",
                         alpha=0.6,
                         )
                
                ax1 = plt.gca()
                ax1.spines['left'].set_color('blue')
                ax1.yaxis.label.set_color('blue')
                ax1.tick_params(axis='y', colors='blue')

                ax2 = plt.twinx(ax1)
                plt.ylabel(ylabel2, fontsize=fontsize)
                plt.scatter(CURRENT[FILTER_START:FILTER_STOP],
                            POWER[FILTER_START:FILTER_STOP], 
                            label=label,
                            color="red",
                            alpha=0.6,
                            edgecolors='black'
                            )
                plt.plot(CURRENT[FILTER_START:FILTER_STOP],
                         POWER[FILTER_START:FILTER_STOP], 
                         label=None,
                         color="red",
                         alpha=0.6,
                         )
                ax2.spines['right'].set_color('red')
                ax2.yaxis.label.set_color('red')
                ax2.tick_params(axis='y', colors='red', labelsize=fontsize-6)
            
            except:
                raise Exception(f"Error while plotting {type_} curve for trial {trial} in {self.name}")
        
        if type_ == "OCP":
            
            # If labels are not provided in method call, use default axis labels for chronopotentiometric data
            if xlabel == None:
                xlabel = "Time (hr)"
            if ylabel == None:
                ylabel = "Voltage (V)"
            plt.xlabel(xlabel, fontsize=fontsize)
            plt.ylabel(ylabel, fontsize=fontsize)
            
            # Iterative treatment of desired trials
            for trial in trials:
                try:
                    
                    # Fetch data series to plot
                    TIME = self.data[trial]["Time"] / 3600
                    VOLTAGE = self.data[trial]["Voltage"]
                        
                    # if no label is provided, use default label, if available
                    _label = label
                    if label is None:
                        _label = self.metadata[trial]["default label"]
                    
                    # Plot using appropriate units and labeling
                    if voltage_units.upper() == "V":
                        plt.plot(TIME,
                                 VOLTAGE[FILTER_START:FILTER_STOP],
                                 label=_label)
                    elif voltage_units.upper() == "MV":
                        plt.ylabel("Voltage (mV)")
                        plt.plot(TIME,
                                 VOLTAGE[FILTER_START:FILTER_STOP]*1000,
                                 label=_label)
                    
                    # Set x-axis time ticks to every 60 minutes
                    '''
                    length_of_longest_trial = 0
                    for line in plt.gca().get_lines():
                        time_line = line.get_xdata()
                        if max(time_line) > length_of_longest_trial:
                            length_of_longest_trial = max(time_line)
                    
                    plt.xticks(np.arange(0, length_of_longest_trial+1, 1.0))
                    '''
                
                except:
                    raise Exception(f"Error while plotting {type_} curve for trial {trial} in {self.name}")
        
        
        # Zoom using x_zoom and y_zoom params
        if (x_zoom != None,None) or (y_zoom != None,None):
            if (x_zoom != None,None):
                plt.xlim(x_zoom)
            if (y_zoom != None,None):
                plt.ylim(y_zoom)
        
        # Plot legend, if applicable
        if legend:
            plt.legend(fontsize=fontsize-6)  
    
    
    def tafel(self, trials=[], voltage_bounds=[None,None], E_eq=0.0):
        '''
        Perform tafel analysis on the specified trial.
        '''

        try:
            # if a list of trials is not passed, default is to perform analysis on all CV data
            if not trials:
                for trial in self.metadata:
                    if self.metadata[trial]["technique"] == "CV":
                        trials.append(trial)
        
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
        
        except:
            raise Exception(f"Error while selecting trials and setting up for Tafel analysis for {self.name}.")
        
        for trial in trials:
            try:
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
            except:
                raise Exception(f"Error while pre-processing and/or filtering data prior to Tafel analysis for {self.name} trial {trial}.")

            try:
                # Calculate overpotential and log(i) from voltage and current data
                overpotential = self.data[trial]["Tafel overpotential"] = tafel_data['Corrected Voltage'] - E_eq
                self.update_history(f"Tafel overpotential calculated and saved for trial {trial}.")
                log_i = self.data[trial]["Tafel current"] = np.log10(abs(tafel_data['Current']))
                self.update_history(f"Tafel current calculated and saved for trial {trial}.")
                
                # Define highest index x value on which to perform Tafel regression (this fudge factor will ideally be eliminated in the future)
                index_limit = round(len(overpotential)/4)
                
                # Extract regression parameters
                tafel_slope, tafel_yint = np.polyfit(overpotential[:index_limit], log_i[:index_limit], 1)
            except:
                raise Exception(f"Error while calculating overpotential, log(i), Tafel slope, or other parameter during Tafel analysis for trial {trial} in {self.name}.")
            
            try:
                # Save extracted values
                # Tafel slope (V/dec)
                self.data[trial].tafel_slope_V = tafel_slope
                self.update_history(f"Tafel slope (V) calculated and saved from Tafel analysis of trial {trial}.")
                # Tafel y-intercept
                self.data[trial].tafel_yint = tafel_yint
                self.update_history(f"Tafel y-intercept calculated and saved from Tafel analysis of trial {trial}.")
                # Tafel slope (mv/dec)
                self.data[trial].tafel_slope_mV = tafel_slope * 1000
                self.update_history(f"Tafel slope (mV) calculated and saved from Tafel analysis of trial {trial}.")
                # exchange current (mA)
                self.data[trial].i_0 = 10**tafel_yint * 1000
                self.update_history(f"Exchange current (mA) calculated and saved from Tafel analysis of trial {trial}.")
            except:
                raise Exception(f"Error while saving extracted values after Tafel analysis for trial {trial} in {self.name}.")
    

    def update_history(self, entry):
        '''
        Append a datetime:entry entry to the object's .history dict object
        '''
        
        try:
            self.history[datetime.today().strftime('%Y-%m-%d %H:%M:%S.%f')] = entry
        except:
            raise Exception(f"Error while updating history with the following entry: {entry}")
    
    
    # Setting the attributes below as 'properties' ensures they are updated (if needed) when other attributes are changed
    # e.g., one important application here is that it allows the object history to auto-update when an attribute is changed

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




# =============================================================================
# Below this needs to be cleaned up and harmonized
# =============================================================================

# This is a duplicate with linear_function above, which is newer; need to compare and downselect
def Regression(overpotential, log_i, step=False):
    if step == False:
        SLOPE, INTERCEPT, RVALUE, PVALUE, STDERR = stats.linregress(overpotential, log_i)
    else:
        SLOPE, INTERCEPT, RVALUE, PVALUE, STDERR = ReverseStepRegression(overpotential, log_i)
    return SLOPE, INTERCEPT
    

# Extract Exchange Current and Exchange Rate Constant (from original pre-chemax Tafel package)
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
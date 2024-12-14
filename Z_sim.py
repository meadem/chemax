import numpy as np
import matplotlib.pyplot as plt

'''
    - To do:
	- Figure out how to make Z-sim more extensible and general use with usage of plt
		- (for example I want to be able to easily custom make a 2x2 image grid on the fly using Z-sim)
    - explore how to make a dashboard UI
    - add a print or export method
'''

class Z_sim:
    '''
    Generate simulated impedance spectroscopy data from a specified circuit model.

    Parameters
    ----------
    circuit : string
        A string used to describe the circuit model; tries to match case and returns error msg if no match.
    frequency_specification : string
        A string that defines whether to use frequency_range and log_spacing to create an array ("range"), 
        or to use the predefined frequency_array ("array").
    frequency_range : tuple
        Describes minimum (frequency_range[0]) and maximum (frequency_range[1]) values for the frequency range.
    frequency_array : numpy 1d array
        Predefined array of frequency values.
    R1 : float
        Resistor 1
    R2 : float
        Resistor 2
    R3 : float
        Resistor 3
    C1 : float
        Capacitor 1
    C2 : float
        Capacitor 2
    C3 : float
        Capacitor 3
    L1 : float
        Inductor 1
    L2 : float
        Inductor 2
    L3 : float
        Inductor 3
    a1 : float
        Coefficient for CPE 1
    a2 : float
        Coefficient for CPE 2
    a3 : float
        Coefficient for CPE 3
    s1 : float
        Coefficient for Warburg element 1
    s2 : float
        Coefficient for Warburg element 2
    s3 : float
        Coefficient for Warburg element 3
    data_points : int
        Number of data points to use in simulation
    log_spacing : bool
        If True, use log spacing for frequency values. If False, linear spacing will be used.

    Returns
    -------
    class object instance
        Calculated impedance and frequency data for the specified circuit model, with methods for plotting.
    
    '''
    
    def __init__(self, 
                 circuit="R1/C1", 
                 frequency_specification="range",
                 frequency_range=[0.000001,1000000], 
                 frequency_array=np.array([0]),
                 R1=1.0, 
                 C1=1.0, 
                 R2=1.0, 
                 C2=1.0, 
                 R3=1.0, 
                 C3=1.0, 
                 L1=1.0, 
                 L2=1.0, 
                 L3=1.0, 
                 a1=1.0, 
                 a2=1.0, 
                 a3=1.0, 
                 s1=1.0, 
                 s2=1.0, 
                 s3=1.0, 
                 data_points=100, 
                 log_spacing=True
                ):
        
        
        if frequency_specification == "range":
            frequency = self.spacing(frequency_range, log_spacing, data_points)
            angular_frequency = 2 * np.pi * frequency
        
        elif frequency_specification == "array":
            frequency = frequency_array
            angular_frequency = 2 * np.pi * frequency
        
        else:
            raise NameError("Frequency specification '" + frequency_specification + "' not found.")
        
        
        Z_R1 = self.resistor_impedance(R1, angular_frequency)
        Z_R2 = self.resistor_impedance(R2, angular_frequency)
        Z_R3 = self.resistor_impedance(R3, angular_frequency)
        Z_C1 = self.capacitor_impedance(C1, angular_frequency)
        Z_C2 = self.capacitor_impedance(C2, angular_frequency)
        Z_C3 = self.capacitor_impedance(C3, angular_frequency)
        Z_L1 = self.inductor_impedance(L1, angular_frequency)
        Z_L2 = self.inductor_impedance(L2, angular_frequency)
        Z_L3 = self.inductor_impedance(L3, angular_frequency)
        Z_Q1 = self.CPE_impedance(C1, a1, angular_frequency)
        Z_Q2 = self.CPE_impedance(C2, a2, angular_frequency)
        Z_Q3 = self.CPE_impedance(C3, a3, angular_frequency)
        Z_W1 = self.Warburg_impedance(s1, angular_frequency)
        Z_W2 = self.Warburg_impedance(s2, angular_frequency)
        Z_W3 = self.Warburg_impedance(s3, angular_frequency)
        
        
        if circuit == "R1":
            self.Z = Z_R1
        
        
        elif circuit == "C1":
            self.Z = Z_C1
        
        
        elif circuit == "Q1":
            self.Z = Z_Q1
        
        
        elif circuit == "L1":
            self.Z = Z_L1
        
        
        elif circuit == "R1+C1":
            self.Z = self.sum_series_impedance(Z_R1, Z_C1)
        
        
        elif circuit == "R1+Q1":
            self.Z = self.sum_series_impedance(Z_R1, Z_Q1)
        
        
        elif circuit == "R1+L1":
            self.Z = self.sum_series_impedance(Z_R1, Z_L1)
        
        
        elif circuit == "R1/C1":
            self.Z = self.sum_parallel_impedance(Z_R1, Z_C1)
        
        
        elif circuit == "R1/Q1":
            self.Z = self.sum_parallel_impedance(Z_R1, Z_Q1)
        
        
        elif circuit == "R1/L1":
            self.Z = self.sum_parallel_impedance(Z_R1, Z_L1)
        
        
        elif circuit == "R1+R2/C2":
            self.Z = self.sum_series_impedance(Z_R1, 
                                               self.sum_parallel_impedance(Z_R2, Z_C2)
                                              )
        
        
        elif circuit == "R1+R2/Q2":
            self.Z = self.sum_series_impedance(Z_R1, 
                                               self.sum_parallel_impedance(Z_R2, Z_Q2)
                                              )
        
        
        elif circuit == "R1/C1+R2/C2":
            self.Z = self.sum_series_impedance(self.sum_parallel_impedance(Z_R1, Z_C1), 
                                               self.sum_parallel_impedance(Z_R2, Z_C2)
                                              )
        
        
        elif circuit == "R1/C1+R2/Q2":
            self.Z = self.sum_series_impedance(self.sum_parallel_impedance(Z_R1, Z_C1), 
                                               self.sum_parallel_impedance(Z_R2, Z_Q2)
                                              )
        
        
        elif circuit == "R1/Q1+R2/Q2":
            self.Z = self.sum_series_impedance(self.sum_parallel_impedance(Z_R1, Z_Q1), 
                                               self.sum_parallel_impedance(Z_R2, Z_Q2)
                                              )
        
        
        elif circuit == "R1+R2/C2+R3/C3":
            # aka Voigt circuit
            self.Z = self.sum_series_impedance(Z_R1,
                                               self.sum_parallel_impedance(Z_R2, Z_C2),
                                               self.sum_parallel_impedance(Z_R3, Z_C3)
                                               )
        
        
        elif circuit == "R1+R2/C2+R3/Q3":
            self.Z = self.sum_series_impedance(Z_R1,
                                               self.sum_parallel_impedance(Z_R2, Z_C2),
                                               self.sum_parallel_impedance(Z_R3, Z_Q3)
                                               )
        
        
        elif circuit == "R1+R2/Q2+R3/Q3":
            self.Z = self.sum_series_impedance(Z_R1,
                                               self.sum_parallel_impedance(Z_R2, Z_Q2),
                                               self.sum_parallel_impedance(Z_R3, Z_Q3)
                                               )
        
        
        elif circuit == "R1/C1+R2/C2+R3/C3":
            self.Z = self.sum_series_impedance(self.sum_parallel_impedance(Z_R1, Z_C1),
                                               self.sum_parallel_impedance(Z_R2, Z_C2),
                                               self.sum_parallel_impedance(Z_R3, Z_C3)
                                               )
        
        
        elif circuit == "R1/Q1+R2/C2+R3/C3":
            self.Z = self.sum_series_impedance(self.sum_parallel_impedance(Z_R1, Z_Q1),
                                               self.sum_parallel_impedance(Z_R2, Z_C2),
                                               self.sum_parallel_impedance(Z_R3, Z_C3)
                                               )
        
        
        elif circuit == "R1/Q1+R2/Q2+R3/C3":
            self.Z = self.sum_series_impedance(self.sum_parallel_impedance(Z_R1, Z_Q1),
                                               self.sum_parallel_impedance(Z_R2, Z_Q2),
                                               self.sum_parallel_impedance(Z_R3, Z_C3)
                                               )
        
        
        elif circuit == "R1/Q1+R2/Q2+R3/Q3":
            self.Z = self.sum_series_impedance(self.sum_parallel_impedance(Z_R1, Z_Q1),
                                               self.sum_parallel_impedance(Z_R2, Z_Q2),
                                               self.sum_parallel_impedance(Z_R3, Z_Q3)
                                               )
        
        
        elif circuit == "R1+C2/(W2+R2)":
            # aka Randles equivalent circuit
            self.Z = self.sum_series_impedance(Z_R1, 
                                               self.sum_parallel_impedance(Z_C2,
                                                                           self.sum_series_impedance(Z_W2, Z_R2)
                                                                          )
                                              )
        
        
        elif circuit == "R1+C2/(R2+(C3/R3))":
            # aka 'ladder'
            self.Z = self.sum_series_impedance(Z_R1, 
                                               self.sum_parallel_impedance(Z_C2,
                                                                           self.sum_series_impedance(Z_R2, 
                                                                                                     self.sum_parallel_impedance(Z_C3, Z_R3))
                                                                          )
                                              )
        
        
        elif circuit == "(R1+(C2/R2))/(R3+C3)":
            # aka 'mixed' circuit
            self.Z = self.sum_parallel_impedance(self.sum_series_impedance(Z_R1, 
                                                                           self.sum_parallel_impedance(Z_C2, Z_R2)),
                                                 self.sum_series_impedance(Z_R3, Z_C3)
                                                )
        
        
        else:
            raise NameError("CIRCUIT MODEL '" + circuit + "' NOT FOUND.")
        
        
        self.circuit = circuit
        self.frequency_range = frequency_range
        self.frequency = frequency
        self.angular_frequency = angular_frequency
        self.data_points = data_points
        self.log_spacing = log_spacing
        self.R1 = R1
        self.C1 = C1
        self.R2 = R2
        self.C2 = C2
        self.R3 = R3
        self.C3 = C3
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.s1 = s1
        self.s2 = s2
        self.s3 = s3
    
    
    def spacing(self, 
                range_boundaries, 
                is_log_spacing, 
                number_of_points
               ):
        
        import numpy as np
        
        low_boundary = range_boundaries[0]
        high_boundary = range_boundaries[1]

        if is_log_spacing:
            log_start = np.log10(low_boundary)
            log_stop = np.log10(high_boundary)
            frequency_exponents = np.linspace(log_start, log_stop, number_of_points)
            frequency = 10**frequency_exponents
            return frequency

        else:
            frequency = np.linspace(low_boundary, high_boundary, number_of_points)
            return frequency
    
    
    def resistor_impedance(self, R, angular_frequency):
        w = angular_frequency
        result = np.full( (1,len(w)), (R + 0j), dtype=complex )
        return result
    
    
    def capacitor_impedance(self, C, angular_frequency):
        w = angular_frequency
        result = 0 - 1j*( w * C )**-1
        return result
    
    
    def inductor_impedance(self, L, angular_frequency):
        w = angular_frequency
        result = 0 + 1j*( w * L )
        return result
    
    
    def CPE_impedance(self, C, a, angular_frequency):
        w = angular_frequency
        result = ( np.cos( a * np.pi / 2 ) / ( C * w**a ) ) - 1j*( np.sin( a * np.pi / 2 ) / ( C * w**a ) )
        return result
    
    
    def Warburg_impedance(self, s, angular_frequency):
        w = angular_frequency
        result = ( s / np.sqrt(w) ) - 1j*( s / np.sqrt(w) )
        return result
    
    
    def sum_series_impedance(self, *args):
        result = np.sum(args, axis=0)
        return result
    
    
    def sum_parallel_impedance(self, Z1, Z2):
        result = ( Z1**-1 + Z2**-1 )**-1
        return result
    
    
    def nyquist(self):
        
        fig = plt.figure(1, figsize=(7,7))
        ax = fig.add_subplot(1,1,1)
        ax.scatter(self.Z.real, -self.Z.imag, linewidth = 3)

        ax.set_xlabel(r'${Z}_{re} \: (\Omega)$', fontsize = 16)
        ax.set_ylabel(r'${-Z}_{im} \: (\Omega)$', fontsize = 16)

        ax.axhline(y=0, color='k')
        ax.axvline(x=0, color='k')
        
        if abs(self.Z.real).max() > abs(self.Z.imag).max():
            ax.set_ylim(ax.get_xlim())
        
        else:
            ax.set_xlim(ax.get_ylim())
        
        ax.set_aspect('equal', 'box')


    def bode(self):

        fig = plt.figure(2, figsize=(7,7))
        ax = fig.add_subplot(1,1,1)
        
        ax.scatter(np.log10(self.frequency), np.log10(abs(self.Z.real)), color='blue')
        ax.set_xlabel(r'$log(freq/Hz)$', fontsize = 16)
        ax.set_ylabel(r'$log(|Z|/\Omega)$', fontsize = 16)
        ax.spines['left'].set_color('blue')
        ax.yaxis.label.set_color('blue')
        ax.tick_params(axis='y', colors='blue')
        
        ax2 = ax.twinx()
        ax2.scatter(np.log10(self.frequency), np.angle(self.Z, deg=True), color='red')
        ax2.set_ylabel(r'$Phase(Z)/deg$', fontsize=16)
        ax2.spines['right'].set_color('red')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')
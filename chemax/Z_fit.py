from test import Z_sim
from scipy.optimize import least_squares, curve_fit

import numpy as np
import matplotlib.pyplot as plt
'''
to do:
 - implement error function and iteration
 - add docstring
 - 
'''
class Z_fit:
    '''
    docstring here
    '''
    
    # need a method to automate adding and evaluating individual circuit elements until a satisfactory fit is found 
    # I can just define a sequence of circuits to try, then use the F-test to justify/discard the addition of a new element
    
    def __init__(self, Z_real, Z_imag, frequency):
        
        # how are the first parameter guesses obtained? must they be provided?
        # do we cycle through all circuit models every time?
        # we could initialize by providing the dataset, then invoke methods to fit against models.
        # the above approach could allow us the memory to compare guesses, as well. For example, is adding a new component justified mathematically?
        
        import numpy as np
        import matplotlib.pyplot as plt
        
        self.Z = Z_real + 1j*Z_imag
        self.frequency = frequency
        self.angular_frequency = 2 * np.pi * frequency
        self.current_circuit_model = ""
        
        self.R1 = 0.0
        self.R2 = 0.0
        # etc?
    
    
    def fit_RC(self):
        
        # I need exit conditions for the minimization function; a percent error and/or max number of guesses
        # I want to minimize error for fit for both the imaginary plane and bode plots; how do I do this?
        
        R1_guess = self.Z.real.min()
        C1_guess = 5e-4    # replace with something else later
        
        best_guess = least_squares(self.get_residuals, [R1_guess, C1_guess], method="lm")
        
        self.current_circuit_model = "R1/C1"
        #result = self.report()
        #return result 
        return best_guess
    
        # what do we want to return when this is invoked? 
        # Without a return we can change the state of the object (I think..), so all we 
        # need to return is some sort of report? Maybe metrics of fit quality, 
        # graphs of residuals, input data plotted against the best guess
        # this would be a good place to implement extensible plotting
        # in order to use Z_sim to generate data for Z_fit, I need to be able to send frequency data to Z_sim instead of it generating its own frequency array every time its invoked
    
    
    def reduced_sum_of_squares(self, guess=[0, 0], weighting="modulus", number_of_model_parameters=0):
        
        x_guess = guess[0]
        y_guess = guess[1]
        
        # Calculate weighting factor(s) using the specified method
        if weighting == "unit":
            real_weighting_factor = imag_weighting_factor = 1

        elif weighting == "modulus":
            real_weighting_factor = imag_weighting_factor = self.modulus_weighting(x_guess, y_guess)
        
        elif weighting == "proportional":
            real_weighting_factor, imag_weighting_factor = self.proportional_weighting(x_guess, y_guess)
        
        else:
            raise NameError("Weighting method '" + weighting + "' not found.")
        
        # S is the weighted sum of squares, and functions as the error signal here
        S = np.sum( real_weighting_factor*( self.Z.real - guess.Z.real )**2 + imag_weighting_factor*( self.Z.imag - guess.Z.imag )**2 )
        
        # dof is degrees of freedom
        dof = 2*(count(self.frequency)) - number_of_model_parameters
        
        # S_v is the reduced sum of squares
        S_v = S / dof
        
        return S_v
    
    
    def reduced_sum_of_squares_test(self, guess=[0, 0]):
        
        x_guess = guess[0]
        y_guess = guess[1]
        
        guess_array = self.send_RC_to_Z_sim(guess[0], guess[1])
        
        real_weighting_factor = imag_weighting_factor = self.modulus_weighting(x_guess, y_guess)
        
        # S is the weighted sum of squares, and functions as the error signal here
        S = np.sum( real_weighting_factor*( self.Z.real - guess_array.Z.real )**2 + imag_weighting_factor*( self.Z.imag - guess_array.Z.imag )**2 )
        
        # dof is degrees of freedom
        number_of_model_parameters = 2
        dof = 2*(len(self.frequency)) - number_of_model_parameters
        
        # S_v is the reduced sum of squares
        S_v = S / dof
        
        return S_v
    
    
    def get_residuals(self, guess=[0, 0]):
        
        x_guess = guess[0]
        y_guess = guess[1]
        
        guess_array = self.send_RC_to_Z_sim(guess[0], guess[1])
        
        residuals = self.Z.imag - guess_array.Z.imag
        
        print([residuals])
        
        return residuals
    
    
    def modulus_weighting(self, x, y):
        result = np.sqrt( x**2 + y**2 )
        return result
    
    
    def proportional_weighting(self, x, y):
        x_result = x**-2
        y_result = y**-2
        return x_result, y_result
    
    
    def report(self):
        
        fig = plt.figure(figsize=(10,10))    # define figure
        ax1 = fig.add_subplot(2,2,1)       # define axes1: Nyquist
        ax2 = fig.add_subplot(2,2,2)       # define axes2: Bode
        #ax3 = fig.add_subplot(2,2,3)       # define axes3: Residuals
        
        
        # Nyquist
        ax1.scatter(self.Z.real, -self.Z.imag, linewidth = 3, label='Experimental')
        # ax1.scatter({REAL_FIT}, {-IMAG_FIT}, linewidth = 3, label='Fit') # Need a distinct marker
        ax1.set_xlabel(r'${Z}_{re} \: (\Omega)$', fontsize = 16)
        ax1.set_ylabel(r'${-Z}_{im} \: (\Omega)$', fontsize = 16)
        ax1.axhline(y=0, color='k')
        ax1.axvline(x=0, color='k')
        
        if abs(self.Z.real).max() > abs(self.Z.imag).max():
            ax1.set_ylim(ax1.get_xlim())
        else:
            ax1.set_xlim(ax1.get_ylim())
        
        ax1.set_aspect('equal', 'box')
        
        
        # Bode
        ax2.scatter(np.log10(self.frequency), np.log10(abs(self.Z.real)), color='blue', label='Experimental') 
        # ax2.scatter(np.log10(self.frequency), np.log10(abs({REAL_FIT})), color='blue', label='Fit') # Need a distinct marker
        ax2.set_xlabel(r'$log(freq/Hz)$', fontsize = 16)
        ax2.set_ylabel(r'$log(|Z|/\Omega)$', fontsize = 16)
        ax2.spines['left'].set_color('blue')
        ax2.yaxis.label.set_color('blue')
        ax2.tick_params(axis='y', colors='blue')
        
        ax2_ = ax2.twinx()
        ax2_.scatter(np.log10(self.frequency), np.angle(self.Z, deg=True), color='red', label='Experimental')
        # ax2_.scatter(np.log10(self.frequency), {PHASEANGLE_FIT}, deg=True), color='red', label='Fit')
        ax2_.set_ylabel(r'$Phase(Z)/deg$', fontsize=16)
        ax2_.spines['right'].set_color('red')
        ax2_.yaxis.label.set_color('red')
        ax2_.tick_params(axis='y', colors='red')
        
        
        # need to capture residuals so I can plot them here
        
        # Can I print out some metrics in quadrant 4?
    
    
    def send_RC_to_Z_sim(self, R1, C1):
        result = Z_sim(circuit="R1/C1", R1=R1, C1=C1, frequency_specification="array", frequency_array=self.frequency)
        return result


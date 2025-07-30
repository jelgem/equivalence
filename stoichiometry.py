import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
from scipy.stats import beta
from scipy import interpolate, integrate


plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11
plt.rcParams['axes.labelsize'] = 14

# = = = = = = = = = = = = = = = = = = = = = = = = = #

class stoichiometry:
    def __init__(self, Fuel_dict, Ox_dict, ylength = 10001):
        # STEP 0: Obtain the stoichiometic properties
        stoich = ct.Solution('gri30.yaml')
        stoich.set_equivalence_ratio(1.0, Fuel_dict, Ox_dict)
        
        # Stoichiometric fuel mass fraction
        self.Yfuel = stoich.Y[stoich.species_index(list(Fuel_dict.keys())[0])] 

        # Stoichiometric Air/Fuel Ratio
        self.s = (1 - self.Yfuel) / self.Yfuel
        
        # Define the random-sample mass fraction and equvialence ratios
        self.Y = np.linspace(0,1.0,ylength)
        self.Phi = self.Phifunc(self.Y)       
    
    #================================================================#
    # Everything you might want to know about a beta distribution
    #================================================================#
    def beta(self, N, verbose = True):
        # Shape Parameters
        a = N * self.Yfuel; b = N * (1 - self.Yfuel)
        # Beta Distribution (of y)
        self.B = beta.pdf(self.Y, a, b)
        # Beta Distribution of Phi
        self.Bphi = self.B * self.s / (self.s + self.Phi)**2
        
        # Medians
        self.medianY    = self.numericalMedian(self.Y, self.B)
        self.medianPhi  = self.numericalMedian(self.Phi, self.Bphi)
        
        # Means
        self.meanY    = self.numericalMean(self.Y, self.B)
        self.meanPhi  = self.numericalMean(self.Phi, self.Bphi)
        
        # Compute the value of the distribution at the medians
        self.medianB     = np.interp(self.medianY, self.Y, self.B)
        self.medianBphi  = np.interp(self.medianPhi, self.Phi, self.Bphi)
        
        if verbose:
            print(f'# - - - - - - - - - - - - - - - - - - - - - #\n\t\t\t\t N = {N}')
            print(f'Median = {self.medianY:1.5f} Mass Fraction, {self.medianPhi:1.5f} ER')
            print(f'Mean   = {self.meanY:1.5f} Mass Fraction, {self.meanPhi:1.5f} ER')
        
    
    #================================================================#
    # Cumulative Distribution Function
    #================================================================#
    
    def CDFunc(self, x, f, low = 0, high = 1e8):
        # Ensure inputs are numpy arrays and sorted
        x = np.asarray(x)
        f = np.asarray(f)
    
        # Interpolate the PDF (linear or cubic for smoothness)
        pdf_interp = interpolate.interp1d(x, f, kind='cubic', fill_value=0.0, bounds_error=False, assume_sorted=True)
    
        # Numerical integration for the CDF between [low, high]
        cdf_val, _ = integrate.quad(pdf_interp, low, high)
        return(cdf_val)
    
    #================================================================#
    # Mean Numerator Integral
    #================================================================#

    def MeanIntegral(self, x, f, low, high):
        """
        Numerically computes the integral of x * f(x) dx from `low` to `high`, given sample points x and PDF values f.
        """            
        # Interpolate x*f(x)
        mean_integrand = interpolate.interp1d(x, x * f, kind='cubic', fill_value=0.0, bounds_error=False, assume_sorted=True)
        
        # Integrate from low to high
        mean_integral, _ = integrate.quad(mean_integrand, low, high)
        return(mean_integral)
    
    #================================================================#
    # System of equations to solve
    #================================================================#

    def SE_phi(self, lowhigh, pcdf=0.9):
        """
        The system of equations to find the roots to identify symmetric points.
        
        lowhigh : List containing the low and high values of the (test) range
        pcdf : the portion of values to expect to find in the range 
        """
        low, high = lowhigh # Split up the list
        # Compute the CDF and the Mean-Integral-Numerator for this range. Return the difference versus the target.
        vex = [self.CDFunc(self.Phi, self.Bphi, low, high) - pcdf, 
                self.MeanIntegral(self.Phi, self.Bphi, low, high) - pcdf]
        # Add a penalty if either value is negative
        for i in range(len(vex)):
            vex[i] *= (1 + 100*(min(lowhigh)**2 * (min(lowhigh)<=0)))
        
        return(vex)
    

        
    # / / / / / / / / / / / / / / / / / / / / / / / / #
    # Functions between non-dimensional metrics
    # \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ #
    def Phifunc(self, y):
        return((y / (1 - y + 1e-17)) / (self.Yfuel / ( 1 - self.Yfuel)))
    
    def Yphifunc(self,phi):
        return(phi / (self.s + phi))
        
    
    
    
    # / / / / / / / / / / / / / / / / / / / / / / / / #
    # Statistical Functions
    # \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ #
    def numericalMedian(self, x, f):
        dx = np.diff(x, append=x[-1])
        cdf = np.cumsum(f * dx)
        cdf /= cdf[-1]
        return(np.interp(0.5, cdf, x))
    
    def numericalMean(self, x, f):
        dx = np.diff(x, append=x[-1])
        mean = np.sum(x * f * dx) / np.sum(f * dx)
        return(np.sum(x * f * dx) / np.sum(f * dx))
    


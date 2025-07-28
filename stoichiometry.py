from scipy.interpolate import interp1d
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cantera as ct
from scipy.stats import beta

plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11
plt.rcParams['axes.labelsize'] = 14

# = = = = = = = = = = = = = = = = = = = = = = = = = #

class stoichiometry:
    def __init__(self, Fuel_dict, Ox_dict, ylength = 5001):
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
        
        
        
    # / / / / / / / / / / / / / / / / / / / / / / / / #
    # Functions between non-dimensional metrics
    # \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ #
    def Phifunc(self, y):
        return((y / (1 - y + 1e-17)) / (self.Yfuel / ( 1 - self.Yfuel)))
        
    
    
    
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
    


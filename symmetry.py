import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from stoichiometry import *
from scipy.optimize import fsolve, curve_fit

# = = = = = = = = = = = = = = = = = = = = = = = = = #

# = = = = = = = = = = = = = = = = = = = = = = = = = #
Fuel_dict   = {'H2': 1.0}
Ox_dict     = {'O2': 0.210085, 'N2': 0.789915}

# Create the stoichiometric mixture
mix = stoichiometry(Fuel_dict, Ox_dict)

# Evaluate the distributions at different samples
Nlist = np.array([100, 250, 500, 1000, 2500, 5000, 10000])

# Use confidence intervals at 10% increments
Plist = np.arange(0.05, 1.0, 0.01)

# Build 2D data structures
Nmat, Pmat = np.meshgrid(Nlist, Plist, indexing='ij')
Phil = Nmat*0.0; Phir = Nmat*0.0

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#%%                     SEARCHING                         #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
for i, N in enumerate(Nlist):
    # Generate the Distribution
    mix.beta(N)
    for j, P in enumerate(Plist):
        # Choose a smart initial guess
        if j > 0:
            guess = [Phil[i, j-1], Phir[i, j-1]]
        elif i > 0:
            guess = [Phil[i-1, j], Phir[i-1, j]]
        else:
            guess = [0.9, 1.1]
        # Solve
        bounds = fsolve(mix.SE_phi, guess, args = P)
        Phil[i, j] = bounds[0]
        Phir[i, j] = bounds[1]


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#%%                  CURVE FITTING                        #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

def f_inverse(x, a, b, c):
    return(c + a / (x + b))

# Fit f_inverse to your data
params, _ = curve_fit(f_inverse, Phil.ravel(), Phir.ravel())
# Define the trendline independent variable
Phil_trend = np.arange(0.01, 1, 0.01)

# Generate the trendline dependent variables using the fitted parameters
Phir_inverse = f_inverse(Phil_trend, *params)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#%%                      PLOTTING                         #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #

fig, axs = plt.subplots(figsize=[6,3.5],layout="constrained")
axs.set_facecolor([1.0]*3)
axs.set_xlabel("Equivalence Ratio")
axs.set_ylabel("Confidence Level")
axs.grid(alpha=0.2, which='both')
axs.set_yticks(np.arange(0.0, 1.01, 0.1))
axs.set_ylim(0, 1)

norm = mcolors.LogNorm(vmin=np.min(Nmat[Nmat>0]), vmax=np.max(Nmat))  # Log

sc0 = axs.scatter(Phil, Pmat, s=14, c = Nmat, norm=norm, cmap = 'rainbow', marker = '<')
sc1 = axs.scatter(Phir, Pmat, s=14, c = Nmat, norm=norm, cmap = 'rainbow', marker = '>')
cbar = plt.colorbar(sc0, ax = axs, label="Sample Size Param.")
# Set colorbar ticks to the values in Nlist (assuming Nlist > 0 and sorted)
cbar.set_ticks(Nlist)
# Format tick labels in usual notation (you can customize formatting as needed)
cbar.ax.set_yticklabels([str(n) for n in Nlist])

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#%%                      PLOTTING                         #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
fig, axs = plt.subplots(figsize=[6,4],layout="constrained")
axs.set_xlabel("Lean Eq. Ratio")
axs.set_ylabel("Rich Eq. Ratio")
axs.grid(alpha=0.2, which='both')
axs.set_xlim([0, 1])

ss0 = axs.scatter(Phil, Phir, s=14, c = Nmat, norm=norm, cmap = 'rainbow', alpha=0.65)

cbar = plt.colorbar(ss0, ax = axs, label="Sample Size Param.")
# Set colorbar ticks to the values in Nlist (assuming Nlist > 0 and sorted)
cbar.set_ticks(Nlist)
# Format tick labels in usual notation (you can customize formatting as needed)
cbar.ax.set_yticklabels([str(n) for n in Nlist])

plt.plot(Phil_trend, Phir_inverse,'k:', alpha = 0.6, 
          label = rf'$\Phi_r = {params[2]:1.4f} + \frac{{{params[0]:1.3f}}}{{\Phi_l + {params[1]:1.3f}}}$')
plt.plot(Phil_trend, 1/Phil_trend, 'k--',alpha=0.4, label = rf'$\Phi_r = \frac{{1}}{{\Phi_l}}$')
plt.legend(fontsize=14)
axs.set_ylim([1,4])


Yfl = mix.Yphifunc(Phil)
Yfr = mix.Yphifunc(Phir)
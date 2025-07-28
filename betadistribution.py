import os, sys
import numpy as np
import matplotlib.pyplot as plt
from stoichiometry import *
# = = = = = = = = = = = = = = = = = = = = = = = = = #

# = = = = = = = = = = = = = = = = = = = = = = = = = #
Fuel_dict   = {'H2': 1.0}
Ox_dict     = {'O2': 0.210085, 'N2': 0.789915}

# Create the stoichiometric mixture
mix = stoichiometry(Fuel_dict, Ox_dict)

# Evaluate the distributions at different samples
Nlist = np.array([100, 250, 500, 1000, 2500, 5000, 10000])



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#%%                      PLOTTING                         #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
xnams = {'Y':'Fuel Mass Fraction', 'Phi':'Equivalence Ratio'}
xvars = list(xnams.keys())
xlims = dict(zip(xvars,    [[0, 0.08], [0, 2]]   ))
figs = {}
axs = {}

for xvar in xvars:
    figs[xvar], axs[xvar] = plt.subplots(figsize=[6,3.3],layout="constrained")
    axs[xvar].grid(alpha = 0.3)
    axs[xvar].set_xlim(*xlims[xvar])
    axs[xvar].set_xlabel(xnams[xvar])
    axs[xvar].set_ylabel("Probability Density")


# = = = = = = = = = = = = = = = = = = = = = = = = = #
colors = plt.get_cmap('rainbow')(mcolors.LogNorm(vmin=Nlist.min(), vmax=Nlist.max())(Nlist))


for i, N in enumerate(Nlist):
    mix.beta(N)
    
    # Preparation of variables
    xvals   = dict(zip(xvars, [mix.Y, mix.Phi]))
    Bdists  = dict(zip(xvars, [mix.B, mix.Bphi]))
    median  = dict(zip(xvars, [mix.medianY, mix.medianPhi]))
    medianB = dict(zip(xvars, [mix.medianB, mix.medianBphi]))
    mean    = dict(zip(xvars, [mix.meanY, mix.meanPhi]))
    
    # Both Y and ER
    for xv in xvars:
        # Plot the distributions
        axs[xv].plot(xvals[xv], Bdists[xv], label = f'N = {N}',color = colors[i], alpha=0.8)
        # Plot the medians
        axs[xv].scatter(median[xv], medianB[xv], s = 17, color = colors[i], alpha = 1, marker = 'v')

idealmean = dict(zip(xvars, [mix.Yfuel, 1.0]))
for xv in xvars:
    axs[xv].legend(fontsize=11)
    axs[xv].set_ylim(0, None) # Freeze the y-range
    axs[xv].plot([idealmean[xv]]*2, list(axs[xv].get_ylim()), 'k:',alpha = 0.3)
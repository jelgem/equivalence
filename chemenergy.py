import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from stoichiometry import *

# = = = = = = = = = = = = = = = = = = = = = = = = = #
# CREATE THE STOICH MIXTURE
# = = = = = = = = = = = = = = = = = = = = = = = = = #
Fuel_dict   = {'C2H2': 1.0}
Ox_dict     = {'O2': 0.210085, 'N2': 0.789915}

# Create the stoichiometric mixture
mix = stoichiometry(Fuel_dict, Ox_dict)

# = = = = = = = = = = = = = = = = = = = = = = = = = #
# RANGES
# = = = = = = = = = = = = = = = = = = = = = = = = = #
Yf  = np.linspace(0, 1, 10001)
Phi = mix.Phifunc(Yf)

Ya = 1.0 - Yf   # undiluted

# = = = = = = = = = = = = = = = = = = = = = = = = = #
# Heat Release
# = = = = = = = = = = = = = = = = = = = = = = = = = #
qn  = (1/mix.Yfuel) * np.minimum(Yf,Ya/mix.s)
qs  = (1 - qn)*np.sign(Phi - 1)

def calculate_f(x, y, C, D, eps=1e-11):
    numerator = D * x / (y+eps) - 1
    denominator = abs(numerator) + eps
    min_part = (x + y / D - abs(x - y / D)) / 2

    f = (1 - min_part / C) * (numerator / denominator)
    return f

#qs = calculate_f(Yf, Ya, mix.Yfuel, mix.s)

xi = (Phi - 1) / (Phi + 1)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#%%                  PLOT: qn and Yf                      #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
fig, axs = plt.subplots(figsize=[6,3.2],layout="constrained")
axs.set_xlabel("Fuel Mass Fraction")
axs.set_ylabel("Chem. Energy (norm)")
axs.grid(alpha=0.2, which='both')
axs.set_xlim(Yf.min(), Yf.max())
axs.set_ylim(0, 1)

axs.plot(Yf, qn, 'k', linewidth =3)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#%%                  PLOT: qn and ER                      #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
fig, axs = plt.subplots(figsize=[6,3.2],layout="constrained")
axs.set_xlabel("Equivalence Ratio")
axs.set_ylabel("Chem. Energy (norm)")
axs.grid(alpha=0.2, which='both')
axs.set_xlim(0, 16)
axs.set_ylim(0, 1)

axs.plot(Phi, qn, 'c', linewidth = 3)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#%%                  PLOT: qn and ER, logscale            #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
fig, axs = plt.subplots(figsize=[6,3.2],layout="constrained")
axs.set_xlabel("Equivalence Ratio")
axs.set_ylabel("Chem. Energy (norm)")
axs.grid(alpha=0.2, which='both')
axs.set_xlim(1e-3, 1e3)
axs.set_ylim(0, 1)

axs.semilogx(Phi, qn, 'c', linewidth = 3)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#%%                  PLOT: qn and qs           #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
fig, axs = plt.subplots(figsize=[6,3.2],layout="constrained")
axs.set_xlabel("Stoichiometry Metric")
axs.set_ylabel("Chem. Energy (norm)")
axs.grid(alpha=0.2, which='both')
axs.set_xlim(qs.min(),qs.max())
axs.set_ylim(0, 1)

axs.plot(qs, qn, 'b', linewidth = 3)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#%%                  PLOT: ER and qs           #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
fig, axs = plt.subplots(figsize=[6,3.2],layout="constrained")
axs.set_xlabel("Stoichiometry Metric")
axs.set_ylabel("Equivalence Ratio")
axs.grid(alpha=0.2, which='both')
axs.set_xlim(xi.min(),xi.max())
axs.set_ylim(1e-3,1e3)

axs.semilogy(qs, Phi, color='b', linewidth = 3, alpha=0.7)

Phi_2 = 2.0 ** np.arange(-9, 10)
qs_2 = calculate_f(Phi_2 / (mix.s + Phi_2), 1 - Phi_2 / (mix.s + Phi_2), mix.Yfuel, mix.s)
axs.scatter(qs_2, Phi_2, marker = 'D', s =15,color='b')

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#%%                  PLOT: ER and EM           #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
fig, axs = plt.subplots(figsize=[6,3.2],layout="constrained")
axs.set_xlabel("Equivalence Metric")
axs.set_ylabel("Equivalence Ratio")
axs.grid(alpha=0.2, which='both')
axs.set_xlim(xi.min(),xi.max())
axs.set_ylim(1e-3,1e3)

axs.semilogy(xi, Phi, color=[0.8, 0.3, 0], linewidth = 3, alpha=0.7)

Phi_2 = 2.0 ** np.arange(-9, 10)
Xi_2 = (Phi_2 - 1) / (Phi_2 + 1)
axs.scatter(Xi_2, Phi_2, marker = 'D', s =15,color=[0.8, 0.3, 0])

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
#%%                  PLOT: qn and  EM           #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = #
fig, axs = plt.subplots(figsize=[6,3.2],layout="constrained")
axs.set_xlabel("Equivalence Metric")
axs.set_ylabel("Chem. Energy (norm)")
axs.grid(alpha=0.2, which='both')
axs.set_xlim(xi.min(),xi.max())
axs.set_ylim(0, 1)

axs.plot(xi, qn, color=[0.8, 0.3, 0], linewidth = 3)

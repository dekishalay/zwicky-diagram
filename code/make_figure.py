""" Make the diagram """

import numpy as np
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from astropy.cosmology import Planck15
# from classes import *

fig, ax = plt.subplots(1,1, figsize=(8,8))

# Stand-in. Using RCF.
search_terms = ['Ia', 'Type II', 'IIb', 'IIn', 'SLSN', 'Ibc', 'Gap', 'ILRT']
for st in search_terms:
    dat = np.loadtxt("%s.txt" %st, dtype=str)
    if np.logical_and(dat.ndim == 1, len(dat) > 0):
        z = dat[1].astype(float)
        m = dat[2].astype(float)
        M = m - Planck15.distmod(z=z).value
        dur = dat[3].astype(float) + dat[4].astype(float)
        ax.scatter(dur, M, label=st)
    if dat.ndim > 1:
        z = dat[:,1].astype(float)
        m = dat[:,2].astype(float)
        M = m - Planck15.distmod(z=z).value
        dur = dat[:,3].astype(float) + dat[:,4].astype(float)
        ax.scatter(dur, M, label=st)
    

ax.set_xlabel("Time Above Half-Max (days)", fontsize=16)
ax.set_ylabel("Peak Luminosity ($M_v$)", fontsize=16)
ax.set_xlim(0.1,100)
ax.set_xscale('log')
ax.set_ylim(-24, -6)
plt.tick_params(axis='both', labelsize=14)
ax.invert_yaxis()

# Axis showing luminosity in erg/s
ax2 = ax.twinx()
# Conversion from M_bol to L_bol, in erg/s
y_f = lambda y_i: 1E7 * 10**((y_i-71.2)/(-2.5))
ymin, ymax = ax.get_ylim()
ax2.set_ylim((y_f(ymin), y_f(ymax)))
ax2.plot([],[])
ax2.set_ylabel("Peak Luminosity (erg/s)", 
        fontsize=16, rotation=270, verticalalignment='bottom')
ax2.tick_params(axis='both', labelsize=14)
ax2.set_yscale('log')

# Display
plt.tight_layout()
plt.show()

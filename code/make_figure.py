""" Make the diagram """

import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from classes import *

fig, ax = plt.subplots(1,1, figsize=(8,8))

ax.set_xlabel("Time Above Half-Max (days)", fontsize=16)
ax.set_ylabel("Peak Bolometric Magnitude", fontsize=16)
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
ax2.set_ylabel("Peak Bolometric Luminosity (erg/s)", 
        fontsize=16, rotation=270, verticalalignment='bottom')
ax2.tick_params(axis='both', labelsize=14)
ax2.set_yscale('log')

# Display
plt.tight_layout()
plt.show()

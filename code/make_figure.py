""" Make the diagram """

import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib.pyplot as plt
plt.rc("font", family="serif")
plt.rc("text", usetex=True)
from astropy.cosmology import Planck15
from classes import snia,relativistic,core_collapse,slsne,fbots,ilrt,lrne,gap
import sys
sys.path.append("/Users/annaho/Dropbox/Projects/Research/Koala/code")
from m_trise import *

fig, ax = plt.subplots(1,1, figsize=(8,6))

snia(ax)
ax.text(6, -20.5, "Thermonuclear\nSupernovae", 
        fontsize=14)#, bbox=dict(facecolor='white', edgecolor='grey'))

core_collapse(ax)
ax.text(
        20, -14.2, "Core-collapse\nSupernovae",
        fontsize=14, color='orange')#,
        #bbox=dict(edgecolor='orange', facecolor='white'))

slsne(ax)
ax.text(20,-23, 'Superluminous\nSupernovae', 
        fontsize=14, color='purple')#, bbox=dict(facecolor='white',
            #edgecolor='purple'))

fbots(ax)
ax.text(2,-20,"FBOTs",fontsize=14,color='k')

ilrt(ax)
ax.text(31, -11, "Intermediate\nLuminosity\nRed Transients",
        horizontalalignment='right', fontsize=14,color='red')

lrne(ax)
ax.text(52, -8, "Luminous\nRed Novae",
        horizontalalignment='right', fontsize=14,color='red')

gap(ax)
ax.text(8, -15, "Ca-rich Gap",
        horizontalalignment='right', fontsize=14,color='black')
ax.text(4, -17, ".Ia Explosions",
        horizontalalignment='right', fontsize=14,color='black')

ax.text(15, -7, "Classical Novae",
        horizontalalignment='left', fontsize=14,color='black')

# Stand-in. Using RCF.
# search_terms = ['Ia', 'Type II', 'IIb', 'IIn', 'SLSN', 'Ibc', 'Gap', 'ILRT']
# for st in search_terms:
#     dat = np.loadtxt("%s.txt" %st, dtype=str)
#     if np.logical_and(dat.ndim == 1, len(dat) > 0):
#         z = dat[1].astype(float)
#         m = dat[2].astype(float)
#         M = m - Planck15.distmod(z=z).value
#         dur = dat[3].astype(float) + dat[4].astype(float)
#         ax.scatter(dur, M, label=st)
#     if dat.ndim > 1:
#         z = dat[:,1].astype(float)
#         m = dat[:,2].astype(float)
#         M = m - Planck15.distmod(z=z).value
#         dur = dat[:,3].astype(float) + dat[:,4].astype(float)
#         ax.scatter(dur, M, label=st)
    
ax.set_xlabel("Time Above Half-Max (rest-frame days)", fontsize=16)
ax.set_xlim(0.01,200)
ax.set_xscale('log')
ax.set_ylim(-5, -30)
ax.set_ylabel("Peak Luminosity ($M_v$)", fontsize=16)

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
ax.tick_params(axis='both', labelsize=14)

# Make an inset for fast transients
# axins = inset_axes(
#         ax, 2, 2, loc=1,
#         bbox_to_anchor=(0.4,0.8),
#         bbox_transform=ax.transAxes)
# axins.set_xlim(1,5)
# axins.set_ylim(-18,-22)
# #mark_inset(ax, axins, loc1=1, loc2=4, fc="none", ec="0.5")
# axins.tick_params(axis='both', labelsize=12)
# plot_ibn(axins)
# fbot(axins)
# kepler(axins)
# ptf09uj(axins)
# cow(axins)
# asu(axins)
# koala(axins)
# gep(axins)

# Display
plt.tight_layout()
#plt.show()
plt.savefig("tau_mv.png", dpi=500)

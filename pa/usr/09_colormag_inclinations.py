# Requires: a file listing inclinations and corresponding visual magnitudes of a star,
#	the same for blue magnitudes;
# Output: a plot of visual magnitude versus color.
# Note: to get the required files run calc_star and calc_spectra for a set of inclinations, then
# 	run filter_spectra once for each filter to get a set of magnitudes.

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc

iodir = '../../' # location of the input/output directory

# the inclinations and the visual magnitudes
incs, V = np.loadtxt(iodir + 'data/vega_V.txt').T

# index of the observed inclination
i = np.searchsorted(incs, 0.08683)

# blue magnitudes
B = np.loadtxt(iodir + 'data/vega_B.txt').T[1]

# colors
color = B - V

# plot
rc('font',**{'family':'serif','serif':['Computer Modern'],'size': 18})
rc('text', usetex=True)
# arrow length
l = 0.01
# dictionary of arrow parameters
d = {'color':'k', 'fill':True, 'linewidth':2, 'length_includes_head':True,\
	'overhang':0.8, 'head_length':l/2, 'head_width':0.03}

max_V = np.max(V)
min_V = np.min(V)
delta_V = max_V - min_V
offset_V = delta_V * 0.1

fig = plt.figure()
ax = plt.axes() 
ax.set_ylim([max_V + offset_V, min_V - offset_V])
ax.scatter(color, V, marker='o', c='b', s=15)
ax.scatter(color[i], V[i], marker='o', c='k', s=15)
ax.set_xlabel('B - V')
ax.set_xticks(0.01 * np.arange(7)[1:])
ax.set_ylabel('V')
ax.arrow(color[i] + 1.2*l, V[i], -l, 0, **d)
fig.tight_layout(pad=0.1)
fig.savefig(iodir + 'vega_colormag.pdf', dpi=200)
plt.close(fig)
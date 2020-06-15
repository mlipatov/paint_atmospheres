# Requires: a file listing inclinations and corresponding visual magnitudes of a star,
#	the same for blue magnitudes;
# Output: a plot of visual magnitude versus color.
# Note: to get the required files run calc_star and calc_spectra for a set of inclinations, then
# 	run filter_spectra once for each filter to get a set of magnitudes.
from pa.lib import util as ut
from pa.lib import star as st
import numpy as np
import pickle
from matplotlib import pyplot as plt
from matplotlib import rc

iodir = '../../' # location of the input/output directory

# unpickle the filtered limb darkening information and get the filter indices
with open(iodir + 'data/limbdark_m01f.pkl', 'rb') as f:
	ld = pickle.load(f)
iV = ld.bands.index('V')
iB = ld.bands.index('B')

# star parameters
omega, luminosity, mass, Req, distance = [0.6151, 40.346, 2.165, 2.815, 2.3694e19] # vega
n_z = 100
incs = np.arccos( np.linspace(1, 0, 50) ) # inclinations, equally spaced in cos(i)
# the observed inclination and record the index
iobs = 0.08683

# create star
star = st.Star(omega, luminosity, mass, Req, distance, n_z, ld) 
# compute its magnitudes at both filters at all inclinations
V = []
B = []
for i in incs:
	mags = star.integrate(i)
	V.append( mags[iV] )
	B.append( mags[iB] )
V = np.array(V)
B = np.array(B)
# colors
color = B - V
# compute the observed values
mags = star.integrate(iobs)
Vobs = mags[iV]
Bobs = mags[iB]
color_obs = Bobs - Vobs

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
ax.scatter(color_obs, Vobs, marker='o', c='k', s=15)
ax.set_xlabel('B - V')
ax.set_xticks(0.01 * np.arange(7)[1:])
ax.set_ylabel('V')
ax.arrow(color_obs + 1.2*l, Vobs, -l, 0, **d)
fig.tight_layout(pad=0.1)
fig.savefig(iodir + 'vega_colormag.pdf', dpi=200)
plt.close(fig)
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

# unpickle the limb darkening information
with open(iodir + 'data/limbdark_m01.pkl', 'rb') as f:
	ld = pickle.load(f)
# the spectrum wavelengths
wl_arr = ld.wl_arr

# visual filter
wlf, T = np.loadtxt(iodir + 'data/Generic_Bessell.V.dat').T
F0 = 3.619e-9 # flux zero point for this filter
fV = ut.Filter(T, wlf, F0)
# blue filter
wlf, T = np.loadtxt(iodir + 'data/Generic_Bessell.B.dat').T
F0 = 6.317e-9 # flux zero point for this filter
fB = ut.Filter(T, wlf, F0)

# star parameters
omega, luminosity, mass, Req = [0.6151, 40.346, 2.165, 2.815] # vega
n_z = 100
distance = 2.3694e19 # cm
incs = np.arccos( np.linspace(1, 0, 50) ) # inclinations, equally spaced in cos(i)
# the observed inclination and record the index
iobs = 0.08683

# create star
star = st.Star(omega, luminosity, mass, Req, n_z, ld) 
# compute its spectrum and flux through the filters at all inclinations
V = []
B = []
for i in incs:
	star_light = star.integrate(i)
	V.append( fV.mag(star_light, wl_arr, distance)[0] )
	B.append( fB.mag(star_light, wl_arr, distance)[0] )
V = np.array(V)
B = np.array(B)
# colors
color = B - V
# compute the observed values
star_light = star.integrate(iobs)
Vobs = fV.mag(star_light, wl_arr, distance)[0]
Bobs = fB.mag(star_light, wl_arr, distance)[0]
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
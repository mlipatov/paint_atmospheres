# Runtime on 2.3GHz MacBook Pro with 8 Gb of RAM: ~7 seconds
# Requires: limb darkening fit information (with fits of intensity w.r.t surface inclination)
#	filter transmission function (from an independent source)
# Output: a plot of one or more transit curves 

from pa.lib import star as st
from pa.lib import util as ut
import numpy as np
import math
import pickle
from matplotlib import pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable

iodir = '../../' # location of the input/output directory

# unpickle the filtered limb darkening information, get the filter index
with open(iodir + 'data/limbdark_m01f.pkl', 'rb') as f:
	ld = pickle.load(f)
iV = ld.bands.index('V')

# star parameters
omega, luminosity, mass, Req, distance = [0.838, 3020, 6.1, 9.16, 1.319e20] # achernar
inclination = np.pi/3
n_z = 100
# transit parameters
b_arr = [-0.3, 0.6]
alpha_arr = [np.pi/3, 0]
radius = 0.01
# transit calculation parameters
n = 200 # number of time points
ns_arr = [7, 7] # number of sight lines per time point

# create star
star = st.Star(omega, luminosity, mass, Req, distance, n_z, ld) 
# compute its magnitude
V = star.integrate(inclination)[iV]
# convert magnitudes to dimension-less flux
star_flux = 10**(V / -2.5)
# # flux
# star_flux = filt.flux(star_light, wl_arr, distance)
# transit flux, normalized
norm_arr = []
# compute normalized transit spectra
for b, alpha, ns in zip(b_arr, alpha_arr, ns_arr):
	# transit info
	tr = st.Transit(b, alpha, radius, n)
	# dimension-less blocked light
	# index 0: location
	blocked_flux = star.transit(inclination, tr, ld, 'V', ns=ns)
	# subtract the transit spectrum from the star spectrum
	transit_flux = star_flux - blocked_flux
	# normalize the fluxes by their maxima
	max_flux = np.max(transit_flux)
	norm_flux = transit_flux / max_flux
	# append to the list of fluxes
	norm_arr.append(norm_flux)

## plot
rc('font',**{'family':'serif','serif':['Computer Modern'],'size': 18})
rc('text', usetex=True)

colors = ['gray', 'black', 'red']
fig = plt.figure()

# main plot
ax = plt.axes() 
ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
for norm_flux, col in zip(norm_arr, colors):
	ax.plot(norm_flux - 1, color=col, linewidth=3)
ax.set(xlabel='Time (arbitrary units)', ylabel=r'$\mathcal{F}/\mathcal{F}_{max} - 1$')

# an inset showing the temperature of the visible surface and the transit line
width, height = fig.get_size_inches() # dimensions of the figure in inches
size_inset = 0.4 # dimension of the inset, in figure heights
# size of the inset in Req
size_req = 2 
# unit conversions
ipr = size_inset * height / size_req # inches per Req
ppi = 72 # points per inch
ppr = ipr * ppi # points per Req
ss = 6 * np.pi * (radius*ppr)**2 # size of scatterplot marker, proportional to the size of the planet
# create the inset axes
ax2 = fig.add_axes([0.32, 0.4, size_inset, size_inset])
# plot the temperature inset
star.plot_temp(ax2, inclination, size_inset * height)
## plot the transit lines
for b, alpha, col in zip(b_arr, alpha_arr, colors):
	# transit
	tr = st.Transit(b, alpha, radius, np.int(n / 10))
	# locations of the planet center
	y, up = tr.locations()
	# scatter plot
	ax2.scatter(y, up, marker='o', c=col, s=ss, alpha=1, zorder=2, clip_on=False)
# save figure
fig.savefig(iodir + "transit.pdf", dpi=200, bbox_inches='tight')
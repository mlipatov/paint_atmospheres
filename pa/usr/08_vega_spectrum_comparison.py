# Requires: two files with synthetic spectra of a star in erg s^-1 Hz^-1 ster^-1
#	a file with the observed spectrum of the star in erg s^-1 Hz^-1 cm^-2
# Output: a plot of all three spectra
# Notes: to get the synthetic spectra, run calc_star and calc_spectra;
#	obtain the observed spectrum from an independent source.

from pa.lib import star as st
from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc

w = 3000 # don't plot above this wavelength
distance = 2.3694e19 # distance to vega in cm

iodir = '../../' # location of the input/output directory

# get the star's parameters
f = open(iodir + 'data/vega/vega0_088418.txt')
for line in f:
	if 'omega' in line:
		omega = np.float( line.split()[-1] )
	elif 'luminosity' in line:
		luminosity = np.float( line.split()[-1] )
	elif 'mass' in line:
		mass = np.float( line.split()[-1] )
	elif 'Req' in line:
		Req = np.float( line.split()[-1] )
	elif 'inclination' in line:
		inclination = np.float( line.split()[-1] )
	if 'z values' in line:
		nz = np.int( line.split()[-1] )
f.close()

# the wavelengths and the synthetic spectrum
wl_syn, I_syn = np.loadtxt(iodir + 'data/vega/vega0_088418.txt').T

m = wl_syn < w
wl_syn = wl_syn[m]
I_syn = I_syn[m] / distance**2 

# a star without limb darkening information for the temperature plot
star = st.Star(omega, luminosity, mass, Req, nz) 

# observed spectrum
wl_obs, I_obs = np.loadtxt(iodir + 'data/vega/vega2.dat').T
m = wl_obs < w
wl_obs = wl_obs[m]
I_obs = I_obs[m]

# spectrum at inclination = pi/2
wl_inc, I_inc = np.loadtxt(iodir + 'data/vega/vega1_570796.txt').T
m = wl_inc < w
wl_inc = wl_inc[m]
I_inc = I_inc[m] / distance**2

## plot
rc('font',**{'family':'serif','serif':['Computer Modern'],'size': 18})
rc('text', usetex=True)

max_I = max(I_syn)
min_I = min(I_syn)
delta_y = max(I_syn) - min(I_syn)
offset_y = delta_y * 0.1

fig = plt.figure()
ax = plt.axes() 
ax.set_ylim([min_I - offset_y, max_I + offset_y])
ax.ticklabel_format(axis='y')
ax.scatter(wl_inc, I_inc, marker='o', c='grey', s=6, alpha=1)
ax.scatter(wl_obs, I_obs, marker='o', c='r', s=3, alpha=1)
ax.scatter(wl_syn, I_syn, marker='o', c='b', s=6, alpha=1)
ax.set_xlabel(r'$\lambda$ (nm)')
ax.set_ylabel(r'$\mathcal{F}_\nu \enspace (\mathrm{erg} \; \mathrm{s}^{-1} \; \mathrm{Hz}^{-1} \; \mathrm{cm}^{-2})$')

## Add an inset showing the temperature of the visible surface
height_inset = 0.2 # height of the inset, in figure heights
width, height = fig.get_size_inches() # size of the current figure in inches
# the width and height of the inset should be the same in inches
width_inset = height_inset * width / height # width of the inset in figure widths
# create the axes with the temperature inset
ax2 = fig.add_axes([0.55, 0.6, width_inset, height_inset])
divider = make_axes_locatable(ax2)
cax = divider.new_vertical(size="5%", pad=0.15, pack_start=True)
fig.add_axes(cax)

# plot the temperature inset
star.plot_temp(ax2, inclination, height_inset * height, cax=cax)

fig.savefig(iodir + 'vega_spectrum_comparison.pdf', dpi=200, bbox_inches='tight')
plt.close(fig)

## plot the Balmer jump
wl1 = 200
wl2 = 550
i1_syn = np.searchsorted(wl_syn, wl1)
i2_syn = np.searchsorted(wl_syn, wl2)
wl_syn = wl_syn[i1_syn:i2_syn]
I_syn = I_syn[i1_syn:i2_syn]
i1_obs = np.searchsorted(wl_obs, wl1)
i2_obs = np.searchsorted(wl_obs, wl2)
I_obs = I_obs[i1_obs:i2_obs]
wl_obs = wl_obs[i1_obs:i2_obs]
i1_inc = np.searchsorted(wl_inc, wl1)
i2_inc = np.searchsorted(wl_inc, wl2)
I_inc = I_inc[i1_inc:i2_inc]
wl_inc = wl_inc[i1_inc:i2_inc]

fig = plt.figure()
ax = plt.axes() 
ax.set_ylim([min_I - offset_y, max_I + offset_y])
ax.ticklabel_format(axis='y')
ax.scatter(wl_inc, I_inc, marker='o', c='grey', s=6, alpha=1)
ax.scatter(wl_obs, I_obs, marker='o', c='r', s=3, alpha=1)
ax.scatter(wl_syn, I_syn, marker='o', c='b', s=6, alpha=1)
ax.set_xlabel(r'$\lambda$ (nm)')
ax.set_ylabel(r'$\mathcal{F}_\nu \enspace (\mathrm{erg} \; \mathrm{s}^{-1} \; \mathrm{Hz}^{-1} \; \mathrm{cm}^{-2})$')

fig.savefig(iodir + 'vega_balmer.pdf', dpi=200, bbox_inches='tight')
plt.close(fig)
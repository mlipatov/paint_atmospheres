# Requires: a limb darkening information file (with fits of intensity w.r.t surface inclination)
# Outputs: a plot that compares different gravity interpolation schemes

import math
import numpy as np
import pickle
import pa.lib.star as st
import pa.lib.util as ut
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import rc

iodir = '../../' # location of the input/output directory

## unpickle the limb darkening information
with open(iodir + 'data/limbdark_m01.pkl', 'rb') as f:
	ld = pickle.load(f)

# constant star parameters
omega, inclination, n = [0, 0, 100]
# temperatures of the stars in Kelvin
temp = np.array([6000, 9000, 12000])
## put the stars on a fictional main sequence where log g = 4.5,
## there is a star with the Sun's mass and luminosity,
## luminosity relates to temperature via a power law,
## and luminosity is proportional to mass**3.5
g = 10**4.5 # gravity in cgs units
logg = 4.5 # log of same
# power and multiplicative constant in the luminosity-temperature dependence,
# with luminosity in solar luminosities and  temperature in Kelvins
m = 28./5 
k = ( 4 * math.pi * ut.sigma * ut.G * ut.Msun / (g * ut.Lsun) )**(7./5)
# the stars' luminosities in solar luminosities
L = k * temp**m
# the stars' masses in solar masses
M = L**(2./7)
# radii of the stars in solar radii
R = L**(1./7) * np.sqrt(ut.G * ut.Msun/ g) / ut.Rsun 

print('Non-rotating stars')
print('Temperatures in Kelvin ' + str(temp))
print('Luminosities in solar luminosities ' + str(L))
print('Masses in solar masses ' + str(M))
print('Radii in solar radii ' + str(R))

# stars with linear temperature interpolaion and full limb darkening information
stars = [st.Star(omega, l, m, r, ut.D10, n, ld=ld) for l, m, r in zip(L, M, R)]
# the light from such stars
light = np.array([s.integrate(inclination) for s in stars])

# remove the limb darkening information at the gravity of the stars
ind_g = np.searchsorted(ld.g, logg)
print('Gravity removed from intensity information: ' + str(ld.g[ind_g]) + \
	'. \nThese should be the same as the star gravity.')
print('Gravities between which we are interpolating: ' + str(ld.g[ind_g - 1]) +\
	' and ' + str(ld.g[ind_g + 1]) )
ld.fit_params = np.delete(ld.fit_params, ind_g, axis=1)
ld.g = np.delete(ld.g, ind_g)

# stars with missing limb darkening information and log gravity interpolation
stars = [st.Star(omega, l, m, r, ut.D10, n, ld=ld) for l, m, r in zip(L, M, R)]
# the light from such stars
light_log = np.array([s.integrate(inclination) for s in stars])

# stars with missing limb darkening information and linear gravity interpolation
stars = [st.Star(omega, l, m, r, ut.D10, n, ld=ld, g_method='lin') for l, m, r in zip(L, M, R)]
# the light from such stars
light_lin = np.array([s.integrate(inclination) for s in stars])

### calculate difference spectra in erg/s/ster/Hz

## proportional spectra
# keep only the light at wavelengths where the precise spectrum is non-zero in all the stars
mask = np.all(light != 0, axis = 0)
diff_log = np.abs(light_log[:, mask] / light[:, mask] - 1)
diff_lin = np.abs(light_lin[:, mask] / light[:, mask] - 1)
light = light[:, mask]

# wavelength array
lam = ld.lam[mask]
# cutoff wavelength in nm for the plot
wl = 2000
ind_wl = max(np.searchsorted(lam, wl) - 1, 1)
ind_ld_wl = max(np.searchsorted(ld.lam, wl) - 1, 1)
# truncate the spectra and the wavelength array at the cutoff wavelength
lam = lam[:ind_wl]
diff_log = diff_log[:, :ind_wl]
diff_lin = diff_lin[:, :ind_wl]
light = light[:, :ind_wl]
diff = np.concatenate((diff_lin, diff_log), axis=1)

# plot results
plt.rcParams.update({'font.size': 18})
rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('text', usetex=True)

ofile = iodir + 'error_g.pdf'
f, axarr = plt.subplots(3, sharex=True)
T_y = [0.17, 0.7, 0.7]
T_x = [0.57, 0.7, 0.7]
min_x = np.min(lam) * 0.9
max_x = np.max(lam) * 1.1
max_y = 6.3
min_y = 6.3e-7
if max_y < np.max(diff):
	print('maximum y on the plot,' + str(max_y) + \
		', should be above maximum y value in the data, ' + str(np.max(diff)))

# instantiate a second array of axes objects that share the same x-axis
axarr2 = [ axarr[i].twinx() for i in range(3) ]
for i in range(3):
	# spectrum axes
	axarr2[i].plot(lam, light[i][:ind_ld_wl], color='grey', linewidth=2, alpha=0.5)
	axarr2[i].tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
	axarr2[i].grid(False)
	axarr2[i].set_zorder(1)

	axarr[i].set_xscale('log')
	axarr[i].set_xlim(left=min_x, right=max_x)
	axarr[i].scatter(lam, diff_lin[i], marker='o', c='b', s=6, alpha=0.5, edgecolors="none")
	axarr[i].scatter(lam, diff_log[i], marker='o', c='g', s=6, alpha=1, edgecolors="none")
	axarr[i].text(T_x[i], T_y[i], 'T = ' + str(temp[i]) + ' K', transform=axarr[i].transAxes)
	axarr[i].tick_params(axis="both")
	axarr[i].set_yscale('log')
	axarr[i].set_ylim(min_y/5, max_y*5)
	axarr[i].set_yticks([1e-6, 1e-3, 1])
	axarr[i].set_zorder(2)
	axarr[i].patch.set_alpha(0)

## add a big axes, hide its frame
ax = f.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
ax.grid(False)
# display the x and y axis labels on the big axes
ax.set_xlabel(r'$\lambda$ (nm)')
ax.set_ylabel(r'$|\delta \mathcal{F}_\nu \,/\, \mathcal{F}_\nu|$', labelpad=20)
## save figure
f.savefig(ofile, dpi=200, bbox_inches='tight')
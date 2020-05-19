# Requires: a file with limb darkening fit information that includes the intensity values
#	a .npy file with locations and values of the lowest I(mu) / I_1 
# 		and the lowest I'(mu) / I_1 for every (T, g, lambda)
# Outputs: a plot with I(mu) fits for the (T, g, lambda) with maximum and median deviations
#	information about these fits
#	information about (T, g, lambda) with the lowest I(mu) / I_1 and the lowest I'(mu) / I_1
#	a plor of the corresponding fits
# Notes: to obtain the required files, run calc_limbdark with the -s option, plus 03a_Imu_fits_min.py

import numpy as np
import pickle
import pa.lib.limbdark as limbdark
import pa.lib.fit as ft
import matplotlib.pyplot as plt
from matplotlib import rc

# given an array of mu, returns an array of p-functions for each mu:
# 0: mu 
# 1: p_i
def pi(mu):
	return np.transpose(np.array([ np.ones_like(mu), mu, mu**2, mu**3, mu**4]))

# points for plots
# inputs: a_ij, I, index of location in I, indices of m_j, p_ij
def Iplot(a, I, ii, im, p):
	aa = a[ tuple(ii) ][im, :]
	Ifit = np.sum(aa * p, axis=-1)
	Igrid = I[tuple(ii)]
	Ifit = Ifit / Igrid[-1]
	Igrid = Igrid / Igrid[-1]
	return Ifit, Igrid

iodir = '../../' # location of the input/output directory

# unpickle the limb darkening information
with open(iodir + 'data/limbdark_m01.pkl', 'rb') as f:
	ld = pickle.load(f)

wl = ld.wl_arr  # 1221 wavelength
g = ld.g_arr	# 11 gravity
T = ld.temp_arr # 61 temperature
bounds = ld.bounds
I = ld.I_arr # (1221, 17, 11, 61) = (wavelength, mu, gravity, temperature)
a = ld.fit_params # (61, 11, 1221, 15) = (temperature, gravity, wavelength, parameter index)
sh = a.shape

# unpickle the minima of intensity fits and their derivatives
Im = np.load(iodir + 'Imin.npy')
# location where the intensity fit is lowest
imin = np.unravel_index(np.argmin(Im[..., 1]), Im[..., 1].shape)
# location where the derivative is lowest
ider = np.unravel_index(np.argmin(Im[..., 3]), Im[..., 3].shape)

# set the mu partition
ft.set_muB(bounds)
# a_ij coefficients: (61, 11, 1221, 3, 5) = (temperature, gravity, wavelength, interval, function)
a = a.reshape( (sh[0], sh[1], sh[2], ft.m, ft.n) ) 

# 17 mu
mu = ft.mu_arr 
# functions at 17 mu: 2D
p = pi(mu)
# intervals where 17 mu are found
i = np.searchsorted(ft.muB_arr, mu, side='right') - 1 
# at each location, in the corresponding interval, for each wavelength, 
# sum up the product of fit parameters and functions of mu
# (61, 11, 1221, 17, 5) -> (61, 11, 1221, 17) = (temperature, gravity, wavelength, mu)
Ifit = np.sum(a[ ..., i, : ] * p[ np.newaxis, np.newaxis, np.newaxis, :, : ], axis=-1)

# permute the dimensions of CK04 intensity grid: (61, 11, 1221, 17) = (temperature, gravity, wavelength, mu)
I = np.transpose(I, axes=[3, 2, 0, 1])
# intensity at mu = 1 from the grid
I1 = I[...,-1][...,np.newaxis]
# relative error in intensity
Ierr = np.abs( (Ifit - I) / (I1 + 1e-300) )

# maximum and median of maximum differences for (temperature, gravity, wavelength) triples
Ierr_max = np.amax(Ierr, axis=-1)

maxerr = np.nanmax(Ierr_max)
imax = np.unravel_index(np.nanargmax(Ierr_max), Ierr_max.shape)

n = np.count_nonzero(~np.isnan(Ierr_max)) # number of non-nan entries
mederr = np.nanmedian(Ierr_max)
dist = np.abs(Ierr_max - mederr)
imed = np.array(np.unravel_index(np.argsort(dist, axis=None), dist.shape))

## points to plot

# mu values to plot
x = np.linspace(0, 1, 100)
# p values to plot
px = pi(x)
# intervals
i = np.searchsorted(ft.muB_arr, x, side='right') - 1 
# intensities
Imaxfit, Imax = Iplot(a, I, imax, i, px)
Imedfit, Imed = Iplot(a, I, imed[:, 0], i, px)
Iminfit, Imin = Iplot(a, I, imin, i, px)
Iderfit, Ider = Iplot(a, I, ider, i, px)

# print relevant values
print('Maximum error is ' + str(maxerr))
print('This is realized at ')
print('T = ' + str(T[imax[0]]))
print('g = ' + str(g[imax[1]])) 
print('lambda = ' + str(wl[imax[2]]))
print()
print('Median error is ' + str(mederr))
print('There are ' + str(n) + ' locations.')
print('One of the locations whose error is closest to this is ')
print('T = ' + str(T[imed[0][0]]))
print('g = ' + str(g[imed[1][0]])) 
print('lambda = ' + str(wl[imed[2][0]]))
print('The error at this location is ' + str(Ierr_max[tuple(imed[:,0])]))
print()
print('Smallest relative value of intensity is ' + str(Im[imin][1]))
print('The location is ')
print('T = ' + str(T[imin[0]]))
print('g = ' + str(g[imin[1]])) 
print('lambda = ' + str(wl[imin[2]]))
print('mu = ' + str(Im[imin][0]))
print()
print('Smallest relative value of the derivative is ' + str(Im[ider][3]))
print('The location is ')
print('T = ' + str(T[ider[0]]))
print('g = ' + str(g[ider[1]])) 
print('lambda = ' + str(wl[ider[2]]))
print('mu = ' + str(Im[ider][2]))
print()

# interval bounds for plotting
boundsp = np.concatenate( ([0], bounds, [1]) )

plt.rcParams.update({'font.size': 18})
rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('text', usetex=True)

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=False,
                                           gridspec_kw = {'height_ratios':[1, 0.34]},
                                           figsize=(6, 6))

ax1.plot(x, Imaxfit, color='b', zorder=1)
ax1.scatter(mu, Imax, color='k', s=20, zorder=2)
ax1.set_ylabel(r'$I_\nu / I_1$')
ax1.set_ylim((-0.05, 1.05))

ax1.plot(x, Imedfit, color='g', zorder=1)
ax1.scatter(mu, Imed, color='k', s=20, facecolors='w', zorder=2, alpha=1)

for b in boundsp:
	ax1.axvline(x=b, color='k', alpha=0.5, lw=0.5)

ax2.plot(mu, Ierr[tuple(imax)], color='b', linewidth=1, linestyle='--', zorder=1)
ax2.scatter(mu, Ierr[tuple(imax)], color='b', s=20, zorder=2)

ax2.plot(mu, Ierr[tuple(imed[:,0])], color='g', linewidth=1, linestyle='--', zorder=1)
ax2.scatter(mu, Ierr[tuple(imed[:,0])], color='g', s=20, facecolors='w', zorder=2, alpha=1)

ax2.margins(y=0.2)
ax2.set_yscale('log')
ax2.set_yticks([1e-9, 1e-6, 1e-3])
ax2.set_ylabel(r'$\left|\delta I_\nu / I_1\right|$', labelpad=10)

ax2.set_xlabel(r'$\mu$')
for b in boundsp:
	ax2.axvline(x=b, color='k', alpha=0.5, lw=0.5)

fig.subplots_adjust(hspace=0.1)
fig.tight_layout()
fig.savefig(iodir + 'Imu.pdf')


# plot the locations with the most negative I and the most negative derivative of I
fig, ax = plt.subplots(1, 1, figsize=(6, 6))

ax.plot(x, Iminfit, color='b', zorder=1)
ax.scatter(mu, Imin, color='k', s=20, zorder=2)
ax.set_ylabel(r'$I_\nu / I_1$')
ax.set_ylim((-0.05, 1.05))

ax.plot(x, Iderfit, color='g', zorder=1)
ax.scatter(mu, Ider, color='k', s=20, facecolors='w', zorder=2, alpha=1)

for b in boundsp:
	ax.axvline(x=b, color='k', alpha=0.5, lw=0.5)

fig.tight_layout()
fig.savefig(iodir + 'Imu_min.pdf')
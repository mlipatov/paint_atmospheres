# Runtime on 2.3GHz MacBook Pro with 8 Gb of RAM: ~10 minutes
# Requires: a file with limb darkening fit information that includes the intensity values
# Outputs: a .npy file with locations and values of the lowest I(mu) / I_1 
# 	and the lowest I'(mu) / I_1 for every (T, g, lambda)
# Notes: to obtain the required file, run calc_limbdark with the -s option

import numpy as np
import pickle
import pa.lib.limbdark as limbdark
import pa.lib.fit as ft
import matplotlib.pyplot as plt
from matplotlib import rc

# real roots of a polynomial within bounds
# input: coefficients, lower bound, upper bound
def rts(a, b1, b2):
	r = np.roots(a)
	r = np.real(r[ np.isreal(r) ])
	r = r[ (b1 <= r) & (r <= b2) ]
	return r

# minimum of a polynomial: (location, value)
# inputs: polynomial's coefficients, locations of local suprema and boundary points
def minim(a, s, b1, b2):
	pts = np.concatenate( (np.array([b1, b2]), s) )
	val = np.polyval(a, pts)
	i = np.argmin(val)
	return pts[i], val[i]

iodir = '../../' # location of the data directory

# unpickle the limb darkening information
with open(iodir + 'data/limbdark_m01.pkl', 'rb') as f:
	ld = pickle.load(f)

wl = ld.lam  # 1221 wavelength
g = ld.g	# 11 gravity
T = ld.T[0:1] # 61 temperature
bounds = ld.bounds
I = ld.I[..., 0:1] # (1221, 17, 11, 61) = (wavelength, mu, gravity, temperature)
a = ld.fit_params[0:1, ...] # (61, 11, 1221, 15) = (temperature, gravity, wavelength, parameter index)
sh = a.shape

# (temperature, gravity, wavelength, 4 values)
# Last dimension gives [smallest-value-mu, smallest-value, smallest_der_mu, smallest-der]
Imin = np.full( (sh[0], sh[1], sh[2], 4), np.nan )

# set the mu partition in the Fit class
ft.set_muB(bounds)
# create a bounds array that includes 0 and 1
bds = np.concatenate( (np.array([0]), bounds, np.array([1])) )
# (61, 11, 1221, 3, 5) = (temperature, gravity, wavelength, interval, function)
a = a.reshape( (sh[0], sh[1], sh[2], ft.m, ft.n) ) 

# permute the dimensions of the stored intensity array to (temperature, gravity, wavelength, mu)
I = np.transpose(I, axes=[3, 2, 0, 1])

print('Computing minima of intensity fits and their derivatives')
# check if any of the fits are negative or have negative derivatives
for iT in range(sh[0]):
	print('T = ' + str(T[iT]))
	for ig in range(sh[1]):
		for iw in range(sh[2]):
			I1 = I[iT, ig, iw,-1] # intensity at mu = 1 from the grid
			Im = [np.nan, np.inf] # location and value of minimum I
			Ipm = [np.nan, np.inf] # location and value of minimum I prime
			for ii in range(ft.m):
				# extract the coefficients on this interval
				aa = a[iT, ig, iw, ii]
				if ~np.isnan(aa[0]): # if the coefficients exist
					# the bounds on this interval
					b1, b2 = bds[ii], bds[ii + 1]
					# coefficients of the derivative
					aap = ( np.arange(5) * aa )[1:]
					# coefficients of the second derivative
					aapp = ( np.arange(4) * aap )[1:]
					# reverse the coefficients for numpy
					aa = np.flip(aa)
					aap = np.flip(aap)
					aapp = np.flip(aapp)
					# real roots of the polynomial within the bounds
					rI = rts(aa, b1, b2)
					# real roots of the derivative within the bounds
					rIp = rts(aap, b1, b2)
					# real roots of the second derivative within the bounds
					rIpp = rts(aapp, b1, b2)
					# location and value of minimum I
					pI, vI = minim(aa, rIp, b1, b2)
					# normalize
					vI /= I1 + 1e-300
					# update
					if vI < Im[1]: Im = [pI, vI] 
					# location and value of minimum I prime
					pIp, vIp = minim(aap, rIpp, b1, b2)
					# normalize
					vIp /= I1 + 1e-300
					# update
					if vIp < Ipm[1]: Ipm = [pIp, vIp]
			Imin[iT, ig, iw] = np.array(Im + Ipm)

np.save(iodir + 'Imin.npy', Imin)


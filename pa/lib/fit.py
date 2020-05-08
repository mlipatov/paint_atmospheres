## This module contains the class Fit as well as the constants that this class uses
# The range of mu, [0, 1], is split into several intervals. 
# The model is linear in the parameters of the fit on each interval.
# The constraints on the function: the intensity at mu = 1 is set, 
#		the function is continuous everywhere in its domain 

import numpy as np
import matplotlib.pyplot as plt

## functions on a given interval
# the non-zero functions on a given interval
f = [lambda _: 1, lambda x: x, lambda x: x**2, lambda x: x**3, lambda x: x**4]
# the number of non-zero functions on a given interval
n = len(f)
# zero functions on a given interval
fz = [lambda _: 0] * n

## the values of mu in Castelli and Kurucz 2004, reversed
mu_arr = np.array([0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.])

# the boundaries for the intervals of the fit: [0, x0], [x0, x1], ... [x(n), 1]
muB_arr = np.array([])
# the number of intervals
m = 0	
# mu array split into a list of arrays according to the boundaries
mu_arr_split = np.array([])
# all the functions on each interval
F = []
# the model for the given complete set of mu values
x = []

# sets the array of boundary values between mu intervals
# then splits the array of mu values accordingly, then computes all the functions on all the intervals
# and the linear model (without coefficients) for all the given mu values
def set_muB(b):
	global muB_arr, m, mu_arr_split, F, x
	# set the array of boundaries, including the left edge of the lowest interval
	muB_arr = np.array([0] + b)
	# set the number of intervals
	m = len(b) + 1
	# split the array of Kurucz's mu values according to the boundaries between intervals
	# don't use the zero that designates the left edge of the lowest interval
	mu_arr_split = np.split(mu_arr, np.searchsorted(mu_arr, muB_arr[1:]))
	# include the boundaries twice, both in the interval on the left and the interval on the right
	for i in range(m - 1):
		mu_arr_split[i] = np.append(mu_arr_split[i], mu_arr_split[i + 1][0])
	# compute all the functions on each interval
	for i in range(m):
		F.append(fz * i + f + fz * (m - i - 1))
	## compute the model for the given complete set of mu values
	# stack the column obtained from applying to each value of mu all the functions 
	# on that value's interval
	xl = []
	for i in range(m): # interval
		Fi = F[i] # all functions on this interval
		for mu in mu_arr_split[i]: # mu on this interval
			xj = [func(mu) for func in Fi]
			xl.append(np.array(xj))
	x = np.array(xl)

# initialize an object of this class with an array of intensities corresponding to the array 
# of mu values in this module and the wavelength, log g and temperature of this fit;
# compute the fit
def fit(I_arr):
	# for the linear algebra fit,
	# duplicate the intensity values at the boundaries of mu intervals
	I_arr_dup = np.copy(I_arr)
	ind = len(I_arr_dup) - 1
	for i in reversed(range(m - 1)):
		sub_ind = len(mu_arr_split[i + 1]) - 1
		ind -= sub_ind
		I_arr_dup = np.insert(I_arr_dup, ind, I_arr_dup[ind])
	# fit
	alpha = np.linalg.lstsq(x, I_arr_dup, rcond=None)[0]
	return alpha

# output: a 2D array of intensity values (location x wavelength)
# units: ergs/cm2(surf)/s/hz/ster
# inputs: a 1D array of mu values (corresponding to locations)
# 	a 3D array of fit parameters (location x wavelength x parameter index)
# note: mu has to be in [0, 1]
def I(mu, p):
	# 2D array of functions evaluated at different locations
	# 0: location
	# 1: function
	f_mu = np.transpose( np.array([ np.ones_like(mu), mu, mu**2, mu**3, mu**4]) )
	# 1D array of lowest interval indices where each mu value is found,
	# each interval of the form [mu1, mu2]
	# 0: location
	i = np.searchsorted(muB_arr, mu, side='right') - 1 
	# reshape the parameter array to distinguish between functions on different intervals
	# sh = p.shape
	# # 0: location
	# # 1: wavelength
	# # 2: interval
	# # 3: function
	# params = p.reshape( (sh[0], sh[1], cls.m, n) )
	params = p.reshape( (m, n) )
	# at each location, in the corresponding interval, for each wavelength, 
	# sum up the product of fit parameters and functions of mu
	output = np.sum(f_mu[ :, np.newaxis, : ] * params[ np.arange(sh[0]), :, i, : ], axis=2)
	return output

# For each function of the fit on each mu interval of the fit, 
# 	computes twice the integral of the function on the intersection of the interval 
# 	and the integration boundaries on mu
# Inputs: coefficients a and b in the relationship between mu and phi: mu = a * cos(phi) + b,
#	a boolean saying whether the integration is *not* for all values of phi in [-pi, pi]
# Output: an array of integrated functions, one for each function on each mu interval
def integrate(belowZ1, a, b):
	# phi in terms of mu
	def phi(mu, a, b):
		cosn = (mu - b) / a
		return np.arccos(cosn)
	# indefinite integrals of the non-zero fit functions as functions of phi,
	# w.r.t. phi on interval i, given the evaluated elliptic integral; 
	# this function should be modified if the fit functions change
	def intgrt(phi, a, b):
		# cosine phi
		cosn = np.cos(phi)
		cos2 = np.cos(2*phi)
		sine = np.sin(phi)
		sin2 = np.sin(2*phi)
		sin3 = np.sin(3*phi)
		sin4 = np.sin(4*phi)
		sin5 = np.sin(5*phi)
		# integral of the polynomial
		integr = [
			b*phi + a*sine,
			(a**2/2. + b**2)*phi + (2*a*b + (a**2*cosn)/2.)*sine,
			((3*a**2*b)/2. + b**3)*phi + ((5*a**3)/6. + 3*a*b**2 + \
				(3*a**2*b*cosn)/2. + (a**3*cos2)/6.)*sine,
			((3*a**4)/8. + 3*a**2*b**2 + b**4)*phi + (3*a**3*b + 4*a*b**3)*sine + \
				(a**4/4. + (3*a**2*b**2)/2.)*sin2 + (a**3*b*sin3)/3. + (a**4*sin4)/32.,
			((15*a**4*b)/8. + 5*a**2*b**3 + b**5)*phi + \
				((5*a**5)/8. + (15*a**3*b**2)/2. + 5*a*b**4)*sine + \
				((5*a**4*b)/4. + (5*a**2*b**3)/2.)*sin2 + \
				(5*a**5*sin3)/48. + (5*a**3*b**2*sin3)/6. + \
				(5*a**4*b*sin4)/32. + (a**5*sin5)/80.
		]
		return np.transpose( np.array(integr) )
	# upper bound on integration w.r.t. mu; mu will decrease
	mu1 = a + b
	# indices of the mu intervals that contain the upper mu integration bound;
	i = np.searchsorted(muB_arr, mu1, side='right') - 1
	# initialize (to zeros) the total integral (divided by 2) for each function on each interval
	# index 0: z
	# index 1: interval
	# index 2: function
	result = np.zeros( a.shape + (m, n) )
	# a mask telling us where a != 0 in the result array
	anz = ~np.array(a == 0)
	## record the integrals at the locations where a = 0
	# if we are looking at the star pole-on, mu doesn't change as phi changes,
	# so we integrate from zero to pi in the current (and only) mu interval
	result_azero = intgrt(np.pi, a[~anz], b[~anz]) - intgrt(0, a[~anz], b[~anz])
	# set the results at a = 0
	result[ ~anz, i[~anz], : ] = result_azero
	### everything below is done for the locations where a != 0
	## when a != 0, mu changes as phi changes, 
	## so we will keep track of whether we enter different mu intervals.
	# lower bound on integration w.r.t. phi; phi will increase
	ph = np.zeros_like(a)
	# the lower bound on integration w.r.t. mu, between 0 and 1;
	# when integrating over phi between -pi and pi, this lower bound is some number above zero,
	# otherwise it is zero
	mu0 = np.zeros_like(a) 
	mu0[ ~belowZ1 & anz] = b[ ~belowZ1 & anz ] - a[ ~belowZ1 & anz ]
	# the corresponding upper bound on integration w.r.t. phi, between 0 and pi
	phi_mu0 = np.full_like(a, np.pi)
	phi_mu0[  belowZ1 & anz ] = phi( 0, a[ belowZ1 & anz ], b[ belowZ1 & anz ] )
	# lower endpoint of the current mu interval
	mu_int = muB_arr[ i ]
	# initialize a mask that filters for the locations
	# where the lower endpoint of the current mu interval 
	# is higher than the lower mu integration bound
	mask = (mu_int > mu0) & anz
	# where the lower endpoints of current mu intervals are higher 
	# than the lower integration bound, integrate to the values of phi corresponding
	# to these lower endpoints
	while np.any(mask):
		# phi corresponding to the lower endpoint of the current mu intervals
		phi_int = phi(mu_int[mask], a[mask], b[mask])
		# compute the integrals on these intervals, from the current values of phi 
		# to those corresponding to the lower endpoints of the current mu intervals
		result[ mask, i[mask], : ] += \
			intgrt(phi_int, a[mask], b[mask]) - intgrt(ph[mask], a[mask], b[mask])
		# update the current values of phi, only leave the values
		# where the lower endpoints of the current mu intervals are higher than the
		# lower mu integration bounds
		ph[mask] = np.copy(phi_int)
		# move one interval to the left
		i[mask] -= 1
		# update the lower endpoints of the current intervals
		mu_int[mask] = muB_arr[ i[mask] ]
		# update the mask telling us where the lower endpoint of the current mu interval 
		# is higher than the lower mu integration bound;
		mask &= (mu_int > mu0)
	# compute the integrals on the remaining parts of the last intervals
	# for all the locations with non-zero a
	result[ anz, i[anz], : ] += intgrt(phi_mu0[anz], a[anz], b[anz]) - intgrt(ph[anz], a[anz], b[anz])
	# flatten the interval and function dimensions into a single dimension
	sh = result.shape
	result = result.reshape( sh[0], sh[1] * sh[2] )
	# return twice the computed integrals
	return 2 * result
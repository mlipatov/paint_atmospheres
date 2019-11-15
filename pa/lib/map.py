from pa.lib import fit as ft
from pa.lib import util as ut
import numpy as np
import math
import sys

np01 = np.array([0, 1]) # a numpy array containing 0 and 1

class Map:
	""" At initialization with a surface, a step in z values and limb darkening information, 
	computes a number of equally spaced discrete z values on (-1, 1), as well as
	wavelength-dependent parameters of the intensity fits and area elements at these z values. 
	Gravity and temperature calculations are based on Espinosa Lara 2011 (EL)."""

	# initialize with the surface shape of the star, the number of z values to use for integration, 
	# the constants needed for temperature and gravity unit conversions, and limb darkening information
	def __init__(self, surf, nz, add_logg, mult_temp, ld, temp_method):

		## initialize the surface
		self.surface = surf # surface shape information
		omega = surf.omega

		# store gravity and temperature modifiers, 
		# as well as the temperature interpolation method
		self.add_logg = add_logg
		self.mult_temp = mult_temp
		self.temp_method = temp_method

		## initialize the array of z values for integration
		## do not use the boundary values
		self.z_arr = np.linspace(-1, 1, nz + 2) # include z = -1, 1
		# store the spacing between the z values
		self.dz = self.z_arr[1] - self.z_arr[0]
		# store the number of z values
		self.nz = nz

		# compute the cylindrical coordinate r
		r_arr = self.surface.R( self.z_arr )

		# compute the area elements
		self.A_arr = surf.A( self.z_arr )

		# compute and store the temperatures and fit parameters
		self.temp_arr, self.params_arr = self.Tp( self.z_arr, r_arr, ld )


	# output: temperatures and intensity fit parameters for a set of locations
	# input: z values, r values and limb darkening information
	def Tp(self, z, r, ld):		
		# compute the spherical coordinate rho
		rho_arr = self.surface.rho( r, z )
		r01 = self.surface.R( np01 ) # a numpy array containing r at z = 0 and r at z = +/-1
		rho0, rho1 = self.surface.rho( r01, np01 ) # rho at z = 0 and +/-1
		
		## compute the effective gravitational acceleration in units of G M / Re**2 
		geff_arr = self.geff(rho_arr, r)
		# convert to log10(gravity in cm / s**2) and to a less memory-intensive data type
		logg_arr = self.add_logg + np.log10(geff_arr)
		logg_arr = logg_arr.astype(np.float32)
		
		## compute the effective temperature in units of ( L / (4 pi sigma Re**2) )**(1/4), 
		## as in EL eqn 31, for each z; then convert it to Kelvin
		[F_arr, F0, F1] = self.F(z, rho_arr, rho0, rho1) # temperature correction values
		# temperatures
		temp_arr = self.mult_temp * ( geff_arr * F_arr )**(1./4)
		# convert temperature values to a less memory-intensive data type
		temp_arr = temp_arr.astype(np.float32)

		# compute the interpolated values of limb darkening fit parameters
		# if given limb darkening information
		if ld is not None:
			params_arr = self.interp(logg_arr, temp_arr, ld, self.temp_method)
		else:
			params_arr = None

		# return
		return [temp_arr, params_arr]
	
	# output: effective gravitational acceleration in units of G M / Re**2 
	# 	as in equations 7 and 31 of EL, for each z
	# inputs: arrays of rho and r coordinates
	def geff(self, rho_arr, r_arr):
		omega = self.surface.omega
		r_arr_sq = np.power(r_arr, 2) # r squared
		geff_arr = np.sqrt( rho_arr**-4 + omega**4 * r_arr_sq - \
			2 * omega**2 * r_arr_sq * rho_arr**-3 )
		return geff_arr

	# output: temperature correction factors (EL equation 26)
	# 	at every value of z; runs the Newton's method algorithm that acts on the entire array 
	# 	of z values that aren't too close to 0 at its every step; for values that 
	# 	are close to 0, returns the temperature correction expected at z = 0
	def F(self, z_arr, rho_arr, rho0, rho1):
		# output: smallest value of x for which to compute 
		#	using full Newton's method step function; a.k.a. x_b
		# inputs: rotational velocity of the star,
		#	temperature correction at x = 0,
		#	order-one factor of proportionality between the error 
		#		in full step function and that in the series expansion
		#	resolution of floating point numbers
		def X1(omega, A, r):
			# a factor that's a little less than 1
			B = 0.9
			# omega below which the approximation yields values greater than 1
			omega_lim = 0.9 * (2./(85*A))**(1./4) * 3**(1./2) * r**(1./4)
			if omega < omega_lim:
				output = 1
			else: # otherwise, evaluate the estimate of x_b
				output = r**(1./6) * (2./(85*A))**(1/6) * \
					( 3**(1./3) * omega**(-2./3) - \
					  3**(-2./3) * (199./255) * omega**(4./3) - \
					  3**(-2./3) * (29123./65025) * omega**(10./3) )
				# in case this estimate exceeds 1 by a little, bring it back to 1
				if output > 1: output = 1
			return output
		# helper function
		# inputs: an array of temperature correction values and 
		#	an array of x = abs(cos(theta)) values
		def G(F_arr, x_arr):
			return np.sqrt(F_arr * (1 - x_arr**2) + x_arr**2)
		# output: full Newton's method step function
		# inputs: F, x, G, rho and omega
		def dF_full(F_arr, x_arr, G_arr, rho_arr, omega):
			mult = -2 * F_arr * G_arr**2
			add1 = (1 - G_arr) / x_arr**2
			add2 = (-1./3) * G_arr * rho_arr**3 * omega**2
			add3 = G_arr * np.log( np.sqrt(F_arr) * (1 + x_arr) / (x_arr + G_arr) ) / x_arr**3
			output = mult * (add1 + add2 + add3)
			return output
		# output: series approximation of Newton's method step function up to third order
		# inputs: F, x, G, rho and omega
		def dF_approx(F, x, omega):
			# helper variables and arrays
			x2 = x**2
			o2 = omega**2
			output = (2*F)/3. + (2*x2)/5. - F**1.5*x2*(1 - o2) - \
				(F**2.5*(10*(1 - o2)**2 + 3*x2*(-3 + 8*o2)))/(15.*(1 - o2))
			return output

		# star / surface parameters
		surf = self.surface
		omega = surf.omega
		o2 = omega**2
		F0 = (1 - o2 * rho0**3)**(-2./3) # F at x = 0
		F1 = math.exp((2./3) * o2 * rho1**3) # F at x = 1
		# absolute value of cosine theta, a.k.a. x
		x_arr = np.abs(z_arr / (surf.f * rho_arr)) 
		# optimal smallest value of x for which to compute using Newton's method
		r = np.finfo(float).eps # resolution of floating point numbers
		A = 1 # a parameter for estimating this value of x
		x1 = X1(omega, A, r)
		# a mask that says which x elements are far enough from zero
		# that we use the full Newton's method step function;
		# for the remaining elements, we use order three series expansion of the 
		# step function instead
		fnm = np.greater(x_arr, x1)
		# initialize the result array (to the half-way point in the possible range)
		F_arr = np.full_like(z_arr, (F0 + F1) / 2)
		# Newton's algorithm; 
		# ten steps is about twice the number we need
		for i in range(10):
			# a helper function evaluated
			G_arr = G(F_arr[fnm], x_arr[fnm])
			# the new values of F at the locations 
			# where we use the full Newton's method step function
			F_arr[ fnm ] = F_arr[fnm] + \
				dF_full(F_arr[fnm], x_arr[fnm], G_arr, rho_arr[fnm], omega)
			# the new values of F at the locations 
			# where we use the series expansion of Newton's method step function
			F_arr[ ~fnm ] = F_arr[~fnm] + \
				dF_approx(F_arr[~fnm], x_arr[~fnm], omega)
			# check if we end up outside the bounds on F 
			# and come back into the bounds if we did
			F_arr[ (F_arr < F1) ] = F1
			F_arr[ (F_arr > F0) ] = F0
		# return
		return (F_arr, F0, F1)

	# Interpolation of fit coefficients w.r.t. gravity and temperature
	# Inputs: array of log gravities at all values of z
	#	array of temperatures at all values of z
	# 	fit coefficients on a grid of temperatures and log gravities
	# 	temperature interpolation method: 'linear', 'log, 'planck'
	# Result: set an attribute array of interpolated fit coefficients with
	# 	index 0: z
	# 	index 1: wavelength
	# 	index 2: parameter index 
	# Note: data types of gravity, temperature and intensity fit parameters 
	# in the limb darkening information should be no more than 6 decimal digits, 
	# to conserve memory
	def interp(self, logg_arr, temp_arr, ld, temp_method):
		
		# The temperature-dependent factor in the Planck function
		# Inputs: temperature and frequency arrays of the same dimensions
		# Output: the factor array, with the same dimensions as each input array
		def planck(T, nu):
			return 1. / ( np.exp( ut.h * nu / (ut.k * T) ) - 1 )
		
		# parameter values of the fit coefficients grid
		wl = ld.wl_arr # wavelengths
		g = ld.g_arr # log gravities
		T = ld.temp_arr # temperatures
		# lengths of the fit coefficients grid
		n_p = ft.Fit.m * ft.n # number of fit coefficients at given physical parameters
		n_wl = len(wl) # number of wavelengths
		# number of z values
		n_z = len(self.z_arr)

		# for each value of z, see where in the limb darkening arrays its values of 
		# log g and temperature are; if the resulting index of an array is ind, 
		# the computed value is greater than or equal to array[ind] and less than array[ind + 1]
		ig = np.searchsorted(g, logg_arr, side='right') - 1
		iT = np.searchsorted(T, temp_arr, side='right') - 1
		# for each value of z, 
		# find the values of gravity and temperature between which we are interpolating
		g1 = g[ig] 
		g2 = g[ig + 1] 
		T1 = T[iT]
		T2 = T[iT + 1]
		# fit parameters for the four nearest neighbors at each wavelength, for each value of z;
		# entries are numpy.NAN when no information is available
		# limb darkening array:
		#	index 1: temperature
		#	index 2: gravity
		#	index 3: wavelength
		#	index 4: parameter index
		# result:
		# 	index 0: z
		# 	index 1: wavelength
		# 	index 2: parameter index
		f11 = ld.fit_params[iT    , ig    ]
		f21 = ld.fit_params[iT + 1, ig    ]
		f12 = ld.fit_params[iT    , ig + 1]
		f22 = ld.fit_params[iT + 1, ig + 1]
		## The following assumes that limb darkening (LD) information for some of the 
		## lowest gravity values at each temperature is missing, as in Castelli and Kurucz 2004
		# boolean arrays, one for each of the neighboring temperatures,
		# saying whether the LD information is missing for the lower gravity value at that temperature
		# use the first wavelength and the first parameter to perform this check
		noinfo1 = np.isnan(f11[:, 0, 0])
		noinfo2 = np.isnan(f21[:, 0, 0])
		noinfo = np.logical_or(noinfo1, noinfo2)
		if np.any(noinfo):
			th = np.get_printoptions().get('threshold')
			np.set_printoptions(threshold=10) # ensure summary printout of large arrays
			print ('We do not have the data to interpolate to find intensity at z = ' + str(self.z_arr[noinfo]) + \
				', where the temperatures are ' + str(temp_arr[noinfo]) + ' and log gravity values are ' +\
				 str(logg_arr[noinfo]) + '. At each of these points, we extrapolate: for each temperature, ' +\
				 'we use the intensity information at the closest gravity where such information is available.')
			np.set_printoptions(threshold=th)
			## at each of the neighboring temperatures, and each z value,
			## keep adding to the gravity index until fit parameters are available
			# gravity indices for lower and upper temperature neighbors, respectively
			ig1, ig2 = [np.copy(ig), np.copy(ig)]
			while np.any(noinfo1):
				## at z values where the lower temperature neighbor's information is missing,
				# increment the gravity limb darkening index 
				ig1[noinfo1] += 1
				# set the fit parameters for both the upper and the lower gravity neighbors
				# to those at the new gravity indices 
				f11[noinfo1] = f12[noinfo1] = ld.fit_params[ iT[noinfo1], ig1[noinfo1] ]
				# update the boolean array saying which lower temperature neighbor info is missing
				noinfo1 = np.isnan(f11[:, 0, 0])
			while np.any(noinfo2):
				## at z values where the upper temperature neighbor information is missing,
				# increment the gravity limb darkening index
				ig2[noinfo2] += 1
				# set the fit parameters for both the upper and the lower gravity neighbors
				# to those at the new gravity indices 
				f21[noinfo2] = f22[noinfo2] = ld.fit_params[ iT[noinfo2] + 1, ig2[noinfo2] ]
				# update the boolean array saying which upper temperature neighbor info is missing
				noinfo2 = np.isnan(f21[:, 0, 0])	

		# if the temperature scale is non-linear, 
		# modify the temperature and gravity arrays accordingly
		if temp_method == 'log':
			T1 = np.log10(T1)
			T2 = np.log10(T2)
			temp_arr = np.log10(temp_arr)
		elif temp_method == 'planck':
			# convert the wavelengths to frequencies
			nu = ut.color_nm_Hz(wl)
			# create grids based on the arrays of temperatures at all values of z and the frequency array
			# index 0: z
			# index 1: frequency
			nn, tt = np.meshgrid(nu, temp_arr)
			NN1, TT1 = np.meshgrid(nu, T1)
			NN2, TT2 = np.meshgrid(nu, T2)
			# replace the 1D temperature with 2D arrays evaluated on the grids
			temp_arr = planck(tt, nn)
			T1 = planck(TT1, NN1)
			T2 = planck(TT2, NN2)
			# add the wavelength dimension to the gravity arrays
			logg_arr = logg_arr[:, np.newaxis]
			g1 = g1[:, np.newaxis]
			g2 = g2[:, np.newaxis]

		## bilinear interpolation (see Wikipedia: Bilinear Interpolation: Algorithm)
		const = (1 / ((g2 - g1) * (T2 - T1)))
		Dg1 = g2 - logg_arr
		Dg2 = logg_arr - g1
		Dt1 = T2 - temp_arr
		Dt2 = temp_arr - T1
		# helper matrix
		w11 = const * Dt1 * Dg1
		w21 = const * Dt2 * Dg1
		w12 = const * Dt1 * Dg2
		w22 = const * Dt2 * Dg2
		# an array of fit function parameters 
		# index 0: z
		# index 1: wavelength
		# index 2: parameter index
		if temp_method == 'planck': 
			# add the parameter index dimension to the helper matrices
			params_arr = \
				f11 * w11[:, :, np.newaxis] + \
				f21 * w21[:, :, np.newaxis] + \
				f12 * w12[:, :, np.newaxis] + \
				f22 * w22[:, :, np.newaxis]
		else:
			# add the wavelength and parameter indices to the helper matrices
			params_arr = \
				f11 * w11[:, np.newaxis, np.newaxis] + \
				f21 * w21[:, np.newaxis, np.newaxis] + \
				f12 * w12[:, np.newaxis, np.newaxis] + \
				f22 * w22[:, np.newaxis, np.newaxis]
		return params_arr
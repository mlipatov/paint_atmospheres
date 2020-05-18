from pa.lib import fit as ft
from pa.lib import util as ut
import numpy as np
import math
import sys

np01 = np.array([0, 1]) # a numpy array containing 0 and 1

class ConvergenceError(Exception):
	""" Error raised during interpolation of intensity coefficients over gravity and temperature """
	def __init__(self, message):
		self.message = message

class InterpolationError(Exception):
	""" Error raised during interpolation of intensity coefficients over gravity and temperature """
	def __init__(self, message):
		self.message = message

class Map:
	""" At initialization with a surface, a step in z values and limb darkening information, 
	computes a number of equally spaced discrete z values on (-1, 1), as well as
	wavelength-dependent parameters of the intensity fits and area elements at these z values. 
	Gravity and temperature calculations are based on Espinosa Lara 2011 (EL)."""

	# initialize with the surface shape of the star, the number of z values to use for integration
	# on the upper half of the star, 
	# the constants needed for temperature and gravity unit conversions, limb darkening information,
	# temperature interpolation method, gravity interpolation method and the number of steps in the
	# Newton's method involved in temperature calculation
	def __init__(self, surf, nz, add_logg, mult_temp, ld, temp_method, g_method, nm):

		## initialize the surface
		self.surface = surf # surface shape information
		omega = surf.omega

		# store gravity and temperature modifiers, 
		# as well as the interpolation methods
		self.add_logg = add_logg
		self.mult_temp = mult_temp
		self.temp_method = temp_method
		self.g_method = g_method

		## initialize the array of z values for integration
		# for the modified trapezoidal "gamma" method, this should be an integral factor of 4 plus 1
		# non-negative z values for the upper half of the star;
		# non-positive values for the lower half are the reversed negatives of this
		self.z_up = np.linspace(0, 1, nz + 1)
		# store the spacing between the z values
		self.dz = self.z_up[1] - self.z_up[0]
		# store the number of z values used for the integration on the upper half of the star
		self.nz = nz

		# cylindrical coordinate r for the upper portion of the star;
		# values for the lower half are a reverse of this
		r_up = self.surface.R( self.z_up )

		# area elements for the upper portion of the star;
		# values for the lower half are a reverse of this
		self.A_up = surf.A( self.z_up )	

		# for the upper half of the star, compute and store the temperatures, fit parameters 
		# and the information about points where we extrapolate instead of interpolating;
		#	this information is a 2D array; each element is an array containing 
		#	[z-coordinate, log gravity, temperature];
		# for the lower half of the star, the corresponding arrays are a reverse of these
		self.temp_up, self.params_up, self.extr_info = self.Tp( self.z_up, r_up, ld, nm )


	# output: temperatures and intensity fit parameters for a set of locations
	# input: z values, r values and limb darkening information
	def Tp(self, z, r, ld, nm):		
		# compute the spherical coordinate rho
		rho = self.surface.rho( r, z )
		r01 = self.surface.R( np01 ) # a numpy array containing r at z = 0 and r at z = +/-1
		rho0, rho1 = self.surface.rho( r01, np01 ) # rho at z = 0 and +/-1
		
		## compute the effective gravitational acceleration in units of G M / Re**2 
		geff = self.geff(rho, r)
		# convert to log10(gravity in cm / s**2) and to a less memory-intensive data type
		logg = self.add_logg + np.log10(geff)
		logg = logg.astype(np.float32)
		
		## compute the effective temperature in units of ( L / (4 pi sigma Re**2) )**(1/4), 
		## as in EL eqn 31, for each z; then convert it to Kelvin
		[F, F0, F1] = self.F(z, rho, rho0, rho1, nm) # temperature correction values
		# temperatures
		temp = self.mult_temp * ( geff * F )**(1./4)
		# convert temperature values to a less memory-intensive data type
		temp = temp.astype(np.float32)

		# compute the interpolated values of limb darkening fit parameters
		# if given limb darkening information
		if ld is not None:
			params, extr_info = self.interp(z, logg, temp, ld, self.temp_method, self.g_method)
		else:
			params = None
			extr_info = None

		# return
		return [temp, params, extr_info]
	
	# output: effective gravitational acceleration in units of G M / Re**2 
	# 	as in equations 7 and 31 of EL, for each z
	# inputs: arrays of rho and r coordinates
	def geff(self, rho, r):
		omega = self.surface.omega
		r_sq = np.power(r, 2) # r squared
		geff = np.sqrt( rho**-4 + omega**4 * r_sq - 2 * omega**2 * r_sq * rho**-3 )
		return geff

	# output: temperature correction factors (EL equation 26)
	# 	at every value of z; runs the Newton's method algorithm that acts on the entire array 
	# 	of z values that aren't too close to 0 at its every step; for values that 
	# 	are close to 0, returns the temperature correction expected at z = 0
	def F(self, z_arr, rho_arr, rho0, rho1, nm):
		# output: smallest value of x for which to compute 
		#	using full Newton's method step function; a.k.a. x_b
		# inputs: rotational velocity of the star,
		#	order-one factor of proportionality between the error 
		#		in full step function and that in the series expansion
		#	resolution of floating point numbers
		def X1(omega, A, r):
			# a factor that's a little less than 1; above 0.95 is enough for r < 1
			B = 0.9
			# omega below which the approximation yields values greater than 1
			omega_lim = B * (2./(85*A))**(1./4) * 3**(1./2) * r**(1./4)
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
		def G(F, x):
			return np.sqrt(F * (1 - x**2) + x**2)
		# output: full Newton's method step function
		# inputs: F, x, G, rho and omega
		def dF_full(F, x, G, rho, omega):
			mult = -2 * F * G**2
			add1 = (1 - G) / x**2
			add2 = (-1./3) * G * rho**3 * omega**2
			add3 = G * np.log( np.sqrt(F) * (1 + x) / (x + G) ) / x**3
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
		A = 1000 # a parameter for estimating this value of x; optimized for the maximum omega
		x1 = X1(omega, A, r)
		# a mask that says which x elements are far enough from zero
		# that we use the full Newton's method step function;
		# for the remaining elements, we use order three series expansion of the 
		# step function instead
		fnm = np.greater(x_arr, x1)
		# initialize the result array (to the half-way point in the possible range)
		F_arr = np.full_like(z_arr, (F0 + F1) / 2)
		# Newton's algorithm; 
		# we never need more than 11 steps
		for i in range(nm):
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

			# # check for convergence
			# if i > 0: 
			# 	diff = np.abs(F_arr - F_prev)
			# 	print(i + 1, diff.max())
			# F_prev = np.copy(F_arr)

		# return
		return (F_arr, F0, F1)

	# Interpolation and extrapolation of fit coefficients w.r.t. gravity and temperature
	# Inputs: array of z
	#	array of log gravities at all values of z
	#	array of temperatures at all values of z
	# 	fit coefficients on a grid of temperatures and log gravities
	# 	temperature interpolation method: 'linear', 'log, 'planck'
	# Outputs: 
	#	3D array of interpolated fit coefficients
	# 		0: z
	# 		1: wavelength
	# 		2: parameter index 
	#	extrapolation info, which is either None or an array of size-3 arrays
	#		[z, gravity, temperature]
	# Note: data types of gravity, temperature and intensity fit parameters 
	# in the limb darkening information should be no more than 6 decimal digits, 
	# to conserve memory
	def interp(self, z_arr, g_arr, temp_arr, ld, temp_method, g_method):
		
		# The temperature-dependent factor in the Planck function
		# Inputs: temperature and frequency arrays of the same dimensions
		# Output: the factor array, with the same dimensions as each input array
		def planck(T, nu):
			return 1. / ( np.exp( ut.h * nu / (ut.k * T) ) - 1 )
		
		# parameter values of the fit coefficients grid
		wl = ld.wl_arr # wavelengths
		g = ld.g_arr # log gravities
		T = ld.temp_arr # temperatures

		# if temperature is outside the bounds of the limb darkening grid at any location, 
		# raise an error
		below = np.any(temp_arr < T[0])
		above = np.any(temp_arr > T[-1]) 
		if below or above:
			message = 'Surface temperature is outside the bounds of the limb darkening information grid:'
			if below: message += ' at some points, it is below ' + str(T[0]) + '.'
			if above: message += ' at some points, it is above ' + str(T[-1]) + '.'		
			raise InterpolationError(message)
		# if gravity is outside the bounds of the limb darkening grid by more than 0.5 dex at any location, 
		# raise an error
		below = np.any(g_arr < g[0] - 0.5)
		above = np.any(g_arr > g[-1] + 0.5) 
		if below or above:
			message = 'Surface gravity is outside the extrapolation bounds of the limb darkening information grid.'
			if below: message += ' Gravity at some points is below ' + str(g[0] - 0.5) + '.'
			if above: message += ' Gravity at some points is above ' + str(g[-1] + 0.5) + '.'
			raise InterpolationError(message)

		# for each z, see where in the limb darkening (LD) arrays the values of 
		# log g and temperature are; if the resulting index of an array is ind, 
		# the computed value is greater than or equal to array[ind] and less than array[ind + 1];
		# here, this is valid only for z values where extrapolation is not required
		ig = np.searchsorted(g, g_arr, side='right') - 1
		iT = np.searchsorted(T, temp_arr, side='right') - 1
		# initialize the indices where the parameter values will be extracted;
		# the initial values are only valid where extrapolation isn't required
		iT11, ig11 = [np.copy(iT), np.copy(ig)] # lower temperature, lower gravity neighbor
		iT12, ig12 = [np.copy(iT), np.copy(ig) + 1] # lower temperature, upper gravity neighbor
		iT21, ig21 = [np.copy(iT) + 1, np.copy(ig)] # upper temperature, lower gravity
		iT22, ig22 = [np.copy(iT) + 1, np.copy(ig) + 1] # upper temperature, upper gravity
		# initialize the values of gravity and temperature between which we are interpolating
		g1 = np.full_like(ig, np.nan, dtype=float)
		g2 = np.full_like(ig, np.nan, dtype=float)
		T1 = np.full_like(iT, np.nan, dtype=float)
		T2 = np.full_like(iT, np.nan, dtype=float)

		# initialize a boolean array of locations where extrapolation is required
		extra = np.full_like(z_arr, False, dtype=bool)
		T1 = T[ iT ] # set the lower temperature values

		# locations where gravity is below the lower bound of the LD grid
		l = (g_arr < g[0])
		extra = np.logical_or(extra, l) # extrapolation is needed here
		ig11[l] = ig12[l] = ig21[l] = ig22[l] = 0 # set gravity indices to the first possible index here
		g1[l] 	= g[0] - 0.5 	# set the lower gravity value to be 0.5 dex below lower bound of LD grid
		g1[~l] 	= g[ ig[~l] ]	# set all other lower gravity values
		# locations where gravity is above the upper bound of the LD grid
		l = (g_arr > g[-1])
		extra = np.logical_or(extra, l) # extrapolation is needed here
		l = np.logical_or(l, g_arr == g[-1]) # add the locations where gravity is at the upper bound
		ig11[l] = ig12[l] = ig21[l] = ig22[l] = len(g) - 1 # set gravity indices to the last index
		g2[l]	= g[-1] + 0.5 # set the upper gravity value to be above the upper bound
		g2[~l]	= g[ ig[~l] + 1 ]  # set all other upper gravity values
		# locations where temperature is exactly at the upper bound of the limb darkening grid
		l = (temp_arr == T[-1])
		iT11[l] = iT12[l] = iT21[l] = iT22[l] = len(T) - 1 # set temperature indices to the last possible index
		T2[l]	= T[-1] + 1000 # set the upper temperature value to be above the upper bound of the LD grid
		T2[~l]	= T[ iT[~l] + 1 ] # set other upper temperature values

		# fit parameters for the four nearest neighbors at each wavelength, for each value of z;
		# entries are numpy.NAN when no information is available; the initial values need to be mofified
		# at z where extrapolation is necessary
		# limb darkening array:
		#	0: temperature
		#	1: gravity
		#	2: wavelength
		#	3: parameter index
		# result:
		# 	index 0: z
		# 	index 1: wavelength
		# 	index 2: parameter index
		f11 = ld.fit_params[iT11, ig11]
		f21 = ld.fit_params[iT21, ig21]
		f12 = ld.fit_params[iT12, ig12]
		f22 = ld.fit_params[iT22, ig22]

		# locations where gravity is within the grid, 
		# but the lower gravity neighbor limb darkening information at the lower temperature is missing
		l = np.logical_and(np.logical_and(g_arr >= g[0], g_arr < g[-1]), np.isnan(f11[:, 0, 0]))
		extra = np.logical_or(extra, l) # extrapolation is needed here
		ig11[l] += 1 # increment the lower gravity / lower temperature indices
		# update the fit parameters
		f11 = ld.fit_params[ iT11, ig11 ]
		f12 = ld.fit_params[ iT12, ig12 ]
		# update the locations
		l = np.logical_and(l, np.isnan(f11[:, 0, 0]))
		if np.any(l): # if there are still parameters missing at lower temperature
			message = 'At some points, surface gravity at the lower-temperature limb darkening (LD) '+\
				'grid neighbor differs by more than 0.5 dex from the closest value with available '+\
				'LD information.'
			raise InterpolationError(message)

		# locations where gravity is within the grid, 
		# but the lower gravity neighbor limb darkening information at the upper temperature is missing
		l = np.logical_and(np.logical_and(g_arr >= g[0], g_arr < g[-1]), np.isnan(f21[:, 0, 0]))
		extra = np.logical_or(extra, l) # update the extrapolation mask
		ig21[l] += 1 # increment the lower gravity / upper temperature indices
		# update the fit parameters
		f21 = ld.fit_params[ iT21, ig21 ]
		f22 = ld.fit_params[ iT22, ig22 ]
		# update the locations
		l = np.logical_and(l, np.isnan(f21[:, 0, 0]))
		if np.any(l): # if there are still parameters missing
			message = 'At some points, surface gravity at the upper-temperature limb darkening (LD) '+\
				'grid neighbor differs by more than 0.5 dex from the closest value with available '+\
				'LD information.'
			raise InterpolationError(message)

		# set the extrapolation info
		if np.any(extra):
			extr_info = np.stack( (self.z_arr[extra], g_arr[extra], temp_arr[extra]), axis=-1)
		else: 
			extr_info = None

		# if the gravity scale is linear,
		# modify the gravity array accordingly
		if g_method == 'lin':
			g1 = np.power(10, g1, dtype=np.float32)
			g2 = np.power(10, g2, dtype=np.float32)
			g_arr = np.power(10, g_arr, dtype=np.float32)

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
			g_arr = g_arr[:, np.newaxis]
			g1 = g1[:, np.newaxis]
			g2 = g2[:, np.newaxis]

		## bilinear interpolation (see Wikipedia: Bilinear Interpolation: Algorithm)
		const = (1. / ((g2 - g1) * (T2 - T1)))
		Dg1 = g2 - g_arr
		Dg2 = g_arr - g1
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
		return [params_arr, extr_info]
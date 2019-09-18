import limbdark.fit as ft
import limbdark.limbdark as limbdark
import star.surface as sf
import util as ut
import numpy as np
import math
import sys, time
from scipy import ndimage
from scipy.interpolate import griddata

class Map:
	""" Contains a map of intensity integrals and the wavelength-dependent parameters of the 
	corresponding intensity fits, area elements, and
	gravity and temperature values across a set of z values on the surface of a rotating star. 
	Gravity and temperature calculations are based on Espinosa Lara 2011 (EL)."""

	# initialize with the surface shape of the star, the step in z, the constants
	# needed for temperature and gravity unit conversions, and limb darkening information
	def __init__(self, surf, z_step, add_logg, mult_temp, ld):
		## initialize surface and limb darkening information
		z1 = surf.z1 # an integration bound on the surface
		self.z_step = z_step
		self.ld = ld # limb darkening information
		self.surf = surf # surface shape information
		omega = self.surf.omega
		## initialize the array of z values for integration
		za = np.arange(-1 + z_step / 2, 1 - z_step / 2, z_step)
		self.z_arr = za[za >= -z1]
		## compute a 2D array of fit function integrals, 
		## one for each z value and each parameter / interval combination of the fit
		self.fitint = np.zeros( (len(self.z_arr), ft.n * ft.Fit.m) ) # initialize the fit functions
		## compute the integrated fit functions
		print ("Integrating the fit functions...")        
		sys.stdout.flush()
		c = 0
		for z in self.z_arr: 
			# at z = -z1, the integral is zero
			if c == 0 and z == -z1:
				self.fitint[0] = np.zeros(ft.n * ft.Fit.m)
			# for z != -z1, compute the integral
			else:
				a, b = surf.ab(z)
				belowZ1 = (z < z1)
				self.fitint[c] = ft.Fit.integrate(belowZ1, a, b)
			c += 1
		## compute area elements for integration
		## as well as cylindrical coordinate r and the spherical coordinate rho
		print ("Computing the area elements, cylindrical coordinates and spherical coordinates...")        
		sys.stdout.flush()
		self.A_arr = surf.A( self.z_arr )
		self.r_arr = np.array([ surf.R(z) for z in self.z_arr ])
		self.rho_arr = surf.rho(self.r_arr, self.z_arr)
		r0, r1 = [ surf.R(0), surf.R(1) ] # r at z = 0 and 1
		self.rho0, self.rho1 = [ surf.rho(r0, 0), surf.rho(r1, 1) ] # rho at z = 0 and 1
		## compute the effective gravitational acceleration in units of G M / Re**2 
		print ("Computing gravities and temperatures...")        
		sys.stdout.flush()
		geff_arr = self.geff(self.rho_arr, self.r_arr)
		g0 = self.geff(self.rho0, r0) # gravity at z = 0 (minimum)
		g1 = self.geff(self.rho1, r1) # gravity at z = 1 (maximum)
		# convert to log10(gravity in cm / s**2)
		self.logg_arr = add_logg + np.log10(geff_arr)
		self.logg0 = add_logg + math.log10(g0)
		self.logg1 = add_logg + math.log10(g1)
		## compute the effective temperature in units of ( L / (4 pi sigma Re**2) )**(1/4), 
		## as in EL eqn 31, for each z; then convert it to Kelvin
		self.F() # compute all the F values
		self.temp_arr = mult_temp * self.Teff( geff_arr, self.F_arr ) # the temperatures
		self.T0 = mult_temp * self.Teff( g0, self.F0 ) # minimum temperature, at z = 0
		self.T1 = mult_temp * self.Teff( g1, self.F1 ) # maximum temperature, at z = 1
		## compute the interpolated values of limb darkening fit parameters
		print ("Interpolating the fit parameters...")        
		sys.stdout.flush()
		self.interp()
	
	# returns the effective gravitational acceleration in units of G M / Re**2 
	# as in equations 7 and 31 of EL, for each z
	# inputs: arrays of rho and r coordinates
	def geff(self, rho_arr, r_arr):
		omega = self.surf.omega
		r_arr_sq = np.power(r_arr, 2) # r squared
		geff_arr = np.sqrt( rho_arr**-4 + omega**4 * r_arr_sq - \
			2 * omega**2 * r_arr_sq * rho_arr**-3 )
		return geff_arr

	# compute the effective temperature in units of ( L / (4 pi sigma Re**2) )**(1/4), 
	# as in EL eqn 31, for each z
	# inputs: arrays of effective gravity and the temperature correction
	def Teff(self, geff_arr, F_arr):
		return ( geff_arr * F_arr )**(1./4)

	# computes the temperature correction factors (EL equation 26)
	# at every value of z; runs the Newton's method algorithm that acts on the entire array 
	# of z values that aren't too close to 0 at its every step; for values that 
	# are close to 0, returns the temperature correction expected at z = 0
	def F(self):
		# helper function
		# inputs: squared rotation speed, 
		#	an array of rho and an array of x = cos(theta)
		def add(o2, rho_arr, x_arr):
			return -x_arr - (1./3) * o2 * rho_arr**3 * x_arr**3
		# helper function
		# inputs: an array of temperature correction values and 
		#	an array of x = cos(theta) values
		def G(F_arr, x_arr):
			return np.sqrt(F_arr * (1 - x_arr**2) + x_arr**2)
		# the function whose root we are looking for, 
		# derived from EL24, with x = cos(theta) and F as defined in EL26
		def f(F_arr, x_arr, G_arr, add_arr): 
			result1 = x_arr / G_arr
			result2 = np.log( np.sqrt(F_arr) * (1 + x_arr) / (x_arr + G_arr) )
			return result1 + result2 + add_arr
		# the derivative of the function we are looking for,
		# with respect to the temperature correction F
		def Df(F_arr, x_arr, G_arr):
			result = (x_arr / G_arr)**3 / (2 * F_arr)
			return result
		# a step of the Newton's method
		# inputs: an array of F values, corresponding arrays of x values and values of 
		#	the additive term that doesn't depend on F
		def step(F_arr, x_arr, add_arr):
			# helper function at the current values of F
			G_arr = G(F_arr, x_arr)
			# the function whose root we are looking for
			f_arr = f(F_arr, x_arr, G_arr, add_arr)
			# the derivative of the function whose root we are looking for
			Df_arr = Df(F_arr, x_arr, G_arr)
			# step
			result = F_arr - f_arr / Df_arr
			return result
		# x is defined as abs(cos(theta))
		# optimal smallest value of x for which to compute using Newton's method
		# assumes float precision = 1e-15
		delta = 0.0015 
		# star / surface parameters
		surf = self.surf
		omega = surf.omega
		o2 = omega**2
		F0 = (1 - o2 * self.rho0**3)**(-2./3) # F at x = 0
		F1 = math.exp((2./3) * o2 * self.rho1**3) # F at x = 1
		z_arr = self.z_arr
		rho_arr = self.rho_arr
		# absolute value of cosine theta, a.k.a. x
		x_arr = np.abs(z_arr / (surf.f * rho_arr)) 
		# a boolean array that says which x elements are so close to zero
		# that we don't want to be using Newton's method
		mask = np.less(x_arr, delta)
		# for Newton method's operation, remove such elements from the arrays
		x_arr = x_arr[~mask]
		rho_arr = rho_arr[~mask] # get the corresponding values of rho
		# an additive term that doesn't depend on F, for every z
		add_arr = add(o2, rho_arr, x_arr) 
		# initial estimates of F
		F_arr = np.full_like(x_arr, (F0 + F1) / 2)
		# Newton's algorithm; 
		# ten steps is about twice the number we need
		for i in range(10):
			# compute the new values of F
			F_arr = step(F_arr, x_arr, add_arr)
			# check if we end up below the lower bound on F 
			# and come back to that lower bound if we did overshoot it
			F_arr[F_arr < F1] = F1
		# reinsert the temperature correction at values of x very close to zero
		ind = np.nonzero(mask)[0]
		F_arr = np.insert(F_arr, ind, F0)
		# save the results and return
		self.F0 = F0
		self.F1 = F1
		self.F_arr = F_arr

	# for each wavelength and intensity fit parameter, interpolate the parameter
	# as a function of temperature and gravity, at each z;
	def interp(self):
		# initialize local variables
		ld = self.ld
		ld.fit_params = np.array(ld.fit_params)
		wl = ld.wl_arr
		fit_params = ld.fit_params
		n_p = ft.Fit.m * ft.n
		n_wl = len(wl)
		n_z = len(self.z_arr)
		# an array of parameters, indexed by the wavelength, value of z and parameter index
		params = np.empty( (n_wl, n_z, n_p) )
		# for each value of z, look to see where in the limb darkening arrays these values of 
		# log g and temperature are; if the resulting index of an array is ind, 
		# the computed value is greater than or equal to array[ind] and less than array[ind + 1]
		ig = np.searchsorted(ld.g_arr, self.logg_arr, side='right') - 1
		iT = np.searchsorted(ld.temp_arr, self.temp_arr, side='right') - 1
		# for each value of z, 
		# find the values of log g and temperature between which we are interpolating
		g1 = ld.g_arr[ig] 
		g2 = ld.g_arr[ig + 1] 
		T1 = ld.temp_arr[iT]
		T2 = ld.temp_arr[iT + 1]
		# fit parameters for the four nearest neighbors at each wavelength, for each value of z;
		# entries are numpy.NAN when no information is available
		# limb darkening array:
		#	index 1: temperature
		#	index 2: gravity
		#	index 3: wavelength
		#	index 4: parameter index
		# result:
		# 	index 1: z
		# 	index 2: wavelength
		# 	index 3: parameter index
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
			print ('We do not have the data to interpolate to find intensity at z = ' + str(self.z_arr[noinfo]) + \
				', where the temperatures are ' + str(self.temp_arr[noinfo]) + ' and log gravity values are ' +\
				 str(self.logg_arr[noinfo]) + '. At each of these points, we extrapolate: for each temperature, ' +\
				 'we use the intensity information at the closest gravity where such information is available.')
			## at each of the neighboring temperatures, 
			## keep adding to the gravity index until fit parameters are available
			# gravity indices for lower and upper temperature neighbors, respectively
			ig1, ig2 = [np.copy(ig), np.copy(ig)]
			while np.any(noinfo1):
				ig1 += 1
				f11 = f12 = ld.fit_params[iT, ig1]
				noinfo1 = np.isnan(f11[:, 0, 0])
			while np.any(noinfo2):
				ig2 += 1
				f21 = f22 = ld.fit_params[iT + 1, ig2]
				noinfo1 = np.isnan(f21[:, 0, 0])				
		## bilinear interpolation (see Wikipedia: Bilinear Interpolation: Algorithm)
		const = (1 / ((g2 - g1) * (T2 - T1)))
		Dg1 = g2 - self.logg_arr
		Dg2 = self.logg_arr - g1
		Dt1 = T2 - self.temp_arr
		Dt2 = self.temp_arr - T1
		self.params_arr = const[:, np.newaxis, np.newaxis] * (\
			f11 * Dt1[:, np.newaxis, np.newaxis] * Dg1[:, np.newaxis, np.newaxis] + \
			f21 * Dt2[:, np.newaxis, np.newaxis] * Dg1[:, np.newaxis, np.newaxis] + \
			f12 * Dt1[:, np.newaxis, np.newaxis] * Dg2[:, np.newaxis, np.newaxis] + \
			f22 * Dt2[:, np.newaxis, np.newaxis] * Dg2[:, np.newaxis, np.newaxis])

	# using the pre-calculated mapped features, integrate to find the light at all wavelengths,
	# in ergs/s/Hz/ster per squared equatorial radius of the star
	def intgrt(self):
		# at each wavelength and z value, obtain the integral of the total fit function over phi;
		# to do this, sum up the products of the fit parameters and the corresponding fit integrals
		# along the fit parameter dimension
		fit_arr = np.sum(self.params_arr * self.fitint[:, np.newaxis, :], axis=2)
		# at each wavelength, sum up the product of the phi integral of the fit function and
		# the dimensionless area element at each z, multiply by the z step
		return self.z_step * np.sum(self.A_arr[:, np.newaxis] * fit_arr, axis=0)

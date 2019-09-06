import limbdark.fit as ft
import limbdark.limbdark as limbdark
import star.surface as sf
import util as ut
import numpy as np
import math
import sys, time
from scipy.interpolate import Rbf

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
		self.ld = ld
		self.f = surf.f # flatness of the surface
		self.omega = surf.omega # rotational velocity of the star with this surface
		omega = self.omega
		## initialize the array of z values
		za = np.arange(-z1, 1, z_step)
		self.z_arr1 = za[za <  0] # negative z values
		self.z_arr2 = za[za >= 0] # non-negative z values
		self.z_arr = za # an array of z values for integration
		## compute a 2D array of fit function integrals, 
		## one for each z value and each parameter / interval combination of the fit
		self.fitint = np.zeros( (len(self.z_arr), ft.n * ft.Fit.m) )
		## compute the integrated fit functions
		print ("Integrating the fit functions...")        
		sys.stdout.flush()
		self.fitint[0] = np.zeros(ft.n * ft.Fit.m) # at z = -z1, the integral is zero
		c = 1
		for z in self.z_arr[1:]: 
			a, b = surf.ab(z)
			belowZ1 = (z < z1)
			self.fitint[c] = ft.Fit.integrate(belowZ1, a, b)
			c += 1
		## compute an array of area elements for integration, one for each value of z
		print ("Computing the area elements, as well as the cylindrical and spherical coordinates...")        
		sys.stdout.flush()
		self.A_arr = np.array([ surf.A(z) for z in self.z_arr ])
		## compute arrays of the cylindrical coordinate r and the spherical coordinate rho for each z
		self.r_arr1 = np.array([ surf.R(z) for z in self.z_arr1 ])
		self.r_arr2 = np.array([ surf.R(z) for z in self.z_arr2 ])
		self.r_arr = np.concatenate(( self.r_arr1, self.r_arr2 ))
		r_arr_sq = np.power(self.r_arr, 2)
		self.rho_arr1 = np.array([ surf.rho(z) for z in self.z_arr1 ])
		self.rho_arr2 = np.array([ surf.rho(z) for z in self.z_arr2 ])
		self.rho_arr = np.concatenate(( self.rho_arr1, self.rho_arr2 ))
		## compute the effective gravitational acceleration in units of G M / Re**2 
		## as in equations 7 and 31 of EL, for each z
		print ("Computing gravities and temperatures...")        
		sys.stdout.flush()
		geff_arr = np.sqrt( self.rho_arr**-4 + omega**4 * r_arr_sq - \
			2 * omega**2 * r_arr_sq * self.rho_arr**-3 )
		# convert to log10(gravity in cm / s**2)
		self.logg_arr = add_logg + np.log10(geff_arr)
		## compute the effective temperature in units of ( L / (4 pi sigma Re**2) )**(1/4), 
		## as in EL eqn 31, for each z
		t_arr = geff_arr**(1./4) * self.Tc()
		# convert the temperature to Kelvin
		self.temp_arr = mult_temp * t_arr
		## compute the interpolated values of limb darkening fit parameters
		print ("Interpolating the fit parameters...")        
		sys.stdout.flush()
		self.params_arr = self.interp()
	

	# returns, as an array, the temperature correction factors (EL equations 31 and 26)
	# at every value of z; runs a bisection algorithm that acts on the entire array 
	# of positive, non-one values of z at its every step
	def Tc(self):
		## variable initializations
		# values of the arrays for positive values of z, excluding z = 1
		z_arr = self.z_arr2[1:-1] 
		rho_arr = self.rho_arr2[1:-1]
		r_arr = self.r_arr2[1:-1]
		# values at z = 1 and 0, for which the correction will be computed separately
		rho1 = self.rho_arr2[-1] # rho at z = 1
		rho0 = self.rho_arr2[0] # rho at z = 0
		r1 = self.r_arr2[-1] # rho at z = 1
		r0 = self.r_arr2[0] # rho at z = 0
		# star / surface parameters
		omega = self.omega
		f = self.f

		## array operations
		cosn_arr = z_arr / (f * rho_arr) # cosine theta for every z
		tan2_arr = (rho_arr - z_arr / f) / r_arr # tan(theta / 2) for every z
		tan_arr = f * r_arr / z_arr # tan theta for every z
		# an array of the function tau evaluated at every value of z, 
		# with tau defined in EL eqn 21, with rho and theta both functions of z
		tau_arr = (1./3) * omega**2 * cosn_arr**3 + cosn_arr + np.log(tan2_arr)

		# at z = 1 (see EL eqn 27), the correction factor is
		Tc1 = math.exp((1./6) * omega**2 * rho1**3)
		# at z = 0 (see EL eqn 28), the correction factor is
		Tc0 = (1 - omega**2 * rho0**3)**(-1./6)

		# initialize to zeros the array of curly theta (EL eqn 24) 
		# for the bisection algorithm
		curly_arr = np.zeros_like(z_arr)
		# size of the range of curly theta
		s = math.pi / 2
		# step through the bisection algorithm the number of times 
		# that guarantees convergence within the precision of floats
		n = int(math.ceil(math.log2( s / np.finfo(1.0).resolution )))
		for i in range(n):
			# turn off RuntimeWarning: divide by zero encountered in log;
			# numpy.log evaluates to -inf when curly theta = 0, which is the 
			# the value we will enter in our array at such values of curly theta
			np.seterr(divide = 'ignore')
			# the function whose root we are looking for, evaluated at the 
			# current values of curly theta 
			f1 = np.where(curly_arr == 0, -np.inf, 
				np.cos(curly_arr) + np.log(np.tan(curly_arr / 2)) - tau_arr)
			# turn on RuntimeWarning: divide by zero encountered in log
			np.seterr(divide = 'warn') 
			# the length of the subintervals into which we are subdividing the range of
			# curly theta at this step of the bisection algorithm
			x = s / 2**(i + 1) 
			# the function evaluated at the sum of the current values of curly theta
			# and the length of a subinterval
			f2 = np.cos(curly_arr + x) + np.log(np.tan((curly_arr + x)/ 2)) - tau_arr
			# compute an array of boolean values for each z; the value says whether the algorithm
			# should add the value of the subinterval to the current estimate of curly theta
			# at a given value of z; then add (or don't add) the value of the subinterval accordingly
			curly_arr = curly_arr + x * np.greater(f1 * f2, 0)
		# the temperature correction for positive values of z
		Tc_arr2 = np.concatenate(( np.sqrt(np.tan(curly_arr) / tan_arr), np.array([Tc1]) ))
		# compute the mirror image of the temperature correction for negative z values
		Tc_arr1 = np.flip( Tc_arr2[:len(self.z_arr1)] )
		# compile the temperature corrections and return
		return np.concatenate(( Tc_arr1, np.array([Tc0]), Tc_arr2 ))


	# for each wavelength and intensity fit parameter, interpolate the parameter
	# as a function of temperature and gravity, then calculate the parameter for each z;
	def interp(self):
		# initialize local variables
		ld = self.ld
		wl_arr = ld.wl_arr
		fit_params = ld.fit_params
		n_params = ft.Fit.m * ft.n
		n_wl = len(wl_arr)
		n_z = len(self.z_arr)
		# an array of parameters, indexed by the wavelength, value of z and parameter index
		params_arr = np.empty( (n_wl, n_z, n_params) )
		# for each value of z, look to see where in the limb darkening arrays these values of 
		# log g and temperature are; if the resulting index of an array is ind, 
		# the computed value is greater than or equal to array[ind] and less than array[ind + 1]
		ind_g_arr = np.searchsorted(ld.g_arr, self.logg_arr, side='right') - 1
		ind_temp_arr = np.searchsorted(ld.temp_arr, self.temp_arr, side='right') - 1
		# for each value of z, 
		# find the values of log g and temperature between which we are interpolating
		g1_arr = ld.g_arr[ind_g_arr] 
		g2_arr = ld.g_arr[ind_g_arr + 1] 
		temp1_arr = ld.temp_arr[ind_temp_arr]
		temp2_arr = ld.temp_arr[ind_temp_arr + 1]
		# fit parameters for the four nearest neighbors at each wavelength, for each value of z;
		# index 1: wavelength
		# index 2: z value
		# index 3: fit parameter index
		f00 = ld.fit_params[ :, ind_g_arr, ind_temp_arr, : ]
		f01 = ld.fit_params[ :, ind_g_arr, ind_temp_arr + 1, : ]
		f10 = ld.fit_params[ :, ind_g_arr + 1, ind_temp_arr, : ]
		f11 = ld.fit_params[ :, ind_g_arr + 1, ind_temp_arr + 1, : ]
		# a boolean array, saying whether for the log g and temperature combination at each z 
		# the fit parameters at either of the four nearest neighbors aren't given;
		# use the first wavelength and the first parameter at each value of z to perform this check
		no_fits = np.logical_or.reduce( (np.isnan(f00[0, :, 0]), np.isnan(f01[0, :, 0]), \
										np.isnan(f10[0, :, 0]), np.isnan(f11[0, :, 0])) )
		if True in no_fits:
			print ('We do not have the data to interpolate to find intensity at z = ' + str(z_arr[no_fits]) + \
				', where the temperatures are ' + str(self.temp_arr[no_fits]) + ' and log gravity values are ' +\
				 str(self.logg_arr[no_fits]) + ".")
		# bilinear interpolation (see Wikipedia: Bilinear Interpolation)
		const = (1 / ((g2_arr - g1_arr) * (temp2_arr - temp1_arr)))
		g0 = g2_arr - self.logg_arr
		g1 = self.logg_arr - g1_arr
		t0 = temp2_arr - self.temp_arr
		t1 = self.temp_arr - temp1_arr
		params_arr = const[np.newaxis, :, np.newaxis] * (\
			f00 * g0[np.newaxis, :, np.newaxis] * t0[np.newaxis, :, np.newaxis] + \
			f10 * g1[np.newaxis, :, np.newaxis] * t0[np.newaxis, :, np.newaxis] + \
			f01 * g0[np.newaxis, :, np.newaxis] * t1[np.newaxis, :, np.newaxis] + \
			f11 * g1[np.newaxis, :, np.newaxis] * t1[np.newaxis, :, np.newaxis])
		return params_arr

	# using the pre-calculated, mapped out features, integrate to find the light at all wavelengths,
	# in ergs/s/Hz/ster per squared equatorial radius of the star
	def intgrt(self):
		# at each wavelength and z value, obtain the integral of the total fit function over phi;
		# to do this, sum up the products of the fit parameters and the corresponding fit integrals
		# along the fit parameter dimension
		fit_arr = np.sum(self.params_arr * self.fitint[np.newaxis, ...], axis=2)
		# at each wavelength, sum up the product of the phi integral of the fit function and
		# the dimensionless area element at each z, multiply by the z step
		return self.z_step * np.sum(self.A_arr[np.newaxis, ...] * fit_arr, axis=1)

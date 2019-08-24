import limbdark.fit as ft
import star.surface as sf
import util as ut
import numpy as np
import math
import sys, time
from scipy.interpolate import Rbf

class Map:
	""" Contains a map of intensity integrals and their wavelength-dependent coefficients 
	(a.k.a. parameters of the corresponding intensity fits), area elements, 
	gravity and temperature values across a set of z values on the surface of a rotating star. 
	Gravity and temperature calculations are based on Espinosa Lara 2011 (EL)."""

	# initialize with the surface shape of the star, the step in z, the constants
	# needed for temperature and gravity unit conversions, wavelengths
	# and an array of intensity parameter fits
	def __init__(self, surf, z_step, add_logg, mult_temp, wl_arr, fit_params):
		z1 = surf.Z1 # an integration bound on the surface
		self.f = surf.f # flatness of the surface
		self.omega = surf.omega # rotational velocity of the star with this surface

		self.z_arr1 = np.arange(-z1, 0, z_step) # negative z values
		self.z_arr2 = np.arange(0, 1, z_step) # positive z values
		self.z_arr = np.concatenate(( self.z_arr1, self.z_arr2 )) # an array of z values for integration

		# a 2D array of integrated fit functions, 
		# one for each z value and each parameter / interval combination of the fit
		self.fitint = np.zeros( (len(self.z_arr), ft.n * ft.Fit.m) )
		## compute the integrated fit functions
		c = 0
		for z in self.z_arr: 
			if z < z1:
				phi1 = surf.phi1(z)
			else:
				phi1 = math.pi
			a, b = surf.ab(z)
			self.fitint[c] = ft.Fit.integrate(phi1, a, b)
			c += 1

		# an array of area elements for integration, one for each value of z
		self.A_arr = np.array([ surf.A(z) for z in self.z_arr ])

		# arrays of the cylindrical coordinate r and the spherical coordinate rho for each z
		self.r_arr1 = np.array([ surf.R(z) for z in self.z_arr1 ])
		self.r_arr2 = np.array([ surf.R(z) for z in self.z_arr2 ])
		self.r_arr = np.concatenate(( self.r_arr1, self.r_arr2 ))
		r_arr_sq = np.power(self.r_arr, 2)
		self.rho_arr1 = np.array([ surf.rho(z) for z in self.z_arr1 ])
		self.rho_arr2 = np.array([ surf.rho(z) for z in self.z_arr2 ])
		self.rho_arr = np.concatenate(( self.rho_arr1, self.rho_arr2 ))
		
		omega = self.omega

		# effective gravitational acceleration in units of G M / Re**2 
		# as in equations 7 and 31 of EL, for each z
		geff_arr = np.sqrt( self.rho_arr**-4 + omega**4 * r_arr_sq - \
			2 * omega**2 * r_arr_sq * self.rho_arr**-3 )
		# log10(gravity in cm / s**2)
		self.logg_arr = add_logg + np.log10(geff_arr)
		
		# effective temperature in units of ( L / (4 pi sigma Re**2) )**(1/4), 
		# as in EL eqn 31, for each z
		t_arr = geff_arr**(1./4) * self.Tc()
		# temperature in Kelvin
		self.temp_arr = mult_temp * t_arr

		# make sure that we are not extrapolating: that all combinations of temperature and gravity
		# on the star are between some four combinations in Kurucz's files
		# np.searchsorted(temp_arr, )

		# for each wavelength and intensity fit parameter, interpolate the parameter
		# as a function of temperature and gravity, then calculate the parameter for each z;
		# results in an array indexed by the wavelength, value of z and parameter index
		n_params = ft.Fit.m * ft.n
		n_wl = len(wl_arr)
		self.params_arr = np.empty( (n_wl, len(self.z_arr), n_params) )
		print ("Interpolating to find the parameters of the fits at all values of z. ")
		sys.stdout.flush()
		start = time.time()
		for ind_wl in range(n_wl):
			if (ind_wl % 100 == 0):
				ut.printf("interpolation for " + str(ind_wl) + " out of " +\
					str(n_wl) + " wavelengths completed.\n")
			sys.stdout.flush()
			for ind_p in range(n_params):
				temp = fit_params[ind_wl, ind_p, :, 0] # temperatures for this parameter and wavelength
				logg = fit_params[ind_wl, ind_p, :, 1] # log g for this parameter and wavelength
				p = fit_params[ind_wl, ind_p, :, 2] # parameter values for this parameter and wavelength
				func = Rbf(temp, logg, p, function='cubic')
				self.params_arr[ind_wl, :, ind_p] = np.array( func(self.temp_arr, self.logg_arr) )
				# [func(t, g)[0] for t, g in zip(self.temp_arr, self.logg_arr)]
		end = time.time()
		print("Done in " + str(end - start) + " seconds")

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

from pa.lib import surface as sf
import pa.lib.map as mp
from pa.lib import fit as ft
from pa.lib import util as ut
import numpy as np
import sys
import math

class Star:
	""" Contains all the information pertaining to a rotating star, including its surface shape,
	a map of physical and geometrical features across its surface, its size, mass and luminosity. 	
	Also performs the 1D integration in the z dimension,
	based on an open-interval formula, equation 4.1.18 of Numerical Recipes, 3rd edition."""
	def __init__(self, omega, luminosity, mass, Req, n_z, ld, temp_method='linear'):
		self.wavelengths = ld.wl_arr # wavelengths
		self.bounds = ld.bounds # the bounds between mu intervals in intensity fits
		self.luminosity = luminosity
		self.mass = mass
		self.Req = Req # equatorial radius, in solar radii
		# information about the surface shape of the star and its inclination
		self.surface = sf.Surface(omega)
		# an additive constant for log g and a multiplicative constant for Teff
		add_logg = math.log10(ut.G * ut.Msun / ut.Rsun**2) + math.log10(mass) - 2 * math.log10(Req)
		mult_temp = ut.Tsun * Req**(-0.5) * luminosity**(0.25)
		# map of gravity, temperature, intensity fit parameters 
		# and other features across the surface of the star
		self.map = mp.Map(self.surface, n_z, add_logg, mult_temp, ld, temp_method)

	# for a given inclination,
	# using the pre-calculated mapped features, integrate to find the light at all wavelengths,
	# in ergs/s/Hz/ster;
	# uses the integration scheme from Numerical Recipes and an additive
	# correction to the scheme that is based on the assumption that the integrand vanishes
	# at the lower integration bound and that the integrand is linear on the interval between
	# the two successive z values around the lower integration bound.
	# also allows a modified Riemann sum and a modified trapezoidal rule.
	def integrate(self, inclination, method='quadratic'):
		# fetch the surface and the map
		surf = self.surface
		mapp = self.map
		z_arr = mapp.z_arr
		dz = mapp.dz
		# set the inclination of the surface
		surf.set_inclination(inclination)
		# get the integration bound for this surface and inclination
		z1 = surf.z1
		
		# produce a mask that says which values of z are strictly above the appropriate integration bound  
		mask = (mapp.z_arr > -z1)
		z = mapp.z_arr[mask] # array of z that are strictly greater than the integration bound
		## compute a 2D array of fit function integrals, 
		## one for each combination of z value, interval and fit function
		a, b = surf.ab(z) # arrays of coefficients needed for the computation of the integrals 
		belowZ1 = np.array( (z < z1) ) # whether only part of the surface at a given z is visible 
		ft.Fit.set_muB(self.bounds) # set the bounds between mu intervals in intensity fits
		fitint = ft.Fit.integrate(belowZ1, a, b)
		# at each z value and wavelength, obtain the integral of the total fit function over phi;
		# to do this, sum up the products of the fit parameters and the corresponding fit integrals
		# along the fit parameter dimension
		fit_arr = np.sum(mapp.params_arr[mask, :, :] * fitint[:, np.newaxis, :], axis=2)
		# obtain the integrand to integrate in the z dimension:
		# at each z and each wavelength, obtain the product of the phi integral of the fit function and
		# the dimensionless area element
		f = mapp.A_arr[mask, np.newaxis] * fit_arr
		## initialize the numerical scheme weights for the set of z values 
		## involved in the integration scheme
		weights = np.ones(f.shape[0], dtype=np.float) # initialize all weights to 1
		# initialize a correction due to the location of the lower integration bound
		# w.r.t. the z values
		corr = 0 

		# depending on the integration method, possibly modify the array of the integrand and set the
		# integration weights
		if method == 'quadratic':
			# find the index of the smallest element of the z values array that is 
			# strictly greater than the lower integration bound
			i = np.searchsorted(z_arr, -z1, side='right')
			# find the difference between this element and the lower integration bound
			d = z_arr[i] + z1
			## based on this difference, possibly modify the set of integrand values involved
			## in the Numerical Recipes integration scheme and compute the correction to the scheme 
			# coefficients of the quadratic approximating the integrand
			a = ( -f[0] / d 			+ f[1] / (d + dz) ) 	/ dz
			b = ( f[0] * (d + dz) / d 	- f[1] * d / (d + dz) ) / dz
			if d >= dz / 2: # if the difference is larger than delta-z
				# correct by the area to the left of the lower integration bound
				d1 = dz - d
				corr = (a / 3) * d1**3 - (b / 2) * d1**2
			else: # the difference should be above zero
				# correct by the area to the right of the lower integration bound
				corr = (a / 3) * d**3 + (b / 2) * d**2
				f = f[1:] # remove the first integrand value from the integration scheme
				weights = weights[1:] # do the same with the weights array
			# set the weights at the boundaries of the open integration interval
			weights[ [ 0,  1,  2] ] = [55./24, -1./6, 11./8]
			weights[ [-1, -2, -3] ]	= [55./24, -1./6, 11./8]
		elif method == 'riemann':
			pass # weights remain ones
		elif method == 'trapezoid':
			weights[-1] = 0.5 # set the upper boundary weight to 1/2
		# sum up the product of the integrand and the weights, add the correction
		result = dz * np.sum(weights[:, np.newaxis] * f, axis=0) + corr
		# return the result in the appropriate units
		return result * (self.Req * ut.Rsun)**2

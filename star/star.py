import star.surface as sf
import star.map as mp
import limbdark.fit as ft
import util as ut
import numpy as np
import sys
import math

class Star:
	""" Contains all the information pertaining to a rotating star, including its surface shape,
	a map of physical and geometrical features across its surface,
	its size, mass and luminosity. """
	def __init__(self, omega, luminosity, mass, Req, z_step, ld):
		self.wavelengths = ld.wl_arr # wavelengths
		self.bounds = ld.bounds # the bounds between mu intervals
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
		self.map = mp.Map(self.surface, z_step, add_logg, mult_temp, ld)


	# for a given inclination,
	# using the pre-calculated mapped features, integrate to find the light at all wavelengths,
	# in ergs/s/Hz/ster
	def integrate(self, inclination):
		# fetch the surface and the map
		surf = self.surface
		mapp = self.map
		# set the inclination of the surface
		surf.set_inclination(inclination)
		# get the integration bound for this surface and inclination
		z1 = surf.z1
		# produce a mask that says which values of z are equal to or 
		# are above the appropriate integration bound
		# (at this lower integration bound the integrand evaluates to zero)
		mask = (mapp.z_arr > -z1)
		z = mapp.z_arr[mask] # array of z at or above the integration bound
		## now that the inclination and the corresponding mask on z values are set,
		## compute a 2D array of fit function integrals, 
		## one for each combination of z value, interval and fit function
		# calculate the integrals
		a, b = surf.ab(z)
		belowZ1 = np.array( (z < z1) )
		ft.Fit.set_muB(self.bounds) # set the bounds between mu intervals in the Fit class
		fitint = ft.Fit.integrate(belowZ1, a, b)

		# at selected z values and each wavelength, obtain the integral of the total fit function over phi;
		# to do this, sum up the products of the fit parameters and the corresponding fit integrals
		# along the fit parameter dimension
		fit_arr = np.sum(mapp.params_arr[mask, :, :] * fitint[:, np.newaxis, :], axis=2)
		# at each wavelength, sum up the product of the phi integral of the fit function and
		# the dimensionless area element at each z, multiply by the z step
		return mapp.z_step * np.sum(mapp.A_arr[mask, np.newaxis] * fit_arr, axis=0) * (self.Req * ut.Rsun)**2

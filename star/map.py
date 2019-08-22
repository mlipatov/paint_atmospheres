import limbdark.fit as ft
import star.surface as sf
import numpy as np
import math

class Map:
	""" Contains a map of intensity integrals, area elements, gravity and temperature across
	a set of z values on the surface of a rotating star. """

	# initialize with the surface shape of the star and the step in z
	def __init__(self, surf, z_step):
		omega = surf.omega # rotational velocity of the star with this surface
		z1 = surf.Z1 # an integration bound on the surface
		self.z_arr = np.arange(-z1, 1, z_step) # an array of z values for integration
		z_arr = self.z_arr

		# a 2D array of integrated fit functions, 
		# one for each z value and each parameter / interval combination of the fit
		self.fitint = np.zeros( (len(self.z_arr), ft.n * ft.Fit.m) )
		## compute the integrated fit functions
		c = 0
		for z in z_arr: 
			if z < z1:
				phi1 = surf.phi1(z)
			else:
				phi1 = math.pi
			a, b = surf.ab(z)
			self.fitint[c] = ft.Fit.integrate(phi1, a, b)
			c += 1

		# an array of area elements for integration, one for each value of z
		self.A_arr = np.array([ surf.A(z) for z in z_arr ])

		# arrays of the cylindrical coordinate r and the spherical coordinate rho for each z
		r_arr = np.array([ surf.R(z) for z in z_arr ])
		r_arr_sq = np.power(r_arr, 2)
		rho_arr = r_arr = np.array([ surf.rho(z) for z in z_arr ])
		# effective gravitational acceleration in units of G*M/Re^2 
		# as in equations 7 and 31 of EL, for each z
		self.geff_arr = np.sqrt( rho_arr**-4 + omega**4 * r_arr_sq - 2 * omega**2 * r_arr_sq * rho_arr**-3 )

	def curly(self):
		pass
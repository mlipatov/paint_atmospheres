import limbdark.fit as ft
import star.surface as sf
import numpy as np
import math

class Star:
	def __init__(self, omega, inclination, luminosity, mass, Req, z_num):
		self.luminosity = luminosity
		self.mass = mass
		self.Req = Req
		
		self.surface = sf.Surface(omega, inclination)
		surf = self.surface

		z1 = surf.Z1 # an integration bound on the surface
		z_step = (1 + z1) / z_num # step size between values of z for integration
		self.z_arr = np.arange(-z1, 1, z_step) # an array of z values for integration
		z_arr = self.z_arr

		# a 2D array of integrated fit functions, 
		# one for each z value and each parameter / interval combination of the fit
		self.fitint = np.zeros( (len(self.z_arr), ft.n * ft.Fit.m) )
		# compute the integrated fit functions
		c = 0
		for z in z_arr[1:]: # all the integrals should be zero at the first z value
			if z < z1:
				phi1 = surf.phi1(z)
			else:
				phi1 = math.pi
			a, b = surf.ab(z)
			self.fitint[c] = ft.Fit.integrate(phi1, a, b)
			c += 1

		# an array of area elements for integration, one for each value of z
		self.A_arr = np.array([ surf.A(z) for z in z_arr ])
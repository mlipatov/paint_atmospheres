import numpy as np
import math

class Star:
	def __init__(self, omega, inclination, luminosity, mass, Req, surface):
		self.omega = omega
		self.inclination = inclination
		self.luminosity = luminosity
		self.mass = mass
		self.Req = Req
		self.surface = surface
		# derived parameters
		self.w = np.inf
		if omega != 0: 
			self.w = 1 + 2 / omega**2
		self.f = 1 + omega**2 / 2
		self.sini = math.sin(inclination)
		self.cosi = math.cos(inclination)



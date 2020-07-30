import numpy as np

class Grid:
	""" Absolute magnitudes for a grid of stars with sun's equatorial radius"""
	def __init__(self, tau, omega, inc, gamma, Z, av, bands):
		# stellar model parameters
		self.tau = tau
		self.omega = omega
		self.inc = inc
		self.gamma = gamma
		self.Z = Z

		# reddenings and bands
		self.av = av
		self.bands = bands

		# magnitudes, indexed in the above order 
		# use a less memory-intensive data type
		shape = tuple(( len(tau), len(omega), len(inc), len(gamma), len(Z), len(av), len(bands) ))
		self.Mag = np.full( shape, np.nan, dtype=np.float32 )
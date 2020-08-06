import numpy as np

class Grid:
	""" Absolute magnitudes for a grid of stars with sun's equatorial radius"""
	def __init__(self, tau, omega, inc, gamma, Z, av, bands, Mag):
		# stellar model parameters
		self.tau = tau
		self.omega = omega
		self.inc = inc
		self.gamma = gamma
		self.Z = Z

		# reddenings
		self.av = av
		# bands
		self.bands = bands

		# dimensions should be in the above order
		self.Mag = Mag
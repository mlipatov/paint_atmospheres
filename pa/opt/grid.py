class Grid:
	""" Absolute magnitudes for a grid of stars """
	# Stores: 
	#	absolute magnitude array
	# 		0: metallicity
	# 		1: gamma = equatorial radius effective gravity (related to mass)
	# 		2: tau = equatorial radius effective temperature (related to luminosity)
	# 		3: omega
	# 		4: inclination
	# 		5: filter
	#	metallicity list
	#	gamma array
	#	tau array
	#	omega array
	#	inclination array
	#	filter filename array
	#	filter zero-point intensity array
	def __init__(self, Mag, Z, gamma, tau, omega, inc, bands):
		# record all information that should be stored in this object
		self.Mag = Mag
		self.Z = Z
		self.gamma = gamma
		self.tau = tau
		self.omega = omega
		self.inc = inc
		self.bands = bands
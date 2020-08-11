import numpy as np

class Grid:
	# dimension parameters
	dims = ['tau', 'omega', 'inc', 'gamma', 'Z', 'av']
	""" Generic grid of stars with sun's equatorial radius """
	def __init__(self, tau, omega, inc, gamma, Z, av):
		# stellar model parameters
		self.tau = tau
		self.omega = omega
		self.inc = inc
		self.gamma = gamma
		self.Z = Z
		# reddenings
		self.av = av

class MagGrid(Grid):
	dims = Grid.dims + ['bands']

	""" Magnitudes grid of stars with sun's equatorial radius """
	def __init__(self, tau, omega, inc, gamma, Z, av, bands, Mag):
		super().__init__(tau, omega, inc, gamma, Z, av)
		# bands
		self.bands = bands
		# dimensions should be in the superclass constructor order
		self.Mag = Mag

	# Inputs:
	# 	minimum tau
	# 	minimum gamma
	#	maximum A_V
	#	approximate factor by which to reduce each dimension 
	#		(making sure to keep first and last values in each dimension)
	def slice(self, taumin=6000, gammin=3.0, avmax=1.0, n=1):
		# magnitudes
		Mag = self.Mag
		# indices to keep in each dimension
		inds = []
		# sliced grid parameters
		params = []
		# for each dimension name
		for p in self.dims:
			# create a local variable with that name and set it to the grid's values
			x = getattr(self, p)
			i = np.arange(len(x))
			if p == 'tau': # pick out only hot stars
				i = i[x >= taumin]
			elif p == 'gamma': # pick out only stars with high surface gravity
				i = i[x >= gammin]
			elif p == 'av':
				i = i[x <= avmax]
			# for all dimensions except bands and omega,
			# create a set of indices that is reduced by a factor, 
			# making sure the first and last elements will be kept;
			# for bands, keep the original indices
			if p != 'bands' and len(i) > 0:
				i = np.concatenate((i[:-1:n], [i[-1]]))
			# add to the indices for magnitudes
			inds.append(i)
			# select the elements for this dimension
			x = x[i]
			params.append(x)
		# select the elements for the magnitudes array
		Mag = Mag[np.ix_(*inds)]
		return MagGrid(*params, Mag)

class LGrid(Grid):
	def __init__(self, tau, omega, inc, gamma, Z, av, Req, lh):
		super().__init__(tau, omega, inc, gamma, Z, av)
		# dimensions should be in the superclass constructor order
		self.Req = Req
		self.lh = lh

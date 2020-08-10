import numpy as np

# dimension parameters
dims = ['tau', 'omega', 'inc', 'gamma', 'Z', 'av', 'bands']

# Inputs:
#	a grid
# 	minimum tau
# 	minimum gamma
#	approximate factor by which to reduce each dimension 
#		(making sure to keep first and last values in each dimension)
def slice(grid, taumin=6000, gammin=3.0, avmax=1.0, n=1):
	# magnitudes
	Mag = grid.Mag
	# indices to keep in each dimension
	inds = []
	# sliced grid parameters
	params = []
	# for each dimension name
	for p in dims:
		# create a local variable with that name and set it to the grid's values
		x = getattr(grid, p)
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
		if p != 'bands':
			i = np.concatenate((i[:-1:n], [i[-1]]))
		# add to the indices for magnitudes
		inds.append(i)
		# select the elements for this dimension
		x = x[i]
		params.append(x)
	# select the elements for the magnitudes array
	Mag = Mag[np.ix_(*inds)]
	return Grid(*params, Mag)

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
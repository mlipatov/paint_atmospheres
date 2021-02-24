import numpy as np
from scipy import interpolate as interpolate

# In this module, Z is the logarithmic relative metallicity, [M/H]
# All models are at 10 parsecs and with solar equatorial radius

# Return the magnitudes, corrected for distance and radii
# Inputs:
#	An array of absolute magnitudes at solar equatorial radius
#	An array of equatorial radii, with same dimensions as the magnitude grid, minus the bands dimension
#		(or an array broadcastable to such dimensions)
#	An array of distance moduli with the same dimension requirements
# Output:
# 	An array of magnitudes corresponding to the magnitude grid, corrected for distance and radius
def correct(Mag, Req, mod):
	Mag += mod
	Mag -= 5 * np.log10( Req )
	return Mag

# Interpolate in a magnitude grid to get magnitudes at a set of points
# Notes:
# 	When one of the neighbors of a point has NAN magnitudes, that point gets NAN magnitudes in the output;
#	when a point is outside the grid, a ValueError is thrown
# Inputs:
#	A magnitude grid
#	An array of points, e.g. [[tau0, omega0, inc0, ...], [tau1, omega1, inc1, ...], ...]
# Output:
#	An array of magnitudes, e.g. [[F435W_0, F555W_0, F814W_0], [F435W_1, F555W_1, F814W_1], ...]
def interp(mg, xi):
	interp_mag = interpolate.interpn((mg.tau, mg.omega, mg.inc, mg.gamma, mg.Z, mg.av), \
		mg.Mag, xi, method='linear', bounds_error=False, fill_value=np.nan)
	return interp_mag

# A version of the above at a particular metallicity and reddening
# Inputs:
#	A magnitude grid
#	An array of points on a 4D grid, e.g. [[tau0, omega0, inc0, gamma0], [tau1, omega1, inc1, gamma0], ...]
#	Metallicity and reddening
# Output:
#	An array of magnitudes, e.g. [[F435W_0, F555W_0, F814W_0], [F435W_1, F555W_1, F814W_1], ...]
def interp4d(mg, xi, Z, AV):
	# find the index of closest metallicity and reddening
	Zi = np.searchsorted(mg.Z, Z, side='right')
	if (Z - mg.Z[Zi - 1]) <= (mg.Z[Zi] - Z): Zi -= 1
	AVi = np.searchsorted(mg.av, AV, side='right')
	if (AV - mg.av[AVi - 1]) <= (mg.av[AVi] - AV): AVi -= 1
	# interpolate in the remaining dimensions
	interp_mag = interpolate.interpn((mg.tau, mg.omega, mg.inc, mg.gamma), \
		mg.Mag[..., Zi, AVi, :], xi, method='linear', bounds_error=False, fill_value=np.nan)
	return interp_mag

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

class Grid:
	""" Absolute magnitudes for a grid of stars """

	# Inputs:
	# 	model parameters, as a list of tuples
	# 		the first element of the tuple is parameter name (e.g. 'tau')
	#		the secont element is an array of parameter values (e.g. array([0, 1, 2, 3]))
	#	band names, as a list (e.g. ['I', 'V', B'])
	#	magnitudes, as a multi-dimensional array
	#		each dimension except the last corresponds to a model parameter
	#		the size of the dimension should be the same as the number of parameter values
	#		the last dimension corresponds to different bands
	def __init__(self, params, bands, Mag):
		self.params = params
		self.bands = bands
		self.Mag = Mag
import star.surface as sf
import star.map as mp

class Star:
	""" Contains all the information pertaining to a rotating star, including its surface shape,
	a map of physical and geometrical features across its surface,
	its size, mass and luminosity. """

	def __init__(self, omega, inclination, luminosity, mass, Req, z_step):
		self.luminosity = luminosity
		self.mass = mass
		self.Req = Req
		self.surface = sf.Surface(omega, inclination)
		self.map = mp.Map(self.surface, z_step)
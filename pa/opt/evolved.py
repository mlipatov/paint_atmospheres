# a class for storing an array of evolved models, including parameters and observables
class Evolved:
	def __init__(self, Z, gamma, tau, omega, inclination, Mag, filt_files):
		self.Z = Z # metallicities of models
		self.gamma = gamma # model gravity parameters
		self.tau = tau # model temperature parameters
		self.omega = omega # model dimensionless rotation velocities
		self.inclination = inclination # model inclinations
		self.Mag = Mag # model magnitudes in several filters
		self.filt_files = filt_files
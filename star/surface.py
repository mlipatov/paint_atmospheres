import numpy as np
import math

class Surface:
	""" Contains all the information pertaining to the surface of a rotating star,
	as defined by a Roche potential that combines gravitational and rotational effects.
	z is defined as the cylindrical coordinate z normalized by the polar radius, r - 
	the cylindrical coordinate r normalized by the equatorial radius, rho - the 
	spherical coordinate rho normalized by the equatorial radius, f - the flatness of the
	star, i.e. the equatorial radius divided by the polar radius. theta is the 
	spherical polar angle. """

	def __init__(self, omega, inclination):
		self.inclination = inclination
		self.omega = omega
		# derived parameters
		if omega != 0: 
			self.w = 1 + 2 / omega**2
		else:
			self.w = np.inf
		self.f = 1 + omega**2 / 2
		self.sini = math.sin(inclination)
		self.cosi = math.cos(inclination)
		self.Z1 = self.Z1()

	## conversions between z and a related variable
	def U(self, z):
		return z / self.f
	def Z(self, u):
		return u * self.f

	## helper functions for computing r(z) and its derivative 
	def T(self, u):
		w = self.w
		return math.acos(\
			(27. - 2*u**6 - 54*w - 6*u**4*w + 27*w**2 - 6*u**2*w**2 - 2*w**3) / \
			(2.*(u**2 + w)**3))
	def V(self, u):
		return math.cos((1./3) * (2 * math.pi + self.T(u))) 
	# s(u)
	def S(self, u):
		w = self.w
		if np.isinf(w):
			return 1 - u**2
		else:
			return (1. / 3) * (2*w - u**2 + 2 * (u**2 + w) * self.V(u))
	# derivative of s(u)
	def Ds(self, u):
		w = self.w
		if np.isinf(w):
			return -2 * u
		else:
			return 2. * u * (1 - 2 * self.V(u)) / (3 + 6 * self.V(u))

	## r(z) and its derivative
	def R(self, z):
		return math.sqrt(self.S(self.U(z)))
	def Drz(self, z):
		if z == 1:
			return np.NINF
		elif z == -1:
			return np.inf
		else:
			return self.Ds(self.U(z)) / (2. * self.f * math.sqrt(self.S(self.U(z))))

	# input: z
	# output: the differential element of area in the units of equatorial radius squared
	#	times the product of the differential element of phi and
	# 	the differential element of z
	def A(self, z):
		return (1./self.f) * math.sqrt(self.S(self.U(z)) + self.Ds(self.U(z))**2 / 4)

	# coefficients in the expression mu = a * cos(phi) + b, returned as a tuple
	def ab(self, z):
		drz = self.f * self.Drz(z)
		sqt = math.sqrt(1 + drz**2)
		a = self.sini / sqt
		b = - self.cosi * drz / sqt
		return [a, b]

	## integration bound on phi for a given z
	# input: z
	# output: integration bound on phi
	def phi1(self, z): 
		if z == -self.Z1:
			return 0
		else:
		  	return math.acos(self.f * self.Drz(z) * self.cosi / self.sini)

	# spherical coordinate rho as a function of z
	def rho(self, z):
		return math.sqrt( self.R(z)**2 + (z / self.f)**2 )

	## integration bound on z
	## this should be computed only once for a given star
	def Z1(self):
		w = self.w
		sini = self.sini
		cosi = self.cosi
		if np.isinf(w):
			return sini
		elif sini == 0:
			return 0
		else:
			## solve for s at the integration bound
			Tsq = (sini / cosi)**2
			# coefficients of the polynomial
			p = np.array([
				-1 - Tsq, \
				6*(1 + Tsq)*w, \
				-15*(1 + Tsq)*w**2, \
				1 - 2*w + w**2 + 20*w**3 + 4*Tsq*(-1 + 2*w - w**2 + 5*w**3), \
				-(w*(4 - 8*w + 4*w**2 + 15*w**3 + 3*Tsq*(-4 + 8*w - 4*w**2 + 5*w**3))), \
				6*w**2*(1 - 2*w + w**2 + w**3 + Tsq*(-2 + 4*w - 2*w**2 + w**3)), \
				-(Tsq*(-2 + 4*w - 2*w**2 + w**3)**2) - w**3*(4 - 8*w + 4*w**2 + w**3), \
				(-1 + w)**2*w**4
			])
			# roots of the polynomial equation
			rts = np.roots(p)
			# find the root that's between zero and one
			condition = ((0 <= rts) & (rts <= 1)) 
			s = np.real(np.extract(condition, rts)[0])
			# now find the value of u at this s
			u = math.sqrt(-(((-1 + s)*(s + s**2 + (-1 + w)**2 - 2*s*w))/(s - w)**2))
			# finally, obtain the bound on z
			return self.Z(u)
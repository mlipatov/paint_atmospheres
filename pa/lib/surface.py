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

	def __init__(self, omega):
		self.omega = omega
		# derived parameters
		if omega != 0: 
			self.w = 1 + 2 / omega**2
		else:
			self.w = np.inf
		self.f = 1 + omega**2 / 2

	# set the inclination of the star's surface
	def set_inclination(self, inclination):
		self.sini = math.sin(inclination)
		self.cosi = math.cos(inclination)
		self.z1 = self.Z1()

	## conversions between z and a related variable
	def U(self, z):
		return z / self.f
	def Z(self, u):
		return u * self.f

	## helper functions for computing r(z) and its derivative 
	def T(self, u):
		w = self.w
		return np.arccos(\
			(27. - 2*u**6 - 54*w - 6*u**4*w + 27*w**2 - 6*u**2*w**2 - 2*w**3) / \
			(2.*(u**2 + w)**3))
	def V(self, u):
		return np.cos((1./3) * (2 * np.pi + self.T(u))) 
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
		result = np.sqrt(self.S(self.U(z)))
		# just in case we tried to take a square root of a negative number
		# because of limited precision, set the result at z = 1 or -1 to zero
		result[np.abs(z) == 1] = 0
		return result

	def Drz(self, z):
		numerator = self.Ds(self.U(z))
		denominator = (2. * self.f * self.R(z))
		return numerator / denominator

	# input: z
	# output: the differential element of area in the units of equatorial radius squared
	#	times the product of the differential element of phi and
	# 	the differential element of z
	def A(self, z):
		return (1./self.f) * np.sqrt(self.S(self.U(z)) + self.Ds(self.U(z))**2 / 4)

	# coefficients in the expression mu = a * cos(phi) + b, returned as a tuple
	def ab(self, z):
		drz = self.f * self.Drz(z)
		sqt = np.sqrt(1 + drz**2)
		a = self.sini / sqt
		b = - self.cosi * drz / sqt
		return [a, b]

	# spherical coordinate rho, given sets of r and z values
	def rho(self, r, z):
		return np.sqrt( r**2 + (z / self.f)**2 )

	# spherical coordinate rho from a set of theta values
	def rh(self, theta):
		w = self.w
		if np.isinf(w):
			result = 1
		else:
			sine = np.sin(theta)
			result = 2 * np.sqrt(w / 3) * np.sin( \
						(1./3) * np.arcsin( \
							(3 * np.sqrt(3) / 2) * (w - 1) * sine / w**(3./2) ) ) / sine
		return result

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
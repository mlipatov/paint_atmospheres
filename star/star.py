import numpy as np
import math

class Star:
	def __init__(self, omega, inclination, luminosity, mass, Req, surface):
		self.omega = omega
		self.inclination = inclination
		self.luminosity = luminosity
		self.mass = mass
		self.Req = Req
		self.surface = surface
		# derived parameters
		self.w = np.inf
		if omega != 0: 
			self.w = 1 + 2 / omega**2
		self.f = 1 + omega**2 / 2
		self.sini = math.sin(inclination)
		self.cosi = math.cos(inclination)

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
	  return math.acos(self.f * self.Drz(z) * self.cosi / self.sini)

	## integration bound on z
	## this should be computed only once for a given star
	def Z1(self):
		w = self.w
		sini = self.sini
		cosi = self.cosi
		if np.isinf(w):
			return sini
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

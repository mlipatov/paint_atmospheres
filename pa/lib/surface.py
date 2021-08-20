import numpy as np
import math
from scipy import interpolate

# functions that provide the volume and the surface area of a star 
# at a given rotation rate in cubed equatorial radii and squared equatorial radii, respectively;
# the associated calculations use vectorized versions of the object methods further down
V = None
A = None
# function that provides the dimensionless omega w.r.t. Omega_K = sqrt(G M / R_eq^3), given 
# the dimensionless omega w.r.t sqrt(G M / R^3), where R is the volume-averaged radius
omega = None
omin = None
omax = None # maximum value of Omega / sqrt(G M / R^3)

def calcVA():
	def T(w, u):
		return np.arccos(\
			(27. - 2*u**6 - 54*w - 6*u**4*w + 27*w**2 - 6*u**2*w**2 - 2*w**3) / \
			(2.*(u**2 + w)**3))

	def vfunc(w, u):
		return np.cos((1./3) * (2 * np.pi + T(w, u)))

	def S(w, u, v):
		return (1. / 3) * (2*w - u**2 + 2 * (u**2 + w) * v)

	def Ds(w, u, v):
	 	return 2. * u * (1 - 2 * v) / (3 + 6 * v)

	def Afunc(f, s, ds): 
		return (1./f) * np.sqrt(s + ds**2 / 4)

	def R(s): 
		return np.sqrt(s)

	# axis 0 : omega, w, f
	# axis 1: z

	omega = np.linspace(0, 1, 1001)
	w = 1 + 2 / omega[1:]**2 # don't use the zero value
	w = w[:, np.newaxis]
	f = 1 + omega[1:]**2 / 2

	z = np.linspace(0, 1, 1001)
	u = z[np.newaxis, :] / f[:, np.newaxis]
	v = vfunc(w, u)
	s = S(w, u, v)
	s[s < 0] = 0
	ds = Ds(w, u, v)
	a = Afunc(f[:, np.newaxis], s, ds)
	r = R(s)

	Va = (1. / f) * 2 * np.pi * np.trapz(r**2, z, axis=-1)
	Aa = 4 * np.pi * np.trapz(a, z, axis=-1)
	# tack on the values at zero rotation
	Va = np.concatenate( (np.array([4*np.pi/3]), Va) )
	Aa = np.concatenate( (np.array([4*np.pi]), Aa) )

	global V, A
	V = interpolate.interp1d(omega, Va, kind='cubic')
	A = interpolate.interp1d(omega, Aa, kind='cubic')

# calculate the Keplerian dimensionless omega as a function of the MESA dimensionless omega
def calcom():
	def f(omega):
		return omega**2 * V(omega)
	
	global omin, omax
	omin = 0
	omax = np.sqrt(V(1) / (4 * np.pi / 3)) # maximum value for Omega / sqrt(G M / R^3)
	otc = np.linspace(omin, omax, 1001)
	ol = np.zeros(len(otc)) # lower bound
	ou = np.ones(len(otc)) # upper bound
	om = (ol + ou) / 2 # midpoint bracket
	i = 0
	while True:
		fl = f(ol) - (4 * np.pi / 3) * otc**2 # function at the lower bound
		fu = f(ou) - (4 * np.pi / 3) * otc**2 # function at the upper bound
		fm = f(om) - (4 * np.pi / 3) * otc**2 # function at the midpoint
		# signs
		sl = np.sign(fl)
		su = np.sign(fu)
		sm = np.sign(fm)
		# set the upper bound to the midpoint when the signs of lower bound and midpoint are opposite
		mask = (sl * sm == -1)
		ou[mask] = om[mask] # new bracket upper bound
		# set the lower bound to the midpoint when the signs of upper bound and midpoint are opposite
		mask = (su * sm == -1)
		ol[mask] = om[mask] # new bracket lower bound
		# new midpoint
		om = (ol + ou) / 2
		# check for convergence
		if i > 0:
			if np.max(np.abs(om - om_prev)) < 0.0001:
				break
		om_prev = np.copy(om) # remember the midpoint for the next iteration
		i+=1
	om[0] = 0
	om[-1] = 1
	global omega
	omega = interpolate.interp1d(otc, om, kind='cubic')

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

	# set the inclination of the star
	def set_inclination(self, inclination):
		self.inclination = inclination
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
		s = self.S(self.U(z))
		# we may get a negative number at z = 1 or -1,
		# due to limited precision; set r = 0 there
		if np.isscalar(z):
			if np.abs(z) == 1:
				s = 0
		else:
			s[ np.abs(z) == 1 ] = 0
		# square root
		output = np.sqrt(s)
		return output

	def Drz(self, z):
		numerator = self.Ds(self.U(z))
		denominator = (2. * self.f * self.R(z))
		return numerator / denominator

	# output: the differential element of area in the units of equatorial radius squared
	#	times the product of the differential element of phi and
	# 	the differential element of z
	# input: z
	def A(self, z):
		return (1./self.f) * np.sqrt(self.S(self.U(z)) + self.Ds(self.U(z))**2 / 4)

	# coefficients in the expression mu = a * cos(phi) + b
	def ab(self, z):
		# mu = cos(i) at the poles
		m = (np.abs(z) != 1)
		a = np.zeros_like(z)
		b = np.full_like(z, self.cosi)
		# calculate coefficients of mu's dependence everywhere except the poles
		drz = self.f * self.Drz(z[m])
		sqt = np.sqrt(1 + drz**2)
		a[m] = self.sini / sqt
		b[m] = - self.cosi * drz / sqt
		
		return [a, b]

	# output: phi corresponding to mu = 0 at a given z
	# inputs: z
	def cos_phi_b(self, z):
		if self.sini == 0:
			# mu = 0 at all phi when z = 0 and no phi when z != 0,
			# so that phi_b is not defined
			return np.nan
		elif z == self.z1:
			return 0
		elif z == -self.z1:
			return 1
		else:
			return self.f * self.Drz(z) * self.cosi / self.sini

	# output: spherical coordinate rho 
	# inputs: sets of r and z values
	def rho(self, r, z):
		return np.sqrt( r**2 + (z / self.f)**2 )

	# output: spherical coordinate rho 
	# input: a set of theta values
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

	## functions related to the projection of the stellar surface
	## onto the view plane, whose coordinates are y and u-prime, both normalized by Req

	# output: u coordinates of the first intersection with the surface for a number of sightlines;
	#		NAN indicates a sightline that doesn't intersect with the surface
	# inputs: u-prime and y coordinates of points where sightlines intersect the viewplane
	# notes: for inclination closer than 0.01 to pi / 2 at standard machine precision, 
	# 		treats inclination as pi / 2
	def sU(self, up_arr, y_arr):
		# u values corresponding to the points on the line
		u_arr = np.full_like(up_arr, np.nan)
		# if inclination is pi / 2, treat separately because cos(i) = 0
		if (self.cosi**6 < 16 * np.finfo(float).eps / 4.) :
			# u and u-prime are the same, so convert u-prime to z
			z_arr = self.Z( up_arr )
			# mask out the values of z that are beyond the star's boundaries
			mask1 = np.abs( z_arr ) <= self.z1 
			# find r for the remaining values
			r_arr = self.R( z_arr[ mask1 ] )
			# mask out the values of y that are greater than r
			mask2 = np.abs( y_arr[ mask1 ] ) <= r_arr
			# since the second indexing of an array refers to a copy produced by the first indexing,
			# we have to do the following on the l.h.s. in order to index only once 
			u_arr[np.nonzero(mask1)[0][mask2]] = up_arr[ mask1 ][ mask2 ]
		else:
			# helper variables for the computation of u
			omega = self.omega
			o2 = omega**2
			o4 = omega**4
			s = 1 / self.cosi
			t = np.tan(self.inclination)
			c2 = np.cos(2 * self.inclination)
			c4 = np.cos(4 * self.inclination)		
			# compute the u value at each point on the line
			for i in range(len(up_arr)):
				up = up_arr[i]
				y = y_arr[i]
				# coefficients of the 6th degree polynomial in u
				p = np.array([
					(o4*t**4*(1 + t**2))/4.,
					-(o4*s*t**3*(2 + 3*t**2)*up)/2.,
					(t**2*(-2*o4 - 2*o4*t**2 - 4*o2*(1 + t**2)))/4. + \
						(t**2*(6*o4*s**2 + 15*o4*s**2*t**2)*up**2)/4. + \
						(t**2*(2*o4 + 3*o4*t**2)*y**2)/4.,
					s*t*(-(o4*s**2) -5*o4*s**2*t**2)*up**3 +\
						up*(s*t*(o4 + 2*o4*t**2 + o2*(2 + 4*t**2)) + s*t*(-o4 - 3*o4*t**2)*y**2),
					1 + o2 + o4/4. + t**2 + o2*t**2 + (o4*t**2)/4. + \
						((o4*s**4)/4. + (15*o4*s**4*t**2)/4.)*up**4 + \
						(-o2 - o4/2. - 2*o2*t**2 - o4*t**2)*y**2 + (o4/4. + (3*o4*t**2)/4.)*y**4 +\
						up**2*(-(o2*s**2) - (o4*s**2)/2. - 6*o2*s**2*t**2 - 3*o4*s**2*t**2 +\
						((o4*s**2)/2. + (9*o4*s**2*t**2)/2.)*y**2),
					(-3*o4*s**5*t*up**5)/2. + up**3*(-(s*(-8*o2*s**2 - 4*o4*s**2)*t)/2. - 3*o4*s**3*t*y**2) +\
						up*(-((4 + 4*o2 + o4)*s*t)/2. - ((-8*o2 - 4*o4)*s*t*y**2)/2. - (3*o4*s*t*y**4)/2.),
					-1 + (o4*s**6*up**6)/4. + ((4 + 4*o2 + o4)*y**2)/4. + ((-4*o2 - 2*o4)*y**4)/4. +\
						(o4*y**6)/4. + up**4*((-4*o2*s**4 - 2*o4*s**4)/4. + (3*o4*s**4*y**2)/4.) +\
						up**2*((4*s**2 + 4*o2*s**2 + o4*s**2)/4. + ((-8*o2*s**2 - 4*o4*s**2)*y**2)/4. +\
						(3*o4*s**2*y**4)/4.)
				])
				# roots of the polynomial equation
				rts = np.roots(p)
				# real roots
				u = rts[np.isreal(rts)]
				# r that solves both the surface equation and the projected ellipse equation
				# is one of the two values that solves the latter
				r = np.sqrt(y**2 + (up - u*self.sini)**2/self.cosi**2)
				# filter the combinations of u and r that are within their domains
				ok = (np.abs(u) < 1 / self.f)*(np.abs(r) < 1) != 0
				u = u[ok]
				r = r[ok]
				# if the line of sight intersects with the star's surface
				if u.size == 2:
					# record the larger u value for use in spectrum calculations
					u_arr[i] = np.max(u)
		return u_arr
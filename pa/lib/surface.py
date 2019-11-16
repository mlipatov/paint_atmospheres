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

	# coefficients in the expression mu = a * cos(phi) + b, returned as a tuple
	def ab(self, z):
		drz = self.f * self.Drz(z)
		sqt = np.sqrt(1 + drz**2)
		a = self.sini / sqt
		b = - self.cosi * drz / sqt
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
	## onto the view plane, whose coordinates are y and u-prime both normalized by Req

	# distance between two points on the view plane
	@staticmethod
	def distance(y1, up1, y2, up2):
		return np.sqrt( (y2 - y1)**2 + (up2 - up1)**2 )

	# output: a mask saying which points are within the transit and are no less than 
	#		a radius away from the projected stellar surface, 
	#	the same for points outside the transit,
	#	u coordinates of the surface corresponding to sightlines through the points;
	#		NAN indicates a sightline that doesn't intersect with the surface
	# inputs: projected impact parameter of planet's orbit, normalized by Req
	#	projected inclination of the planet's orbit w.r.t. the star's rotation axis (between 0 and pi/2)
	#	planet's radius in Req
	#	u-prime and y coordinates of points on a straight line through the viewplane,
	#		with first and last sets outside the star
	def transit_locations(self, b, alpha, radius, up_arr, y_arr):
		# output: u coordinates of two points where the line of sight intersects with the surface
		# inputs: y and u-prime coordinates of a line of sight
		#	helper variables
		def U(y, up, hvars):
			o2, o4, s, t, c2, c4 = hvars
			# coefficients of the 6th degree polynomial in u
			p = np.array([
				(o4*t**4*(1 + t**2))/4.,
				-(o4*s*t**3*(2 + 3*t**2)*up)/2.,
				(t**2*(-4*o2*(1 + t**2) + 2*o4*(-1 + 3*s**2*up**2 + y**2) +\
					o4*t**2*(-2 + 15*s**2*up**2 + 3*y**2)))/4.,
				s*t*up*(o2*(2 + 4*t**2) - o4*(-1 + s**2*up**2 + y**2 +\
					t**2*(-2 + 5*s**2*up**2 + 3*y**2))),
				1 + t**2 - o2*(-1 + s**2*up**2 + y**2 + t**2*(-1 + 6*s**2*up**2 + 2*y**2)) +\
					(o4*((-1 + s**2*up**2 + y**2)**2 +\
					t**2*(1 + 15*s**4*up**4 - 4*y**2 + 3*y**4 + 6*s**2*up**2*(-2 + 3*y**2))))/4.,
				-(s*t*up*(4 + o2*(4 - 8*s**2*up**2 - 8*y**2) + \
					o4*(1 + 3*s**4*up**4 - 4*y**2 + 3*y**4 + 2*s**2*up**2*(-2 + 3*y**2))))/2.,
				((-1 + s**2*up**2 + y**2)*(4 - 4*o2*(s**2*up**2 + y**2) + \
					o4*(s**4*up**4 + y**2*(-1 + y**2) + s**2*up**2*(-1 + 2*y**2))))/4.
			])
			# roots of the polynomial equation
			rts = np.roots(p)
			# the point on the visible surface is
			# the largest real root between negative one and one
			condition = ((-1 <= rts) & (rts <= 1) & np.isreal(rts))
			u = np.real(np.extract(condition, rts))
			return u

		# output: |y| corresponding to a sightline tangent to the surface
		# input: y and u-prime coordinates of a sightline where there is no intersection
		#	y and u-prime coordinates of a nearby sighline with an intersection
		def Y(y_n, up_n, y, up):
			# number of points on a grid where to look for the intersection
			# closest to the tangent sightline
			n = 10
			# y and u-prime coordinates of the points
			y_arr = np.linspace(y_n, y, n)
			up_arr = np.linspace(up_n, up, n)
			# go through the points until there is an intersection
			for y, up in zip(y_arr, up_arr):
				u = U(y, up, hvars)
				# if the line of sight intersects with the star's surface
				if u.size > 0:
					break
			# average of the z values between the two intersection points
			z = self.Z( np.average(u) )
			# if this z value is outside the boundaries where mu can be zero,
			# bring it back within the boundaries
			if z < -self.z1:
				z = -self.z1
			elif z > self.z1:
				z = self.z1
			# cos(phi) corresponding to the point where mu = 0 for this z value
			cos_phi = self.cos_phi_b( z )
			# abs(y) corresponding to this point
			output = np.sqrt(1 - cos_phi**2)
			# return
			return output

		# slope and u-prime intersect of the line
		m = np.tan(alpha)
		up0 = b / np.cos(alpha)
		# r corresponding to the upper-most z
		r1 = self.R( self.z1 )
		# u-prime corresponding to the upper-most z
		up1 = self.U(self.z1) * self.sini + r1 * self.cosi

		# the sign of the y coordinate at the first and last points 
		# where the transit line intersects the projection of the surface
		ys_first = -1 # when alpha = 0
		ys_last = 1 # when alpha = 0
		if alpha > 0:
			ys_first = np.sign(-up1 - up0)
			ys_last = np.sign(up1 - up0)
		if alpha < 0:
			ys_first = np.sign( -1 * (up1 - up0) )
			ys_last = np.sign( -1 * (-up1 - up0) )

		# helper variables for the computation of u
		omega = self.omega
		o2 = omega**2
		o4 = omega**4
		s = 1 / self.cosi
		t = np.tan(self.inclination)
		c2 = np.cos(2 * self.inclination)
		c4 = np.cos(4 * self.inclination)
		hvars = [o2, o4, s, t, c2, c4]

		# u values corresponding to the points on the line
		u_arr = np.full_like(up_arr, np.nan)
		# whether the first location where the line of sight
		# through the planet's center intersects with the star's surface
		# hasn't been seen
		first = True
		# compute the u value at each point on the line
		for i in range(len(up_arr)):
			up = up_arr[i]
			y = y_arr[i]
			u = U(y, up, hvars)
			# if the line of sight intersects with the star's surface
			if u.size > 0:
				# record the larger u value for potential use in spectrum calculations
				u_arr[i] = np.max(u)
				# if this is the first intersection
				if first:
					# approximate the point on the transit line where
					# a sightline is tangent to the surface
					if self.inclination == 0:
						y_first = ( -m*up0 - np.sqrt(1 + m**2 - b**2) ) / (1 + m**2)
					else:
						y_first = ys_first * Y(y_arr[i - 1], up_arr[i - 1], y, up, ys_first)
					# note that the first point of intersection has passed
					first = False
				# update the index of the last seen intersection
				i_last = i
		# if at least one line of sight through the planet's center 
		# intersects with the surface
		if not first:
			# approximate the second point 
			# where the transit line intersects with the projected surface
			if self.inclination == 0:
				y_last = ( -m*up0 + np.sqrt(1 + m**2 - b**2) ) / (1 + m**2)
			else:
				y_last = ys_last * Y(y_arr[i_last + 1], up_arr[i_last + 1], y_arr[i_last], up_arr[i_last], ys_last)
			## u-prime values of the transit intersections
			up_first = y_first * m + up0
			up_last  = y_last  * m + up0
			# mask saying which values of u-prime and y arrays are between the first 
			# and the second intersection and are no less than a radius away from 
			# these intersections
			within = np.logical_and.reduce( ( y_arr > y_first, y_arr < y_last, 
				self.distance(y_arr, up_arr, y_first, up_first) >= radius,
				self.distance(y_arr, up_arr, y_last, up_last) >= radius ) )
			# mask saying which values of u-prime and y arrays are outside the first 
			# and the second intersection and are no less than a radius away from 
			# these intersections
			outside = np.logical_or( 
				np.logical_and( y_arr < y_first, 
					self.distance(y_arr, up_arr, y_first, up_first) >= radius ),
				np.logical_and( y_arr > y_last, 
					self.distance(y_arr, up_arr, y_last, up_last) >= radius ) )
		return [within, outside, u_arr]
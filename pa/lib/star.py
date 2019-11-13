from pa.lib import surface as sf
import pa.lib.map as mp # use this form of import for map
from pa.lib import fit as ft
from pa.lib import util as ut

import matplotlib as mpl # for the temperature color bar
from matplotlib import pyplot as plt

import scipy.optimize as optimization

import numpy as np
import sys
import math

class Star:
	""" Contains all the information pertaining to a rotating star, including its surface shape,
	a map of physical and geometrical features across its surface, its size, mass and luminosity. 	
	Performs the 1D integration in the z dimension,
	based on an open-interval formula, equation 4.1.18 of Numerical Recipes, 3rd edition.
	Paints the temperature across the surface of the star."""
	def __init__(self, omega, luminosity, mass, Req, nz, ld=None, temp_method='planck'):
		if ld is not None:
			self.wavelengths = ld.wl_arr # wavelengths
			self.bounds = ld.bounds # the bounds between mu intervals in intensity fits
		self.luminosity = luminosity
		self.mass = mass
		self.Req = Req # equatorial radius, in solar radii
		# information about the surface shape of the star and its inclination
		self.surface = sf.Surface(omega)
		# an additive constant for log g and a multiplicative constant for Teff
		add_logg = math.log10(ut.G * ut.Msun / ut.Rsun**2) + math.log10(mass) - 2 * math.log10(Req)
		mult_temp = ut.Tsun * Req**(-0.5) * luminosity**(0.25)
		# map of gravity, temperature, intensity fit parameters 
		# and other features across the surface of the star
		self.map = mp.Map(self.surface, nz, add_logg, mult_temp, ld, temp_method)

	# for a given inclination,
	# using the pre-calculated mapped features, integrate to find the light at all wavelengths,
	# in ergs/s/Hz/ster;
	# uses the integration scheme from Numerical Recipes and an additive
	# correction to the scheme that is based on the assumption that the integrand vanishes
	# at the lower integration bound and that the integrand is quadratic on the interval between
	# the two successive z values around the lower integration bound.
	# also allows a modified midpoint rule.
	def integrate(self, inclination, method='quadratic'):
		# fetch the surface and the map
		surf = self.surface
		mapp = self.map
		# get the data from the map, don't use z = +/- 1
		z_arr = mapp.z_arr[ 1:-1 ] 
		params_arr = mapp.params_arr[ 1:-1]
		A_arr = mapp.A_arr[ 1:-1 ]
		dz = mapp.dz
		# set the inclination of the surface
		surf.set_inclination(inclination)
		# get the integration bound for this surface and inclination
		z1 = surf.z1
		
		# produce a mask that says which values of z are strictly above the appropriate integration bound  
		mask = (z_arr > -z1)
		z = z_arr[mask] # array of z that are strictly greater than the integration bound
		## compute a 2D array of fit function integrals, 
		## one for each combination of z value, interval and fit function
		a, b = surf.ab(z) # arrays of coefficients needed for the computation of the integrals 
		belowZ1 = np.array( (z < z1) ) # whether only part of the surface at a given z is visible 
		ft.Fit.set_muB(self.bounds) # set the bounds between mu intervals in intensity fits
		fitint = ft.Fit.integrate(belowZ1, a, b)
		# at each z value and wavelength, obtain the integral of the total fit function over phi;
		# to do this, sum up the products of the fit parameters and the corresponding fit integrals
		# along the fit parameter dimension
		fit_arr = np.sum(params_arr[mask, :, :] * fitint[:, np.newaxis, :], axis=2)
		# obtain the integrand to integrate in the z dimension:
		# at each z and each wavelength, obtain the product of the phi integral of the fit function and
		# the dimensionless area element
		f = A_arr[mask, np.newaxis] * fit_arr

		## initialize the numerical scheme weights for the set of z values 
		## involved in the integration scheme
		weights = np.ones(f.shape[0], dtype=np.float) # initialize all weights to 1
		# initialize a correction due to the location of the lower integration bound
		# w.r.t. the z values
		corr = 0 

		# depending on the integration method, possibly modify the array of the integrand and set the
		# integration weights
		if method == 'quadratic':
			# find the index of the smallest element of the z values array that is 
			# strictly greater than the lower integration bound
			i = np.searchsorted(z_arr, -z1, side='right')
			# find the difference between this element and the lower integration bound
			d = z_arr[i] + z1
			## based on this difference, possibly modify the set of integrand values involved
			## in the Numerical Recipes integration scheme and compute the correction to the scheme 
			# coefficients of the quadratic approximating the integrand
			a = ( -f[0] / d 			+ f[1] / (d + dz) ) 	/ dz
			b = ( f[0] * (d + dz) / d 	- f[1] * d / (d + dz) ) / dz
			if d >= dz / 2: # if the difference is larger than delta-z / 2
				# correct by the area to the left of the lower integration bound
				d1 = dz - d
				corr = (a / 3) * d1**3 - (b / 2) * d1**2
			else: # the difference should be above zero
				# correct by the area to the right of the lower integration bound
				corr = (a / 3) * d**3 + (b / 2) * d**2
				f = f[1:] # remove the first integrand value from the integration scheme
				weights = weights[1:] # do the same with the weights array
			# set the weights at the boundaries of the open integration interval
			weights[ [ 0,  1,  2] ] = [55./24, -1./6, 11./8]
			weights[ [-1, -2, -3] ]	= [55./24, -1./6, 11./8]
		elif method == 'midpoint':
			weights[-1] = 1.5 # set the upper boundary weight to 1.5

		# sum up the product of the integrand and the weights, add the correction
		result = dz * np.sum(weights[:, np.newaxis] * f, axis=0) + corr
		# return the result in the appropriate units
		return result * (self.Req * ut.Rsun)**2

	# paint the temperature of the visible surface of the star
	# on a set of axes
	# Inputs: axes, inclination, size of the axes in inches, 
	#	an optional axes where to draw a horizontal color bar
	def plot_temp(self, ax, inclination, size_inches, cax=None):

		sine = np.sin(inclination)
		cosn = np.cos(inclination)

		# extract the z values
		z = self.map.z_arr
		# get the r values
		r = self.surface.R( z )
		# convert the z values to the same scale as the r values
		z = z / self.surface.f
		## thicknesses of the ellipses 
		r_diff = np.diff(r)
		z_diff = np.diff(z)
		# line thickness in units of Req
		# double it, so that successive ellipses draw over half the previous line
		l = 2 * (r_diff * cosn + z_diff * sine)
		## don't draw the ellipse at z = +/-1
		r = r[1:]
		z = z[1:]
		# minor axes of the ellipses
		b = r * cosn
		# ellipses' locations in the direction of their minor axes
		c = z * sine
		# temperatures of the ellipses
		T = self.map.temp_arr[1:]
		T_min = np.min(T)
		T_max = np.max(T)
		T_range = T_max - T_min
		# colors = max_col * (T_max - T) / T_range + (1 - max_col)
		colors = (T - T_min) / T_range

		# size of the image in Req
		size_req = 2 
		# unit conversions
		ipr = size_inches / size_req # inches per Req
		ppi = 72 # points per inch
		ppr = ipr * ppi # points per Req

		## draw the star
		min_col = 0.25 # minimum color in the color map
		cmapBig = mpl.cm.get_cmap('YlOrBr', 512)
		cmap = mpl.colors.ListedColormap(cmapBig(np.linspace(1, min_col, 256)))
		ax.set_aspect(1)
		ax.axis('off')
		ax.set_frame_on(False)
		ax.set_xlim([-1, 1])
		ax.set_ylim([-1, 1])
		for i in range(len(z)):
			ellipse = mpl.patches.Ellipse(xy=(0, c[i]), width=2*r[i], height=2*b[i], edgecolor=cmap(colors[i]),\
				fc=cmap(colors[i]), fill=True, lw=ppr * l[i])
			ax.add_patch(ellipse)
		if cax is not None:
			norm = mpl.colors.Normalize(vmin=T_min/1000, vmax=T_max/1000)
			cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
			cb.set_label(r'Temperature, $10^3$ K')

	# output: intensity in erg/s/Hz/ster from a point on the surface of a star
	#	covered by a planet
	# inputs: inclination of the star's rotation axis,
	#	projected impact parameter of planet's orbit,
	#	projected inclination of the planet's orbit w.r.t. the star's rotation axis
	def transit(self, inclination, b, alpha):

		# z values, normalized by the polar radius
		z_tilde = self.map.z_arr
		# z values normalized by the equatorial radius
		u = z_tilde / self.surface.f

		### parametrize the line of transit
		# the point where the line is closest to the center of the star
		x0 = -b * np.sin(alpha)
		y0 = b * np.cos(alpha)

		# helper variables
		omega = self.surface.omega
		o2 = omega**2
		o4 = omega**4
		s = 1 / np.cos(inclination)
		t = np.tan(inclination)
		c2 = np.cos(2 * inclination)
		c4 = np.cos(4 * inclination)

		# number of time points
		n = 200
		# locations of points along the line;
		# the star cannot be outside one equatorial radius of the line's zero point
		loc = np.linspace(-1, 1, n)
		# x and y values of the points on the line
		x_arr = np.cos(alpha) * loc
		y_arr = np.sin(alpha) * loc
		for x, y in zip(x_arr, y_arr):
			# coefficients of the 6th degree polynomial in u
			p = np.array([
				(o4*t**4*(1 + t**2))/4.,
				-(o4*s*t**3*(2 + 3*t**2)*x)/2.,
				(t**2*(-4*o2*(1 + t**2) + 2*o4*(-1 + 3*s**2*x**2 + y**2) +\
					o4*t**2*(-2 + 15*s**2*x**2 + 3*y**2)))/4.,
				s*t*x*(o2*(2 + 4*t**2) - o4*(-1 + s**2*x**2 + y**2 +\
					t**2*(-2 + 5*s**2*x**2 + 3*y**2))),
				1 + t**2 - o2*(-1 + s**2*x**2 + y**2 + t**2*(-1 + 6*s**2*x**2 + 2*y**2)) +\
					(o4*((-1 + s**2*x**2 + y**2)**2 +\
					t**2*(1 + 15*s**4*x**4 - 4*y**2 + 3*y**4 + 6*s**2*x**2*(-2 + 3*y**2))))/4.,
				-(s*t*x*(4 + o2*(4 - 8*s**2*x**2 - 8*y**2) + \
					o4*(1 + 3*s**4*x**4 - 4*y**2 + 3*y**4 + 2*s**2*x**2*(-2 + 3*y**2))))/2.,
				((-1 + s**2*x**2 + y**2)*(4 - 4*o2*(s**2*x**2 + y**2) + \
					o4*(s**4*x**4 + y**2*(-1 + y**2) + s**2*x**2*(-1 + 2*y**2))))/4.
			])
			# roots of the polynomial equation
			r = np.roots(p)
			# find the root that's between negative one and one
			condition = ((-1 <= r) & (r <= 1) & np.isreal(r)) 
			rts = np.extract(condition, r)
			print(x, y, rts)

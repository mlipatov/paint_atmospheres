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
	def __init__(self, omega, luminosity, mass, Req, nz, ld=None, temp_method='planck', g_method='log'):
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
		self.map = mp.Map(self.surface, nz, add_logg, mult_temp, ld, temp_method, g_method)

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
	# 	on a set of axes
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

	# output: flux in zero flux units from a series of points on the surface of a star
	#	covered by a planet
	# inputs: inclination of the star's rotation axis,
	#	planetary transit information
	#	limb darkening information
	#	filter filter information
	#	number of points in the shadow of a planet where to calculate intensities (1, 7, 19 or >19)
	# requres: -pi / 2 < alpha < pi / 2
	# notes: the computation of intensity for each point is much more time-consuming than
	#	the computation of the point itself
	def transit(self, inclination, distance, tr, ld, filt, ns=7):
		# radius of the planet
		radius = tr.radius
		## notes: we work in the viewplane coordinates 
		##	the star cannot be outside one equatorial radius of the line's zero point
		# set the star's inclination
		self.surface.set_inclination(inclination)
		# locations of the planet center
		yc, upc = tr.locations()
		# locations of points within the shadow
		if ns <= 19: # pack smaller circles within the shadow
			# set the ratio of the shadow's radius to that of a smaller circle
			if ns == 1:		a = 1
			elif ns == 7: 	a = 3
			elif ns == 19: 	a = 1 + np.sqrt(2) + np.sqrt(6)
			# coordinates of the points in the shadow, relative to the shadow's center
			y_pts = (radius / a) * np.array([ 0 ])
			up_pts = (radius / a) * np.array([ 0 ])
			if ns == 7 or ns == 19: # if 7 or 19 circles, calculate at least 6 more points
				sqrt3 = np.sqrt(3) # helper variable
				y_pts2 = (radius / a) * np.array([ 2, -2, 1, -1, 1, -1 ])
				up_pts2 = (radius / a) * np.array([ 0, 0, sqrt3, sqrt3, -sqrt3, -sqrt3 ])
				y_pts = np.concatenate( (y_pts, y_pts2) )
				up_pts = np.concatenate( (up_pts, up_pts2) )
				# if 19 circles, calculate 12 more points
				if ns == 19:
					y_pts2 = (radius / a) * np.array([ 1, -1, 1, -1, \
						1 + sqrt3, -(1 + sqrt3), 1 + sqrt3, -(1 + sqrt3), \
						2 + sqrt3, -(2 + sqrt3), 2 + sqrt3, -(2 + sqrt3) ])
					up_pts2 = (radius / a) * np.array([ 2 + sqrt3, 2 + sqrt3, -(2 + sqrt3), -(2 + sqrt3), \
						1 + sqrt3, 1 + sqrt3, -(1 + sqrt3), -(1 + sqrt3), \
						1, 1, -1, -1 ])
					y_pts = np.concatenate( (y_pts, y_pts2) )
					up_pts = np.concatenate( (up_pts, up_pts2) )
		elif ns > 19:
			# generate a number of points within a unit square centered on the origin that will result in
			# approximately the requested number of points within a circle inscribed in it
			N = np.int(np.ceil( ns * 4 / np.pi ))
			y_pts = np.random.uniform(low=-0.5, high=0.5, size=N )
			up_pts = np.random.uniform(low=-0.5, high=0.5, size=N )
			mask = np.sqrt(y_pts**2 + up_pts**2) < 1
			y_pts = radius * y_pts[ mask ]
			up_pts = radius * up_pts[ mask ]
			ns = len(y_pts)

		# flat arrays of y and u-prime values of the representative locations on the planet shadow
		# these include the center of the shadow and several other points within it
		y_arr = np.tile(y_pts, tr.n) + np.repeat(yc, ns)
		up_arr = np.tile(up_pts, tr.n) + np.repeat(upc, ns)

		# u coordinates of the surface corresponding to sightlines through the points
		u_arr = self.surface.sU(up_arr, y_arr)
		# flux at sightlines
		flux_arr = np.zeros_like(u_arr)
		# a mask saying at which points the sightlines intersect the surface
		mask = ~np.isnan(u_arr) 

		# z, r arrays for the points where intensity is calculated
		z_arr = self.surface.Z( u_arr[ mask ] )
		r_arr = self.surface.R( z_arr )
		# mask saying where x < 0 at the intersection of sightline and surface
		xneg = up_arr[ mask ] > u_arr[ mask ] * self.surface.sini
		# abs(cos(phi)) array for the points where intensity is calculated
		cos_phi = np.sqrt( 1 - (y_arr[ mask ] / r_arr)**2 )
		# multiply by -1 where x < 0
		cos_phi[ xneg ] = -1 * cos_phi[ xneg ]
		# a, b, mu arrays for the points where intensity is calculated
		a_arr, b_arr = self.surface.ab(z_arr)
		mu_arr = a_arr * cos_phi + b_arr	
		## intensity fit parameters
		# 0: point index
		# 1: wavelength index
		# 2: parameter index
		params_arr = self.map.Tp( z_arr, r_arr, ld )[1]
		sh = np.shape(params_arr)
		ft.Fit.set_muB(self.bounds) # set the bounds between mu intervals in intensity fits
		## intensities in ergs/s/Hz/ster
		# 0: point index
		# 1: wavelength index
		I_arr = ft.Fit.I(mu_arr, params_arr) * np.pi * radius**2 * (self.Req * ut.Rsun)**2

		# fluxes
		# 0: point index
		flux_arr[ mask ] = filt.flux(I_arr, ld.wl_arr, distance)
		# reshape the flux array according to belonging to different planet shadows
		# and sum the fluxes from different sightlines for a given shadow
		flux = np.sum( np.reshape(flux_arr, (tr.n, ns)), axis=1 )
		# divide the flux by the number of points in a planet's shadow
		flux = flux / ns 
		# return the fluxes
		return flux

class Transit:
	""" Information pertaining to a planetary transit:
	projected impact parameter of planet's orbit, normalized by Req
	projected inclination of the planet's orbit w.r.t. the star's rotation axis (between 0 and pi/2)
	radius of the planet in Req
	number of time points across two equatorial radii of the star and four planetary radii """
	def __init__(self, b, alpha, radius, n):
		self.b = b
		self.alpha = alpha
		self.radius = radius
		self.n = n
		# helper variables
		self.sina = np.sin(alpha)
		self.cosa = np.cos(alpha)

	# outputs: y coordinates of points along the transit line
	#	up coordinates of these points
	# inputs: self
	def locations(self):
		# locations of planet center along the transit line
		loc = np.linspace(-1 - 2*self.radius, 1 + 2*self.radius, self.n)
		# locations of the planet center
		y  = -self.b * self.sina + self.cosa * loc
		up =  self.b * self.cosa + self.sina * loc
		return [y, up]
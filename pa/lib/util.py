# Generally useful functions. Uses cgs units unless specified otherwise.
import math
import sys
import numpy as np
from scipy.interpolate import interp1d

### physical constants
Lsun = 3.839e33 # solar luminosity in erg/s
G = 6.674e-8 # gravitational constant in cm**3 g**-1 s**-2
Msun = 1.989e33 # mass of the sun in grams
Rsun = 6.9551e10 # solar radius in cm
sigma = 5.6704e-5 # Stefan-Boltzmann constant in erg*cm**-2*s**-1*K**-4
c = 2.99792458e10 # speed of light in cm/s
h = 6.62606885e-27 # planck's constant in erg*s
k = 1.3806504e-16 # Boltzmann constant in erg/K
Tsun = (Lsun / (4*math.pi*sigma*Rsun**2))**(0.25) # temperature of the sun in Kelvins
Zsun = 0.017 # from Grevesse, N., & Sauval, A. J., 1998, Space Sci. Rev., 85, 161 (used by Castelli and Kurucz 2004)

# printf() function from O'Reilly's Python Cookbook
def printf(format, *args):
    sys.stdout.write(format % args)

# inputs: starting and ending time in seconds
# output: a string with the time in hours, minutes and seconds
def timef(atime):
	hours, rem = divmod(atime, 3600)
	minutes, seconds = divmod(rem, 60)
	res = "{:0>2}:{:0>2}:{:05.2f}".format(int(hours), int(minutes), seconds)
	return res

# approximate the bolometric luminosity of a star in erg/s/ster
# 	using the trapezoidal rule
# input: light from the star at many wavelengths in erg/s/ster/Hz
#	wavelengths in nm
def bolometric(light, wl):
	# convert intensity per Hz of frequency to per nm of wavelength 
	# this is the integrand in our integral w.r.t. wavelength
	f = ut.Hz_to_nm(light, wl)
	# calculate the differences between wavelengths in nm
	diff = np.diff(wl)
	## estimate the integral using the trapezoidal rule with variable argument differentials
	# array of averaged differentials
	d = 0.5 * ( np.append(diff, 0) + np.insert(diff, 0, 0) )
	return np.sum(d * f)

### Wavelength / frequency conversions

# inputs: an array of wavlengths in nanometers, 
# output: an array of wavelengths in Angstroms
def color_nm_A(wl_arr):
	return wl_arr * 10

# inputs: an array of wavelengths in nanometers
# output: an array of frequencies in Hz
def color_nm_Hz(wl):
	c_nm = 1.e7 * c # speed of light in nm per second
	return c_nm / wl

### Flux conversions

# inputs: an array of some quantity per Angstrom of wavelength, 
# 	an array of corresponding wavelengths in nanometers
# output: an array of the same quantity per Hertz of frequency
def A_to_Hz(f_arr, wl_arr):
	cA = 1.e8 * c # speed of light in angstroms per second
	wl_A = 1.e1 * wl_arr # wavelengths in angstroms
	return f_arr * wl_A**2 / cA

# inputs: an array of some quantity per Hz of frequency, 
# 	an array of corresponding wavelengths in nanometers
# output: an array of the same quantity per Angstrom of wavelength
def Hz_to_A(f_arr, wl_arr):
	cA = 1.e8 * c # speed of light in angstroms per second
	wl_A = 1.e1 * wl_arr # wavelengths in angstroms
	return f_arr * cA / wl_A**2

# inputs: an array of some quantity per nanometer of wavelength, 
# 	an array of corresponding wavelengths in nanometers
# output: an array of the same quantity per Hertz of frequency
def nm_to_Hz(f_arr, wl_arr):
	c_nm = 1.e7 * c # speed of light in nm per second
	return f_arr * wl_arr**2 / c_nm

# inputs: an array of some quantity per Hz of frequency, 
# 	an array of corresponding wavelengths in nanometers
# output: an array of the same quantity per nm of wavelength
def Hz_to_nm(f_arr, wl_arr):
	c_nm = 1.e7 * c # speed of light in nm per second
	return f_arr * c_nm / wl_arr**2

# input: an array of intensity per square centimeter of photodetector, 
# 	distance to the star in centimeters
# output: an array of intensity per steradian
def cm2_to_ster(I_arr, distance):
	return I_arr * distance**2

# input: an array of intensity per steradian, 
# 	distance to the star in centimeters
# output: an array of intensity per square centimeter of photoreceptor
def ster_to_cm2(I_arr, distance):
	return I_arr / distance**2

## Filter
class Filter:
	""" Information pertaining to a filter """
	# inputs: filter multipliers, filter wavelengths, flux zero point in erg/s/cm2/A
	def __init__(self, filt, wlf, f0):
		self.filt = filt
		self.f0 = f0
		# a cubic spline based on the filter
		self.f = interp1d(wlf, filt, kind='cubic', bounds_error=False, fill_value=0)

	# output: 1D array (by location) of flux through a filter 
	#	in units of the filter's flux zero point
	# inputs: 2D array of fluxes in erg/s/ster/Hz (location x wavelength) or
	#		1D array (wavelength)
	#	wavelengths for the light in nm 
	#	distance to the object
	def flux(self, flux_arr, wll, distance):
		# convert flux from per Hz to per angstrom
		flux_arr = Hz_to_A(flux_arr, wll)
		# convert flux from per steradian to per cm2 of photodetector
		flux_arr = ster_to_cm2(flux_arr, distance)
		# convert wavelength to angstroms
		wll_A = color_nm_A(wll)
		# filter evaluated at the light's wavelengths
		fil = self.f(wll_A) 
		# multiply the flux by the filter function
		integrand = np.multiply(flux_arr, fil[np.newaxis, :])
		# calculate the differences between light's wavelengths in A
		diff = np.diff(wll_A)
		## estimate the integral using the trapezoidal rule with variable argument differentials
		# array of averaged differentials
		d = 0.5 * ( np.append(diff, 0) + np.insert(diff, 0, 0) )
		# approximation of the integral, flux integrated over wavelengths
		flux = np.sum( d[np.newaxis, :] * integrand, axis=1 )
		# zero point of the integrated flux
		integrand = fil * self.f0
		flux_zero = np.sum(d * integrand)
		# output
		output = np.empty_like(flux)
		mask = np.logical_or( flux == 0, np.isnan(flux) )
		output[ mask ] = flux[ mask ]
		output[ ~mask ] = flux[ ~mask ] / flux_zero
		return output

	# output: magnitude of star's light through a filter
	# inputs: same as those of flux() in this module
	def mag(self, light, wll, distance):
		return -2.5 * np.log10( self.flux(light, wll, distance) )
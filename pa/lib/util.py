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
D10 = 3.085678e+19 # ten parsecs in cm
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

# pseudo effective temperature
def tau(L, Req):
	return ( L * Lsun / (4 * np.pi * sigma * (Req * Rsun)**2) )**(1./4)
# luminosity in solar luminosities
def L(tau, Req):
	return 4 * np.pi * (Req * Rsun)**2 * sigma * tau**4 / Lsun

# log pseudo effective gravity
def gamma(M, Req):
	return np.log10( G * M * Msun / (Req * Rsun)**2 )
# mass in solar masses
def M(gamma, Req):
	return 10**gamma * (Req * Rsun)**2 / (G * Msun)

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

# inputs: an array of wavelengths in nanometers
# output: an array of frequencies in Hz
def color_nm_Hz(wl):
	c_nm = 1.e7 * c # speed of light in nm per second
	return c_nm / wl

### intensity / flux conversions

# inputs: an array of some quantity per Angstrom of wavelength, 
# 	an array of corresponding wavelengths in Angstroms
# output: an array of the same quantity per Hertz of frequency
def A_to_Hz(f_arr, wl_arr):
	cA = 1.e8 * c # speed of light in angstroms per second
	return f_arr * wl_arr**2 / cA

# inputs: an array of some quantity per Hz of frequency 
#	(last dimension is wavelength), 
# 	an array of corresponding wavelengths in nanometers
# output: an array of the same quantity per Angstrom of wavelength
def Hz_to_A(x, wl):
	cA = 1.e8 * c # speed of light in angstroms per second
	wl_A = 1.e1 * wl # wavelengths in angstroms
	return x * cA / wl_A**2

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

# integrate intensity (last dimension is wavelength) 
# convolved with the transmission curve, normalize
# by the integral of the transmission curve
# Inputs: intensity array in ergs cm^-2 s^-1 Hz^-1 ster^-1 (wavelength dimension should be last)
#	intensity wavelengths in nm
# 	transmission curve
#	filter wavelengths in nm
# Output: intensities in erg cm^-2 s^-1 nm^-1 ster^-1, with the wavelength dimension filtered out
def filter(I, wll, trans, wlf):
	# a cubic spline based on the filter
	Tfunc = interp1d(wlf, trans, kind='cubic', bounds_error=False, fill_value=0)
	# evaluate the transmission curve at the light's wavelengths
	T = Tfunc(wll)
	# convert intensity from per Hz to per nm
	I = Hz_to_nm(I, wll)
	# calculate the differences between light's wavelengths
	diff = np.diff(wll)
	## estimate the integral using the trapezoidal rule with variable argument differentials
	# array of averaged differentials
	d = 0.5 * ( np.append(diff, 0) + np.insert(diff, 0, 0) )
	# approximation of the integral, flux integrated over wavelengths
	intensity = np.sum( d * I * T, axis=-1 )
	# where intensity is not zero or nan, normalize by the integral of the transmission curve
	output = np.empty_like(intensity)
	mask = np.logical_or( intensity == 0, np.isnan(intensity) )
	output[ mask ] = intensity[ mask ]
	output[ ~mask ] = intensity[ ~mask ] / np.sum( d * T )
	return output
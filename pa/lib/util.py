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

# output: star's light through a filter
# inputs: light from the star in erg/s/cm2/A at many wavelengths
#	wavelengths for the light in A
#	filter 
# 	wavelengths for the filter in A
#	filter's intensity zero point in erg/s/cm2/A
def mag(light, wll, filt, wlf, I0):
	# compute a cubic spline based on the filter
	f = interp1d(wlf, filt, kind='cubic', bounds_error=False, fill_value=0)
	# filter evaluated at the light's wavelengths
	fil = f(wll) 
	# and multiply by light
	integrand = np.multiply(fil, light)
	# calculate the differences between light's wavelengths in A
	diff = np.diff(wll)
	## estimate the integral using the trapezoidal rule with variable argument differentials
	# array of averaged differentials
	d = 0.5 * ( np.append(diff, 0) + np.insert(diff, 0, 0) )
	# approximation of the integral
	flux = np.sum(d * integrand)
	# approximate the flux zero point
	integrand = fil * I0
	flux_zero = np.sum(d * integrand)
	# return the magnitude
	return -2.5 * np.log10(flux / flux_zero)

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

### Intensity conversions

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
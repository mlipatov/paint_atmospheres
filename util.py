# Generally useful functions and classes. Uses cgs units
import star.surface as sf
import math
import sys
import numpy as np
import mpmath as mp
import scipy.interpolate as interp

### physical constants
Lsun = 3.839e33 # solar luminosity in erg/s
G = 6.674e-8 # gravitational constant in cm**3 g**-1 s**-2
Msun = 1.989e33 # mass of the sun in grams
Rsun = 6.9551e10 # solar radius in cm
sigma = 5.6704e-5 # Stefan-Boltzmann constant in erg*cm**-2*s**-1*K**-4
c = 2.99792458e10 # speed of light in cm/s
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

# input: wavelength in nm
# output: frequency in Hz
def wl_to_freq(wl):
	return 1.e7 * c / wl
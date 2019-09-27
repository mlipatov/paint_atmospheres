# Generally useful functions and classes. Uses cgs units
import math
import sys
import numpy as np

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

# not needed
# modification of a function from Python documentation on bisect
# def index(a, x):
#     'Locate the leftmost value exactly equal to x'
#     i = bisect_left(a, x)
#     if i != len(a) and a[i] == x:
#         return i
#     return None # this is the line that was modified

# inputs: an array of intensity per Angstrom of wavelength, 
# 	an array of corresponding wavelengths in nanometers
# output: an array of intensity per Hertz of frequency
def convert_from_A(I_arr, wl_arr):
	cA = 1.e8 * c # speed of light in angstroms per second
	wl_A = 1.e1 * wl_arr # wavelengths in angstroms
	return I_arr * wl_A**2 / cA

# inputs: an array of intensity per nanometer of wavelength, 
# 	an array of corresponding wavelengths in nanometers
# output: an array of intensity per Hertz of frequency
def convert_from_nm(I_arr, wl_arr):
	c_nm = 1.e7 * c # speed of light in nm per second
	return I_arr * wl_arr**2 / c_nm

# input: an array of intensity per square centimeter of photoreceptor, 
# 	distance to the star in centimeters
# output: an array of intensity per steradian
def convert_to_ster(I_arr, distance):
	return I_arr * distance**2
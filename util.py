# Generally useful functions and classes. Uses cgs units
import star.surface as sf
import math
import sys
import numpy as np
import mpmath as mp
import scipy.interpolate as interp

c = 2.99792458e10 # speed of light in cm/s

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
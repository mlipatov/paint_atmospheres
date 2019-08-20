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
	res = "{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds)
	return res

# input: wavelength in nm
# output: frequency in Hz
def wl_to_freq(wl):
	return 1.e7 * c / wl

class Ellip:
	""" Implementation of elliptic integrals involved in calculating the integral of sqrt(mu).
	Specifically, given m, the second argument to the incomplete elliptic integral of the 
	second kind E, computes E( arcsin(1/sqrt(m)), m )"""

	# low and high values of inclination and rotation speed;
	# minimum and maximum values of the argument to the elliptic integral are the corresponding 
	# minimum and maximum obtained from the four surfaces that result from the four 
	# combinations of these parameters
	i_bounds = [ 0.01, math.pi / 2 ] # i != 0 because logm is undefined there
	omega_bounds = [ 0, 0.9 ]

	def __init__(self, m_num, z_step, check):
		# argument of the elliptic integral in terms of the parameters of the dependence of
		# mu on cos(phi)
		def m(a, b):
			return 2 / (1 + b / a)
		### Calculate the elliptic integral on an evenly spaced grid of log(m),
		### which corresponds to evenly spaced values of z
		## find the bounds on log(m)
		# surfaces with extreme values of inclination and omega
		surf_arr = [ sf.Surface(omega, i) for omega in self.omega_bounds for i in self.i_bounds ]
		# lower value of z for each surface
		z_low_arr = [ -surf.Z1 + z_step for surf in surf_arr ]
		# parameters of the dependence of mu on cos(phi) at the lower value of z for each surface
		ab_low_arr = [ surf.ab(z_low) for surf, z_low in zip(surf_arr, z_low_arr) ]
		# maximum values of log(m), which occur at the lower values of z
		logm_max_arr = []
		for a, b in ab_low_arr: 
			logm_max_arr.append(math.log(m(a, b)))
		# maximum value of log(m)
		self.logm_max = max(logm_max_arr) 
		# upper value of z
		z_up = 1 - z_step
		# parameters of the dependence of mu on cos(phi) at maximum value of z for each surface
		ab_up_arr = [ surf.ab(z_up) for surf in surf_arr ]
		# minimum values of log(m)
		logm_min_arr = []
		for a, b in ab_up_arr: 
			logm_min_arr.append(math.log(m(a, b)))
		# minimum value of log(m)
		self.logm_min = min(logm_max_arr) 
		## an array of logm values
		self.logm_step = (self.logm_max - self.logm_min) / m_num
		self.logm_arr = np.arange(self.logm_min, self.logm_max + self.logm_step, self.logm_step)
		# an array of m values
		self.m_arr = np.exp(self.logm_arr)
		## integration
		self.int_arr = np.array([ float(mp.ellipe(math.asin(1/math.sqrt(m)), m).real) for m in self.m_arr ])
		## interpolation
		self.func = interp.interp1d(self.logm_arr, self.int_arr, kind='cubic')
		if check:
			self.check()
		self.z_step = z_step

	# return the interpolated value of the elliptic integral for a given m
	def ellip(m):
		return self.func(m)

	def check(self):
		# stagger the values of logm
		logm_arr_2 = np.arange(self.logm_min + self.logm_step / 2, self.logm_max, self.logm_step)
		# an array of m values
		m_arr_2 = np.exp(logm_arr_2)
		# integrate
		int_arr_2 = np.array([ mp.ellipe(math.asin(1/math.sqrt(m)), m) for m in m_arr_2 ])
		# evaluate the interpolated function on these new values of logm,
		# subtract from exact values and calculate the maximum deviation
		max_dev = np.amax(np.abs(int_arr_2 - self.func(logm_arr_2)))
		printf("Maximum deviation of the interpolated elliptic integral from true value of the integral is" + \
			str(max_dev))
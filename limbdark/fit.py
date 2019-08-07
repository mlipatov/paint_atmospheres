## This module contains the class Fit as well as the constants that this class uses
# The range of mu, [0, 1], is split into several intervals. 
# The model is linear in the parameters of the fit on each interval.
# The constraints on the function: the intensity at mu = 1 is set, 
#		the function is continuous everywhere in its domain 

import numpy as np
import math
import matplotlib.pyplot as plt
import util as ut

## functions on a given interval
# the non-zero functions on a given interval
f = [lambda _: 1, lambda x: x, lambda x: x**2, lambda x: x**3, lambda x: math.sqrt(x)]
# the number of non-zero functions on a given interval
n = len(f)
# zero functions on a given interval
fz = [lambda _: 0] * n

## the values of mu in Castelli and Kurucz 2004, reversed
mu_arr = np.array([0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.])
# the boundaries [x0, x1, ... x(n)] between the intervals of the fit: [0, x0), [x0, x1), ... [x(n), 1]
mu0_arr = np.array([0.05, 0.25])
# the number of intervals
m = len(mu0_arr) + 1
# mu array split into a list of arrays according to the boundaries
mu_arr_split = np.split(mu_arr, np.searchsorted(mu_arr, mu0_arr))

## all the functions on each interval
F = []
for i in range(m):
	F.append(fz * i + f + fz * (m - i - 1))
## the model for the given complete set of mu values
# stack the column obtained from applying to each value of mu all the functions 
# on that value's interval
xl = []
for i in range(m): # interval
	Fi = F[i] # all functions on this interval
	for mu in mu_arr_split[i]: # mu on this interval
		xj = [func(mu) for func in Fi]
		xl.append(np.array(xj))
x = np.array(xl)

# number of steps in mu for monotonicity checks
mu_step = 0.0001
# array of mu values for which the fitted intensity functions should be discretely monotonic
mu_check_arr = np.arange(0, 1, mu_step)
# maximum allowed deviation of the given intensity values from computed fit values,
# divided by the given intensity value at mu = 1
max_dev = 0.001


class Fit:
	""" Class containing intensity vs mu for some temperature and gravity, 
		along with the corresponding fit """

	# initialize with an array of intensities corresponding to the array of mu values in this class
	# and the wavelength, log g and temperature of this fit
	# computes the fit
	def __init__(self, I_arr, wl, g, temp, ch):
		self.I_arr = I_arr
		self.wl = wl
		self.g = g
		self.temp = temp
		# fit
		alpha = np.linalg.lstsq(x, I_arr, rcond=None)[0]
		self.params = alpha
		if ch:
			self.check()

	# intensity vs mu in terms of computed parameters
	# mu has to be in [0, 1]
	def I(self, mu):
		i = np.searchsorted(mu0_arr, mu, side='right') # interval where this mu value is found
		Fi = F[i] # all the functions on this interval
		return sum(self.params * [func(mu) for func in Fi])

	# check that the fit is good
	def check(self):
		vI = np.vectorize(self.I) # vectorized version of the intensity function
		# check the zero value
		I0 = self.I(0)
		if I0 < 0:
			ut.printf ("Intensity is negative (" + str(I0) + ") for mu = 0 at wavelength " + str(self.wl) +\
				", g " + str(self.g) + ", temperature " + str(self.temp) + "\n")		
		# check monotonicity
		I_check_arr = vI(mu_check_arr)
		ch = np.all(np.diff(I_check_arr) >= 0)
		if not ch:
			ut.printf ("Intensity is not monotonically increasing at wavelength " + str(self.wl) +\
				", g " + str(self.g) + ", temperature " + str(self.temp) + "\n")
		# check that the fit is good
		I1 = self.I_arr[-1] # given intensity value at mu = 1
		diff = np.abs(np.amax(vI(mu_arr) - self.I_arr))
		if diff == 0:
			dev = 0
		elif I1 == 0:
			dev = np.inf
		else:
			dev = diff / I1
		if dev > max_dev:
			ut.printf ("Maximum deviation of fitted intensity from given values divided by given intensity " +\
				"at mu = 1 is too high: " + str(dev) + " at wavelength " + str(self.wl) + ", g " +\
				str(self.g) + ", temperature " + str(self.temp) + "\n")


	# plot the information in this object
	def plot(self):
		I_arr = self.I_arr
		wl = self.wl
		g = self.g
		temp = self.temp
		params = self.params
		# construct the label string
		lab = ""
		c = 0 # count
		i = 0 # interval
		for a in params:
			if c % n == 0:
				i += 1
				if i != 1:
					lab = lab + '\n'
				lab = lab + 'a' + str(i) + ': '
			lab = lab + '%.2E '
			c += 1
		# construct the offsets of data from the edges of the plot
		delta_y = np.max(I_arr) - np.min(I_arr)
		offset_y = delta_y * 0.1
		# vectorized version of the intesity vs mu function
		vI = np.vectorize(self.I)

		fig = plt.figure()
		plt.axes().set_ylim([np.min(I_arr) - offset_y, np.max(I_arr) + offset_y])
		plt.scatter(mu_arr, I_arr, marker='o', c='b', s=6)
		plt.plot(mu_arr, vI(mu_arr), 'g--', label=lab % tuple(params))
		plt.title('I vs mu wl='+str(wl)+' T='+str(temp)+' g='+str(g))
		plt.xlabel('mu')
		plt.ylabel('intensity')
		plt.legend()
		fig.savefig('I_vs_mu_wl'+str(wl)+'T'+str(temp)+'g'+str(g)+'.png')
# Adapted from the code in pa.lib.map.F()
# Output: a plot of the relative error in temperature correction for omega in [0, 0.999]

from pa.lib import fit as ft
from pa.lib import util as ut
import numpy as np
from mpmath import mp
import math
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
import time

iodir = '../../'

def rho(x, omega):
	output = np.empty_like(omega)
	m = x == 1
	output[m] = 1. / f(omega[m])
	x = x[~m]
	omega = omega[~m]
	output[~m] = (2*np.sqrt(2 + omega**2) * \
		np.sin(np.arcsin((3*np.sqrt(3 - 3*x**2)*omega) / (2 + omega**2)**1.5)/3.)) / \
		(np.sqrt(3 - 3*x**2)*omega)
	return output

def f(omega):
	return 1 + omega**2 / 2

def F_0(omega):
	return (1 - omega**2) ** (-2./3)

def F_1(omega): 
		return np.exp((2./3) * omega**2 * (1 / f(omega))**3)

# output: smallest value of x for which to compute 
#	using full Newton's method step function; a.k.a. x_b
# inputs: a grid of rotational velocities omega
#	order-one factor of proportionality between the error 
#		in full step function and that in the series expansion
#	resolution of floating point numbers
def X1(omega, k, q):
	# a factor that's a little less than 1
	B = 0.9
	# omega below which the approximation yields values greater than 1
	omega_lim = B * (2./(85*k))**(1./4) * 3**(1./2) * q**(1./4)
	output = np.empty_like(omega)
	mask = (omega < omega_lim)
	output[mask] = 1
	output[~mask] = q**(1./6) * (2./(85*k))**(1/6) * \
			( 3**(1./3) * omega[~mask]**(-2./3) - \
			  3**(-2./3) * (199./255) * omega[~mask]**(4./3) - \
			  3**(-2./3) * (29123./65025) * omega[~mask]**(10./3) )
	# in case this estimate exceeds 1 by a little, bring it back to 1
	output [output > 1] = 1
	return output

# helper function
# inputs: an array of temperature correction values and 
#	an array of x = abs(cos(theta)) values
def G(F, x):
	return np.sqrt(F * (1 - x**2) + x**2)
# output: full Newton's method step function
# inputs: F, x, G, rho and omega
def dF_full(F, x, G, rho, omega):
	mult = -2 * F * G**2
	add1 = (1 - G) / x**2
	add2 = (-1./3) * G * rho**3 * omega**2
	add3 = G * np.log( np.sqrt(F) * (1 + x) / (x + G) ) / x**3
	output = mult * (add1 + add2 + add3)
	return output
# same as above, with higher precision; uses mpmath
def dF_full_prec(F, x, G, rho, omega):
	mult = -2 * F * G**2
	add1 = (1 - G) / x ** 2
	add2 = (-1./3) * G * rho**3 * omega**2
	logarg = np.array([mp.sqrt(y) for y in F]) * (1 + x) / (x + G)
	add3 = G * np.array([mp.log(y) for y in logarg]) / x ** 3 
	output = mult * (add1 + add3 + add2)
	return output
# output: series approximation of Newton's method step function up to third order
# inputs: F, x, G, rho and omega
def dF_approx(F, x, omega):
	# helper variables and arrays
	x2 = x**2
	o2 = omega**2
	output = (2*F)/3. + (2*x2)/5. - F**1.5*x2*(1 - o2) - \
		(F**2.5*(10*(1 - o2)**2 + 3*x2*(-3 + 8*o2)))/(15.*(1 - o2))
	return output

nm = 15 # number of steps to run the double-precision versions of the algorithm
nmp = 20 # number of steps to run the higher precision version
# omega_max = 0.999
delta = np.logspace(-3, 0, num=400, base=10)
omega = np.flip(1 - delta)
# omega = np.linspace(0, omega_max, 200)
o2 = omega**2
F0 = F_0(omega) # F at x = 0
F1 = F_1(omega) # F at x = 1
k = 100 # a parameter for estimating this value of x
q = np.finfo(float).eps # resolution of floating point numbers
# optimal smallest value of x for which to compute using Newton's method
xb = X1(omega, k, q)
# rho at these values of x
rho_b = rho(xb, omega)

# initialize the result arrays (to the half-way point in the possible range)
F_full = np.full_like(omega, (F0 + F1) / 2)
F_approx = np.full_like(omega, (F0 + F1) / 2)
F_etalon = mp.mpf(1) * np.full_like(omega, (F0 + F1) / 2)
# Newton's algorithm using the two variants of double precision 
start = time.time()
for i in range(nm):
	# helper function
	G_full = G(F_full, xb)
	G_approx = G(F_approx, xb)
	# the new values of F at the locations 
	# where we use the full Newton's method step function
	F_full = F_full + dF_full(F_full, xb, G_full, rho_b, omega)
	# the new values of F at the locations 
	# where we use the series expansion of Newton's method step function
	F_approx = F_approx + dF_approx(F_approx, xb, omega)
	# check if we end up outside the bounds on F 
	# and come back into the bounds if we did
	m = (F_full < F1);		F_full[ m ] = F1[ m ]
	m = (F_full > F0);		F_full[ m ] = F0[ m ]
	m = (F_approx < F1); 	F_approx[ m ] = F1[ m ]
	m = (F_approx > F0);	F_approx[ m ] = F0[ m ]
end = time.time()
print('Time for the two sets of double precision evaluations in seconds: ' + str(end - start), flush=True)
# Newton's algorithm using the full expression method with higher precision
xb = mp.mpf(1) * xb
F0 = mp.mpf(1) * F0
F1 = mp.mpf(1) * F1
mp.dps = 100 # number of digits after decimal point in higher precision calculations
start = time.time()
for i in range(nm):
	# helper function
	G_etalon = G(F_etalon, xb)
	# the new values of F at the locations 
	# where we use the etalon Newton's method step function
	F_etalon = F_etalon + dF_full_prec(F_etalon, xb, G_etalon, rho_b, omega)
	# check if we end up outside the bounds on F 
	# and come back into the bounds if we did
	m = (F_etalon < F1);	F_etalon[ m ] = F1[ m ]
	m = (F_etalon > F0);	F_etalon[ m ] = F0[ m ]
	# # uncomment the following four lines to see that the etalon values converge
	# if i > 0: 
	# 	diff = np.abs(F_etalon - F_prev)
	# 	print(i + 1, float(diff.max()))
	# F_prev = np.copy(F_etalon)
end = time.time()
print('Time for the high precision evaluations: ' + str(end - start), flush=True)

dfull = np.abs(F_full/F_etalon - 1).astype(float)
dapprox = np.abs(F_approx/F_etalon - 1).astype(float)
print('k = A*B ' + str(k))
print('Maximum error using full formula: ' + str(dfull.max()))
print('Maximum error using series approximation: ' + str(dapprox.max()))

diff = np.concatenate((dfull, dapprox))
max_diff = np.max(diff)
min_diff = np.min(diff)



plt.rcParams.update({'font.size': 18})
rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('text', usetex=True)

# convergence plot figure
fig = plt.figure()
# axes
ax = plt.axes()

ax.set_yscale('log')
ax.set_xscale('log')
ax.invert_xaxis()
ax.set_ylim(q / 1e3, max_diff * 4)

ax.scatter(1 - omega, dapprox, marker='o', facecolors='none', edgecolors='g', s=6)
ax.scatter(1 - omega, dfull, marker='o', facecolors='b', edgecolors='b', s=6)

om_label = [0, 0.9, 0.99, 0.999]
ax.set_xticks(1 - np.asarray(om_label))
ax.set_xticklabels(['%g' % x for x in om_label])
ax.set_xlim(1.2, 1e-3 * 0.8)

ax.set_xlabel(r'$\omega$')
ax.set_ylabel(r'$\left|\delta F(x_b) \,/\, F(x_b)\right|$', labelpad=5)
fig.savefig(iodir + 'error_F.pdf', dpi=200, bbox_inches='tight')

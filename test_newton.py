import numpy as np
import math
import star.surface as sf
import scipy.optimize as optimization # contains bisect

# helper function
# inputs: squared rotation speed, 
#	an array of rho and an array of x = cos(theta)
def add(o2, rho_arr, x_arr):
	return -x_arr - (1/3) * o2 * rho_arr**3 * x_arr**3

# helper function
# inputs: an array of temperature correction values and 
#	an array of x = cos(theta) values
def G(F_arr, x_arr):
	return np.sqrt(F_arr * (1 - x_arr**2) + x_arr**2)

# the function whose root we are looking for, 
# derived from EL24, with x = cos(theta) and F as defined in EL26
def f(F_arr, x_arr, G_arr, add_arr): 
	result1 = x_arr / G_arr
	result2 = np.log( np.sqrt(F_arr) * (1 + x_arr) / (x_arr + G_arr) )
	return result1 + result2 + add_arr

# wrapper function to pass to scipy optimization routines
def f_wrap(F, x, add):
	return f(F, x, G(F, x), add)

# the derivative of the function we are looking for,
# with respect to the temperature correction F
def Df(F_arr, x_arr, G_arr):
	result = (x_arr / G_arr)**3 / (2 * F_arr)
	return result

# x is defined as cos(theta)
tiny = np.finfo(1.0).resolution # precision of floats
delta = 0.0015 # optimal smallest value of x for precision = 1e-15
omega = 0.99
surf = sf.Surface(omega, math.pi/2)
rho0 = 1 # rho at x = 0
rho1 = 1 / surf.f # rho at x = 1
o2 = omega**2
F0 = (1 - o2 * rho0**3)**(-2/3) # F at x = 0
F1 = math.exp((2/3) * o2 * rho1**3) # F at x = 1
Fi = (F0 + F1) / 2 # initial guess

# z_arr = np.arange(delta*1.1, 1, delta)
z_arr = np.array([delta*1.1, delta*10, 0.5, 1])
r_arr = np.array([ surf.R(z) for z in z_arr ])
rho_arr = surf.rho(r_arr, z_arr)
x_arr = z_arr / (surf.f * rho_arr) # cosine theta, same as x
# an additive factor that doesn't depend on F, for every z
add_arr = add(o2, rho_arr, x_arr) 

# initial estimates of F
F_arr = np.full_like(x_arr, Fi)

# F_exact_arr = np.array([
# 	optimization.newton(f_wrap, Fi, args=(x, a))
# 	for Fi, x, a in zip(F_arr, x_arr, add_arr)
# 	])
# F_exact_arr[-1] = F1
print(x_arr)
print(F1, F0)
# print(F_exact_arr)
print()

for i in range(3):
	# print(i)
	# print(i, F_arr)
	# helper function at the current values of F
	G_arr = G(F_arr, x_arr)
	# the function whose root we are looking for, 
	# evaluated at the current values of F
	f_arr = f(F_arr, x_arr, G_arr, add_arr)
	# print(f_arr)
	# the derivative of the function whose root we are looking for,
	# evaluated at the current values of F
	Df_arr = Df(F_arr, x_arr, G_arr)
	# print(Df_arr)
	# compute the new values of F
	F_arr = F_arr - f_arr / Df_arr
	pass
print(F_arr)
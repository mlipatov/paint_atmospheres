import numpy as np
import math
import star.surface as sf
import scipy.optimize as optimization # contains bisect

def correction(curly, rhs): # equation 24 of EL
    return math.cos(curly) + math.log(math.tan(curly/2)) - rhs
def tau(omega, rho, cosn, tan2): # right-hand side of equation 24 of EL (see below)
    return (1.0/3) * omega**2 * rho**3 * cosn**3 + cosn + math.log(tan2)

tiny = 1e-15
omega = 0.9
## Estimate the optimal number of iterations and the resulting precision of the bisection algorithm 
## that obtains curly theta
# size of the range of curly theta
s = math.pi
# looking of curly theta where theta = pi/2, so that tau = 0
surf = sf.Surface(omega, math.pi/2)
z_arr = np.array([-1, 0, 1, 0.5])
r_arr = np.array([ surf.R(z) for z in z_arr ])
rho_arr = surf.rho(r_arr, z_arr)
cosn_arr = z_arr / (surf.f * rho_arr) # cosine theta
tan2_arr = (rho_arr - z_arr / surf.f) / r_arr # tan(theta / 2)
tau_arr = (1./3) * omega**2 * rho_arr**3 * cosn_arr**3 + cosn_arr + np.log(tan2_arr)
curly_exact_arr = np.array([
	optimization.brentq(correction, tiny, math.pi - tiny, args=(tau(omega, rho, cosn, tan2)))
	for rho, cosn, tan2 in zip(rho_arr, cosn_arr, tan2_arr)
	])
curly_exact_arr[0] = math.pi
curly_exact_arr[1] = math.pi/2
curly_exact_arr[2] = 0
print(z_arr)
print(tau_arr)
print(curly_exact_arr)
print()
# the bisection algorithm is only accurate up to n = 19 iterations, 
# down to precision in curly theta equal to 3e-6
curly_arr = np.zeros_like(curly_exact_arr)
for i in range(54):
	print(i, curly_arr - curly_exact_arr)
	# turn off RuntimeWarning: divide by zero encountered in log;
	# numpy.log evaluates to -inf when curly theta = 0, which is the 
	# the value we want in our array at such values of curly theta
	np.seterr(divide = 'ignore')
	# the function whose root we are looking for, evaluated at the 
	# current values of curly theta 
	f1 = np.where(curly_arr == 0, -np.inf, 
		np.cos(curly_arr) + np.log(np.tan(curly_arr / 2)) - tau_arr)
	# turn on RuntimeWarning: divide by zero encountered in log
	np.seterr(divide = 'warn') 
	# the length of the subintervals into which we are subdividing the range of
	# curly theta at this step of the bisection algorithm
	x = s / 2**(i + 1) 
	# the function evaluated at the sum of the current values of curly theta
	# and the length of a subinterval
	f2 = np.cos(curly_arr + x) + np.log(np.tan((curly_arr + x)/ 2)) - tau_arr
	print('  ', f1)
	print('  ',	f2)
	# compute an array of boolean values for each z; the value says whether the algorithm
	# should add the value of the subinterval to the current estimate of curly theta
	# at a given value of z; then add (or don't add) the value of the subinterval accordingly
	curly_arr = curly_arr + x * np.greater(f1 * f2, 0)
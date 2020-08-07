# Computes magnitudes for stars on the grid
# Requires: a list of limbdarkening (LD) objects at a number of metallicites, each for a number of reddenings
# Result: a grid of magnitudes

from pa.lib import util as ut
from pa.lib import star
import pa.lib.map as mp
from pa.opt import grid as gd

import numpy as np
import glob, os, sys, time
import pickle
from more_itertools import sort_together
import multiprocessing as mlp
import psutil
import gc

iodir = '../../'

# compute magnitudes for stars on the grid
# Inputs:
#	metallicity index
#	luminosity index
#	list of limbdarkening (LD) files
#	lists of luminosities, omegas, inclinations, masses, metallicities, reddenings, bands
# Output: magnitudes indexed by [omega][inclination][mass][band][reddening]
def computemag(z, l, limb, L, tau, omega, inc, M, Z, av, bands):
	# array of magnitudes
	shape = tuple(( len(omega), len(inc), len(M), len(bands), len(av) ))
	# use a less memory intensive data type
	result = np.full( shape, np.nan, dtype=np.float32 )
	for m in np.arange(len(M)):
		for o in np.arange(len(omega)):
			try:
				# pre-calculate inclination-independent quantities
				st = star.Star(omega[o], L[l], M[m], 1, ut.D10, 100, ld=limb)
				for i in np.arange(len(inc)):
					# calculate the magnitudes at this inclination;
					# these should be a 1D array, corresponding to bands and reddenings as follows:
					# (b0, av0), (b0, av1), ..., (b0, avm), (b1, av0), ..., (bn, avm)
					mags = st.integrate(inc[i])
					result[o, i, m, :, :] = mags.reshape(len(bands), len(av))
			except mp.InterpolationError as err:
				pass
	print('.', end='', flush=True)
	return result

sockets = 2 # number of chips on the machine
num_cpu = sockets * psutil.cpu_count(logical = False) # number of cores in each parallelized run

# limb darkening
# 0: Z
# 1: A_V
with open(iodir + 'data/ldlist.pkl', 'rb') as fl:
	ldlist = pickle.load(fl)
# bands
bands = ldlist[0].bands

## model parameters

# metallicities and reddening coefficients
Z = []
for iz in range(len(ldlist)):
	l = ldlist[iz]
	Z.append( l.Z )
	if iz == 0:
		av = l.a_v
Z = np.array(Z)

# tau and gamma, based on T and g in the limb darkening information
i = np.argwhere(Z == 0.0)[0][0]
ld = ldlist[i]
# effective temperatures at the equator, related to luminosity:
t = np.sort(np.copy(ld.T)) 
tref = 4 # base refinement in tau, for tau >= 6250
t = ut.refine(t, tref) # increase the size of the tau grid by about this factor
# for cooler stars, refine the tau grid more
m = t < 6250
t = np.concatenate(( ut.refine(t[m], 2), t[~m] ))
m = t < 4500
tau = np.concatenate(( ut.refine(t[m], 2), t[~m] ))

# log effective gravities at the equator, related to mass
gamma = np.sort(np.copy(ld.g)).astype(np.float32)
# corresponding luminosities and masses of stars (in solar luminosities and masses)
# whose equatorial radius is equal to that of the sun
L = ut.L(tau, 1)
M = ut.M(gamma, 1)

# omegas
omega = np.append( np.linspace(0, 1, 20)[:-1], np.array([0.99]) )
# inclinations
inc = np.linspace(0, np.pi/2, 20)

n = len(av) * len(Z) * len(gamma) * len(L) * len(omega) * len(inc)
print('Dimensions of the grid:')
print('A_V\tZ\tgamma\ttau\tomega\tinclination')
print( len(av), '\t', len(Z), '\t', len(gamma), '\t', len(L), '\t', len(omega), '\t', len(inc) )
print('The grid will contain ' + format(n / 1e6, '.0f') + ' million stellar models')
print('It will take up ' + format(n * len(bands) * 4 / 1e9, '.2f') + ' Gb of hard disk space')
print('At every metallicity, this corresponds to ' + format(len(L), '.0f') + ' luminocities / dots below' )
sys.stdout.flush()

# split the computation into pool runs by metallicity, 
# to reduce total RAM usage during every run;
# this necessitates storing the results of previous runs on hard disk
results = [] # this will be indexed by [metallicity][luminosity][omega][inclination][mass][band][reddening]
# pickle the results
with open(iodir + 'data/temp.pkl', 'wb') as f:
	pickle.dump(results, f)

for z in range(len(Z)): 
	# limb darkening information for this metallicity
	limb = ldlist[z]
	# initialize the pool object
	pool = mlp.Pool( num_cpu )

	start = time.time()
	print('Start Z = ' + str(Z[z]))
	sys.stdout.flush()

	# run the code in parallel; the result is a list of arrays;
	# the list is indexed by [luminosity]; the arrays are indexed by [omega][inclination][mass][band][reddening] 
	result = pool.starmap( computemag, [(z, l, limb, L, tau, omega, inc, M, Z, av, bands) for l in range(len(L))] )
	pool.close() # close the pool
	result = np.array(result) # convert the result of this pool run into numpy array

	# unpickle the total list of results
	with open(iodir + 'data/temp.pkl', 'rb') as f:
		results = pickle.load(f)
	# append to the total list of results
	results.append(result) 
	# pickle the results
	with open(iodir + 'data/temp.pkl', 'wb') as f:
		pickle.dump(results, f)
	# dereference all the results arrays
	del results
	del result
	# garbage collect
	gc.collect() 

	end = time.time()
	print('\nZ = ' + str(Z[z]) + ' finished in ' + ut.timef(end - start))
	sys.stdout.flush()

# unpickle the total array of results
with open(iodir + 'data/temp.pkl', 'rb') as f:
	results = pickle.load(f)
# delete the temporary file
os.remove(iodir + 'data/temp.pkl')
# convert the results to numpy array
results = np.array(results) 
# permute from 	[metallicity][luminosity][omega][inclination][mass][band][reddening]
# to 			[luminosity][omega][inclination][mass][metallicity][reddening][band]
results = np.transpose(results, axes=[1, 2, 3, 4, 0, 6, 5])

# absolute magnitudes array for stars with sun's equatorial radius
# move the reddening dimension from the starting to the penultimate position
grid = gd.Grid( tau, omega, inc, gamma, Z, av, bands, results )

nvalid = np.count_nonzero(~np.isnan(grid.Mag)) / len(bands) 
print(format(nvalid / 1e6, '.0f') + ' millions valid stellar models')

### Pickle the grid
with open(iodir + 'data/mag_grid.pkl', 'wb') as f:
	pickle.dump(grid, f)

from pa.lib import util as ut
from pa.lib import star
import pa.lib.map as mp
import numpy as np
import os, sys, time
import pickle
from more_itertools import sort_together

class Grid:
	""" Absolute magnitudes for a grid of stars """
	# Computes absolute magnitudes in several filters, on a grid of
	#	metallicity
	#	gamma = a gravitational parameter (twice as fine as the gravity grid in LD, 21 values)
	#	tau = a thermal parameter (as fine as the temperature grid in LD, 75 values)
	#	omega = a dimensionless rotational parameter (11 values)
	#	inclination (10 values)
	# Inputs: files with filter information, 
	#	zero-point intensities for filters,
	#	limb darkening (LD) information files,
	#	index of the LD file with temperature and gravity corresponding to tau and gamma
	#	an array of omegas
	#	an array of inclinations
	# Stores: 
	#	absolute magnitude array
	# 		0: metallicity
	# 		1: gamma = equatorial radius effective gravity (related to mass)
	# 		2: tau = equatorial radius effective temperature (related to luminosity)
	# 		3: omega
	# 		4: inclination
	# 		5: filter
	#	metallicity list
	#	gamma array
	#	tau array
	#	omega array
	#	inclination array
	#	filter filename array
	#	filter zero-point intensity array
	def __init__(self, filt_files, filt_I0, ld_files, ld_index, omega, inc):
		
		# ten parsecs in cm
		distance = 3.085678e+19 

		# filter information
		filters = []
		count = 0
		for ffile in filt_files:
			f = open(ffile)
			filt = []
			wl = []
			for line in f:
				if line.strip():
					data = line.split()
					wl.append(float(data[0]))
					filt.append(float(data[1]))
			f.close()
			filt = ut.Filter(filt, wl, filt_I0[count])
			filters.append(filt)
			count += 1
		
		# metallicities relative to solar on log scale
		Z = []
		for f in ld_files:
			# obtain Z in the form -01
			z_str = os.path.splitext(os.path.basename(f))[0].split('_')[1].replace('m', '-').replace('p', '')
			# insert the decimal point and convert to float
			z = float(z_str[:-1] + '.' + z_str[-1])
			Z.append(z)
		Z, ld_files = sort_together([Z, ld_files])

		# temperatures and gravities from limb darkening information
		with open(ld_files[ld_index], 'rb') as f:
			ld = pickle.load(f)
		# equatorial radius effective temperatures and log equatorial effective gravities
		# both sorted from smallest to largest; gravities converted to a more memory-intensive data type
		tau = np.sort(np.copy(ld.temp_arr)) # related to luminosity
		# insert additional values into the gravity array
		g_arr = np.copy(ld.g_arr)
		g_arr2 = (g_arr + 0.25)
		gamma = np.sort(np.concatenate( (g_arr, g_arr2) )).astype(np.float32)
		# corresponding luminosities and masses of stars (in solar luminosities and masses)
		# whose equatorial radius is equal to that of the sun
		L = 4 * np.pi * ut.Rsun**2 * ut.sigma * tau**4 / ut.Lsun
		M = 10**gamma * ut.Rsun**2 / (ut.G * ut.Msun)

		# absolute magnitudes array for stars with sun's equatorial radius
		# use a less memory-intensive data type
		Mag = np.full( (len(Z), len(gamma), len(tau), len(omega), len(inc), len(filters)), np.nan, dtype=np.float32 )
		tot = len(L) * len(omega) # total number of stars (not counting inclinations) at a given mass
		print('Calculating magnitudes on a grid of stars.')
		# compute magnitudes for stars on the grid at each metallicity
		for z in np.arange(len(Z)):
			sys.stdout.flush()
			start = time.time()
			# load the limb darkening information
			with open(ld_files[z], 'rb') as f:
				ld = pickle.load(f)
			print('Computing magnitudes for [M/H] = ' + str(Z[z]))
			for m in np.arange(len(M)):
				interp = 0 # number of stars computed using only interpolation
				extrap = 0 # number of stars where extrapolation was used
				print('    ', end='')
				for l in np.arange(len(L)):
					print('.', end='', flush=True)
					for o in np.arange(len(omega)):
						try:
							# pre-calculate inclination-independent quantities
							st = star.Star(omega[o], L[l], M[m], 1, 200, ld=ld)
							if st.map.extr_info is None:
								interp += 1
							else:
								extrap += 1
							for i in np.arange(len(inc)):
								# calculate the spectrum at this inclination
								light = st.integrate(inc[i])
								# calculate the filter magnitudes
								mags = []
								for filt in filters:
									mags.append( filt.mag(light, ld.wl_arr, distance)[0] )
								Mag[z, m, l, o, i, :] = np.array(mags)
								# print(z, m, l, o, i)
								# print(Z[z], M[m], L[l], omega[o], inc[i])
								# print(gamma[m], tau[l])
								# print(mags)
								# sys.exit()
						except mp.InterpolationError as err:
							pass
				print()
				print('    gamma = ' + str(gamma[m]))
				print('    ' + str(interp) + '/' + str(tot) + ' stars computed with interpolation only.')
				print('    ' + str(extrap) + '/' + str(tot) + ' stars computed with extrapolation.')
				print('    ' + str(tot - extrap - interp) + '/' + str(tot) + ' stars not computed.')
			end = time.time()
			print('Time: for [M/H] = ' + str(Z[z]) + ': ' + str(end - start) + ' seconds.')
			sys.stdout.flush()
		# record all information that should be stored in this object
		self.Mag = Mag
		self.filt_files = filt_files
		self.filt_I0 = np.array(filt_I0)
		self.Z = np.array(Z)
		self.gamma = gamma
		self.tau = tau
		self.omega = omega
		self.inc = inc
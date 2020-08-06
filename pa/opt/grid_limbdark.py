# Given an array of metallicities and an array of reddening coefficients,
# constructs an array of limbdarkening (LD) objects, one for each metallicity and all the reddenings.

from pa.lib import limbdark as ld
from pa.lib import util as ut
import numpy as np
import pickle
import glob
import os
from tqdm import tqdm
from more_itertools import sort_together

iodir = "../../"

# bounds on mu intervals
bounds = [0.1, 0.4]

# A_V values
av = np.linspace(0, 0.5, 20)

# filter files
filtfiles = glob.glob(iodir + "data/filters/*")
filtfiles.sort()

# limb darkening files from Castelli and Kurucz
ldfiles = glob.glob(iodir + "data/limbdark/*.pck*")
# metallicities of these files
Z = []
for f in ldfiles:
	zstr = os.path.basename(f)[1:4].replace('m', '-').replace('p', '')
	Z.append( float(zstr[:-1] + '.' + zstr[-1]) )
Z, ldfiles = sort_together([Z, ldfiles])

lam, _, _ = ld.getlamTg(ldfiles[0]) # wavelengths in the Kurucz files
ut.setlam( lam ) # set the wavelength array for filtering and reddening

ldlist = [] # list of limbdarkening objects 
for i in tqdm(range(len(ldfiles))):
	lam, T, g = ld.getlamTg(ldfiles[i])
	# get the intensity array from this file
	I = ld.getI(ldfiles[i], lam, T, g) 
	# initialize a limb darkening object 
	l = ld.LimbDark(lam, T, g, Z[i]) 
	# filter and redden the object, return band x reddening intensities
	Ib = l.filter(I, filtfiles, a_v=av) 
	# fit polynomials for I(mu)
	l.fit(Ib, bounds) 
	# record the object
	ldlist.append(l)

with open(iodir + 'data/ldlist.pkl', 'wb') as f:
	pickle.dump(ldlist, f)
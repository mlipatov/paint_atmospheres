# calculates spectra of a given star at different inclinations
from pa.utils import limbdark
from pa.utils import fit as ft
from pa.utils import star
from pa.utils import util as ut
import numpy as np
import sys
import time
import argparse
import pickle

def run():
	parser = argparse.ArgumentParser(description="Example: \n" +\
		"calc_spectra \'data/vega.pkl\' \'data/vega.txt\' -i 0 1.5707963267948966 0.01; " +\
		"calc_spectra \'data/vega.pkl\' \'data/vega.txt\' -i 0.08683")
	parser.add_argument("pkl_sfile", help="the pickled star file")
	parser.add_argument("output", help="an output spectrum text file to create")
	parser.add_argument('-i', type=float, nargs='+', help='either a single inclination in radians ' +
		'or a range specified by minimum, maximum and step', required=True)
	args = parser.parse_args()

	## inputs
	txt_sfile = args.output # spectrum text file
	pkl_sfile = args.pkl_sfile # pickled star
	i = args.i # inclinations
	li = len(i)
	if li not in [1, 3]:
		sys.exit("Please specify either a single inclination in radians (one number) " +\
			"or a range specified by minimum, maximum and step (three numbers).")
	elif li == 1:
		inclinations = i
	elif li == 3:
		inclinations = np.arange(*i)
	# unpickle the star
	with open(pkl_sfile, 'rb') as f:
		st = pickle.load(f)
	# get the wavelengths at which we see light from this star
	wl = st.wavelengths

	## write the spectrum of the star in text format
	# create this file if it doesn't exist, open it for writing
	f = open(txt_sfile,'w+') 
	# write the header
	f.write('# omega: ' + str(st.surface.omega) + '\n')
	f.write('# luminosity: ' + str(st.luminosity) + '\n')
	f.write('# mass: ' + str(st.mass) + '\n')
	f.write('# Req: ' + str(st.Req) + '\n')
	f.write('# z resolution: ' + str(st.map.dz) + '\n')
	f.close() # close the file

	# open the file for appending
	f = open(txt_sfile, 'a') 
	# calculate the spectra
	leni = len(inclinations)
	light = np.empty( (leni, len(wl)) )
	for i, inc in np.ndenumerate(inclinations):
		light[i] = st.integrate(inc)
		if i[0] % 10 == 0:		
			print(str(i[0]) + " out of " + str(leni) + " inclinations calculated.")        
			sys.stdout.flush()

	# write the spectra to file
	f.write('\n')
	f.write('# intensity(ergs/s/Hz/ster) \n')
	f.write('# wavelength(nm)\\inclination(rad)') 
	for i, inc in np.ndenumerate(inclinations):
		f.write( ' \t' + str(inc) )
	f.write('\n')
	for j, w in np.ndenumerate(wl):
		f.write( str(w) )
		for i, inc in np.ndenumerate(inclinations):
			f.write('\t %.5E' % light[i, j])
		f.write('\n')
	f.close()
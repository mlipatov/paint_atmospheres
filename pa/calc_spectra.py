# calculates spectra of a given star at different inclinations
from pa.lib import limbdark
from pa.lib import fit as ft
from pa.lib import star
from pa.lib import util as ut
import numpy as np
from numpy.core import defchararray as ch
import sys
import time
import argparse
import pickle
import os

def run():
	parser = argparse.ArgumentParser(description="Example: \n" +\
		"calc_spectra \'data/vega.pkl\' \'data/vega/\' -i 0 1.5707963267948966 0.01; " +\
		"calc_spectra \'data/vega.pkl\' \'data/vega/\' -i 0.08683")
	parser.add_argument("pkl_sfile", help="the pickled star file")
	parser.add_argument("output", help="an output directory")
	parser.add_argument('-i', type=float, nargs='+', help='either a single inclination in radians ' +
		'or a range specified by minimum, maximum and step', required=True)
	args = parser.parse_args()

	## inputs
	pkl_sfile = args.pkl_sfile # pickled star file
	output = args.output # output location
	i = args.i # inclinations
	li = len(i)
	if li not in [1, 3]:
		sys.exit("Please specify either a single inclination in radians (one number) " +\
			"or a range specified by minimum, maximum and step (three numbers).")
	elif li == 1:
		inclinations = i
	elif li == 3:
		inclinations = np.arange(*i)
		# decimal precision of inclination
		prec = np.int( np.ceil( -np.log10( i[2] ) ) )
	leni = len(inclinations)
	# unpickle the star
	with open(pkl_sfile, 'rb') as f:
		st = pickle.load(f)
	# get the wavelengths at which we see light from this star
	wl = st.wavelengths

	## write the spectra of the star in text format
	# create the directory
	if not os.path.exists(output):
		os.mkdir(output)
	# filenames
	if not output.endswith('/'):
		output += '/'
	filename = os.path.splitext(os.path.basename(pkl_sfile))[0]
	ofiles = ch.add(output + filename, np.round(inclinations, decimals=prec).astype(str))
	ofiles = ch.replace(ofiles, '.', '_')
	ofiles = ch.add(ofiles, '.txt')

	for i, ofile in np.ndenumerate(ofiles):
		# message
		if i[0] % 10 == 0:		
			print(str(i[0]) + " out of " + str(leni) + " inclinations calculated.")        
			sys.stdout.flush()
		# current inclination
		inc = inclinations[i] 
		# calculate the spectrum
		light = st.integrate(inc)

		# create this file if it doesn't exist, open it for writing
		f = open(ofile,'w+') 
		# write the header
		f.write('# omega: ' + str(st.surface.omega) + '\n')
		f.write('# luminosity: ' + str(st.luminosity) + '\n')
		f.write('# mass: ' + str(st.mass) + '\n')
		f.write('# Req: ' + str(st.Req) + '\n')
		f.write('# z resolution: ' + str(st.map.dz) + '\n')
		f.write('# inclination(rad): ' + str(inclinations[i]) + '\n')
		# write the spectrum to the file
		f.write('\n')
		f.write('# wavelength(nm)\tintensity(ergs/s/Hz/ster)\n') 
		for j, w in np.ndenumerate(wl):
			f.write( str(w) )
			f.write('\t %.5E' % light[j])
			f.write('\n')
		f.close()
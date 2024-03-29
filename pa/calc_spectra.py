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
	parser = argparse.ArgumentParser(description="Examples: \n" +\
		"calc_spectra data/vega.pkl data/vega/ -i 0.000 1.5707963267948966 150; " +\
		"calc_spectra data/vega.pkl data/vega/ -i 0.088418; " +\
		"calc_spectra data/altair.pkl data/altair/ -i 0.8840; " +\
		"calc_spectra data/achernar.pkl data/achernar/ -i 1.0577")
	parser.add_argument("pkl_sfile", help="the pickled star file")
	parser.add_argument("output", help="the output directory")
	parser.add_argument('-i', type=float, nargs='+', help='either a single inclination in radians ' +
		'or a equally spaced values specified by minimum, maximum and number', required=True)
	parser.add_argument("-m", help="longitudinal integration method: 0=cubic(default), 1=trapezoidal", type=int, \
			default=0)
	args = parser.parse_args()

	## inputs
	pkl_sfile = args.pkl_sfile # pickled star file
	output = args.output # output location
	
	# integration method
	if args.m == 0: 
		m = 'cubic'
	elif args.m == 1:
		m = 'trapezoid'
	else:
		sys.exit("Longitudinal integration method should be either 0 (cubic) or 1 (trapezoidal).")

	# inclinations
	i = args.i 
	li = len(i)
	if li not in [1, 3]:
		sys.exit("Please specify either a single inclination in radians (one number) " +\
			"or a range specified by minimum, maximum and step (three numbers).")
	elif li == 1:
		inclinations = np.array( i )
		# decimal precision of inclination for printout
		prec = 6
	elif li == 3:
		mi, ma, num = i
		inclinations = np.linspace( mi, ma, num=int(num) )
		# decimal precision of inclination for printout
		prec = np.int( np.ceil( -np.log10( (ma - mi) / num ) ) )
	leni = len(inclinations)
	
	# unpickle the star
	with open(pkl_sfile, 'rb') as f:
		st = pickle.load(f)
	# get the wavelengths at which we see light from this star
	wl = st.wavelengths

	## write the spectra of the star in text format
	# create the directory if it doesn't exist
	if not os.path.exists(output):
		os.mkdir(output)
	# filenames
	if not output.endswith('/'):
		output += '/'
	filename = os.path.splitext(os.path.basename(pkl_sfile))[0]
	inc_str = np.array([("%." + str(prec) + "f") % x for x in np.round(inclinations, decimals=prec)])
	ofiles = ch.add(output + filename + '_', ch.replace(inc_str, '.', '_'))
	ofiles = ch.add(ofiles, '.txt')

	for i, ofile in np.ndenumerate(ofiles):
		# message
		if i[0] % 10 == 0:		
			print(str(i[0]) + " out of " + str(leni) + " inclinations calculated.")        
			sys.stdout.flush()
		# current inclination
		inc = inclinations[i] 
		# calculate the spectrum or the magnitudes
		light = st.integrate(inc, method=m)

		# create this file if it doesn't exist, open it for writing
		f = open(ofile,'w+') 
		# write the header
		f.write('# luminosity: ' + str(st.luminosity) + '\n')
		f.write('# omega: ' + str(st.surface.omega) + '\n')
		f.write('# inclination(rad): ' + str(inclinations[i]) + '\n')
		f.write('# mass: ' + str(st.mass) + '\n')
		f.write('# Req: ' + str(st.Req) + '\n')
		f.write('# distance: ' + format(st.distance, '.2e') + ' cm\n')
		f.write('# A_V: ' + format(*(st.a_v), '.2f') + '\n')
		f.write('# number of upper half z values: ' + str(st.map.nz) + '\n')
		# write the spectrum to the file
		f.write('\n')
		if st.bands is None: # spectrum mode
			f.write('# wavelength(nm)\tflux(ergs/s/Hz/ster)\n') 
			for j, w in np.ndenumerate(wl):
				f.write( str(w) )
				f.write('\t %.5E' % light[j])
				f.write('\n')
		else: # photometry mode
			f.write('# filter\twavelength(nm)\tmagnitude\n') 
			for j, w in enumerate(wl):
				f.write( st.bands[j] )
				f.write('\t %.6g' % w )
				f.write('\t %.8f' % light[j])
				f.write('\n')
		f.close()
		
# in case we are running this file as the main program
if __name__ == "__main__":
	run()
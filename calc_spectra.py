# calculates spectra of a given star at different inclinations
import limbdark.limbdark as limbdark
import limbdark.fit as ft
import star.star as star
import util as ut
import numpy as np
import sys
import time
import argparse
import pickle

parser = argparse.ArgumentParser(description="Example: \n" +\
	"python calc_spectra.py \'vega.txt\' \'vega.pkl\' -i 0 1.5707963267948966 0.01; " +\
	"python calc_spectra.py \'vega.txt\' \'vega.pkl\' -i 0.08683")
parser.add_argument("output", help="an output spectrum text file to create")
parser.add_argument("pkl_sfile", help="the pickled star file")
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
	inclinations = np.concatenate(( np.arange(*i), i[-2:-1] ))
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
f.close() # close the file
# open the file for appending
f = open(txt_sfile, 'a') 
# calculate and write the spectra to the file
leni = len(inclinations)
for i, inc in np.ndenumerate(inclinations):
	if i[0] % 10 == 0:		
		print(str(i[0]) + " out of " + str(leni) + " inclinations calculated and printed.")        
		sys.stdout.flush()
	f.write('\n')
	f.write('# inclination(rad): ' + str(inc) + '\n')
	f.write('# wavelength(nm)\tintensity(ergs/s/Hz/ster)\n') 
	light = st.integrate(inc)
	ind = 0
	while (ind < len(wl)):
		f.write(str(wl[ind]) + '\t %.5E\n' % light[ind])
		ind += 1
f.close()
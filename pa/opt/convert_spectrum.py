# Plot a spectrum
from pa.lib import util as ut
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pickle

def run():
	parser = argparse.ArgumentParser(description="Examples: \n" +\
		"convert_spectrum \'data/vega1.dat\' \'data/vega.dat\' -d 2.37e19 -s A -w A,    \
		 convert_spectrum \'data/sun1.dat\' \'data/sun.dat\' -d 1.496e13 -p W -a m**2 -s nm")
	parser.add_argument("ofile", help="output spectrum text file")
	parser.add_argument("ifile", help="input spectrum text file")
	parser.add_argument("-d", help="distance to the star in centimeters \
		(if the intensity is per area of photoreceptor, convert to be per steradian)", 
		type=float)
	parser.add_argument("-p", help="units of power (e.g. W; convert to erg/s)")
	parser.add_argument("-a", help="units of photoreceptor area (e.g. m**2; convert to cm**2)")
	parser.add_argument("-s", help="units of specific intensity denominator (e.g. A or nm; convert to Hz)")
	parser.add_argument("-w", help="units of wavelength (e.g. A; convert to nm)")
	args = parser.parse_args()

	ifile = args.ifile # input spectrum file
	ofile = args.ofile # output spectrum file
	distance = args.d # distance to the star in centimeters
	power = args.p # units of power
	area = args.a # units of photoreceptor area
	specific = args.s # units of specific intensity denominator
	wavelength = args.w # units of wavelength

	# obtain the wavelengths and the intensities from the file
	f = open(ifile)
	wl_arr = []
	I_arr = []
	count = 0;
	for line in f:
	    if (not '#' in line) and line.strip():
	        data = line.split()
	        wl_arr.append(float(data[0]))
	        I_arr.append(float(data[1]))
	f.close()
	I_arr = np.array(I_arr)
	wl_arr = np.array(wl_arr)

	# convert units of wavelength to nm
	if wavelength is not None:
		if wavelength == 'A':
			wl_arr = 1.e-1 * wl_arr
	# convert units of power to erg/s
	if power is not None:
		if power == 'W':
			I_arr = 1.e7 * I_arr
	# convert units of photoreceptor area to cm**2
	if area is not None:
		if area == 'm**2':
			I_arr = 1.e-4 * I_arr
	# convert units of the specific intensity denominator to Hz
	if specific is not None:
		if specific == 'A':
			I_arr = ut.convert_from_A(I_arr, wl_arr)
		elif specific == 'nm':
			I_arr = ut.convert_from_nm(I_arr, wl_arr)
	# convert to intensity per steradian from intensity per unit photoreceptor area
	if distance is not None:
		I_arr = ut.convert_to_ster(I_arr, distance)

	# write to the output file
	f = open(ofile,'w+') 
	# write the header
	f.write('# wavelength(nm)\tintensity(ergs/s/Hz/ster)\n') 
	f.close() # close the file
	f = open(ofile, 'a') # open the file for appending
	# write the spectrum to the file
	ind = 0
	while (ind < len(wl_arr)):
		f.write(str(wl_arr[ind]) + '\t %.5E\n' % I_arr[ind])
		ind += 1
	f.close()
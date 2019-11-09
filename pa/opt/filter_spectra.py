# Plot a spectrum
from pa.lib import util as ut
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pickle
import os

def run():
	parser = argparse.ArgumentParser(description="Examples: \n" +\
		"filter_spectra \'data/Generic_Bessell.V.dat\' 3.619e-9 \'data/vega/\' 2.3694e19 \'data/vega_V.txt\'; " +\
		"filter_spectra \'data/Generic_Bessell.B.dat\' 6.317e-9 \'data/vega/\' 2.3694e19 \'data/vega_B.txt\'")
	parser.add_argument("ffile", help="text file with filter in Angstroms")
	parser.add_argument("zero", help="intensity zero point in erg/cm2/s/A", type=float)
	parser.add_argument("idir", help="directory with spectrum text files in nm and ergs/s/Hz/ster")
	parser.add_argument("dist", help="distance to the star in centimeters", type=float)
	parser.add_argument("ofile", help="output file with inclination and magnitude")
	args = parser.parse_args()

	ffile = args.ffile # input filter file
	idir = args.idir # directory with spectrum files
	ofile = args.ofile # output file
	distance = args.dist # distance to the star in centimeters
	I0 = args.zero # intensity zero point

	# get filter information
	f = open(ffile)
	wlf = []
	filt = []
	for line in f:
		if line.strip():
			data = line.split()
			wlf.append(float(data[0]))
			filt.append(float(data[1]))
	f.close()

	# add slashes to directory names if necessary
	if not idir.endswith('/'):
		idir += '/'

	## get input file names
	# ignore hidden files
	input_filenames = (f.name for f in os.scandir(idir) if not f.name.startswith('.'))
	# add directory names
	infiles = sorted( [idir + i for i in input_filenames] )

	mags = []
	incs = []
	## filter spectra
	for infile in infiles:
		# get spectrum and inclination
		f = open(infile)
		wl_arr = []
		I_arr = []
		for line in f:
			if 'inclination' in line:
				inclination = np.float( line.split()[-1] )
				incs.append(inclination)
			if (not '#' in line) and line.strip():
				data = line.split()
				wl_arr.append(float(data[0]))
				I_arr.append(float(data[1]))
		f.close()
		I_arr = np.array(I_arr)
		wl_arr = np.array(wl_arr)
		# convert intensity from per Hz to per angstrom
		I_arr = ut.Hz_to_A(I_arr, wl_arr)
		# convert intensity from per steradian to per cm2 of photodetector
		I_arr = ut.ster_to_cm2(I_arr, distance)
		# convert wavelength to angstroms
		wl_arr = ut.color_nm_A(wl_arr)
		# compute the magnitude
		mags.append(ut.mag(I_arr, wl_arr, filt, wlf, I0))


	# write to the output file
	f = open(ofile,'w+') 
	# write the header
	f.write('# filter: ' + os.path.basename(ffile) + '\n')
	f.write('# intensity zero point: %.5e \n' % I0)
	f.write('# distance: %.5e cm \n' % distance)
	f.write('# inclination(rad)\tmagnitude\n') 
	# write the spectrum to the file
	ind = 0
	while (ind < len(incs)):
		f.write(str(incs[ind]) + '\t %.5E\n' % mags[ind])
		ind += 1
	f.close() # close the file

# in case we are running this file as the main program
if __name__ == "__main__":
	run()
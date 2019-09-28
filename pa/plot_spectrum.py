# Plot a spectrum
import matplotlib.pyplot as plt
import numpy as np
import argparse

def run():
	parser = argparse.ArgumentParser(description="Examples: \n" +\
		"plot_spectrum \'data/vega.txt\' \'data/vega.png\' 3000,\
		 plot_spectrum \'data/sun.txt\' \'data/sun.png\' 4000")
	parser.add_argument("ifile", help="input: spectrum text file")
	parser.add_argument("ofile", help="output: spectrum .png file")
	parser.add_argument("wl", help="upper cutoff for wavelength, in nanometers", type=int)
	args = parser.parse_args()

	ofile = args.ofile # output file
	ifile = args.ifile # input file
	wl = args.wl # cutoff wavelength

	# obtain the wavelengths and the intensities from the file
	f = open(ifile)
	wl_arr = []
	I_arr = []
	count = 0;
	for line in f:
	    if (not '#' in line) and line.strip():
	        data = line.split()
	        if wl < float(data[0]):
	        	break
	        wl_arr.append(float(data[0]))
	        I_arr.append(float(data[1]))
	f.close()
	I_arr = np.array(I_arr)
	wl_arr = np.array(wl_arr)

	# plot
	delta_y = max(I_arr) - min(I_arr)
	offset_y = delta_y * 0.1
	fig = plt.figure()
	plt.axes().set_ylim([min(I_arr) - offset_y, max(I_arr) + offset_y])
	plt.scatter(wl_arr, I_arr, marker='o', c='b', s=3)
	plt.title('intensity vs wavelength from ' + ifile)
	plt.xlabel('wavelength (nm)')
	plt.ylabel('intensity (erg/s/Hz/ster)')
	fig.savefig(ofile)
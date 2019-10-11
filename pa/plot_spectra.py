import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def run():
	parser = argparse.ArgumentParser(description="Examples: \n" +\
		"plot_spectra \'data/vega/\' \'data/vega_img/\' 3000,\
		 plot_spectra \'data/sun/\' \'data/sun_img/\' 4000")
	parser.add_argument("input", help="input directory with text files")
	parser.add_argument("output", help="output directory with .png files")
	parser.add_argument("wl", help="upper cutoff for wavelength, in nanometers", type=int)
	args = parser.parse_args()

	output_dir = args.output # output directory
	input_dir = args.input # input directory
	wl = args.wl # cutoff wavelength

	# add slashes to directory names if necessary
	if not input_dir.endswith('/'):
		input_dir += '/'
	if not output_dir.endswith('/'):
		output_dir += '/'

	# create the output directory if it doesn't exist
	if not os.path.exists(output_dir):
		os.mkdir(output_dir)

	infiles = sorted( [input_dir + i for i in os.listdir(input_dir)] )
	for infile in infiles:
		# create the output file path
		outfile = output_dir + os.path.splitext(os.path.basename(infile))[0] + '.png'

		f = open(infile)
		wl_arr = []
		I_arr = []
		count = 0;
		for line in f:
			if 'inclination' in line:
				inc = np.round( np.float( line.split()[-1] ), decimals=2 )
			if (not '#' in line) and line.strip():
				data = line.split()
				if wl < float(data[0]):
					break
				wl_arr.append(float(data[0]))
				I_arr.append(float(data[1]))
		f.close()
		I_arr = np.array(I_arr)
		wl_arr = np.array(wl_arr)

		## plot
		if inc == 0: # temporary
			max_I = max(I_arr)
			min_I = min(I_arr)
			delta_y = max(I_arr) - min(I_arr)
			offset_y = delta_y * 0.1
		fig = plt.figure()
		plt.axes().set_ylim([min_I - offset_y, max_I + offset_y])
		plt.scatter(wl_arr, I_arr, marker='o', c='b', s=3)
		plt.title( 'Spectrum at i = %5.2f' % inc )
		plt.xlabel(r'$\lambda$, nm')
		plt.ylabel(r'Intensity,  erg / s $\times$ Hz $\times$ ster')
		fig.savefig(outfile, dpi=200)
		plt.close(fig)
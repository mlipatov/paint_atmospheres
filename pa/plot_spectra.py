from pa.lib import star as st
from matplotlib import pyplot as plt
import numpy as np
import argparse
import os

def run():
	parser = argparse.ArgumentParser(description="Examples: \n" +\
		"plot_spectra \'data/vega/\' \'data/vega_img/\' 3000 -t,\
		 plot_spectra \'data/sun/\' \'data/sun_img/\' 4000")
	parser.add_argument("input", help="input directory with text files")
	parser.add_argument("output", help="output directory with .png files")
	parser.add_argument("wl", help="upper cutoff for wavelength, in nanometers", type=int)
	parser.add_argument("-t", help="include visible surface temperature", action="store_true")
	args = parser.parse_args()

	output_dir = args.output # output directory
	input_dir = args.input # input directory
	wl = args.wl # cutoff wavelength
	temp = args.t # whether to include surface temperature

	# add slashes to directory names if necessary
	if not input_dir.endswith('/'):
		input_dir += '/'
	if not output_dir.endswith('/'):
		output_dir += '/'

	# create the output directory if it doesn't exist
	if not os.path.exists(output_dir):
		os.mkdir(output_dir)

	## get input file names
	# ignore hidden files
	input_filenames = (f.name for f in os.scandir(input_dir) if not f.name.startswith('.'))
	# add directory names
	infiles = sorted( [input_dir + i for i in input_filenames] )

	# if the surface temperature picture is needed, create a star
	# with the parameters found in the first input file
	if temp:
		f = open(infiles[0])
		for line in f:
			if 'omega' in line:
				omega = np.float( line.split()[-1] )
			elif 'luminosity' in line:
				luminosity = np.float( line.split()[-1] )
			elif 'mass' in line:
				mass = np.float( line.split()[-1] )
			elif 'Req' in line:
				Req = np.float( line.split()[-1] )
			if 'z values' in line:
				nz = np.int( line.split()[-1] )
		f.close()
		# create a star without limb darkening information
		star = st.Star(omega, luminosity, mass, Req, nz) 


	## plot spectra
	for infile in infiles:
		# create the output file path
		outfile = output_dir + os.path.splitext(os.path.basename(infile))[0] + '.jpeg'

		f = open(infile)
		wl_arr = []
		I_arr = []
		count = 0;
		for line in f:
			if 'inclination' in line:
				inclination = np.float( line.split()[-1] )
				inc = np.round( inclination, decimals=2 )
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
		# if inc == 0: # uncomment if we want to plot all plots with the same y-range
		max_I = max(I_arr)
		min_I = min(I_arr)
		delta_y = max(I_arr) - min(I_arr)
		offset_y = delta_y * 0.1

		fig = plt.figure()
		ax = plt.axes() 
		ax.set_ylim([min_I - offset_y, max_I + offset_y])
		ax.set_ylim((-2.98405e+18, 3.282455e+19)) # temporary
		ax.scatter(wl_arr, I_arr, marker='o', c='b', s=3)
		ax.set_title( 'Spectrum at i = %5.2f' % inc )
		ax.set_xlabel(r'$\lambda$, nm')
		ax.set_ylabel(r'Intensity,  erg / s $\times$ Hz $\times$ ster')

		## Add an inset showing the temperature of the visible surface
		if temp: 
			height_inset = 0.2 # height of the inset, in figure heights
			width, height = fig.get_size_inches() # size of the current figure in inches
			# the width and height of the inset should be the same in inches
			width_inset = height_inset * width / height # width of the inset in figure widths
			# create the new axes
			ax2 = fig.add_axes([0.6, 0.6, width_inset, height_inset])
			# plot the temperature inset
			star.plot_temp(ax2, inclination, height_inset * height)

		fig.savefig(outfile, dpi=200)
		plt.close(fig)

# in case we are running this file as the main program
if __name__ == "__main__":
	run()
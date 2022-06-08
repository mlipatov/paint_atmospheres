from pa.lib import star as st
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import argparse
import os

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif", 
    "font.serif": "Computer Modern",
    "font.size": 20,
    "figure.figsize": (10, 6)
})
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

parser = argparse.ArgumentParser(description="Examples: \n" +\
	"python plot_bd_inclinations.py \'data/vega/\' \'data/vega_img/\' 3000 -t")
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
			omega = float( line.split()[-1] )
		elif 'luminosity' in line:
			luminosity = float( line.split()[-1] )
		elif 'mass' in line:
			mass = float( line.split()[-1] )
		elif 'Req' in line:
			Req = float( line.split()[-1] )
		if 'z values' in line:
			nz = int( line.split()[-1] )
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
			inclination = float( line.split()[-1] )
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
	ds = 10 # downsampling factor
	if inc == 0: # temporary
		max_I = max(I_arr)
		min_I = min(I_arr)
		delta_y = max(I_arr) - min(I_arr)
		offset_y = delta_y * 0.1
	fig = plt.figure()
	plt.axes().set_ylim([min_I - offset_y, max_I + offset_y])
	# plt.scatter(wl_arr, I_arr, marker='o', c='b', s=3)
	plt.plot(wl_arr[::ds], I_arr[::ds], c=colors[0])
	inc_deg = np.around(inc * 180 / np.pi, 0)
	title = r'Inclination = $' + str(int(inc_deg)) + '^{\circ} $'
	if inc_deg < 10: title += ' '
	plt.title( r'Inclination = $' + str(int(inc_deg)) + '^{\circ} $' )
	plt.xlabel(r'$\lambda$ (nm)')
	plt.ylabel(r'Flux  (erg / s $\times$ Hz $\times$ ${\rm cm}^2$)')

	## Add an inset showing the temperature of the visible surface
	if temp: 
		height_inset = 0.2 # height of the inset, in figure heights
		width, height = fig.get_size_inches() # size of the current figure in inches
		# the width and height of the inset should be the same in inches
		width_inset = height_inset * width / height # width of the inset in figure widths
		# create the new axes
		ax2 = fig.add_axes([0.3, 0.6, width_inset, height_inset])
		# plot the temperature inset
		divider = make_axes_locatable(ax2)
		cax = divider.new_vertical(size="5%", pad=0.15, pack_start=True)
		fig.add_axes(cax)
		# plot the temperature inset
		star.plot_temp(ax2, inclination, height_inset * height, cax=cax, exp=0)

	fig.subplots_adjust(bottom=0.15)
	fig.savefig(outfile, dpi=200)
	plt.close(fig)
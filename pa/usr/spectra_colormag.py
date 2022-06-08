# Requires: a directory of spectra of a star at different inclinations,
#	two files with filter transmission functions
# Output: plots of the star's spectrum and position on a CMD for different inclinations
# Note: these can be assembled into a movie, for example with the following command:
# 	ffmpeg -framerate 10 -pattern_type glob -i '*.jpeg' -c:v libx264 -pix_fmt yuv420p -vf pad="width=ceil(iw/2)*2:height=ceil(ih/2)*2"

from pa.lib import star as st
from pa.lib import util as ut
from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np
import os

ddir = '../../' # location of the data directory
idir = ddir + 'data/vega/' # spectrum input directory
odir = ddir + 'data/vega_img/' # output directory
Vfile = ddir + 'data/Generic_Bessell.V.dat' # filter file
Bfile = ddir + 'data/Generic_Bessell.B.dat' # filter file
w = 3000 # don't plot above this wavelength
distance = 2.3694e19 # distance to the star in cm
V0 = 3.619e-9 # visual flux zero point
B0 = 6.317e-9 # blue flux zero point

Vwl, Vt = np.loadtxt(Vfile).T
Vfilt = ut.Filter(Vt, Vwl, V0)
Bwl, Bt = np.loadtxt(Bfile).T
Bfilt = ut.Filter(Bt, Bwl, B0)

# create the output directory if it doesn't exist
if not os.path.exists(odir):
	os.mkdir(odir)

## get spectrum input file names
# ignore hidden files
input_filenames = (f.name for f in os.scandir(idir) if not f.name.startswith('.'))
# add directory names
infiles = sorted( [idir + i for i in input_filenames] )

# create a star with the parameters found in the first input file used for the temperature inset
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

# get spectra
inclinations = []
incs = []
wl = []
flux = []
V = []
B = []
for infile in infiles:
	# create the output file path
	outfile = odir + os.path.splitext(os.path.basename(infile))[0] + '.jpeg'
	# get inclination
	f = open(infile)
	for line in f:
		if 'inclination' in line:
			incl = np.float( line.split()[-1] )
			inclinations.append(incl)
			incs.append(np.round( incl, decimals=2 ))
			break
	f.close()
	if len(wl) == 0:
		wl, fl = np.loadtxt(infile).T
		m = wl < w
		wl = wl[m]
	else:
		fl = np.loadtxt(infile)[:, 1]
	fl = fl[m] 
	flux.append(fl)
	# compute the magnitudes
	V.append( Vfilt.mag(fl, wl, distance) )
	B.append( Bfilt.mag(fl, wl, distance) )
V = np.array(V)
B = np.array(B)
flux = np.array(flux)
wl = np.array(wl)
inclinations = np.array(inclinations)
incs = np.array(incs)
color = B - V

## plot spectra
plt.rcParams.update({'font.size': 18})
rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('text', usetex=True)
dpi=200

for i in range(len(inclinations)):
	# output file
	outfile = odir + os.path.splitext(os.path.basename(infiles[i]))[0] + '.jpeg'

	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
	fig.suptitle('Inclination = %5.2f radians' % incs[i])

	# spectrum plot
	max_I = max(flux[i])
	min_I = min(flux[i])
	delta_y = max(flux[i]) - min(flux[i])
	offset_y = delta_y * 0.1

	ax1.set_ylim([min_I - offset_y, max_I + offset_y])
	ax1.set_ylim((-2.98405e+18, 3.282455e+19)) # temporarily fix the flux range
	ax1.scatter(wl, flux[i], marker='o', c='b', s=3)
	ax1.set_xlabel(r'$\lambda$ (nm)')
	ax1.set_ylabel(r'$\mathcal{F}_\nu \enspace (\mathrm{erg} \; \mathrm{s}^{-1} \; \mathrm{Hz}^{-1} \; \mathrm{cm}^{-2})$')
	
	# an inset showing the temperature of the visible surface
	height_inset = 0.2 # height of the inset, in figure heights
	width, height = fig.get_size_inches() # size of the current figure in inches
	# the width and height of the inset should be the same in inches
	width_inset = height_inset * width / height # width of the inset in figure widths
	# create the new axes
	ax1b = fig.add_axes([0.15, 0.6, width_inset, height_inset])
	# plot the temperature inset
	star.plot_temp(ax1b, inclinations[i], height_inset * height)
	
	# magnitudes plot
	max_V = np.max(V)
	min_V = np.min(V)
	delta_V = max_V - min_V
	offset_V = delta_V * 0.1
 
	ax2.set_ylim([max_V + offset_V, min_V - offset_V])
	ax2.scatter(color, V, marker='o', c='b', s=15)
	ax2.scatter(color[i], V[i], marker='o', edgecolors='darkorange', facecolors='none', s=60, linewidths=3)
	ax2.set_xlabel('B - V')
	ax2.set_xticks(0.01 * np.arange(7)[1::2])
	ax2.set_ylabel('V')

	# save figure
	fig.subplots_adjust(wspace=0.25)
	fig.savefig(outfile, dpi=dpi, bbox_inches='tight')
	plt.close(fig)
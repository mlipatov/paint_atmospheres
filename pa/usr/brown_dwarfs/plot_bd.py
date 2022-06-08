# collect all the spectra of rotating and non-rotating brown dwarfs and make figures for the paper.
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import os, glob, argparse
from pa.lib import util as ut

## for each inclination of the rotating dwarf, 
## find the non-rotating dwarf whose spectrum is closest
# given an independent spectrum x (e.g. a rotating dwarf's spectrum)
# and a dependent spectrum y (e.g. a non-rotating dwarf's spectrum)
# conduct linear regression: model y(x) as a straight line, minimize the rms deviation
# between y(x) and y, output the value 
# (Wikipedia:Simple_linear_regression, Wikipedia:Root-mean-square_deviation)
def rmsd_func(x, y):
	n = len(x)
	x_norm = x - x.mean() 
	y_norm = y - y.mean()
	beta = np.sum(x_norm * y_norm) / np.sum(x_norm**2)
	return np.sqrt( np.sum((y_norm - beta*x_norm)**2 / n) )

parser = argparse.ArgumentParser(description="Examples: \n" +\
		"python plot_bd.py ../../../data/bd_spectra/betapicb/ --min 950 --max 1125 --ratio; " +\
		"python plot_bd.py ../../../data/bd_spectra/J0348-6022_few/ -d 10 -i 45 --ratio")
parser.add_argument("directory", help="the directory with spectra and a directory called plots")
parser.add_argument("-i", help="inclination of the rotating dwarf in the spectrum ratio", type=int, default=0)
parser.add_argument("-d", help="downsample factor d, plot every dth wavelength; default is 1", type=int, default=1)
parser.add_argument("--min", help="minimum wavelength to plot in nm; default is 0", type=float, default=0)
parser.add_argument("--max", help="maximum wavelength to plot in nm; default is infinity", type=float, default=-1)
parser.add_argument("--ratio", default=False, action="store_true", help="Include ratios with non-rotating dwarf spectra")
parser.add_argument("--temploc", default='upper_right', help="Location of temperature color bar: e.g., lower_right")
args = parser.parse_args()

# wavelength limits
lammin = args.min
if args.max < 0:
	lammax = np.inf
	strmax = 'inf'
else:
	lammax = args.max
	strmax = str(int(lammax))
ds = args.d # downsample parameter
inc = args.i # reference inclination

# directory and filename
iodir = args.directory.strip("/")
name = os.path.basename(iodir)
if name == 'betapicb':
	fancy_name = r'$\beta$ Pictoris b'
elif 'J0348-6022' in name:
	fancy_name = r'J0348-6022, $\bar{T}_{\rm eff} = 880\,\mathrm{K}$'
else:
	fancy_name = name
os.makedirs(iodir + '/plots', exist_ok=True) # make the directory if it does not exist
outfile = iodir + '/plots/' + name + '_' + str(int(lammin)) + '-' + strmax +\
	'_ds' + str(ds) + '_i' + str(inc) + '.pdf' 

# plotting parameters
rat_yrange = 2.0 # y range for the ratio plot; optimized for the 880 K model
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif", 
    "font.serif": "Computer Modern",
    "font.size": 20
})
cols = plt.rcParams['axes.prop_cycle'].by_key()['color']
xlabel = r'$\lambda$ (nm)'
ylabel = r'${\cal F}$ (erg / s $\times$ Hz $\times$ ${\rm cm}^2$)'
mask = None
lam = None

# parameters of rotating dwarfs
incs_rot = []
labels_rot = []
fluxes_rot = []
colors_rot = []
lum_rot = None
omega_rot = None

# parameters of non-rotating dwarfs
lums_norot = []
labels_norot = []
fluxes_norot = []
radius_norot = None

# if working with J0348-6022, include only certain rotational speeds and inclinations
if 'J0348-6022' in name:
	omega_rot = 0.42
	incs_rot = [0, 45, 90]

# input data: spectra of rotating BDs at different inclinations
# and those of non-rotating objects at different luminosities
i = 0 # rotating count
for file in sorted(glob.glob(iodir + '/**/*.txt', recursive=True)):
	# get the parameters
	comments = open(file).readlines()
	luminosity = float(comments[0].split()[2])
	omega = np.around(float(comments[1].split()[2]), 3)
	inclination = np.around(float(comments[2].split()[2]) * 180 / np.pi, 3)
	radius = np.around(float(comments[4].split()[2]), 4)
	if not ('J0348-6022' in name) or \
			( omega in [0.0, omega_rot] and inclination in incs_rot ):
		# load intensity versus wavelength
		data = np.loadtxt(file)
		# set the wavelengths and the mask
		if mask is None: mask = (data[:, 0] <= lammax) & (data[:, 0] >= lammin)
		if lam is None: lam = data[:, 0][mask][::ds]
		# flux
		flux = data[:,1][mask][::ds]
		if omega > 0: # rotating
			if omega_rot is None: omega_rot = omega
			if lum_rot is None: lum_rot = luminosity
			lb = r'$\omega_{\rm R}$ = ' + str(omega)
			# if the inclination is the one to compare, record the index
			if inclination == inc: i_inc = len(fluxes_rot)
			# label includes inclination
			lb = lb + ', $i$ = ' + str(int(inclination)) + '$^{\circ} $'
			labels_rot.append(lb)
			colors_rot.append(cols[i])
			fluxes_rot.append(flux)
			incs_rot.append(inclination)
			i += 1
		else: # non-rotating
			if radius_norot is None: radius_norot = radius
			lb = r'$\omega$ = 0'
			# label includes luminosity
			lb = lb + ', $L = L_' + str(len(lums_norot)) + '$'
			labels_norot.append(lb)
			fluxes_norot.append(flux)
			lums_norot.append(luminosity)
# sort the non-rotating lists and remove extras if working with the T7 dwarf
if ('J0348-6022' in name):
	lums_norot = np.array(lums_norot)
	inds = lums_norot.argsort()
	lums_norot = lums_norot[inds] # [::4]
	fluxes_norot = np.array(fluxes_norot)[inds] # [::4]
	labels_norot = np.array(labels_norot)[inds] # [::4]

# bolometric and frequency-averaged fluxes
freq = 1/lam
freq_diff = np.diff(freq)
weights = np.append(freq_diff, 0) + np.insert(freq_diff, 0, 0)
# a dictionary indexed by rotational speed and inclination, just like fluxes
bolo_rot = [np.sum(f*weights) for f in fluxes_rot]
# frequency-averaged flux
avg_rot = bolo_rot / (freq[-1] - freq[0])

# index of the non-rotating model that fits the rotating model at a given inclination best
min_ind = []
for i in range(len(fluxes_rot)): # for each inclination of the rotating dwarf
	frot = fluxes_rot[i]
	rmsd = [rmsd_func(frot, fnrt) for fnrt in fluxes_norot]
	minr = min(rmsd)
	mini = rmsd.index(minr)
	min_ind.append(mini)

# effective temperature in K
T = ( lums_norot / (4 * np.pi * radius_norot**2 * ut.sigma_sun) )**(1/4)
Tmin = T[min_ind[i_inc]] # temperature of the spectrum with smallest deviation

## plot the fluxes of rotating models
# figure and axes
if args.ratio:
	fig, (ax0, ax1) = plt.subplots(2, 1, sharex=True, figsize=(10, 10))
else:
	fig, (ax0) = plt.subplots(1, 1, sharex=True, figsize=(10, 5))
# set the labels for the regular spectrum plots
ax0.set_ylabel(ylabel)
# plot the rotators
for i in range(len(fluxes_rot)):
	ax0.plot(lam, fluxes_rot[i], label=labels_rot[i], c=colors_rot[i])
# for full-range plot, draw a horizontal line at the average flux of the focal inclination
if lammin == 0 and lammax == np.inf:
	ax0.hlines(y=avg_rot[i_inc], xmin=lam[0], xmax=lam[-1], colors='grey',\
		linestyles='--', label=r'$\bar{\cal F}_{i=' + str(inc) + '^\circ}$')
ax0.legend(title=fancy_name, fontsize=16, title_fontsize=16)
# plot the flux ratios
if args.ratio:
	cmapBig = mpl.cm.get_cmap('cividis', 512)
	cmap = mpl.colors.ListedColormap(cmapBig(np.linspace(0, 1, 256)))
	norm = mpl.colors.Normalize(vmin=T[0], vmax=T[-1])
	## plot the ratios of fluxes
	ax1.set_ylabel(r'${\cal F}_{\omega_{\rm R} = ' + str(omega_rot) + ',\,i=' + str(inc) + \
		'^{\circ}} / {\cal F}_{\omega = 0}$') 
	ax1.set_xlabel(xlabel)
	# ratio of one inclination to non-rotators
	for i in range(len(fluxes_norot)):
		if i == min_ind[i_inc]: 
			col = 'crimson' # color for the best-matching non-rotating spectrum
		else: 
			col = cmap((T[i] - T[0])/(T[-1] - T[0]))
		lb = r'${\cal F}_{i=' + str(inc) + '^{\circ}} / {\cal F}_{\omega = 0,\, L = L_' + str(i) +'}$'
		ax1.plot(lam, fluxes_rot[i_inc]/fluxes_norot[i], '-', label=lb, c=col )
	# adjust the top y-limit so that y range is the same in all ratio plots
	ylim = ax1.get_ylim()
	print('y limits and range: ' + str(ylim) + ", " + str(ylim[-1] - ylim[0]))
	ymin = ax1.get_ylim()[0]
	ax1.set_ylim(ymax = ymin + rat_yrange)
	## color bar to indicate the effective temperature of the non-rotating dwarfs
	# the colorbar axes
	if args.temploc == 'upper_right': cax = fig.add_axes([0.65, 0.4, 0.2, 0.015]) 
	elif args.temploc == 'lower_right':	cax = fig.add_axes([0.65, 0.2, 0.2, 0.015])
	elif args.temploc == 'upper_center': cax = fig.add_axes([0.4, 0.4, 0.2, 0.015])	
	elif args.temploc == 'lower_center': cax = fig.add_axes([0.4, 0.2, 0.2, 0.015])
	ticks = ticker.LinearLocator(2)
	cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal', ticks=ticks, format='%.0f')
	cb.set_label(r'$T_{{\rm eff},\,\omega = 0}$ (K)')
	cb.ax.plot([Tmin]*2, [-0.5, 1.5], 'crimson', lw=5) # my data is between 0 and 1
# save figure file
fig.savefig(outfile, bbox_inches='tight')
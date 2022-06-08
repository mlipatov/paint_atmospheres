# make plots that analyze the spectra of all models - rotating at various speeds and inclinations
# and non-rotating - of the T7 dwarf
import numpy as np
import glob, os, argparse, pickle
# from sigfig import round
from matplotlib import pyplot as plt

parser = argparse.ArgumentParser(description="Examples: \n" +\
		"python plot_T7.py 400" )
parser.add_argument("temp", help="average temperature", type=int)
args = parser.parse_args()
tempstr = '_' + str(args.temp) +'/'

# plotting parameters
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif", 
    "font.serif": "Computer Modern",
    "font.size": 20
})
cols = plt.rcParams['axes.prop_cycle'].by_key()['color']

# given an independent spectrum x (e.g. a rotating dwarf's spectrum)
# and a dependent spectrum y (e.g. a non-rotating dwarf's spectrum)
# conduct linear regression: model y(x) as a straight line, minimize the rms deviation
# between y(x) and y, output the value 
# (Wikipedia:Simple_linear_regression, Wikipedia:Root-mean-square_deviation)
# Inputs: 
# 	independent spectrum x
# 	dependent spectrum y
#	weights of x and y due to differences between neighboring wavelengths
def rmsd_func(x, y, x_mean):
	n = len(x)
	x_norm = x - x.mean() # subtract the mean of the independent spectrum
	y_norm = y - y.mean() # same with the dependent spectrum
	# slope of the best-fitting linear model for the dependent spectrum
	beta = np.sum(x_norm * y_norm) / np.sum(x_norm**2)  
	# root mean squared deviation of the model from the dependent spectrum
	rmsd = np.sqrt( np.sum((y_norm - beta*x_norm)**2 / n) )
	# normalize the minimum RMSD by the average flux of the independent, observed spectrum
	rmsd /= x_mean
	return rmsd

# convert Roche-model dimensionless velocities to actual dimensionless velocities 
# we expect for the rotating dwarfs. To do this, calculate the dimensionless moment of inertia
# from the Darwin-Radau relation, the observed oblateness and the observed velocity at this oblateness.
# Then assume that the moment of inertia remains the same at lower rotational velocities.
# Inputs:
#	list of Roche velocities
#	observed oblateness
#	observed rotational velocity
# Output: list of real dimensionless velocities
def oms(omsRoche, f_obs=0.08, om_obs=0.36050):
	# oblatenesses corresponding to the given Roche velocities
	fR = 1/(1 + 2/omsRoche**2) 
	# moment of inertia from observed oblateness and velocity
	C = (2/3)*(1 - (2/5)*((5/2)*(om_obs**2/f_obs)-1)**(1/2)) 
	# velocities from oblatenesses and moment of inertia
	om = np.sqrt(fR*(116 - 300*C + 225*C**2)/40)
	return om 

def cos2deg(cos): # works only for angles between 0 and 90
	cos[cos > 1] = 1
	cos[cos < 0] = 0
	return np.arccos(cos) * (180 / np.pi)
def deg2cos(deg):
	return np.cos(deg * np.pi / 180)

lum_rot = None # luminosity of rotating models
fluxes_rot = dict() # rotating model fluxes, indexed by rotational speed
incs_rot = [] # inclinations of rotating models, should be the same for each rotational speed
lam = None # wavelengths
radius_norot = None # radius of non-rotating models
fluxes_norot = [] # non-rotating model fluxes
lums_norot = [] # non-rotating model luminosities

# input data: spectra of rotating BDs at different inclinations and rotational speeds
# and those of non-rotating objects at different luminosities
ratlumdir = '../../../data/bd_spectra/J0348-6022_ratlum/'
os.makedirs(ratlumdir, exist_ok=True)
iodir = '../../../data/bd_spectra/J0348-6022' + tempstr
for file in sorted(glob.glob(iodir + '**/*.txt', recursive=True)):
	# get the parameters
	comments = open(file).readlines()
	luminosity = float(comments[0].split()[2])
	omega = np.around(float(comments[1].split()[2]), 3)
	inclination = np.around(float(comments[2].split()[2]) * 180 / np.pi, 3)
	radius = np.around(float(comments[4].split()[2]), 4)
	dist = float(comments[5].split()[2])
	# load intensity versus wavelength
	data = np.loadtxt(file)
	# set the wavelengths and the flux
	if lam is None: lam = data[:, 0]
	flux = data[:,1]
	if omega > 0: # rotating
		if lum_rot is None: lum_rot = luminosity
		if omega not in fluxes_rot: fluxes_rot[omega] = []
		fluxes_rot[omega].append(flux)
		# inclinations should all be the same, so record only for one rotational speed
		if omega == 0.5:
			incs_rot.append(inclination)
	else: # non-rotating
		if radius_norot is None: radius_norot = radius
		fluxes_norot.append(flux)
		lums_norot.append(luminosity)
# cosine of inclination, rounded
cosi = np.cos(np.array(incs_rot) * np.pi / 180)
cosi = np.around(cosi, 10)
# rotational speeds
omegas = list(fluxes_rot.keys())

# bolometric fluxes
freq = 1/lam
freq_diff = np.diff(freq)
weights = np.append(freq_diff, 0) + np.insert(freq_diff, 0, 0)
bolo_nrt = [np.sum(f*weights) for f in fluxes_norot]
# a dictionary indexed by rotational speed and inclination, just like fluxes
bolo_rot = {omega: [np.sum(f*weights) for f in fluxes_rot[omega]] for omega in fluxes_rot}
# frequency-averaged flux
avg_rot = {omega: [f / (freq[-1] - freq[0]) for f in bolo_rot[omega]] for omega in bolo_rot}

# minimum rms deviation of the linear model of a non-rotating spectrum from the rotating dwarf's spectrum,
# indexed by rotational speed and inclination
min_rmsd = dict()
# factor in the ratio of the best-fit model luminosity to the actual luminosity of the dwarf
# that only deviates from 1 because of the spatial asymmetry of the observed flux due to rotational effects,
# as opposed to the difference between the dwarf's flux and the flux of the best-matching model; here,
# all flux is bolometric. Indexed by rotational speed and inclination.
rat_lum = dict()
for omega in omegas: # for each rotational speed
	min_rmsd[omega] = []
	rat_lum[omega] = []
	for i in range(len(incs_rot)): # for each inclination of the rotating dwarf
		frot = fluxes_rot[omega][i]
		arot = avg_rot[omega][i]
		rmsd = [rmsd_func(frot, fnrt, arot) for fnrt in fluxes_norot]
		minr = min(rmsd)		
		min_rmsd[omega].append(minr)
		mini = rmsd.index(minr)
		# ratio of the luminosity of the non-rotating model to that of the rotating dwarf, 
		# corrected by the corresponding inverse ratio of bolometric fluxes
		rat_lum[omega].append( 1 / ((lums_norot[mini] / lum_rot) * (bolo_rot[omega][i] / bolo_nrt[mini])) )

# save the anisotropy ratio 
file = ratlumdir + 'rat_lum_J0348-6022' + tempstr.replace('/','') + '.pkl'
with open(file, 'wb') as f:
    pickle.dump(rat_lum, f)

# coefficients for quadratic fits (highest power first) 
# to the dependence of luminosity ratio correction on cosine of inclination
pr = []
for omega in omegas: 
	pr.append(np.polyfit(cosi, rat_lum[omega], 3))
pr = np.array(pr)
# coefficients for meta-fits - quadratic fits of the coefficients as functions of omega
ppr = []
for i in range(pr.shape[1]): # for each power of the fits, starting from highest
	p = np.polyfit(omegas, pr[:, i], 2)
	ppr.append(p)
ppr = np.array(ppr)
ppr = np.around(ppr, 4) # round to 3 digits after the decimal
# do the calculation backwards from the meta-coefficients of the fits against omega
# to get coefficients of the fits of luminosity ratio versus cosine of inclination for each omega
pr2 = np.zeros_like(pr)
for i in range(ppr.shape[0]): # for each power of the original fits, starting from highest
	pr2[:, i] = np.polyval(ppr[i, :], omegas) # compute the coefficients for all omegas from the meta-fits

## plots
xlabel = r'$\cos{i}$'
plotdir = iodir + 'plots/'
os.makedirs(plotdir, exist_ok=True)
legend_title = r'$\bar{T}_{\rm eff} = ' + tempstr.replace('_','').replace('/','') + '\,\mathrm{K}$'
legend_fontsize = 14
incs = [0, 30, 60, 90]

# plot the rms deviation versus inclination for different rotational velocities
if '880' in tempstr: om = omegas[::2]
else: om = omegas
for omega in om:
	plt.plot(cosi, min_rmsd[omega], label=r'$\omega = ' + str(np.around(oms(omega), 2)) + '$')
plt.legend(title=legend_title, fontsize=legend_fontsize, title_fontsize=legend_fontsize)
plt.xlabel(xlabel)
plt.ylabel(r'$\min{{\rm RMSD}({\hat{\cal F}}_{\omega = 0})/\bar{{\cal F}}}$')
plt.gca().invert_xaxis()
secax = plt.gca().secondary_xaxis('top', functions=(cos2deg, deg2cos))
secax.set_xlabel(r'$i$ (degrees)', labelpad=10)
secax.set_xticks(incs)
plt.tight_layout()
plt.savefig(plotdir + 'min_rmsd_J0348-6022' + tempstr.replace('/','') + '.pdf', dpi=200)
plt.clf()

# plot the rotational asymmetry factor in luminosity determination,
# the discrete points - for each object temperature
xarr = np.linspace(cosi[0], cosi[-1], 100)
glb = glob.glob(ratlumdir + '*.pkl', recursive=True) # glob of anisotropy ratios at all BD temperatures
# check if this is the main temperature
if '880' in tempstr: start = 0; skip = 2
else: start = 0; skip = 1
if len(glb) > 0:
	rat_diff = [] # all anisotropy ratios
	for file in glb:
		# load the anisotropy ratio 
		with open(file, 'rb') as f: rat_lum = pickle.load(f)
		j = 0
		rd = []	# small ratio difference
		for i in range(start, len(omegas), skip):
			omega = omegas[i]
			# evaluate the fit at the locations of discrete data points
			rl_fit = np.polyval(pr2[i, :], cosi)
			# plot
			plt.scatter(cosi, rat_lum[omega], c=cols[j], alpha=0.5, lw=0) #, s=3, marker='.')
			# store the anisotropy ratio difference between the fit and the data points
			rd.append(rl_fit - rat_lum[omega])
			# append the color counter
			j+=1
		rat_diff.append(rd)
	# dimensions: average temperature, omega, inclination
	rat_diff = np.array(rat_diff)
	# index of the maximum absolute difference: this is probably the largest inclination, omega and temperature 
	rat_ind = np.unravel_index(np.argmax(np.abs(rat_diff)), rat_diff.shape)
	print('This should be (4, 3, 10) - indices of largest temperature, omega and inclination: ' + str(rat_ind))
	print('The largest fit minus data points difference at these parameter values is ' + str(rat_diff[rat_ind]))
	# rmsd for the shown data points
	rd_flat = rat_diff.flatten()
	print('The RMSD for the shown data points is ' + str(np.sum(rd_flat**2 / rd_flat.shape)))
	# plot the fit lines
	j = 0
	for i in range(start, len(omegas), skip):
		label = r'$\omega = ' + str(np.around(oms(omegas[i]), 2)) +	'$, $\omega_\mathrm{R} = ' + str(omegas[i]) + '$'	
		plt.plot(xarr, np.polyval(pr2[i, :], xarr), label=label, c=cols[j])
		j+=1
	# plot details
	plt.hlines(y=1, xmin=cosi[0], xmax=cosi[-1], color='grey', linestyle='--', zorder=-10)
	plt.legend(title=legend_title, fontsize=legend_fontsize, title_fontsize=legend_fontsize)
	plt.xlabel(xlabel)
	plt.ylabel(r'$L\,p / (4\pi r^2{\cal F})$')
	plt.gca().invert_xaxis()
	secax = plt.gca().secondary_xaxis('top', functions=(cos2deg, deg2cos))
	secax.set_xlabel(r'$i$ (degrees)', labelpad=10)
	secax.set_xticks(incs)
	plt.tight_layout()
	plt.savefig(plotdir + 'rat_lum_J0348-6022' + tempstr.replace('/','') + '.pdf', dpi=200)

print('coefficients of 2nd degree polynomials in omega, ' +
	'which give coefficients of 3rd degree polynomials in inclination, ' + 
	'highest power first.')
np.set_printoptions(suppress=True)
print(ppr)
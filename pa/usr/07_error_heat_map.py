# Runtime on 2.3GHz MacBook Pro with 8 Gb of RAM: ~16 minutes
# Requires: a directory of (~1000) spectrum files that use nz = 10,000 at a set of equally spaced inclinations,
#   a directory of spectrum files that use nz = 100, corresponding to the same inclinations;
#   np.sort(filenames) in each directory should sort them in inclination, from zero to pi/2.
# Notes: to get the required files, run calc_star and calc_spectra.

import matplotlib.pyplot as plt
from matplotlib import gridspec, rc, cm
from astropy.io import fits
import numpy as np
from scipy import interpolate
import glob
import os

nz = '100'
om = '999'

vmax = -2
vmin = -5

rc('font',**{'family':'serif','serif':['Computer Modern'],'size': 18})
rc('text', usetex=True)

iodir = '../../' # location of the input/output directory

# These should give sorted file lists.  Note that you have to choose
# the names carefully to make sure that they are sorted by np.sort.
# I used a rename shell script, but there are other ways.  It might 
# be worth modifying calc_spectra.py to print zeros, e.g. 1.100 rather
# than 1.1.

filelist1 = list(np.sort(glob.glob(iodir + 'data/vega' + str(om) + '_100/vega*.txt'))) # nz = 100
# for x in filelist1: print(os.path.basename(x)) 
filelist2 = list(np.sort(glob.glob(iodir + 'data/vega' + str(om) + '_10K/vega*.txt'))) # nz = 10,000

if len(filelist1) != len(filelist2):
    print("Mismatched lists of files!")
    exit()
n_inc = len(filelist1)

# Compute the normalized difference.  The 1e-300 prevents division by zero.

spec_max = 0

diffarr = np.zeros((n_inc, np.loadtxt(filelist1[0]).shape[0]))
for i in range(n_inc):
    spec1 = np.loadtxt(filelist1[i])[:, 1]
    spec2 = np.loadtxt(filelist2[i])[:, 1]
    diffarr[i] = np.abs(spec1 - spec2)/(spec2 + 1e-300)
    spec_max = np.max( ( (spec2[:131]/spec2.max()).max(), spec_max) )

print('Below 100 nm, spectra relative to maximum do not rise above ' + str(spec_max))
print('Above 100 nm, greatest relative error is ' + str(diffarr[:,131:].max()))

# Read in spectra to plot on the top of the axes.  We'll read in
# the lowest, middle, and highest inclinations.  We will normalize
# each one by the spectrum of the pole-on star.

lam, spec_mid = np.loadtxt(filelist2[n_inc//2]).T
spec_pole = np.loadtxt(filelist2[0]).T[1]
spec_eq = np.loadtxt(filelist2[-1]).T[1]
norm = np.amax(spec_pole)
spec_mid /= norm
spec_pole /= norm
spec_eq /= norm


# Two stacked plots
fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=False,
                                           gridspec_kw = {'height_ratios':[0.34, 1]},
                                           figsize=(6*0.9, 6))

# The reference spectrum 
ax1.plot(np.linspace(0, 1, len(spec_pole)), spec_pole, color='grey')
ax1.set_yticks([])

# The heat map plot.  1e-300 prevents log of zero.
img = ax2.imshow(np.log10(diffarr + 1e-300), origin='lower', cmap='hot_r',
                 vmin=vmin, vmax=vmax, extent=[0, 1, 0, 1], aspect=0.81)

if om == '615':
    ax1.text(0.7, 0.7, r'$\omega = 0.6151$', fontsize=18)
elif om == '999':
    ax1.text(0.7, 0.7, r'$\omega = 0.9990$', fontsize=18)

# Label the wavelength axis
lam_label = [10, 100, 250, 500, 1000, 3000, 10000]
f = interpolate.interp1d(lam, np.linspace(0, 1, len(lam)), kind='cubic')
ax2.set_xticks(f(np.asarray(lam_label)))
ax2.set_xticklabels(['%d' % x for x in lam_label])
ax2.set_xlabel(r'$\lambda$ (nm)')

# Label the inclination axis
inc_label = [0, 30, 60, 90]
ax2.set_yticks(np.asarray(inc_label)/90.)
ax2.set_yticklabels(['%d' % x for x in inc_label])
ax2.set_ylabel('Inclination (degrees)')

# Remove excess spacing and save
plt.subplots_adjust(hspace=.0)
plt.subplots_adjust(wspace=.0)
plt.savefig(iodir + 'accuracy_nz' + nz + '_om0p' + om + '.pdf')



# The color bar figure
fig, (a, ax) = plt.subplots(2, 1, gridspec_kw = {'height_ratios':[0.34, 1]}, 
                                  figsize=(1, 6))
fig.subplots_adjust(left=0, right=0.2)
a.axis('off')

ticks = np.arange(vmin, vmax + 1e-10, 1)
labs = ['$%d$' % (tick) for tick in ticks]

# Color map
cmap = cm.hot_r
norm = plt.Normalize(vmin=vmin, vmax=vmax)

# Color bar
cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ticks=ticks,
             cax=ax, orientation='vertical')
cb.ax.set_xticklabels(labs)
cb.ax.set_ylabel(r'$\log_{10}{|\delta \mathcal{F}_\nu \,/\, \mathcal{F}_\nu|}$',
                  fontsize='large', labelpad=6)

# Save figure
fig.savefig(iodir + 'accuracy_colorbar.pdf')
# Runtime on 2.3GHz MacBook Pro with 8 Gb of RAM: ~1 minute
# Requires: a limb darkening information file (with fits of intensity w.r.t surface inclination)
# Output: convergence of the longitudinal integral at a specific wavelength

import math
import numpy as np
import pickle
import pa.lib.star as star
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
from matplotlib import rc

iodir = '../../' # location of the input/output directory

## unpickle the limb darkening information
with open(iodir + 'data/limbdark_m01.pkl', 'rb') as f:
	ld = pickle.load(f)

lam = ld.lam

# choose a wavelength that is close to the star's peak at all inclinations
wl = 511. 
ind = np.where(lam == wl)[0][0]

# number of z values: log2 scale
low, high = [3, 13]
m = high - low + 1 
n_z = np.logspace(low, high, num=m, base=2.)

diff_trap = []
diff_cubic = []
omega, luminosity, mass, Req, distance, inclination = [0.6151, 40.346, 2.165, 2.815, 2.3694e19, math.pi/4]

print("calculating the etalon spectrum with N = 10,000...")
st = star.Star(omega, luminosity, mass, Req, distance, 1e4, ld=ld)
ref = st.integrate(inclination, method='cubic')
mask = ref > 0

print("calculating the flux at the indicator wavelength")
for i, n in np.ndenumerate(n_z):
	st = star.Star(omega, luminosity, mass, Req, distance, n, ld=ld)

	light = st.integrate(inclination, method='cubic')
	sd = np.abs(1 - light[ind] / ref[ind])
	diff_cubic.append( sd )

	light = st.integrate(inclination, method='trapezoid')
	sd = np.abs(1 - light[ind] / ref[ind])
	diff_trap.append( sd )

plt.rcParams.update({'font.size': 18})
rc('font',**{'family':'serif','serif':['Computer Modern']})
rc('text', usetex=True)

# convergence plot figure
fig = plt.figure()
# axes
ax = plt.axes()
ax.scatter(n_z, diff_trap, marker='o', c='b', s=40)
ax.scatter(n_z, diff_cubic, marker='o', c='g', s=40)
ax.loglog(n_z, diff_trap, c='b', linewidth=1)
ax.loglog(n_z, diff_cubic, c='g', linewidth=1)
ax.set_xlabel(r'${\cal N}$')
ax.set_ylabel(r'$|\delta \mathcal{F}_\nu \,/\, \mathcal{F}_\nu|$', labelpad=5)
diff = np.concatenate((diff_trap, diff_cubic))
max_diff = np.max(diff)
min_diff = np.min(diff)
max_n = max(n_z)
min_n = min(n_z)
ax.set_ylim(min_diff * 0.5, max_diff * 2)
ax.set_yticks([1e-9, 1e-7, 1e-5, 1e-3])
ax.set_xlim(min_n * 0.9, max_n * 1.1)
fig.savefig(iodir + 'vega6151_45deg511nm.pdf', dpi=200, bbox_inches='tight')
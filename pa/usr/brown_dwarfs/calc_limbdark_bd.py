# Given a set of brown dwarf atmosphere files, 
# interpolate to find intensities at wavelengths where they are missing,
# interpolate to find same at the temperature and gravity where they are missing,
# create a limbdarkening object that PARS can use.

import os
import numpy as np
from scipy.interpolate import griddata
from scipy.interpolate import interpn
import pickle
# PARS imports
from pa.lib import limbdark as ld
from pa.lib import util as ut
from pa.lib import fit as ft

di = './bd_atmospheres/'
files = os.listdir(di)
files.sort()
# viewing angles, count 16
mu = np.array([ float(m) for m in open(di + files[0], 'r').readline().split(',')[1:] ])
# wavenumbers in cm^-1, count 9831
data = np.loadtxt(di + files[0], skiprows=1, delimiter=',')
wno = data[:, 0]
# temperatures and gravities
lis = []
for file in files:
	Tstr, gstr = file.split('_')[0:2]
	T1 = int(Tstr[1:-1]); g1 = int(gstr[1:-3])
	lis.append([T1, g1])
lis.sort()
ar = np.array(lis)
Tlist = ar[:, 0] # temperatures in filenames
glist = ar[:, 1] # gravities in filenames
T = np.unique(Tlist) # unique temperatures
g = np.unique(glist) # unique gravities
grid = [glist[Tlist == t] for t in T] # gravities at each temperature

# intensity array should have dimensions (wavelengths, viewing angles, gravities, temperatures)
# missing entries should be numpy.nan
i = 0
I = np.full((len(wno), len(mu), len(g), len(T)), np.nan, dtype=np.float64)
for file in files:
	Tstr, gstr = file.split('_')[0:2]
	T1 = int(Tstr[1:-1]); g1 = int(gstr[1:-3])
	Ti = np.where(T == T1)[0][0]; gi = np.where(g == g1)[0][0]
	data = np.genfromtxt(di + file, skip_header=1, delimiter=',', filling_values=np.nan)[:, 1:]
	miss = np.where(np.isnan(data)) # locations where there is no data
	# interpolate to fill at the wavenumbers where intensity is missing
	# assume all mu at a wno are either present or not
	if len(miss[0]) > 0:
		wmiss = np.unique(miss[0]) # missing wavenumber indices
		points = ( np.delete(wno, wmiss), mu )
		values = np.delete(data, wmiss, axis=0)
		xi = np.array([wno[miss[0]], mu[miss[1]]]).T
		data[miss] = interpn(points, values, xi, method='linear')
	I[:, :, gi, Ti] = data # intensities on the grid of wavenumbers and viewing angles
	i+=1
	if i % 10 == 0: print(i)
# interpolate to get intensities at the temperature and gravity where they are missing
# this may be fastest to do "by hand"
mi = np.where(np.isnan(I)) # locations in the intensity array where data is missing
gi = np.unique(mi[2])[0]; Ti = np.unique(mi[3])[0] # indices 
# See Wikipedia: Bilinear Interpolation
x1 = g[gi - 1]; x2 = g[gi + 1]; x = g[gi]
y1 = T[Ti - 1]; y2 = T[Ti + 1]; y = T[Ti]
f11 = I[:, :, gi - 1, Ti - 1]
f12 = I[:, :, gi - 1, Ti + 1]
f21 = I[:, :, gi + 1, Ti - 1]
f22 = I[:, :, gi + 1, Ti + 1]
Iint = float((x2 - x1)*(y2 - y1))**-1 * \
	( f11*(x2 - x)*(y2 - y) + f21*(x - x1)*(y2 - y) + f12*(x2 - x)*(y - y1) + f22*(x - x1)*(y - y1) )
# insert the interpolated values
I[:, :, gi, Ti] = Iint
# convert wavenumber to wavelength in nanometers
lam = 1e7/wno
# flip arrays so that wavelength is increasing
wno = np.flip(wno); lam = np.flip(lam); I = np.flip(I, axis=0)
# convert gravity to cm / s^2 and take the log
logg = np.log10(g * 100)
# metallicity
Z = 0.0
# whether to save discrete intensities
save = False
# boundaries between viewing angle intervals
bounds = [0.1, 0.4]

# PARS calls
l = ld.LimbDark(lam, T, logg, Z) # initialize a limb darkening object
ft.mu_arr = mu # set the viewing angles
l.fit(I, bounds, save) # fit polynomials for I(mu)
### Pickle the limb darkening information
with open('./limbdark_BD_00.pkl', 'wb') as f:
	pickle.dump(l, f)
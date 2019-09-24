### explores the convergence properties of the synthetic star spectra
import math
import numpy as np
import pickle
import star.star as star
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def expon(x, a, b, c):
    return a*np.exp(b*x)+c

# get the information from the precise spectrum (z step = 5e-5)
# obtain the wavelengths and the intensities from the file
f = open('vega_40K_zvalues.txt')
wl_arr = []
I_arr = []
count = 0;
for line in f:
    if (not '#' in line) and line.strip():
        data = line.split()
        wl_arr.append(float(data[0]))
        I_arr.append(float(data[1]))
f.close()
I_arr = np.array(I_arr)
wl_arr = np.array(wl_arr)
mask = (I_arr != 0)

# z steps of the approaching spectra
z_step = np.arange(2e-3, 2e-4, -1e-4)
diff = []
omega, luminosity, mass, Req, inclination = [0.8760, 40.124, 2.135, 2.818, 0.08683]

## unpickle the limb darkening information
with open('limbdark_m01.pkl', 'rb') as f:
	ld = pickle.load(f)

for i, zs in np.ndenumerate(z_step):
	st = star.Star(omega, luminosity, mass, Req, zs, ld)
	light = st.integrate(inclination)
	d = np.abs(1 - light[mask] / I_arr[mask])
	sd = np.sum(d) / len(d)
	print ( zs, sd )
	diff.append( sd )

# fit an exponential
popt, pcov = curve_fit(expon, z_step, diff, p0=(1, 1e-6, 1))
zst = np.linspace(0, 2e-4, 1000)
diff_fit = expon(zst, *popt)
print(popt)

# plot
fig = plt.figure()
plt.scatter(z_step, diff, marker='o', c='g', s=6)
plt.title('difference with the previous step value')
plt.xlabel('z step')
plt.ylabel('difference')
plt.ylim(-0.0001, max(diff)+0.0001)
plt.xlim(0, 0.0025)
plt.plot(zst, diff_fit)
plt.show()
fig.savefig('plot.png')
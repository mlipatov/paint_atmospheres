# Plot the dependence of the proportion of Keplerian limit versus the MESA dimensionless velocity

from pa.lib import surface as sf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc

iodir = '../../'

plt.rcParams.update({
    "text.usetex": True,
    "font.serif": ["Computer Modern"],
    "font.size": 20
})

sf.calcVA()
sf.calcom()
oM = np.linspace(sf.omin, sf.omax, 1000)

fig = plt.figure()
ax = plt.axes()
ax.plot(oM, sf.omega(oM), lw=3)

ax.set_xlabel(r'$\omega_{\rm M}\sqrt{1 - \frac{L}{L_{\rm Edd}}}$')
ax.set_ylabel(r'$\omega$')
fig.savefig(iodir + 'omega_conversion.pdf', dpi=200, bbox_inches='tight')

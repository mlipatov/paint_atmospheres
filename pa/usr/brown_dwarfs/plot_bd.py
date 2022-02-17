# plot brown dwarf spectrum
from matplotlib import pyplot as plt
import numpy as np
import glob

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif", 
    "font.serif": "Computer Modern",
    "font.size": 20,
    "figure.figsize": (10, 6)
})
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
colors[2] = 'grey'
xlabel = r'$\lambda$ (nm)'
ylabel = r'Flux  (erg / s $\times$ Hz $\times$ ${\rm cm}^2$)'

iodir = '../../../data/J0348-6022_3inc/'

ds = 10 # downsample
data = []
labels = []
mask = None
i = 0
# rotating BD at different inclinations
for file in glob.glob(iodir + '*.txt'):
	# get the rotation rate and inclination
	comments = open(file).readlines()
	omega = np.around(float(comments[1].split()[2]), 3)
	inclination = np.around(float(comments[2].split()[2]) * 180 / np.pi, 3)
	# intensity versus wavelength
	dt = np.loadtxt(file)
	if mask is None:
		mask = (dt[:, 0] < 1350) & (dt[:, 0] > 1150)
	lb = r'$\omega = ' + str(omega) + ', i = ' + str(int(inclination)) + '^{\circ} $'
	data.append(dt) # save for later
	labels.append(lb)
	plt.plot(dt[:,0][mask], dt[:,1][mask], label=lb, c=colors[i])
	i += 1
# non-rotating version
dt = np.loadtxt('../../../data/J0348-6022_omega0/J0348-6022_omega01_413720.txt')
lb = r'$\omega = 0, R_{\rm p} < R < R_{\rm e}$'
data.append(dt)
labels.append(lb)
plt.plot(dt[:,0][mask], dt[:,1][mask], c=colors[2], label=lb)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend()
plt.tight_layout()
plt.savefig(iodir + 'plots/J0348-6022_zoomin.pdf')
# plt.show() # show the plot with all wavelengths

# make a downsampled plot
plt.clf()
data = np.array(data)[:, ::ds, :]
for i in range(data.shape[0]):
	plt.plot(data[i,:,0], data[i,:,1], label=labels[i], c=colors[i])
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.legend()
plt.tight_layout()
plt.savefig(iodir + 'plots/J0348-6022_ds10.pdf')
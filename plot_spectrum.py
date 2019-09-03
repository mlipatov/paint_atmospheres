# Plot a spectrum
import util as ut
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pickle

parser = argparse.ArgumentParser(description="Examples: \n" +\
	"python plot_spectrum.py \'vega1.dat\' 3000,\
	 python plot_spectrum.py \'sun1.dat\' 4000")
parser.add_argument("txt_file", help="spectrum text file")
parser.add_argument("wl", help="upper cutoff for wavelength, in nanometers", type=int)
args = parser.parse_args()

txt_file = args.txt_file # spectrum file
wl = args.wl # cutoff wavelength

# obtain the wavelengths and the intensities from the file
f = open(txt_file)
wl_arr = []
I_arr = []
count = 0;
for line in f:
    if (not '#' in line) and line.strip():
        data = line.split()
        if wl < float(data[0]):
        	break
        wl_arr.append(float(data[0]))
        I_arr.append(float(data[1]))
f.close()
I_arr = np.array(I_arr)
wl_arr = np.array(wl_arr)

# plot
delta_y = max(I_arr) - min(I_arr)
offset_y = delta_y * 0.1
fig = plt.figure()
plt.axes().set_ylim([min(I_arr) - offset_y, max(I_arr) + offset_y])
plt.scatter(wl_arr, I_arr, marker='o', c='b', s=3)
plt.title('intensity vs wavelength from ' + txt_file)
plt.xlabel('wavelength (nm)')
plt.ylabel('intensity (erg/s/Hz/ster)')
fig.savefig('I_vs_wl.png')
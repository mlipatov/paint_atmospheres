# Plot a spectrum
import util as ut
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pickle

parser = argparse.ArgumentParser(description="Examples: \n" +\
	"python plot_spectrum.py \'vega.dat\' 4000 -d 2.37e19 -s A -w A,    \
	 python plot_spectrum.py \'sun.dat\' 4000 -d 1.496e13 -p W -a m**2 -s nm")
parser.add_argument("txt_file", help="input spectrum text file")
parser.add_argument("wl", help="upper cutoff for wavelength, in nanometers", type=int)
parser.add_argument("-d", help="distance to the star in centimeters \
	(if the intensity is per area of photoreceptor, default is per steradian)", 
	type=float)
parser.add_argument("-p", help="units of power (e.g. W; default is erg/s)")
parser.add_argument("-a", help="units of photoreceptor area (e.g. m**2; default is cm**2)")
parser.add_argument("-s", help="units of specific intensity denominator (e.g. A or nm; default is Hz)")
parser.add_argument("-w", help="units of wavelength (e.g. A; default is nm)")
args = parser.parse_args()

txt_file = args.txt_file # spectrum file
wl = args.wl # cutoff wavelength
distance = args.d # distance to the star in centimeters
power = args.p # units of power
area = args.a # units of photoreceptor area
specific = args.s # units of specific intensity denominator
wavelength = args.w # units of wavelength

# convert units of cutoff wavelength to the units of the input file
if wavelength is not None:
	if wavelength == 'A':
		wl = 1.e1 * wl

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


# convert units of wavelength to nm
if wavelength is not None:
	if wavelength == 'A':
		wl_arr = 1.e-1 * wl_arr
# convert units of power to erg/s
if power is not None:
	if power == 'W':
		I_arr = 1.e7 * I_arr
# convert units of photoreceptor area to cm**2
if area is not None:
	if area == 'm**2':
		I_arr = 1.e-4 * I_arr
# convert units of the specific intensity denominator to Hz
if specific is not None:
	if specific == 'A':
		I_arr = ut.convert_from_A(I_arr, wl_arr)
	elif specific == 'nm':
		I_arr = ut.convert_from_nm(I_arr, wl_arr)
# convert to intensity per steradian from intensity per unit photoreceptor area
if distance is not None:
	I_arr = ut.convert_to_ster(I_arr, distance)

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
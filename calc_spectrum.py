# Sets up the spectrum calculation for a star
import limbdark.limbdark as limbdark
import limbdark.fit as ft
import star.star as star
import util as ut
import numpy as np
import sys
import time
import argparse
import pickle

parser = argparse.ArgumentParser(description="Example: \n" +\
	"python calc_spectrum.py \'vega.txt\' \'limbdark.pkl\' 0.0005 " +\
	"0.8760 0.08683 40.124 2.135 2.818")
parser.add_argument("output", help="an output spectrum text file to create")
parser.add_argument("pkl_lfile", help="the limb darkening .pkl file to access")
parser.add_argument("z_step", help="integration step for the normalized z coordinate", type=float)
parser.add_argument("omega", help="rotation speed divided by its Keplerian value at the equator", type=float)
parser.add_argument("inclination", help="angle between line of sight and rotation axis", type=float)
parser.add_argument("luminosity", help="luminosity in solar luminosities", type=float)
parser.add_argument("mass", help="mass in solar masses", type=float)
parser.add_argument("Req", help="equatorial radius in solar radii", type=float)
args = parser.parse_args()

## inputs
txt_sfile = args.output # spectrum text file
pkl_lfile = args.pkl_lfile # limb darkening file
z_step = args.z_step
## star parameter inputs 
omega = args.omega # dimensionless rotation speed
inclination = args.inclination # the angle between the axis of rotation and the line of sight
luminosity = args.luminosity # luminosity of the star in solar luminosities
mass = args.mass # mass of the star in solar masses
Req = args.Req # equatorial radius of the star in solar radii

## unpickle the limb darkening information
with open(pkl_lfile, 'rb') as f:
	ld = pickle.load(f)
# set the bounds between mu intervals in the Fit class
# with the bounds in the limb darkening information
ft.Fit.set_muB(ld.bounds)

# print a message
ut.printf ("Mass: " + str(mass) + " sun(s)\nLuminosity: "+ str(luminosity) + " sun(s)\n" +\
	"Equatorial radius: " + str(Req) + " solar radi(i/us)\nRotation speed: " + str(omega) +\
	" of the Keplerian angular velocity at the equator\nInclination: " +\
	str(inclination) + " radians.\n")
## Integrate light from a star with given physical parameters, resolution of map, 
## and limb darkening information
st = star.Star(omega, inclination, luminosity, mass, Req, z_step, ld)

# wl = 393
# ind_wl = np.argwhere(ld.wl_arr == wl)[0][0]
# print(st.light_arr[ind_wl])
# print(st.map.params_arr[ind_wl, :, -5:])
# print(st.map.fitint[:, -5:])
# print(st.map.logg_arr)
# print(st.map.temp_arr)

# write the spectrum of the star in text format
# create this file if it doesn't exist, open it for writing
f = open(txt_sfile,'w+') 
# write the header
f.write('# omega: ' + str(omega) + '\n')
f.write('# inclination: ' + str(inclination) + '\n')
f.write('# luminosity: ' + str(luminosity) + '\n')
f.write('# mass: ' + str(mass) + '\n')
f.write('# Req: ' + str(Req) + '\n')
f.write('# wavelength(nm)\tintensity(ergs/s/Hz/ster)\n') 
f.close() # close the file
f = open(txt_sfile, 'a') # open the file for appending
# write the spectrum to the file
ind = 0
while (ind < len(ld.wl_arr)):
	f.write(str(ld.wl_arr[ind]) + '\t %.5E\n' % st.light_arr[ind])
	ind += 1
f.close()
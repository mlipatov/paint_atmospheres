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
	"python calc_spectrum.py \'star.pkl\' \'limbdark.pkl\' 0.0005 " +\
	"0.8760 0.08683 40.124 2.135 2.818")
parser.add_argument("pkl_sfile", help="a name for a .pkl spectrum file to create")
parser.add_argument("pkl_lfile", help="name of the limb darkening .pkl file to access")
parser.add_argument("z_step", help="step for normalized z coordinate", type=float)
parser.add_argument("omega", help="rotation speed divided by its Keplerian value at the equator", type=float)
parser.add_argument("inclination", help="angle between line of sight and rotation axis", type=float)
parser.add_argument("luminosity", help="luminosity in solar luminosities", type=float)
parser.add_argument("mass", help="mass in solar masses", type=float)
parser.add_argument("Req", help="equatorial radius in solar radii", type=float)
args = parser.parse_args()

## inputs
pkl_sfile = args.pkl_sfile # spectrum file
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

# for Vega at 401 nm, 
# 4.918108154053354e+19 ergs/s/nm/ster for zstep = 0.0005 (4.4 sec),
# 4.92012754246987e+19 for z = 0.0001 (45.6 sec, 0.04% more) 
wl = 401
ind_wl = np.argwhere(ld.wl_arr == wl)[0][0]
# print(st.light_arr[ind_wl])
# print(st.map.params_arr[ind_wl])
# print(st.map.fitint)
print(st.map.logg_arr)
# print(st.map.temp_arr)

# pickle the star
with open(pkl_sfile, 'wb') as f:
	pickle.dump(st, f)
# Sets up the spectrum calculation for a star
import limbdark.fit as ft
import star.star as star
import util as ut
import argparse
import pickle

parser = argparse.ArgumentParser(description="Example: \n" +\
	"python calc_star.py \'vega.pkl\' \'limbdark_m01.pkl\' 4000 " +\
	"0.8760 40.124 2.135 2.818")
parser.add_argument("output", help="an output file containing the pickled star")
parser.add_argument("ld_file", help="the limb darkening .pkl file to access")
parser.add_argument("n_z", help="number of values for the normalized z coordinate", type=int)
parser.add_argument("omega", help="rotation speed divided by its Keplerian value at the equator", type=float)
parser.add_argument("luminosity", help="luminosity in solar luminosities", type=float)
parser.add_argument("mass", help="mass in solar masses", type=float)
parser.add_argument("Req", help="equatorial radius in solar radii", type=float)
args = parser.parse_args()

## inputs
pkl_sfile = args.output # star pickle file
pkl_lfile = args.ld_file # limb darkening file
n_z = args.n_z
## star parameter inputs 
omega = args.omega # dimensionless rotation speed
luminosity = args.luminosity # luminosity of the star in solar luminosities
mass = args.mass # mass of the star in solar masses
Req = args.Req # equatorial radius of the star in solar radii

## unpickle the limb darkening information
with open(pkl_lfile, 'rb') as f:
	ld = pickle.load(f)

# print a message
ut.printf ("Mass: " + str(mass) + " sun(s)\nLuminosity: "+ str(luminosity) + " sun(s)\n" +\
	"Equatorial radius: " + str(Req) + " solar radi(i/us)\nRotation speed: " + str(omega) +\
	" of the Keplerian angular velocity at the equator\n")
## For a star with given physical parameters, resolution of map, 
## and limb darkening information, pre-compute quantities that are needed for 
## integrating the starlight and are independent of the inclination
st = star.Star(omega, luminosity, mass, Req, n_z, ld)

### Pickle the star
with open(pkl_sfile, 'wb') as f:
	pickle.dump(st, f)
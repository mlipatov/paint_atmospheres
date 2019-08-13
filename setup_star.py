# Sets up the spectrum calculation for a star
import limbdark.limbdark as limbdark
import star.star as star
import sys
import time
import argparse
import pickle

parser = argparse.ArgumentParser(description="Example: \n" +\
	"python setup_spectrum.py \'spectrum.pkl\' \'limbdark.pkl\' \'roche\' " +\
	"0.8760 0.08683 40.124 2.135 2.818")
parser.add_argument("pkl_sfile", help="a name for a .pkl spectrum file to create")
parser.add_argument("pkl_lfile", help="a name for a limb darkening .pkl file from which to obtain wavelengths")
parser.add_argument("surface", help="surface of the star: either roche or ellipsoid")
parser.add_argument("omega", help="rotation speed divided by its Keplerian value at the equator", type=float)
parser.add_argument("inclination", help="angle between line of sight and rotation axis", type=float)
parser.add_argument("luminosity", help="luminosity in solar luminosities", type=float)
parser.add_argument("mass", help="mass in solar masses", type=float)
parser.add_argument("Req", help="equatorial radius in solar radii", type=float)
args = parser.parse_args()

### Inputs
pkl_sfile = args.pkl_sfile # spectrum file
pkl_lfile = args.pkl_lfile # limb darkening file
## star parameters 
surface = args.surface # functional approximation of the surface of the star
omega = args.omega # dimensionless rotation speed
inclination = args.inclination # the angle between the axis of rotation and the line of sight
luminosity = args.luminosity # luminosity of the star in solar luminosities
mass = args.mass # mass of the star in solar masses
Req = args.Req # equatorial radius of the star in solar radii

## unpickle the limb darkening information
with open(pkl_lfile, 'rb') as f:
	ld = pickle.load(f)

### Pickle the star
st = star.Star(omega, inclination, luminosity, mass, Req, surface)
with open(pkl_sfile, 'wb') as f:
	pickle.dump(st, f)

# print an initialization message
if (surface == "ellipsoid"):
	surfmes = " an ellipsoid"
elif (surface == "roche"):
	surfmes = "obtained from a Roche potential"
ct.printf ("Spectrum initialized for a star with a mass of " + str(mass) + " sun(s), a luminosity of "\
	+ str(luminosity) + " sun(s), an equatorial radius equal to " + str(Req) + \
	" solar radii, a rotation speed " + str(omega) +\
	" of the Keplerian angular velocity at the equator, and an inclination of " +\
	str(inclination) + " radians. The star's surface is " + surfmes + ".\n")
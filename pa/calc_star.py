# Sets up the spectrum calculation for a star
from pa.lib import fit as ft
from pa.lib import star
from pa.lib import util as ut
import argparse
import pickle

def run():
	parser = argparse.ArgumentParser(description="Example: \n" +\
		"calc_star \'data/limbdark_m01.pkl\' \'data/vega.pkl\' " +\
		"0.6151 40.346 2.165 2.815 100")
	parser.add_argument("ld_file", help="the limb darkening .pkl file to access")
	parser.add_argument("output", help="an output file containing the pickled star")
	parser.add_argument("omega", help="rotation speed divided by its Keplerian value at the equator", type=float)
	parser.add_argument("luminosity", help="luminosity in solar luminosities", type=float)
	parser.add_argument("mass", help="mass in solar masses", type=float)
	parser.add_argument("Req", help="equatorial radius in solar radii", type=float)
	parser.add_argument("nz", help="number of z coordinates in the upper half (normally 100)", type=int)
	parser.add_argument("-t", help="temperature interpolation: 0=planck(default), 1=linear, 2=log", type=int, \
			default=0)
	parser.add_argument("-n", help="number of steps in the Newton's method for temperature calculation (default is 15)", type=int, \
			default=15)
	args = parser.parse_args()

	## input files
	pkl_sfile = args.output # star pickle file
	pkl_lfile = args.ld_file # limb darkening file
	## star parameter inputs 
	omega = args.omega # dimensionless rotation speed
	luminosity = args.luminosity # luminosity of the star in solar luminosities
	mass = args.mass # mass of the star in solar masses
	Req = args.Req # equatorial radius of the star in solar radii
	# integration parameter
	nz = args.nz
	# temperature interpolation parameter
	if args.t == 0:
		tm = 'planck'
	if args.t == 1: 
		tm ='linear'
	elif args.t == 2:
		tm ='log'
	# Newton's method parameter
	nm = args.n 

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
	st = star.Star(omega, luminosity, mass, Req, nz=nz, ld=ld, temp_method=tm, nm=nm)

	### Pickle the star
	with open(pkl_sfile, 'wb') as f:
		pickle.dump(st, f)
		
# in case we are running this file as the main program
if __name__ == "__main__":
	run()
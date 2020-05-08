# This file calculates the fits of Intensity versus mu
# from the limb darkening information of Castelli and Kurucz 2004

from pa.lib import limbdark
import argparse
import pickle

def run():
	parser = argparse.ArgumentParser(description="Example: \n" +\
		"calc_limbdark \'data/im01k2.pck\' \'data/limbdark_m01.pkl\' -b 0.1 0.4 -s")
	parser.add_argument("ldfile", help="an existing file with limb darkening information")
	parser.add_argument("pkl_lfile", help="name for a .pkl file with limb darkening information to create")
	parser.add_argument('-b', type=float, nargs='+', help='list of boundaries between intervals', required=True)
	parser.add_argument("-s", help="save the original discrete intensities in the .pkl file", 
						action="store_true")
	args = parser.parse_args()

	pkl_lfile = args.pkl_lfile # limb darkening pickle file
	ldfile = args.ldfile # file with limb darkening information
	bounds = args.b # list of boundaries of separate intervals for the fits
	save = args.s # save the intensity data

	## Load the limb darkening information, calculate the fits and perform checks
	## (about 4 - 6 seconds to calculate the fits)
	ld = limbdark.LimbDark(ldfile, bounds, save)

	### Pickle the limb darkening information
	with open(pkl_lfile, 'wb') as f:
		pickle.dump(ld, f)
		
# in case we are running this file as the main program
if __name__ == "__main__":
	run()
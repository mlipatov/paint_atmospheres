# This file calculates the fits of Intensity versus mu
# from the limb darkening information of Castelli and Kurucz 2004

from pa.lib import limbdark
import argparse
import pickle
import glob

def run():
	parser = argparse.ArgumentParser(description="Examples: \n" +\
		"calc_limbdark data/im01k2.pck data/limbdark_m01f.pkl 0.1 0.4 -f data/filters/; " +\
		"calc_limbdark data/im01k2.pck data/limbdark_m01.pkl 0.1 0.4 -s")
	parser.add_argument("ldfile", help="an existing file with intensities on an atmosphere grid")
	parser.add_argument("pkl_lfile", help="name for an output .pkl file with limb darkening fits")
	parser.add_argument('bounds', type=float, nargs='+', help='list of boundaries between intervals')
	parser.add_argument("-f", help="directory with filter files")
	parser.add_argument("-s", help="save the discrete intensities in the .pkl file", 
						action="store_true")
	args = parser.parse_args()

	pkl_lfile = args.pkl_lfile # limb darkening pickle file
	ldfile = args.ldfile # file with limb darkening information
	bounds = args.bounds # list of boundaries of separate intervals for the fits
	save = args.s # whether to save the intensity data
	if args.f is not None:
		filtfiles = glob.glob(args.f + "*")
		filtfiles.sort()
	else:
		filtfiles = []

	## Load the limb darkening information, calculate the fits and perform checks
	## (about 4 - 6 seconds to calculate the fits)
	ld = limbdark.LimbDark(ldfile, filtfiles, bounds, save)

	### Pickle the limb darkening information
	with open(pkl_lfile, 'wb') as f:
		pickle.dump(ld, f)
		
# in case we are running this file as the main program
if __name__ == "__main__":
	run()
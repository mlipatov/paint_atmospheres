# This file calculates the fits of Intensity versus mu
# from the limb darkening information of Castelli and Kurucz 2004

from pa.lib import limbdark as ld
from pa.lib import util as ut
import argparse
import pickle
import glob, os

def run():
	parser = argparse.ArgumentParser(description="Examples: \n" +\
		"calc_limbdark data/im01k2.pck data/limbdark_m01f.pkl 0.1 0.4 -f data/filters/ -a 0.1; " +\
		"calc_limbdark data/im01k2.pck data/limbdark_m01.pkl 0.1 0.4 -s")
	parser.add_argument("ldfile", help="an existing file with intensities on an atmosphere grid")
	parser.add_argument("pkl_lfile", help="name for an output .pkl file with limb darkening fits")
	parser.add_argument('bounds', type=float, nargs='+', help='list of boundaries between intervals')
	parser.add_argument("-f", help="directory with filter files for photometry mode")
	parser.add_argument("-a", type=float, default=0, help="A_V reddening coefficient, " +\
		"currently only implemented in photometry mode; default is zero; R_V is 3.1")
	parser.add_argument("-s", help="save the discrete intensities in the .pkl file", 
						action="store_true")
	args = parser.parse_args()

	pkl_lfile = args.pkl_lfile # limb darkening pickle file
	ldfile = args.ldfile # file with limb darkening information
	bounds = args.bounds # list of boundaries of separate intervals for the fits
	save = args.s # whether to save the intensity data
	a = args.a

	# metallicity
	zstr = os.path.basename(ldfile)[1:4].replace('m', '-').replace('p', '')
	Z = float(zstr[:-1] + '.' + zstr[-1])

	lam, T, g = ld.getlamTg(ldfile) # wavelengths in the Kurucz file
	I = ld.getI(ldfile, lam, T, g) # intensity array from the file
	l = ld.LimbDark(lam, T, g, Z) # initialize a limb darkening object

	if args.f is not None:
		filtfiles = glob.glob(args.f + "*")
		filtfiles.sort()
		ut.setlam( lam ) # set the wavelength array for filtering
		I = l.filter(I, filtfiles, a_v=a) # filter and redden the object, return band intensities

	# fit polynomials for I(mu)
	l.fit(I, bounds, save) 

	### Pickle the limb darkening information
	with open(pkl_lfile, 'wb') as f:
		pickle.dump(l, f)
		
# in case we are running this file as the main program
if __name__ == "__main__":
	run()
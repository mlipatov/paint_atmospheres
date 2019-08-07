# This file calculates the fits of Intensity versus mu
# from the limb darkening information of Castelli and Kurucz 2004

import limbdark.limbdark as limbdark
import argparse
import pickle

parser = argparse.ArgumentParser(description="Example: \n" +\
	"python calc_limbdark.py \'im01k2.pck\' \'limbdark.pkl\' -c")
parser.add_argument("ldfile", help="an existing file with limb darkening information")
parser.add_argument("pkl_lfile", help="a name for a .pkl file with limb darkening information to create")
parser.add_argument("-c", help="check each fit for non-negativity, monotonicity and goodness", 
					action="store_true")
args = parser.parse_args()

pkl_lfile = args.pkl_lfile # limb darkening pickle file
ldfile = args.ldfile # file with limb darkening information
check = args.c # check the fit: true or false 

## Load the limb darkening information and calculate the fits
## (about 4 - 6 seconds)
ld = limbdark.LimbDark(ldfile, check)

### Pickle the limb darkening information
with open(pkl_lfile, 'wb') as f:
	pickle.dump(ld, f)
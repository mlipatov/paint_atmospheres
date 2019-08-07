import limbdark.limbdark as limbdark
import argparse
import pickle

parser = argparse.ArgumentParser(description="Example: \n" +\
	"python plot_Imu.py \'limbdark.pkl\' 401 0 3500")
parser.add_argument("pkl_lfile", help="a name of a .pkl file with limb darkening information")
parser.add_argument("wavelength", help="wavelength in nm", type=float)
parser.add_argument("logg", help="log base 10 of g in cgs units", type=float)
parser.add_argument("temp", help="temperature in Kelvin", type=float)
args = parser.parse_args()

pkl_lfile = args.pkl_lfile # limb darkening pickled file
wl = args.wavelength
g = args.logg
temp = args.temp

## unpickle the fits
with open(pkl_lfile, 'rb') as f:
	ld = pickle.load(f)

## Plot a dependence of intensity on mu
ld.plotFit(wl, g, temp)
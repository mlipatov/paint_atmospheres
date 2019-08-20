# Sets up integration of spectra of rotating stars
import util as ut
import sys
import time
import argparse
import pickle

parser = argparse.ArgumentParser(description="Example: \n" +\
	"python setup_integration.py \'ellip.pkl\' 1000 0.0005 -c")
parser.add_argument("pkl_efile", help="a name for a .pkl elliptic integration file to create")
parser.add_argument("m_num", help="number of m values for interpolation", type=int)
parser.add_argument("z_step", help="step for normalized z coordinate", type=float)
parser.add_argument("-c", help="check the elliptic integration object for accuracy", 
					action="store_true")
args = parser.parse_args()

pkl_efile = args.pkl_efile
m_num = args.m_num
z_step = args.z_step
check = args.c # check the fit: true or false 
### Initialize and pickle the elliptic integration object
el = ut.Ellip(m_num, z_step, check)
with open(pkl_efile, 'wb') as f:
	pickle.dump(el, f)
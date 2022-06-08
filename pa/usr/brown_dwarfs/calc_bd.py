# start with a rotating brown dwarf at some luminosity and equatorial radius
# compute the average of its equatorial and polar radii
# for each rotational velocity omega in a set of decreasing omegas
# 	reduce the equatorial radius to keep the average radius constant
# 	compute the spectra of a rotating brown dwarf at this velocity and a number of inclinations
# compute non-rotating brown dwarf spectra at luminosities between those corresponding to
# the equatorial and the polar temperatures of the fastest-rotating brown dwarf;
# keep the radius of the non-rotating dwarf equal to the constant average radius

import numpy as np
import os, pickle, glob, argparse, shutil
from sigfig import round as sf_round
from pa.lib import util as ut 

datadir = '../../../data/'

parser = argparse.ArgumentParser(description="Examples: \n" +\
		"python calc_bd.py bd_pkl/J0348-6022_880/ bd_spectra/J0348-6022_880/ 4.28643e-6 11 21 -o 0.1 0.5 9; " + 
		"python calc_bd.py bd_pkl/J0348-6022_few/ bd_spectra/J0348-6022_few/ 4.28643e-6 3 9 -o 0.417;")
parser.add_argument("pkl_dir", help="directory with pickle files")
parser.add_argument("spec_dir", help="directory with spectrum files")
parser.add_argument("luminosity", help="constant luminosity of the rotating dwarf, in solar luminosities", type=float)
parser.add_argument("inum", help="number of inclinations", type=int)
parser.add_argument("nrt", help="number of non-rotating dwarfs at different luminosities", type=int)
parser.add_argument('-o', type=float, nargs='+', help='either a single Roche-equivalent omega ' +
	'or a equally spaced values specified by minimum, maximum and number', required=True)
args = parser.parse_args()

# omegas
omegas = args.o 
lo = len(omegas)
if lo not in [1, 3]:
	sys.exit("Please specify either a single omega (one number) " +\
		"or a range specified by minimum, maximum and step (three numbers).")
elif lo == 1:
	omegas = np.array( omegas )
elif lo == 3:
	mi, ma, num = omegas
	omegas = np.linspace( mi, ma, num=int(num) )

mass = 0.041 # constant mass of all dwarfs, in solar masses
dist = 3.09e19 # constant distance of 10 pc, in cm
znum = 100 # constant number of z values in the upper half of the dwarf

### compute pickle files and spectra of rotating brown dwarfs at several rotational velocities and inclinations
lum_rot = args.luminosity # constant luminosity of rotating dwarfs, in solar luminosities
inum = args.inum # number of inclinations
pkl_dir = args.pkl_dir
spec_dir = args.spec_dir
# clear the directories
shutil.rmtree(datadir + pkl_dir)
shutil.rmtree(datadir + spec_dir)

# omax = args.omax # maximum omega
Remax = 0.093 # maximum equatorial radius in solar radii
f = 1 / (1 + 2/omegas[-1]**2) # oblateness according to the Roche model
Rp = Remax * (1 - f) # polar radius of the maximum omega model, according to the Roche model
# average radius will be the radius of the non-rotating models
R = np.around( (Rp + Remax) / 2, 4) # precision at four decimal places
# omegas = np.around(np.linspace(0, omax, 11), 3)[2:] # use this line for nine omegas
Req = 2 * R * (2 + omegas**2) / (4 + omegas**2) # corresponding equatorial radii that keep average radius constant
omstr = [str(o).replace('.', '')[:3] for o in omegas] # strings that tag omegas
pkls = [] # pickle files of rotating brown dwarfs
for i in range(len(omegas)):
	omdir = datadir + pkl_dir + 'rotating/'
	os.makedirs(omdir, exist_ok=True)
	pkl = omdir + 'J0348-6022_om' + omstr[i] + '.pkl' # pickle file name
	pkls.append(pkl)
	# compute the dwarf pickle file
	command = 'calc_star ' + datadir + 'limbdark_BD_00.pkl ' + pkl + ' ' + str(omegas[i]) + ' ' + str(lum_rot) + \
		' ' + str(mass) + ' ' + str(Req[i]) + ' ' + str(dist) + ' ' + str(znum)
	os.system(command)
	# compute the spectra at different inclinations
	omdir = datadir + spec_dir + 'rotating/omega' + omstr[i] + '/'
	os.makedirs(omdir, exist_ok=True) # make the directory for this omega if it does not exist
	for i in glob.glob(omdir + '*'): os.remove(i) # remove files from that directory if there are any
	command = "calc_spectra " + pkl + " " + omdir + " -i 0.0 1.5707963267948966 " + str(inum)
	os.system(command)

# fastest rotating model
with open(pkls[-1], 'rb') as f:
	st = pickle.load(f)
T = st.map.temp_up # temperatures across the upper portion of the star

# luminosities that correspond to the effective temperatures 
# at the equator and the poles, in solar luminosities
Le, Lp = T.take([0, -1])**4 * (4 * np.pi * R**2 * ut.sigma_sun)

# luminosities of the non-rotating models, slightly rounded
L = np.linspace(Le, Lp, args.nrt)
L = [float(sf_round(l, sigfigs=6, notation='sci')) for l in L]
# corresponding strings for filenames
Lstr = [str("{:e}".format(l)).replace('.', '')[:4] for l in L]

# compute the non-rotating brown dwarf models and their spectra
for i in range(len(L)):
	pkldir = datadir + pkl_dir + "non_rotating/"
	os.makedirs(pkldir, exist_ok=True)
	pkl = "J0348-6022_L" + Lstr[i] + ".pkl"
	command = "calc_star " + datadir + "limbdark_BD_00.pkl " + \
		pkldir + pkl + " 0 " + str(L[i]) + " " + str(mass) + " " + str(R) + " " + str(dist) + " 100"
	os.system(command)
	specdir = datadir + spec_dir + "non_rotating/"
	os.makedirs(specdir, exist_ok=True)
	command = "calc_spectra " + pkldir + pkl + " " + specdir + " -i 0.0"
	os.system(command)
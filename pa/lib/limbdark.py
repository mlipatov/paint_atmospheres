# from a file provided by Castelli and Kurucz 2004,
# this module obtains wavelengths, log g values, temperatures and intensities
from pa.lib import util as ut
from pa.lib import fit as ft
import numpy as np
import sys, time
import re

# given the name of a file provided by Castelli and Kurucz 2004,
# get the wavelengths in nm, sorted log g and temperature arrays
def getlamTg(filename):

    f = open(filename,'U')
    lam = np.empty(0, dtype=np.float64)
    count = 0;
    for line in f:
        if ('TEFF' in line):
            count += 1
        if count > 1: break
        if (not '#' in line) and line.strip():
            data = line.split()
            lam = np.append(lam, np.float(data[0]))
    f.close()

    # arrays of discrete temperatures and gravities
    T = np.empty(0, dtype=np.float16)
    g = np.empty(0, dtype=np.float16)

    f = open(filename,'U')
    for line in f:
        if ('TEFF' in line):        
            data = line.split()
            T = np.append(T, np.float32(data[1]))
            g = np.append(g, np.float16(data[3]))
    f.close()
    # make the temperatures and values of log g sorted and unique
    T = np.unique(T)
    T = np.sort(T)
    g = np.unique(g)
    g = np.sort(g)
    return [lam, T, g]
    
# For every combination of wavelength, mu, log g (called g below) and temperature, record the intensity
# Intensity indices: [wavelength][mu][log gravity][temperature]
# Intensity units: erg cm^-2 s^-1 Hz^-1 ster^-1 on a grid
# Where the intensity is not given for some temperature, a NAN
def getI(filename, lam, T, g):
    I = np.full((len(lam), len(ft.mu_arr), len(g), len(T)), np.nan, dtype=np.float64)
    il = 0; ig = 0; iT = 0
    lam = -1; mu = -1; gv = -1; Tv = -1;
    f = open(filename,'U')
    for line in f:
        if ('TEFF' in line):
            il = 0
            data = line.split()
            Tv = float(data[1])
            iT = np.where(T == Tv)[0][0]
            gv = float(data[3])
            ig = np.where(g == gv)[0][0]
        if (not '#' in line) and line.strip():
            # split the data line according to the number of characters alotted to each field
            # wavelength:   9
            # I(mu = 1):    10
            # I(mu < 1):    6
            line = line.rstrip() # remove trailing spaces
            data = []
            data.extend([ line[0:9], line[9:19] ])
            line = line[19:]
            while len(line) > 5:
                data.append(line[:6])
                line = line[6:]
            data = np.array(data)
            # intensity at mu = 1
            I1 = float(data[1]) 
            I[il, -1, ig, iT] = I1
            # intensities at other values of mu
            I_rest = I1 * np.array([float(x) for x in data[2:]])/100000.0 
            # flip the intensity array because the mu values are flipped
            I[il, 0:-1, ig, iT] = np.flip(I_rest)
            il += 1
    f.close()
    return I


class LimbDark:
    """ 
    Limb darkening fits.

    Uses python's numpy.float32, which is precise to 6 decimal digits; 
    Kurucz gives 5 decimal digits for each intensity, so that the derived 
    coefficients of intensity fits are only meaningful to 5 decimal digits.
    
    Admits two configurations:
    1. Spectrum. Contains fits of specific intensities to viewing angle, at every T, g and lambda.
    2. Photometry. Contains fits of filtered intensities to viewing angle, at every T, g and band.
        Also contains the additive magnitude offset for every band.
    Optionally, contains the original intensity grid at every T, g, lambda/band and mu.

     """

    # initialize an object
    # Inputs:
    #   wavelengths
    #   metallicity
    # Fields modified:
    #   gravity and temperature arrays
    #   metallicity
    def __init__(self, lam, T, g, Z):
        self.g = g
        self.T = T
        self.lam = lam
        self.Z = Z
        self.a_v = np.array([0]) # default reddening
        self.bands = None
        self.F0 = None

    # filter and redden the intensities (photometry mode)
    # Inputs:
    #   array of intensities on a grid of [wavelength][mu][log gravity][temperature]
    #   list of files with band transmission curves and flux zero points (photometry mode)
    #   optional set of reddening coefficients A_V, equal to zero if no reddening
    # Fields modified:
    #   reddening coefficients: (av0, av1, ..., avm)
    #   band names:             (b0, b1, ..., bn)
    #   magnitude offset for each band: (F0 x m), (F1 x m), ... (Fn x m)
    #   characteristic wavelength for each band: (lam0 x m), (lam1 x m), ... (lamn x m)
    # Output:
    #   intensity array, on a grid of [band x reddening][mu][log gravity][temperature],
    #       where the first dimension corresponds to pairs 
    #       (b0, av0), (b0, av1), ..., (b0, avm), (b1, av0), ..., (bn, avm)
    def filter(self, I, filtfiles, a_v=[0]):
        anum = len(a_v) # number of A_V values
        bands = [] # band names
        F0 = [] # flux zero points

        # array with band intensities
        Iband = []
        # make the wavelength the last dimension in the intensity array
        I = np.moveaxis(I, 0, -1)
        # band parameters
        wl_band = [] # wavelengths
        for file in filtfiles: # for every band
            f = open(file)
            filetext = f.read()
            f.close()
            # band name
            bands.append( re.search("# name.*", filetext)[0].split(':')[-1].strip() )
            # zero pt flux, converted from erg cm^-2 s^-1 A^-1 to erg cm^-2 s^-1 nm^-1
            F0.append( float(re.search("# flux zero point.*", filetext)[0].split(':')[-1].strip()) * 10. )
            # characteristic wavelength, converted from A to nm
            wl_band.append( float(re.search("# mean wavelength.*", filetext)[0].split(':')[-1].strip()) / 10. )
            # transmission curve, wavelengths converted from A to nm
            wlf, T = np.loadtxt(file).T
            wlf = wlf / 10.
            # intensities in erg cm^-2 s^-1 nm^-1 ster^-1, 
            # with the wavelength dimension filtered out; reddened
            for a in a_v:
                Iband.append( ut.filter(I, T, wlf, a) )
        I = np.array( Iband ) # a filtered intensity at each band

        self.a_v = np.array(a_v) # record the reddening coefficients
        self.bands = np.array(bands) # record band names
        self.lam = np.repeat(wl_band, anum) # characteristic band wavelengths, repeated for each reddening
        self.F0 = np.repeat(F0, anum) # flux magnitude zero points, repeated for each reddening
        return I # return the resulting discrete intensities

    # fit intensities to polynomials vs mu
    # Inputs:
    #   bounds that define a partition of mu's range
    #   optional: whether to save the discrete intensities
    # Fields modified:
    #   mu partition bounds
    #   fit parameters
    def fit(self, I, bounds, save=False):
        self.bounds = bounds
        if save:
            self.I = I
        # initialize the module that computes the fits
        ft.set_muB(bounds)

        n_temp = len(self.T)
        n_g = len(self.g)
        n_wl = len(self.lam)
        n_param = ft.n * ft.m
        # initialize a list of fit coefficients for interpolation to NAN values for each 
        # combination of gravity, temperature, wavelength and the fit parameter index: 
        # index 1 : temperature;
        # index 2 : log g;
        # index 3 : wavelength (band);
        # index 4 : fit parameter index (if applicable).
        self.fit_params = \
            np.full( (n_temp, n_g, n_wl, n_param), np.nan, dtype=np.float32 ) 
        
        # for each combination of wavelength, gravity and temperature, calculate the fit at these values; 
        # record the fit parameters in the array needed for interpolation
        
        print ("Computing fits of intensity versus mu. ")
        sys.stdout.flush()
        start = time.time()
        for ind_g in np.arange(n_g):
            ut.printf(str(ind_g) + " out of " + str(n_g) + " gravity values completed.\n")
            sys.stdout.flush()
            for ind_temp in np.arange(n_temp):
                if not np.isnan(I[0, 0, ind_g, ind_temp]):
                    fp = np.empty( (n_wl, n_param) )
                    for ind_wl in np.arange(n_wl):
                        I_slice = I[ind_wl, :, ind_g, ind_temp] # get the intensities at different mus
                        wl = self.lam[ind_wl]
                        g = self.g[ind_g]
                        temp = self.T[ind_temp]
                        # fit and record the fit parameters in the array that is later used by interpolation
                        fp[ind_wl] = ft.fit(I_slice)
                    self.fit_params[ind_temp][ind_g] = fp
        end = time.time()
        print("Done in " + str(end - start) + " seconds")
        sys.stdout.flush()
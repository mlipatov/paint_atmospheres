from pa.lib import util as ut
from pa.lib import fit as ft
import numpy as np
import sys, time

# given the name of a file provided by Castelli and Kurucz 2004,
# output arrays of wavelengths, mu, values of log g, temperatures, intensities and 
#     an array specifying which values of temperature each log g has 
#    (this contains a -1 for each combination of temperature and log g that doesn't exist)
# all arrays other than the intensity array should be sorted
# Intensity indices: [wavelength][mu][log gravity][temperature]
# Where the intensity is not given for some temperature, the record is -1
def getdata(filename):
    # get the wavelengths
    f = open(filename,'U')
    wl_arr = np.empty(0, dtype=np.float64)
    count = 0;
    for line in f:
        if ('TEFF' in line):
            count += 1
        if count > 1: break
        if (not '#' in line) and line.strip():
            data = line.split()
            wl_arr = np.append(wl_arr, np.float(data[0]))
    f.close()
    # get the values of log g and temperature
    f = open(filename,'U')
    g_arr = np.empty(0, dtype=np.float16) 
    temp_arr = np.empty(0, dtype=np.float32)
    for line in f:
        if ('TEFF' in line):        
            data = line.split()
            temp = np.float32(data[1])
            temp_arr = np.append(temp_arr, temp)
            g = np.float16(data[3])
            g_arr = np.append(g_arr, g)
    f.close()
    # make the temperatures and values of log g sorted and unique
    temp_arr = np.unique(temp_arr)
    temp_arr = np.sort(temp_arr)
    g_arr = np.unique(g_arr)
    g_arr = np.sort(g_arr)

    # For every combination of wavelength, mu, log g (called g below) and temperature, record the intensity
    # Intensity indices: [wavelength][mu][log gravity][temperature]
    # Where the intensity is not given for some temperature, a NAN
    I_arr = np.full((len(wl_arr), len(ft.mu_arr), len(g_arr), len(temp_arr)), np.nan, dtype=np.float64)
    ind_wl = 0; ind_g = 0; ind_temp = 0
    wl = -1; mu = -1; g = -1; temp = -1;
    f = open(filename,'U')
    for line in f:
        if ('TEFF' in line):
            ind_wl = 0
            data = line.split()
            temp = float(data[1])
            ind_temp = np.where(temp_arr == temp)[0][0]
            g = float(data[3])
            ind_g = np.where(g_arr == g)[0][0]
        if (not '#' in line) and line.strip():
            data = line.split()
            I1 = float(data[1]) # intensity at mu = 1
            I_rest = I1 * np.array([float(x) for x in data[2:]])/100000.0 # intensities at other values of mu
            # flip the intensity array because the mu values are flipped
            I_arr[ind_wl, -1, ind_g, ind_temp] = I1
            I_arr[ind_wl, 0:-1, ind_g, ind_temp] = np.flip(I_rest)
            ind_wl += 1
    f.close()

    return [wl_arr, g_arr, temp_arr, I_arr]


class LimbDark:
    """ Class containing all limbdarkening information.
    Uses single-precision floating point to conserve space. Specifically, uses python's numpy.float32, which
    is precise to 6 decimal digits; Kurucz gives 5 decimal digits for each intensity, so that the derived 
    coefficients of intensity fits are only meaningful to 5 decimal digits. """

    # initialize with a file containing the limb darkening information from Castelli and Kurucz 2004
    def __init__(self, datafile, bounds, check, save):
        wl_arr, g_arr, temp_arr, I_arr = getdata(datafile)
        self.wl_arr = wl_arr
        self.g_arr = g_arr
        self.temp_arr = temp_arr
        self.bounds = bounds
        if save:
            self.I_arr = I_arr
        ft.Fit.set_muB(bounds)

        n_temp = len(temp_arr)
        n_g = len(g_arr)
        n_wl = len(wl_arr)
        n_param = ft.n * ft.Fit.m
        # initialize a list of fit coefficients for interpolation to NAN values for each 
        # combination of gravity, temperature, wavelength and the fit parameter index: 
        # index 1 : temperature;
        # index 2 : log g;
        # index 3 : wavelength; (if applicable)
        # index 4 : fit parameter index (if applicable).
        self.fit_params = \
            np.full( (n_temp, n_g, n_wl, n_param), np.nan, dtype=np.float32 ) 
        # for each combination of wavelength, gravity and temperature, initialize a fit object, 
        # thus calculating the fit at these values; record the fit parameters in the array needed for 
        # interpolation
        print ("Computing fits of intensity versus mu. ")
        sys.stdout.flush()
        start = time.time()
        for ind_g in np.arange(n_g):
            ut.printf(str(ind_g) + " out of " + str(n_g) + " gravity values completed.\n")
            sys.stdout.flush()
            for ind_temp in np.arange(n_temp):
                if not np.isnan(I_arr[0, 0, ind_g, ind_temp]):
                    fp = np.empty( (n_wl, n_param) )
                    for ind_wl in np.arange(n_wl):
                        I_slice = I_arr[ind_wl, :, ind_g, ind_temp] # get the intensities at different mus
                        wl = wl_arr[ind_wl]
                        g = g_arr[ind_g]
                        temp = temp_arr[ind_temp]
                        # initialize and fit
                        fit = ft.Fit(I_slice, wl, g, temp, check)
                        # record the fit parameters in the array that is later used by interpolation
                        fp[ind_wl] = fit.p
                    self.fit_params[ind_temp][ind_g] = fp
            if check:
                print (ft.Fit.I0_min, ft.Fit.min_step, ft.Fit.max_dev)
        end = time.time()
        print("Done in " + str(end - start) + " seconds")
        sys.stdout.flush()


    # plots the information corresponding to given wavelength, log g and temperature
    def plotFit(self, wl, g, temp):
        # create a fit object
        ind_wl = np.where(self.wl_arr == wl)[0][0]
        ind_g = np.where(self.g_arr == g)[0][0]
        ind_temp = np.where(self.temp_arr == temp)[0][0]
        I_slice = self.I_arr[ind_wl, :, ind_g, ind_temp]
        check = False
        fit = ft.Fit(I_slice, wl, g, temp, check)
        # plot it
        fit.plot()
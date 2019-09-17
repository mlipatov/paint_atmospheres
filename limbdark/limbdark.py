import util as ut
import limbdark.fit as ft
import numpy as np
import sys, time
import math

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
    wl_arr = np.empty(0)
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
    g_arr = np.empty(0)
    temp_arr = np.empty(0)
    for line in f:
        if ('TEFF' in line):        
            data = line.split()
            temp = np.float(data[1])
            temp_arr = np.append(temp_arr, temp)
            g = np.float(data[3])
            g_arr = np.append(g_arr, g)
    f.close()
    # make the temperatures and values of log g sorted and unique
    temp_arr = np.unique(temp_arr)
    temp_arr = np.sort(temp_arr)
    g_arr = np.unique(g_arr)
    g_arr = np.sort(g_arr)
    ## create an array that records which values of temperature each log g has
    # initialize the array elements to -1
    g_temp_arr = -1*np.ones((len(g_arr), len(temp_arr)))
    # populate the array
    f = open(filename,'U')
    for line in f:
        if ('TEFF' in line):
            data = line.split()
            temp = np.float(data[1])
            g = np.float(data[3])
            # find out which index of the g array has this g
            ind_g = np.where(g_arr == g)
            # find out which index of the temperature array has this temperature
            ind_t = np.where(temp_arr == temp)
            if (not ind_g[0].size==0) and (not ind_t[0].size==0): # if there are such values of g and temp
                ind_g = ind_g[0][0] # get the indices
                ind_t = ind_t[0][0]
                # modify the corresponding element of the array being populated
                g_temp_arr[ind_g][ind_t] = temp
    f.close()

    # For every combination of wavelength, mu, log g (called g below) and temperature, record the intensity
    # Intensity indices: [wavelength][mu][log gravity][temperature]
    # Where the intensity is not given for some temperature, record -1
    I_arr = -1*np.ones((len(wl_arr), len(ft.mu_arr), len(g_arr), len(temp_arr)))
    ind_wl = 0; ind_mu = 0; ind_g = 0; ind_temp = 0
    wl = -1; mu = -1; g = -1; temp = -1;
    f = open(filename,'U')
    for line in f:
        if ('TEFF' in line):
            ind_wl = 0
            data = line.split()
            temp = np.float(data[1])
            ind_temp = np.where(temp_arr == temp)[0][0]
            g = np.float(data[3])
            ind_g = np.where(g_arr == g)[0][0]
        if (not '#' in line) and line.strip():
            data = line.split()
            I1 = np.float(data[1]) # intensity at mu = 1
            I_rest = I1 * np.array([float(x) for x in data[2:]])/100000.0 # intensities at other values of mu
            # flip the intensity array because the mu values are flipped
            I_arr[ind_wl, -1, ind_g, ind_temp] = I1
            I_arr[ind_wl, 0:-1, ind_g, ind_temp] = np.flip(I_rest)
            ind_wl += 1
    f.close()

    return [wl_arr, g_arr, temp_arr, I_arr, g_temp_arr]


class LimbDark:
    """ Class containing all limbdarkening information """

    # initialize with a file containing the limb darkening information from Castelli and Kurucz 2004
    def __init__(self, datafile, bounds, check, save):
        wl_arr, g_arr, temp_arr, I_arr, g_temp_arr = getdata(datafile)
        self.wl_arr = wl_arr
        self.g_arr = g_arr
        self.temp_arr = temp_arr
        # self.g_temp_arr = g_temp_arr # deprecated, use None in the parameter array instead
        self.bounds = bounds
        if save:
            self.I_arr = I_arr
        ft.Fit.set_muB(bounds)

        g_temp_shape = np.shape(g_temp_arr)
        n_wl = len(wl_arr)
        n_param = ft.n * ft.Fit.m
        # initialize a list of fit coefficients for interpolation to NAN values for each 
        # combination of gravity and temperature; for the combinations where coefficients exist,
        # the last two indices will stand for the wavelength and the fit parameter index: 
        # index 1 : log g;
        # index 2 : temperature;
        # index 3 : wavelength; (if applicable)
        # index 4 : fit parameter index (if applicable).
        fit_params = [[None for j in range(g_temp_shape[1])] for i in range(g_temp_shape[0])]
        # for each combination of wavelength, gravity and temperature, initialize a fit object, 
        # thus calculating the fit at these values; record the fit parameters in the array needed for 
        # interpolation
        print ("Computing fits of intensity versus mu. ")
        sys.stdout.flush()
        start = time.time()
        for ind_g in np.arange(g_temp_shape[0]):
            ut.printf(str(ind_g) + " out of " + str(g_temp_shape[0]) + " gravity values completed.\n")
            sys.stdout.flush()
            for ind_temp in np.arange(g_temp_shape[1]):
                if g_temp_arr[ind_g][ind_temp] != -1:
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
                    fit_params[ind_g][ind_temp] = fp
            if check:
                print (ft.Fit.I0_min, ft.Fit.min_step, ft.Fit.max_dev)
        self.fit_params = fit_params
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
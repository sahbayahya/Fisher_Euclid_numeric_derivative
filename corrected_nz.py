

import numpy as np
import pylab as P
import scipy.integrate
import scipy.interpolate
import scipy.optimize
import sys
from fit_Euclid import dvdz


DEBUG_PLOT = True # Whether to plot fitting functions or not
NU_21 = 1.420 # HI emxission line freq. in GHz
FULLSKY = (4.*np.pi * (180./np.pi)**2.) # deg^2 in the full sky
NGAL_MIN = 1e3 # Min. no. of galaxies to tolerate in a redshift bin
CBM = 1. #np.sqrt(1.57) # Correction factor due to effective beam for MID/MK (OBSOLETE)
CTH = 0.5 # Correction factor due to taking 5 sigma (not 10 sigma) cuts for SKA1
SBIG = 500. # Flux rms to extrapolate dn/dz out to (constrains behaviour at large Srms)
C = 3e5 # Speed of light, km/s

# Define fitting coefficients from Mario's note (HI_specs.pdf)
Srms = np.array([0., 1., 3., 5., 6., 7.3, 10., 23., 40., 70., 100., 150., 200.,])
c1 = [6.23, 7.33, 6.91, 6.77, 6.84, 6.76, 6.64, 6.02, 5.74, 5.62, 5.63, 5.48, 5.00]
c2 = [1.82, 3.02, 2.38, 2.17, 2.23, 2.14, 2.01, 1.43, 1.22, 1.11, 1.41, 1.33, 1.04]
c3 = [0.98, 5.34, 5.84, 6.63, 7.13, 7.36, 7.95, 9.03, 10.58, 13.03, 15.49, 16.62, 17.52]
c4 = [0.8695, 0.5863, 0.4780, 0.5884, 0.5908, 0.5088, 0.4489, 0.5751, 0.5125,
      0.6193, 0.6212, 1., 1., 1.]
c5 = [0.2338, 0.6410, 0.9181, 0.8076, 0.8455, 1.0222, 1.2069, 0.9993, 1.1842,
      1.0179, 1.0759, 0., 0., 0.]
c1 = np.array(c1); c2 = np.array(c2); c3 = np.array(c3)
c4 = np.array(c4); c5 = np.array(c5)
Smax = np.max(Srms)



def flux_redshift(z):
    """
    Flux rms as a function of redshift.
    """
    z = np.atleast_1d(z)
    nu = NU_21 / (1. + z)
    if nucrit[ID] is not None:
        Sz = fluxrms[ID] * Scorr[ID] * np.ones(nu.size)
        idxs = np.where(nu*1e3 > nucrit[ID])
        Sz[idxs] *= (nu[idxs]*1e3 / nucrit[ID])
    else:
        Sz = fluxrms[ID] * Scorr[ID]
        Sz = nu * Sz if not Sconst[ID] else Sz * np.ones(nu.size)
    return Sz
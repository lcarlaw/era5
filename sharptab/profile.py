''' Create the Sounding (Profile) Object '''
from numba import njit
from numba.experimental import jitclass
from numba import float32, int32, float64

import numpy as np

from sharptab import thermo
from sharptab import utils
from .constants import *

###############################################################################
# Jitted form of sharppy.sharptab.profile.
#
# This code has been heavily re-written and subsequently simplitied to work
# with numba/jit.
#
# For example, **kwargs are not allowed and had to make further refinements
# to remove calls to masked numpy arrays. This should be okay for our purposes
# since there should be no 'missing' data.
#
# In addition, there is some odd behavior in the original code with try/catch
# blocks. These don't necessarily function as anticipated with numba, so
# several "hacks" had to be put in place to specifically call sections of
# various functions...the main one being the wetlift->wobf->satlift sections.
###############################################################################

@njit
def create_profile(pres, tmpc, dwpc, wspd, wdir, hght):
    return Profile(pres, tmpc, dwpc, wspd, wdir, hght)

spec = [
        ('pres', float64[:]),
        ('tmpc', float64[:]),
        ('dwpc', float64[:]),
        ('wspd', float64[:]),
        ('wdir', float64[:]),
        ('hght', float64[:]),
        ('logp', float64[:]),
        ('vtmp', float64[:]),
        ('sfc', int32),
        ('top', int32),
        ('wetbulb', float64[:]),
        ('theta_', float64[:]),
        ('thetae', float64[:]),
        ('wvmr', float64[:]),
        ('relh', float64[:]),
        ('u', float64[:]),
        ('v', float64[:])
]
@jitclass(spec)
class Profile(object):
    def __init__(self, pres, tmpc, dwpc, wspd, wdir, hght):
    #def __init__(self, pres):
        ## get the data and turn them into arrays
        self.pres = pres
        self.hght = hght
        self.tmpc = tmpc
        self.dwpc = dwpc
        self.wdir = wdir
        self.wspd = wspd

        self.u, self.v = utils.vec2comp(self.wdir, self.wspd)

        self.logp = np.log10(self.pres.copy())
        self.vtmp = thermo.virtemp( self.pres, self.tmpc, self.dwpc )


        #idx = np.where(self.pres > 0)[0]
        #self.vtmp[self.dwpc.mask[idx]] = self.tmpc[self.dwpc.mask[idx]] # Masking any virtual temperature

        ## get the index of the top and bottom of the profile
        self.sfc = np.where(self.tmpc)[0].min()
        self.top = np.where(self.tmpc)[0].max()

        ## generate the wetbulb profile
        self.wetbulb = self.get_wetbulb_profile()
        ## generate theta-e profile
        self.thetae = self.get_thetae_profile()
        ## generate theta profile
        self.theta_ = self.get_theta_profile()
        ## generate water vapor mixing ratio profile
        self.wvmr = self.get_wvmr_profile()
        ## generate rh profile
        self.relh = self.get_rh_profile()

    def get_rh_profile(self):
        '''
        Function to calculate the relative humidity profile

        Parameters
        ----------
        None

        Returns
        -------
        Array of the relative humidity profile
        '''

        rh = thermo.relh(self.pres, self.tmpc, self.dwpc)
        return rh

    def get_wvmr_profile(self):
        '''
            Function to calculate the water vapor mixing ratio profile.

            Parameters
            ----------
            None

            Returns
            -------
            Array of water vapor mixing ratio profile
            '''

        wvmr = thermo.calc_mixratio( self.pres, self.dwpc )
        return wvmr

    def get_thetae_profile(self):
        '''
            Function to calculate the theta-e profile.

            Parameters
            ----------
            None

            Returns
            -------
            Array of theta-e profile
            '''
        thetae = np.empty(self.pres.shape[0])
        for i in range(len(self.pres)):
            thetae[i] = thermo.ctok( thermo.calc_thetae(self.pres[i], self.tmpc[i], self.dwpc[i]) )
        return thetae

    def get_theta_profile(self):
        '''
            Function to calculate the theta profile.

            Parameters
            ----------
            None

            Returns
            -------
            Array of theta profile
            '''
        theta_ = np.empty(self.pres.shape[0])
        for i in range(len(self.pres)):
            theta_[i] = thermo.theta(self.pres[i], self.tmpc[i])
        theta_ = thermo.ctok(theta_)
        return theta_

    def get_wetbulb_profile(self):
        '''
            Function to calculate the wetbulb profile.

            Parameters
            ----------
            None

            Returns
            -------
            Array of wet bulb profile
            '''

        wetbulb = np.empty(self.pres.shape[0])
        for i in range(len(self.pres)):
            wetbulb[i] = thermo.calc_wetbulb( self.pres[i], self.tmpc[i], self.dwpc[i] )
        return wetbulb
"""
###############################################################################
#
# THERMODYNAMIC FUNCTIONS SCRAPED FROM SHARPY.SHARPTAB.THERMO
#
###############################################################################
@njit
def relh(p, t, td):
    '''
    Returns the virtual temperature (C) of a parcel.

    Parameters
    ----------
    p : number
        The pressure of the parcel (hPa)
    t : number
        Temperature of the parcel (C)
    td : number
        Dew point of parcel (C)

    Returns
    -------
    Relative humidity (%) of a parcel

    '''
    return 100. * vappres(td)/vappres(t) #$mixratio(p, td) / mixratio(p, t)

@njit
def calc_mixratio(p, t):
    '''
    Returns the mixing ratio (g/kg) of a parcel

    Parameters
    ----------
    p : number, numpy array
        Pressure of the parcel (hPa)
    t : number, numpy array
        Temperature of the parcel (hPa)

    Returns
    -------
    Mixing Ratio (g/kg) of the given parcel

    '''
    x = 0.02 * (t - 12.5 + (7500. / p))
    wfw = 1. + (0.0000045 * p) + (0.0014 * x * x)
    fwesw = wfw * vappres(t)
    return 621.97 * (fwesw / (p - fwesw))

@njit
def calc_thetae(p, t, td):
    '''
    Returns the equivalent potential temperature (C) of a parcel.

    Parameters
    ----------
    p : number
        The pressure of the parcel (hPa)
    t : number
        Temperature of the parcel (C)
    td : number
        Dew point of parcel (C)

    Returns
    -------
    Equivalent potential temperature (C)

    '''
    p2, t2 = drylift(p, t, td)
    return theta(100., wetlift(p2, t2, 100.))

@njit
def ctok(t):
    '''
    Convert temperature from Celsius to Kelvin

    Parameters
    ----------
    t : number, numpy array
        The temperature in Celsius

    Returns
    -------
    Temperature in Kelvin (number or numpy array)

    '''
    return t + ZEROCNK

@njit
def calc_wetbulb(p, t, td):
    '''
    Calculates the wetbulb temperature (C) for the given parcel

    Parameters
    ----------
    p : number
        Pressure of parcel (hPa)
    t : number
        Temperature of parcel (C)
    td : number
        Dew Point of parcel (C)

    Returns
    -------
    Wetbulb temperature (C)
    '''
    p2, t2 = drylift(p, t, td)
    return wetlift(p2, t2, p)

@njit
def wetlift(p, t, p2):
    '''
    Lifts a parcel moist adiabatically to its new level.

    Parameters
    -----------
    p : number
        Pressure of initial parcel (hPa)
    t : number
        Temperature of initial parcel (C)
    p2 : number
        Pressure of final level (hPa)

    Returns
    -------
    Temperature (C)

    '''
    thta = theta(p, t)
    thetam = thta - wobf(thta) + wobf(t)
    return satlift(p2, thetam)

@njit
def wobf(t):
    '''
    Implementation of the Wobus Function for computing the moist adiabats.

    .. caution::
        The Wobus function has been found to have a slight
        pressure dependency (Davies-Jones 2008).  This dependency
        is not included in this implementation.

    Parameters
    ----------
    t : number, numpy array
        Temperature (C)

    Returns
    -------
    Correction to theta (C) for calculation of saturated potential temperature.

    '''
    t = t - 20
    if t <= 0:
        npol = 1. + t * (-8.841660499999999e-3 + t * ( 1.4714143e-4 + t * (-9.671989000000001e-7 + t * (-3.2607217e-8 + t * (-3.8598073e-10)))))
        npol = 15.13 / (np.power(npol,4))
        return npol
    else:
        ppol = t * (4.9618922e-07 + t * (-6.1059365e-09 + t * (3.9401551e-11 + t * (-1.2588129e-13 + t * (1.6688280e-16)))))
        ppol = 1 + t * (3.6182989e-03 + t * (-1.3603273e-05 + ppol))
        ppol = (29.93 / np.power(ppol,4)) + (0.96 * t) - 14.8
        return ppol

@njit
def satlift(p, thetam):
    '''
    Returns the temperature (C) of a saturated parcel (thm) when lifted to a
    new pressure level (hPa)

    .. caution::
        Testing of the SHARPpy parcel lifting routines has revealed that the
        convergence criteria used the SHARP version (and implemented here) may cause
        drifting the pseudoadiabat to occasionally "drift" when high-resolution
        radiosonde data is used.  While a stricter convergence criteria (e.g. 0.01) has shown
        to resolve this problem, it creates a noticable departure from the SPC CAPE values and therefore
        may decalibrate the other SHARPpy functions (e.g. SARS).

    Parameters
    ----------
    p : number
        Pressure to which parcel is raised (hPa)
    thetam : number
        Saturated Potential Temperature of parcel (C)
    conv : number
        Convergence criteria for satlift() (C)

    Returns
    -------
    Temperature (C) of saturated parcel at new level

    '''
    #try:
    # If p and thetam are scalars
    if np.fabs(p - 1000.) - 0.001 <= 0:
        return thetam
    eor = 999
    t2 = None       # All vars must be pre-defined for njit
    e2 = None       # All vars must be pre-defined for njit
    while np.fabs(eor) - 0.1 > 0:
        if eor == 999:                  # First Pass
            pwrp = np.power((p / 1000.),ROCP)
            t1 = (thetam + ZEROCNK) * pwrp - ZEROCNK
            e1 = wobf(t1) - wobf(thetam)
            rate = 1
        else:                           # Successive Passes
            rate = (t2 - t1) / (e2 - e1)
            t1 = t2
            e1 = e2
        t2 = t1 - (e1 * rate)
        e2 = (t2 + ZEROCNK) / pwrp - ZEROCNK
        e2 += wobf(t2) - wobf(e2) - thetam
        eor = e2 * rate
    return t2 - eor

@njit
def drylift(p, t, td):
    '''
    Lifts a parcel to the LCL and returns its new level and temperature.

    Parameters
    ----------
    p : number, numpy array
        Pressure of initial parcel in hPa
    t : number, numpy array
        Temperature of inital parcel in C
    td : number, numpy array
        Dew Point of initial parcel in C

    Returns
    -------
    p2 : number, numpy array
        LCL pressure in hPa
    t2 : number, numpy array
        LCL Temperature in C

    '''
    t2 = lcltemp(t, td)
    p2 = thalvl(theta(p, t), t2)
    return p2, t2

@njit
def lcltemp(t, td):
    '''
    Returns the temperature (C) of a parcel when raised to its LCL.

    Parameters
    ----------
    t : number, numpy array
        Temperature of the parcel (C)
    td : number, numpy array
        Dewpoint temperature of the parcel (C)

    Returns
    -------
    Temperature (C) of the parcel at it's LCL.

    '''
    s = t - td
    dlt = s * (1.2185 + 0.001278 * t + s * (-0.00219 + 1.173e-5 * s -
        0.0000052 * t))
    return t - dlt

@njit
def thalvl(theta, t):
    '''
    Returns the level (hPa) of a parcel.

    Parameters
    ----------
    theta : number, numpy array
        Potential temperature of the parcel (C)
    t : number, numpy array
        Temperature of the parcel (C)

    Returns
    -------
    Pressure Level (hPa [float]) of the parcel
    '''

    t = t + ZEROCNK
    theta = theta + ZEROCNK
    return 1000. / (np.power((theta / t),(1./ROCP)))

@njit
def theta(p, t):
    '''
    Returns the potential temperature (C) of a parcel.

    Parameters
    ----------
    p : number, numpy array
        The pressure of the parcel (hPa)
    t : number, numpy array
        Temperature of the parcel (C)
    p2 : number, numpy array (default 1000.)
        Reference pressure level (hPa)

    Returns
    -------
    Potential temperature (C)

    '''
    return ((t + ZEROCNK) * np.power((1000. / p),ROCP)) - ZEROCNK

@njit
def virtemp(p, t, td):
    '''
    Returns the virtual temperature (C) of a parcel.  If
    td is masked, then it returns the temperature passed to the
    function.

    Parameters
    ----------
    p : number
        The pressure of the parcel (hPa)
    t : number
        Temperature of the parcel (C)
    td : number
        Dew point of parcel (C)

    Returns
    -------
    Virtual temperature (C)

    '''

    tk = t + ZEROCNK
    w = 0.001 * mixratio(p, td)
    vt = (tk * (1. + w / eps) / (1. + w)) - ZEROCNK
    return vt

@njit
def mixratio(p, t):
    '''
    Returns the mixing ratio (g/kg) of a parcel

    Parameters
    ----------
    p : number, numpy array
        Pressure of the parcel (hPa)
    t : number, numpy array
        Temperature of the parcel (hPa)

    Returns
    -------
    Mixing Ratio (g/kg) of the given parcel

    '''
    x = 0.02 * (t - 12.5 + (7500. / p))
    wfw = 1. + (0.0000045 * p) + (0.0014 * x * x)
    fwesw = wfw * vappres(t)
    return 621.97 * (fwesw / (p - fwesw))

@njit
def vappres(t):
    '''
    Returns the vapor pressure of dry air at given temperature

    Parameters
    ------
    t : number, numpy array
        Temperature of the parcel (C)

    Returns
    -------
    Vapor Pressure of dry air

    '''
    pol = t * (1.1112018e-17 + (t * -3.0994571e-20))
    pol = t * (2.1874425e-13 + (t * (-1.789232e-15 + pol)))
    pol = t * (4.3884180e-09 + (t * (-2.988388e-11 + pol)))
    pol = t * (7.8736169e-05 + (t * (-6.111796e-07 + pol)))
    pol = 0.99999683 + (t * (-9.082695e-03 + pol))
    return 6.1078 / pol**8
"""

from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import numpy as np

from sharptab.constants import *
import sharptab.thermo as thermo
import sharptab.profile as profile
import sharptab.params as params
import sharptab.interp as interp
import sharptab.winds as winds
import sharptab.utils as utils

sigma = 1.5
def sharppy_calcs(**kwargs):
    tmpc = kwargs.get('tmpc')
    dwpc = kwargs.get('dwpc')
    hght = kwargs.get('hght')
    wdir = kwargs.get('wdir')
    wspd = kwargs.get('wspd')
    pres = kwargs.get('pres')

    mucape = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    mlcape = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    mlcin = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    mulpl = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    ebot = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    etop = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    eshr = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    esrh = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    estp = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    ebwd_u = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    ebwd_v = np.zeros((tmpc.shape[1], tmpc.shape[2]))

    for j in range(tmpc.shape[1]):
        for i in range(tmpc.shape[2]):
            prof = profile.create_profile(pres=pres[:,j,i], tmpc=tmpc[:,j,i],
                                          hght=hght[:,j,i], dwpc=dwpc[:,j,i],
                                          wspd=wspd[:,j,i], wdir=wdir[:,j,i])

            # Effective inflow and shear calculations
            eff_inflow = params.effective_inflow_layer(prof)
            ebot[j,i] = interp.to_agl(prof, interp.hght(prof, eff_inflow[0]))
            etop[j,i] = interp.to_agl(prof, interp.hght(prof, eff_inflow[1]))

            # This isn't quite right...need to find midpoint between eff Inflow
            # bottom and EL
            ebwd_u[j,i], ebwd_v[j,i] = winds.wind_shear(prof, pbot=eff_inflow[0],
                                                        ptop=500)
            eshr[j,i] = utils.mag(ebwd_u[j,i], ebwd_v[j,i])

            # Bunkers storm motion function not implemented yet.
            #srwind = params.bunkers_storm_motion(prof)
            #esrh[j,i] = winds.helicity(prof, ebot[j,i], etop[j,i], stu=srwind[0],
            #                           stv = srwind[1])[0]
            esrh[j,i] = winds.helicity(prof, ebot[j,i], etop[j,i])[0]

            # Parcel buoyancy calculations
            mupcl = params.parcelx(prof, flag=3)
            mlpcl = params.parcelx(prof, flag=4)
            mucape[j,i] = mupcl.bplus
            mulpl[j,i] = mupcl.lplhght
            mlcape[j,i] = mlpcl.bplus
            mlcin[j,i] = mlpcl.bminus
            estp[j,i] = params.stp_cin(mlpcl.bplus, esrh[j,i], eshr[j,i],
                                     mlpcl.lclhght, mlpcl.bminus)

    # Apply some data masks
    mlcin = np.where(mlcape < 25., 0, mlcin)
    mulpl = np.where(mlcape < 25., 0, mulpl)
    mucape = gaussian_filter(mucape, sigma=sigma)
    mlcape = gaussian_filter(mlcape, sigma=sigma)
    mlcin = gaussian_filter(mlcin, sigma=sigma)
    mulpl = gaussian_filter(mulpl, sigma=sigma)

    eshr = np.where(np.isnan(eshr), 0, eshr)
    eshr = gaussian_filter(eshr, sigma=sigma)
    eshr = np.where(eshr < 24., np.nan, eshr)
    ebwd_u = np.where(eshr < 24., np.nan, ebwd_u)
    ebwd_v = np.where(eshr < 24., np.nan, ebwd_v)

    esrh = np.where(np.isnan(esrh), 0, esrh)
    esrh = gaussian_filter(esrh, sigma=sigma)
    estp = np.where(np.isnan(estp), 0, estp)
    estp = gaussian_filter(estp, sigma=sigma)

    idx = np.where(np.isnan(ebot))
    ebot = np.where(np.isnan(ebot), 0, ebot)
    ebot = gaussian_filter(ebot, sigma=sigma)
    ebot[idx] = np.nan

    ret = {
        'mucape' : mucape,
        'mlcape' : mlcape,
        'mlcin' : mlcin,
        'mulpl' : mulpl,
        'ebot' : ebot,
        'etop': etop,
        'ebwd_u': ebwd_u,
        'ebwd_v': ebwd_v,
        'eshr': eshr,
        'esrh': esrh,
        'estp': estp
    }
    return ret

def routine_calcs(**kwargs):
    tmpc = kwargs.get('tmpc')
    dwpc = kwargs.get('dwpc')
    hght = kwargs.get('hght')
    uwnd = kwargs.get('uwnd')
    vwnd = kwargs.get('vwnd')
    pres = kwargs.get('pres')
    sp = kwargs.get('sp')
    dx = kwargs.get('dx')
    dy = kwargs.get('dy')

    # Temperature
    t_925 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    t_850 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    t_700 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    t_500 = np.zeros((tmpc.shape[1], tmpc.shape[2]))

    # moisture
    td_925 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    td_850 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    rh = np.zeros((tmpc.shape[1], tmpc.shape[2]))

    # Heights
    z_925 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    z_850 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    z_700 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    z_500 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    z_300 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    z_250 = np.zeros((tmpc.shape[1], tmpc.shape[2]))

    # Kinematics
    u_925 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    v_925 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    u_850 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    v_850 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    u_700 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    v_700 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    u_500 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    v_500 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    u_300 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    v_300 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    u_250 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    v_250 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    wspd_500 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    wspd_300 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    div_300 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    wspd_250 = np.zeros((tmpc.shape[1], tmpc.shape[2]))
    div_250 = np.zeros((tmpc.shape[1], tmpc.shape[2]))

    for j in range(tmpc.shape[1]):
        for i in range(tmpc.shape[2]):
            # 700-500 mb mean RH
            idx = np.where(np.logical_and(pres[:,j,i]>=500, pres[:,j,i]<=700))
            t_mean = np.mean(tmpc[idx[0],j,i])
            td_mean = np.mean(dwpc[idx[0],j,i])
            rh[j,i] = thermo.rh_from_dewpoint(t_mean+ZEROCNK, td_mean+ZEROCNK)

            # Geopotential Heights
            z_925[j,i], t_925[j,i] = calc_hght_and_t(tmpc[:,j,i], pres[:,j,i],
                                                sp[j,i], hght[:,j,i]*G, 925.)
            z_850[j,i], t_850[j,i] = calc_hght_and_t(tmpc[:,j,i], pres[:,j,i],
                                                sp[j,i], hght[:,j,i]*G, 850.)
            z_700[j,i], t_700[j,i] = calc_hght_and_t(tmpc[:,j,i], pres[:,j,i],
                                                sp[j,i], hght[:,j,i]*G, 700.)
            z_500[j,i] = np.interp(500, pres[::-1,j,i], hght[::-1,j,i])
            z_300[j,i] = np.interp(300, pres[::-1,j,i], hght[::-1,j,i])
            z_250[j,i] = np.interp(250, pres[::-1,j,i], hght[::-1,j,i])

            # U- and V- winds
            u_925[j,i] = np.interp(925, pres[::-1,j,i], uwnd[::-1,j,i])
            v_925[j,i] = np.interp(925, pres[::-1,j,i], vwnd[::-1,j,i])
            u_850[j,i] = np.interp(850, pres[::-1,j,i], uwnd[::-1,j,i])
            v_850[j,i] = np.interp(850, pres[::-1,j,i], vwnd[::-1,j,i])
            u_700[j,i] = np.interp(700, pres[::-1,j,i], uwnd[::-1,j,i])
            v_700[j,i] = np.interp(700, pres[::-1,j,i], vwnd[::-1,j,i])
            u_500[j,i] = np.interp(500, pres[::-1,j,i], uwnd[::-1,j,i])
            v_500[j,i] = np.interp(500, pres[::-1,j,i], vwnd[::-1,j,i])
            u_300[j,i] = np.interp(300, pres[::-1,j,i], uwnd[::-1,j,i])
            v_300[j,i] = np.interp(300, pres[::-1,j,i], vwnd[::-1,j,i])
            u_250[j,i] = np.interp(250, pres[::-1,j,i], uwnd[::-1,j,i])
            v_250[j,i] = np.interp(250, pres[::-1,j,i], vwnd[::-1,j,i])

            # Temperatures
            #t_850[j,i] = np.interp(850., pres[::-1,j,i], tmpc[::-1,j,i])
            #t_700[j,i] = np.interp(700., pres[::-1,j,i], tmpc[::-1,j,i])
            t_500[j,i] = np.interp(500., pres[::-1,j,i], tmpc[::-1,j,i])

            # Dewpoints
            td_850[j,i] = np.interp(850., pres[::-1,j,i], dwpc[::-1,j,i])
            td_925[j,i] = np.interp(925., pres[::-1,j,i], dwpc[::-1,j,i])

    # Wind speeds and post-processing
    wspd_500 = np.sqrt(u_500**2 + v_500**2) * MS2KTS
    wspd_300 = np.sqrt(u_300**2 + v_300**2) * MS2KTS
    wspd_250 = np.sqrt(u_250**2 + v_250**2) * MS2KTS
    div_300 = thermo.divergence(u_300, v_300, dx, dy) * 10e4
    div_300 = gaussian_filter(div_300, sigma=sigma)
    div_250 = thermo.divergence(u_250, v_250, dx, dy) * 10e4
    div_250 = gaussian_filter(div_250, sigma=sigma)
    ret = {
        'rh_700-500': rh*100,
        'z_925': z_925,
        'z_850': z_850,
        'z_700': z_700,
        'z_500': z_500,
        'z_300': z_300,
        'z_250': z_250,
        'u_925': u_925,
        'v_925': v_925,
        'u_850': u_850,
        'v_850': v_850,
        'u_700': u_700,
        'v_700': v_700,
        'u_500': u_500,
        'v_500': v_500,
        'u_300': u_300,
        'v_300': v_300,
        'u_250': u_250,
        'v_250': v_250,
        't_925': t_925,
        't_850': t_850,
        't_700': t_700,
        't_500': t_500,
        'td_925': td_925,
        'td_850': td_850,
        'wspd_500': wspd_500,
        'wspd_300': wspd_300,
        'div_300': div_300,
        'wspd_250': wspd_250,
        'div_250': div_250
    }
    return ret

def calc_hght_and_t(t_arr, p_arr, p_s, phi_arr, p):
    """Following equations 5, 6, 15, and 16-19b in _[1], computes the
    geopotential height and temperatures for regions where the geopotential
    is below ground level.

    _[1]: Trenberth, K. E., J. C. Berry, and L. E. Buja, 1993: Vertical
        Interpolation and Truncation of Model-Coordinate Data. NCAR Tech. Note
        NCAR/TN-396+STR, 54 pp, http://dx.doi.org/10.5065/D6HX19NH

    Parameters
    ----------
    t_arr : numpy array
        Array of temperatures on model levels [C]
    p_arr : numpy array
        Array of pressures on model levels [hPa]
    p_s : float
        Surface pressure [hPa]
    phi_s : float
        Surface geopotential [m**2/s**2]
    p : float or int
        Isobaric surface we're interpolating to

    Returns
    -------
    z : float
        Geopotential height at the requested p level [m]
    t : float
        Temperature at the requested p level [C]
    """
    alpha = 0.0065 * R_d / G
    t_arr = t_arr + ZEROCNK
    t_NL = t_arr[0]
    p_NL = p_arr[0]
    t_star = t_NL + (alpha * t_NL * ((p_s/p_NL) - 1))

    inc = R_d * t_star * np.log(p/p_s) * (1 + 0.5 * alpha * np.log(p/p_s) +
          (1./6.) * (alpha * np.log(p/p_s)) ** 2)
    phi_s = phi_arr[0]
    phi = phi_s - inc
    z = phi / G

    to = t_star + (0.0065 * phi_s / G)
    t_pl = min(to, 298)
    if 2000. <= z <= 2500.:
        to_prime = 0.002 * ((2500. - (phi_s/G)) * to + (phi_s/G-2000.)*t_pl)
    elif z > 2500.:
        to_prime = t_pl

    if z > 2000.:
        if to_prime < t_star:
            alpha = 0.
        else:
            alpha = (R_d * (to_prime - t_star)) / phi_s
    T = t_star * (1 + alpha * np.log(p/p_s) + 0.5 * (alpha * np.log(p/p_s))**2 +
                  (1./6.) * (alpha * np.log(p/p_s)) ** 3)

    if p > p_s: # Below ground
        z = phi / G
        t = T
    else:       # Above ground level. Simple vertical interpolation
        z = np.interp(p, p_arr[::-1], phi_arr[::-1]/G)
        t = np.interp(p, p_arr[::-1], t_arr[::-1])
    return z, t-ZEROCNK

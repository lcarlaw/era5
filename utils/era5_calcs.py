"""The ERA5 vertical dimension is defined by an eta coordinate system where
pure pressure coordinates at the model top transition to a terrain-following
sigma coordinate system at the surface. As a result, geopotential and
pressure values need to be computed on each model level based on the surface
pressure and geopotential height.

You don't need all model levels, but this script does require a contiguous set
of levels all the way to the bottom.
"""

import xarray as xr
from numba import njit
import numpy as np

# ------------------------------------------------------------------------------
# Defining some constants
# ------------------------------------------------------------------------------
R_d = 287.06    # dry air,     J K**-1 kg**-1
lvlmax = 137    # ERA5 vertical levels

@njit
def calc_z_and_p(lnsp, tk, q, z, levels):
    """
    This function to calculate geopotential and pressure on ERA5 model levels
    and is based on Peter Ukkonen's Julia module *calc_geopotential_ERA5.jl*
    but re-written to run with numba's jit decorator.

    Parameters
    ----------
    lnsp: numpy array [time, latitude, longitude]
        Logarithm of surface pressure (Pa)
    tk: numpy array [time, level, latitude, longitude]
        Temperature (Kelvin)
    q: numpy array [time, level, latitude, longitude]
        Specific Humidity (kg kg**-1)
    z: numpy array [time, latitude, longitude]
        Geopotential (m**2 s**-2)
    levels: numpy array
        Listing of the model levels.

    Returns
    -------
    hght : numpy array [level]
        Geopotential on model levels (m**2 s**-2)
    pres: numpy array [level]
        Pressure on model levels (Pa)
    """
    nlevs = levels.shape[0]
    sp = np.exp(lnsp)
    z_h = 0.
    hght = np.zeros(nlevs)
    pres = np.zeros(nlevs)

    # No good way to handle the coefficients. Can't pass a pandas array (csv
    # read-in) or define this as a global variable in the jitted function, so
    # have to declare here.
    a_coeff = [
        0.0,2.000365,3.102241,4.666084,6.827977,9.746966,13.605424,18.608931,
        24.985718,32.985710,42.879242,54.955463,69.520576,86.895882,107.415741,
        131.425507,159.279404,191.338562,227.968948,269.539581,316.420746,
	    368.982361,427.592499,492.616028,564.413452,643.339905,729.744141,
        823.967834,926.344910,1037.201172,1156.853638,1285.610352,1423.770142,
        1571.622925,1729.448975,1897.519287,2076.095947,2265.431641,2465.770508,
        2677.348145,2900.391357,3135.119385,3381.743652,3640.468262,3911.490479,
        4194.930664,4490.817383,4799.149414,5119.895020,5452.990723,5798.344727,
        6156.074219,6526.946777,6911.870605,7311.869141,7727.412109,8159.354004,
        8608.525391,9076.400391,9562.682617,10065.978516,10584.631836,11116.662109,
        11660.067383,12211.547852,12766.873047,13324.668945,13881.331055,
        14432.139648,14975.615234,15508.256836,16026.115234,16527.322266,
        17008.789063,17467.613281,17901.621094,18308.433594,18685.718750,
        19031.289063,19343.511719,19620.042969,19859.390625,20059.931641,
        20219.664063,20337.863281,20412.308594,20442.078125,20425.718750,
        20361.816406,20249.511719,20087.085938,19874.025391,19608.572266,
        19290.226563,18917.460938,18489.707031,18006.925781,17471.839844,
	    16888.687500,16262.046875,15596.695313,14898.453125,14173.324219,
        13427.769531,12668.257813,11901.339844,11133.304688,10370.175781,
        9617.515625,8880.453125,8163.375000,7470.343750,6804.421875,6168.531250,
        5564.382813,4993.796875,4457.375000,3955.960938,3489.234375,3057.265625,
        2659.140625,2294.242188,1961.500000,1659.476563,1387.546875,1143.250000,
        926.507813,734.992188,568.062500,424.414063,302.476563,202.484375,
        122.101563,62.781250,22.835938,3.757813,0.0,0.0
    ]

    b_coeff = [
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
        0.0,0.000007,0.000024,0.000059,0.000112,0.000199,0.000340,0.000562,
        0.000890,0.001353,0.001992,0.002857,0.003971,0.005378,0.007133,
        0.009261,0.011806,0.014816,0.018318,0.022355,0.026964,0.032176,
        0.038026,0.044548,0.051773,0.059728,0.068448,0.077958,0.088286,
        0.099462,0.111505,0.124448,0.138313,0.153125,0.168910,0.185689,0.203491,
    	0.222333,0.242244,0.263242,0.285354,0.308598,0.332939,0.358254,0.384363,
        0.411125,0.438391,0.466003,0.493800,0.521619,0.549301,0.576692,0.603648,
        0.630036,0.655736,0.680643,0.704669,0.727739,0.749797,0.770798,0.790717,
        0.809536,0.827256,0.843881,0.859432,0.873929,0.887408,0.899900,0.911448,
        0.922096,0.931881,0.940860,0.949064,0.956550,0.963352,0.969513,
    	0.975078,0.980072,0.984542,0.988500,0.991984,0.995003,0.997630,1.0
    ]

    Ph_levplusone = a_coeff[lvlmax] + (b_coeff[lvlmax] * sp)
    ilevel = nlevs - 1
    for lev in levels[::-1]:
        t_level = tk[ilevel]
        q_level = q[ilevel]

        # Moist temperature
        t_level = t_level * (1. + 0.609133 * q_level)

        # Pressures on half-levels
        levint = int(lev)
        Ph_lev = a_coeff[levint-1] + (b_coeff[levint-1] * sp)
        pres[ilevel] = (Ph_lev + Ph_levplusone) / 2.

        if lev == 1:
            dlogP = np.log(Ph_levplusone / 0.1)
            alpha = np.log(2)
        else:
            dlogP = np.log(Ph_levplusone / Ph_lev)
            dP = Ph_levplusone - Ph_lev
            alpha = 1. - ((Ph_lev / dP) * dlogP)
        TRd = t_level * R_d

        # z_f is the geopotential of this full level
        # Integrate from previous (lower) half-level z_h to the full level
        z_f = z_h + (TRd * alpha)

        # Geopotential (add in surface geopotential)
        hght[ilevel] = z_f + z

        # z_h is the geopotential of 'half-levels'
        # Integrate z_h to next half level
        z_h = z_h + (TRd * dlogP)
        Ph_levplusone = Ph_lev

        # Up we go!
        ilevel -= 1

    return hght, pres

def exec(sfc_file, lvl_file):
    """Execute the calculation of geopotential and pressure on the ERA5 native
    model levels. Uses numba to massively speed up serial calculations.

    Parameters
    ----------
    sfc_file : string
        Path to the ERA5 file containing the surface variables on level "1"
    lvl_file: string
        Path to the ERA5 file containing the model variables

    Returns
    -------
    Outputs a netcdf4 file containing the geopotential and pressure on model
    levels.
    """
    ds_sfc = xr.open_dataset(sfc_file)
    ds_lev = xr.open_dataset(lvl_file, chunks={'time':20})

    q = ds_lev['q'].values
    tk = ds_lev['t'].values
    levels = ds_lev['level'].values
    z = ds_sfc['z'].values
    lnsp = ds_sfc['lnsp'].values

    hght = np.zeros((tk.shape))
    pres = np.zeros((tk.shape))
    ntime, nlevel, nlat, nlon = tk.shape
    for k in range(ntime):
        print(k)
        for j in range(nlat):
            for i in range(nlon):
                hght[k,:,j,i], pres[k,:,j,i] = calc_z_and_p(lnsp[k,j,i],
                                                            tk[k,:,j,i],
                                                            q[k,:,j,i],
                                                            z[k,j,i],
                                                            levels)
    output_data = {'hght':hght, 'pres':pres}
    ds = xr.Dataset(
        {'hght': (('time','level','latitude','longitude'), hght),
         'pres': (('time','level','latitude','longitude'), pres),
         'lnsp': (('time','latitude','longitude'), lnsp)
        },
        coords = {
            'time': ds_lev.time.values,
            'longitude': ds_lev.longitude.values,
            'latitude': ds_lev.latitude.values,
            'level': levels
        }
    )
    ds.to_netcdf('derived.nc')
    return output_data

#exec("./data/era5_surface_reduced.nc", "./data/era5_reduced.nc")

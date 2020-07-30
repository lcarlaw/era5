
##############################################################################*
# Default data
#
# The DEFAULT_DATA dictionary contains base information about field names and
# pressure levels associated with the ERA5 dataset. These specific variables
# are eventually passed to the cdsapi module for downloading off the MARS
# servers. As a result, **DO NOT EDIT** the strings in the DEFAULT_DATA
# dictionary. Full documentation on the various parameters, levels, and fields
# available for download can be found at the following page:
#
# https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation
###############################################################################
pressure_vars = 'pressure_vars'
surface_vars = 'surface_vars'
pressure_levs = 'pressure_levs'

DEFAULT_DATA = {

    pressure_vars: [
                    'geopotential',
                    'vorticity',
                    'temperature',
                    'specific_humidity',
                    'u_component_of_wind',
                    'v_component_of_wind'
    ],

    surface_vars: [
                   '10m_u_component_of_wind',
                   '10m_v_component_of_wind',
                   '2m_temperature',
                   '2m_dewpoint_temperature',
                   'mean_sea_level_pressure',
                   'k_index',
                   'total_totals_index',
                   'surface_pressure',
                   'total_column_water_vapour',
                   'geopotential'
    ],

    pressure_levs : [
                     '1000', '975', '950',
                     '925',  '900', '875',
                     '850',  '825', '800',
                     '775',  '750', '700',
                     '650',  '600', '550',
                     '500',  '450', '400',
                     '350',  '300', '250',
                     '200',  '150', '125',
                     '100', '70', '50', '30'
    ]
}

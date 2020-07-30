import numpy as np
from collections import defaultdict
import argparse
from datetime import datetime, timedelta
import os
import cdsapi

from utils.era5_config import *
from utils.mapinfo import domains

pwd = os.environ['PWD']
def download_era5(start, end, data_path, domain=None):
    c = cdsapi.Client()

    if domain is not None:
        plot_bounds = domains[domain]
    else:
        plot_bounds = domains['US']

    times = defaultdict(list)
    while start <= end:
        time_string = datetime.strftime(start, '%Y%m%d_%H')
        times['time-str'].append(time_string)
        times['years'].append(str(start.year))
        times['months'].append(str(start.month).zfill(2))
        times['days'].append(str(start.day).zfill(2))
        times['hours'].append(str(start.hour).zfill(2) + ':00')
        start += timedelta(hours=1)

    years = sorted(set(times['years']))
    months = sorted(set(times['months']))
    days = sorted(set(times['days']))
    hours = sorted(set(times['hours']))

    outfile = '%s-%s.nc' % (times['time-str'][0], times['time-str'][-1])

    '''
    # Surface Level variables (need for Pressure and Geopotential calculations)
    c.retrieve(
        'reanalysis-era5-complete',
        {
            'class': 'ea',
            'date': '2020-04-01/to/2020-04-30',
            'expver': '1',
            'levelist': '1',
            'levtype': 'ml',
            'param': '129/152',
            'stream': 'oper',
            'time': '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00',
            'type': 'an',
            'area' : '%s/%s/%s/%s' % (plot_bounds[-1]+5,plot_bounds[0]-8,plot_bounds[1]-3,plot_bounds[2]+3),
            'grid': '0.25/0.25',
            'format': 'netcdf'
    },
    '%s/%s' % (data_path, outfile))
    '''


    # Model Level variables
    c.retrieve(
        'reanalysis-era5-complete',
        {
            'class': 'ea',
            'date': '2020-05-01/to/2020-05-31',
            'expver': '1',
            'levelist': '1/57/58/59/60/61/62/63/64/65/66/67/68/69/70/71/72/73/74/75/76/77/78/79/80/81/82/83/84/85/86/87/88/89/90/91/92/93/94/95/96/97/98/99/100/101/102/103/104/105/106/107/108/109/110/111/112/113/114/115/116/117/118/119/120/121/122/123/124/125/126/127/128/129/130/131/132/133/134/135/136/137',
            'levtype': 'ml',
            'param': '129/130/131/132/133/152',
            'stream': 'oper',
            'time': '00:00:00/01:00:00/02:00:00/03:00:00/04:00:00/05:00:00/06:00:00/07:00:00/08:00:00/09:00:00/10:00:00/11:00:00/12:00:00/13:00:00/14:00:00/15:00:00/16:00:00/17:00:00/18:00:00/19:00:00/20:00:00/21:00:00/22:00:00/23:00:00',
            'type': 'an',
            'area' : '%s/%s/%s/%s' % (plot_bounds[-1]+5,plot_bounds[0]-8,plot_bounds[1],plot_bounds[2]),
            'grid': '0.25/0.25',
            'format': 'netcdf'
    },
    '%s/%s' % (data_path, outfile))
    





    """
    if not os.path.exists(data_path + '/' + outfile):
        c.retrieve(
            'reanalysis-era5-pressure-levels',
            {
                'product_type': 'reanalysis',
                'variable': DEFAULT_DATA[pressure_vars],
                'pressure_level': DEFAULT_DATA[pressure_levs],
                'year': years,
                'month': months,
                'day': days,
                'time': hours,
                'format': 'netcdf',
                'area' : [plot_bounds[-1]+5,plot_bounds[0]-8,plot_bounds[1],plot_bounds[2]],
                'grid': [0.25, 0.25],
            },
            '%s/%s' % (data_path, outfile))

    # Surface and other single-level variables
    if not os.path.exists(data_path + '/sfc_' + outfile):
        c.retrieve(
            'reanalysis-era5-single-levels',
            {
                'product_type': 'reanalysis',
                'variable': DEFAULT_DATA[surface_vars],
                'year': times['years'],
                'month': times['months'],
                'day': times['days'],
                'time': times['hours'],
                'format': 'netcdf',
                'area' : [plot_bounds[-1]+1,plot_bounds[0]-1,plot_bounds[1]-1,plot_bounds[2]+1],
                'grid': [0.25, 0.25],
            },
            '%s/sfc_%s' % (data_path, outfile))
    """
    return

def convert_gaussian_to_latlon(infile, outfile):
    """Convert the ERA5 reduced gaussian grid to lat-lon for easier reading in
    xarray
    """
    arg = "cdo -s -f nc4 -t ecmwf setgridtype,regular %s %s.nc" % (infile, outfile)
    os.subprocess(arg, shell=True)

def parse_time(time_str):
    plot_time = datetime.strptime(time_str, '%Y-%m-%d/%H')
    return plot_time

def download(start_time, end_time, domain):
    data_path = pwd + '/data/'
    start = datetime.strptime(start_time, '%Y-%m-%d/%H')
    end = datetime.strptime(end_time, '%Y-%m-%d/%H')
    _file = download_era5(start, end, data_path, domain=None)
    #_radar_file = download_radar(plot_time, data_path)
    return

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-s', '--start-time', dest='start_time', help="YYYY-MM-DD/HH The desired valid time of the archived RAP data to view.")
    ap.add_argument('-e', '--end-time', dest='end_time', help="YYYY-MM-DD/HH The desired valid time of the archived RAP data to view.")
    ap.add_argument('-d', '--domain', dest='domain', help="[MW|SGP|CGP|NGP|GL|SE|MA|NE|GL|NW|GBSN|SW|CONUS] Plotting domain string. Set to None for no plot.")
    args = ap.parse_args()
    np.seterr(all='ignore')

    download(args.start_time,
             args.end_time,
             args.domain,
            )

if __name__ == "__main__":
    main()

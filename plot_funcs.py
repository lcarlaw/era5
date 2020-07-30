"""plot_funcs.py

    Primary function controlling the actual plot creation. Plot parameters are
    specified in the plotparms.py file within this top-level directory.

"""

import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime
import pylab
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.ndimage import gaussian_filter
from scipy import spatial

import metpy.calc as mpcalc
from plotparms import *
from sharptab.constants import G, ZEROCNK, MS2KTS
import utils.data_calcs as data_calcs
import sharptab.thermo as thermo
from utils.mapinfo import *
import utils.era5_calcs as era5_calcs
from airports import metadata

def make_basemap(bounds):
    """Create a basemap object so we don't have to re-create this over again.

    Parameters
    ----------
    bounds : list
        Defined in the ./utils/mapinfo.py file. Format for the list convention:
        [Left, Bottom, Right, Top]

    Returns
    -------
    m_ : Basemap plotting object
    """
    m_ = Basemap(projection='stere', llcrnrlon=bounds[0], llcrnrlat=bounds[1],
                 urcrnrlon=bounds[2], urcrnrlat=bounds[3], lat_ts=25,
                 lat_0=25, lon_0=-95, resolution='i')
    m_.drawcoastlines(color='#a09ea0', linewidth=2)
    m_.drawstates(color='#a09ea0', linewidth=2)
    m_.drawcountries(color='#a09ea0', linewidth=2)
    m_.drawcounties(color="#e6e6e6", linewidth=1)
    return m_

def process_files(filename, domain, realtime=False):
    """Read in the presure-level and surface-level data files, merge them
    based on the surface pressure values, and convert units.

    Parameters
    ----------
    filename : string
        Either a full path to the main pressure-level fields file, or path to
        a directory containing real-time downloaded grib data. If the latter,
        realtime must = True
    domain : string
        Domain plot area, specified in the command line arguments at main.py

    Returns
    -------
    ret : dictionary
        Dictionary containing the merged datasets and additional variables
        necessary for sharppy and plotting routines.
    """
    bounds = domains[domain]

    # We don't need the full dataset. Figure out what data lies in our specified
    # domain and subset it. Give a litte buffer for the map projection.
    # Latitudes are reversed from ERA in NCEP (realtime) data.
    if not realtime:
        data = xr.open_dataset(filename)

        # Is this data on pressure levels or full-resolution native coordinates?
        # Data entries and preprocessing will change for each.
        dtype = data.level.values.dtype
        if dtype == 'float64':
            data_format = 'native'
            data['longitude'] = data.longitude.values - 360.
            lons = data['longitude']
        elif dtype == 'int32':
            data_format = 'isobaric'
            lons = data.longitude.values

        lats = data.latitude.values
        ds = data.sel(latitude=slice(bounds[3]+5, bounds[1]),
                      longitude=slice(bounds[0]-8, bounds[2]))
    else:
        data = xr.open_dataset(filename, engine='cfgrib',
                               backend_kwargs={'filter_by_keys':
                                              {'typeOfLevel': 'hybrid'}})
        data = data.assign_coords(time=data.valid_time.values)
        lons = data.longitude.values - 360.
        lats = data.latitude.values
        idx_lon = np.where(np.logical_and(lons>=bounds[0]-1, lons<=bounds[2]+1),
                           1, np.nan)
        idx_lat = np.where(np.logical_and(lats>=bounds[1]-1, lats<=bounds[3]+1),
                           1, np.nan)
        idx = np.where(np.logical_and(idx_lon==1, idx_lat==1))
        ds = data.sel(x=slice(idx[1].min(), idx[1].max()),
                      y=slice(idx[0].min(), idx[0].max()))

    print("Reading variables")
    if realtime:
        pres = ds.pres.values / 100.
        hght = ds.gh.values
    else:
        if data_format == 'native':
            # We have to compute the Z and P values for raw ERA5 data
            era5_calcs.exec("./data/era5_surface_reduced.nc",
                                  "./data/era5_reduced.nc")
            #pres = derived['pres'] / 100.
            #hght = derived['hght'] / 9.80665
            ds_tmp = xr.open_dataset('/Users/leecarlaw/scripts/era5-to-spc/derived.nc')
            ds_tmp['longitude'] = ds_tmp.longitude.values - 360.
            data_tmp = ds_tmp.sel(latitude=slice(bounds[3]+5, bounds[1]),
                                  longitude=slice(bounds[0]-8, bounds[2]))
            pres = data_tmp['pres'].values / 100.
            hght = data_tmp['hght'].values / G
            sp = np.exp(data_tmp['lnsp'].values) / 100.
        elif data_format == 'isobaric':
            hght = ds.z.values / G
            tmp = np.tile(ds.level.values,
                          (hght.shape[0],hght.shape[2],hght.shape[3]))
            pres = np.reshape(tmp.ravel(), hght.shape, order='F')

    uwnd = ds.u.values
    vwnd = ds.v.values
    tmp = ds.t.values

    spfh = ds.q.values
    dwp = thermo.dewpoint_from_specific_humidity(spfh, tmp, pres)

    print("outputting dictionary")
    lats = ds.latitude.values
    lons = ds.longitude.values
    lons, lats = np.meshgrid(lons, lats)

    # Define the new arrays and convert units
    tmpc = tmp - ZEROCNK
    dwpc = dwp
    wdir = thermo.wind_direction(uwnd, vwnd)
    wspd = thermo.wind_speed(uwnd, vwnd) * MS2KTS

    if realtime: sp = -99.

    # For our jitted sharppy routines, need to be REALLY careful with units or
    # things will break. If this isn't realtime data, we'll need to reverse
    # the vertical orientation since native ERA5 data is oriented top to
    # bottom.
    order = 1
    if not realtime and data_format == 'native': order = -1
    tmpc = np.array(tmpc[:,::order], dtype='float64')
    dwpc = np.array(dwpc[:,::order], dtype='float64')
    hght = np.array(hght[:,::order], dtype='float64')
    wdir = np.array(wdir[:,::order], dtype='float64')
    wspd = np.array(wspd[:,::order], dtype='float64')
    pres = np.array(pres[:,::order], dtype='float64')
    uwnd = uwnd[:,::order]
    vwnd = vwnd[:,::order]

    ret = {
        'tmpc': tmpc,
        'dwpc': dwpc,
        'hght': hght,
        'wdir': wdir,
        'wspd': wspd,
        'u': uwnd,
        'v': vwnd,
        'pres': pres,
        'lons': lons,
        'lats': lats,
        'sp': sp,
        'time': ds.time.values
    }
    return ret

def make_plots(filename, domain, realtime):
    """Primary driving function to render output plots.
    """
    bounds = domains[domain]
    data = process_files(filename, domain, realtime)

    fig_aspect = 1.5
    fig_wid = 15
    fig_hght = fig_wid / fig_aspect

    fig = pylab.figure(figsize=(fig_wid, fig_hght), dpi=100)
    ax_left = 0.01
    ax_bot = 0.0475
    ax_hght = 0.935
    ax_wid = ax_hght / fig_aspect
    ax = pylab.axes([ax_left, ax_bot, 0.975, ax_hght], xticks=[], yticks=[])

    # Create the re-usable basemap plotting object
    m = make_basemap(bounds)
    x, y = m(data['lons'], data['lats'])
    dx, dy = thermo.lat_lon_grid_deltas(data['lons'], data['lats'])

    for t in range(data['tmpc'].shape[0]):
        dt = pd.to_datetime(data['time'][t])
        valid_date = datetime.strftime(dt, '%Y%m%d/%H00')
        save_date = datetime.strftime(dt, '%Y%m%d_%H')

        prof_data = {'pres':data['pres'][t], 'tmpc':data['tmpc'][t],
                     'dwpc':data['dwpc'][t], 'hght':data['hght'][t],
                     'wdir':data['wdir'][t], 'wspd':data['wspd'][t],
                     'uwnd':data['u'][t], 'vwnd':data['v'][t],
                     'dx':dx, 'dy':dy, 'sp':data['sp'][t]}

        # Routine plots and routines
        routine_data = data_calcs.routine_calcs(**prof_data)

        # SHARPpy plots and routines. Includes thermodynamics and any wind shear
        # calculations.
        arrs = data_calcs.sharppy_calcs(**prof_data)

        # Merge the datasets
        arrs = {**arrs, **data, **routine_data}
        #arrs = {**data, **routine_data}
        for parm in PLOT_DICT.keys():
            # Reset the various plotting objects to start the new loop
            cf_data = None; cf = None; cb = None
            c1_data = None; c1 = None; clab1 = None
            c2_data = None; c2 = None; clab2 = None
            c3_data = None; c3 = None; clab3 = None
            barb = None
            meta = PLOT_DICT[parm]

            if meta['cf_data']:
                #cf_data = gaussian_filter(arrs[meta['cf_data']], sigma=sigma)
                cf_data = arrs[meta['cf_data']]
                cf = m.contourf(x, y, cf_data, meta['plot_levs'][0],
                                **PLOT_KWARGS[parm][0])

            if meta['c1_data']:
                c1_data = arrs[meta['c1_data']]
                c1 = m.contour(x, y, c1_data, meta['plot_levs'][1],
                               **PLOT_KWARGS[parm][1])
                clab1 = pylab.clabel(c1, fmt="%d", inline_spacing=8)

            if meta['c2_data']:
                c2_data = arrs[meta['c2_data']]
                c2 = m.contour(x, y, c2_data, meta['plot_levs'][2],
                               **PLOT_KWARGS[parm][2])
                clab2 = pylab.clabel(c2, fmt="%d", inline_spacing=8)

            if meta['c3_data']:
                c3_data = arrs[meta['c3_data']]
                c3 = m.contour(x, y, c3_data, meta['plot_levs'][3],
                               **PLOT_KWARGS[parm][3])
                clab3 = pylab.clabel(c3, fmt="%d", inline_spacing=8)

            # Wind barbs
            if meta['plot_barbs'][0]:
                u_name = meta['plot_barbs'][1][0]
                v_name = meta['plot_barbs'][1][1]
                if len(meta['plot_barbs']) == 2:
                    u = arrs[u_name]
                    v = arrs[v_name]
                elif len(meta['plot_barbs']) == 3:
                    level = meta['plot_barbs'][2]
                    u = arrs[u_name][t,level]
                    v = arrs[v_name][t,level]
                if parm in ['500mb', '700mb', '300mb', '250mb']:
                    barb_color='k'
                else:
                    barb_color = '#c99546'
                if domain == 'US':
                    skip = 7
                else:
                    skip = 4
                barb = m.barbs(x[::skip,::skip], y[::skip,::skip], u[::skip,::skip],
                               v[::skip,::skip], length=6.25, linewidth=0.85,
                               zorder=99, sizes=dict(width=0.15,emptybarb=0.15,
                               spacing=0.1375, height=0.3), flagcolor=None,
                               color=barb_color)

            # Colorbar
            if meta['cf_data']:
                cax = inset_axes(ax, width="18%", height="2.5%", loc="lower left",
                                 bbox_to_anchor=(-0.002,-0.055,1,1),
                                 bbox_transform=ax.transAxes)
                cb = pylab.colorbar(cf, cax=cax, orientation='horizontal',
                                    extend='max')
                cax.xaxis.set_ticks_position('top')
                cax.xaxis.set_label_position('top')
                cax.xaxis.set_tick_params(pad=0.001)
                cax.tick_params(labelsize=10)

            # Image text and information headers
            img_info = "%s %s" % (valid_date, meta['plot_info'])
            t1 = ax.annotate(img_info, xy=(0.5, -0.031), va='top', ha='center',
                             xycoords='axes fraction', fontsize=11)
            t2 = ax.annotate('ERA5 0.25 x 0.25 deg Reanalysis-SPC Meld Project',
                             xy=(0,1), va='bottom', xycoords='axes fraction',
                             fontsize=12, color='b')
            save_name = "%s/%s_%s_%s.png" % (PLOT_DIR, domain, parm, save_date)
            #pylab.savefig(save_name, dpi=pylab.gcf().dpi, transparent=True)
            pylab.savefig(save_name, dpi=pylab.gcf().dpi)
            clean_objects(cf, cb, c1, c2, c3, clab1, clab2, clab3, t2, t1, barb)
    pylab.close()
    return data

def clean_objects(cf, cb, c1, c2, c3, clab1, clab2, clab3, t2, t1, barb):
    """Remove previous plotting objects
    """
    if cf:
        for member in cf.collections:
            member.remove()
    if cb: cb.remove()
    if c1:
        for member in c1.collections:
            member.remove()
    if c2:
        for member in c2.collections:
            member.remove()
    if c3:
        for member in c3.collections:
            member.remove()
    if barb:
        barb[0].remove()
        barb[1].remove()
    if clab1:
        for label in clab1:
            label.remove()
    if clab2:
        for label in clab2:
            label.remove()
    if clab3:
        for label in clab3:
            label.remove()
    t1.remove()
    t2.remove()

'''
.. caution::
    For some reason, merely having this function uncommented results in nothing
    plotting after the 1st plot. Not sure at all why that would be the case...
def clean_objects(*args):
    """Remove previous plotting objects. This is so we don't have to instantiate
    a whole new plot figure and basemap object
    """
    for arg in args:
        attrs = dir(arg)
        # QuadContourSet
        if 'collections' in attrs:
            for item in arg.collections: item.remove()
        # Wind barbs and contour labels
        elif '__len__' in attrs:
            for item in arg:
                item.remove()
        # Text strings and everything else
        elif 'get_text' in attrs: arg.remove()
        else: pass
'''

def nearest_idx(points, lon, lat):
    """Search for the nearest grid point using a KDTree

    Parameters
    ----------
    points : list
        List of lists indicating lon/lat pais: [[LON1, LAT1], [LON2, LAT2]]
    lat : np.array
        2-d array of gridded latitudes
    lon : np.arrat
        2-d array of gridded longitudes
    """
    #if tree is None:
    lonlat = np.column_stack((lon.ravel(), lat.ravel()))
    tree = spatial.cKDTree(lonlat)
    dist, idx = tree.query(points, k=1)
    ind = np.column_stack(np.unravel_index(idx, lon.shape))
    return [(j,i) for j,i in ind]

def write_to_sounding(data, idx_locs, ids, date):
    """Writes data to sounding files.

    Parameters
    ----------
    data : dictionary
        Python dictionary (likely created by function process_files) containing
        the necessary meteorological data
    idx_locs : list of tuples
        Each tuple contains the (j,i) reference points to each ID location in
        the specified lat /lon grid
    ids : list
        List of airport IDs. See the airports.py file for acceptable ID names:
        [ORD, MDW, etc]

    Returns
    -------
    Output sharppy-readable text file. Output is into the SOUNDING_DIR which is
    specified in the plotparms.py file.
    """
    knt = 0
    for idx in idx_locs:
        levs = data['pres'][0,:,idx[0],idx[1]]
        t_out = data['tmpc'][0,:,idx[0],idx[1]]
        td_out = data['dwpc'][0,:,idx[0],idx[1]]
        wdir_out = data['wdir'][0,:,idx[0],idx[1]]
        wspd_out = data['wspd'][0,:,idx[0],idx[1]]
        hght_out = data['hght'][0,:,idx[0],idx[1]]

        out_time = "%s%s%s/%s00"%(date[2:4], date[4:6], date[6:8], date[8:10])
        out_file = "%s/%s.%s" % (SOUNDING_DIR, date, ids[knt])
        f = open(out_file, 'w')

        f.write("%TITLE%\n")
        f.write(" %s   %s" % (ids[knt], out_time))
        f.write("\n\n")
        f.write('   LEVEL     HGHT     TEMP     DWPT     WDIR     WSPD\n')
        f.write('------------------------------------------------------\n')
        f.write('%RAW%\n')

        #start_, end_, inc_, pres_incr = 0, t_out.shape[0], 1, 1
        for row in range(t_out.shape[0]):
            if hght_out[row] > 0:
                out_line = "%s,%s,%s,%s,%s,%s" % (levs[row],
                                                  hght_out[row],
                                                  t_out[row],
                                                  td_out[row],
                                                  wdir_out[row],
                                                  wspd_out[row])
                f.write(out_line +'\n')
        f.write('%END%')
    knt += 1

def get_soundings(data, ids):
    """Main entry function which controls the creation of SHARPpy-readable
    sounding data files by gathering the necessary data for writing.

    Parameters
    ----------
    ids : list
        List of airport IDs. See the airports.py file for acceptable ID names:
        [ORD, MDW, ..., etc.]
    """
    lats = data['lats']
    lons = data['lons']
    points = []
    elevs = []
    for id_ in ids:
        entry = metadata[id_]
        points.append([entry[1], entry[0]])
    idx_locs = nearest_idx(points, lons, lats)

    try:
        date = pd.to_datetime(data['time'][0])
        date = datetime.strftime(date, "%Y%m%d%H")
    except:
        date = data['time'][0]
    write_to_sounding(data, idx_locs, ids, date)

def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))

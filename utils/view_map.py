'''Simple script to plot a domain specified in the mapinfo.py file'''
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mapinfo import domains

import argparse
import xarray as xr
import numpy as np
import numpy.ma as ma

skip = 3

def _make_basemap(bounds):
    m_ = Basemap(projection='stere', llcrnrlon=bounds[0], llcrnrlat=bounds[1],
                 urcrnrlon=bounds[2], urcrnrlat=bounds[3], lat_ts=25,
                 lat_0=25, lon_0=-95, resolution='l')
    m_.drawcoastlines(color='#a09ea0', linewidth=2)
    m_.drawstates(color='#a09ea0', linewidth=2)
    m_.drawcountries(color='#a09ea0', linewidth=2)
    #m_.drawcounties(color="#e6e6e6", linewidth=1)
    return m_

def _clean_objects(*args):
    """Remove previous plotting objects"""

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


def plotmap(domain):
    bounds = domains[domain]

    # Plotting figure
    fig_aspect = 1.5
    fig_wid = 15
    fig_hght = fig_wid / fig_aspect

    fig = pylab.figure(figsize=(fig_wid, fig_hght), dpi=100)
    axes_left = 0.01
    axes_bot = 0.045
    axes_hght = 0.935
    axes_wid = axes_hght / fig_aspect
    ax = pylab.axes([axes_left, axes_bot, 0.975, axes_hght], xticks=[], yticks=[])
    m = _make_basemap(bounds)

    '''
    fname = '/Users/leecarlaw/scripts/era5-to-spc/data/20200701_00-20200701_02.nc'
    ds = xr.open_dataset(fname)
    tmpc = ds.t.values[0,0]
    u = ds.t.values[0,0]
    v = ds.t.values[0,0]
    lons = ds.longitude.values
    lats = ds.latitude.values
    lons, lats = np.meshgrid(lons, lats)
    x, y = m(lons, lats)

    #cax = fig.add_axes([0.01, axes_bot-0.05, 0.15, 0.025])
    cax = inset_axes(ax, width="25%", height="2.5%", loc="lower left",
                     bbox_to_anchor=(-0.002,-0.058,1,1), bbox_transform=ax.transAxes)
    cb = pylab.colorbar(cf, cax=cax, orientation='horizontal', extend='max')
    cax.xaxis.set_ticks_position('top')
    cax.xaxis.set_label_position('top')
    cax.xaxis.set_tick_params(pad=0.001)
    cax.tick_params(labelsize=9)
    t1 = ax.annotate('BOTTOM', xy=(0.5,-0.04), va='top', xycoords='axes fraction')
    t2 = ax.annotate('TEST ERA', xy=(0., 1), va='bottom', xycoords='axes fraction')
    pylab.savefig('test1.png')
    _clean_objects(c1, t1, t, t2, barb, cf, cb)
    pylab.savefig('test2.png')
    '''
    plt.show()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(dest='domain', help='Domain key from mapinfo config file')
    args = ap.parse_args()
    plotmap(args.domain)

if __name__ == "__main__": main()

import argparse
from plotparms import *
from plot_funcs import make_plots, get_soundings
#from utils.mapinfo import *

def plotter(filename, domain=None, realtime=False, ids=None):
    """
    """
    if not domain: domain='US'
    data = make_plots(filename, domain, realtime)

    if ids:
        try:
            ids = ids.strip().split(',')
            get_soundings(data, ids)
        except:
            raise ValueError("Improperly formatted sounding locations: ORD,MDW,DPA")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(dest='filename', help="Full path to netCDF4 file")
    ap.add_argument('-r', '--realtime', dest='realtime', help="default=False. If True, will attempt to download and parse realtime files")
    ap.add_argument('-d', '--domain', dest='domain', help="[MW|SGP|CGP|NGP|GL|SE|MA|NE|GL|NW|GBSN|SW|CONUS] Plotting domain string. Set to None for no plot.")
    ap.add_argument('-s', '--ids', dest='ids', help="Plot soundings: ORD,VPZ,DFW")
    args = ap.parse_args()

    plotter(args.filename,
            domain=args.domain,
            realtime=args.realtime,
            ids=args.ids)

if __name__ == "__main__": main()

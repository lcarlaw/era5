# era5-to-spc
Scripts in this repo attempt to re-create SPC's [Mesoanalysis Archive](https://www.spc.noaa.gov/exper/mesoanalysis/) using the ECMWF's [ERA5 Dataset](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5). Initially, this is aimed at producing static images, but an extension to output various SPC-derived parameters into netCDF files could be tested, as well as apply the code to real time NWP output.

The main overhaul here was to massively accelerate the CPU-intensive thermodynamic calculations performed by SHARPpy (mainly due to constant lifting of parcels) using the [Python Numba](http://numba.pydata.org/) module. This is a non-trivial task, as several standard Python modules and code are not supported by Numba. In addition, the nature of the "Just-in-time" compilation requires explicit `type` declarations within Python `Classes`. As a result, the original SHARPpy code had to be parsed out line-by-line to allow it to work with Numba and the jit decorator, and some flexibility has certainly been lost here. The biggest issues were the lack of `**kwarg` support and numpy masked arrays. In this current iteration of code, it's assumed that the input meteorological arrays are full without any missing/masked data.

## Code Execution Time Improvements
The initial "Just-in-time compiling" of the associated SHARPpy code (creating the Profile and Parcel classes) takes about 30 seconds, and this is the main overhead here. Each subsequent time slice [for a full CONUS domain] takes about 10-20 seconds, and about 4 seconds for a smaller "Midwest" domain test. Compare this with the straight-out-of-the-box SHARPpy code (running in serial) where each time slice takes upwards of 25 to 30 minutes!

* Simple `tqdm` output for a time slice (~10 seconds):
![](https://github.com/lcarlaw/era5-to-spc/raw/master/readme_images/numba.png?raw=true)

* Same for straight serial out of the box run (~30 minutes):
![](https://github.com/lcarlaw/era5-to-spc/raw/master/readme_images/serial.png?raw=true)

## ERA5 Dataset Information
### From the ERA5 online documentation:

ERA5 was produced using 4D-Var data assimilation in CY41R2 of ECMWF’s Integrated Forecast System (IFS), with 137 hybrid sigma/pressure (model) levels in the vertical, with the top level at 0.01 hPa. Atmospheric data are available on these levels and they are also interpolated to 37 pressure, 16 potential temperature and 1 potential vorticity level(s). "Surface or single level" data are also available, containing 2D parameters such as precipitation, 2m temperature, top of atmosphere radiation and vertical integrals over the entire atmosphere. The IFS is coupled to a soil model, the parameters of which are also designated as surface parameters, and an ocean wave model.

The ERA5 dataset contains one (hourly, 31 km) high resolution realisation (referred to as "reanalysis" or "HRES") and a reduced resolution ten member ensemble (referred to as "ensemble" or "EDA"). Generally, the data are available at a sub-daily and monthly frequency and consist of analyses and short (18 hour) forecasts, initialised twice daily from analyses at 06 and 18 UTC. Most analysed parameters are also available from the forecasts. There are a number of forecast parameters, e.g. mean rates and accumulations, that are not available from the analyses.

### Available Parameters
A full listing of the downloadable parameters can be [found here](https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation). Tables 2 and 9 provide information on the most likely parameters you'll want to download.

### A note about derived parameters
The ERA5 vertical dimension is defined by an eta coordinate system where pure pressure coordinates at the model top transition to a terrain-following sigma coordinate system at the surface. As a result, geopotential and pressure values need to be computed on each model level based on the surface pressure and geopotential height.

While we could instead utilize the ERA5 dataset on isobaric surfaces, the benefit of using the native vertical coordinate dataset is the much higher vertical resolution (especially near the ground) which aids in accurate thermodynamic calculations. Another issue to contend with is this requires deriving quantities that fall below the "ground surface" (pressure levels, temperature, dewpoint, etc.). An old, but useful, reference for this derivation can be found in this [NCAR Technical Note](http://dx.doi.org/10.5065/D6HX19NH) from 1993 by Kevin Trenberth, Jeffery Berry, and Lawrence Buja [[1]](#1).

## Basic Setup Notes
You must first create an account on the [CDS Site](https://cds.climate.copernicus.eu/#!/home). Once you register, you'll be directed to copy off an API key and a url string that you'll save into a file: `$HOME/.cdsapirc`.

### Creating the base environment
The setup here proceeds using Anaconda, as well as assuming a completely vanilla Python3 install.  I've also edited my `~/.condarc` file to add conda-forge to the default channels. The order is important here for versioning control as newer versions of matplotlib and basemap cause issues with some of the map features.

```
conda create --name era5 python=3.7
conda activate era5

conda install matplotlib=2.2.4
conda install basemap=1.2.0 basemap-data-hires=1.2.0
conda install xarray netcdf4
conda install cdsapi
conda install numba
conda install metpy
conda install cdo zarr: (optional: only needed for realtime RAP applications)
```

#### Install Julia:
Download the Julia installer binary from https://julialang.org/downloads/. Add the `julia` binary install location to your `PATH` variable. For my install this was:

```
vi ~/.bash_profile
export PATH="/Applications/Julia-1.4.app/Contents/Resources/julia/bin/:$PATH"
source ~/.bash_profile

python3 -m pip install julia
```

Next, install the required packages by opening a julia REPL:

```
julia
>>>import Pkg
>>>Pkg.add("NCDatasets")
```

#### Exporting .yml information:

Export the conda `.yml` file with:

```
conda env export --from-history -f environment.yml
```

## Useage
General usage information.

### Downloading ERA5 Data
```
python get_era5.py -s START_TIME -e END_TIME [-d DOMAIN]
```
`START_TIME`: Initial (earliest) slice of data request. Form is `YYYY-MM-DD/HH`
`END_TIME`: Latest (last) slice of data request. Form is `YYYY-MM-DD/HH`
`DOMAIN`: Optional argument to download a specific domain subset. If left off, this will default to the CONUS region. The acceptable values are located in the `./utils/mapinfo.py` file. There is a helper script called `./utils/view_map.py` to visualize the domains that will be downloaded and plotted.

NetCDF4 files will be saved into a folder called `data` in your present working directory.

Note that depending on the size of your file request (length of time, domain size), download/queuing times can become quite long, potentially on the order of days. You can see the live [queue here](https://cds.climate.copernicus.eu/live/queue), and also actively monitor [your requests here](https://cds.climate.copernicus.eu/cdsapp#!/yourrequests).

### Image Plotting
Plotting is controlled at the top level by `main.py` which accepts several command-line arguments and passes them to the various plotting functions. The general use is:
```
python main.py FILENAME [-d DOMAIN]
```
`FILENAME`: Required argument specifying the full path to the ERA5 netCDF file.
`DOMAIN`: Optional argument specifying the domain to be plotted. These are defined in the `./utils/mapinfo.py` file and can be edited.

The `plotparms.py` file controls most of the plotting configuration, from linestyles and colors, to wind barbs and text labels. It's a goofy file since it has to control so much, and there are undoubtedly much better ways to handle this...

## References
<a id="1">[1]</a>
Trenberth, K. E., J. C. Berry, and L. E. Buja, 1993: Vertical Interpolation and Truncation of Model-Coordinate Data. NCAR Tech. Note NCAR/TN-396+STR, 54 pp, http://dx.doi.org/10.5065/D6HX19NH

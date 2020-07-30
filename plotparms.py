import numpy as np
import os
from utils.era5_config import DEFAULT_DATA
pwd = os.environ['PWD']
PLOT_DIR = pwd + '/images'
SOUNDING_DIR = pwd + '/soundings'

alpha=0.85
sigma = 1.5

pressure_levels = DEFAULT_DATA['pressure_levs']

levels = {
          'mlcape': [100,250,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,
                     5500,6000,6500,7000,7500,8000,9000,10000],
          'mucape': [100,250,500,1000,2000,3000,4000,5000,6000,7000,8000,10000],
          'vort': [10,12,16,20,24,28,32,36,40]
}

colors = {
         'mlcape': ['#ec7c7a','#ec7c7a','#d84138','#d84138','#d84138','#d84138',
                    '#bc271c','#bc271c','#bc271c','#bc271c','#801811','#801811',
                    '#801811','#801811','#801811','#801811','#801811','#801811'],

         'mucape': ['#ec7c7a','#ec7c7a','#d84138','#d84138','#bc271c','#bc271c',
                    '#bc271c','#801811','#801811','#801811','#801811','#801811'],

         'vort': ['#397e23','#5ac039','#9ef94e','#edec4f','#f7d549', '#ed8433',
                  '#e63524', '#e63524']
}
################################################################################
#
# PLOTTING CONFIGURATIONS
#
#   Parameter reference name. Will also controll the save name of the plot.
#      |
#      |
#   'mucape' : {'cf_data': Dictionary key of contour filled parameter
#               'c1_data': Dictionary key of contoured parameter
#               'c2_data': Dictionary key of 2nd contoured parameter
#               'plot_barbs': [True|False, (u name, v name), level]
#               'plot_info': Plot information text string at bottom of image
#               'plot_levs':[[500,1000,2000,4000],
#                            np.arange(250, 2250, 250), # Lifted Parcel Height (m; contour)
#                            levels['cape']]
#    },
################################################################################
# Controls the plotting routines requiring the use of SHARPpy and its lifting
# routines.
SHARPPY_DICT = {
    'mucape' : {'cf_data': 'mulpl', # Lifted Parcel Height (m)
                'c1_data': 'mulpl', # Lifted Parcel Height (m; contour)
                'c2_data': 'mucape', # MUCAPE
                'c3_data': None,
                'plot_barbs': [False, ()],
                'plot_info': 'MUCAPE (J/kg) and lifted parcel level (m AGL, fill)',
                'plot_levs':[[500,1000,2000,4000],
                             np.arange(250, 2250, 250), # Lifted Parcel Height (m; contour)
                             levels['mucape']]
    },

    'mlcape' : {'cf_data': 'mlcin', # MLCIN
               'c1_data': 'mlcin', # MLCIN (contour)
               'c2_data': 'mlcape',
               'c3_data': None,
               'plot_barbs': [True, ('u', 'v'), 0],
               'plot_info': 'MLCAPE (contour) and MLCIN (J/kg, shaded)',
               'plot_levs': [[-500, -100, -25],
                             [-1000,-500,-400,-300,-250,-200,-150,-100,-50,-25],
                             levels['mlcape']]
    },

    'eshr' : {'cf_data': None, # No contour fill
               'c1_data': 'eshr', # EBWD
               'c2_data': None,
               'c3_data': None,
               'plot_barbs': [True, ('ebwd_u', 'ebwd_v')],
               'plot_info': 'Effective bulk shear (kt)',
               'plot_levs': [None,
                            [25,30,35,40,45,50,55,60,65,70,75,80],
                             None]
    },

    'effh' : {'cf_data': 'ebot', # Effective inflow base (m)
               'c1_data': 'esrh', # Effective SRH
               'c2_data': None,
               'c3_data': None,
               #'plot_barbs': [True, ('storm_u', 'storm_v')],
               'plot_barbs': [False, ()],
               'plot_info': 'Eff. Inflow Base (fill, m AGL), ESHR (m2/s2) and storm motion (kt)',
               'plot_levs': [[10, 250, 1000, 3000],
                             [25,50,100,200,300,400,500,600,700,800,900,1000],
                              None]
    },

    'stpc' : {'cf_data': 'mlcin', # Effective inflow base (m)
               'c1_data': 'estp', # Effective SRH
               'c2_data': None,
               'c3_data': None,
               'plot_barbs': [False, ()],
               'plot_info': 'Significant Tornado Parameter (eff layer) and MLCIN (J/kg, shaded at 25 and 100)',
               'plot_levs': [[-500, -100, -25],
                             [.25,.5,1,2,3,4,5,6,7,8,9,10],
                              None]
    }
}

ROUTINE_DICT = {
    '925mb': {'cf_data': 'td_925',
              'c1_data': 't_925',
              'c2_data': 't_925',
              'c3_data': 'z_925',
              'plot_barbs': [True, ('u_925', 'v_925')],
              'plot_info': '925mb height (m MSL, black), temp (C, red), dewpoint, and wind (kt)',
              'plot_levs': [[8,16,18,20,22,24,26,28],
                            np.arange(2,90,2),
                            np.arange(-60,2,2),
                            np.arange(0, 5000, 30)]
    },

    '850mb': {'cf_data': 'td_850',
              'c1_data': 't_850',
              'c2_data': 't_850',
              'c3_data': 'z_850',
              'plot_barbs': [True, ('u_850', 'v_850')],
              'plot_info': '850mb height (m MSL, black), temp (C, red), dewpoint, and wind (kt)',
              'plot_levs': [[6,12,14,16,18,20,22,24],
                            np.arange(2,90,2),
                            np.arange(-60,2,2),
                            np.arange(900, 5000, 30)]
    },

    '700mb': {'cf_data': 'rh_700-500',
              'c1_data': 't_700',
              'c2_data': 't_700',
              'c3_data': 'z_700',
              'plot_barbs': [True, ('u_700', 'v_700')],
              'plot_info': '700mb height (m, MSL), wind (kt), temp (C, red), and 700-500 mb mean RH > 70',
              'plot_levs': [[70,105],
                            np.arange(2,32,2),
                            np.arange(-50,2,2),
                            np.arange(1110, 5000, 30)]
    },

    '500mb': {'cf_data': 'wspd_500',
              'c1_data': 't_500',
              'c2_data': 'z_500',
              'c3_data': None,
              'plot_barbs': [True, ('u_500', 'v_500')],
              'plot_info': '700mb height (m, MSL), wind (kt), temp (C, red), and wind (kt, hatched >= 40 kt)',
              'plot_levs': [[40,60,70,100,120,140,160],
                             np.arange(-90,20,1),
                             np.arange(3000, 9000, 60)]
    },

    '300mb': {'cf_data': 'wspd_300',
              'c1_data': 'z_300',
              'c2_data': 'div_300',
              'c3_data': None,
              'plot_barbs': [True, ('u_300', 'v_300')],
              'plot_info': '300mb height (m, MSL), divergence (10**-5/s, magenta), and wind (kt, hatched >= 60kt)',
              'plot_levs': [[60,80,100,120,140,160,180],
                            np.arange(9000, 15000, 60),
                            np.arange(2,52,2)]

    },

    '250mb': {'cf_data': 'wspd_250',
              'c1_data': 'z_250',
              'c2_data': 'div_250',
              'c3_data': None,
              'plot_barbs': [True, ('u_250', 'v_250')],
              'plot_info': '250mb height (m, MSL), divergence (10**-5/s, magenta), and wind (kt, hatched >= 60kt)',
              'plot_levs': [[60,80,100,120,140,160,180],
                            np.arange(9000, 15000, 60),
                            np.arange(2,52,2)]
}

    #'250mb': [
    #          [60,80,100,120,140,160,180],
    #          np.arange(9000, 15000, 60),
    #          np.arange(60,310,10),
    #          np.arange(2,102,2),
    #          '250mb height (m, MSL, black), divergence (10**-5/s, magenta), and wind (kt, hatched >= 60 kt)'
    #],

    #'500mb': [
    #          [40,60,70,100,120,140,160],
    #          np.arange(3000, 9000, 60),
    #          [40,60,70,100,120,140,160,180,200],
    #          np.arange(-90,20,1),
    #          '500mb height (m, MSL, black), temp (C, red), and wind (kt, hatched >= 40 kt)'
    #]

    #'vadv': [
    #          [10,12,16,20,24,28,32,36,40],
    #          np.arange(3000, 9000, 60),
    #          np.arange(10,100,2),
    #          np.arange(-50,0,5),
    #          np.arange(5,100,5),
    #          '500 mb height, abs. vorticity (fill, 10**5/s-1), 700-400 mb differential vorticity advection'
    #]
}




PLOT_KWARGS = {
    'mucape': [dict(colors=['#9ef94e','#f7d549','#c3882f'], alpha=alpha),
               dict(colors='k', linewidths=1, linestyles=[(0, (15,9))]),
               dict(colors=colors['mucape'], linewidths=[1,1,2,2,2,2,2,2,2,2,2,2]),
               dict()],

    'mlcape': [dict(colors=['#4caee3', '#69e7e7'], alpha=alpha),
               dict(colors='#212dd0', linewidths=1, linestyles=[(0,(15,9))]),
               dict(colors=colors['mlcape'], linewidths=[1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]),
               dict()],

    'eshr':   [dict(),
               dict(colors='#244f88', linewidths=[1,1,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]),
               dict(), dict()],

    'effh':   [dict(colors=['#f1b1bc','#eb7370','#da4238'], alpha=alpha),
               dict(colors=['#559df5','#559df5','#559df5','#244587','#244587','#244587',
                            '#244587','#244587','#244587','#244587','#244587','#244587',
                            '#244587','#244587','#244587','#244587','#244587','#244587'],
                            linewidths=[0.5,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2]),
               dict(), dict()],

    'stpc':   [dict(colors=['#4caee3', '#69e7e7'], alpha=alpha),
               dict(colors=['#ed8433','#ed8433','#ed8433','#e83224','#e83224','#e83224',
                            '#801811','#801811','#801811','#801811','#801811','#801811'],
                    linewidths=[1,1.5,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],
                    linestyles=['--','--','-','-','-','-','-','-','-','-','-','-','-']),
               dict(), dict()],

    '925mb':  [dict(colors=['#3d8825','#5dc83c','#a1fa4e','#6bebec',
                           '#4db1e9', '#4191f7','#efec4e'], alpha=alpha),
               dict(colors='#e83224', linewidths=1, linestyles=[(0,(15,9))]),
               dict(colors='#1b34c0', linewidths=1, linestyles=[(0,(15,9))]),
               dict(colors='k', linewidths=2)],

    '850mb':  [dict(colors=['#3d8825','#5dc83c','#a1fa4e','#6bebec',
                           '#4db1e9', '#4191f7','#efec4e'], alpha=alpha),
               dict(colors='#e83224', linewidths=1, linestyles=[(0,(15,9))]),
               dict(colors='#1b34c0', linewidths=1, linestyles=[(0,(15,9))]),
               dict(colors='k', linewidths=2)],

    '700mb':  [dict(colors=['#3d8726'], alpha=alpha),
               dict(colors='#e83224', linewidths=[1,1,1,1,2,2,2,2,2,2,2,2,2,2,2],
                    linestyles=[(0,(15,9))]),
               dict(colors='#1b34c0', linewidths=1., linestyles=[(0,(15,9))]),
               dict(colors='k', linewidths=2)],

    '500mb':  [dict(colors=['#4293f5','#4db1e7','#6bebeb','#866fc6','#8840ef',
                          '#802088'], alpha=alpha),
               dict(colors='#e83224', linewidths=1., linestyles=[(0,(15,9))]),
               dict(colors='k', linewidths=2),dict()],

    '300mb': [dict(colors=['#4293f5','#4db1e7','#6bebeb','#866fc6','#8840ef',
                           '#802088'], alpha=alpha),
              dict(colors='k', linewidths=2),
              dict(colors='#e740e7', linewidths=2.5)],

    '250mb': [dict(colors=['#4293f5','#4db1e7','#6bebeb','#866fc6','#8840ef',
                           '#802088'], alpha=alpha),
              dict(colors='k', linewidths=2),
              dict(colors='#e740e7', linewidths=2.5)]
}



"""
plot_kwargs = {
    'mucape':[dict(colors=['#9ef94e','#f7d549','#c3882f'], alpha=alpha),
              dict(colors='k', linewidths=1, linestyles=[(0,(15,9))]),
              dict(colors=cape_cols, linewidths=[1,1,2,2,2,2,2,2,2,2,2,2]),
              dict()],

    '250mb': [dict(colors=['#4293f5','#4db1e7','#6bebeb','#866fc6','#8840ef',
                           '#802088'], alpha=alpha),
              dict(colors='k', linewidths=1.7),
              dict(colors='#255284', linewidths=1.),
              dict(colors='#e740e7', linewidths=1.7)],

    '500mb': [dict(colors=['#4293f5','#4db1e7','#6bebeb','#866fc6','#8840ef',
                          '#802088'], alpha=alpha),
              dict(colors='k', linewidths=1.7),
              dict(colors='#2f67a3', linewidths=1.),
              dict(colors='#e83224', linewidths=1., linestyles=[(0,(12,4))])],

    '700mb': [dict(colors=['#3d8726'], alpha=alpha),
              dict(colors=['#5ec93c'], linewidths=1.7),
              dict(colors='k', linewidths=1.7),
              dict(colors='#1b34c0', linewidths=1., linestyles=[(0, (12,4))]),
              dict(colors='#e83224', linewidths=[1,1,1,1,2,2,2,2,2,2,2,2,2,2,2],
                   linestyles=[(0, (12,4))])],

    'vadv': [dict(colors=['#397e23','#5ac039','#9ef94e','#edec4f','#f7d549',
                          '#ed8433', '#e63524', '#e63524'], alpha=alpha),
              dict(colors='k', linewidths=2.3),
              dict(colors='k', linewidths=1.),
              dict(colors='r', linewidths=1.5, linestyles=[(0, (15,9))]),
              dict(colors='b', linewidths=1.5, linestyles=[(0, (15,9))])],

    '700mb2': [dict(colors=['#4293f5','#4db1e7','#6bebeb','#866fc6','#8840ef',
                          '#802088'], alpha=alpha),
              dict(colors='k', linewidths=1.7),
              dict(colors='#255284', linewidths=1.),
              dict()],

    '850mb': [dict(colors=['#3d8825','#5dc83c','#a1fa4e','#6bebec',
                           '#4db1e9', '#4191f7','#efec4e'], alpha=alpha),
              dict(colors='k', linewidths=1.7),
              dict(colors='#1b34c0', linewidths=1., linestyles=[(0,(15,9))]),
              dict(colors='r', linewidths=1, linestyles=[(0,(15,9))])],

    '850mb2': [dict(colors=['#4293f5','#4db1e7','#6bebeb','#866fc6','#8840ef',
                           '#802088'], alpha=alpha),
              dict(colors='k', linewidths=1.7),
              dict(colors='#1b34c0', linewidths=1., linestyles=[(0,(15,5))]),
              dict(colors='r', linewidths=1, linestyles=[(0,(15,9))])],

    '925mb': [dict(colors=['#3d8825','#5dc83c','#a1fa4e','#6bebec',
                           '#4db1e9', '#4191f7','#efec4e'], alpha=alpha),
              dict(colors='k', linewidths=1.7),
              dict(colors='#1b34c0', linewidths=1., linestyles=[(0,(15,9))]),
              dict(colors='r', linewidths=1, linestyles=[(0,(15,9))])],

    '700mb_vort': [dict(colors=vort_cols, alpha=alpha),
                   dict(colors='k', linewidths=1.7),
                   dict(colors='#1b34c0', linewidths=1., linestyles=[(0,(15,5))]),
                   dict(colors='r', linewidths=[1,1,1,1,2,2,2,2,2,2,2,2,2,2,2],
                         linestyles=[(0,(15,9))])],

    'mnwd': [dict(colors=['#4293f5','#4db1e7','#6bebeb','#866fc6','#8840ef',
                           '#802088'], alpha=alpha),
              dict(colors='k', linewidths=1.7), dict(), dict()],

    'prop': [dict(colors=['#9df54d','#5dc63c','#3c8526'], alpha=alpha),
             dict(colors='k', linewidths=0.5),
             dict(colors=['r', 'b'], linewidths=[2,2]), dict(),
             dict(colors='k', linewidths=3)],

    'mlcp': [dict(colors=['#4caee3', '#69e7e7'], alpha=alpha),
             dict(colors=cape_cols, linewidths=[1,1,2,2,2,2,2,2,2,2,2,2]),
             dict(colors='b', linewidths=0.5, linestyles=[(0,(15,9))]), dict(), dict()],

    'laps': [dict(colors=['#f7d0a7'], alpha=alpha), dict(colors=['#669c4e',
             '#ed8433','#ed8433','#db453a','#db453a','#db453a','#db453a','#db453a'],
             linewidths=[1,2,2,2,2,2,2,2,2,2]),dict(colors=['#669c4e',
             '#ed8433','#ed8433','#db453a','#db453a','#db453a','#db453a','#db453a'],
             linewidths=[1,2,2,2,2,2,2,2,2,2]), dict(), dict()],

    'shr6': [dict(colors=['#4caee3', '#69e7e7'], alpha=0),
             dict(colors='#244f88', linewidths=[1,1,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]),
             dict(colors='b', linewidths=0.5, linestyles=[(0,(15,9))]), dict(), dict()],

    'qlcs1': [dict(colors=['#81472b','#c3882f','#f5da4a','#c5c542','#f1b2bd'], alpha=alpha),
              dict(colors=cape_cols[::-1], linewidths=[1,1,2,2,2,2,2,2,2,2,2,2][::-1]),
              dict(),dict()]
}
"""

"""
qpe_levs = [.01,.02,.05,.08,.1,.15,.2,.25,.3,.4,.5,.6,.7,.8,.9,1,1.2,1.4,1.6,1.8,
        2,2.5,3,3.5,4,5,6,7,8,9,10,11,12,13,14,15,16,17,
        18,19,20,21,22]
qpe_cols = ['#D5D5D5','#BEBEBD','#A5A5A5','#838383','#CAFEC1','#B7FAAE','#7DF378',
        '#3FCD45','#29B12B','#1D9D1F','#1C65CD','#2E83E8','#53A4EE','#97D2F6',
        '#B4EDF6','#FEFBD4','#FFFAAF','#FEE880','#FDC04B','#FDA028','#FC6120',
        '#FC351D','#DE1917','#BF0712','#A3050D','#850309','#633C33','#8C645B',
        '#B48D84','#C69F96','#EFDFD7','#CFC9DE','#BEB5D3','#9B8ABB','#8774AD',
        '#725EA1','#750D74','#8C118B','#AF19AF','#BF1CBF','#CB1ECA','#DF22DE']


mrms_levs = np.arange(5,80,5)
mrms_cols = ['#9b9b9a','#cccbcc','#4f9af8','#0026f5','#75f94c',
             '#59c139','#fffd55','#f8cc46','#ee6f2e','#eb3323','#bd401e','#8d1a10',
             '#ea40f7','#8d41c5']
"""

# Merge the plotting dictionaries
PLOT_DICT = {**SHARPPY_DICT, **ROUTINE_DICT}
#PLOT_DICT = {**ROUTINE_DICT}

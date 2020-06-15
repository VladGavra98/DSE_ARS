'''
	AE3200 Design Synthesis Exercise
	Group 09 - Autonomous Environmental Sensing

	Authors: Widmann Sebastian
	Created: 12.06.2020
	
	Project Supervisors:
		- Dr. Irene C. Dedoussi
		- Dr. Ir. Mirjam Snellen
		- Ir. Lorenzo Pasqualetto Cassinis
		- Mark Schelbergen

	This script is used for the data processing of the AQ data around Schiphol airport	
'''

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from math import pi, exp, sqrt
from pathlib import Path
from sympy.solvers import solve
from sympy.abc import xi, Q
import datetime
from scipy.stats import zscore

# ==============================================================================================
# Set global plotting parameters
# ==============================================================================================
# texpsize= [18,20,22]

# plt.rc('font', size=texpsize[1], family='serif')                        # controls default text sizes
# plt.rc('axes', titlesize=texpsize[1])                                   # fontsize of the axes title
# plt.rc('axes', labelsize=texpsize[1])                                   # fontsize of the x and y labels
# plt.rc('xtick', labelsize=texpsize[0])                                  # fontsize of the tick labels
# plt.rc('ytick', labelsize=texpsize[0])                                  # fontsize of the tick labels
# plt.rc('legend', fontsize=texpsize[0])                                  # legend fontsize
# plt.rc('figure', titlesize=texpsize[2])                                 # fontsize of the figure title
# matplotlib.rcParams['lines.linewidth']  = 1.5
# matplotlib.rcParams['figure.facecolor'] = 'white'
# matplotlib.rcParams['axes.facecolor']   = 'white'
# matplotlib.rcParams["legend.fancybox"]  = False

# ==============================================================================================
# Read Measurement Data from Pickle into DataFrame
# ==============================================================================================
colors = ['tab:blue', 'tab:orange', 'tab:grey']

schiphol_data = {'CO'  : {'s1': 'Badhoevedorp'},
				 'NO'  : {'s1': 'Badhoevedorp', 's2': 'Hoofddorp', 's3': 'OudeMeer'},
				 'NO2' : {'s1': 'Badhoevedorp', 's2': 'Hoofddorp', 's3': 'OudeMeer'},
				 'O3'  : {'s2': 'Hoofddorp'},
				 'PM25': {'s1': 'Badhoevedorp'},
				 'PM10': {'s1': 'Badhoevedorp', 's2': 'Hoofddorp', 's3': 'OudeMeer'},
				 'PS'  : {'s1': 'Badhoevedorp', 's2': 'Hoofddorp', 's3': 'OudeMeer'}}

pollutants = {}
times = [datetime.time(i,0) for i in range(6,24)]

for species, stations in schiphol_data.items():
	pollutant = pd.DataFrame() # initialise dataframe

	for key, station in stations.items():
		filedir = 'AQ_SchipholData'
		file = station + '_' + species + '.pickle'
		path = Path(filedir, file) # merge file direction and file name

		df = pd.read_pickle(path) # read data from pickle into dataframe
		df.value = df.value.astype('float64') # convert data type for concentration back from object to flaot64
		df.value = df.value * 1E-6 # convert concentrations from [microgram] to [gram]

		try:
			dropindex = df.index[df.timestamp_measured.str.startswith('2020')][0] # get start index for 2020 measurements
			df = df[df.index < dropindex] # remove measurements for 2020
		except IndexError:
			pass

		df.timestamp_measured = pd.to_datetime(df.timestamp_measured, format='%Y-%m-%dT%H:%M:%S+00:00')
		df = df[df['timestamp_measured'].dt.time.isin(times)] # remove times between 23:00:00 - 06:00:00

		# print(df.value.max())

		df = df[(np.abs(zscore(df['value'])) < 4)] # remove outliers from dataset which are below/above 3 std. deviations
		# df = df[df.value != 0] # remove rows which measured no concentration; negative concentrations are not excluded
		df = df.reset_index() # reset indices of dataframe
		df = df.drop(['index', 'formula'], axis=1) # remove unimportant columns

		pollutant[station] = df.value

	pollutants[species] = pollutant

# ==============================================================================================
# Plot Data in Boxplots
# ==============================================================================================
exportfolder = 'AQ_Boxplots'

# for key, value in pollutants.items():
	# fig1, ax1 = plt.subplots(nrows=1, ncols=1, squeeze=False, figsize=(8,6)) # create subplots

	# ylabel = 'Concentration ' + (r'[particles/$m^3$]' if key == 'PS' else r'[g/$m^3$]')

	# value.boxplot(ax=ax1[0,0]) # plot boxplots
	# ax1[0,0].set_ylabel(ylabel)

	# exportname = str(key) + '_boxplot.png'
	# exportpath = Path(exportfolder, exportname)
	# fig1.savefig(exportpath, bbox_inches='tight', dpi=300)

# ==============================================================================================
# Gaussian Plume Model
# ==============================================================================================
Q = pd.Series([24.912, # CH4
			  27.84, # CO
			  237.778, # NO
			  237.778, # NO2
			  float(pollutants['O3'].max()), # CO
			  0.373, #PM2.5
			  0.746, # PM10
			  max(pollutants['PS'].max()), # UFP
			  1.05]) #SO2

unit1 = pd.Series(['g/s','g/s','g/s','g/s','g/s','g/s','g/s','p/s', 'g/s'])

unit2 = pd.Series(['g/m3','g/m3','g/m3','g/m3','g/m3','g/m3','g/m3','p/m3','g/m3'])

emissions = pd.DataFrame()
emissions['Q'], emissions['Unit of Q'] = Q, unit1

V = 12 # [m/s] freestream velocity
x = 106 # [m] minimum distance between the aircraft and the drone
y, z, h = 0, 0, 0

sigma_y = 0.08 * x * (1 + 0.0001 * x)**(-0.5) # dispersion coefficient based on literature
sigma_z = 0.06 * x * (1 + 0.0015 * x)**(-0.5) # dispersion coefficient based on literature
sigma_y = sigma_y if sigma_y > sigma_z else sigma_z # set sigma_y = sigma_z if sigma_y < sigma_Z

emissions['Xi'] = emissions.Q / (2 * pi * sigma_y * sigma_z * V) * exp( - y**2 / (2 * sigma_y**2) - (z - h)**2 / (2 * sigma_z**2))
emissions['Unit of Xi'] = unit2

R = 0.082057366080960 # [L*atm/(K*mol)]
T = 283 # [K] averaged yearly temperature around Schiphol airport
p = 1 # [atm] ambient pressure

molarmasses = pd.Series([16.043, 28.010, 30.006, 46.006, 47.996, np.nan, np.nan, np.nan, 64.064])

emissions['Xi [ppmm]'] = emissions['Xi'] * R * T * 1E3 / (p * molarmasses) # convert gram/m3 to ppm by weight

print(emissions) # print dataframe
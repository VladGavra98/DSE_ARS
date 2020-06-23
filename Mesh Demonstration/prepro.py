# -*- coding: utf-8 -*-
'''
	AE3200 Design Synthesis Exercise
	Group 09 - Autonomous Environmental Sensing

	@author: vgavra
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
from math import pi, exp
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
		filedir = r"C:\Users\vladg\OneDrive\Documents\GitHub\DSE_ARS\PlumeModelling\AQ_SchipholData"
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

		df.timestamp_measured = pd.to_datetime(df.timestamp_measured)
		df = df[df['timestamp_measured'].dt.time.isin(times)] # remove times between 23:00:00 - 06:00:00

		df = df[(np.abs(zscore(df['value'])) < 3)] # remove outliers from dataset which are below/above 3 std. deviations
		# df = df[df.value != 0] # remove rows which measured no concentration; negative concentrations are not excluded
		df = df.reset_index() # reset indices of dataframe
		df = df.drop(['index', 'formula'], axis=1) # remove unimportant columns

		pollutant[station] = df.value

	pollutants[species] = pollutant

g= open("NO_data.txt",'w')

out = pollutants['NO'].values
np.savetxt(g, out)


g.close()
# ==============================================================================================
# Plot Data in Boxplots
# ==============================================================================================
# concentration = pd.DataFrame(columns=['Species', 'Station', 'Max', 'Unit']) # initialise dataframe for plotting
# index = 0

# exportfolder = 'AQ_Boxplots'

# for key, value in pollutants.items():
#  	fig1, ax1 = plt.subplots(nrows=1, ncols=1, squeeze=False, figsize=(8,6)) # create subplots

#  	ylabel = 'Concentration ' + (r'particles/$m^3$' if key == 'PS' else r'g/$m^3$')

#  	value.boxplot(ax=ax1[0,0]) # plot boxplots
#  	ax1[0,0].set_ylabel(ylabel)

#  	exportname = str(key) + '_boxplot.png'
#  	exportpath = Path(exportfolder, exportname)
#  	fig1.savefig(exportpath, bbox_inches='tight', dpi=300)

#  	xi_max, pos = value.max(), value.max().idxmax(axis=1)
#  	unit = 'p/m^3' if key == 'PS' else 'g/m^3' # print string of unit into dataframe
#  	concentration.loc[index] = [key, pos, max(xi_max), unit]
#  	index += 1


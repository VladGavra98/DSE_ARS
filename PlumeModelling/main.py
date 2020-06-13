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

	This script is used for the propeller and motor sizing of the coaxial quadcopter		
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import pi, exp
from pathlib import Path
from scipy.stats import norm



# ==============================================================================================
# Gaussian Plume Model
# ==============================================================================================
V = 12 # [m/s] freestream velocity
xi = lambda x, y, z, h, Q, Ky, Kz: Q / (2 * pi * V) * (2 * Ky * x / V)**(-0.5) * (2 * Kz * x / V)**(-0.5) * exp( - y**2 / 2 * (2 * Ky * x / V)**(-1) - (z - h)**2 / 2 * (2 * Kz * x / V)**(-1)) #[g/m^3] pollutant concentration

# ==============================================================================================
# Read Measurement Data from Pickle into DataFrame
# ==============================================================================================
# CO_stations 	= ['Badhoevedorp']
# NO_stations 	= ['Badhoevedorp', 'Hoofddorp', 'Oude_Meer']
# NO2_stations 	= ['Badhoevedorp', 'Hoofddorp', 'Oude_Meer']
# O3_stations 	= ['Hoofddorp']
# PM25_stations 	= ['Badhoevedorp']
# PM10_stations 	= ['Badhoevedorp', 'Hoofddorp', 'Oude_Meer']
PS_stations 	= ['Badhoevedorp', 'Hoofddorp', 'Oude_Meer']

species = ['CO', 'NO', 'NO2', 'O3', 'PM25', 'PM10', 'PS']
species = ['PS']
# stations_species = [CO_stations, NO_stations, NO2_stations, O3_stations, PM25_stations, PM10_stations, PS_stations]
stations_species = [PS_stations]

df_CO 	= pd.DataFrame(columns=['concentration', 'timestamp'])
df_NO 	= pd.DataFrame(columns=['concentration', 'timestamp'])
df_NO2 	= pd.DataFrame(columns=['concentration', 'timestamp'])
df_O3 	= pd.DataFrame(columns=['concentration', 'timestamp'])
df_PM25 = pd.DataFrame(columns=['concentration', 'timestamp'])
df_PM10 = pd.DataFrame(columns=['concentration', 'timestamp'])
df_PS 	= pd.DataFrame()

fig1, ax1 = plt.subplots(nrows=3, ncols=1, squeeze=False, figsize=(30,10))

for index1, stations in enumerate(stations_species):
	for index2, station in enumerate(stations):
		filedir = 'AQ_SchipholData'
		file = str(station) + '_' + species[index1] + '.pickle'
		path = Path(filedir, file)

		df = pd.read_pickle(path) # import measurement data from pickle
		df = df.drop(['formula'], axis=1) # remove column indicating species
		df.columns = ['concentration', 'timestamp'] # rename column headers
		df.concentration = df.concentration / 1000 # convert from milligram to gram

		# Get index where measurements for 2020 start - if no measurements taken in 2020, skip command
		try:
			dropindex = df.index[df.timestamp.str.startswith('2020')][0]
			df = df[df.index < dropindex]
		except IndexError:
			pass

		mu 		= df.concentration.mean()
		sigma 	= df.concentration.std()
		x 		= np.linspace(mu - 3 * sigma, mu + 3 * sigma, 1000)

		ax1[index2, 0].plot(x, norm.pdf(x, mu, sigma))
		# print(index2)


plt.show()


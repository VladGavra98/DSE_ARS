'''
	AE3200 Design Synthesis Exercise
	Group 09 - Autonomous Environmental Sensing

	Authors: Widmann Sebastian
	Created: 15.05.2020
	
	Project Supervisors:
		- Dr. Irene C. Dedoussi
		- Dr. Ir. Mirjam Snellen
		- Ir. Lorenzo Pasqualetto Cassinis
		- Mark Schelbergen

	This script is used for the propeller and motor sizing of the coaxial quadcopter		
'''

# Import modules
from math import pi, sqrt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define constants
m 			= 8.856 # [kg] Maximum take-off mass
g 			= 9.80665 # [m/s^2] Gravitational Acceleration
n_prop 		= 8 # [-] Number of propellers
rho 		= 1.225 # [kg/m^3] sea-level density
mu 			= 1.81E-5 # [Ns/m^2] Dynamic visocity of air
a 			= 343 # [m/s] Speed of sound of air
D 			= 0.0254 * np.array([22, 24, 26, 28, 30, 32, 40]) # [m] Propeller Diameter
eta_coaxial = 0.85 # [-] Coaxial efficiency loss
eta_mech 	= 0.8 # [-] Total mechanical efficiency loss
TW 			= 2 # [-] Thrust-to-Weight Ratio

# Define mission constraints
T_hover 	= (m*g) / 8 # [N] Required thrust per propeller during hovering
T_manoeuvre = TW * T_hover # [N] Required thrust per propeller during manoeuvres

# Import propeller data from excel file in the following format
# ['Motor Speed [1/s', 'Thrust [N]', 'Torque [Nm]', 'Voltage[V]', 'Current [A]']
xls 	= pd.ExcelFile("PropellerData_Meizlik.xlsx")
blade22 = pd.read_excel(xls, '22')
blade24 = pd.read_excel(xls, '24')
blade26 = pd.read_excel(xls, '26')
blade28 = pd.read_excel(xls, '28')
blade30 = pd.read_excel(xls, '30')
blade32 = pd.read_excel(xls, '32')

# Store blade dataframes into list for easy iterating
blades 	= [blade22, blade24, blade26, blade28, blade30, blade32]
bladetype = ['22" x 7.4"', '24" x 8.1"', '26" x"8.7"', '28" x 9.4"', '30" x 10"', '32" x 10.6"']

# Initialise dataframe for hovering point which stores data of all propellers
hover 	= pd.DataFrame(columns=['RPM [1/s]', 'Power [W]', 'FOM [-]', 'Disk Loading [N/m^2]'])

 # Initialise figure with 2 rows and 2 columns
fig1, ax1 = plt.subplots(nrows=2, ncols=2, squeeze=False, figsize=(10,8))

# Iterate over dataframes to add columns
# Mechanical power = 2 * pi * n * torque
# Electrical power = voltage * current
# Thrust coefficient = thrust / (rho * n**2 * D**4)
# Power coefficient = power / (rho * n**3 * D**5)
#
for index, df in enumerate(blades):
	P_mech 	= 2 * pi * df['Motor Speed [1/s]'] * df['Torque [Nm]'] / eta_coaxial
	P_elec 	= df['Voltage [V]'] * df['Current [A]']
	Ct 		= df['Thrust [N]'] / (rho * df['Motor Speed [1/s]']**2 * D[index]**4)
	Cp 		= P_mech / (rho * df['Motor Speed [1/s]']**3 * D[index]**5)
	FOM 	= Ct**(3/2) / (sqrt(2) * Cp)

	# c 		= D[index]/2 * 0.121867779 # [m] chord at r/R = 0.75
	# J 		= (2 * pi * df['Motor Speed [1/s]'] * c) / (2 * pi * df['Motor Speed [1/s]'] * D[index])
	# eta 	= J * Ct / Cp
	# Re 		= rho * (2 * pi * df['Motor Speed [1/s]']) * c / mu # [-] Reynolds number at given rpm
	# M 		= (2 * pi * df['Motor Speed [1/s]'] * D[index] * 0.5) / a # [-] Mach number at given rpm for tip speeds

	# Insert columns to dataframe
	df['Mechanical Power [W]'] = P_mech
	df['Electrical Power [W]'] = P_elec
	df['Thrust Coefficient [-]'] = Ct
	df['Power Coefficient [-]'] = Cp
	df['Figure of Merit [-]'] = FOM
	# df['Advance Ratio [-]'] = J
	# df['Efficiency [-]'] = eta
	# df['Reynolds Number [-]'] = Re
	# df['Mach Number [-]'] = M

	# Remove first row from dataframe
	df = df.iloc[1:]

	# Calculate operating point for hovering flight
	rpm_hover 			= np.interp(T_hover, df['Thrust [N]'], df['Motor Speed [1/s]'])
	P_hover 			= np.interp(rpm_hover, df['Motor Speed [1/s]'], df['Mechanical Power [W]'])
	P_hover 			= P_hover / eta_mech
	FOM_hover 			= np.interp(rpm_hover, df['Motor Speed [1/s]'], df['Figure of Merit [-]'])
	DL_hover 			= T_hover / (pi / 4 * D[index])
	hover.loc[index] 	= [rpm_hover*60, P_hover, FOM_hover, DL_hover] # add row to dataframe

	# Polynomial fit to C_T vs. FOM curve
	pfit = np.polyfit(df['Thrust Coefficient [-]'].to_numpy(), df['Figure of Merit [-]'].to_numpy(), deg=3)
	trendline = np.poly1d(pfit)

# 	# Plot advance ratio against efficiency for each propeller
	ax1[0,0].plot(df['Motor Speed [1/s]'], df['Thrust [N]'], marker='o', label=bladetype[index])
	ax1[0,0].set_xlabel('Rotational Speed Motor $\Omega$ [1/s]')
	ax1[0,0].set_ylabel('Thrust [N]')
	ax1[0,0].minorticks_on() # set minor ticks
	ax1[0,0].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
	ax1[0,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid

	ax1[0,1].plot(df['Motor Speed [1/s]'], df['Mechanical Power [W]'], marker='o')
	ax1[0,1].set_xlabel('Rotational Speed Motor $\Omega$ [1/s]')
	ax1[0,1].set_ylabel('Mechanical Power [W]')
	ax1[0,1].minorticks_on() # set minor ticks
	ax1[0,1].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
	ax1[0,1].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid

	# ax1[1,0].scatter(df['Thrust Coefficient [-]'], df['Figure of Merit [-]'])
	ax1[1,0].plot(df['Thrust Coefficient [-]'], trendline(df['Thrust Coefficient [-]']))
	ax1[1,0].set_xlabel('Thrust Coefficient $c_T$ [-]')
	ax1[1,0].set_ylabel('Figure of Merit [-]')
	ax1[1,0].minorticks_on() # set minor ticks
	ax1[1,0].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
	ax1[1,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid

fig1.legend()

plt.show()
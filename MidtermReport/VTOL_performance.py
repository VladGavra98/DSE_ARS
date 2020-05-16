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

	This script is used for the initial sizing of the propellers / aerodynamic surfaces for the concept selection		
'''

# Import modules
from math import sqrt, pi, radians, cos, sin
import numpy as np
from scipy.optimize import fsolve

# ---------------------------------------------------------------------------------------------
# VTOL Flying Wing
# ---------------------------------------------------------------------------------------------
# Battery assumptions for Li-S technology
n 		= 6											# [-] number of cells in series
C_cell 	= 19 										# [Ah] cell capacity
V_cell 	= 2.1 										# [V] nominal cell voltage
C_batt  = n * C_cell 								# [Ah] battery capacity
V_batt 	= n * V_cell 								# [V] battery voltage
C_batt 	= C_batt * V_batt 							# [Wh] battery capacity converted from Ah to Wh

# Mass breakdown
g 		= 9.81 										# [m/s^2] gravitational acceleration
eta_p	= 0.8 										# [-] propeller efficiency
# E 		= 10800										# [s] endurance
# R 		= 1000 										# [m] range
CLCD 	= 14 										# [-] lift-to-drag ratio (taken from https://pubs.acs.org/doi/abs/10.1021/acsenergylett.8b02195)
W_PL	= 1.935										# [kg] payload mass
W_B 	= n * 0.141 								# [kg] battery weight for 6S-configuration

# f_batt 	= 1.3 * (g / (eta_p * E) * R / (CLCD)) 		# [-] battery weight fraction

func	= lambda W_TO: -W_TO + W_PL * 2.20462 / (1 - (W_B * 2.20462 / W_TO) - (5.1E-6*W_TO + 0.42))
W_TO 	= float(0.453592 * fsolve(func, 10))		# [kg] take-off weight
W_E 	= 0.453592 * W_TO * (5.1E-6*W_TO + 0.42) 	# [kg] empty weight

# Power calculation
f 		= 1.03 										# [-] fuselage downwash correction
FoM 	= 0.7 										# [-] figure of merit 
rho 	= 1.225 									# [kg/m^3] sea-level density	
V_climb	= 0.00508 * 500 							# [m/s] required minimum climb speed
eta_m 	= 0.9 										# [-] electromechancial efficiency
D 		= 0.0254 * 6 								# [m] propeller diameter for VTOL; input diameter in inches
A 		= 3* pi/4 * D **2 							# [m^2] propeller disk area for 3 propellers
V 		= 25 		 								# [m/s] Cruise speed; minimum required max. speed equal to 50kts / 25m/s
sigma	= 1.0 										# [-] duct expansion ratio 

P_vert 	= 1/eta_m * (f * W_TO / FoM * sqrt(f * (W_TO / A) / (2 * rho)) + 0.5 * W_TO * V_climb)		# [W] required power for vertical flight
P_hori	= 1/(eta_m * eta_p) * (W_TO * V / CLCD)														# [W] required power for horizontal flight
P_hover = sqrt(W_TO ** 3 / (4* sigma *rho * A))														# [W] required power to hover

t_vert 	= C_batt / P_vert * 60						# [s] endurance time in pure vertical flight
t_hori 	= C_batt / P_hori * 60						# [s] endurance time in pure horizontal flight
t_hover	= C_batt / P_hover * 60						# [s] endurance time in pure hover flight

# Print important variables
print('VTOL Flying Wing Caluclations')
print('----------------------------------------------')
print('Initial battery weight =', round(W_B,3), 'kg')
print('Initial empty weight =', round(W_E,3), 'kg')
print('Initial takeoff weight =', round(W_TO,3), 'kg')
print('Required power for vertical flight =', round(P_vert,3), 'W')
print('Required power for cruise flight =', round(P_hori,3), 'W')
print('Required power for hover flight =', round(P_hover,3), 'W')
print('Maximum time in vertical flight =', round(t_vert,3), 's')
print('Maximum time in cruise flight =', round(t_hori,3), 's')
print('Maximum time in hover flight =', round(t_hover,3), 's')
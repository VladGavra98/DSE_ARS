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
# Battery assumptions for LiPo technology (https://www.beslist.nl/speelgoed_spelletjes/d1029505946/ZOP_Power_111_V_7500_mAh_35C_3S_Lipo_Batterij_XT60_Plug_voor_RC_Quadcopter_Auto.html)
C_batt  = 10000E-3 									# [Ah] battery capacity
V_batt 	= 11.1 										# [V] battery voltage
C_batt 	= C_batt * V_batt							# [Wh] battery capacity converted from Ah to Wh

# Mass breakdown
g 		= 9.81 										# [m/s^2] gravitational acceleration
eta_p	= 0.8 										# [-] propeller efficiency
# E 		= 10800										# [s] endurance
# R 		= 1000 										# [m] range
CLCD 	= 23 										# [-] lift-to-drag ratio (taken from https://pubs.acs.org/doi/abs/10.1021/acsenergylett.8b02195)
W_PL	= 2.394										# [kg] payload mass
W_B 	= 0.448 									# [kg] battery weight
# f_batt 	= 1.3 * (g / (eta_p * E) * R / (CLCD)) 		# [-] battery weight fraction

func	= lambda W_TO: -W_TO + W_PL * 2.20462 / (1 - (W_B * 2.20462 / W_TO) - (-0.00296*W_TO + 0.87))
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

P_vert 	= 1/(eta_m * eta_p) * (f * W_TO / FoM * sqrt(f * (W_TO / A) / (2 * rho)) + 0.5 * W_TO * V_climb)	# [W] required power for vertical flight
P_hori	= 1/eta_p * (W_TO * V / CLCD) + 100																# [W] required power for horizontal flight
P_hover = 1/(eta_m * eta_p) * sqrt((W_TO*g) ** 3 / (4* sigma * rho * A))											# [W] required power to hover

P_PL 	= 100										# [W] required power for payload
P_hover = P_hover + P_PL 							# [W] required hover power with payload power included

t_vert 	= C_batt / P_vert * 60						# [min] endurance time in pure vertical flight
t_hori 	= C_batt / P_hori * 60						# [min] endurance time in pure horizontal flight
t_hover	= C_batt / P_hover * 60						# [min] endurance time in pure hover flight

# Print important variables
print('VTOL Flying Wing Calculations')
print('----------------------------------------------')
print('Required payload weight =', W_PL, 'kg')
print('Initial battery weight =', round(W_B,3), 'kg')
print('Initial empty weight =', round(W_E,3), 'kg')
print('Initial takeoff weight =', round(W_TO,3), 'kg')
print('Initial battery capacity =', C_batt, 'Wh')
print('Required power for vertical flight =', round(P_vert,3), 'W')
print('Required power for horizontal flight =', round(P_hori,3), 'W')
print('Required power for hover flight =', round(P_hover,3), 'W')
print('Maximum time in vertical flight =', round(t_vert,3), 'min')
print('Maximum time in horizontal flight =', round(t_hori,3), 'min')
print('Maximum time in hover flight =', round(t_hover,3), 'min')
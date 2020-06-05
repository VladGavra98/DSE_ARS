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
from math import sqrt, pi, radians, cos, sin, atan
import numpy as np
from scipy.optimize import fsolve

# ---------------------------------------------------------------------------------------------
# Co-axial Octocopter
# ---------------------------------------------------------------------------------------------
# Battery assumptions for LiPo technology (https://www.beslist.nl/speelgoed_spelletjes/d1029505946/ZOP_Power_111_V_7500_mAh_35C_3S_Lipo_Batterij_XT60_Plug_voor_RC_Quadcopter_Auto.html)
C_batt  = 10000E-3 									# [Ah] battery capacity
V_batt 	= 11.1 										# [V] battery voltage
C_batt 	= C_batt * V_batt							# [Wh] battery capacity converted from Ah to Wh

# Mass breakdown
g 		= 9.81 										# [m/s^2] gravitational acceleration
eta_p	= 0.8 										# [-] propeller efficiency
m_PL  	= 0.75										# [kg] payload mass
m_B 	= 0.448 									# [kg] battery weight

func	= lambda m_TO: -m_TO + m_PL * 2.20462 / (1 - (m_B * 2.20462 / m_TO) - (-4.6E-5*m_TO + 0.68)) # coefficients of quad copter is used as structurally more similar to quadcopter than an octocopter
m_TO 	= float(0.453592 * fsolve(func, 10))		# [kg] take-off weight
m_E  	= 0.453592 * m_TO * (-4.6E-5*m_TO + 0.68) 	# [kg] empty weight

# Power calculation
CD_top   = 1.05                                     # assumed as a cube
CD_side  = 0.0613                                   # for i = 15. this is to be iterated
S_top    = 0.5**2                                   # [m2]
S_side   = 0.5*0.25                                 # [m2] thickness of 5cm
r_prop   = 0.35                                     # [m] maximum radius of propeller for quadcopter config with 2cm space between propellers
A        = 8 * (pi * r_prop**2)                     # [m2] total area of the 8 propellers
rho  	 = 1.225 									# [kg/m^3] sea-level density	
V_climb	 = 0.00508 * 500 							# [m/s] required minimum climb speed
eta_m    = 0.9 								 		# [-] electromechancial efficiency
V 	     = 25 		 								# [m/s] Cruise speed; minimum required max. speed equal to 50kts / 25m/s
eta_prop = 0.8
eta_mech = 0.9
eta_coaxial = 0.78125                               # efficiency of coaxial rotors

P_PL 	= 100										# [W] required power for payload

P_vert 	= 1/(eta_prop*eta_mech*eta_coaxial) * (m_TO * g + CD_top * 0.5 * rho * V_climb**2 * S_top) * V_climb		# [W] required power for vertical flight

T_hover = m_TO * g 	 					            # [N] required thrust to hover
P_hover = 1/(eta_prop*eta_mech*eta_coaxial) * sqrt((m_TO*g) ** 3 / (2*rho*A)) # [W] required power to hover

D 		= 0.5*rho*V**2*S_side * CD_side
W 		= m_TO * g
T 		= W

theta 	= atan(D/W) 								# [rad] inclination angle

func2 	= lambda Vi: -Vi + T/(2*rho*A)*1/sqrt((V*cos(theta))**2 + (V*sin(theta) + Vi)**2)
Vi 		= float(fsolve(func2, 0))

P_hori 	= P_PL + 1/(eta_prop*eta_mech*eta_coaxial)*(T*Vi + D*V)	# [W] required power for horizontal flight

P_hover = P_hover + P_PL 							# [W] required hover power with payload power included

t_vert 	= C_batt / P_vert * 60						# [min] endurance time in pure vertical flight
t_hori 	= C_batt / P_hori * 60						# [min] endurance time in pure horizontal flight
t_hover	= C_batt / P_hover * 60						# [min] endurance time in pure hover flight

# Print important variables
print('Coaxial Octocopter Calculations')
print('----------------------------------------------')
print('Required payload weight =', m_PL, 'kg')
print('Initial battery weight =', round(m_B,3), 'kg')
print('Initial empty weight =', round(m_E,3), 'kg')
print('Initial takeoff weight =', round(m_TO,3), 'kg')
print('Initial battery capacity =', C_batt, 'Wh')
print('Required power for vertical flight =', round(P_vert,3), 'W')
print('Required power for horizontal flight =', round(P_hori,3), 'W')
print('Required power for hover flight =', round(P_hover,3), 'W')
print('Maximum time in vertical flight =', round(t_vert,3), 'min')
print('Maximum time in horizontal flight =', round(t_hori,3), 'min')
print('Maximum time in hover flight =', round(t_hover,3), 'min')
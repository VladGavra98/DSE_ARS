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
from math import sqrt, pi, radians, cos, sin, atan, degrees
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

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
W_PL  	= 2.394										# [kg] payload mass
W_B 	= 0.448 									# [kg] battery weight

func	= lambda W_TO: -W_TO + W_PL * 2.20462 / (1 - (W_B * 2.20462 / W_TO) - (-4.6E-5*W_TO + 0.68)) # coefficients of quad copter is used as structurally more similar to quadcopter than an octocopter
W_TO 	= float(0.453592 * fsolve(func, 10))		# [kg] take-off weight
W_E  	= 0.453592 * W_TO * (-4.6E-5*W_TO + 0.68) 	# [kg] empty weight

# Power calculation
CD_top   = 1.05                                     # assumed as a cube
# CD_side  = 0.0613  
CD_side = CD_top                                    # for i = 15. this is to be iterated
S_top    = 0.5**2                                   # [m2]
# S_side   = 0.5*0.25                                 # [m2] 
S_side = 0.3                                        # TO BE UPDATED 
r_prop   = 0.381                                    # [m] maximum radius of propeller for quadcopter config with 2cm space between propellers
A        = 8 * (pi * r_prop**2)                     # [m2] total area of the 8 propellers
rho  	 = 1.225 									# [kg/m^3] sea-level density	
V_climb	 = 0.00508 * 500 							# [m/s] required minimum climb speed
eta_m    = 0.9 								 		# [-] electromechancial efficiency
V 	     = 10 		 								# [m/s] Cruise speed; minimum required max. speed equal to 50kts / 25m/s
eta_prop = 0.8
eta_mech = 0.9
eta_coaxial = 0.78125                               # efficiency of coaxial rotors

P_PL 	= 100										# [W] required power for payload


P_vert 	= (W_TO * g + CD_top * 0.5 * rho * V_climb**2 * S_top) * V_climb / eta_prop / eta_mech / eta_coaxial		# [W] required power for vertical flight

T_hover = W_TO * g 	 					             	# [N] required thrust to hover
P_hover = sqrt((W_TO*g)**3 / (2*rho*A)) / eta_prop / eta_mech / eta_coaxial	# [W] required power to hover											# [W] required power to hover

D 		= 0.5 * rho * V**2 * S_side * CD_side
W 		= W_TO * g

# 
theta 	= atan(D/W)
# theta = radians(15)
T 		= W /cos(theta)

func2 	= lambda Vi: -Vi + T / (2*rho*A) / sqrt((V*cos(theta))**2 + (V*sin(theta) + Vi)**2)
Vi 		= float(fsolve(func2, 0))

P_hori 	= P_PL + 1/(eta_prop*eta_mech*eta_coaxial) * (T*Vi + D*V)	# [W] required power for horizontal flight


P_hover = P_hover + P_PL 							# [W] required hover power with payload power included

t_vert 	= C_batt / P_vert * 60						# [min] endurance time in pure vertical flight
t_hori 	= C_batt / P_hori * 60						# [min] endurance time in pure horizontal flight
t_hover	= C_batt / P_hover * 60						# [min] endurance time in pure hover flight

# Print important variables
print('Coaxial Octocopter Calculations')
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


##### Predefined/Required Specifications
grav = 9.80665 # [m/s2]
vel_cruise = 10 # [m/2] - based on autonomous sensor 
s_2pts = 200 # [m]
s_2pts_max = sqrt(2 * 200**2) # diagonal dist between 2 pts



##### Mass Estimation
M_PL  =  1.5
M_Batt  = 2 * 3.58 # 2 Tattu 32Ah Batteries
M_Frame = 3.8 #  Estimated based on data of YANGDA YD6-1600L Hexacopter
M_tot = (M_PL + M_Batt + M_Frame) * 1.1 # Added contingency
M_TO = M_tot

##### Thrust Calculation
TWratio_original = 2.5 # TO BE VARIED 
Throttle_ratio = 0.50
TWratio = TWratio_original * Throttle_ratio
Thrust = TWratio * M_tot * grav

##### Acceleration Calculation
dt = 0.001
acc_cruise_list = []
vel = 0
vel_data = [0]
t_acc_horz = 0
tilt_angle = radians(15) # TO BE VARIED
Thrust_cruise = Thrust * cos(pi/2 - tilt_angle)
while vel < 10:
    acc_cruise = (Thrust_cruise - 0.5 * CD_side * rho * S_side * vel**2) / M_tot
    vel += acc_cruise * dt
    vel_data.append(vel)
    t_acc_horz += dt
    acc_cruise_list.append(acc_cruise) # Check if acceleration in realistic range
# print(t_acc_horz) # time taken to reach vel_cruise from 0 m/s
s_acc_horz = np.trapz(np.array(vel_data), dx = dt)
# print(s_acc_horz)

##### Time Calculation
t_2pts = 2 * t_acc_horz + (s_2pts_max - 2 * s_acc_horz) / vel_cruise

print(t_2pts)
print("Take-off Mass:", M_TO)
print("Total Thrust during cruise:", Thrust)
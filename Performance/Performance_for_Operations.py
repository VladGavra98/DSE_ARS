"""
@author: Sebastian & Changkyu
"""
# Import my slaves
from math import sqrt, pi, radians, cos, sin, atan, degrees
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# Assumptions
# 1. Acceleration = -1 * Decceleration
# 2. Complete still-hovering during measurements, v_cruise = 0 [m/s]
# 3. The central body of the octocopter is a cuboid


""" Fixed Parameters """

# General Parameters
grav = 9.80665 # [m/s2]
rho_air = 1.225 # [km/m3]


# Fixed Performance Parameters
vel_cruise = 10 # [m/s] - limited by the autonomy sensor
vel_cruise_max = 26 # [m/s] - EASA's requirement
vel_climb_max = 0.00508 * 500 # [m/s] - EASA's requirement
N_prop = 8
P_PL = 25 # [W] - sum of power req by payloads
M_PL = 1.47 # [kg] - sum of mass of payloads


# Fixed Operational Parameters
dist_2pts = 200 # [m]
dist_2pts_max = sqrt(2 * 200**2) # [m] - diagonal dist between 2 pts


""" VARIABLE Parameters - Should be adjusted accordingly with finalised values """

# Efficiencies
# eta_prop = 0.8 # included in EM
eta_EM = 0.8 # Electromechanical Efficiency
eta_coaxial = 0.85 # THIS SHOULD BE CORRECTED
# eta_motor = 0.85 # included in EM


# Battery Parameters
Cap_batt_Ah = 2 * 32 # [Ah] - 2 Tattu 32Ah Batteries
V_batt = 11.1 # [V]
Cap_batt = Cap_batt_Ah * V_batt # [Wh] 
M_Batt = 2 * 2.58 # [kg]


# Structural Parameters
M_Frame = 0.38 # [kg] - Estimated based on data of YANGDA YD6-1600L Hexacopter
M_TO = (M_PL + M_Batt + M_Frame) * 1.2 # [kg] - added 20% contingency
W_TO = M_TO * grav # [N]
CD_top = 1.05 # Assumed as a cube
CD_side = 1.05 # Assumed as a cube
S_top = 0.5*0.5 # [m2]
S_side = 0.3 # [m2]


# Aero & Prop Parameters
r_prop = 0.381 # [m] - propeller with diameter of 30 inches
A_disk = N_prop * pi * r_prop**2 # [m] - Total disk area of all propellers


""" Calculations  """

# Power Calculations
P_vert = (M_TO * grav + CD_top * 0.5 * rho_air * vel_climb_max**2 * S_top) * vel_climb_max / eta_EM / eta_coaxial # [W] required power for vertical flight
P_hover = sqrt((M_TO * grav)**3 / (2 * rho_air * A_disk))  / eta_EM / eta_coaxial	# [W] required power to hover WITHOUT PAYLOAD										# [W] required power to hover

D_cruise	= 0.5 * rho_air * vel_cruise**2 * S_side * CD_side
theta = atan(D_cruise / W_TO)
T_cruise	= W_TO /cos(theta)
func_vel_i 	= lambda vel_i: -vel_i + T_cruise / (2 * rho_air * A_disk) / sqrt((vel_cruise * cos(theta))**2 + (vel_cruise * sin(theta) + vel_i)**2)
vel_i 		= float(fsolve(func_vel_i, 0))
P_cruise 	= (T_cruise * vel_i + D_cruise * vel_cruise) / (eta_EM * eta_coaxial) 	# [W] required power for horizontal flight

T_hover = W_TO

# Endurance Calculations
t_vert 	= Cap_batt / P_vert * 60						# [min] endurance time in pure vertical flight
t_cruise 	= Cap_batt / P_cruise * 60					# [min] endurance time in pure horizontal flight
t_hover	= Cap_batt / P_hover * 60						# [min] endurance time in pure hover flight


""" # Bunch of Print Statements for Power and Endurance results """
# print('Required payload mass =', M_PL, 'kg')
# print('Initial battery mass =', round(M_Batt,3), 'kg')
# print('Initial takeoff mass =', round(M_TO,3), 'kg')
# print('Initial battery capacity =', Cap_batt, 'Wh')
print('Required power for vertical flight =', round(P_vert,3), 'W')
print('Required power for horizontal flight =', round(P_cruise,3), 'W')
print('Required power for hover flight =', round(P_hover,3), 'W')
print('Maximum time in vertical flight =', round(t_vert,3), 'min')
print('Maximum time in horizontal flight =', round(t_cruise,3), 'min')
print('Maximum time in hover flight =', round(t_hover,3), 'min')


# Calculation for Max Thrust
D_cruise_max	= 0.5 * rho_air * vel_cruise_max**2 * S_side * CD_side
theta_cruise_max = atan(D_cruise_max / W_TO)
T_cruise_max	= W_TO /cos(theta_cruise_max)
print("Thrust at 25m/s =", T_cruise_max)

acc_max = 6 # [m/s2] - typical drone's acceleration is 4; max = 6
Thrust_acc_max = T_hover + M_TO * acc_max
print("Thrust to max acc =", Thrust_acc_max)

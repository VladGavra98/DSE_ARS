# Import my slaves
from math import sqrt, pi, radians, cos, sin, atan, degrees, ceil, floor
import pandas as pd
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import matplotlib

################ SENSITIVITY STUFF #################
jump = 0.01
# percent = np.arange(-1.0, 5.0+jump, jump) # When plotting, rmb to subtract 1 and multiply by 100 to shown % change
# percent_dod = np.arange(0.1, 1.0+jump, jump)
# percent_area = np.arange(1.0, 5.0+jump, jump)
percent_op = np.arange(0.0, 6.0+jump, jump)


# # For the table
# percent = np.arange(0.8, 1.5+jump, jump) # When plotting, rmb to subtract 1 and multiply by 100 to shown % change
# percent_dod = percent
# percent_area = percent
# percent_op = percent


#outcome is all the flight time 
struc_out = []
DOD_out = []
SA_out = []
op_out = []

# hover_time = 7
# g 			= 9.80665 # [m/s^2] Gravitational Acceleration
# rho 		= 1.225 # [kg/m^3] sea-level density
# mu 			= 1.81E-5 # [Ns/m^2] Dynamic visocity of air
# a 			= 343 # [m/s] Speed of sound of air
# D 			= 0.0254 * np.array([28, 30, 32]) # [m] Propeller Diameter
# eta_coaxial = 0.85 # [-] Coaxial efficiency loss
# eta_mech 	= 0.85 # [-] Total mechanical efficiency loss

# V_cruise	= 10 # [m/s] cruise speed
# V_max 		= 26 # [m/s] maximum speed
# CD 			= 1.05 # [-] assumed to be a cube
# S 			= 1.9596966763292372 # [m^2] frontal area

# D_cruise 	= 0.5 * rho * V_cruise ** 2 * CD * S # [N] drag force during cruise flight
# D_max 		= 0.5 * rho * V_max ** 2 * CD * S # [N] drag force during cruise flight

# # Mass budget
# n_propellers	= 8 # [-] number of propellers
# m_sensors 		= 1.124 # [kg] noise and air pollution sensors mass
# m_autonomy 		= 0.219 # [kg] autonomy hardware mass
# m_motors 		= n_propellers * 0.630 #[kg] brushless motor mass
# m_propellers 	= n_propellers * 91.46E-3 # [kg] propeller mass

# m_structures 	= (2.289 + 1.4) * percent # [kg] operational empty mass
# for m_structures_sens in m_structures:

#     # Frame data
#     w_frame = 0.230 # [m] frame width
#     h_frame = 0.260 # [m] frame height
    
#     # Motor data
#     V_motor = 36 # [V] nominal voltage of motor
    
#     # Battery configuration - Chosen cell: EEMB LP623454
#     DoD = 0.85 # [-] Depth of Discharge
#     C_cell = 1200E-3 # [Ah] cell capacity
#     m_cell = 22.5E-3 # [kg] cell mass
#     V_cell = 3.7 # [V] nominal cell voltage
    
#     h_cell = 69E-3 # [m] cell height
#     w_cell = 34.5E-3 # [m] cell width 
#     t_cell = 6.2E-3 # [m] cell thickness
    
#     n_series = ceil(V_motor / V_cell)
    
#     C_battery = C_cell # [Ah] battery capacity
#     m_battery = n_series * m_cell # [kg] battery mass
#     V_battery = n_series * V_cell # [V] battery voltage
#     E_battery = V_battery * C_battery # [Wh] battery capacity
    
#     A_frame = w_frame * h_frame
#     A_battery = h_cell * n_series * t_cell
#     A_battery = 1.05 * A_battery # contingency margin for cooling / packaging
    
#     n_max = floor(A_frame / A_battery) - 1
#     n_max = 16
    
#     m_parallel = n_max * m_battery
#     C_parallel = n_max * E_battery
    
    
#     xls 	= pd.ExcelFile("PropellerData_Meizlik.xlsx")
#     blade32 = pd.read_excel(xls, '32')
    
#     m 	= m_sensors + m_autonomy + m_parallel + m_structures_sens + m_motors + m_propellers
    
#     theta_cruise 	= atan(D_cruise / (m * g)) # [rad] inclination angle during cruise
#     theta_max 		= atan(D_max / (m * g)) # [rad] inclination angle during max speed
    
#     # Define mission constraints
#     T_hover 	= (m*g) / n_propellers # [N] Required thrust per propeller during hovering
#     T_cruise 	= (m*g) / cos(theta_cruise) # [N] Required thrust during cruise
    
#     P_mech 	= 2 * pi * blade32['Motor Speed [1/s]'] * blade32['Torque [Nm]'] / eta_coaxial
    
#     blade32['Mechanical Power [W]'] = P_mech
#     pfit_thrust = np.polyfit(blade32['Motor Speed [1/s]'], blade32['Thrust [N]'], deg=3)
#     tl_thrust	= np.poly1d(pfit_thrust)
    
#     pfit_power 	= np.polyfit(blade32['Motor Speed [1/s]'], blade32['Mechanical Power [W]'], deg=3)
#     tl_power	= np.poly1d(pfit_power)
    
#     pfit_torque	= np.polyfit(blade32['Motor Speed [1/s]'], blade32['Torque [Nm]'], deg=3)
#     tl_torque	= np.poly1d(pfit_torque)
    
#     """ Fixed Parameters """
    
#     # General Parameters
#     grav = 9.80665 # [m/s2]
#     rho_air = 1.225 # [km/m3]
    
#     # Fixed Performance Parameters
#     vel_cruise = 10 # [m/s] - limited by the autonomy sensor
#     vel_cruise_max = 25 # [m/s] - EASA's requirement
#     vel_climb_max = 0.00508 * 500 # [m/s] - EASA's requirement
#     N_prop = 8
#     P_PL = 25 # [W] - sum of power req by payloads
#     P_max = 360 * N_prop # [W] - max power output of chosen motor 
    
#     # Fixed Operational Parameters
#     dist_2pts = 150 # [m]
#     dist_2pts_max = sqrt(2 * dist_2pts**2) # [m] - diagonal dist between 2 pts
    
#     # Propeller Parameters
#     D_prop = 0.8128
#     r_prop = D_prop / 2
#     t_prop = 0.04
#     A_disk = N_prop * pi * r_prop**2
    
#     """ VARIABLE Parameters - Should be adjusted accordingly with finalised values """
    
#     # Efficiencies
#     # eta_prop = 0.8 # included in EM
#     eta_EM = 0.85 # Electromechanical Efficiency
#     eta_coaxial = 0.85 # THIS SHOULD BE CORRECTED
    
#     # Battery Parameters
#     Cap_batt_Ah = 192 # [Ah]
#     V_batt = 37 # [V]
#     Cap_batt = Cap_batt_Ah * V_batt / 10 * DoD # [Wh] 
#     M_Batt =  3.6 # [kg]
    
#     # Motor Parameters
#     N_motor = 8
#     D_motor = 0.09
#     t_motor = (27.4 + 12.5) / 1000 
    
#     # Structural Parameters
#     N_arm = 4
#     l_arm = 0.564
#     r_arm = 0.024
#     D_arm = r_arm * 2
#     l_body = 0.351
#     w_body = 0.274
#     h_body = 0.002
    
#     M_TO = m
#     W_TO = M_TO * grav # [N]
#     CD_top = 1.05 # Assumed as a cube
#     CD_side = 1.05 # Assumed as a cube
#     S_top = (0.274 * 0.351) * 1.1  # [m2]
#     S_side_0deg = (N_motor * D_motor * t_motor + N_arm * D_arm * l_arm + h_body*w_body) * 1.1
    
#     """ Calculations for cruise power and tilt angle """
#     # Power Calculations
#     P_PL = 17.467 + 14.660
#     P_hover = sqrt((W_TO)**3 / (2 * rho_air * A_disk))  / eta_EM /eta_coaxial # [W] required power to hover WITHOUT PAYLOAD
#     P_vert = P_hover + CD_top * 0.5 * rho_air * vel_climb_max**2 * S_top * vel_climb_max / eta_EM / eta_coaxial
    
#     # Calculate tilt angle at cruise velocity
#     S_side_guess = 0 # [m2]
#     S_side_outcome = S_top * sin(radians(5)) + S_side_0deg / cos(radians(5))
#     while np.abs(S_side_guess - S_side_outcome) > 1E-5:
        
#         S_side_guess = S_side_outcome
#         D_cruise	= 0.5 * rho_air * vel_cruise**2 * S_side_guess * CD_side
#         theta       = atan(D_cruise / W_TO)
#         T_cruise	= W_TO /cos(theta)
#         S_side_outcome = S_top * sin(theta) +S_side_0deg / cos(theta)
        
#     S_side_cruise = S_side_outcome
#     S_side = S_side_cruise
    
#     rpm_dp2 = np.interp(T_cruise / n_propellers, blade32['Thrust [N]'], blade32['Motor Speed [1/s]'])
#     T_dp2 = tl_thrust(rpm_dp2)
#     P_dp2 = tl_power(rpm_dp2)
#     torque_dp2 = tl_torque(rpm_dp2)
#     P_cruise = (n_propellers * P_dp2 + D_cruise * V_cruise) / eta_mech 
        
#     # Calculte tilt angle at max cruise velocity
    
#     S_side_guess = 0 # [m2]
#     S_side_outcome2 = S_top * sin(radians(10)) + S_side_0deg / cos(radians(10))
#     while np.abs(S_side_guess - S_side_outcome2) > 1E-5:
        
#         S_side_guess = S_side_outcome2
        
#         D_cruise_max	= 0.5 * rho_air * vel_cruise_max**2 * S_side_guess * CD_side
#         theta_max = atan(D_cruise_max / W_TO)
#         T_cruise_max	= W_TO /cos(theta_max)
#         S_side_outcome = S_top * sin(theta_max) + S_side_0deg / cos(theta_max)
#         # print(S_side_outcome)
    
#     S_side_maxcruise = S_side_outcome
#     rpm_dp2 = np.interp(T_cruise_max / n_propellers, blade32['Thrust [N]'], blade32['Motor Speed [1/s]'])
#     P_cruise_max = (n_propellers * P_dp2 + D_cruise * V_cruise) / eta_mech 
#     T_hover = W_TO
    
#     # Endurance Calculations
#     t_vert 	= Cap_batt / P_vert * 60						# [min] endurance time in pure vertical flight
#     t_cruise = Cap_batt / P_cruise * 60					# [min] endurance time in pure horizontal flight
#     t_hover	= Cap_batt / P_hover * 60						# [min] endurance time in pure hover flight
    
#     # Calculation for Max Thrust PREM VALUS
#     D_cruise_max	= 0.5 * rho_air * vel_cruise_max**2 * S_side_maxcruise * CD_side
#     theta_cruise_max = atan(D_cruise_max / W_TO)
#     T_cruise_max	= W_TO /cos(theta_cruise_max)
#     acc_max = 4 # [m/s2] - typical drone's acceleration is 4; max = 6
#     Thrust_acc_max = sqrt(T_hover**2 + M_TO**2* acc_max**2)
    
#     # Iteration Calculation
#     dt = 0.1
#     acc_to_cruise = np.arange(4,20+dt,dt)
#     Energy_data = []
#     t_2pts = []
#     for i in acc_to_cruise:
#         Thrust_acc_to_cruise = np.sqrt(T_hover**2 + (M_TO * i)**2) # Const Thrust level for the specified acceleration
#         t_acc_to_cruise = (vel_cruise - 0) / i # Simple a = (v - u) / t formula
#         dist_acc_to_cruise = 0.5 * i * t_acc_to_cruise**2 # Simple s = ut + 0.5at^2 formula
#         num_time_steps = int(round(t_acc_to_cruise / dt, 0))
#         vel_range = np.linspace(0, vel_cruise, num_time_steps) # array of vel with same as its corresponding time steps array
#         Power_acc_to_cruise = Thrust_acc_to_cruise * vel_range / eta_coaxial / eta_mech # Power = Thrust * velocity
#         Energy_acc_to_cruise = 0.5 * Power_acc_to_cruise[-1] * t_acc_to_cruise # area of triangle
        
#         dist_vel_cruise = dist_2pts - dist_acc_to_cruise * 2
#         t_vel_cruise = dist_vel_cruise / vel_cruise
#         Energy_vel_cruise = P_cruise * t_vel_cruise
        
#         Energy_total = (Energy_acc_to_cruise * 2 + Energy_vel_cruise  + (P_hover +P_PL) * hover_time) / 60**2
#         Energy_data.append(Energy_total)
        
#         t_total = t_acc_to_cruise * 2 + t_vel_cruise
#         t_2pts.append(t_total)
        
#     opt_hori_acc = acc_to_cruise[Energy_data.index(np.min(Energy_data))]
#     opt_hori_acc_energy = np.min(Energy_data) # [Wh]
#     min_energy_plot = np.ones((np.size(acc_to_cruise))) * opt_hori_acc_energy
#     opt_acc_plot = np.ones((np.size(Energy_data))) * opt_hori_acc
    
#     tot_energy_1measure = opt_hori_acc_energy
#     num_measurements = Cap_batt / tot_energy_1measure
#     total_flight_duration = num_measurements * (hover_time + t_total)
#     struc_out.append(total_flight_duration)
    


# percent = np.round(percent, 2)
# percent = percent.tolist()
# index1_percent = percent.index(1.0)
# struc_out = (np.array(struc_out) - struc_out[index1_percent]) / struc_out[index1_percent] * 100
# percent = (np.array(percent)-1)*100


# texpsize= [26,28,30]
# SMALL_SIZE  = texpsize[0]
# MEDIUM_SIZE = texpsize[1]
# BIGGER_SIZE = texpsize[2]
# plt.rc('font', size=MEDIUM_SIZE)                    ## controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)                ## fontsize of the axes title
# plt.rc('axes', labelsize=SMALL_SIZE)                ## fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)               ## legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)             ## fontsize of the figure title
# matplotlib.rcParams['lines.linewidth']  = 1.5
# matplotlib.rcParams['figure.facecolor'] = 'white'
# matplotlib.rcParams['axes.facecolor']   = 'white'
# matplotlib.rcParams["legend.fancybox"]  = False

# fig, ax = plt.subplots(1,1,squeeze=False,figsize=(12,9))
# ax[0,0].plot(percent, struc_out, 'b--', label='Extrapolated')
# ax[0,0].plot(percent[100:300], struc_out[100:300], 'r-')
# ax[0,0].set_xlabel('% change in Structural Mass')
# ax[0,0].set_ylabel('% change in Flight Duration')

# ax[0,0].axvline(0,color=(0,0,0),linewidth=1.3)
# ax[0,0].axhline(0,color=(0,0,0),linewidth=1.3)
# ax[0,0].minorticks_on() # set minor ticks
# ax[0,0].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
# ax[0,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid

# plt.legend()
# fig.savefig("sensitivity_structure.png", bbox_inches='tight')            

# plt.show()












# hover_time = 7
# g 			= 9.80665 # [m/s^2] Gravitational Acceleration
# rho 		= 1.225 # [kg/m^3] sea-level density
# mu 			= 1.81E-5 # [Ns/m^2] Dynamic visocity of air
# a 			= 343 # [m/s] Speed of sound of air
# D 			= 0.0254 * np.array([28, 30, 32]) # [m] Propeller Diameter
# eta_coaxial = 0.85 # [-] Coaxial efficiency loss
# eta_mech 	= 0.85 # [-] Total mechanical efficiency loss

# V_cruise	= 10 # [m/s] cruise speed
# V_max 		= 26 # [m/s] maximum speed
# CD 			= 1.05 # [-] assumed to be a cube
# S 			= 1.9596966763292372 # [m^2] frontal area

# D_cruise 	= 0.5 * rho * V_cruise ** 2 * CD * S # [N] drag force during cruise flight
# D_max 		= 0.5 * rho * V_max ** 2 * CD * S # [N] drag force during cruise flight

# # Mass budget
# n_propellers	= 8 # [-] number of propellers
# m_sensors 		= 1.124 # [kg] noise and air pollution sensors mass
# m_autonomy 		= 0.219 # [kg] autonomy hardware mass
# m_motors 		= n_propellers * 0.630 #[kg] brushless motor mass
# m_propellers 	= n_propellers * 91.46E-3 # [kg] propeller mass

# m_structures 	= (2.289 + 1.4) # [kg] operational empty mass


# # Frame data
# w_frame = 0.230 # [m] frame width
# h_frame = 0.260 # [m] frame height

# # Motor data
# V_motor = 36 # [V] nominal voltage of motor

# # Battery configuration - Chosen cell: EEMB LP623454

# # DoD_sens = 0.85 * percent_dod # [-] Depth of Discharge
# DoD_sens = percent_dod
# for DoD in DoD_sens:
    
#     C_cell = 1200E-3 # [Ah] cell capacity
#     m_cell = 22.5E-3 # [kg] cell mass
#     V_cell = 3.7 # [V] nominal cell voltage
    
#     h_cell = 69E-3 # [m] cell height
#     w_cell = 34.5E-3 # [m] cell width 
#     t_cell = 6.2E-3 # [m] cell thickness
    
#     n_series = ceil(V_motor / V_cell)
    
#     C_battery = C_cell # [Ah] battery capacity
#     m_battery = n_series * m_cell # [kg] battery mass
#     V_battery = n_series * V_cell # [V] battery voltage
#     E_battery = V_battery * C_battery * DoD # [Wh] battery capacity
    
#     A_frame = w_frame * h_frame
#     A_battery = h_cell * n_series * t_cell
#     A_battery = 1.05 * A_battery # contingency margin for cooling / packaging
    
#     n_max = floor(A_frame / A_battery) - 1
#     n_max = 16
    
#     m_parallel = n_max * m_battery
#     C_parallel = n_max * E_battery
    
    
#     xls 	= pd.ExcelFile("PropellerData_Meizlik.xlsx")
#     blade32 = pd.read_excel(xls, '32')
    
#     m 	= m_sensors + m_autonomy + m_parallel + m_structures + m_motors + m_propellers
    
#     theta_cruise 	= atan(D_cruise / (m * g)) # [rad] inclination angle during cruise
#     theta_max 		= atan(D_max / (m * g)) # [rad] inclination angle during max speed
    
#     # Define mission constraints
#     T_hover 	= (m*g) / n_propellers # [N] Required thrust per propeller during hovering
#     T_cruise 	= (m*g) / cos(theta_cruise) # [N] Required thrust during cruise
    
#     P_mech 	= 2 * pi * blade32['Motor Speed [1/s]'] * blade32['Torque [Nm]'] / eta_coaxial
    
#     blade32['Mechanical Power [W]'] = P_mech
#     pfit_thrust = np.polyfit(blade32['Motor Speed [1/s]'], blade32['Thrust [N]'], deg=3)
#     tl_thrust	= np.poly1d(pfit_thrust)
    
#     pfit_power 	= np.polyfit(blade32['Motor Speed [1/s]'], blade32['Mechanical Power [W]'], deg=3)
#     tl_power	= np.poly1d(pfit_power)
    
#     pfit_torque	= np.polyfit(blade32['Motor Speed [1/s]'], blade32['Torque [Nm]'], deg=3)
#     tl_torque	= np.poly1d(pfit_torque)
    
#     """ Fixed Parameters """
    
#     # General Parameters
#     grav = 9.80665 # [m/s2]
#     rho_air = 1.225 # [km/m3]
    
#     # Fixed Performance Parameters
#     vel_cruise = 10 # [m/s] - limited by the autonomy sensor
#     vel_cruise_max = 25 # [m/s] - EASA's requirement
#     vel_climb_max = 0.00508 * 500 # [m/s] - EASA's requirement
#     N_prop = 8
#     P_PL = 25 # [W] - sum of power req by payloads
#     P_max = 360 * N_prop # [W] - max power output of chosen motor 
    
#     # Fixed Operational Parameters
#     dist_2pts = 150 # [m]
#     dist_2pts_max = sqrt(2 * dist_2pts**2) # [m] - diagonal dist between 2 pts
    
#     # Propeller Parameters
#     D_prop = 0.8128
#     r_prop = D_prop / 2
#     t_prop = 0.04
#     A_disk = N_prop * pi * r_prop**2
    
#     """ VARIABLE Parameters - Should be adjusted accordingly with finalised values """
    
#     # Efficiencies
#     # eta_prop = 0.8 # included in EM
#     eta_EM = 0.85 # Electromechanical Efficiency
#     eta_coaxial = 0.85 # THIS SHOULD BE CORRECTED
    
#     # Battery Parameters
#     Cap_batt_Ah = 192 # [Ah]
#     V_batt = 37 # [V]
#     Cap_batt = Cap_batt_Ah * V_batt / 10 * DoD # [Wh] 
#     M_Batt =  3.6 # [kg]
    
#     # Motor Parameters
#     N_motor = 8
#     D_motor = 0.09
#     t_motor = (27.4 + 12.5) / 1000 
    
#     # Structural Parameters
#     N_arm = 4
#     l_arm = 0.564
#     r_arm = 0.024
#     D_arm = r_arm * 2
#     l_body = 0.351
#     w_body = 0.274
#     h_body = 0.002
    
#     M_TO = m
#     W_TO = M_TO * grav # [N]
#     CD_top = 1.05 # Assumed as a cube
#     CD_side = 1.05 # Assumed as a cube
#     S_top = (0.274 * 0.351) * 1.1  # [m2]
#     S_side_0deg = (N_motor * D_motor * t_motor + N_arm * D_arm * l_arm + h_body*w_body) * 1.1
    
#     """ Calculations for cruise power and tilt angle """
#     # Power Calculations
#     P_PL = 17.467 + 14.660
#     P_hover = sqrt((W_TO)**3 / (2 * rho_air * A_disk))  / eta_EM /eta_coaxial # [W] required power to hover WITHOUT PAYLOAD
#     P_vert = P_hover + CD_top * 0.5 * rho_air * vel_climb_max**2 * S_top * vel_climb_max / eta_EM / eta_coaxial
    
#     # Calculate tilt angle at cruise velocity
#     S_side_guess = 0 # [m2]
#     S_side_outcome = S_top * sin(radians(5)) + S_side_0deg / cos(radians(5))
#     while np.abs(S_side_guess - S_side_outcome) > 1E-5:
        
#         S_side_guess = S_side_outcome
#         D_cruise	= 0.5 * rho_air * vel_cruise**2 * S_side_guess * CD_side
#         theta       = atan(D_cruise / W_TO)
#         T_cruise	= W_TO /cos(theta)
#         S_side_outcome = S_top * sin(theta) +S_side_0deg / cos(theta)
        
#     S_side_cruise = S_side_outcome
#     S_side = S_side_cruise
    
#     rpm_dp2 = np.interp(T_cruise / n_propellers, blade32['Thrust [N]'], blade32['Motor Speed [1/s]'])
#     T_dp2 = tl_thrust(rpm_dp2)
#     P_dp2 = tl_power(rpm_dp2)
#     torque_dp2 = tl_torque(rpm_dp2)
#     P_cruise = (n_propellers * P_dp2 + D_cruise * V_cruise) / eta_mech 
        
#     # Calculte tilt angle at max cruise velocity
    
#     S_side_guess = 0 # [m2]
#     S_side_outcome2 = S_top * sin(radians(10)) + S_side_0deg / cos(radians(10))
#     while np.abs(S_side_guess - S_side_outcome2) > 1E-5:
        
#         S_side_guess = S_side_outcome2
        
#         D_cruise_max	= 0.5 * rho_air * vel_cruise_max**2 * S_side_guess * CD_side
#         theta_max = atan(D_cruise_max / W_TO)
#         T_cruise_max	= W_TO /cos(theta_max)
#         S_side_outcome = S_top * sin(theta_max) + S_side_0deg / cos(theta_max)
#         # print(S_side_outcome)
    
#     S_side_maxcruise = S_side_outcome
#     rpm_dp2 = np.interp(T_cruise_max / n_propellers, blade32['Thrust [N]'], blade32['Motor Speed [1/s]'])
#     P_cruise_max = (n_propellers * P_dp2 + D_cruise * V_cruise) / eta_mech 
#     T_hover = W_TO
    
#     # Endurance Calculations
#     t_vert 	= Cap_batt / P_vert * 60						# [min] endurance time in pure vertical flight
#     t_cruise = Cap_batt  / P_cruise * 60					# [min] endurance time in pure horizontal flight
#     t_hover	= Cap_batt  / P_hover * 60						# [min] endurance time in pure hover flight
    
#     # Calculation for Max Thrust PREM VALUS
#     D_cruise_max	= 0.5 * rho_air * vel_cruise_max**2 * S_side_maxcruise * CD_side
#     theta_cruise_max = atan(D_cruise_max / W_TO)
#     T_cruise_max	= W_TO /cos(theta_cruise_max)
#     acc_max = 4 # [m/s2] - typical drone's acceleration is 4; max = 6
#     Thrust_acc_max = sqrt(T_hover**2 + M_TO**2* acc_max**2)
    
#     # Iteration Calculation
#     dt = 0.1
#     acc_to_cruise = np.arange(4,20+dt,dt)
#     Energy_data = []
#     t_2pts = []
#     for i in acc_to_cruise:
#         Thrust_acc_to_cruise = np.sqrt(T_hover**2 + (M_TO * i)**2) # Const Thrust level for the specified acceleration
#         t_acc_to_cruise = (vel_cruise - 0) / i # Simple a = (v - u) / t formula
#         dist_acc_to_cruise = 0.5 * i * t_acc_to_cruise**2 # Simple s = ut + 0.5at^2 formula
#         num_time_steps = int(round(t_acc_to_cruise / dt, 0))
#         vel_range = np.linspace(0, vel_cruise, num_time_steps) # array of vel with same as its corresponding time steps array
#         Power_acc_to_cruise = Thrust_acc_to_cruise * vel_range / eta_coaxial / eta_mech # Power = Thrust * velocity
#         Energy_acc_to_cruise = 0.5 * Power_acc_to_cruise[-1] * t_acc_to_cruise # area of triangle
        
#         dist_vel_cruise = dist_2pts - dist_acc_to_cruise * 2
#         t_vel_cruise = dist_vel_cruise / vel_cruise
#         Energy_vel_cruise = P_cruise * t_vel_cruise
        
#         Energy_total = (Energy_acc_to_cruise * 2 + Energy_vel_cruise  + (P_hover +P_PL) * hover_time) / 60**2
#         Energy_data.append(Energy_total)
        
#         t_total = t_acc_to_cruise * 2 + t_vel_cruise
#         t_2pts.append(t_total)
        
#     opt_hori_acc = acc_to_cruise[Energy_data.index(np.min(Energy_data))]
#     opt_hori_acc_energy = np.min(Energy_data) # [Wh]
#     min_energy_plot = np.ones((np.size(acc_to_cruise))) * opt_hori_acc_energy
#     opt_acc_plot = np.ones((np.size(Energy_data))) * opt_hori_acc
    
#     tot_energy_1measure = opt_hori_acc_energy
#     num_measurements = Cap_batt / tot_energy_1measure
#     total_flight_duration = num_measurements * (hover_time + t_total)
#     DOD_out.append(total_flight_duration)
    
# percent_dod = np.round(percent_dod, 2)
# percent_dod = percent_dod.tolist()
# index1_percent_dod =percent_dod.index(1.0)
# DOD_out = (np.array(DOD_out) - DOD_out[index1_percent_dod]) / DOD_out[index1_percent_dod] * 100
# percent_dod = (np.array(percent_dod)-1)*100


# texpsize= [26,28,30]
# SMALL_SIZE  = texpsize[0]
# MEDIUM_SIZE = texpsize[1]
# BIGGER_SIZE = texpsize[2]
# plt.rc('font', size=MEDIUM_SIZE)                    ## controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)                ## fontsize of the axes title
# plt.rc('axes', labelsize=SMALL_SIZE)                ## fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)               ## legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)             ## fontsize of the figure title
# matplotlib.rcParams['lines.linewidth']  = 1.5
# matplotlib.rcParams['figure.facecolor'] = 'white'
# matplotlib.rcParams['axes.facecolor']   = 'white'
# matplotlib.rcParams["legend.fancybox"]  = False

# fig, ax = plt.subplots(1,1,squeeze=False,figsize=(12,9))
# ax[0,0].plot(percent_dod, DOD_out, 'r-')
# # ax[0,0].plot(percent[100:300], struc_out[100:300], 'r-')
# ax[0,0].set_xlabel('DoD [%]')
# ax[0,0].set_ylabel('% change in Flight Duration')

# ax[0,0].axvline(0,color=(0,0,0),linewidth=1.3)
# ax[0,0].axhline(0,color=(0,0,0),linewidth=1.3)
# ax[0,0].minorticks_on() # set minor ticks
# ax[0,0].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
# ax[0,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid

# # plt.legend()
# fig.savefig("sensitivity_DOD.png", bbox_inches='tight')            

# plt.show()
    
    
    
    
    
    
    
    
    
# ## FOR SURFACE AREA CONTINGENCY, INSTEAD OF DOING SA * 1.1 * PERCENTAGE, SA * PERCENTAGE WAS DONE
    
# hover_time = 7
# g 			= 9.80665 # [m/s^2] Gravitational Acceleration
# rho 		= 1.225 # [kg/m^3] sea-level density
# mu 			= 1.81E-5 # [Ns/m^2] Dynamic visocity of air
# a 			= 343 # [m/s] Speed of sound of air
# D 			= 0.0254 * np.array([28, 30, 32]) # [m] Propeller Diameter
# eta_coaxial = 0.85 # [-] Coaxial efficiency loss
# eta_mech 	= 0.85 # [-] Total mechanical efficiency loss

# V_cruise	= 10 # [m/s] cruise speed
# V_max 		= 26 # [m/s] maximum speed
# CD 			= 1.05 # [-] assumed to be a cube
# S 			= 1.9596966763292372 # [m^2] frontal area

# D_cruise 	= 0.5 * rho * V_cruise ** 2 * CD * S # [N] drag force during cruise flight
# D_max 		= 0.5 * rho * V_max ** 2 * CD * S # [N] drag force during cruise flight

# # Mass budget
# n_propellers	= 8 # [-] number of propellers
# m_sensors 		= 1.124 # [kg] noise and air pollution sensors mass
# m_autonomy 		= 0.219 # [kg] autonomy hardware mass
# m_motors 		= n_propellers * 0.630 #[kg] brushless motor mass
# m_propellers 	= n_propellers * 91.46E-3 # [kg] propeller mass

# m_structures 	= (2.289 + 1.4) # [kg] operational empty mass


# # Frame data
# w_frame = 0.230 # [m] frame width
# h_frame = 0.260 # [m] frame height

# # Motor data
# V_motor = 36 # [V] nominal voltage of motor

# # Battery configuration - Chosen cell: EEMB LP623454
# DoD= 0.85 
    
# C_cell = 1200E-3 # [Ah] cell capacity
# m_cell = 22.5E-3 # [kg] cell mass
# V_cell = 3.7 # [V] nominal cell voltage

# h_cell = 69E-3 # [m] cell height
# w_cell = 34.5E-3 # [m] cell width 
# t_cell = 6.2E-3 # [m] cell thickness

# n_series = ceil(V_motor / V_cell)

# C_battery = C_cell # [Ah] battery capacity
# m_battery = n_series * m_cell # [kg] battery mass
# V_battery = n_series * V_cell # [V] battery voltage
# E_battery = V_battery * C_battery * DoD # [Wh] battery capacity

# A_frame = w_frame * h_frame
# A_battery = h_cell * n_series * t_cell
# A_battery = 1.05 * A_battery # contingency margin for cooling / packaging

# n_max = floor(A_frame / A_battery) - 1
# n_max = 16

# m_parallel = n_max * m_battery
# C_parallel = n_max * E_battery


# xls 	= pd.ExcelFile("PropellerData_Meizlik.xlsx")
# blade32 = pd.read_excel(xls, '32')

# m 	= m_sensors + m_autonomy + m_parallel + m_structures + m_motors + m_propellers

# theta_cruise 	= atan(D_cruise / (m * g)) # [rad] inclination angle during cruise
# theta_max 		= atan(D_max / (m * g)) # [rad] inclination angle during max speed

# # Define mission constraints
# T_hover 	= (m*g) / n_propellers # [N] Required thrust per propeller during hovering
# T_cruise 	= (m*g) / cos(theta_cruise) # [N] Required thrust during cruise

# P_mech 	= 2 * pi * blade32['Motor Speed [1/s]'] * blade32['Torque [Nm]'] / eta_coaxial

# blade32['Mechanical Power [W]'] = P_mech
# pfit_thrust = np.polyfit(blade32['Motor Speed [1/s]'], blade32['Thrust [N]'], deg=3)
# tl_thrust	= np.poly1d(pfit_thrust)

# pfit_power 	= np.polyfit(blade32['Motor Speed [1/s]'], blade32['Mechanical Power [W]'], deg=3)
# tl_power	= np.poly1d(pfit_power)

# pfit_torque	= np.polyfit(blade32['Motor Speed [1/s]'], blade32['Torque [Nm]'], deg=3)
# tl_torque	= np.poly1d(pfit_torque)

# """ Fixed Parameters """

# # General Parameters
# grav = 9.80665 # [m/s2]
# rho_air = 1.225 # [km/m3]

# # Fixed Performance Parameters
# vel_cruise = 10 # [m/s] - limited by the autonomy sensor
# vel_cruise_max = 25 # [m/s] - EASA's requirement
# vel_climb_max = 0.00508 * 500 # [m/s] - EASA's requirement
# N_prop = 8
# P_PL = 25 # [W] - sum of power req by payloads
# P_max = 360 * N_prop # [W] - max power output of chosen motor 

# # Fixed Operational Parameters
# dist_2pts = 150 # [m]
# dist_2pts_max = sqrt(2 * dist_2pts**2) # [m] - diagonal dist between 2 pts

# # Propeller Parameters
# D_prop = 0.8128
# r_prop = D_prop / 2
# t_prop = 0.04
# A_disk = N_prop * pi * r_prop**2

# """ VARIABLE Parameters - Should be adjusted accordingly with finalised values """

# # Efficiencies
# # eta_prop = 0.8 # included in EM
# eta_EM = 0.85 # Electromechanical Efficiency
# eta_coaxial = 0.85 # THIS SHOULD BE CORRECTED

# # Battery Parameters
# Cap_batt_Ah = 192 # [Ah]
# V_batt = 37 # [V]
# Cap_batt = Cap_batt_Ah * V_batt / 10 * DoD # [Wh] 
# M_Batt =  3.6 # [kg]

# # Motor Parameters
# N_motor = 8
# D_motor = 0.09
# t_motor = (27.4 + 12.5) / 1000 

# # Structural Parameters
# N_arm = 4
# l_arm = 0.564
# r_arm = 0.024
# D_arm = r_arm * 2
# l_body = 0.351
# w_body = 0.274
# h_body = 0.002

# M_TO = m
# W_TO = M_TO * grav # [N]
# CD_top = 1.05 # Assumed as a cube
# CD_side = 1.05 # Assumed as a cube


# for area_cont in percent_area:
#     S_top = (0.274 * 0.351) * area_cont  # [m2]
#     S_side_0deg = (N_motor * D_motor * t_motor + N_arm * D_arm * l_arm + h_body*w_body) * area_cont
    
#     """ Calculations for cruise power and tilt angle """
#     # Power Calculations
#     P_PL = 17.467 + 14.660
#     P_hover = sqrt((W_TO)**3 / (2 * rho_air * A_disk))  / eta_EM /eta_coaxial # [W] required power to hover WITHOUT PAYLOAD
#     P_vert = P_hover + CD_top * 0.5 * rho_air * vel_climb_max**2 * S_top * vel_climb_max / eta_EM / eta_coaxial
    
#     # Calculate tilt angle at cruise velocity
#     S_side_guess = 0 # [m2]
#     S_side_outcome = S_top * sin(radians(5)) + S_side_0deg / cos(radians(5))
#     while np.abs(S_side_guess - S_side_outcome) > 1E-5:
        
#         S_side_guess = S_side_outcome
#         D_cruise	= 0.5 * rho_air * vel_cruise**2 * S_side_guess * CD_side
#         theta       = atan(D_cruise / W_TO)
#         T_cruise	= W_TO /cos(theta)
#         S_side_outcome = S_top * sin(theta) +S_side_0deg / cos(theta)
        
#     S_side_cruise = S_side_outcome
#     S_side = S_side_cruise
    
#     rpm_dp2 = np.interp(T_cruise / n_propellers, blade32['Thrust [N]'], blade32['Motor Speed [1/s]'])
#     T_dp2 = tl_thrust(rpm_dp2)
#     P_dp2 = tl_power(rpm_dp2)
#     torque_dp2 = tl_torque(rpm_dp2)
#     P_cruise = (n_propellers * P_dp2 + D_cruise * V_cruise) / eta_mech 
        
#     # Calculte tilt angle at max cruise velocity
    
#     S_side_guess = 0 # [m2]
#     S_side_outcome2 = S_top * sin(radians(10)) + S_side_0deg / cos(radians(10))
#     while np.abs(S_side_guess - S_side_outcome2) > 1E-5:
        
#         S_side_guess = S_side_outcome2
        
#         D_cruise_max	= 0.5 * rho_air * vel_cruise_max**2 * S_side_guess * CD_side
#         theta_max = atan(D_cruise_max / W_TO)
#         T_cruise_max	= W_TO /cos(theta_max)
#         S_side_outcome = S_top * sin(theta_max) + S_side_0deg / cos(theta_max)
#         # print(S_side_outcome)
    
#     S_side_maxcruise = S_side_outcome
#     rpm_dp2 = np.interp(T_cruise_max / n_propellers, blade32['Thrust [N]'], blade32['Motor Speed [1/s]'])
#     P_cruise_max = (n_propellers * P_dp2 + D_cruise * V_cruise) / eta_mech 
#     T_hover = W_TO
    
#     # Endurance Calculations
#     t_vert 	= Cap_batt / P_vert * 60						# [min] endurance time in pure vertical flight
#     t_cruise = Cap_batt  / P_cruise * 60					# [min] endurance time in pure horizontal flight
#     t_hover	= Cap_batt  / P_hover * 60						# [min] endurance time in pure hover flight
    
#     # Calculation for Max Thrust PREM VALUS
#     D_cruise_max	= 0.5 * rho_air * vel_cruise_max**2 * S_side_maxcruise * CD_side
#     theta_cruise_max = atan(D_cruise_max / W_TO)
#     T_cruise_max	= W_TO /cos(theta_cruise_max)
#     acc_max = 4 # [m/s2] - typical drone's acceleration is 4; max = 6
#     Thrust_acc_max = sqrt(T_hover**2 + M_TO**2* acc_max**2)
    
#     # Iteration Calculation
#     dt = 0.1
#     acc_to_cruise = np.arange(4,20+dt,dt)
#     Energy_data = []
#     t_2pts = []
#     for i in acc_to_cruise:
#         Thrust_acc_to_cruise = np.sqrt(T_hover**2 + (M_TO * i)**2) # Const Thrust level for the specified acceleration
#         t_acc_to_cruise = (vel_cruise - 0) / i # Simple a = (v - u) / t formula
#         dist_acc_to_cruise = 0.5 * i * t_acc_to_cruise**2 # Simple s = ut + 0.5at^2 formula
#         num_time_steps = int(round(t_acc_to_cruise / dt, 0))
#         vel_range = np.linspace(0, vel_cruise, num_time_steps) # array of vel with same as its corresponding time steps array
#         Power_acc_to_cruise = Thrust_acc_to_cruise * vel_range / eta_coaxial / eta_mech # Power = Thrust * velocity
#         Energy_acc_to_cruise = 0.5 * Power_acc_to_cruise[-1] * t_acc_to_cruise # area of triangle
        
#         dist_vel_cruise = dist_2pts - dist_acc_to_cruise * 2
#         t_vel_cruise = dist_vel_cruise / vel_cruise
#         Energy_vel_cruise = P_cruise * t_vel_cruise
        
#         Energy_total = (Energy_acc_to_cruise * 2 + Energy_vel_cruise  + (P_hover +P_PL) * hover_time) / 60**2
#         Energy_data.append(Energy_total)
        
#         t_total = t_acc_to_cruise * 2 + t_vel_cruise
#         t_2pts.append(t_total)
        
#     opt_hori_acc = acc_to_cruise[Energy_data.index(np.min(Energy_data))]
#     opt_hori_acc_energy = np.min(Energy_data) # [Wh]
#     min_energy_plot = np.ones((np.size(acc_to_cruise))) * opt_hori_acc_energy
#     opt_acc_plot = np.ones((np.size(Energy_data))) * opt_hori_acc
    
#     tot_energy_1measure = opt_hori_acc_energy
#     num_measurements = Cap_batt / tot_energy_1measure
#     total_flight_duration = num_measurements * (hover_time + t_total)
#     SA_out.append(total_flight_duration)    


# percent_area = np.round(percent_area, 2)
# percent_area = percent_area.tolist()
# index1_percent_area = percent_area.index(1.0)
# SA_out = (np.array(SA_out) - SA_out[index1_percent_area]) / SA_out[index1_percent_area] * 100
# percent_area = (np.array(percent_area)-1) * 100

# texpsize= [26,28,30]
# SMALL_SIZE  = texpsize[0]
# MEDIUM_SIZE = texpsize[1]
# BIGGER_SIZE = texpsize[2]
# plt.rc('font', size=MEDIUM_SIZE)                    ## controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)                ## fontsize of the axes title
# plt.rc('axes', labelsize=SMALL_SIZE)                ## fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)               ## legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)             ## fontsize of the figure title
# matplotlib.rcParams['lines.linewidth']  = 1.5
# matplotlib.rcParams['figure.facecolor'] = 'white'
# matplotlib.rcParams['axes.facecolor']   = 'white'
# matplotlib.rcParams["legend.fancybox"]  = False

# fig, ax = plt.subplots(1,1,squeeze=False,figsize=(12,9))
# ax[0,0].plot(percent_area, SA_out, 'r-')
# # ax[0,0].plot(percent[100:300], struc_out[100:300], 'r-')
# ax[0,0].set_xlabel('Contingency added to Surface Area [%]')
# ax[0,0].set_ylabel('% change in Flight Duration')

# ax[0,0].axvline(0,color=(0,0,0),linewidth=1.3)
# ax[0,0].axhline(0,color=(0,0,0),linewidth=1.3)
# ax[0,0].minorticks_on() # set minor ticks
# ax[0,0].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
# ax[0,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid

# # plt.legend()
# fig.savefig("sensitivity_SAcontingency.png", bbox_inches='tight')            

# plt.show()
















hover_time_sens = 7 * percent_op
for hover_time in hover_time_sens:
    g 			= 9.80665 # [m/s^2] Gravitational Acceleration
    rho 		= 1.225 # [kg/m^3] sea-level density
    mu 			= 1.81E-5 # [Ns/m^2] Dynamic visocity of air
    a 			= 343 # [m/s] Speed of sound of air
    D 			= 0.0254 * np.array([28, 30, 32]) # [m] Propeller Diameter
    eta_coaxial = 0.85 # [-] Coaxial efficiency loss
    eta_mech 	= 0.85 # [-] Total mechanical efficiency loss
    
    V_cruise	= 10 # [m/s] cruise speed
    V_max 		= 26 # [m/s] maximum speed
    CD 			= 1.05 # [-] assumed to be a cube
    S 			= 1.9596966763292372 # [m^2] frontal area
    
    D_cruise 	= 0.5 * rho * V_cruise ** 2 * CD * S # [N] drag force during cruise flight
    D_max 		= 0.5 * rho * V_max ** 2 * CD * S # [N] drag force during cruise flight
    
    # Mass budget
    n_propellers	= 8 # [-] number of propellers
    m_sensors 		= 1.124 # [kg] noise and air pollution sensors mass
    m_autonomy 		= 0.219 # [kg] autonomy hardware mass
    m_motors 		= n_propellers * 0.630 #[kg] brushless motor mass
    m_propellers 	= n_propellers * 91.46E-3 # [kg] propeller mass
    
    m_structures 	= (2.289 + 1.4) # [kg] operational empty mass
    
    
    # Frame data
    w_frame = 0.230 # [m] frame width
    h_frame = 0.260 # [m] frame height
    
    # Motor data
    V_motor = 36 # [V] nominal voltage of motor
    
    # Battery configuration - Chosen cell: EEMB LP623454
    DoD= 0.85 
        
    C_cell = 1200E-3 # [Ah] cell capacity
    m_cell = 22.5E-3 # [kg] cell mass
    V_cell = 3.7 # [V] nominal cell voltage
    
    h_cell = 69E-3 # [m] cell height
    w_cell = 34.5E-3 # [m] cell width 
    t_cell = 6.2E-3 # [m] cell thickness
    
    n_series = ceil(V_motor / V_cell)
    
    C_battery = C_cell # [Ah] battery capacity
    m_battery = n_series * m_cell # [kg] battery mass
    V_battery = n_series * V_cell # [V] battery voltage
    E_battery = V_battery * C_battery * DoD # [Wh] battery capacity
    
    A_frame = w_frame * h_frame
    A_battery = h_cell * n_series * t_cell
    A_battery = 1.05 * A_battery # contingency margin for cooling / packaging
    
    n_max = floor(A_frame / A_battery) - 1
    n_max = 16
    
    m_parallel = n_max * m_battery
    C_parallel = n_max * E_battery
    
    
    xls 	= pd.ExcelFile("PropellerData_Meizlik.xlsx")
    blade32 = pd.read_excel(xls, '32')
    
    m 	= m_sensors + m_autonomy + m_parallel + m_structures + m_motors + m_propellers
    
    theta_cruise 	= atan(D_cruise / (m * g)) # [rad] inclination angle during cruise
    theta_max 		= atan(D_max / (m * g)) # [rad] inclination angle during max speed
    
    # Define mission constraints
    T_hover 	= (m*g) / n_propellers # [N] Required thrust per propeller during hovering
    T_cruise 	= (m*g) / cos(theta_cruise) # [N] Required thrust during cruise
    
    P_mech 	= 2 * pi * blade32['Motor Speed [1/s]'] * blade32['Torque [Nm]'] / eta_coaxial
    
    blade32['Mechanical Power [W]'] = P_mech
    pfit_thrust = np.polyfit(blade32['Motor Speed [1/s]'], blade32['Thrust [N]'], deg=3)
    tl_thrust	= np.poly1d(pfit_thrust)
    
    pfit_power 	= np.polyfit(blade32['Motor Speed [1/s]'], blade32['Mechanical Power [W]'], deg=3)
    tl_power	= np.poly1d(pfit_power)
    
    pfit_torque	= np.polyfit(blade32['Motor Speed [1/s]'], blade32['Torque [Nm]'], deg=3)
    tl_torque	= np.poly1d(pfit_torque)
    
    """ Fixed Parameters """
    
    # General Parameters
    grav = 9.80665 # [m/s2]
    rho_air = 1.225 # [km/m3]
    
    # Fixed Performance Parameters
    vel_cruise = 10 # [m/s] - limited by the autonomy sensor
    vel_cruise_max = 25 # [m/s] - EASA's requirement
    vel_climb_max = 0.00508 * 500 # [m/s] - EASA's requirement
    N_prop = 8
    P_PL = 25 # [W] - sum of power req by payloads
    P_max = 360 * N_prop # [W] - max power output of chosen motor 
    
    # Fixed Operational Parameters
    dist_2pts = 150 # [m]
    dist_2pts_max = sqrt(2 * dist_2pts**2) # [m] - diagonal dist between 2 pts
    
    # Propeller Parameters
    D_prop = 0.8128
    r_prop = D_prop / 2
    t_prop = 0.04
    A_disk = N_prop * pi * r_prop**2
    
    """ VARIABLE Parameters - Should be adjusted accordingly with finalised values """
    
    # Efficiencies
    # eta_prop = 0.8 # included in EM
    eta_EM = 0.85 # Electromechanical Efficiency
    eta_coaxial = 0.85 # THIS SHOULD BE CORRECTED
    
    # Battery Parameters
    Cap_batt_Ah = 192 # [Ah]
    V_batt = 37 # [V]
    Cap_batt = Cap_batt_Ah * V_batt / 10 * DoD # [Wh] 
    M_Batt =  3.6 # [kg]
    
    # Motor Parameters
    N_motor = 8
    D_motor = 0.09
    t_motor = (27.4 + 12.5) / 1000 
    
    # Structural Parameters
    N_arm = 4
    l_arm = 0.564
    r_arm = 0.024
    D_arm = r_arm * 2
    l_body = 0.351
    w_body = 0.274
    h_body = 0.002
    
    M_TO = m
    W_TO = M_TO * grav # [N]
    CD_top = 1.05 # Assumed as a cube
    CD_side = 1.05 # Assumed as a cube


    area_cont = 1.1
    S_top = (0.274 * 0.351) * area_cont  # [m2]
    S_side_0deg = (N_motor * D_motor * t_motor + N_arm * D_arm * l_arm + h_body*w_body) * area_cont
    
    """ Calculations for cruise power and tilt angle """
    # Power Calculations
    P_PL = 17.467 + 14.660
    P_hover = sqrt((W_TO)**3 / (2 * rho_air * A_disk))  / eta_EM /eta_coaxial # [W] required power to hover WITHOUT PAYLOAD
    P_vert = P_hover + CD_top * 0.5 * rho_air * vel_climb_max**2 * S_top * vel_climb_max / eta_EM / eta_coaxial
    
    # Calculate tilt angle at cruise velocity
    S_side_guess = 0 # [m2]
    S_side_outcome = S_top * sin(radians(5)) + S_side_0deg / cos(radians(5))
    while np.abs(S_side_guess - S_side_outcome) > 1E-5:
        
        S_side_guess = S_side_outcome
        D_cruise	= 0.5 * rho_air * vel_cruise**2 * S_side_guess * CD_side
        theta       = atan(D_cruise / W_TO)
        T_cruise	= W_TO /cos(theta)
        S_side_outcome = S_top * sin(theta) +S_side_0deg / cos(theta)
        
    S_side_cruise = S_side_outcome
    S_side = S_side_cruise
    
    rpm_dp2 = np.interp(T_cruise / n_propellers, blade32['Thrust [N]'], blade32['Motor Speed [1/s]'])
    T_dp2 = tl_thrust(rpm_dp2)
    P_dp2 = tl_power(rpm_dp2)
    torque_dp2 = tl_torque(rpm_dp2)
    P_cruise = (n_propellers * P_dp2 + D_cruise * V_cruise) / eta_mech 
        
    # Calculte tilt angle at max cruise velocity
    
    S_side_guess = 0 # [m2]
    S_side_outcome2 = S_top * sin(radians(10)) + S_side_0deg / cos(radians(10))
    while np.abs(S_side_guess - S_side_outcome2) > 1E-5:
        
        S_side_guess = S_side_outcome2
        
        D_cruise_max	= 0.5 * rho_air * vel_cruise_max**2 * S_side_guess * CD_side
        theta_max = atan(D_cruise_max / W_TO)
        T_cruise_max	= W_TO /cos(theta_max)
        S_side_outcome = S_top * sin(theta_max) + S_side_0deg / cos(theta_max)
        # print(S_side_outcome)
    
    S_side_maxcruise = S_side_outcome
    rpm_dp2 = np.interp(T_cruise_max / n_propellers, blade32['Thrust [N]'], blade32['Motor Speed [1/s]'])
    P_cruise_max = (n_propellers * P_dp2 + D_cruise * V_cruise) / eta_mech 
    T_hover = W_TO
    
    # Endurance Calculations
    t_vert 	= Cap_batt / P_vert * 60						# [min] endurance time in pure vertical flight
    t_cruise = Cap_batt  / P_cruise * 60					# [min] endurance time in pure horizontal flight
    t_hover	= Cap_batt  / P_hover * 60						# [min] endurance time in pure hover flight
    
    # Calculation for Max Thrust PREM VALUS
    D_cruise_max	= 0.5 * rho_air * vel_cruise_max**2 * S_side_maxcruise * CD_side
    theta_cruise_max = atan(D_cruise_max / W_TO)
    T_cruise_max	= W_TO /cos(theta_cruise_max)
    acc_max = 4 # [m/s2] - typical drone's acceleration is 4; max = 6
    Thrust_acc_max = sqrt(T_hover**2 + M_TO**2* acc_max**2)
    
    # Iteration Calculation
    dt = 0.1
    acc_to_cruise = np.arange(4,20+dt,dt)
    Energy_data = []
    t_2pts = []
    for i in acc_to_cruise:
        Thrust_acc_to_cruise = np.sqrt(T_hover**2 + (M_TO * i)**2) # Const Thrust level for the specified acceleration
        t_acc_to_cruise = (vel_cruise - 0) / i # Simple a = (v - u) / t formula
        dist_acc_to_cruise = 0.5 * i * t_acc_to_cruise**2 # Simple s = ut + 0.5at^2 formula
        num_time_steps = int(round(t_acc_to_cruise / dt, 0))
        vel_range = np.linspace(0, vel_cruise, num_time_steps) # array of vel with same as its corresponding time steps array
        Power_acc_to_cruise = Thrust_acc_to_cruise * vel_range / eta_coaxial / eta_mech # Power = Thrust * velocity
        Energy_acc_to_cruise = 0.5 * Power_acc_to_cruise[-1] * t_acc_to_cruise # area of triangle
        
        dist_vel_cruise = dist_2pts - dist_acc_to_cruise * 2
        t_vel_cruise = dist_vel_cruise / vel_cruise
        Energy_vel_cruise = P_cruise * t_vel_cruise
        
        Energy_total = (Energy_acc_to_cruise * 2 + Energy_vel_cruise  + (P_hover +P_PL) * hover_time) / 60**2
        Energy_data.append(Energy_total)
        
        t_total = t_acc_to_cruise * 2 + t_vel_cruise
        t_2pts.append(t_total)
        
    opt_hori_acc = acc_to_cruise[Energy_data.index(np.min(Energy_data))]
    opt_hori_acc_energy = np.min(Energy_data) # [Wh]
    min_energy_plot = np.ones((np.size(acc_to_cruise))) * opt_hori_acc_energy
    opt_acc_plot = np.ones((np.size(Energy_data))) * opt_hori_acc
    
    tot_energy_1measure = opt_hori_acc_energy
    num_measurements = Cap_batt / tot_energy_1measure
    total_flight_duration = num_measurements * (hover_time + t_total)
    op_out.append(total_flight_duration)    
    


percent_op = np.round(percent_op, 2)
percent_op = percent_op.tolist()
index1_percent_op = percent_op.index(1.0)
op_out = (np.array(op_out) - op_out[index1_percent_op]) / op_out[index1_percent_op] * 100
percent_op = (np.array(percent_op) - 1) * 100

texpsize= [26,28,30]
SMALL_SIZE  = texpsize[0]
MEDIUM_SIZE = texpsize[1]
BIGGER_SIZE = texpsize[2]
plt.rc('font', size=MEDIUM_SIZE)                    ## controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)                ## fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)                ## fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)               ## legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)             ## fontsize of the figure title
matplotlib.rcParams['lines.linewidth']  = 1.5
matplotlib.rcParams['figure.facecolor'] = 'white'
matplotlib.rcParams['axes.facecolor']   = 'white'
matplotlib.rcParams["legend.fancybox"]  = False

fig, ax = plt.subplots(1,1,squeeze=False,figsize=(12,9))
ax[0,0].plot(percent_op, op_out, 'r-')
# ax[0,0].plot(percent[100:300], struc_out[100:300], 'r-')
ax[0,0].set_xlabel('% change in Hover Time')
ax[0,0].set_ylabel('% change in Flight Duration')

ax[0,0].axvline(0,color=(0,0,0),linewidth=1.3)
ax[0,0].axhline(0,color=(0,0,0),linewidth=1.3)
ax[0,0].minorticks_on() # set minor ticks
ax[0,0].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
ax[0,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid

# plt.legend()
fig.savefig("sensitivity_hovertime.png", bbox_inches='tight')            

plt.show()











# plt.figure()
# plt.plot(percent, struc_out, label='structure')
# plt.plot(percent_dod, DOD_out, label='DoD')
# plt.plot(percent_area, SA_out, label='SA')
# plt.plot(percent_op, op_out, label = 'hover time')
# plt.legend()
# plt.show()
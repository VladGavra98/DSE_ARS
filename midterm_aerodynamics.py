'''
	AE3200 Design Synthesis Exercise
	Group 09 - Autonomous Environmental Sensing

	Authors: Park ChangKyu, Widmann Sebastian
	Created: 12.05.2020
	
	Project Supervisors:
		- Dr. Irene C. Dedoussi
		- Dr. Ir. Mirjam Snellen
		- Ir. Lorenzo Pasqualetto Cassinis
		- Mark Schelbergen

	This script is used for the initial sizing of the propellers / aerodynamic surfaces for the concept selection		
'''
from math import pi, sqrt
import numpy as np

# Common Parameters
grav = 9.8065 # [m/s2]
op_alt = 50 # [m]
safety_factor = 1.5
horz_acc_req = 1 * grav # [m/s2]
vert_acc_req = 1 * grav 
rho_air = 1.225 # [kg/m3]
dim = 2
op_alt = 50

# Masses [kg]
noise_sensor_mass = 0.1
airpol_sensor_mass = 10
Payload_mass = noise_sensor_mass + airpol_sensor_mass
Battery_mass = 20 


# Propeller function
def prop_power(N, rho_air, area, rotor_vel):
    '''
    Parameters
    ----------
    N : Number of propellers
    rho_air : Real time air density
    area : Area swept by propellers
    rotor_vel : Instantaneous peripheral velocity of rotors    

    Returns
    -------
    Thrust : Total thrust produced by all propellers

    '''
    Thrust = N * rho_air * area * rotor_vel**2
    return Thrust

########### VTOL: critical Motion = vertical
VTOL_prop_no = 3
VTOL_struc_mass = 50
VTOL_tot_mass = VTOL_struc_mass + Payload_mass + Battery_mass
VTOL_top_area = 1**2


Cd_VTOL_top = 4 # TBU
max_vert_vel = sqrt(0 + 2 * vert_acc_req * op_alt)
max_thrust_req_VTOL = VTOL_tot_mass * (grav + vert_acc_req) + Cd_VTOL_top * 0.5 * rho_air * max_vert_vel**2 * VTOL_top_area # at alt = op_alt
vel_thru_disk = 50 # TBU
max_power_req_VTOL = max_thrust_req_VTOL * (max_vert_vel + vel_thru_disk)

########### Drone: critical motion = horizontal 
Drone_prop_no = 4
Drone_struc_mass = 50
Drone_tot_mass = Drone_struc_mass + Payload_mass + Battery_mass


Cd_Drone_top = 2.56 # TBU
Drone_top_area = 1**2 # TBU
max_thrust_req_Drone_vert = Drone_tot_mass * (grav + vert_acc_req) + Cd_Drone_top * 0.5 * rho_air * max_vert_vel**2 * Drone_top_area
max_power_req_Drone_vert = max_thrust_req_Drone_vert * (max_vert_vel + vel_thru_disk)

Cd_Drone_side = (0.0464 + 0.0895) / 2 # side
Drone_side_area = 1**2 # TBU
max_horz_vel = sqrt(0 + 2 * horz_acc_req * 50)
inc_angle = np.radians(35)  # rad
max_thrust_req_Drone_horz =  Drone_tot_mass * horz_acc_req + Cd_Drone_side * 0.5 * rho_air * max_horz_vel**2 * Drone_side_area
max_power_req_Drone_horz = max_thrust_req_Drone_horz  * (max_horz_vel + vel_thru_disk) / np.sin(inc_angle)


############ Airship: critical motion = horizontal
Airship_prop_no = 2
Airship_struc_mass = 200
Airship_tot_mass = Airship_struc_mass + Payload_mass + Battery_mass
# Vol_airship = 1000
Cd_Airship = 0.1 # TBU
Airship_frontarea = 1.5**2

max_thrust_req_Airship_horz =  Airship_tot_mass * horz_acc_req + Cd_Airship * 0.5 * rho_air * max_horz_vel**2 * Airship_frontarea

max_power_req_Airship_horz = max_thrust_req_Airship_horz  * (max_horz_vel + vel_thru_disk)

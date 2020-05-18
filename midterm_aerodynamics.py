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
horz_acc_req = 5
vert_acc_req = 0.15
rho_air = 1.225 # [kg/m3]
dim = 2
op_alt = 30
volt = 2.1 # [V]
Battery_cap = 160 * 60 * 60 # [Wh] 2 cells = 76Ah
trip_dist = 1500 # [m]
max_vert_vel = sqrt(0 + 2 * vert_acc_req * op_alt)
max_horz_vel = 25

# Vi calculations
ang_vel_disk = 660 # [rad/s]
radius_disk = 23e-3 # [m]
vel_thru_disk = 50 # TBU

# Masses [kg]
noise_sensor_mass = 0.1
airpol_sensor_mass = 5
Payload_mass = noise_sensor_mass + airpol_sensor_mass
Battery_mass = 5 


# Propeller function
def vel_thru_disk(angular_vel_disk, radius_disk, thrust, rho_air):
    '''
    Parameters
    ----------
    angular_vel_disk : angular velocity of the disk(propeller) [rad/s]
    radius_disk : radius of the disk [m]
    thrust : thrust [N]
    rho_air : air density [kg/m3]

    Returns
    -------
    vi : velocity of air through the disk [m/s]

    '''
    import numpy as np
    area_disk = np.pi * radius_disk**2
    thrust_coeff = thrust * 2 / (radius_disk * angular_vel_disk)**2 / area_disk
    vi = 0.5 * np.sqrt(thrust_coeff) * radius_disk * angular_vel_disk
    
    return vi


########### VTOL: critical Motion = vertical
VTOL_prop_no = 3
VTOL_struc_mass = 10
VTOL_tot_mass = VTOL_struc_mass + Payload_mass + Battery_mass
VTOL_top_area = 0.2
Cd_VTOL_top = 1.28 
VTOL_side_area = 0.05
Cd_VTOL_side = 0.05

max_thrust_req_VTOL_vert = VTOL_tot_mass * (grav + vert_acc_req) + Cd_VTOL_top * 0.5 * rho_air * max_vert_vel**2 * VTOL_top_area # at alt = op_alt
max_power_req_VTOL_vert = max_thrust_req_VTOL_vert * (max_vert_vel + vel_thru_disk(ang_vel_disk, radius_disk, max_thrust_req_VTOL_vert, rho_air))


max_thrust_req_VTOL_horz =  VTOL_tot_mass * horz_acc_req + Cd_VTOL_side * 0.5 * rho_air * max_horz_vel**2 * VTOL_side_area
max_power_req_VTOL_horz = max_thrust_req_VTOL_horz  * (max_horz_vel + vel_thru_disk(ang_vel_disk, radius_disk, max_thrust_req_VTOL_horz, rho_air)) 

Energy_per_trip_VTOL = max_power_req_VTOL_vert * 20 # the time req to reach operational height is 5seconds
Num_of_trips_VTOL = Battery_cap / Energy_per_trip_VTOL

########### Drone: critical motion = horizontal 
Drone_prop_no = 4
Drone_struc_mass = 50
Drone_tot_mass = Drone_struc_mass + Payload_mass + Battery_mass
Cd_Drone_top = 1.28 
Drone_top_area = 1**2 # TBU

max_thrust_req_Drone_vert = Drone_tot_mass * (grav + vert_acc_req) + Cd_Drone_top * 0.5 * rho_air * max_vert_vel**2 * Drone_top_area
max_power_req_Drone_vert = max_thrust_req_Drone_vert * (max_vert_vel + vel_thru_disk(ang_vel_disk, radius_disk, max_thrust_req_Drone_vert, rho_air))

Cd_Drone_side = (0.0464 + 0.0895) / 2 # side
Drone_side_area = 1**2 # TBU
max_horz_vel = 25
inc_angle = np.radians(35)  # rad

max_thrust_req_Drone_horz =  Drone_tot_mass * horz_acc_req + Cd_Drone_side * 0.5 * rho_air * max_horz_vel**2 * Drone_side_area
max_power_req_Drone_horz = max_thrust_req_Drone_horz  * (max_horz_vel + vel_thru_disk(ang_vel_disk, radius_disk, max_thrust_req_Drone_horz, rho_air)) / np.sin(inc_angle)

Energy_per_trip_Drone = max_power_req_VTOL_vert * 20 + max_power_req_Drone_horz * 5 # the time req to reach operational height is 5seconds, 3 seconds to reach max top horz speed
Num_of_trips_Drone = Battery_cap / Energy_per_trip_Drone

############ Airship: critical motion = horizontal
Airship_prop_no = 2
Airship_struc_mass = 200
Airship_tot_mass = Airship_struc_mass + Payload_mass + Battery_mass
# Vol_airship = 1000
Cd_Airship = 0.1 # TBU
Airship_frontarea = 1.5**2

max_thrust_req_Airship_horz =  Airship_tot_mass * horz_acc_req + Cd_Airship * 0.5 * rho_air * max_horz_vel**2 * Airship_frontarea
max_power_req_Airship_horz = max_thrust_req_Airship_horz  * (max_horz_vel + vel_thru_disk(ang_vel_disk, radius_disk, max_thrust_req_Airship_horz, rho_air))

Energy_per_trip_Airship = max_power_req_Airship_horz * 5
Num_of_trips_Airhsip = Battery_cap / Energy_per_trip_Airship

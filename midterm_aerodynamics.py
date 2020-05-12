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
from math import pi
# Other Parameters
grav = 9.8065 # [m/s2]
op_alt = 50 # [m]
safety_factor = 1.5
dim = 2 [m]

# Masses [kg]
Payload_mass = 50   
Battery_mass = 20 
Structural_mass = 20
Total_mass = Payload_mass + Battery_mass + Structural_mass


# Motion - vertical (VTOL)
vert_acc_needed = 0.5 # [m/s2]
Lift_req_vert = Total_mass *  (vert_acc_needed + grav) # [N]

# Motion - horizontal (Airship, Drone)
hor_acc_needed = 50 # [m/s2]

# Airship calculations
a_hull = 0.03
rho_He = 0.164 # [km/m3]
Vol_airship = 4/3 * pi * (dim/2)**3
Mass_airship = a_hull * rho_He * Vol_airship

def 
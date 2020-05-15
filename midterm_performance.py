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

# Import modules
from math import sqrt, pi, radians, cos, sin

# ---------------------------------------------------------------------------------------------
# VTOL Flying Wing
# ---------------------------------------------------------------------------------------------
# Define constants
SF 		= 1.5  								# [-] safety factor for uncertainty
g 		= 9.81 								# [m/s^2] gravitational acceleration
rho 	= 1.225 							# [kg/m^3] sea-level density
D_prop 	= 0.152 							# [m] propeller diameter - equivalent to 6 inch propeller
A 		= 3 * (pi/4 * D_prop**2)			# [m^2] total propeller disk area - minimum of 3 propellers required for static stability
eta 	= 0.85 								# [-] preliminary estimation for propulsion efficiency
sigma 	= 1.2 								# [-] expansion ratio of duct

# Battery assumptions for Li-S technology
costpb 	= 200 								# [$/kWh] estimated battery costs per capacity
C_cell 	= 19 								# [Ah] cell capacity
V_cell 	= 2.1 								# [V] nominal cell voltage
n 		= 4									# [-] number of cells in series
C_batt  = n * C_cell 						# [Ah] battery capacity
V_batt 	= n * V_cell 						# [V] battery voltage

C_batt 	= C_batt * V_batt 					# [Wh] battery capacity converted from Ah to Wh
price 	= C_batt / (costpb / 1000) 			# [$] estimated battery costs

# Mass breakdown
m_PL 	= 3.7 								# [kg] payload mass
m_batt 	= 1.0								# [kg] battery mass
m_stru 	= 1.2 * 3.802						# [kg] structural mass - 3.802kg without contingency

m 		= m_PL + m_batt + m_stru 			# [kg] maximum take-off mass
m 		= 1.2 * m 							# [kg] add contingency of 20%

# Hover calculations
T_hover = m * g 	 						# [N] required thrust to hover
P_hover = sqrt(T_hover ** 3 / (4*sigma*rho*A))	# [W] required power to hover
P_hover = P_hover / eta * SF				# [W] corrected power requirement with safety factor and propulsive efficiency
t_hover = C_batt / P_hover * 60				# [s] maximum time hovering

# Forward flight caluclations
V 		= 30 								# [m/s] estimated maximum speed
S 		= 0.5 								# [m^2] estimated surface area based on reference aircraft (https://scholars.direct/Articles/aerospace-engineering-and-mechanics/jaem-4-020.php?jid=aerospace-engineering-and-mechanics)
b 		= 2									# [m] maximum wingspan
MAC 	= 0.2 								# [m] mean aerodynamic chord based on reference aircraft
A 		= b **2 / S 						# [-] aspect ratio
mu 		= 1.81e-5							# [kg/m*s] dynamic visocity of air
Re 		= rho * V * MAC / mu 				# [-] estimated Reynolds number at max. speed 
CD0 	= 4.98 / sqrt(Re) 					# [-] estimated parasite drag coefficient (http://apmonitor.com/me575/uploads/Main/2013_Flying_Wing.pdf)
e 		= 0.7 								# [-] estimated efficiency factor (https://www.grc.nasa.gov/www/k-12/airplane/induced.html)

P_req 	= 0.5 * CD0 * rho * V **3 * S + 2 * m * g / (pi * e * A * rho * V * S)
P_req 	= SF * P_req						# [W] required power with safety factor added
t_req	= C_batt / P_req * 60				# [s] maximum time flying

# Print important variables
print('VTOL Flying Wing Caluclations')
print('----------------------------------------------')
print('Cost of Li-S battery =', round(price,1), '$')
print('Required hovering power =', round(P_hover,1), 'W')
print('Maximum estimated hover time =', round(t_hover,1), 'min')
print('Required flight power =', round(P_req,1), 'W')
print('Maximum estimated flight time =', round(t_req,1), 'min')


# ---------------------------------------------------------------------------------------------
# Drone 
# ---------------------------------------------------------------------------------------------
# Define constants
# D_prop 	= 0.406 							# [m] propeller diameter - equivalent to 16 inch propeller
A 		= 8 * (pi/4 * D_prop**2)			# [m^2] total propeller disk area - minimum of 8 propellers required for static stability

# Mass breakdown
m_PL 	= 3.84 								# [kg] payload mass
m_batt 	= 1.0								# [kg] battery mass
m_stru 	= 1.2 * 1.423						# [kg] structural mass - 1.423kg without contingency

m 		= m_PL + m_batt + m_stru 			# [kg] maximum take-off mass
m 		= 1.2 * m 							# [kg] add contingency of 20%

# Hover calculations
T_hover = m * g 	 						# [N] required thrust to hover
P_hover = sqrt(T_hover ** 3 / (2*rho*A))	# [W] required power to hover
P_hover = P_hover / eta * SF				# [W] corrected power requirement with safety factor and propulsive efficiency
t_hover = C_batt / P_hover * 60				# [s] maximum time hovering

# Forward flight caluclations
CD 		= 0.15				 				# [-] estimated drag coefficient at 30deg inclincation (https://aip.scitation.org/doi/pdf/10.1063/1.4981989)
i 		= radians(30)						# [rad] inclincation angle (https://aip.scitation.org/doi/pdf/10.1063/1.4981989)
S 		= 1.325 * sin(i) 					# [m^2] frontal area based on reference values (https://freeflysystems.com/alta-8/specs)

T_req 	= T_hover / cos(i)					# [N] required additional thrust for horizontal flight due to inclincation
D 		= 0.5 * rho * V**2 * S * CD 		# [N] experienced drag during horizontal flight
T_tot 	= T_req + D 						# [N] total thrust required

P_req 	= sqrt(T_tot ** 3 / (2*rho*A))		# [W] required power during horizontal flight
P_req 	= SF * P_req						# [W] required power during horizontal flight with safety factor
t_req	= C_batt / P_req * 60				# [s] maximum time flying

# Print important variables
print()
print()
print('Drone Caluclations')
print('----------------------------------------------')
print('Cost of Li-S battery =', round(price,1), '$')
print('Required hovering power =', round(P_hover,1), 'W')
print('Maximum estimated hover time =', round(t_hover,1), 'min')
print('Required flight power =', round(P_req,1), 'W')
print('Maximum estimated flight time =', round(t_req,1), 'min')
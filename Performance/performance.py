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

	This script is used for the propeller and motor sizing of the coaxial quadcopter		
'''

# Import modules
from math import pi, sqrt, atan, cos, sin, ceil, floor
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, fsolve
from sklearn.metrics import r2_score
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
from PIL import Image
import cv2
from itertools import zip_longest
# ==============================================================================================
# Set global plotting parameters
# ==============================================================================================
texpsize= [18,20,22]

plt.rc('font', size=texpsize[1], family='serif')                        # controls default text sizes
plt.rc('axes', titlesize=texpsize[1])                                   # fontsize of the axes title
plt.rc('axes', labelsize=texpsize[1])                                   # fontsize of the x and y labels
plt.rc('xtick', labelsize=texpsize[0])                                  # fontsize of the tick labels
plt.rc('ytick', labelsize=texpsize[0])                                  # fontsize of the tick labels
plt.rc('legend', fontsize=texpsize[0])                                  # legend fontsize
plt.rc('figure', titlesize=texpsize[2])                                 # fontsize of the figure title
matplotlib.rcParams['lines.linewidth']  = 1.5
matplotlib.rcParams['figure.facecolor'] = 'white'
matplotlib.rcParams['axes.facecolor']   = 'white'
matplotlib.rcParams["legend.fancybox"]  = False

# Define constants
g 			= 9.80665 # [m/s^2] Gravitational Acceleration
rho 		= 1.225 # [kg/m^3] sea-level density
mu 			= 1.81E-5 # [Ns/m^2] Dynamic visocity of air
a 			= 343 # [m/s] Speed of sound of air
D 			= 0.0254 * np.array([28, 30, 32]) # [m] Propeller Diameter
eta_coaxial = 0.85 # [-] Coaxial efficiency loss
eta_mech 	= 0.8 # [-] Total mechanical efficiency loss

V_cruise	= 10 # [m/s] cruise speed
V_max 		= 26 # [m/s] maximum speed
CD 			= 1.05 # [-] assumed to be a cube
S 			= 0.3 # [m^2] frontal area

D_cruise 	= 0.5 * rho * V_cruise ** 2 * CD * S # [N] drag force during cruise flight
D_max 		= 0.5 * rho * V_max ** 2 * CD * S # [N] drag force during cruise flight

# Mass budget
n_propellers	= 8 # [-] number of propellers
m_sensors 		= 2.386 # [kg] noise and air pollution sensors mass
m_autonomy 		= 0.219 # [kg] autonomy hardware mass
m_motors 		= n_propellers * 0.630 #[kg] brushless motor mass
m_propellers 	= n_propellers * 91.46E-3 # [kg] propeller mass
m_structures 	= 2.289 # [kg] operational empty mass

# Frame data
w_frame = 0.236 # [m] frame width
h_frame = 0.313 # [m] frame height

# Motor data
V_motor = 36 # [V] nominal voltage of motor

# Battery configuration - Chosen cell: EEMB LP623454
DoD = 0.8 # [-] Depth of Discharge
C_cell = 1200E-3 # [Ah] cell capacity
m_cell = 22.5E-3 # [kg] cell mass
V_cell = 3.7 # [V] nominal cell voltage

h_cell = 69E-3 # [m] cell height
w_cell = 34.5E-3 # [m] cell width 
t_cell = 6.2E-3 # [m] cell thickness

n_series = ceil(V_motor / V_cell)

C_battery = n_series * C_cell # [Ah] battery capacity
m_battery = n_series * m_cell # [kg] battery mass
V_battery = n_series * V_cell # [V] battery voltage
E_battery = V_battery * C_battery # [Wh] battery capacity

A_frame = w_frame * h_frame
A_battery = h_cell * n_series * t_cell
A_battery = 1.05 * A_battery # contingency margin for cooling / packaging

n_max = floor(A_frame / A_battery) - 1

# Import propeller data from excel file in the following format
# ['Motor Speed [1/s', 'Thrust [N]', 'Torque [Nm]', 'Voltage[V]', 'Current [A]']
xls 	= pd.ExcelFile("PropellerData_Meizlik.xlsx")
blade28 = pd.read_excel(xls, '28')
blade30 = pd.read_excel(xls, '30')
blade32 = pd.read_excel(xls, '32')

# Store blade dataframes into list for easy iterating
blades 	= [blade28, blade30, blade32]
bladetype = ['28" x 9.4"', '30" x 10"', '32" x 10.6"']
exportname = ['28', '30', '32']

# Initialise figures
fig1, ax1 = plt.subplots(nrows=2, ncols=2, squeeze=False, figsize=(20,20))

# Colors used for iterative plotting
colors = ['tab:blue', 'tab:orange', 'tab:grey']

# ==============================================================================================
# ITERATE OVER PROPELLER DIAMETERS
# ==============================================================================================
for index, df in enumerate(blades):
	m_parallel = n_max * m_battery 

	# Mass budget
	m 				= m_sensors + m_autonomy + m_parallel + m_structures + m_motors + m_propellers
	m 				= 1.05 * m # [kg] contingency added

	theta_cruise 	= atan(D_cruise / (m * g)) # [rad] inclination angle during cruise
	theta_max 		= atan(D_max / (m * g)) # [rad] inclination angle during max speed

	# Define mission constraints
	T_hover 	= (m*g) / n_propellers # [N] Required thrust per propeller during hovering
	T_cruise 	= (m*g) / cos(theta_cruise) # [N] Required thrust during cruise
	T_max 		= (m*g) / cos(theta_max) # [N] Required thrust during max speed

	P_mech 	= 2 * pi * df['Motor Speed [1/s]'] * df['Torque [Nm]'] / eta_coaxial
	P_elec 	= df['Voltage [V]'] * df['Current [A]']
	Ct 		= df['Thrust [N]'] / (rho * df['Motor Speed [1/s]']**2 * D[index]**4)
	Cp 		= P_mech / (rho * df['Motor Speed [1/s]']**3 * D[index]**5)
	FOM 	= Ct**(3/2) / (sqrt(2) * Cp)

	c 		= 0.5* D[index] * 0.121867779 # [m] chord at r/R = 0.75
	Re 		= rho * (df['Motor Speed [1/s]'] * c) / mu # [-] Reynolds number at given rpm
	M 		= (df['Motor Speed [1/s]'] * D[index] * 0.5) / a # [-] Mach number at given rpm for tip speeds

	df['Mechanical Power [W]'] = P_mech
	df['Electrical Power [W]'] = P_elec
	df['Thrust Coefficient [-]'] = Ct
	df['Power Coefficient [-]'] = Cp
	df['Figure of Merit [-]'] = FOM
	df['Reynolds Number [-]'] = Re
	df['Mach Number [-]'] = M
	df = df.iloc[1:] # Remove first row from dataframe

	# ==============================================================================================
	# REGRESSION ANALYSIS OF EXPERIMENTAL DATA
	# ==============================================================================================
	pfit_thrust = np.polyfit(df['Motor Speed [1/s]'], df['Thrust [N]'], deg=3)
	tl_thrust	= np.poly1d(pfit_thrust)
	R2_thrust 	= r2_score(df['Thrust [N]'], tl_thrust(df['Motor Speed [1/s]']))

	pfit_power 	= np.polyfit(df['Motor Speed [1/s]'], df['Mechanical Power [W]'], deg=3)
	tl_power	= np.poly1d(pfit_power)
	R2_power 	= r2_score(df['Mechanical Power [W]'], tl_power(df['Motor Speed [1/s]']))

	pfit_torque	= np.polyfit(df['Motor Speed [1/s]'], df['Torque [Nm]'], deg=3)
	tl_torque	= np.poly1d(pfit_torque)
	R2_torque 	= r2_score(df['Torque [Nm]'], tl_torque(df['Motor Speed [1/s]']))

	# ==============================================================================================
	# PLOT ROTATIONAL SPEED VS. THRUST
	# ==============================================================================================
	ax1[0,0].scatter(df['Motor Speed [1/s]'], df['Thrust [N]'], color=colors[index], marker='o')
	ax1[0,0].plot(df['Motor Speed [1/s]'], tl_thrust(df['Motor Speed [1/s]']), color=colors[index], label=bladetype[index])
	ax1[0,0].set_xlabel('Rotational Speed Motor $n$ [rev/s]')
	ax1[0,0].set_ylabel('Thrust [N]')
	ax1[0,0].set_xlim(0,50)
	ax1[0,0].set_ylim(0,60)
	ax1[0,0].minorticks_on() # set minor ticks
	ax1[0,0].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
	ax1[0,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid
	ax1[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k') # set legend for subplot

	# ==============================================================================================
	# PLOT ROTATIONAL SPEED VS. MECHANICAL POWER
	# ==============================================================================================
	ax1[0,1].scatter(df['Motor Speed [1/s]'], df['Mechanical Power [W]'], color=colors[index], marker='o')
	ax1[0,1].plot(df['Motor Speed [1/s]'], tl_power(df['Motor Speed [1/s]']), color=colors[index], label=bladetype[index])
	ax1[0,1].set_xlabel('Rotational Speed Motor $n$ [rev/s]')
	ax1[0,1].set_ylabel('Mechanical Power [W]')
	ax1[0,1].set_xlim(0,50)
	ax1[0,1].set_ylim(0,600)
	ax1[0,1].minorticks_on() # set minor ticks
	ax1[0,1].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
	ax1[0,1].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid

	# ==============================================================================================
	# PLOT ROTATIONAL SPEED VS. TORQUE
	# ==============================================================================================
	ax1[1,0].scatter(df['Motor Speed [1/s]'], df['Torque [Nm]'], color=colors[index], marker='o')
	ax1[1,0].plot(df['Motor Speed [1/s]'], tl_torque(df['Motor Speed [1/s]']), color=colors[index], label=bladetype[index])
	ax1[1,0].set_xlabel('Rotational Speed Motor $n$ [rev/s]')
	ax1[1,0].set_ylabel('Torque [Nm]')
	ax1[1,0].set_xlim(0,50)
	ax1[1,0].set_ylim(0,2)
	ax1[1,0].minorticks_on() # set minor ticks
	ax1[1,0].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
	ax1[1,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid

	# ==============================================================================================
	# PLOT DISK LOADING VS. POWER LOADING
	# ==============================================================================================
	DL = df['Thrust [N]'] / ((pi / 4) * D[index]**2)
	PL = df['Thrust [N]'] / df['Mechanical Power [W]']

	pfit_pl	= np.polyfit(DL, PL, deg=10)
	tl_pl	= np.poly1d(pfit_pl)
	R2_pl 	= r2_score(PL, tl_pl(DL))

	ax1[1,1].scatter(DL, PL, marker='o', color=colors[index])
	ax1[1,1].plot(DL, tl_pl(DL), color=colors[index], label=bladetype[index])
	ax1[1,1].set_xlabel(r'Disk Loading [N/$\mathrm{m}^2$]')
	ax1[1,1].set_ylabel('Power Loading [N/W]')
	ax1[1,1].set_xlim(0,100)
	ax1[1,1].set_ylim(0,0.4)
	ax1[1,1].minorticks_on() # set minor ticks
	ax1[1,1].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
	ax1[1,1].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid

# ==============================================================================================
# CALCULATE DESIGN POINTS
# ==============================================================================================
rpm_nom 	= np.interp(T_cruise / n_propellers, df['Thrust [N]'], df['Motor Speed [1/s]'])
T_nom 		= tl_thrust(rpm_nom)
P_nom 		= tl_power(rpm_nom)
torque_nom 	= tl_torque(rpm_nom)
PL_nom 		= T_nom / P_nom

rpm_max 	= np.interp(T_max / n_propellers, df['Thrust [N]'], df['Motor Speed [1/s]'])
T_max 		= tl_thrust(rpm_max)
P_max 		= tl_power(rpm_max)
torque_max 	= tl_torque(rpm_max)
PL_max 		= T_max / P_max

# ==============================================================================================
# PLOT DESIGN POINTS
# ==============================================================================================
ax1[0,0].axhline(y=T_nom, xmin=0, xmax=1, linestyle='--', color='black')
ax1[0,0].text(1, T_nom, r'T$_\mathrm{nom}$', fontsize=18, horizontalalignment='left', verticalalignment='bottom')
ax1[0,0].axhline(y=T_max, xmin=0, xmax=(1), linestyle='--', color='black')
ax1[0,0].text(1, T_max, r'T$_\mathrm{max}$', fontsize=18, horizontalalignment='left', verticalalignment='bottom')

ax1[0,1].axhline(y=P_nom, xmin=0, xmax=1, linestyle='--', color='black')
ax1[0,1].text(1, P_nom, r'P$_\mathrm{nom}$', fontsize=18, horizontalalignment='left', verticalalignment='bottom')
ax1[0,1].axhline(y=P_max, xmin=0, xmax=1, linestyle='--', color='black')
ax1[0,1].text(1, P_max, r'P$_\mathrm{max}$', fontsize=18, horizontalalignment='left', verticalalignment='bottom')

ax1[1,0].axhline(y=torque_nom, xmin=0, xmax=1, linestyle='--', color='black')
ax1[1,0].text(1, torque_nom, r'$\tau_\mathrm{nom}$', fontsize=18, horizontalalignment='left', verticalalignment='bottom')
ax1[1,0].axhline(y=torque_max, xmin=0, xmax=1, linestyle='--', color='black')
ax1[1,0].text(1, torque_max, r'$\tau_\mathrm{max}$', fontsize=18, horizontalalignment='left', verticalalignment='bottom')

ax1[1,1].axhline(y=PL_nom, xmin=0, xmax=1, linestyle='--', color='black')
ax1[1,1].text(1, PL_nom, r'PL$_\mathrm{nom}$', fontsize=18, horizontalalignment='left', verticalalignment='bottom')

fig1.savefig('performance_propellers.png', bbox_inches='tight', dpi=300)


# ==============================================================================================
# PLOT MOTOR DATA
# ==============================================================================================
def image_resize(image, width = None, height = None, inter = cv2.INTER_AREA):
    # initialize the dimensions of the image to be resized and
    # grab the image size
    dim = None
    (h, w) = image.shape[:2]

    # if both the width and height are None, then return the
    # original image
    if width is None and height is None:
        return image

    # check to see if the width is None
    if width is None:
        # calculate the ratio of the height and construct the
        # dimensions
        r = height / float(h)
        dim = (int(w * r), height)

    # otherwise, the height is None
    else:
        # calculate the ratio of the width and construct the
        # dimensions
        r = width / float(w)
        dim = (width, int(h * r))

    # resize the image
    resized = cv2.resize(image, dim, interpolation = inter)

    # return the resized image
    return resized

img1 = cv2.imread('Maxon_360.jpg') # read image
# img1 = image_resize(img1, height=1000)
height, width = img1.shape[:2] # determine size of image

dx, dy = 127, 73 # determine grid spacing
nx = width / (dx) + 1  # calculate number of grid lines on x-axis
ny = height / (dy) # calculate number of grid lines on y-axis

xmax = 1.5
ymax = 5000 / 60
xlabels = np.linspace(0, xmax, int(nx))
ylabels = np.linspace(0, ymax, int(ny))

# draw horizontal lines
for i in range(0, int(ny)):
	cv2.line(img1, (0, height - i*dy), (width, height - i*dy), (169,169,169), thickness=2) # draw major gridlines
	dy_minor = dy / 10
	for j in range(0, 12):
		dh = int(j * dy_minor)
		cv2.line(img1, (0, height - i*dy - dh), (width, height - i*dy - dh), (169,169,169), thickness=1) # draw minor gridlines

# draw vertical lines
for i in range(1, int(nx)):
	cv2.line(img1, (i*dx, height), (i*dx, 0), (169,169,169), thickness=2) # draw major gridlines
	cv2.putText(img1, str(round(xlabels[i],2)), (i*dx, height - 5), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,0,0), thickness=2) # add ticks on y-axis inside image
	dx_minor = dx / (10)
	for j in range(0, 10):
		dw = int(j * dx_minor)
		cv2.line(img1, ((i-1)*dx + dw, height), ((i-1)*dx + dw, 0), (169,169,169), thickness=1) # draw minor gridlines

for i in range(0, int(ny)):
	cv2.putText(img1, str(round(ylabels[i])), (0, height - i*dy -5), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,0,0), thickness=2) # add ticks on y-axis inside image

for i in range(1, int(nx)):
	cv2.putText(img1, str(round(xlabels[i],2)), (i*dx, height - 5), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,0,0), thickness=2) # add ticks on y-axis inside image

cv2.putText(img1, 'Motor Speed [1/s]', (5, 20), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,0,0), thickness=2) # add ticks on y-axis inside image
cv2.putText(img1, 'T[Nm]', (width - 50, height - 10), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,0,0), thickness=2) # add ticks on y-axis inside image

x1, y1 = torque_nom, rpm_nom
x2, y2 = torque_max, rpm_max

cv2.circle(img1, ((int(x1 / xmax * width), int(height - y1 / ymax * height))), radius=3, color=(0,0,0), thickness=6)
cv2.circle(img1, ((int(x2 / xmax * width), int(height - y2 / ymax * height))), radius=3, color=(0,0,0), thickness=6)
cv2.putText(img1, '32" Nom', ((int(x1 / xmax * width) + 20, int(height - y1 / ymax * height))), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,0,0), thickness=2) # add ticks on y-axis inside image
cv2.putText(img1, '32" Max', ((int(x2 / xmax * width) + 20, int(height - y2 / ymax * height))), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,0,0), thickness=2) # add ticks on y-axis inside image
cv2.putText(img1, 'Pmax = 360W', ((int(1.02 / xmax * width), int(height - 60 / ymax * height))), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0,0,0), thickness=2) # add ticks on y-axis inside image

cv2.imwrite('M360W.PNG', img1)


# ==============================================================================================
# ITERATE OVER DIFFERENT BATTERY CONFIGURATION FOR EACH PROPELLER DIAMETER
# ==============================================================================================
# Initialise figures
fig2, ax2 = plt.subplots(nrows=1, ncols=2, squeeze=False, figsize=(20,10))

# Colors used for iterative plotting
colors = ['tab:blue', 'tab:orange', 'tab:grey']

for index, df in enumerate(blades):

	# Initialise dataframe for hovering point which stores data of all propellers
	hover 	= pd.DataFrame(columns=['RPM [1/s]', 'Thrust [N]', 'Power [W]', 'Torque [Nm]', 'Mass [kg]', 'Endurance [min]'])
	cruise 	= pd.DataFrame(columns=['RPM [1/s]', 'Thrust [N]', 'Power [W]', 'Torque [Nm]', 'Mass [kg]', 'Endurance [min]'])
	maxop 	= pd.DataFrame(columns=['RPM [1/s]', 'Thrust [N]', 'Power [W]', 'Torque [Nm]', 'Mass [kg]', 'Endurance [min]'])

	for n in range(1, 200):
		C_parallel = n * E_battery
		m_parallel = n * m_battery 

		# Mass budget
		m_structures 	= 2.73 # [kg] operational empty mass
		m 				= m_sensors + m_autonomy + m_parallel + m_structures + m_motors + m_propellers
		m 				= 1.05 * m # [kg] contingency added

		theta_cruise 	= atan(D_cruise / (m * g)) # [rad] inclination angle during cruise
		theta_max 		= atan(D_max / (m * g)) # [rad] inclination angle during max speed

		# Define mission constraints
		T_hover 	= (m*g) / n_propellers # [N] Required thrust per propeller during hovering
		T_cruise 	= (m*g) / cos(theta_cruise) # [N] Required thrust during cruise
		T_max 		= (m*g) / cos(theta_max) # [N] Required thrust during max speed

		# ==============================================================================================
		# REGRESSION ANALYSIS OF EXPERIMENTAL DATA
		# ==============================================================================================
		pfit_thrust = np.polyfit(df['Motor Speed [1/s]'], df['Thrust [N]'], deg=3)
		tl_thrust	= np.poly1d(pfit_thrust)
		R2_thrust 	= r2_score(df['Thrust [N]'], tl_thrust(df['Motor Speed [1/s]']))

		pfit_power 	= np.polyfit(df['Motor Speed [1/s]'], df['Mechanical Power [W]'], deg=3)
		tl_power	= np.poly1d(pfit_power)
		R2_power 	= r2_score(df['Mechanical Power [W]'], tl_power(df['Motor Speed [1/s]']))

		pfit_torque	= np.polyfit(df['Motor Speed [1/s]'], df['Torque [Nm]'], deg=3)
		tl_torque	= np.poly1d(pfit_torque)
		R2_torque 	= r2_score(df['Torque [Nm]'], tl_torque(df['Motor Speed [1/s]']))

		# ==============================================================================================
		# CALCULATE DESIGN POINTS
		# ==============================================================================================
		rpm_dp1	= np.interp(T_hover, df['Thrust [N]'], df['Motor Speed [1/s]'])
		T_dp1	= tl_thrust(rpm_dp1)
		P_dp1 	= tl_power(rpm_dp1)
		torque_dp1 = np.interp(rpm_dp1, df['Motor Speed [1/s]'], df['Torque [Nm]'])
		E_dp1 	= (DoD * C_parallel) / ((P_dp1 / eta_mech) * n_propellers) * 60

		rpm_dp2 = np.interp(T_cruise / n_propellers, df['Thrust [N]'], df['Motor Speed [1/s]'])
		T_dp2 = tl_thrust(rpm_dp2)
		P_dp2 = tl_power(rpm_dp2)
		torque_dp2 = tl_torque(rpm_dp2)

		P_cruise = (n_propellers * P_dp2 + D_cruise * V_cruise)
		P_dp2 = P_cruise / n_propellers
		E_dp2 	= (DoD * C_parallel) / (P_cruise / eta_mech) * 60 

		rpm_dp3 = np.interp(T_max / n_propellers, df['Thrust [N]'], df['Motor Speed [1/s]'])
		T_dp3 = tl_thrust(rpm_dp3)
		P_dp3 = tl_power(rpm_dp3)
		torque_dp3 = tl_torque(rpm_dp3)

		P_max = (n_propellers * P_dp3 + D_max * V_max)
		P_dp3 = P_max / n_propellers
		E_dp3 = (DoD * C_parallel) / (P_max / eta_mech) * 60 

		hover.loc[n] 	= [rpm_dp1, T_dp1, P_dp1, torque_dp1, m_parallel, E_dp1] # add row to dataframe
		cruise.loc[n] 	= [rpm_dp2, T_dp2, P_dp2, torque_dp2, m_parallel, E_dp2] # add row to dataframe
		maxop.loc[n] 	= [rpm_dp3, T_dp3, P_dp3, torque_dp3, m_parallel, E_dp3] # add row to dataframe

	# ==============================================================================================
	# PLOT NUMBER OF BATTERY MODULES AGAINST ENDURANCE
	# ==============================================================================================
	ax2[0,0].scatter((cruise['Mass [kg]'] / m_battery), cruise['Endurance [min]'], marker='.', color=colors[index], label=bladetype[index])
	ax2[0,0].set_xlabel('Number of Battery Modules [-]')
	ax2[0,0].set_ylabel('Endurance at Nominal Power [min]')
	ax2[0,0].set_xlim(0,200)
	ax2[0,0].set_ylim(0,600)
	ax2[0,0].minorticks_on() # set minor ticks
	ax2[0,0].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
	ax2[0,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid
	ax2[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k') # set legend for subplot

	ax2[0,1].scatter((maxop['Mass [kg]'] / m_battery), maxop['Endurance [min]'], marker='.', color=colors[index], label=bladetype[index])
	ax2[0,1].set_xlabel('Number of Battery Modules [-]')
	ax2[0,1].set_ylabel('Endurance at Maximum Power [min]')
	ax2[0,1].set_xlim(0,200)
	ax2[0,1].set_ylim(0,350)
	ax2[0,1].minorticks_on() # set minor ticks
	ax2[0,1].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
	ax2[0,1].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid

ax2[0,0].scatter((cruise['Mass [kg]'].iloc[n_max-1] / m_battery), (cruise['Endurance [min]'].iloc[n_max-1]), s=100, marker='x', color='black')

fig2.savefig('endurance_propellers.png', bbox_inches='tight', dpi=300)


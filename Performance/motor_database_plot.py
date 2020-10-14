import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


motors_data = pd.read_excel("components_database.xlsx")
power = pd.DataFrame.to_numpy(motors_data.Power)
nom_volt = pd.DataFrame.to_numpy(motors_data.Nominal_Voltage)
torque = pd.DataFrame.to_numpy(motors_data.Maximum_Torque) / 1000

P_nom = 133.155
Torque_nom = 0.735

power_req = np.ones((50)) * P_nom
power_req_corres_torque = np.linspace(0, Torque_nom, 50)
torque_req = np.ones((50)) * Torque_nom
torque_req_corres_power = np.linspace(0, P_nom, 50)

max_volt = np.ones((50)) * 44.4
max_volt_corres_power = np.linspace(0, P_nom, 50)
power_req_corres_volt = np.linspace(0, Torque_nom, 50)

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


ax[0,0].set_xlim(left=0, right=np.max(torque)+0.08)
ax[0,0].set_ylim(0, np.max(power)+20)

ax[0,0].fill_between(np.linspace(0, np.max(torque)+0.08,50), 0, power_req, facecolor = 'red')
ax[0,0].fill_between(np.linspace(0, Torque_nom+0.005,50), 0, np.max(power)+20, facecolor = 'red', label='Infeasible')

ax[0,0].scatter(torque,power, s=60, marker = 'x', c='black' )
ax[0,0].set_xlabel('Nominal Torque [Nm]')
ax[0,0].set_ylabel('Power Output [W]')

# ax[0,0].set_xlim(left=0, right=np.max(nom_volt)+2)
# ax[0,0].set_ylim(0, np.max(power)+20)

# ax[0,0].fill_between(np.linspace(0, np.max(nom_volt)+2),0,power_req, facecolor = 'red')
# ax[0,0].fill_between(np.linspace(44.4, np.max(nom_volt)+2, 50), 0,np.max(power)+20, facecolor = 'red', label='Infeasible')
# ax[0,0].scatter(nom_volt,power, s=60, marker = 'x', c='black' )
# ax[0,0].set_xlabel('Nominal Voltage [V]')
# ax[0,0].set_ylabel('Power Output [W]')


ax[0,0].axvline(0,color=(0,0,0),linewidth=1.3)
ax[0,0].axhline(0,color=(0,0,0),linewidth=1.3)
ax[0,0].minorticks_on() # set minor ticks
ax[0,0].grid(which='major', linestyle='-', linewidth='0.5', color='black') # customise major grid
ax[0,0].grid(which='minor', linestyle=':', linewidth='0.5', color='grey') # customise minor grid

plt.legend()
fig.savefig("power_vs_power.png", bbox_inches='tight')            

plt.show()


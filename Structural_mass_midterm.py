"""Imports"""
import math
from math import pi
from decimal import *
import numpy as np
from itertools import product



#Structural Mass Airship Midterm
V_air = 229.85      #[m3]
r_air = 3.8         #[m]

S_air = 4*pi*r_air*r_air    #[m2]

ro_hull_list = [0.224653, 0.074, 0.18, 0.196, 0.208, 0.17] #[kg/m2]
material_list = ["Average","Tedlar","Zylon Polyurethane and Tedlar", "Vectran and Polyethane", "Zylon Polyurethane", "Vectran"]
m_struc_air = []

for i in range(0,len(ro_hull_list)):
    m_struc = round(ro_hull_list[i]*S_air,4)
    m_struc_air.append(m_struc)
    print("M_struc_air for is", m_struc_air[i], "kg considering", material_list[i])

#print(np.average(m_struc_air))
    

#Structural Mass Drone Midterm
r_prop = 0.4 #[m]
r_extra = 0.05 #[m]
w_arm = 0.04 #[m]
t_arm = 0.004 #[mm]
A_extra = 0.2*0.2 #[m2]
amount_rot = [8]

r_arm  = r_prop+ r_extra


material_drone = ["Aluminium Alloy 6060", "Aramid", "Carbon Fiber"]
density_drone = [2700, 1400, 1700]
m_struc_drone = []

for i in range(0, len(material_drone)):
    for j in range(0, len(amount_rot)):
        V_frame = amount_rot[j]*r_arm*w_arm*t_arm + A_extra*t_arm  #[m3]
        m_struc= round(density_drone[i]*V_frame,4)                      #[kg]
        m_struc_drone.append(m_struc)                           
        print("M_struc_drone for is", m_struc_drone[-1], "kg considering", material_drone[i], "and", amount_rot[j], "rotors")
print(V_frame*1000000)
#Structural Mass VTOL
S_w = np.arange(21.527821 ,43.055642, 1)                     #Wingarea [ft2]
W_fw = 3.840/0.453592                                        #weight fuel wing [lb]
A = np.arange(5,10,0.05)                                     #Aspect ratio [-]
taper = np.arange(0.2, 0.5, 0.01)                            #taper [-]
q = 2222#1993.8                                                   #dynamic pressure at cruise [lb/ft2]
tc = np.arange(0.1, 0.2, 0.005)                               #thickness over chord [-]
Nz = 1.5*1.5                                                 # loadfactor [-]
Wdg = 12.07/0.453592                                         #design gross weight [lb]


W_wing_list = []

#for s in S_w:
#    for a in A:
#        for t in taper:
#                for c in tc:
#                        r1 = 0.036*s**0.758
#                        r2 = W_fw**0.0035
#                        r3 = (a/(math.cos(0))**2)**0.6
#                        r4 = q**0.006
#                        r5 = t**0.04
#                        r6 = (100*c/(math.cos(0)))**(-0.3)
#                        r7 = (Nz*Wdg)**0.59
#                        W_wing = r1*r2*r3*r4*r5*r6*r7*0.453592 #[kg]
#                        W_wing_list.append(W_wing)
#                       
#print(np.amin(W_wing_list),np.amax(W_wing_list))
#print(np.average(W_wing_list))


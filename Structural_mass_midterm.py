"""Imports"""
from math import pi



#Structural Mass Airship Midterm
V_air = 229.85      #[m3]
r_air = 3.8         #[m]

S_air = 4*pi*r_air*r_air    #[m2]

ro_hull_list = [0.224653, 0.074, 0.18, 0.196, 0.208, 0.17] #[kg/m2]
material_list = ["Average","Tedlar","Zylon Polyurethane and Tedlar", "Vectran and Polyethane", "Zylon Polyurethane", "Vectran"]
m_struc_air = []

for i in range(0,len(ro_hull_list)):
    m_struc = ro_hull_list[i]*S_air
    m_struc_air.append(m_struc)
    #print("M_struc_air for is", m_struc_air[i], "kg considering", material_list[i])
    

#Structural Mass Drone Midterm
r_prop = 0.4 #[m]
r_extra = 0.05 #[m]
w_arm = 0.04 #[m]
t_arm = 0.004 #[mm]
A_extra = 0.25*0.25 #[m2]
amount_rot = 4

r_arm  = r_prop+ r_extra
V_frame = amount_rot*r_arm*w_arm*t_arm + A_extra*t_arm#[m3]

material_drone = ["Aluminium Alloy 6060", "Aramid", "Carbon Fiber"]
density_drone = [2700, 1400, 1700]
m_struc_drone = []

for i in range(0, len(material_drone)):
        m_struc= density_drone[i]*V_frame
        m_struc_drone.append(m_struc)
        print("M_struc_drone for is", m_struc_drone[i], "kg considering", material_drone[i])
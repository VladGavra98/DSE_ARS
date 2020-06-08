import numpy as np


concepts = 3
names = ["Airship", "Drone", "E-VTOL"]
#    Average weights and standard deviation
w_avg = np.array([14,29,22,19,16])
w_SD  = np.array([4.9	,5.83	,5.45	,4.84	,3.5])

#rankings:
R1 = 10
R2 = 7
R3 = 4

#Scores given to each concept:
#Concept order: Airship, Drone, E-VTOL (in this order)
#only risks will stay like this!!

#Operational Cost
CFR = np.array([R2,R1,R3]) #Cost and Freq of Repair
CC = np.array([R3,R2,R1]) #Cost of Consumables
CBM = np.array([R3,R1,R2])  #Cost Between Mission

#Flight Performance
MVFT = np.array([5,7.4,10]) #max vertical flight time
MHFT = np.array([0.7,8.57,10])  #max horizontal flight time
MHT = np.array([10,9.81,4.38])  #max hover time
Mo = np.array([R1,R2,R3])  #Mobility
St = np.array([R3,R1,R2])  #Stability

#Payload Application
SC = np.array([R1,R2,R3])  #Sensing Capability
SIA = np.array([R1,R2,R2]) #Signal Isolation Ability

#Risk
RISK = np.array([8,7.6,6.5])
#Sustainability
SUST = np.array([5.71,7.71,6.57])

OpCost = (CFR + CC + CBM) / 3
FlPerf1 = (MVFT + MHFT + MHT) / 3
FlPerf =  (FlPerf1*0.5 + Mo*0.25 + St*0.25)
PaAppl = (SC + SIA) / 2



S = np.zeros([5,3])
S[0,] = OpCost
S[1,] = FlPerf
S[2,] = PaAppl
S[3,] = RISK
S[4,] = SUST

# -*- coding: utf-8 -*-
"""
Created on Thu May 14 08:38:32 2020

@author: vladg,maxor
"""
import numpy as np
import random
import matplotlib.pyplot as plt
import itertools as it

concepts = 3
names = ["Airship", "Drone", "E-VTOL"]
#    Average weights and standard deviation
w_avg = np.array([14,29,22,19,16])
w_SD  = np.array([4.9	,5.83	,5.45	,4.84	,3.5])
percadj = 10 #adjustment of variations, =10 means 10%       #set to zero for no change
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
SUST = np.array([5.71,7.71,6.57]) #ADDed

OpCost = (CFR + CC + CBM) / 3
FlPerf1 = (MVFT + MHFT + MHT) / 3
FlPerf =  (FlPerf1 + Mo + St) / 3
PaAppl = (SC + SIA) / 2



S = np.zeros([5,3])
S[0,] = OpCost
S[1,] = FlPerf
S[2,] = PaAppl
S[3,] = RISK
S[4,] = SUST

grades = np.zeros((concepts))
top = []
sec = []
trd = []

f = open("sensitivity_scores.txt",'w')
f.writelines("Scores with all possible 1*SD combinations (now)\n\n\n")

for i,order in enumerate(list(set(list(it.permutations([1,1,1,-1,-1,-1,0,0,0],3))))):
    S2 = S
    for p in range(5): #runs through each row and tries the combination for each row (no combination of rows: too much)
        S2[p,] = S2[p,]+(percadj/10)*np.array(order)
        #print(S2)

        #print(order)
        for k in range(concepts):
            grades[k] = np.dot(w_avg,S2[:,k])/100
            # print(names[k],": ",grades[k])
        #print(names[np.argsort(grades)[2]],names[np.argsort(grades)[1]],names[np.argsort(grades)[0]])
        top.append(names[np.argsort(grades)[2]])
        sec.append(names[np.argsort(grades)[1]])
        trd.append(names[np.argsort(grades)[0]])
        f.writelines("\n"+str(order)+"\n")
        f.writelines(str(names[np.argsort(grades)[2]]+': '+str(round(grades[np.argsort(grades)[2]],2))+'  '+names[np.argsort(grades)[1]]+': '+str(round(grades[np.argsort(grades)[1]],2))+'  '+names[np.argsort(grades)[0]]+': '+str(round(grades[np.argsort(grades)[0]],2))+'\n'))
        #f.writelines(str(names[np.argsort(grades)[2]]))  uncomment for just first

print('\ntotal combinations: '+str(5*len((list(set(list(it.permutations([1,1,1,-1,-1,-1,0,0,0],3))))))))

print('\nFirst Count:')
print('Airship: '+str(top.count('Airship')))
print('Drone: '+str(top.count('Drone')))
print('E-VTOL: '+str(top.count('E-VTOL')))

print('\nSecond Count:')
print('Airship: '+str(sec.count('Airship')))
print('Drone: '+str(sec.count('Drone')))
print('E-VTOL: '+str(sec.count('E-VTOL')))

print('\nThird Count:')
print('Airship: '+str(trd.count('Airship')))
print('Drone: '+str(trd.count('Drone')))
print('E-VTOL: '+str(trd.count('E-VTOL')))

f.close()
print('\n\nDone')




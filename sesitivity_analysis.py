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
SDadj = 2 #adjustment of SD, =1 means examines one SD deviation

#Scores given to each concept:
#Concept order: Airship, Drone, E-VTOL (in this order)
#only risks will stay like this!!
S = np.array([[6,6.33,6],
              [6, 5.67, 6.33],
              [9,5,4],
              [8, 7.6, 6.5],
              [7,6,5]])

grades = np.zeros((concepts))
top = []
sec = []
trd = []

f = open("sensitivity.txt",'w')
f.writelines("Scores with all possible 1*SD combinations\n\n\n")


for i,order in enumerate(list(set(list(it.permutations([1,1,1,1,1,-1,-1,-1,-1,-1,0,0,0,0,0],5))))):
    # aux = w_avg
    delta = np.dot(SDadj*w_SD,np.array(order))
    # aux[n] = aux[n] + order
    aux = delta + w_avg
    aux = aux/sum(aux)*100

    print(order)
    for k in range(concepts):
        grades[k] = np.dot(aux,S[:,k])/100
        # print(names[k],": ",grades[k])
    print(names[np.argsort(grades)[2]],names[np.argsort(grades)[1]],names[np.argsort(grades)[0]])
    top.append(names[np.argsort(grades)[2]])
    sec.append(names[np.argsort(grades)[1]])
    trd.append(names[np.argsort(grades)[0]])
    f.writelines("\n"+str(order)+"\n")
    f.writelines(str(names[np.argsort(grades)[2]]+'  '+names[np.argsort(grades)[1]]+'  '+names[np.argsort(grades)[0]]+'\n'))
    #f.writelines(str(names[np.argsort(grades)[2]]))  uncomment for just first

print('\ntotal combinations: '+str(len((list(set(list(it.permutations([1,1,1,1,1,-1,-1,-1,-1,-1,0,0,0,0,0],5))))))))

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


print('\n\ndone')















# print("\nAlll weights + SD:")

# for i in range(len(w_avg)):
#     aux = w_avg
#     aux[i] = aux[i] + 1*w_SD[i]
#     aux = aux/sum(aux)*100
#     # print(aux,sum(aux))
#     for k in range(concepts):
#         grades[k] = np.dot(aux,S[:,k])/100
#         # print(names[k],": ",grades[k])
#     print(names[np.argsort(grades)[2]],names[np.argsort(grades)[1]],names[np.argsort(grades)[0]])


# print("\nAll weights - SD:")

# for i in range(len(w_avg)):
#     aux = w_avg
#     aux[i] = aux[i] - 1*w_SD[i]
#     aux = aux/sum(aux)*100
#     # print(aux,sum(aux))
#     for k in range(concepts):
#         grades[k] = np.dot(aux,S[:,k])/100
#         # print(names[k],": ",grades[k])
#     print(names[np.argsort(grades)[2]],names[np.argsort(grades)[1]],names[np.argsort(grades)[0]])

# for n in range(len(w_avg)):
#     aux = w_avg
#     aux[n] = aux[n] - 1*w_SD[n]
#     for l in range(len(w_avg)):
#         if l!=n:
#             for i in range(len(w_avg)):
#                 if i!=n and i!=l:
#                     aux[i] = aux[i] - 1*w_SD[i]
#                     for j in range(len(w_avg)):
#                         if j != i and j!=l and j!=n:
#                             # print(n,l,i,j)
#                             aux = aux/sum(aux)*100
#                             for k in range(concepts):
#                                 grades[k] = np.dot(aux,S[:,k])/100
#                                 # print(names[k],": ",grades[k])
#                             print(names[np.argsort(grades)[2]],names[np.argsort(grades)[1]],names[np.argsort(grades)[0]])

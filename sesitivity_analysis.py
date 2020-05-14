# -*- coding: utf-8 -*-
"""
Created on Thu May 14 08:38:32 2020

@author: vladg,maxr
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

#Scores given to each concept:
#Concept order: Airship, Drone, E-VTOL (in this order)
#only risks will stay like this!!
S = np.array([[5.33, 4.67, 3.67],
              [6.67, 8, 4.33],
              [5.5, 7.5, 1.5],
              [8, 7.6, 6.5],
              [6, 5, 2]])

grades = np.zeros((concepts))


f = open("sensitivity.txt",'w')
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


for i,order in enumerate(list(it.combinations([1,1,1,1,1,-1,-1,-1,-1,-1,0,0,0,0,0],5))):
    # aux = w_avg
    delta = np.dot(w_SD,np.array(order))
    # aux[n] = aux[n] + order
    aux = delta + w_avg
    aux = aux/sum(aux)*100

    print(order)
    for k in range(concepts):
        grades[k] = np.dot(aux,S[:,k])/100
        # print(names[k],": ",grades[k])
    print(names[np.argsort(grades)[2]],names[np.argsort(grades)[1]],names[np.argsort(grades)[0]])
    f.writelines(str(names[np.argsort(grades)[2]]+names[np.argsort(grades)[1]]+names[np.argsort(grades)[0]]+'\n'))

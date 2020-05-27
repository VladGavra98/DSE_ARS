# -*- coding: utf-8 -*-
"""
Created on Thu May 14 08:38:32 2020

@author: vladg,maxor
"""
import numpy as np
import random
import matplotlib.pyplot as plt
import itertools as it

concepts = 4
names = ["Li-Po", "Li-ion", "LiSu","Solar"]
#    Average weights and standard deviation
des = 'VTOL'  #airship, drone or VTOL

if des == 'airship':
    w_avg = np.array([30,10,10,24.5,10.5,15])
if des == 'drone':
    w_avg = np.array([16,24,10,24.5,10.5,15])
if des == 'VTOL':
    w_avg = np.array([12,28,10,24.5,10.5,15])

scr = 90 #  adjustment score (abs mult) 1 = 100%, 1.50=50% increase




S = np.zeros([6,4])
S[0,] = [4.6,5.8,9.1,10]
S[1,] = [8.57,6.14,10,2.73]
S[2,] = [8,10,7,5]
S[3,] = [9,9,0,1]
S[4,] = [9,6,10,9]
S[5,] = [6,7,5,9]

grades = np.zeros((concepts))
top = []
sec = []
trd = []
frt = []


f = open("sensitivity_weights.txt",'w')
f.writelines("Scores with all possible 1*SD combinations\n\n\n")


for i,order in enumerate(list(set(list(it.permutations([1,1,1,1,-1,-1,-1,-1,0,0,0,0],4))))):

    for p in range(4):
        delta = scr*(np.array(order))
        # aux[n] = aux[n] + order
        S[p,] = delta + S[p,]
        #aux = aux/sum(aux)*100
        #print(order)

        for k in range(concepts):
            grades[k] = np.dot(w_avg,S[:,k])/100


        top.append(names[np.argsort(grades)[3]])
        sec.append(names[np.argsort(grades)[2]])
        trd.append(names[np.argsort(grades)[1]])
        frt.append(names[np.argsort(grades)[0]])

        f.writelines("\n"+str(order)+"\n")
        f.writelines(str(names[np.argsort(grades)[3]])+': '+str(round(grades[np.argsort(grades)[3]],2))+'  '+str(names[np.argsort(grades)[2]])+': '+str(round(grades[np.argsort(grades)[2]],2))+'  '+names[np.argsort(grades)[1]]+': '+str(round(grades[np.argsort(grades)[1]],2))+'  '+names[np.argsort(grades)[0]]+': '+str(round(grades[np.argsort(grades)[0]],2))+'\n')
        #f.writelines(str(names[np.argsort(grades)[2]]))  uncomment for just first

ll = len(4*(list(set(list(it.permutations([1,1,1,1,-1,-1,-1,-1,0,0,0,0],4))))))

print('\n'+str(des)+' concept:  total combinations: '+str(ll))

lst = [top,sec,trd,frt]

for k in range(2): #CHANGE TO FOUR FOR ALL POSITIONS
    print('\nResult: '+str(k+1))
    for j in range(4):
        print(str(names[j])+'  '+str(lst[k].count(names[j]))+'   ('+str(round(100*lst[k].count(names[j])/ll,2))+'%)')


f.close()
print('\n\ndone')







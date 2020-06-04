# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 14:37:40 2020
EM = (total) emission map

This serves as the main script with the needed functions for the assement
@author: vladg
"""

import numpy as np
import matplotlib.pyplot as plt

a = 1
b=  1
c = 1
d = 10

class Point(np.ndarray):

    def __init__(self,location,distance,aicraft,noise,concentrantions):
        """ Initiliase the Point class with:
                 location = [x,y,z] np array
                 distance = time stamp at which the measurment was taken
                 aircraft = name of the aircraft -- str
                 noise = noise intensity value in AdB
                 concetrantions = [O3,CO,NO,NO2,SO2,PM,UFP,CH4]  -- KEEP THIS ORDER!

        """
        if(len(concentrantions)!=8):
            raise "Wrong number of recorded values!"

        self.x = location[0]
        self.y = location[1]
        self.z = location[2]
        self.aicraft  =  str(aicraft)
        self.concentrations = concentrantions
        self.O3      =  concentrantions[0]
        self.CO     =  concentrantions[1]
        self.NO      =  concentrantions[2]
        self.NO2      =  concentrantions[3]
        self.SO2      =  concentrantions[4]
        self.PM      =  concentrantions[5]
        self.UFP      = concentrantions[6]
        self.CH4     = concentrantions[7]
        self.noise   = noise


        def __str__(self):
            return str("Point at ") + str(self.location) + str(': ') +str(max(self.noise,max(self.concentrations)))

def calcMap(point):

    #unpack point
    x,y,z = point.location
    return a*x + b*y + c*z + d


points_lst = []
cons= np.ones((7))

point1= Point([1,1,1],100,'A320',80,cons)
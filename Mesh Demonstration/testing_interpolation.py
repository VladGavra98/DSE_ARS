# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 22:56:31 2020

@author: vladg
"""
import numpy as np
import numpy.linalg
import scipy as sp
import scipy.interpolate
import scipy.ndimage
import matplotlib.pyplot as plt
from interpolator import Spline_RBF

def foo(x):
    return np.sin(x)



xx = np.arange(0,10,0.1)
y = foo(xx)
xx_long = np.arange(-2,12,0.1)

grid = np.arange(2,8,0.1)   #much sparser grid
data = foo(grid)


epsilon_lst = np.arange(0.1,10,0.1)
error_lst = []

for epsilon_val in epsilon_lst:
    S = Spline_RBF([grid],data,ndim=1,epsilon=epsilon_val)
    y_hat = S.interpolate([xx])

    error_lst.append(np.sqrt(sum((y_hat-y)**2)))



# Plotting:
plt.plot(xx,y,'-b',label='real')
plt.plot(xx,y_hat,'-r', label='spline')
plt.legend(loc='best')
plt.figure("Error Plot")
plt.plot(epsilon_lst,error_lst)
plt.show()


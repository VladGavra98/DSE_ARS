# -*- coding: utf-8 -*-
'''

	AE3200 Design Synthesis Exercise
	Group 09 - Autonomous Environmental Sensing

   Created on Tue Jun  9 22:56:31 2020
	@author: vladg

	Project Supervisors:
		- Dr. Irene C. Dedoussi
		- Dr. Ir. Mirjam Snellen
		- Ir. Lorenzo Pasqualetto Cassinis
		- Mark Schelbergen

	This script is used for the initial sizing of the propellers / aerodynamic surfaces for the concept selection
'''

import numpy as np
import numpy.linalg
import scipy as sp
import scipy.interpolate
import scipy.ndimage
import matplotlib
import matplotlib.pyplot as plt
from interpolator import Spline_RBF, Point

def foo(x,y,z,):
    return 5 + np.exp(-x*x) + np.exp(-z*z) + np.exp(-y*y)


refine_factor = 2

def main(refine_factor):
    R = 5000
    h =150
    N_points = int(2 * R / h) #in one direction only



    # Domain to simulate:
    xx = yy = np.linspace(-R,R,refine_factor *  N_points)
    zz = np.linspace(0,50,refine_factor * N_points)
    print("Domain size:", len(xx)**2 * len(zz))
    # Calulcate real value of the metric:
    psi = foo(xx,yy,zz)

    # Perform measurements:
    gridx = gridy = np.linspace(-R,R,N_points)   #much sparser grid
    gridz = np.zeros(N_points)
    gridz[int(N_points/3):int(2*N_points/3)] = 15
    gridz[int(2*N_points/3):] = 30

    data = foo(gridx,gridy,gridz)

    # Store points in class:
        # (not working now but for later)
    points = Point([gridx,gridy,gridz],0,data)

    epsilon_lst = np.arange(120,170,1)
    error_lst = []

    # # Interpolate:
    # S = Spline_RBF([gridx,gridy,gridz],data,ndim=3,basis_name='gaussian')
    # psi_hat,weights = S.interpolate([xx,yy,zz])

    # error_lst.append(np.sqrt(sum((psi_hat-psi)**2)))

    # print("Total RMS error:", error_lst)
    # print("Relative error [%]:", sum(abs((psi_hat-psi)/psi)) *100)

    for epsilon_val in epsilon_lst:
        # Interpolate:
        S = Spline_RBF([gridx,gridy,gridz],data,ndim=3,epsilon=epsilon_val, basis_name="air_quality")
        psi_hat,weights = S.interpolate([xx,yy,zz])

        error_lst.append(sum((psi_hat-psi)**2)/len(psi))



    # Plotting:
    plt.figure("Error Plot")

    texpsize = [15, 18, 22]

    plt.rc('font', size=texpsize[1], family='serif')  # controls default text sizes
    plt.rc('axes', titlesize=texpsize[1])  # fontsize of the axes title
    plt.rc('axes', labelsize=texpsize[1])  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=texpsize[0])  # fontsize of the tick labels
    plt.rc('ytick', labelsize=texpsize[0])  # fontsize of the tick labels
    plt.rc('legend', fontsize=texpsize[0])  # legend fontsize
    plt.rc('figure', titlesize=texpsize[2])  # fontsize of the figure title
    matplotlib.rcParams['lines.linewidth'] = 1.5
    matplotlib.rcParams['figure.facecolor'] = 'white'
    matplotlib.rcParams['axes.facecolor'] = 'white'
    matplotlib.rcParams["legend.fancybox"] = False




    plt.plot(epsilon_lst,error_lst,label=str(refine_factor))
    plt.xlabel(r'$\epsilon$ [-]')
    plt.ylabel('MSE [-]')
    plt.title(r'MSE vs. $\epsilon$')
    plt.axvline(0,color=(0,0,0),linewidth=0.8) #comment out for no axes line
    plt.axhline(0,color=(0,0,0),linewidth=0.8) #comment out for no axes line
    plt.minorticks_on() # set minor ticks
    plt.grid(which='major', linestyle='-', linewidth='0.3', color='black') # customise major grid
    plt.grid(which='minor', linestyle=':', linewidth='0.3', color='grey') # customise minor grid
    plt.legend()
    plt.show()

    # plt.plot(xx,y,'-b',label='real')
    # plt.plot(xx,y_hat,'-r', label='spline')
    # plt.legend(loc='best')

plt.show()

for k in range(5,6):
    main(k)
# -*- coding: utf-8 -*-
'''

	AE3200 Design Synthesis Exercise
	Group 09 - Autonomous Environmental Sensing

   Created on Tue Jun  14  2020
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
import time
from interpolator import Spline_RBF, Point
import matplotlib.pyplot as plt
import matplotlib

# Data frame format:
#  'NO'  : {'s1': 'Badhoevedorp', 's2': 'Hoofddorp', 's3': 'OudeMeer'}

def air_qaulity():
    # LOad and flatten data:
    N       = 100 # points
    data_NO = np.genfromtxt("NO_data.txt",max_rows=N)
    psi     = data_NO.flatten() * 10**6
    Nfactor = 2

    # quick & dirty check:
    for item in list(psi):
        if item=='nan':
            print(item)


    """
    Coordinates of the stations (x,y):

        Badhoevedorp : 52.333, 4.7738
        Hoofddorp    : 52.327, 4.716
        OudeMeer     : 57.279, 4.770

        Order: B,H,O

        OudeMeer is the orginin of the locally defined coordinate system

    """
    xB,yB = 52.333, 4.7738
    xH,yH =  52.327, 4.716

    xO,yO =  57.279, 4.770

    xB,yB = 52.333-xO, 4.7738- yO
    xH,yH =  52.327 - xO, 4.716 - yO

    xtab = ytab= ztab = np.zeros(len(psi))

    xtab[:N] = xB
    xtab[N:2*N] = xH
    xtab[2*N:] = xO

    ytab[:N] = xB
    ytab[N:2*N] = xH
    ytab[2*N:] = xO


    timetab = range(len(psi))  #time in hours from the first measurement

    epsilon_lst = np.arange(0.8, 1, 0.01)
    error_lst = []

    # # Interpolate:
    # S = Spline_RBF([gridx,gridy,gridz],data,ndim=3,basis_name='gaussian')
    # psi_hat,weights = S.interpolate([xx,yy,zz])

    for epsilon_val in epsilon_lst:
            # Interpolate:
            S = Spline_RBF([xtab[::Nfactor],ytab[::Nfactor],ztab[::Nfactor],timetab[::Nfactor]],psi[::Nfactor],ndim=4,epsilon=epsilon_val,basis_name="air_quality")
            psi_hat,weights = S.interpolate([xtab,ytab,ztab,timetab])

            error_lst.append(sum((psi_hat-psi)**2)/len(psi))



    # Plotting:
    plt.figure("Error Plot Air Quality")

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




    plt.scatter([0.9],[8.588],color='r')
    plt.plot(epsilon_lst,error_lst,label="Error")
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

    print("Total RMS error:", error_lst)


# %% Noise modelling validation:
def noise():
    # LOad and flatten data:
    N       = 500 # points
    data = np.genfromtxt("output.txt",max_rows=N)

    data = data[:-1]
    psi     = []

    # for i in range(5,len(data),5):
    #     psi.append(np.average(data[i-1:i]))

    psi = np.array(data)
    psi = psi[psi!=0]
    Nfactor = 2 #skip each Nfactor points


    timetab = range(len(psi))  #time in hours from the first measurement

    epsilon_lst = np.arange(2, 5, 0.1)
    error_lst = []

    # Interpolate:
    # S = Spline_RBF([timetab[::Nfactor]],psi[::Nfactor],ndim=1,basis_name='gaussian')
    # psi_hat,weights = S.interpolate([timetab])

    for epsilon_val in epsilon_lst:
            # Interpolate:
            S = Spline_RBF([timetab[::Nfactor]],psi[::Nfactor],ndim=1,epsilon=epsilon_val,basis_name="noise")
            psi_hat,weights = S.interpolate([timetab])

            error_lst.append(sum((psi_hat-psi)**2)/len(psi))

    plt.figure("Noise data (testing)")
    plt.scatter(timetab,psi)
    print(psi,data)



    # Plotting:
    plt.figure("Error Plot Noise")

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





    plt.plot(epsilon_lst,error_lst,label="Error")
    plt.scatter([3.5],[9.6675],color='r')
    plt.xlabel(r'$\epsilon$ [-]')
    plt.ylabel('MSE [%]')
    plt.title(r'MSE vs. $\epsilon$')
    plt.axvline(0,color=(0,0,0),linewidth=0.8) #comment out for no axes line
    plt.axhline(0,color=(0,0,0),linewidth=0.8) #comment out for no axes line
    plt.minorticks_on() # set minor ticks
    plt.grid(which='major', linestyle='-', linewidth='0.3', color='black') # customise major grid
    plt.grid(which='minor', linestyle=':', linewidth='0.3', color='grey') # customise minor grid
    plt.legend()
    plt.show()

    # plt.plot(timetab,psi,'-b',label='real')
    # plt.plot(timetab,psi_hat,'-r', label='spline')
    # plt.legend(loc='best')

    # plt.show()



air_qaulity()
noise()
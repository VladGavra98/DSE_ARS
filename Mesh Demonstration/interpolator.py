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
import time



#++++++++++++++++++++++++++++++++++++++++++++ Data Point Class +++++++++++++++++++++++++++++++++++
class Point:
    """Data point self-defined class, use if numpy arrays are not sufficient"""
    def __init__(self,position,time, psi):

        self.position =position
        self.time = time
        self.psi = psi



#+++++++++++++++++++++++++++++++++++++ Cubic Spline Class +++++++++++++++++++++++++++++++++++++++++
class Spline (sp.interpolate.CubicSpline):

    '''Class that generates cubic (natural) splines'''

    def __init__(self,data,pos):
        """
        Interpoalte unknown function function from data, it takes as arguments:
                s(x) = a*(x-x0)**3 + b*(x-x1)**2 + c*(x-x2) +d

            data = Psi value at a grid point
            pos =  numpy ndarray containing the grid locations

            Returns:
                Cubic spline coffieints in numpy array: [a,b,c,d]


        """
        self.data   = np.array(data)
        self.grid   = np.array(pos)
        self.n      = len(self.data)
        self.h      = self.grid[1:] - self.grid[:-1]

        if self.data.shape != self.grid.shape:

            raise TypeError ("Data has wrong dimension")




class interpolate:

    ''' see example below for how to get the coefficients. Keep in mind its ai(x-xi)^3 etc'''
    def __init__(self, data, pos):
        self.data   = data
        self.n      = len(data)

        #Create a list of all the distances between points.
        self.h      = np.zeros(self.n-1)
        for i in range(self.n-1):
            self.h[i]  = pos[i+1]-pos[i]


        #Setup the matrixes
        co_matrix   = np.zeros((self.n-2,self.n-2))
        f_matrix    = np.zeros((self.n-2,1))

        for i in range(self.n-2):
            co_matrix[i][i]         = (self.h[i]+self.h[i+1])/3

            if i != 0:
                co_matrix[i][i-1]   = self.h[i]/6
                co_matrix[i-1][i]   = self.h[i]/6

            #f matrix is shifted because i is defined weird due to m0 and mn being 0 and not part of the comatrixes


            f_matrix[i][0]          = (data[i+2] - data[i+1])/self.h[i+1] - (data[i+1] - data[i])/self.h[i]


        m_matrix    = np.linalg.solve(co_matrix,f_matrix)





        m_matrix    = np.linalg.solve(co_matrix,f_matrix)

        #Add boundary conditions for m_matrix
        m_matrix    = np.vstack((np.array([0]),m_matrix))
        m_matrix    = np.vstack((m_matrix,np.array([0])))

        self.abcd   = np.zeros((self.n - 1 , 4))



        for i in range(self.n-1):

            ai  = (m_matrix[i+1] - m_matrix[i] )/ (6 * self.h[i])
            bi  = m_matrix[i] / 2
            ci  = (self.data[i+1] - self.data[i]) / self.h[i] - self.h[i] / 3 * m_matrix[i] - self.h[i] / 6 * m_matrix[i+1]
            di  = self.data[i]

            self.abcd[i] = [ai,bi,ci,di]

#++++++++++++++++++++++++++++++++++ RBF Class +++++++++++++++++++++++++++++++++++++++++++++++++++++

class Spline_RBF(sp.interpolate.Rbf):
    """
        Spline class based on scipy radial basis functions, generates multidimensional splines
    """

    def __init__(self,grid,data,ndim,epsilon=None,basis_name=None):

        """
        Interpoalte unknown function from data, it takes as arguments:


            data = Psi value at a grid point
            pos =  numpy ndarray containing the grid locations  -- list of grid values
            ndim = numer of dimensions of the gird  -- int
            epsilon = valuea for epsilon optimisation
            basis_name = string or callable to be used inside the RBF constructor -- str/callable
                        - noise
                        - air_quality
                        - None (default)
            Returns:
                rbfObj

        """
        self.data   = np.array(data)
        self.grid   = np.array(grid)

        if len(self.grid) != ndim:

            raise TypeError ("Grid has wrong dimension, use a list around the grid")

        self.dimensions = ndim

        # if isinstance(basis, str):
        #     self.basis = basis.lower()
        # elif callable(basis):
        #     self.basis = basis
        # else:
        #     raise ValueError("Basis given is not valid!")

        self.epsilon = epsilon
        self.basis_name = basis_name
        print(self.basis_name)


    def basis_airquality(self,r):
        return np.exp(-r*r/(self.epsilon * self.epsilon))

    def basis_noise(self,r):

         return  1/np.sqrt((self.epsilon + r*r))


    def interpolate(self,target_grid):
        """
        Returns the interpoalted values at the target grid values:
            Input:  target_grid  -- numpy nd array
            Output: interpoalted value(s) at grid location   -- numpy ndarray
                    weights -- of each basis function -- numpy ndarray,
                                same length as data grid dimensions
        """



        if self.basis_name !=None:
            if self.basis_name == "air_quality" or self.basis_name == "air quality":
                # no name specified, use the self defined one
                self.rbfObj = sp.interpolate.Rbf(*self.grid,self.data,function=self.basis_airquality)

            elif self.basis_name == "noise":
                # no name specified, use the self defined one
                self.rbfObj = sp.interpolate.Rbf(*self.grid,self.data,function=self.basis_noise)

            else:
                # Basis name spceified, let Scipy find the correct base in its list
                # Epsilon parameter cannot be used anymore

                self.rbfObj = sp.interpolate.Rbf(*self.grid,self.data,function=self.basis_name)
        else:
            raise Warning("No function is given so default scipy will be used!")
            self.rbfObj = sp.interpolate.Rbf(*self.grid,self.data)


        weights = self.rbfObj.nodes

        # tested with self defined function:
        # RMS error depends on epsilon so 1d optimisation can be performed
        # tested with gaussian basis: low error is achieved and epislon,  indeed, does not matter


        return self.rbfObj(*target_grid),weights



# ++++++++++++++++++++++++++++++++ Helper Functions ++++++++++++++++++++++++++++++++++++++++++++++++
def normalise(x):

    """
    Helper function to normlise the a vector
    Can be used to remove the units of physical dimension
    """
    x =np.array(x)
    maxim = np.max(x)

    if maxim :
        return x/maxim
    return x

def main():

    # #Example 4.2 from the ANA reader
    # data = [0,0.2624,0.6419,1.0296,]
    # pos = [0,0.1,0.3,0.6]

    # S_hat = Spline(data,pos).abcd  #coefficients

    # # Scipy univariate cubic interpolation:
    # S =  sp.interpolate.CubicSpline(pos,data,bc_type= 'natural')
    # x = 0.5
    # print("Spline value at %f is %f \n" %(x, S(x)))


    # Testing....
    # Geberate RBF spline:
    x = y= z= t= np.arange(-30,30,1)
    y = np.zeros(len(x)) +1
    x = normalise(x)
    y = normalise(y)
    z = normalise(z)
    t = normalise(t)

    data = np.linspace(10,20,len(x))

    print("Grid dimesnion: %i points " %(len(x)**4))

    # spline2 = Spline_RBF([x,y,z,t],data,4)
    spline2 = Spline_RBF([x,y,z,t],data,4)    #works with a given new basis function (user defined)


    xi = yi= zi = ti= np.linspace(0,1,30)
    # data_hat = spline2.interpolate([xi,yi,zi]) #should return 20 numbers -> correct!
    # data_hat = spline2.interpolate([1,1,1])  # should return one number  -> correct!
    data_hat, basis_weights = spline2.interpolate([xi,yi,zi,ti])  # should return one number, outside of domain (extrapolate)  -> correct!

    print(data_hat)

def test():
    x, y, z, d = np.random.rand(4, 50)
    rbfi = sp.interpolate.Rbf(x, y, z, d)  # radial basis function interpolator instance
    xi = yi = zi = np.linspace(0, 1, 20)
    di = rbfi(xi, yi, zi)   # interpolated values
    print(di.shape)

if __name__=='__main__':

    # test()
    start = time.time()
    main()
    end= time.time()

    print("Elapsed time [s]: ", end-start)



#Example 4.2 from the ANA reader
#data = [0,0.2624,0.6419,1.0296]
#pos = [0,0.1,0.3,0.6]
#print(interpolate(data,pos).abcd)

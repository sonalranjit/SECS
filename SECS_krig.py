__author__ = 'sonal'
from geostatsmodels import utilities, variograms, model, kriging, geoplot
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *

def lla2ecef(lla):
    #source: http://github.com/timduly4/pyglow/blob/master/pyglow/coord.py
    #Constats WGS84
    a = 6378137.
    b = 6356752.3142
    e2 = 1-(b**2)/(a**2)

    # check for 1D case:
    dim = len(lla.shape)
    if dim == 1:
        lla = np.reshape(lla,(1,3))

    # convert lat and on to radians
    lat = lla[:,0]/180.*pi
    lon = lla[:,1]/180.*pi

    xyz = np.array(np.zeros(lla.shape))

    N = a/np.sqrt(1-e2*np.sin(lat)*np.sin(lat))
    # Calculate the X-coordinate
    xyz[:,0] = (N)*np.cos(lat)*np.cos(lon)

    # Calculate the Y-coordinate
    xyz[:,1] = N*np.sin(lon)*np.cos(lat)

    # Keep the SECS data as it is
    xyz[:,2] = lla[:,2]

    return np.array(xyz)

data = np.loadtxt('/home/sonal/SECS_EICS/SECS/SECS20110301/01/SECS20110301_120000.dat')
xyz = lla2ecef(data)
pt = lla2ecef(np.array([43.7, 280.6,0]))
tolerance = 100000
lags =np.arange(tolerance, 3000000, tolerance)
sill =np.var(xyz[:,2])

#svm = model.semivariance(model.exponential, (900000, sill))
#geoplot.semivariogram(xyz, lags, tolerance, model=svm)

covfct = model.covariance(model.exponential, (900000, sill))
ptz = kriging.simple(xyz,covfct,pt[:,:2],N=10)

'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2],c='g',marker = '^')
ax.scatter(pt[0,0],pt[0,1],ptz[0],c='r')
plt.show()'''
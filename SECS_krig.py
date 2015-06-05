__author__ = 'sonal'
from geostatsmodels import utilities, variograms, model, kriging, geoplot
import numpy as np
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os
from math import *

def lla2ecef(lla):
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

sat_data = np.loadtxt('/home/sonal/SECS/sat_data_march.txt')
zero_col = np.zeros((len(sat_data),1))
sat_data = np.column_stack((sat_data,zero_col))
print sat_data.shape
for i in range(0,100):
#for i in range(len(sat_data)):
    print "Processing file "+str(i)+" of "+str(len(sat_data))
    secs_path = '/home/sonal/SECS_EICS/SECS/'
    sat_y = str(int(sat_data[i,0]))
    sat_m = str(int(sat_data[i,1])).zfill(2)
    sat_d = str(int(sat_data[i,2])).zfill(2)
    sat_ymd = sat_y+sat_m+sat_d
    sat_h = str(int(sat_data[i,3])).zfill(2)
    sat_mins = str(int(sat_data[i,4])).zfill(2)
    sat_secs = str(int(floor(sat_data[i,5]))).zfill(2)
    sat_hms = sat_h+sat_mins+sat_secs

    SEC_file = secs_path+'SECS'+sat_ymd+'/'+sat_d+'/'+'SECS'+sat_ymd+'_'+sat_hms+'.dat'

    if os.path.exists(SEC_file):

        sec_grid = np.loadtxt(SEC_file)
        grid_xyz = lla2ecef(sec_grid)

        sat_latlon = np.zeros((1,3))
        sat_latlon[:,(0,1)] = sat_data[i,(6,7)]
        sat_xyz = lla2ecef(sat_latlon)

        sill = np.var(grid_xyz[:,2])
        covfct = model.covariance(model.exponential,(900000, sill))

        ptz = kriging.simple(grid_xyz,covfct,sat_xyz[:,:2],N=10)

        sat_data[i,8] = ptz[0]

np.savetxt('test.txt', sat_data)
'''
data = np.loadtxt('/home/sonal/SECS_EICS/SECS/SECS20110301/01/SECS20110301_000000.dat')
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

'''#Plotting
lons, lats = np.meshgrid(data[:,1],data[:,0])
Z = data[:,2]
m = Basemap(width=12000000, height=8000000, resolution='l', projection='laea',\
            lat_ts=min(data[:,0]), lat_0=np.median(data[:,0]),lon_0=np.median(data[:,1]))

m.drawcoastlines()
m.drawparallels(np.arange(-80.,81.,20.))
m.drawmeridians(np.arange(-180.,181.,20.))
x,y =m(data[:,1],data[:,0])
cs = m.contour(lons,lats,Z,8000,linewidth=0.5,colors='k',latlon=True)
#m.scatter(x,y)
plt.show()'''
'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2],c='g',marker = '^')
ax.scatter(pt[0,0],pt[0,1],ptz[0],c='r')
plt.show()'''
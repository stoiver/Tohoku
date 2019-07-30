#!/usr/bin/env python

import os
import sys
from mpl_toolkits.basemap import Basemap,interp
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp2d
#from anuga.shallow_water.shallow_water_domain import Domain
#from anuga.abstract_2d_finite_volumes.mesh_factory import rectangular_cross
#from anuga.abstract_2d_finite_volumes.quantity import Quantity
#from matplotlib.mlab import griddata
sys.path.append("/usr/lib64/python2.6/dist-packages/")
from scipy.interpolate import griddata
from Scientific.IO.NetCDF import NetCDFFile
from pyproj import Proj


fileroot = sys.argv[1]
print 'plotting %s.dem' % fileroot

def plot_netcdf(filename):
    ncf = NetCDFFile(filename,'r')
    # Check if not ANUGA dem
    if not hasattr(ncf,'cellsize'):
        print 'File %s is not an ANUGA DEM (no cellsize attribute)' % filename
        return(-1)
    # Set UTM grid params from netcdf header
    cellsz = ncf.cellsize[0]
    nrows  = ncf.nrows[0]
    ncols  = ncf.ncols[0]
    xll    = ncf.xllcorner[0]   # Easting of lower left corner
    yll    = ncf.yllcorner[0]  # Northing of lower left corner
    xur    = xll + (ncols-1)*cellsz
    yur    = yll + (nrows-1)*cellsz
    zone   = ncf.zone[0]
    zdat   = np.flipud(ncf.variables['elevation'][:])

    # Transform corners to lat,lon
    p = Proj(proj='utm',zone=zone,ellps='WGS84')
    lnll,ltll = p(xll,yll,inverse=True) 
    lnur,ltur = p(xur,yur,inverse=True)
    print lnll,ltll,lnur,ltur 
    # Transform from utm to lat,lon grid
    X_utm = np.linspace(xll,xll+(ncols-1)*cellsz,ncols)
    Y_utm = np.linspace(yll,yll+(nrows-1)*cellsz,nrows)
    X_utm,Y_utm = np.meshgrid(X_utm,Y_utm)
    lns,lts = p(X_utm,Y_utm,inverse=True)
    lnsi = np.linspace(lnll,lnur,1000)
    ltsi = np.linspace(ltll,ltur,1000)
    zi = griddata((lns.flatten(), lts.flatten()), zdat.flatten(), \
                  (lnsi[None,:], ltsi[:,None]), method='linear', \
                      fill_value=0.)
    # 
    map = Basemap(lnll,ltll,lnur,ltur,projection='merc',resolution='h')

    x, y = map(lnsi,ltsi)
    xx, yy = np.meshgrid(x,y)
    CS = map.contour(xx,yy,zi,15,linewidths=0.5,colors='k')
    CS = map.contourf(xx,yy,zi,15,cmap=plt.cm.jet)
#    plt.colorbar() # draw colorbarsys.exit(0)
    map.drawcoastlines()
    map.drawcountries()
    map.fillcontinents(color = 'coral')
    map.drawmapboundary()
    map.drawmeridians(np.arange(0, 360, 1),labels=[1,0,0,1])
    map.drawparallels(np.arange(-90, 90, 1),labels=[1,0,0,1])

Fig = plt.figure()
for i in range(0,1):
#    plt.subplot(320+i+1)
    plot_netcdf(fileroot+'.dem')
plt.show()

sys.exit(0)

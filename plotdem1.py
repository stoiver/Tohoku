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



def read_dem(filename):
    ncf = NetCDFFile(filename,'r')
    # Set ANUGA file - UTM grid params from netcdf header
    if hasattr(ncf,'cellsize'):
        cellsz = ncf.cellsize[0]
        nrows  = ncf.nrows[0]
        ncols  = ncf.ncols[0]
        xll    = ncf.xllcorner[0]   # Easting of lower left corner
        yll    = ncf.yllcorner[0]  # Northing of lower left corner
        xur    = xll + (ncols-1)*cellsz
        yur    = yll + (nrows-1)*cellsz
        x = np.linspace(xll,xur,ncols)
        y = np.linspace(yll,yur,ncols)
        zone   = ncf.zone[0]
        zone = 51
        zdat   = np.flipud(ncf.variables['elevation'][:])
    # Made from GMT?
    elif 'x_range' in ncf.variables:
        xrng = ncf.variables['x_range']
        yrng = ncf.variables['y_range']
        zrng = ncf.variables['z_range']
        dxdy = ncf.variables['spacing']
        lnll = xrng[0]
        lnur = xrng[1]
        ltll = yrng[0]
        ltur = yrng[1]
        # Pixel registration
        if ncf.variables['z'].node_offset[0] == 1:
            lnll += 0.5*dxdy[0]
            lnur -= 0.5*dxdy[0]
            ltll += 0.5*dxdy[1]
            ltur -= 0.5*dxdy[1]
        nx = 1+int(0.1+(lnur-lnll)/dxdy[0])
        ny = 1+int(0.1+(ltur-ltll)/dxdy[1])
        zdat = ncf.variables['z'][:].reshape(ny,nx)
        x = np.linspace(lnll,lnur,nx)
        y = np.linspace(ltur,ltll,ny)
    elif ncf.Conventions == 'COARDS/CF-1.0':
        x = ncf.variables['x'][:]
        y = ncf.variables['y'][:]
        zdat=ncf.variables['z'][:]
    else:
        print 'Not GDAL/GMT netcdf file: %s' % filename
        return(-1)
    return (x,y,zdat)

def plot_netcdf(xll,yll,xur,yur,zdat,zone,cellsz):
    (nrows,ncols) = zdat.shape
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
    map = Basemap(lnll,ltll,lnur,ltur,projection='merc',resolution='i')

    x, y = map(lnsi,ltsi)
    xx, yy = np.meshgrid(x,y)
    CS = map.contour(xx,yy,zi,15,linewidths=0.5,colors='k')
    CS = map.contourf(xx,yy,zi,15,cmap=plt.cm.jet)
    plt.colorbar() # draw colorbarsys.exit(0)
    map.drawcoastlines()
    map.drawcountries()
#    map.fillcontinents(color = 'coral')
    map.drawmapboundary()
    map.drawmeridians(np.arange(0, 360, 1),labels=[1,0,0,1])
    map.drawparallels(np.arange(-90, 90, 1),labels=[1,0,0,1])

filename = sys.argv[1]
if len(sys.argv) < 3:
    fileref = ''
    print 'plotting %s' % (filename)
else:
    fileref = sys.argv[2]
    print 'plotting %s - %s' % (filename,fileref)
Fig = plt.figure()
zdat0 = 0.
for i in range(0,1):
#    plt.subplot(320+i+1)
    (x,y,zdat)  = read_dem(filename)
    if fileref != '':
        (x,y,zdat0) = read_dem(fileref)
        zdat -= zdat0
    xll = x.min()
    xur = x.max()
    yll = y.min()
    yur = y.max()
    cellsize = x[1]-x[0]
    print xll,yll,xur,yur,cellsize
    plot_netcdf(xll,yll,xur,yur,zdat,54,cellsize) 
plt.savefig('dem3600.pdf')
plt.show()


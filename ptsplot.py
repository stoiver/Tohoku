#!/usr/bin/env python

import os
import sys
import numpy as np
from netCDF4 import Dataset as Dataset
import matplotlib.pyplot as plt
import matplotlib

cdict = {'red'  : ((0.0, 0.0, 0.0),(0.5, 0.0, 0.7),(1.0, 1.0, 1.0)),
         'green': ((0.0, 0.0, 0.0),(0.5, 1.0, 0.0),(1.0, 1.0, 1.0)),
         'blue': ((0.0, 0.0, 1.0),(0.5, 0.0, 0.0),(1.0, 0.5, 1.0))}
my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)


nci = Dataset(sys.argv[1],'r',format='NETCDF3_CLASSIC')
varname = 'elevation'
pts = nci.variables['points'][:]
elev = nci.variables['elevation'][:]
#elev = nci.variables['0.010'][:]
if len(sys.argv) > 2:
    zmnx = float(sys.argv[2])
    elev = np.where(elev>zmnx, zmnx,elev)
    elev = np.where(elev<-zmnx, -zmnx,elev)
fig = plt.figure(figsize=(12.,12.))
ax = fig.add_subplot(111)
sc = ax.scatter(pts[:,0],pts[:,1],c=elev,edgecolor='none',cmap=plt.cm.jet)
plt.colorbar(sc)
plt.show()

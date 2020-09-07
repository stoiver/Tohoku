
import matplotlib.pyplot as plt

import anuga


from scipy.io import loadmat


from numpy import loadtxt, array2string

infile = 'Elevation.txt'
elev = loadtxt(infile)

print elev


nrows, ncols = elev.shape
xllcorner = 419669.750000000
yllcorner = 4042429.17000000
cellsize = 4040.4
NODATA_value =  -9999


filename = 'okada.asc'
fid = open(filename, 'w')

msg = """ncols         %g
nrows         %g
xllcorner     %g
yllcorner     %g
cellsize      %g
NODATA_value  %g
""" % (ncols, nrows, xllcorner, yllcorner, cellsize, NODATA_value)

print msg

fid.write(msg)

#Create linear function
ref_points = []
ref_elevation = []
x0 = 2000
y = 3010
yvec = range(nrows)
xvec = range(ncols)

for i in range(nrows):
    srow = ' '.join([ str(x) for x in elev[i, :] ])
    fid.write(srow)
    fid.write('\n')

fid.close()




metafilename = 'okada.prj'
fid = open(metafilename, 'w')
fid.write("""Projection UTM
Zone 54
Datum WGS84
Zunits NO
Units METERS
Spheroid WGS84
Xshift 0.0000000000
Yshift 10000000.0000000000
Parameters
""")
fid.close()

import anuga

anuga.asc2dem(filename)
anuga.dem2pts(filename)

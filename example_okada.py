import numpy as num

#from anuga.tsunami_source.tsunami_okada import earthquake_tsunami

from anuga.tsunami_source.okada_tsunami import earthquake_source

import anuga

"""
Pick the test you want to do; T= 0 test a point source,
T= 1  test single rectangular source, T= 2 test multiple
rectangular sources
"""


# Choose what test to proceed
T = 1

if T == 1:
    # Initial condition of earthquake for multiple source
    x0 = 40000.0
    y0 = 40000.0
    depth = 15.0
    strike = 0.0
    dip = 15.0
    length = 20.0
    width = 6.0
    slip = 10.0
    rake = 90.0
    ns = 1
elif T == 2:
    # Initial condition of earthquake for multiple source
    x0 = [40000.0, 40000.0]
    y0 = [45000.0, 35000.0]
    depth = [15.0, 15.0]
    strike = [0.0, 0.0]
    dip = [15.0, 15.0]
    length = [10.0, 10.0]
    width = [6.0, 6.0]
    slip = [10.0, 10.0]
    rake = [90.0, 90.0]
    ns = 2

source = num.array([x0, y0, depth, strike, dip, length, width,
                  slip, rake, num.zeros((ns,))]).transpose()

print source

# Create domain
dx = dy = 1000
L = 80000
W = 80000


# Create topography
def topography(x, y):
    el = -1000
    return el


domain = anuga.rectangular_cross_domain(int(L/dx), int(W/dy), len1=L, len2=W)

domain.set_name('test')
domain.set_quantity('elevation', function=topography, location='centroids')

# Ts = earthquake_tsunami(ns=ns, NSMAX=NSMAX, length=length, width=width,
#                         strike=strike, depth=depth, dip=dip,
#                         xi=x0, yi=y0, z0=0, slip=slip, rake=rake,
#                         domain=domain, verbose=False)


Ts = earthquake_source(
             source=source,
             domain=domain,
             lon_lat_degrees=False,
             lon_before_lat=True,
             utm_zone=None,
             #proj4string=None,
             verbose=False)

domain.set_quantity('stage', function=Ts, location='centroids')


dplotter = anuga.Domain_plotter(domain)

print num.max(dplotter.stage), num.min(dplotter.stage)
print num.max(dplotter.elev), num.min(dplotter.elev)
print num.max(dplotter.depth), num.min(dplotter.depth)


dplotter.plot_stage_frame(vmin=-3.0, vmax=8.0)

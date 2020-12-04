import numpy as num

#from anuga.tsunami_source.tsunami_okada import earthquake_tsunami

#from anuga.tsunami_source.okada_tsunami import earthquake_source

import anuga

"""
Pick the test you want to do; T= 0 test a point source,
T= 1  test single rectangular source, T= 2 test multiple
rectangular sources
"""


# Choose what test to proceed
T = 2

if T == 1:
    # Initial condition of earthquake for multiple source
    params = dict(
        xoff   = 400000, 
        yoff   = 400000,
        depth  = 15000.0, 
        length = 100000.0, 
        width  = 60000.0, 
        slip   = 10.0, 
        opening= 0.0, 
        strike = 180.0, 
        dip    = 15.0, 
        rake   = 90.0,
        nu     = 0.25)
elif T == 2:
    # Initial condition of earthquake for multiple source
    params = dict(
        xoff   =[400000.0, 400000.0], 
        yoff   =[450000.0, 350000.0],
        depth  =[15000.0, 15000.0], 
        length =[10000.0, 10000.0], 
        width  =[6000.0, 6000.0], 
        slip   =[10.0, 10.0], 
        opening=[0.0, 0.0], 
        strike =[0.0, 0.0], 
        dip    =[15.0, 15.0], 
        rake   =[90.0, 90.0],
        nu     =0.25)


# Create domain
dx = dy = 10000
L = 800000
W = 800000


# Create topography
def topography(x, y):
    el = -1000
    return el


domain = anuga.rectangular_cross_domain(int(L/dx), int(W/dy), len1=L, len2=W)

domain.set_name('test')
domain.set_quantity('elevation', function=topography, location='centroids')


def tsunami_function(x,y):
    import okada

    params['x'] = x
    params['y'] = y

    #import pprint
    #pprint.pprint(params)

    ue,un,uz = okada.forward(**params)

    return uz

domain.set_quantity('stage', function=tsunami_function, location='centroids')


dplotter = anuga.Domain_plotter(domain)

print (num.max(dplotter.stage), num.min(dplotter.stage))
print (num.max(dplotter.elev), num.min(dplotter.elev))
print (num.max(dplotter.depth), num.min(dplotter.depth))


dplotter.plot_stage_frame(vmin=-0.2, vmax=0.2)

#!/usr/bin/env python

import os
import sys
import project

import anuga
import numpy

scenario = project.scenario
name = 'Tohoku_Fujii_35299'

print 'output dir:', name
which_var = 0

if which_var == 0:    # Stage
    outname = name + '_stage'
    quantityname = 'stage'

if which_var == 1:    # Absolute Momentum
    outname = name + '_momentum'
    quantityname = '(xmomentum**2 + ymomentum**2)**0.5'    #Absolute momentum

if which_var == 2:    # Depth
    outname = name + '_depth'
    quantityname = 'stage-elevation'  #Depth

if which_var == 3:    # Speed
    outname = name + '_speed'
    quantityname = '(xmomentum**2 + ymomentum**2)**0.5/(stage-elevation+1.e-30)'  #Speed

if which_var == 4:    # Elevation
    outname = name + '_elevation'
    quantityname = 'elevation'  #Elevation

print 'start sww2array'
#origin=[54, ]
#origin=[54, ]
print quantityname
for i in range(0,61):
#    print 'not here!'
#    continue
    
    index = i
    time = 120 * index
    outname = '%s_%s_%04d' % (name,quantityname,time)
    print 80*'='
    print outname,index

    cellsize = 1000.0
    a = anuga.sww2array(name+'.sww',
                  quantity=quantityname,
                  cellsize=cellsize,
                  reduction=index,
                  verbose=True)

    import ipdb as pdb
    pdb.set_trace()

    # produce png file

    import pylab
    import numpy
    a = numpy.where(a == -9999, numpy.nan, a)
    #a = numpy.where(a > 10.0, numpy.nan, a)


    
    a = a[::-1,:]
    nrows = a.shape[0]
    ncols = a.shape[1]

    ratio = float(nrows)/float(ncols)
    #print ratio
    y = numpy.arange(nrows)*cellsize
    x = numpy.arange(ncols)*cellsize

    #Setup fig size to correpond to array size
    fig = pylab.figure(figsize=(10, 10*ratio))

    levels = numpy.arange(-12, 20, 1.0)
    CF = pylab.contourf(x,y,a, levels=levels)
    CB = pylab.colorbar(CF, shrink=0.8, extend='both')
    CC = pylab.contour(x,y,a, levels=[0.0])

    pylab.title('Tohoku %4d secs'% time )
    pylab.savefig(outname+'.png', dpi=600)





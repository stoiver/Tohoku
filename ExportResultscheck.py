#!/usr/bin/env python

import os
import sys
import project

import anuga

scenario = project.scenario
name = 'Fukushima_' + scenario

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

print 'start sww2dem'

print quantityname
## Extract data for different time steps

for i in range(0,2):
#    continue
    index = i*7
    if quantityname=='stage':
       outname = '%s_stage%04d' % (name,10.*index)
    else:
       outname = '%s_elev%04d' % (name,10.*index)

    print name, outname, index
#
for i in range(1,21):
    jndex = 30*i
    if quantityname=='stage':
       outname = '%s_stage%04d' % (name,10.*jndex)
    else:
       outname = '%s_elev%04d' % (name,10.*jndex)

    if jndex < 270:
        namein = name
        index = jndex
    elif jndex < 546:
        namein = name+'_time_2700_0'
        index = jndex-269
    else:
        namein = name+'_time_5460_0'
        index = jndex-545
    
    print namein,outname, index


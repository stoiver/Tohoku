""" Common filenames and locations for topographic data, meshes and outputs.
    This file defines the parameters of the scenario you wish to run.
"""

from __future__ import print_function
import os, sys
from os.path import join
import anuga
from anuga.geometry.polygon import plot_polygons



#------------------------------------------------------------------------------
# Define scenario 
#------------------------------------------------------------------------------
#scenario = 'okada'  # bad okada source
#scenario = 'Fujii' # good fit to dart 21418
#scenario = 'Ammon'
#scenario = 'Hayes' # good fit to dart 21418
#scenario = 'UCSB3'
scenario = 'Caltech'
#------------------------------------------------------------------------------
# Filenames
#------------------------------------------------------------------------------
name_stem = 'Tohoku'
combined_elevation_basename = 'Tohoku'
meshname = name_stem + '_'+scenario+ '_.msh'
topographies_folder = join(os.getcwd(),'topo') 
ascii_grid_filenames = []
point_filenames = ['etopo4min_corrected.txt',
		   'etopo1minNan.txt',
                   'bathy_utmNan.txt']

verbose= True

# Filename for locations where timeseries are to be produced

#gauge_filename = 'gauges.csv'

#------------------------------------------------------------------------------
# Domain definitions
#------------------------------------------------------------------------------
# bounding polygon for study area
#bounding_polygon = anuga.read_polygon('polygons/bounding_n.csv')
bounding_polygon = anuga.read_polygon('polygons/bounding_extended.csv')

A = anuga.polygon_area(bounding_polygon) / 1000000.0
print ('Area of bounding polygon = %.2f km^2' % A)

#------------------------------------------------------------------------------
# Interior region definitions
#------------------------------------------------------------------------------
# Read interior polygons
#poly_1min = anuga.read_polygon('polygons/bounding_1min.csv')
poly_level1 = anuga.read_polygon('polygons/polygon_new1.csv')
poly_level2 = anuga.read_polygon('polygons/bounding_gps.csv')
poly_level3 = anuga.read_polygon('polygons/bounding_inundation.csv')

# Optionally plot points making up these polygons

#plot_polygons([bounding_polygon, poly_level1, poly_level2,poly_level3],figname="plot_polyl.png")

# Define resolutions (max area per triangle) for each polygon

rfact = 30
#rfact = 100
res_whole   =   500000*rfact 
res_level1  =   500000*rfact  
res_level2  =   120000*rfact 
res_level3  =    20000*rfact 

# Define list of interior regions with associated resolutions
build_regions = [poly_level1, poly_level2]
interior_regions = [[poly_level1, res_level1],[poly_level2, res_level2], [poly_level3, res_level3]]
number_triangles=anuga.geometry.polygon.number_mesh_triangles(interior_regions,bounding_polygon,res_whole)

print ('Number of triangles: %s'%number_triangles)

# Define directory

output_folder = os.getcwd() 
output_run = join(output_folder, '_output')+'_'+ scenario
print ('output directory: %s'%output_run) 
source_folder=join(os.getcwd(),'sources')
source_file=join(source_folder,scenario+'.pts')
print ('source directory: %s'%source_file)



lores_cellsize = 20000

#!/usr/bin/env python
"""
Build the elevation data to run a tsunami inundation scenario
for Darwin, NT, Australia.
Input: elevation data from project.py
utput: pts file stored in project.topographies_dir
The run_model.py is reliant on the output of this script.
"""
#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------
# Standard modules
import sys
import os
from os.path import join
import anuga
import project 

# Fine pts file to be clipped to area of interest
#------------------------------------------------------------------------------
print 'project.bounding_polygon', project.bounding_polygon
print 'project.combined_elevation_basename', project.combined_elevation_basename
# Create Geospatial data from ASCII files
geospatial_data = {}
if not project.ascii_grid_filenames == []:
    for filename in project.ascii_grid_filenames:
       	absolute_filename = join(project.topographies_folder, filename)
        # convert_dem_from_ascii2netcdf(absolute_filename,                                    
				# basename_out=absolute_filename,                                      
				# use_cache=True,                                     
				# verbose=True)        
	# dem2pts(absolute_filename, use_cache=True, verbose=True)	
        anuga.asc2dem(absolute_filename+'.asc',
			use_cache=False, 
			verbose=True)
	
        anuga.dem2pts(absolute_filename+'.dem', 
			use_cache=False, 
			verbose=True)
        G_grid = anuga.geospatial_data.Geospatial_data(file_name=absolute_filename+'.pts',verbose=True)
 	print 'Clip geospatial object'
        geospatial_data[filename] = G_grid.clip(project.bounding_polygon)
# Create Geospatial data from TXT files
if not project.point_filenames == []:
   for i in range(0,len(project.point_filenames)):
       filename = project.point_filenames[i]
#    for filename in project.point_filenames:        
       print 'filename:',filename,i
       absolute_filename = project.topographies_folder+'/'+str(filename)       
       G_points = anuga.geospatial_data.Geospatial_data(file_name=absolute_filename, verbose=True)
       print 'Clip %s inside polygon' % filename
       if i==0:
            G_points = G_points.clip(project.bounding_polygon)
       else:
            G_points = G_points.clip(project.interior_regions[i-1][0])
       if filename != project.point_filenames[-1]:
            print 'Clip %s to outside interior polygon' % filename
            print project.interior_regions[i][0]
            G_points = G_points.clip_outside(project.interior_regions[i][0])
       geospatial_data[filename] = G_points
#-------------------------------------------------------------------------------
# Combine, clip and export dataset 
#-------------------------------------------------------------------------------
print 'Add geospatial objects' 
# except', project.offshore_name5
G = None
for key in geospatial_data:
    #if key == project.point_filenames[4] or key == project.ascii_filenames[1]:    
	#    # Skip these files for now    
	#    continue	
	G += geospatial_data[key]
#print 'Clip combined geospatial data'
##G_clip = G.clip_outside(project.poly_aoi1)
##G_all = G_clip + geospatial_data[project.point_filenames[4]] 
#G_clipped = G_all.clip(project.poly_all)
#G_clip = G.clip(project.bounding_polygon)
print 'Export combined DEM file'
G.export_points_file(project.combined_elevation_basename + '.pts')
print 'Do txt version too'
# Use for comparision in ARC
G.export_points_file(project.combined_elevation_basename + '.txt')


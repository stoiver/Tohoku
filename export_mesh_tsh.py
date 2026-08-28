"""Write the notebook mesh to ANUGA's ASCII .tsh format, georeferenced.

The mesh is the one every entry point builds: notebook_tohoku_open_elevation.ipynb,
setup_simulation.create_domain() and calibrate_deterministic.py all call
create_domain_from_regions() with the polygons and resolutions in project.py, so
this reproduces it exactly (354 946 triangles at rfact = 30).

    python export_mesh_tsh.py [outfile.tsh]

The point of the script is the georeference.  Left to itself,
create_mesh_from_regions() resolves the zone to DEFAULT_ZONE (-1, undefined)
because no poly_geo_reference is supplied, and writes that into the file --
the corner is right but nothing records that these are UTM 54N coordinates.
Passing an explicit mesh_geo_reference puts the zone in.

Note what .tsh *cannot* carry.  Its georeference block is exactly three values --
zone, xllcorner, yllcorner (see Geo_reference.write_ASCII).  Hemisphere, EPSG and
false_easting have no slot, so a reader gets zone 54 and must supply the rest.
That matters because Geo_reference defaults the hemisphere to *southern* when it
is undefined: on the way back in, set it explicitly, and prefer
`domain.set_epsg(32654)` over `set_zone`/`set_hemisphere` (see CLAUDE.md,
*Coordinate system*).  Use .msh (NetCDF) if the full georeference must survive.
"""
import sys

import numpy as np

import project
from anuga.coordinate_transforms.geo_reference import Geo_reference
from anuga.pmesh.mesh_interface import create_mesh_from_regions

ZONE = 54                       # UTM zone 54N, EPSG:32654

outfile = sys.argv[1] if len(sys.argv) > 1 else 'Tohoku_notebook_mesh.tsh'

# Same corner create_mesh_from_regions() would derive on its own -- the lower
# left of the bounding polygon -- so only the zone changes.
poly = np.asarray(project.bounding_polygon)
geo = Geo_reference(xllcorner=float(poly[:, 0].min()),
                    yllcorner=float(poly[:, 1].min()),
                    zone=ZONE)

create_mesh_from_regions(
    project.bounding_polygon,
    boundary_tags={'bottom': [0], 'ocean_east': [1], 'top': [2], 'onshore': [3]},
    maximum_triangle_area=project.res_whole,
    interior_regions=project.interior_regions,
    mesh_geo_reference=geo,
    filename=outfile,
    use_cache=False,
    verbose=False,
)

print(f'wrote {outfile}  (zone {ZONE}, '
      f'xll {geo.xllcorner:.1f}, yll {geo.yllcorner:.1f})')

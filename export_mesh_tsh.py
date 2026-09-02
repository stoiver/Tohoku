"""Export the mesh used by notebook_tohoku_open_elevation.ipynb to .tsh (or .msh).

    python export_mesh_tsh.py [outfile.tsh]

The domain is built with the notebook's own call -- same polygons, boundary tags,
resolutions and `set_epsg(32654)` -- and the mesh is then pulled straight off it
with `domain.mesh.save_to_file()`.  Going through the domain rather than calling
create_mesh_from_regions() a second time matters: it is the mesh the runs
actually use, so triangle indices line up with the cells in the `.sww` output.
(`save_mesh_to_tsh()` does the same thing but is deprecated.)

The georeference comes along by itself -- `set_epsg(32654)` puts zone 54 on the
domain, and that is what gets written.

Two things .tsh does not carry:

* **The outline.**  Bounding polygon, interior regions, holes and max areas
  belong to the mesh generator, not the runtime mesh, so that section is empty
  (see Mesh.save_to_file).  The file is a triangulation, not a re-generatable
  recipe -- `project.py` remains the recipe.
* **Anything but the zone.**  Geo_reference.write_ASCII emits exactly three
  values: zone, xllcorner, yllcorner.  Hemisphere, EPSG and false_easting have
  no slot, so the object read back reports `hemisphere='undefined'` -- and
  Geo_reference treats undefined as *southern*.  Set it explicitly on the way
  back in (prefer `set_epsg`, see CLAUDE.md *Coordinate system*), or write
  `.msh`, which is NetCDF and keeps the rest.
"""
import sys

import anuga

import project

outfile = sys.argv[1] if len(sys.argv) > 1 else 'Tohoku_notebook_mesh.tsh'

# --- notebook_tohoku_open_elevation.ipynb, cell 10 -------------------------
domain = anuga.create_domain_from_regions(
    project.bounding_polygon,
    boundary_tags={'bottom': [0], 'ocean_east': [1], 'top': [2], 'onshore': [3]},
    maximum_triangle_area=project.res_whole,
    interior_regions=project.interior_regions,
    use_cache=False,
    verbose=False,
)
domain.set_epsg(32654)  # WGS 84 / UTM zone 54N
# ---------------------------------------------------------------------------

domain.mesh.save_to_file(outfile)

print(f'triangles : {len(domain)}')
print(f'geo ref   : {domain.geo_reference}')
print(f'wrote     : {outfile}')
print('note: .tsh keeps the zone only -- hemisphere/EPSG are not persisted, '
      'and the outline section is empty')

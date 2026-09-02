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

**Prefer .msh unless you specifically need ASCII.**  Same mesh, but 11.4 MB
against 24.8 MB (NetCDF binary rather than one text line per vertex and per
triangle, the latter at full 16-digit float repr for a mesh whose cells are
~250 m), ANUGA reads it natively, and it keeps the whole georeference where
.tsh keeps only the zone.  Measured on this mesh: .tsh 24.8 MB, .tsh.gz 9.2,
.tsh.xz 7.2, .msh 11.4, .msh.xz 5.2.  Just pass a .msh filename.

Two things .tsh does not carry:

* **The outline.**  Bounding polygon, interior regions, holes and max areas
  belong to the mesh generator, not the runtime mesh, so that section is empty
  (see Mesh.save_to_file).  The file is a triangulation, not a re-generatable
  recipe -- `project.py` remains the recipe.
* **Anything but the zone.**  Geo_reference.write_ASCII emits exactly three
  values: zone, xllcorner, yllcorner.  Hemisphere, EPSG and false_easting have
  no slot, so the object read back reports `hemisphere='undefined'` -- and
  Geo_reference treats undefined as *southern*.  Set it explicitly on the way
  back in (prefer `set_epsg`, see CLAUDE.md *Coordinate system*).  A .msh of
  this same mesh reads back as `(zone=54, ..., hemisphere=northern,
  epsg=32654)` -- verified, not assumed.
"""
import gzip
import os
import shutil
import sys

import anuga

import project

outfile = sys.argv[1] if len(sys.argv) > 1 else 'Tohoku_notebook_mesh.tsh'

# A trailing .gz is handled here rather than by ANUGA: save_to_file() dispatches
# on the last four characters, so it will not take 'x.tsh.gz'.  Write the plain
# file, compress it, drop the original.  Note nothing in ANUGA reads .gz back --
# gunzip first, or use .msh, which is smaller than .tsh.gz anyway and keeps the
# full georeference.
gz = outfile.endswith('.gz')
inner = outfile[:-3] if gz else outfile

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

domain.mesh.save_to_file(inner)

if gz:
    with open(inner, 'rb') as fin, gzip.open(outfile, 'wb', compresslevel=9) as fout:
        shutil.copyfileobj(fin, fout)
    raw = os.path.getsize(inner)
    os.remove(inner)
    print(f'triangles : {len(domain)}')
    print(f'geo ref   : {domain.geo_reference}')
    print(f'wrote     : {outfile}  '
          f'({os.path.getsize(outfile)/1e6:.1f} MB, from {raw/1e6:.1f} MB raw)')
else:
    print(f'triangles : {len(domain)}')
    print(f'geo ref   : {domain.geo_reference}')
    print(f'wrote     : {outfile}  ({os.path.getsize(outfile)/1e6:.1f} MB)')
print('note: .tsh keeps the zone only -- hemisphere/EPSG are not persisted, '
      'and the outline section is empty')

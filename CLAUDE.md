# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository simulates the 2011 Tohoku tsunami using the [ANUGA](https://github.com/anuga-community/anuga_core) hydrodynamic solver. It models earthquake-generated sea-floor deformation via Okada fault models and propagates the resulting tsunami wave over a triangulated mesh of the Japanese coast.

## Environment

ANUGA is installed editable from `~/anuga_core` into the conda env `anuga_env_3.14`, which `.bashrc` activates by default (currently `anuga 3.3.8.dev`). This repo is a *scenario* repo — it has no build step, no package, and no test suite of its own.

**Always run scripts and notebooks from the repository root.** `project.py` resolves `topo/`, `sources/` and the output directory with `os.getcwd()`, so running from elsewhere silently produces wrong paths. The notebooks begin with a `try: os.chdir('Tohoku')` guard for the case where Jupyter was started one level up.

## Running Simulations

**Okada KL source** — the script is MPI-aware (`distribute`/`finalize`), so the same file runs either way:
```bash
python run_Tohoku_okada.py                # single process
mpirun -n 4 python run_Tohoku_okada.py    # parallel, for fine meshes
```
Parallel runs write per-process `.sww` files that `domain.sww_merge(delete_old=True)` combines at the end.

**Pre-computed `.pts` source:**
```bash
python run_Tohoku.py                      # uses project.scenario, e.g. Caltech
```

**Notebook workflow** (preferred for interactive work):
```bash
jupyter notebook notebook_tohoku_source_example.ipynb
```

**Key notebooks:**
- `notebook_tohoku_source_example.ipynb` — end-to-end simulation with KL slip field, driven by `setup_simulation.py`
- `notebook_tohoku_open_elevation.ipynb` — same simulation but with the freely downloadable NOAA DEM in place of `Tohoku.pts` (see *Open elevation data*); the only fully reproducible path without the proprietary data
- `notebook_tohoku_okada_example.ipynb` — Okada source model example
- `notebook_tohoku_okada_kl_example.ipynb` — Okada + KL slip field
- `notebook_random_slipfields.ipynb` / `notebook_okada_kl_test.ipynb` — KL slip-field generation and testing in isolation (no ANUGA run)
- `notebook_plot_sww.ipynb` — visualise `.sww` output files
- `notebook_simple_runup_example.ipynb` — basic ANUGA runup demo
- `notebook_okada_tutorial.ipynb` — tutorial on Okada deformation model

## Architecture

### Configuration (`project.py`)
Central config module imported by all simulation scripts. There is no CLI or config file — **you change a scenario by editing the module-level variables and re-running**, and `project.py` prints polygon area / triangle count as an import side effect. Controls:
- `scenario` — which earthquake source to use (`Caltech`, `Fujii`, `Ammon`, `Hayes`, `UCSB3`, `okada`); selected by uncommenting one of the assignments
- `rfact` — inverse mesh-resolution factor scaling all four `res_*` triangle areas. Larger `rfact` = coarser. Triangle counts depend on the `res_level3` divisor as well: with the historical `res_level3 = 20000*rfact`, `rfact=30`→~87 000, `rfact=60`→~44 000, `rfact=100`→~24 800. `res_level3` is now `2000*rfact` (10× finer inundation region) with `rfact = 30` (current) → ~355 000 triangles, so a full run takes minutes rather than seconds; `rfact = 100` is still there commented out for a fast path
- `interior_regions` — nested refinement polygons `poly_level1/2/3` read from `polygons/*.csv`
- `meshname` = `Tohoku_<scenario>_.msh`, `output_run` = `_output_<scenario>`, `source_file` = `sources/<scenario>.pts`

### Simulation pipeline

1. **Mesh creation** — `anuga.create_domain_from_regions()` builds a triangular mesh from boundary polygons in `polygons/`. Mesh resolution is refined in nested interior polygons (`poly_level1/2/3`).

2. **Elevation setup** — `Tohoku.pts` (combined bathymetry/topography point cloud) is interpolated onto mesh centroids.

3. **Earthquake source** — Two approaches:
   - *Pre-computed `.pts` sources* (`run_Tohoku.py`): reads deformation grids from `sources/<scenario>.pts` and adds them directly to the stage.
   - *KL random slip field* (`run_Tohoku_okada.py`, `setup_simulation.py`): generates a stochastic slip field using Karhunen–Loève expansion then computes Okada surface deformation.

4. **Evolution** — ANUGA integrates shallow-water equations; results saved as `.sww` binary files.

### Key modules

| File | Purpose |
|------|---------|
| `okada.py` | Okada (1985) elastic dislocation model — computes surface displacements (uE, uN, uZ) for a rectangular fault |
| `okada_subfaults.py` | Divides a fault into sub-faults and sums their Okada displacements |
| `okada_kl_subfaults.py` | Extends `okada_subfaults` with a KL-based random slip field for uncertainty quantification; `kl_deformation()` is the main entry point |
| `setup_simulation.py` | Higher-level helpers: `create_domain()`, `apply_deformation()`, `evolve_domain()`, and gauge recording classes |
| `project.py` | Scenario parameters, mesh polygons, resolution settings |
| `tsunami_observations.py` | Loads the TTJS surveyed inundation/run-up heights for validation; subsetting by UTM extent or by gauge |

### Coordinate system
All simulations use **UTM Zone 54N**, set with `domain.set_epsg(32654)`. This replaced the older `domain.set_hemisphere('northern')` + `domain.set_zone(54)` pair — use `set_epsg` in new code. The two are equivalent for zone/hemisphere/EPSG, with one difference: `set_epsg` also populates `false_easting` from pyproj (0 → 500000, the correct UTM value), which is metadata written into `.sww` and surfaces as `Xshift` in `.prj` files produced by `sww2dem` (used by `ExportResults.py`). Note the ordering rule if you ever go back to the old calls — `set_zone` defaults the hemisphere to *southern* when it is still undefined, so `set_hemisphere` must come first.

ANUGA stores centroid coordinates *relative* to the domain lower-left corner, so anything expressed in absolute UTM (fault epicentre, gauge positions) must be offset by `domain.geo_reference.xllcorner`/`yllcorner` before use — this is what `setup_simulation.apply_deformation()` does with `xoff`/`yoff`. Gauge lat/lon is converted with `utm.from_latlon(..., force_zone_number=54)`.

### Output
- `.sww` files — ANUGA's binary output format (NetCDF), visualised with `anuga.SWW_plotter`, the plot notebooks, or `anuga_sww_gui` (below).
- `_output_<scenario>/` directories — per-run outputs including `.sww`, slip field PNGs, `slips.txt`, and stage time-series.

## Visualising results

`~/anuga_core/scripts/anuga_sww_gui.py` is a Tk GUI/renderer for `.sww` files, driven by a TOML config. `tohoku_source_example_001.toml` is the worked example for this repo:

```bash
python ~/anuga_core/scripts/anuga_sww_gui.py --config tohoku_source_example_001.toml
```

It renders `stage` over an OpenStreetMap/satellite basemap (`epsg = 32654` matching UTM 54N) and writes numbered PNG frames to `_plot/`. The `*_osm.jpg` / `*_satellite.jpg` files with their `.georef` sidecars are exported basemap tiles; a `.georef` line is `xll yll xres yres width height`.

## Gotchas

- **`iseed` in `run_Tohoku_okada.py` does not select the slip realisation.** Despite its comment, the top-level `iseed` only feeds the run name (`Okada_<iseed>` → `_output_Okada_<iseed>`); the actual KL draw is fixed by the hard-coded `iseed=1001` in the `okl.kl_deformation(...)` call. `setup_simulation.apply_deformation()` similarly hard-codes `iseed=1234`. Change both if you want a genuinely different realisation.
- **Mesh caching differs between entry points.** `run_Tohoku_okada.py` calls `create_domain_from_regions(..., use_cache=True)` while `setup_simulation.create_domain()` uses `use_cache=False`; both write to the same `project.meshname`. After changing `rfact` or the polygons, delete `Tohoku_<scenario>_.msh` to be sure the mesh is rebuilt.
- **Different scripts run different durations.** `run_Tohoku_okada.py` evolves 4 hours at a 5-minute yieldstep with `tide = -0.45`; `setup_simulation.evolve_domain()` evolves 2 hours at 2 minutes with `tide = 0.0`. Don't compare their outputs without accounting for this.
- **`evolve_domain()` expects specific gauge keys** — it iterates `[21418, 0, 1, 2]`, so the `gauges` dict passed in must contain all four.
- **Generated artefacts are gitignored, input data is not.** `.gitignore` covers `*.sww`, `*.msh`, `anuga_*.log`, `_output_*/`, `_plot/`, `screenshots/`, `*.tif`, `*.georef` and the `tohoku_open_dem*.jpg` basemaps. It deliberately does **not** ignore `*.pts` — `Tohoku.pts` and `sources/*.pts` are tracked inputs. Don't add `*.pts` to it.

## Dependencies

`anuga>=3.3`, `numpy`, `scipy`, `utm`, `matplotlib`. The open-elevation path additionally needs `rasterio`, `pyproj` and `requests`.

The ANUGA source is at `~/anuga_core`; its `CLAUDE.md` documents the build system, test commands and full architecture, and is the place to look when behaviour needs tracing into the solver. Rebuild after editing ANUGA's C/Cython:

```bash
cd ~/anuga_core
pip install --no-build-isolation -v -e .
pytest --pyargs anuga --run-fast   # ~40 s, skips MPI/slow tests
pytest --pyargs anuga              # full suite, ~3 min
```

## Data Files

- `Tohoku.pts` — combined bathymetry/topography, required by every script except the open-elevation notebook. Tracked in git at **31 MB** (NetCDF binary, so it barely compresses) and by far the largest object in the repo; it dominates clone time. Leave it alone — it is deliberately *not* in `.gitignore`, and removing it now would need a history rewrite.
- `polygons/*.csv` — boundary and interior mesh polygons in UTM coordinates. Only four are live in `project.py` (`bounding_extended`, `polygon_new1`, `bounding_gps`, `bounding_inundation`); the rest are historical alternatives.
- `sources/*.pts` — pre-computed earthquake source deformation fields, one per scenario
- `21418_notide.txt` — de-tided DART buoy 21418 observations used to validate the modelled wave
- `observations/ttjs_survey_20121229.csv` — the 2011 Tohoku Earthquake Tsunami Joint Survey (TTJS) nationwide field survey: 5907 measured inundation/run-up heights, tide-corrected, 31.5-43.7 N. Downloaded from <http://www.coastal.jp/ttjt/> and converted Shift-JIS -> UTF-8. Loaded by `tsunami_observations.py` (`load_ttjs`, `in_extent`, `near`), which also re-downloads it and the NOAA NCEI run-up records. The survey group's terms require citing them and the release date — see `observations/README.md`. Used by the validation section of `notebook_tohoku_open_elevation.ipynb`.
- `build_elevation.py` — legacy Python 2 script for building `Tohoku.pts` from raw topo files. **Does not run** (`SyntaxError` on Python 2 `print` statements) and its docstring still refers to a Darwin, NT scenario.

### Open elevation data

`download_elevation.py` downloads freely available global bathymetry/topography as a GeoTIFF that can substitute for the proprietary `Tohoku.pts`:

```python
from download_elevation import download_noaa_dem

# Download GEBCO-quality DEM for the project domain (no API key needed)
download_noaa_dem('Tohoku_dem.tif')                    # bounds auto-derived from project.py
download_noaa_dem('Tohoku_dem.tif', resolution_arcsec=30)  # coarser/faster download

# Use in a simulation instead of Tohoku.pts:
domain.set_quantity('elevation', filename='Tohoku_dem.tif', location='centroids')
```

Source: NOAA NCEI Global DEM Mosaic (GEBCO-based), EPSG:4326, 15 arc-second native resolution, land + ocean. An OpenTopography alternative (free API key required) is also provided in the same file.

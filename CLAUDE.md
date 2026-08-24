# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository simulates the 2011 Tohoku tsunami using the [ANUGA](https://github.com/anuga-community/anuga_core) hydrodynamic solver. It models earthquake-generated sea-floor deformation via Okada fault models and propagates the resulting tsunami wave over a triangulated mesh of the Japanese coast.

## Environment

ANUGA is installed editable from `~/anuga_core` into the conda env `anuga_env_3.14`, which `.bashrc` activates by default (currently `anuga 3.3.8.dev`). This repo is a *scenario* repo — it has no build step and no package.

It does have a small test suite for the source-model code, which needs only `numpy`, `scipy`, `matplotlib` and `pytest` — **no ANUGA** — so it runs anywhere in a couple of seconds:

```bash
python -m pytest tests -v
```

`tests/test_okada.py` covers the Okada dislocation model (linearity in slip, zero-slip, decay with distance, subfault summation against the single-fault solution); `tests/test_okada_kl.py` covers the KL slip field, including regression guards for two bugs that were live in the repo — see *Gotchas*. `.github/workflows/tests.yml` runs them on push and PR across ubuntu-latest and macos-latest × Python 3.11/3.12/3.13. **The Python spread is what gives the complex-eigenvalue guard its teeth**, because that bug tracks the numpy version rather than the OS: a CI probe with the bug reinstated failed on numpy 2.5.2 under both OpenBLAS and Accelerate and passed on 2.4.6 under both. Keep 3.12 in the matrix — that is where it was first hit in the wild. Every job prints `numpy.show_config()`, which is the first thing worth reading when a numerical failure appears on one leg only.

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

## Calibration against the TTJS survey

Scored with Aida's (1978) *K* and *&kappa;* over the ~1700 survey points in the inundation close-up box, alongside the DART 21418 peak (observed 1.87 m at 33 min). *K* is the geometric mean of observed/modelled, so **K = 1 is unbiased and K < 1 means the model runs high**; *&kappa;* is the geometric standard deviation, i.e. the typical scatter factor. `tsunami_observations.py` loads the survey; the notebook's validation section does the scoring.

### Notebook defaults

`notebook_tohoku_open_elevation.ipynb` ships `flow_algorithm = 'DE_ader2'`, `slip = 84`, `friction = 0.033`, `elevation_source = 'open'`: DART peak 1.88 m against 1.87 m observed, **K = 1.00**, κ = 1.70, bias +0.20 m, 68 of ~1700 survey points left dry. K = 0.95 is inside the Japanese guideline band (0.95 < K < 1.05); κ is not, and no lever tested so far moves it — see the κ floor below.

**Calibrate in two stages, in this order — the knobs are orthogonal.** The modelled DART peak is linear in slip (0.023 m of peak per metre of slip) and friction does not touch it at all: at slip 81 the peak is 1.86 m for n = 0.03, 0.035 and 0.04 alike. So

1. set `slip` from the DART 21418 peak, then
2. set `friction` from Aida K against the survey.

Under DE_ader2 at slip 84: n = 0.03 → K = 0.93; n = 0.033 → K = 1.00, κ = 1.70, bias +0.20 m; n = 0.035 → K = 1.06. (Under DE0 at slip 81 the same sweep gave n = 0.03 → K = 0.95, n = 0.035 → K = 1.08, n = 0.04 → K = 1.24.) **The calibrated pair is solver-specific** — changing the flow algorithm shifts slip by a few per cent and friction by ~10%, so recalibrate if you change it. K and the arithmetic bias optimise at slightly different n because K is a geometric mean of ratios (weighting the many small-height points) and the bias is an arithmetic mean of differences (weighting the few large ones). K is the conventional measure, so 0.03 is the shipped choice.

These values are tied to the fixed KL code (`eigh`, normal Sobol deviates) — see *Gotchas*. Before those fixes the notebook used slip 60 with n = 0.04; that pairing is meaningless now, since the old uniform-deviate bug inflated the slip field and the same settings give K = 1.54 with a DART peak of 1.38 m against the corrected code.

### Best configuration found

**UCSB3 source + `Tohoku.pts` elevation + Manning n = 0.05**, on the full mesh:

| run | DART | K | κ | bias | RMS | survey pts left dry |
|-----|------|------|------|-------|------|----|
| UCSB3, open DEM, n = 0.04 | 1.86 | 0.95 | 1.72 | +0.47 | 3.47 | 64 |
| UCSB3, `Tohoku.pts`, n = 0.04 | 1.85 | 0.84 | 1.70 | +0.97 | 3.92 | **3** |
| **UCSB3, `Tohoku.pts`, n = 0.05** | **1.85** | **1.06** | 1.71 | **−0.04** | **3.74** | 5 |
| KL slip 81, open DEM, n = 0.03 (what the notebook ships) | 1.86 | 0.95 | 1.69 | +0.45 | 4.06 | 68 |

The notebook's own configuration is comparable on K and κ, but misses 68 surveyed points entirely; the best configuration misses 5. The gap is inundation extent, not height accuracy — which is a property of the DEM, not the source.

The UCSB3 rows predate the KL fixes but are unaffected by them: those runs take their deformation from `sources/UCSB3.pts` and never touch the KL code. The KL row is post-fix.

Two things to carry away:

- **Friction and DEM are coupled — tune them together.** n = 0.04 is right for the 450 m open DEM and n = 0.05 for the 150 m `Tohoku.pts`: the finer bathymetry lets more water through, so it wants more roughness. Retuning n on the finer DEM moved K from 0.84 to 1.06 and the bias from +0.97 m to −0.04 m while costing only two extra dry points.
- **κ ≈ 1.7 is a floor, and it is *not* the mesh.** Every configuration tried — three sources, two DEMs at 450 m and 150 m, friction from 0 to 0.05, one and two fault segments, before and after the KL fixes — bottoms out at κ = 1.67–1.73, and only ever goes *up* from there (to 2.0 at badly chosen friction). This was assumed to be discretisation until it was measured directly: see *Mesh resolution* below. It is not — nor is it the solver (*Flow algorithm*), nor missing sea defences (*Coastal defences*), each tested and eliminated. Do not expect a finer mesh to buy you the κ < 1.45 guideline target.

### Flow algorithm

`DE0` is the ANUGA default and what this repo used until now.  Measured at fixed source, DEM, friction and mesh (UCSB3 + `Tohoku.pts` + n = 0.05, ~250 m):

| algorithm | timestepping | runtime | DART | K | κ | RMS |
|-----------|--------------|---------|------|------|------|------|
| `DE0` | euler | 88 s | 1.85 | 1.06 | 1.71 | 3.74 |
| `DE1` | rk2 | 129 s | 2.02 | 0.93 | 1.73 | 5.09 |
| **`DE_ader2`** | ader2 | **88 s** | 1.81 | 1.07 | **1.66** | **3.63** |

`DE_ader2` is a free improvement over the default: better κ and RMS at *identical* runtime, because its higher-order reconstruction permits a larger stable timestep that pays for the extra work per step.  The notebook now uses it.

**`DE1` is the odd one out and worth treating with suspicion** — second-order time integration made the fit markedly *worse* (DART 2.02 m against 1.87 observed, RMS 5.09) at 1.5× the cost.  That may be an interaction with the limiter betas `_set_DE1_defaults()` installs, which differ from DE0's; it has not been chased down.

All three run under `set_compute_mode('unified')` with no fallback to legacy.  `set_flow_algorithm()` must be called *before* `set_compute_mode()` and before any quantities are set, since it resets the limiter betas and the timestepping method.

Note the κ column: 1.66–1.73 across all three.  The solver is not what sets the scatter either — see *Mesh resolution*.

### Coastal defences (riverwalls)

`Tohoku.pts` cannot resolve the sea defences, and adding them as riverwalls does not help.

**The DEM has no defences in it.** At 147 x 159 m spacing against a Sendai seawall ~30-40 m wide at the base, there is roughly *0.27 of a sample* across the structure, so a wall does not get smoothed — it falls between sample points and vanishes. Measured directly, the maximum elevation within 300 m inland of the 0 m line is 2.30 m at Arahama, 2.24 m at Yuriage and 2.11 m at Sendai airport, against a real crest of ~6.2 m T.P. This is what `friction` is standing in for, and it is why refining the mesh changed nothing (see *Mesh resolution*): the missing feature is in the data, not the discretisation, and a 62 m mesh still cannot see a 30 m wall.

**ANUGA's riverwalls are the right tool and they work here.** A riverwall is a sub-grid barrier carried on mesh edges with a weir flux, so it represents a 30 m wall on a 250 m mesh. They are implemented in the GPU kernels (`gpu_domain_core.c`, `core_kernels.c`), they run under `set_compute_mode('unified')` with `DE_ader2`, and they require a DE algorithm. The wall alignment must coincide *exactly* with mesh edges: pass it to `breaklines` in `create_domain_from_regions`, then call `domain.riverwallData.create_riverwalls()`, which handles the georeference itself (pass absolute UTM).

**But the intact-wall idealisation is wrong.** A 6.2 m crest along 33.6 km of the Sendai coast (alignment derived from the DEM's 0 m contour, offset 150 m landward, 230 edges instantiated), paired against an identical run on the identical breakline mesh:

| Sendai box, 953 points | K | κ | bias | RMS | onshore wet cells |
|------------------------|------|------|-------|------|-----|
| no wall | 1.17 | 1.61 | −0.78 | 1.97 | 4545 |
| 6.2 m riverwall | 1.94 | 1.97 | −2.02 | 3.04 | 3519 |

Every measure gets worse: 23% less onshore wetting, K to 1.94, κ *up* from 1.61 to 1.97. The real Sendai levees were overtopped within minutes and largely destroyed, so the survey heights record what happened *after* they failed; a wall that stands for four hours blocks water that actually went inland. This is also evidence against the natural hypothesis that missing defences cause the κ scatter — adding them makes the scatter worse.

A fair test needs **failure**: crest held at 6.2 m until overtopping, then dropped. ANUGA riverwalls are static, but an operator could lower `riverwall_elevation` on a time or overtopping criterion. That needs per-segment failure information which does not appear to exist at the required resolution.

One trap worth recording: the first attempt at this pairing inherited `n = 0.033` (calibrated on the *open* DEM) onto the `pts` DEM, leaving it badly under-frictioned at κ = 2.81 rather than ~1.6. At that operating point the wall's effect was partly buried — κ appeared flat (Δ = −0.002) and the K shift was milder. **Pair experiments at a calibrated operating point**, or the comparison measures the mis-calibration instead.

### Mesh resolution

Refining the mesh does **not** reduce κ.  Same source, DEM and friction; only `res_level3` (the inundation region) refined, from ~250 m cells to ~62 m:

| mesh | cells | DART | K | κ | bias | RMS | dry |
|------|-------|------|------|------|-------|------|-----|
| ~250 m | 354 946 | 1.85 | 1.06 | 1.71 | −0.04 | 3.74 | 5 |
| ~62 m | 4 735 001 | 1.82 | 1.02 | **1.73** | +0.37 | 4.15 | 2 |

**13× the cells and κ does not move** (1.71 → 1.73, marginally worse).  That disproves the earlier assumption — recorded here as fact for some time — that the ~250 m discretisation was what set the scatter.  It is not.

What remains: at 62 m cells the model is reading a 150 m DEM, so **the bathymetry is now the binding constraint** rather than the mesh; or the residual is genuine bay-to-bay variability in the survey that no model at this scale reproduces.  The first is testable with a finer DEM over a small area; the second is not fixable at all.

The secondary effects are small and in the expected directions: K improves 1.06 → 1.02 and the dry count falls 5 → 2 as finer cells resolve more low ground, while bias and RMS drift positive because n = 0.05 was calibrated on the 250 m mesh — the friction coupling again.

Practical notes for repeating it: the run was 4.7 M triangles, 91 minutes on the GPU, 1.9 GB of `.sww`.  To keep the file that size the run stored only `stage` and `elevation` (no momenta) and used a 2-minute yieldstep for the first 90 simulated minutes then 10 minutes thereafter — K and κ are read from the `max_*_c` fields, which the operator writes every timestep regardless of yieldstep.  Two traps in that:

- `set_quantities_to_be_stored()` must be called **before** `set_collect_max_quantities()`.  The operator *augments* the store list at construction, so the other order silently drops every `max_*` field and leaves a file that cannot be scored.
- `anuga.SWW_plotter` requires `xmomentum_c` and will raise `KeyError` on such a file.  Read those runs with `netCDF4` directly, or keep the momenta.

### Friction

All rows below are the same run — KL single plane, slip 81, open DEM, full mesh, fixed KL code — with only Manning *n* changed:

| n | DART peak | K | κ | bias | RMS | survey pts left dry |
|------|------|------|------|--------|------|-----|
| 0     | 1.86 | 0.63 | 1.87 | +2.61 | 5.31 | 59 |
| 0.02  | 1.86 | 0.79 | 2.01 | +1.46 | 4.43 | 62 |
| 0.025 | 1.86 | 0.85 | 1.72 | +0.97 | 4.20 | 65 |
| **0.03** | 1.86 | **0.95** | **1.69** | +0.45 | 4.06 | 68 |
| 0.035 | 1.86 | 1.08 | 1.72 | −0.07 | 3.99 | 68 |
| 0.04  | 1.86 | 1.24 | 1.79 | −0.56 | 3.98 | 69 |
| 0.05  | 1.86 | 1.62 | 1.92 | −1.42 | 4.13 | 78 |

Reading it:

- **The DART peak is 1.86 m at every value, including n = 0.** Friction is inert in the far field, so slip and friction can be calibrated independently and in that order.
- **n = 0.03 is the K optimum** (K = 0.95, inside the guideline band); the arithmetic bias zeroes near 0.035 and RMS bottoms at 0.04. The three criteria disagree across a 0.01 window, which is about the tuning precision this model supports — don't read more than ±0.005 into it.
- **κ bottoms out at n = 0.03 (1.69), the same place as K.** The κ = 2.01 at n = 0.02 breaks the pattern, sitting above both its neighbours; unexplained, possibly a partial-wetting effect at that damping, and left in rather than smoothed away.
- **`dry` rises monotonically with n** — height accuracy is bought with inundation extent. That trades directly against the DEM choice, where `Tohoku.pts` takes dry from 64 to 3, so pair the two rather than tuning either alone.

An earlier version of this table was measured at slip 60 before the KL fixes (`eigh`, normal Sobol deviates). The shape was the same but the optimum sat at n = 0.04; those numbers are superseded.

### Sources

Scored on a coarse mesh (all `res_*` × 8) unless noted, with the deformation added to both `stage` and `elevation`. Note `run_Tohoku.py` adds to `stage` only, so its results will not line up with these:

| source | uZ max | DART | K | note |
|--------|--------|------|------|------|
| Caltech | 6.9 m | 0.87 | 0.97 | best coastal fit of the .pts sources, DART half-size |
| Fujii | 12.5 m | 0.90 | 0.88 | |
| Ammon | 9.5 m | 0.83 | 0.51 | |
| Hayes | 8.6 m | 0.88 | 2.00 | under-predicts the coast 2× |
| UCSB3 | 16.0 m | **1.51** | 0.77 | best joint fit; at full mesh with n = 0.04: DART 1.86, K 0.95 |
| okada.pts | 30.5 m | 4.13 | 0.20 | unusable, as `project.py` already notes |

Fault-parameter sweeps (depth 9–20 km, width 50–120 km, length 200–500 km, and a two-segment source with a shallow near-trench strip) did **not** decouple the far field from the coast: every configuration strong enough to match DART over-predicted the survey ~2×, and enlarging the fault raised coastal height faster than the DART peak. The single-plane 200 × 50 km at 20 km depth remains the best joint geometry. Friction, not source geometry, is what resolves the tension — which is why the friction table above matters more than any of the source variants.

## Gotchas

- **`iseed` in `run_Tohoku_okada.py` does not select the slip realisation.** Despite its comment, the top-level `iseed` only feeds the run name (`Okada_<iseed>` → `_output_Okada_<iseed>`); the actual KL draw is fixed by the hard-coded `iseed=1001` in the `okl.kl_deformation(...)` call. `setup_simulation.apply_deformation()` similarly hard-codes `iseed=1234`. Change both if you want a genuinely different realisation.
- **Mesh caching differs between entry points.** `run_Tohoku_okada.py` calls `create_domain_from_regions(..., use_cache=True)` while `setup_simulation.create_domain()` uses `use_cache=False`; both write to the same `project.meshname`. After changing `rfact` or the polygons, delete `Tohoku_<scenario>_.msh` to be sure the mesh is rebuilt.
- **Different scripts run different durations.** `run_Tohoku_okada.py` evolves 4 hours at a 5-minute yieldstep with `tide = -0.45`; `setup_simulation.evolve_domain()` evolves 2 hours at 2 minutes with `tide = 0.0`. Don't compare their outputs without accounting for this.
- **`evolve_domain()` expects specific gauge keys** — it iterates `[21418, 0, 1, 2]`, so the `gauges` dict passed in must contain all four.
- **ANUGA's default Manning friction is 0.0, and most scripts here never set it.** That was the single largest error in the model: with nothing dissipating the wave between the shelf break and the shore, *every* source tried — the KL plane, a two-segment KL source, and all five published inversions in `sources/` — over-predicted the surveyed inundation by roughly a factor of two once it was strong enough to match DART 21418. `notebook_tohoku_open_elevation.ipynb` now sets `domain.set_quantity('friction', 0.04)`, which takes Aida K from 0.55 to 1.08 and leaves the DART peak unchanged to 0.01 m (friction is negligible in 5700 m of water). Treat 0.04 as an effective value standing in for unresolved bed roughness and coastal defences on a 450 m DEM, not a measured one — and note the calibration is coupled to the source: at n = 0.04, `slip = 60` gives K = 1.08 while `slip = 35` under-predicts at K = 1.58.
- **The KL slip field had two bugs, now fixed — don't reintroduce them.** `kl_correlation_matrices` used `np.linalg.eig` on a symmetric covariance matrix; that is the general non-symmetric LAPACK routine, which may return complex-conjugate eigenpairs, and the complex values propagate through `sqrtD` into the slip field and `okada()`. **It depends on the numpy version, not the platform**: with `eig` restored, CI fails on numpy 2.5.2 under both OpenBLAS and Accelerate, and passes on numpy 2.4.6 under both (run 32253704090). It was first hit on macOS with Python 3.12, which made it look platform-specific; it is not. Use `np.linalg.eigh`, and clip eigenvalues at 0 before the sqrt. Separately, `sample='sobol'` fed raw Sobol points — uniform on [0,1) — into an expansion that wants standard normals, giving every mode a coefficient with mean 0.5 instead of 0; this inflated the slip field, which is why the KL source needed `slip = 60` to match DART while published inversions peak at 7–16 m of uplift. Both are covered by `tests/test_okada_kl.py`.
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

#### Choosing between the two DEMs

`notebook_tohoku_open_elevation.ipynb` takes either, via `elevation_source` in the elevation cell:

- `'open'` — the downloaded GeoTIFF alone (15 arcsec, ~450 m). The reproducible path, and the default.
- `'pts'` — `Tohoku.pts` (1.34 M points, ~150 m) interpolated onto the mesh with `LinearNDInterpolator` and overlaid where it has coverage, with the GeoTIFF filling the rest. `Tohoku.pts` spans x 422–899 km, y 3954–4518 km while the domain runs east to x 1300 km, so it *cannot* supply the deep ocean on its own — fitting it across the whole domain extrapolates to nonsense (this is also why `add_quantity('stage', filename='sources/<scenario>.pts')` in `run_Tohoku.py` needs care: those source fields cover an even smaller box). The run writes `tohoku_pts_dem.sww` rather than `tohoku_open_dem.sww`.

The two agree on the mean surface to +0.03 m over the 330 918 shared cells, but differ by 4.9 m RMS on the shelf, 7.3 m in the −20…0 m band and 12.2 m on coastal land — i.e. wherever the elevation actually decides the inundation. The measured effect, holding source and friction fixed, is almost entirely on *extent* rather than height: heights barely move (K 0.95 → 0.84, κ 1.72 → **1.70**) but the number of surveyed points the model leaves dry falls from **64 to 3**. The 450 m DEM cannot resolve the low coastal strip the water crossed.

Two things worth knowing before tuning further: κ is essentially unchanged by tripling the bathymetric resolution, and equally unchanged by refining the mesh (see *Mesh resolution*), so neither is what sets the point-to-point scatter; and friction and DEM are coupled — n = 0.04 was calibrated on the open DEM, so the finer bathymetry lets more water through and wants a slightly larger n to recover K ≈ 1.

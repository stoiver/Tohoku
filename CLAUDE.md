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
| `okada_kl_subfaults.py` | Extends `okada_subfaults` with two slip fields: a KL random draw (`kl_deformation()`) and a smooth moment-normalised taper (`deterministic_slip()` / `deterministic_deformation()`). Both share `sum_subfault_deformation()` for the Okada summation |
| `calibrate_deterministic.py` | Scriptable one-run calibration driver for either source (`--source det|kl`): builds the notebook's domain, applies the source, evolves, and appends DART peak / Aida K / &kappa; / bias / RMS / dry count to `calibration_results.jsonl` |
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

`notebook_tohoku_open_elevation.ipynb` ships `flow_algorithm = 'DE_ader2'`, `slip = 71`, `friction = 0.038`, `elevation_source = 'open'`. It previously shipped `slip = 84` / `friction = 0.033`, which were calibrated when the KL field's `alpha` was 0.75, which inflated the realised slip (see *Gotchas*). At the corrected `alpha = 0.4` the same nominal `slip = 84` raises the realised mean from 50.1 m to 65.9 m, a 1.32x stronger source.

**The recalibrated pair is `slip = 71`, `friction = 0.038`** — DART 1.87 m against 1.87 m observed, K = 1.01, κ = 1.76, bias +0.39 m, 67 of ~1700 survey points left dry. Measured on the full mesh under DE_ader2 with the open DEM:

| slip | n | DART | K | κ | bias | RMS | dry |
|------|-------|------|------|------|-------|------|-----|
| 64 | 0.033 | 1.68 | 0.95 | 1.70 | +0.50 | 4.05 | 68 |
| 71 | 0.033 | 1.87 | 0.86 | **1.69** | +1.02 | 4.51 | 66 |
| 71 | 0.036 | 1.87 | 0.95 | 1.72 | +0.63 | 4.41 | 66 |
| **71** | **0.038** | **1.87** | **1.01** | 1.76 | +0.39 | 4.36 | 67 |
| 71 | 0.040 | 1.87 | 1.08 | 1.88 | +0.15 | 4.33 | 69 |

Note that **K and κ no longer co-optimise**. Before the `alpha` fix they coincided at n = 0.033 (K 1.00, κ 1.70); now κ bottoms at 1.69 at n = 0.033 where K is only 0.86, and reaching K = 1 costs 0.07 of κ. The deterministic source splits the same way (see *Deterministic source*), so the old coincidence looks like an artefact of the negative-slip field rather than a property worth chasing. K is the conventional measure, so 0.038 is the recommended choice; take n = 0.036 (K 0.95, κ 1.72) if κ matters more for your purpose — both sit inside the Japanese guideline band 0.95 < K < 1.05.

Also note the KL slip scaling is **not** simply the ratio of realised mean slips: predicting from that gave slip 64 and a DART peak of 1.68 m, because `alpha` reshapes the field and DART responds to peak uplift as well as mean strength. Measure the sensitivity instead — it is 0.0263 m of DART peak per metre of nominal slip here.

**Calibrate in two stages, in this order — the knobs are orthogonal.** The modelled DART peak is linear in slip (0.023 m of peak per metre of slip) and friction does not touch it at all: at slip 81 the peak is 1.86 m for n = 0.03, 0.035 and 0.04 alike. So

1. set `slip` from the DART 21418 peak, then
2. set `friction` from Aida K against the survey.

Under DE_ader2 at slip 84: n = 0.03 → K = 0.93; n = 0.033 → K = 1.00, κ = 1.70, bias +0.20 m; n = 0.035 → K = 1.06. (Under DE0 at slip 81 the same sweep gave n = 0.03 → K = 0.95, n = 0.035 → K = 1.08, n = 0.04 → K = 1.24.) **The calibrated pair is solver-specific** — changing the flow algorithm shifts slip by a few per cent and friction by ~10%, so recalibrate if you change it. K and the arithmetic bias optimise at slightly different n because K is a geometric mean of ratios (weighting the many small-height points) and the bias is an arithmetic mean of differences (weighting the few large ones). K is the conventional measure, so 0.03 is the shipped choice. **All the numbers in this subsection predate the `alpha` fix and are superseded** by the recalibration table above; the reasoning about stage ordering and solver-specificity still stands.

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
- **κ ≈ 1.7 is a floor, and it is *not* the mesh.** *(The mechanism is now known: κ is mathematically invariant under any uniform rescaling of the modelled heights, so no strength knob can move it — see* &kappa; *is scale-invariant, below.)* Every configuration tried — three sources, two DEMs at 450 m and 150 m, friction from 0 to 0.05, one and two fault segments, before and after the KL fixes — bottoms out at κ = 1.67–1.73, and only ever goes *up* from there (to 2.0 at badly chosen friction). This was assumed to be discretisation until it was measured directly: see *Mesh resolution* below. It is not — nor is it the solver (*Flow algorithm*), nor missing sea defences (*Coastal defences*), each tested and eliminated. Do not expect a finer mesh to buy you the κ < 1.45 guideline target.

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

An earlier version of this table was measured at slip 60 before the KL fixes (`eigh`, normal Sobol deviates). The shape was the same but the optimum sat at n = 0.04; those numbers are superseded. **This table is itself now superseded** — it predates the `alpha` 0.75 -> 0.4 fix. Its shape survives (κ bottoming near the K optimum, `dry` rising monotonically with n, DART flat in n), but the operating point has moved: see the recalibration table under *Notebook defaults*.

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

Those sweeps varied the *plane* with uniform slip, which changes moment and slip distribution together. The distribution *within* the plane is a separate lever, tested later — see *Deterministic source* below. It does decouple the two (narrowing the width raised the DART peak 36% while K moved 4%), but not usefully: every configuration that improved the far-field/coast ratio degraded &kappa;.

### Deterministic source

`deterministic_slip()` replaces the KL random draw with a smooth reproducible
slip field: a separable Gaussian taper over the plane, scaled so the fault
carries a prescribed seismic moment.

```python
uE, uN, uZ, slips = okl.deterministic_deformation(
    x, y, xoff=x0 - xll, yoff=y0 - yll,
    M0=3.8e22, u0=0.5, sig_u=0.15, v0=0.25, sig_v=0.30,
    depth=23*km, length=400*km, width=150*km)
```

`u` is the along-strike coordinate and `v` the down-dip one, both normalised to
(0, 1) with **v = 0 at the trench**; `u0`/`v0` place the slip peak and
`sig_u`/`sig_v` set how tightly it concentrates. Amplitude is *not* a free
parameter — it follows from `M0 = rigidity * area * mean(slip)`, so retapering
redistributes slip without touching the moment.

**Subfault indexing is easy to get wrong.** `slips[i, j]` has i along strike and
j down dip, with **j = 0 the shallow trench edge**. With `strike = 195` the
along-strike sense flips, so **i = 0 is the *south* end**.

The point of it is checkability. The KL source's amplitude knob (`slip`) has no
physical meaning, and the shipped calibration sits at `slip = 84` on a
200 x 50 km plane — Mw 8.95 packed onto a tenth of the real rupture area, an
effective source rather than a physical one. The deterministic source can be
checked against published quantities, and at M0 = 4.41e22 it hits all of them:
Mw 9.03, peak slip 52 m, peak sea-bed uplift 13.8 m, ~1.2 m of coastal
subsidence, slip above 20 m spanning ~240 km along strike.

Sweep it with `calibrate_deterministic.py`, which mirrors the notebook exactly
(DE_ader2, open DEM, Flather boundaries, 2 h evolve) and scores one run in
**64 s** on the GPU at the full ~355 000-triangle mesh:

```bash
python calibrate_deterministic.py --M0 3.8e22 --sig-u 0.15 --v0 0.25 --friction 0.045
```

**Run these one at a time** — one GPU, and concurrent runs contend for it.

#### &kappa; is scale-invariant, which explains the &kappa; floor

*K* is a geometric mean of `obs/model` and &kappa; its geometric standard
deviation. Scaling every modelled height by a constant `c` shifts each
log-ratio by `-log(c)` — a constant — leaving the standard deviation **exactly**
unchanged. So

> **no uniform amplitude knob can move &kappa;.** Not `slip`, not `M0`, not
> anything that rescales the whole modelled field. Such knobs move *K* alone.

Measured, at fixed shape and friction: M0 = 4.20e22 gives &kappa; = 1.718 and
M0 = 4.41e22 gives &kappa; = 1.715, a 5% amplitude change moving &kappa; by 0.003.

This is the mechanism behind the &kappa; ~ 1.7 floor recorded throughout this
file. Only levers that change the *shape* of the modelled height distribution —
friction, DEM, mesh, solver, source *geometry* — can touch &kappa; at all, and
measurably they move it by ~0.1 at best. Any future attempt on &kappa; must
change shape; tuning strength is provably wasted effort.

#### Calibration strategy

Three stages, but **only friction is truly orthogonal to DART**:

1. **`M0` from the DART 21418 peak.** Linear and friction-independent — one
   rescale by 1.87/1.78 moved the peak from 1.78 m to 1.87 m at 32 min and it
   stayed there across every friction tested.
2. **`v0` / `sig_v` from the coastal subsidence** (GEONET measured ~1.0-1.2 m at
   Oshika in 2011; the model's mean over the survey box is the proxy). This is
   an observable the repo had not been using, and it is what pins the down-dip
   slip position.
3. **`friction` from Aida *K*.**

**Stages 1 and 2 are not independent of each other.** Any geometry change —
`sig_u`, `width`, `length` — moves the DART peak, because DART responds to peak
uplift as well as displaced volume. Displaced volume *is* nearly invariant at
fixed moment (174-245 km3 across every distribution tried), but that is not
enough to pin DART. Re-run stage 1 after any geometry change.

#### Measured results

Friction curve at the best deterministic configuration found
(M0 = 3.80e22, `sig_u` 0.15, `v0` 0.25, 400 x 150 km at 23 km, open DEM):

| n | DART | K | &kappa; | bias | RMS | dry |
|-------|------|------|------|-------|------|-----|
| 0.033 | 1.85 | 0.72 | **1.66** | +1.89 | 4.29 | 62 |
| 0.040 | 1.85 | 0.88 | 1.73 | +0.92 | 3.95 | 64 |
| **0.045** | 1.85 | **1.03** | 1.81 | +0.26 | 3.83 | 65 |
| 0.050 | 1.85 | 1.22 | 1.92 | &minus;0.36 | 3.70 | 73 |

The source-shape survey, all at n = 0.033 unless noted:

| configuration | DART | K | &kappa; | bias | dry |
|---|------|------|------|-------|-----|
| `v0` 0.50, `sig_u` 0.20, M0 4.41e22 | 1.87 | 0.62 | 1.72 | +2.87 | 57 |
| `v0` 0.35 | 1.77 | 0.67 | 1.67 | +2.31 | 62 |
| `v0` 0.25 | 1.86 | 0.70 | **1.67** | +2.04 | 60 |
| `sig_u` 0.15 (`v0` 0.50) | 2.17 | 0.56 | 1.74 | +3.81 | **20** |
| `sig_u` 0.12, M0 3.50e22 | 1.96 | 0.72 | 1.66 | +1.89 | 62 |
| width 100 km, depth 17 km | 2.40 | 0.75 | 1.92 | +2.13 | 60 |
| width 75 km, depth 14 km | 2.52 | 0.73 | 1.79 | +2.22 | 64 |

**Verdict: the deterministic source is comparable to the shipped KL config, not
better.** Calibrated at n = 0.045 it gives K = 1.03, &kappa; = 1.81, bias +0.26,
RMS 3.83, 65 dry, DART 1.85 — against the KL's K = 1.00, &kappa; = 1.70,
bias +0.20, RMS 4.06, 68 dry. It wins slightly on RMS and dry count and loses
clearly on &kappa;. Use it for reproducibility and physical defensibility, not
for a better fit.

Four things to carry away:

- **For this source family *K* and &kappa; do not co-optimise.** &kappa; bottoms at
  **1.66** at n = 0.033 — the lowest value anywhere in this file — but *K* is
  0.72 there, and reaching *K* = 1 costs 0.15 of &kappa;. Under the KL source the
  two coincide at n = 0.033. The deterministic source carries a
  friction-independent excess of coastal loading relative to its DART
  amplitude, so friction gets pushed past its own optimum to compensate.
- **Width is not the fix.** Narrowing the plane raised the DART peak steeply
  (1.85 -> 2.40 -> 2.52 m) while *K* barely moved (0.72 -> 0.75 -> 0.73), so it
  does decouple far field from coast — but it degrades &kappa; (1.66 -> 1.92 at
  100 km), and since renormalising M0 to recover DART cannot move &kappa; at all
  (see above), that degradation is permanent. Tested and eliminated.
- **Compactness buys inundation *extent*.** `sig_u` = 0.15 with `v0` = 0.50 cut
  the dry count from 57 to **20**, the best figure on the open DEM by a wide
  margin (previous best 64; beating it before required the proprietary 150 m
  `Tohoku.pts`). It costs *K* and &kappa;, but if extent is what matters, a
  shorter sharper source is the lever.
- **&kappa; ~ 1.7 survives another class of test.** Source *distribution* — taper
  width, down-dip position, plane width, moment — is a lever never previously
  tried here, and it lands in 1.66-1.92 like everything else. Sources, DEMs,
  friction, mesh, solver, riverwalls and now slip distribution all bottom out
  in the same place.


### Split fault

`split_fault_deformation()` hangs several segments off a common trench anchor,
each abutting the previous in both horizontal distance and depth, with its own
dip, width, moment and taper. `segment_placement()` does the geometry; it is
covered by `test_split_fault_reduces_to_single_plane`, which pins a
one-segment split against `deterministic_deformation()` — worth having, since a
wrong placement shifts the source tens of km while still producing an entirely
plausible-looking deformation field.

The motivation is a real limitation of the single plane. Matching the ~1.2 m of
surveyed coastal subsidence forces slip **down dip**, which is exactly where it
over-loads the coast; trench-peaked slip fixes the coast and gives near-zero
subsidence. A listric source can do both, and does: a narrow split (shallow
75 km at 9 deg from 5 km depth, deep 90 km at 20 deg, half the moment each) is
the **only source tested here that satisfies all four physical observables at
once** — Mw 9.00, peak slip 57 m, peak uplift 13.9 m, subsidence &minus;1.36 m.
Every single-plane variant misses at least one; the best of them is 55% high on
uplift and produces a third of the observed subsidence.

Measured, all at n = 0.033 on the open DEM:

| geometry | M0 | DART | K | &kappa; | bias | RMS | dry |
|---|------|------|------|------|-------|------|-----|
| narrow, f = 0.5 | 4.0e22 | 0.86 | 0.97 | **1.68** | &minus;0.16 | **2.74** | 69 |
| narrow, f = 0.6 | 4.0e22 | 1.00 | 1.12 | 1.69 | &minus;0.91 | 2.71 | 76 |
| wide, f = 0.85 | 5.2e22 | 1.59 | 0.71 | 1.79 | +1.77 | 3.51 | 63 |
| wide, f = 0.85 | 6.2e22 | 1.90 | 0.60 | 1.78 | +2.95 | 4.54 | 57 |

`f` is the fraction of the moment on the shallow segment; "wide" is shallow
150 km at 9 deg with a 60 km deep segment.

**The split fault does decouple the far field from the coast** — the first
thing tested here that does. Moving moment from the deep segment to the shallow
one (f 0.5 -> 0.6) raised the DART peak 16% *and* dropped coastal heights 15%,
in opposite directions. Width, taper position and compactness all failed to do
that on a single plane.

**It still cannot match both.** Made strong enough to reach DART 1.87 it
over-predicts the coast worse than the single plane (K 0.60 against 0.72),
because keeping peak slip physical requires widening the shallow segment, which
broadens the uplift and destroys the decoupling. Extrapolating the narrow
geometry to hit both targets lands at f = 0.77, M0 = 6.05e22 and ~130 m of peak
slip. Confirmed from both directions: **within physical slip the split cannot
satisfy DART and the survey simultaneously.**

#### The unresolved tension

The f = 0.5 row is the important one. It is the **best coastal fit measured in
this repo** — K 0.97, &kappa; 1.68, RMS 2.74, bias &minus;0.16 m, with *no
friction retuning at all* — and it carries **46% of the required DART
amplitude**. For comparison, at K ~ 1 the best single plane manages &kappa; 1.81
/ RMS 3.83 and the recalibrated KL &kappa; 1.76 / RMS 4.36.

Lining the DART-matched sources up by how badly they then over-predict the
coast:

| source at DART ~ 1.87 | K at n = 0.033 |
|---|---|
| KL 200 x 50, slip 71 | 0.86 |
| deterministic 400 x 150 | 0.72 |
| split wide, f = 0.85 | 0.60 |
| *split narrow f = 0.5, at DART 0.86* | *0.97* |

So the coast prefers a source **roughly half** the far-field strength that DART
implies, and the discrepancy is a property of the model rather than of any
source: no distribution tested closes it. Candidates, none tested:

- the coast is over-amplified — 450 m bathymetry over the shelf, the effective
  friction standing in for unresolved roughness, or the shallow-water physics
  at the shelf break;
- the far field is under-propagated — numerical dissipation over the 500 km
  path to 21418, again on 450 m bathymetry. **Not the solver**: `DE1` on the
  identical split f = 0.5 case raised DART 10% (0.86 -> 0.95 m) but raised
  coastal heights ~20% (K 0.97 -> 0.81), taking the DART/coast ratio from
  0.834 to 0.769 — the gap widens — at 1.66x the runtime (106 s against 64 s).
  &kappa; moved 1.68 -> 1.69. Tested and eliminated;
- the source is mislocated in a way that costs far-field amplitude. Note the
  `x0` shift of 40 km west, which was tuned for *arrival time* on the old
  200 x 50 km plane and has never been revisited against amplitude.

That last one is the cheapest to check and should be first.


## Gotchas

- **`iseed` in `run_Tohoku_okada.py` does not select the slip realisation.** Despite its comment, the top-level `iseed` only feeds the run name (`Okada_<iseed>` → `_output_Okada_<iseed>`); the actual KL draw is fixed by the hard-coded `iseed=1001` in the `okl.kl_deformation(...)` call. `setup_simulation.apply_deformation()` similarly hard-codes `iseed=1234`. Change both if you want a genuinely different realisation.
- **Mesh caching differs between entry points.** `run_Tohoku_okada.py` calls `create_domain_from_regions(..., use_cache=True)` while `setup_simulation.create_domain()` uses `use_cache=False`; both write to the same `project.meshname`. After changing `rfact` or the polygons, delete `Tohoku_<scenario>_.msh` to be sure the mesh is rebuilt.
- **Different scripts run different durations.** `run_Tohoku_okada.py` evolves 4 hours at a 5-minute yieldstep with `tide = -0.45`; `setup_simulation.evolve_domain()` evolves 2 hours at 2 minutes with `tide = 0.0`. Don't compare their outputs without accounting for this.
- **`evolve_domain()` expects specific gauge keys** — it iterates `[21418, 0, 1, 2]`, so the `gauges` dict passed in must contain all four.
- **ANUGA's default Manning friction is 0.0, and most scripts here never set it.** That was the single largest error in the model: with nothing dissipating the wave between the shelf break and the shore, *every* source tried — the KL plane, a two-segment KL source, and all five published inversions in `sources/` — over-predicted the surveyed inundation by roughly a factor of two once it was strong enough to match DART 21418. `notebook_tohoku_open_elevation.ipynb` now sets `domain.set_quantity('friction', 0.04)`, which takes Aida K from 0.55 to 1.08 and leaves the DART peak unchanged to 0.01 m (friction is negligible in 5700 m of water). Treat 0.04 as an effective value standing in for unresolved bed roughness and coastal defences on a 450 m DEM, not a measured one — and note the calibration is coupled to the source: at n = 0.04, `slip = 60` gives K = 1.08 while `slip = 35` under-predicts at K = 1.58.
- **`nu` was silently ignored by the KL path, now fixed.** `kl_deformation()`
  accepted a `nu` argument but the inner `okada.forward()` call hard-coded
  `nu=0.25`, so changing it in a notebook did nothing. It is plumbed through
  now. The value 0.25 (Poisson solid) remains the right default, so no
  previously recorded result is affected.
- **The KL slip field's `alpha` was 0.75 and produced negative slip.** `alpha`
  in `kl_correlation_matrices` is the coefficient of variation (`sigma = alpha*mu`),
  and the expansion is Gaussian with no positivity constraint. At 0.75 a
  400 x 100 km / 10x10 draw gave slip from &minus;26.9 m to 100.5 m with **21 of
  100 subfaults slipping backwards**, and dragged the realised mean from the
  nominal 40 m down to 23.8 m. It is now **0.4**, the value at which negative
  slip disappears (CoV 0.47, min 4.3 m, peak 72.3 m). Note this raises the
  realised mean slip for a given nominal `slip`, so **the pre-existing KL
  calibration at `slip = 84` no longer holds** — recalibrate stage 1 before
  trusting a KL run. A second knob, `r0 = 0.2*width`, is the correlation
  length; at width 100 km it is 20 km against 40 km along-strike subfaults, so
  the field is nearly uncorrelated along strike. Longer `r0` smooths it but
  *increases* the spread at fixed `alpha`, so tune the two together.
  Recalibrated: the KL path now wants `slip = 71` with `friction = 0.038`
  (was 84 / 0.033) — see *Notebook defaults*. `notebook_tohoku_open_elevation.ipynb`
  carries the new pair. The other KL notebooks were **not** touched, because
  the calibration is specific to the 200 x 50 km plane on the open DEM under
  DE_ader2 and does not transfer to their geometry — but their sources are
  still 1.32x stronger than before the fix, so any tuning done in them is
  stale.
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

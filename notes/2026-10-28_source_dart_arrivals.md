# Source, DART 21418 signal and coastal arrival times

Discussion notes, 2026-10-28.  All model numbers are from
`notebook_tohoku_open_elevation.ipynb` as committed in `52b72ef` (open DEM,
DE_ader2 beta 1.0, n = 0.042).

## 1. The fault in the open-elevation notebook

A **two-segment split (listric) fault**.  Its geometry comes from the literature
rather than being fitted, the moment is the Global CMT value, and Manning
friction is the only parameter fitted to data.

**Placement.** Anchored by its **up-dip (trench) edge** at `(750290, 4181170)`
UTM 54N, not by its centroid.  Strike 195 deg, rake 87 deg (nearly pure thrust),
400 km along strike.  The anchor came from a sweep scored on the DART 21418
phase.  It leaves the model ~1.6 min early on purpose, because matching the
phase exactly needed 107 m of peak slip and a moment 9% above gCMT.

**Segments.** Each starts where the previous one ends, in horizontal position
and in depth:

| | near-trench | deep |
|---|---|---|
| dip | 5 deg | 15 deg |
| down-dip width | 80 km | 120 km |
| depth | 2-9 km | 9-40 km |
| share of moment | 90% | 10% |
| along-strike taper `sig_u` | 0.30 (sigma ~120 km) | 0.20 |
| down-dip taper | `v0` 0.45, `sig_v` 0.35 | same |

The steepening dip follows the Japan Trench interface (3-5 deg at the trench,
~12-15 deg inboard).  The 90% shallow share follows tsunami inversions, which
put the large slip above ~10 km depth.

**Slip.** A smooth, reproducible Gaussian taper on each segment (not a KL
random draw), with amplitude set by M0 = rigidity x area x mean slip.

**Strength and timing.** M0 = 5.3e22 N m (gCMT, Mw 9.08).  Released 69 s after
the origin time (the gCMT centroid time); the delay is applied when comparing
with DART.

**Applied to the water.** Okada deformation on 10 x 10 subfaults per segment.
The water surface gets the uplift plus the Tanioka & Satake horizontal-push
term, which is comparable to the uplift for a source this shallow.  Dropping it
costs 64% of the DART peak.

| quantity | model | observed / published |
|---|---|---|
| peak slip | 73 m | 50-70 m near the trench |
| peak sea-bed uplift | 14.6 m | |
| moment centroid | 37.87 N 143.32 E, 7.1 km | gCMT 37.52 N 143.05 E, 20 km |
| DART 21418 peak | 1.88 m | 1.87 m |
| coast (open DEM, n = 0.042) | K 1.00, kappa 1.58 | |

The model centroid sits ~44 km NNE of gCMT's and much shallower, because the
fault is placed by its trench edge.

## 2. The early signal in the DART 21418 record

The spikes in the first ~10 min are **seismic waves shaking the sea floor**,
not water.

A DART station measures bottom pressure (~5.7 km down) and reports
eta = dp / (rho g).  A vertical sea-floor acceleration *a* changes the pressure
by rho h a, so it reads as

    eta ~ h a / g

At h = 5700 m, an acceleration of only ~1.7 mm/s^2 (0.0002 g) reads as 1 m.

Epicentre to buoy is 552 km.  From the origin time:

| arrival | expected | record |
|---|---|---|
| P (~8 km/s) | ~1.1 min | up to +3 cm in the samples stamped 96 s |
| S (~4.5 km/s) | ~2.0 min | first large values at 2.6 min (+0.83 m) |
| Rayleigh (~3.7 km/s) | ~2.5 min | largest swings at 3.6 min (-0.96 m) |
| tsunami | | rise from ~27 min, peak 1.87 m at ~33 min |

Smaller oscillations follow to ~10 min, likely sound waves bouncing between the
sea surface and the sea floor plus seismic coda.  The record is then quiet
(0.018 m standard deviation over 15-25 min).  The buoy is far outside the
permanent deformation zone, so only the pulse at ~30 min is real sea surface.

The record looks spiky because these waves have periods of seconds to tens of
seconds, far below its sampling interval.  At 1.6, 2.6 and 3.6 min there are
**four values per timestamp**, which looks like 15 s event-mode samples stamped
to the whole minute.  That explains the repeated timestamps in
`21418_notide.txt`, but it is inferred from the file, not documented.

The model has none of this signal: ANUGA solves the shallow-water equations
with an instantaneous source, so it has no elastic Earth, no seismic waves and
no compressible water.  This is why the DART phase fit in
`calibrate_deterministic.py` uses only the 20-60 min window.

## 3. Why the P wave barely shows

It is probably there.  The samples stamped 96 s, just after the expected
arrival at ~66 s, wobble by up to +3 cm against ~0.6 cm before.  Two reasons it
stays small:

1. **P waves are the weakest arrivals.** At 550 km, P moves the ground much
   less than S, and S much less than surface waves.  Surface waves spread over a
   surface rather than a volume, so they lose amplitude more slowly with distance
   and dominate at regional range.
2. **The recorder averages.** DART reports averaged pressure (15 s in event
   mode).  P waves are ~1 Hz, so dozens of cycles fall in each average and
   cancel.  S and surface waves (periods ~10-100 s) survive.  The averaging is a
   low-pass filter that suppresses the P wave most.

Caveats: the interval a timestamp covers is inferred, not documented.  The
pressure-to-height conversion is only an analogy for seismic motion: it shows
*that* the floor shook, not true seismic amplitudes.

## 4. Arrival times along the coast

From `tohoku_open_dem.sww`.  "First rise" is the earliest 10% of nearshore
(0-10 m deep) cells per band rising +0.5 m above the initial stage.  Minutes are
after the origin time (14:46 JST), including the 69 s source delay.  Frames are
every 2 min, so times are good to about +-1-2 min.

| coast | first +0.5 m | first +2 m | peak (median) |
|---|---|---|---|
| 39.6-39.8 N | 27 | 31 | 33 |
| Sanriku rias, 39.0-39.6 N | 17-19 | 27-29 | 35-37 |
| 38.6-39.0 N | 23-25 | 29 | 43-45 |
| Oshika / Onagawa, 38.4-38.6 N | 25 | 29 | 61 |
| 38.2-38.4 N (incl. Oshika tip) | 35 | 35 | 69 |
| Sendai plain, 38.0-38.2 N | 59 | 63 | 73 |
| 37.8-38.0 N | 49 | 55 | 69 |
| Fukushima, 37.4-37.8 N | 37-41 | 43-47 | 55-59 |

| place | ground | first +0.5 m | max stage |
|---|---|---|---|
| Sendai Bay offshore | -37.2 m | 39 min | 7.3 m at 53 min |
| Arahama (nearshore) | -8.8 m | 59 min | 7.6 m at 71 min |
| Yuriage | +0.5 m | 69 min | 6.3 m at 75 min |
| Sendai airport | +1.8 m | 75 min | 4.7 m at 77 min |

**Why the Sendai plain is latest despite being nearest the epicentre:**

- **A wide, shallow shelf.** Speed is sqrt(g h): ~850 km/h at 5.7 km depth,
  ~300 km/h at 700 m, ~80 km/h at 50 m.  Crossing ~50 km of shallow Sendai Bay
  adds most of half an hour.
- **The north is near the source in deep water.** The Sanriku coast faces the
  northern half of the rupture across a narrow shelf.
- **The source starts far offshore.** The uplift is near the trench, about
  200 km out, so even the fastest arrivals take ~20 min.

**Against what happened** (from memory, not checked against primary sources):
first large waves reached the Iwate ria coast ~25-30 min after the earthquake,
and the Sendai plain ~15:50-16:00 JST (~65-75 min), with Sendai airport flooded
at ~70 min.  The model is in that range, consistent with its ~1.6 min early
arrival at DART 21418.  It was never fitted to coastal arrival times, so this
is a check, not a calibration.

### To reproduce section 4

```python
import numpy as np, anuga
from pyproj import Transformer
sp = anuga.SWW_plotter('tohoku_open_dem.sww', absolute=True)
t, st, el = np.asarray(sp.time), np.asarray(sp.stage), np.asarray(sp.elev)
lon, lat = Transformer.from_crs('EPSG:32654', 'EPSG:4326', always_xy=True).transform(sp.xc, sp.yc)
e0 = el[0] if el.ndim == 2 else el
near = (e0 < -0.5) & (e0 > -10) & (lat > 37.4) & (lat < 39.8) & (lon < 142.2)
rise = st[:, near] - st[0, near]
hit = rise > 0.5
first = np.where(hit.any(0), t[hit.argmax(0)], np.nan)   # s after release; add 69 s for origin time
```

"""Scriptable calibration driver for the deterministic Okada source.

Mirrors the configuration of notebook_tohoku_open_elevation.ipynb (DE_ader2,
open DEM, Flather boundaries, 2 hour evolve) but takes the source and friction
from the command line and prints the scores, so the three calibration stages
can be swept without re-running a notebook by hand.

    python calibrate_deterministic.py --M0 4.2e22 --friction 0.033
"""
import argparse, json, os, time
import numpy as np
import anuga
import okada_kl_subfaults as okl
import project

km = 1000.0
closeup_extent = [475_000, 580_000, 4_150_000, 4_280_000]
md = 0.01

p = argparse.ArgumentParser()
p.add_argument('--flow-algorithm', default='DE_ader2')
p.add_argument('--elevation-source', choices=['open', 'pts'], default='open')
p.add_argument('--no-horizontal-push', dest='push', action='store_false',
               help='drop the Tanioka & Satake (1996) horizontal-push term '
                    '(what every run before 2026-09 did)')
p.set_defaults(push=True)
p.add_argument('--source', choices=['det', 'kl', 'split'], default='det')
p.add_argument('--split-frac', type=float, default=0.6)   # moment fraction, shallow segment
p.add_argument('--dip1', type=float, default=8.0)
p.add_argument('--width1', type=float, default=75000.0)
p.add_argument('--top1', type=float, default=5000.0)
p.add_argument('--dip2', type=float, default=18.0)
p.add_argument('--width2', type=float, default=90000.0)
p.add_argument('--xt', type=float, default=750290.0)   # trench edge, UTM 54N
p.add_argument('--yt', type=float, default=4181170.0)
p.add_argument('--slip', type=float, default=84.0)   # kl only: nominal mean slip
p.add_argument('--iseed', type=int, default=1234)    # kl only
p.add_argument('--M0', type=float, default=4.2e22)
p.add_argument('--sig-u', type=float, default=0.20)
p.add_argument('--v0', type=float, default=0.5)
p.add_argument('--sig-v', type=float, default=0.30)
p.add_argument('--friction', type=float, default=0.033)
p.add_argument('--length', type=float, default=400*km)
p.add_argument('--width', type=float, default=150*km)
p.add_argument('--depth', type=float, default=23*km)
p.add_argument('--x0', type=float, default=680000.0)
p.add_argument('--y0', type=float, default=4200000.0)
p.add_argument('--finaltime', type=float, default=2.0*3600)
p.add_argument('--tag', default='det')
p.add_argument('--results', default='calibration_results.jsonl')
a = p.parse_args()

t_start = time.time()

domain = anuga.create_domain_from_regions(
    project.bounding_polygon,
    boundary_tags={'bottom': [0], 'ocean_east': [1], 'top': [2], 'onshore': [3]},
    maximum_triangle_area=project.res_whole,
    interior_regions=project.interior_regions,
    use_cache=False, verbose=False)

domain.set_flow_algorithm(a.flow_algorithm)
domain.set_epsg(32654)
domain.set_collect_max_quantities()
domain.set_compute_mode('unified')
domain.set_name(f'calib_{a.tag}')

xll = domain.geo_reference.xllcorner
yll = domain.geo_reference.yllcorner
print(f'Triangles: {len(domain)}')

domain.set_quantity('elevation', filename='Tohoku_dem.tif', location='centroids')

if a.elevation_source == 'pts':
    # Tohoku.pts (1.34 M points, ~150 m) overlaid where it has coverage, with
    # the GeoTIFF filling the rest -- it cannot supply the deep ocean east of
    # x ~ 899 km on its own.  The LinearNDInterpolator build is a Delaunay
    # triangulation of 1.34 M points, so cache the result on the mesh: a
    # friction sweep reuses the same mesh every run.
    cache = f'_pts_elev_{len(domain)}.npy'
    Elev0 = domain.quantities['elevation'].centroid_values
    if os.path.exists(cache):
        elev_pts = np.load(cache)
        print(f'pts elevation: loaded {cache}')
    else:
        from scipy.interpolate import LinearNDInterpolator
        from anuga.geospatial_data.geospatial_data import Geospatial_data
        t_i = time.time()
        pts_dem = Geospatial_data('Tohoku.pts')
        interp = LinearNDInterpolator(pts_dem.get_data_points(absolute=True),
                                      np.asarray(pts_dem.get_attributes()))
        elev_pts = interp(domain.centroid_coordinates[:, 0] + xll,
                          domain.centroid_coordinates[:, 1] + yll)
        np.save(cache, elev_pts)
        print(f'pts elevation: built in {time.time()-t_i:.0f} s -> {cache}')
    covered = np.isfinite(elev_pts)
    Elev0[:] = np.where(covered, elev_pts, Elev0)
    print(f'Tohoku.pts covers {covered.sum()} of {covered.size} cells')
tide = 0.0
domain.set_quantity('stage', tide)
Elevation = domain.quantities['elevation'].centroid_values
Stage     = domain.quantities['stage'].centroid_values
Stage[:]  = np.maximum(Elevation, Stage)
domain.set_quantity('friction', a.friction)

# --- deterministic source -------------------------------------------------
x = domain.centroid_coordinates[:, 0]
y = domain.centroid_coordinates[:, 1]
if a.source == 'det':
    uE, uN, uZ, slips = okl.deterministic_deformation(
        x, y, xoff=a.x0 - xll, yoff=a.y0 - yll,
        M0=a.M0, u0=0.5, sig_u=a.sig_u, v0=a.v0, sig_v=a.sig_v,
        depth=a.depth, length=a.length, width=a.width)
elif a.source == 'split':
    f = a.split_frac
    segs = [dict(dip=a.dip1, width=a.width1, top_depth=a.top1, M0=f*a.M0,
                 sig_u=a.sig_u, v0=0.45, sig_v=0.35),
            dict(dip=a.dip2, width=a.width2, M0=(1.0 - f)*a.M0,
                 sig_u=0.20, v0=0.45, sig_v=0.35)]
    uE, uN, uZ, slips_list = okl.split_fault_deformation(
        x, y, a.xt - xll, a.yt - yll, segs, length=a.length)
    slips = np.concatenate([s_.ravel() for s_ in slips_list])
else:
    uE, uN, uZ, slips = okl.kl_deformation(
        x, y, xoff=a.x0 - xll, yoff=a.y0 - yll,
        E_subfault=10, N_subfault=10, sample='sobol', iseed=a.iseed,
        depth=a.depth, length=a.length, width=a.width,
        strike=195.0, dip=14.0, rake=87.0, nu=0.25,
        slip=a.slip, opening=0.0)

# Tanioka & Satake (1996): where the sea bed slopes, horizontal motion
# displaces water too.  Moving the bed east by uE replaces the bed at a point
# with bed that was previously to the west, so the effective vertical change is
# -uE*dz/dx - uN*dz/dy on top of uZ.  The term scales with horizontal
# displacement and bathymetric slope, so it is largest for shallow near-trench
# slip -- exactly where uZ (which scales as slip*sin(dip)) is weakest.  It is
# applied to the water surface only: the sea bed itself moves by uZ.
if a.push:
    elev_q = domain.quantities['elevation']
    elev_q.compute_gradients()                 # pre-earthquake bathymetry
    dzdx = elev_q.x_gradient
    dzdy = elev_q.y_gradient
    push = -(uE*dzdx + uN*dzdy)
    push = np.where(Elevation < 0.0, push, 0.0)   # only where there is water
    print(f'horizontal push: {push.min():.2f} .. {push.max():.2f} m '
          f'(uZ {uZ.min():.2f} .. {uZ.max():.2f})')
else:
    push = 0.0

Elevation[:] += uZ
Stage[:]     += uZ + push
if a.source == 'split':
    Mw = (np.log10(a.M0) - 9.1)/1.5
else:
    Mw = (np.log10(40e9*a.length*a.width*slips.mean()) - 9.1)/1.5
print(f'source: Mw {Mw:.2f}  mean slip {slips.mean():.1f} m  peak {slips.max():.1f} m  '
      f'uZ {uZ.min():.2f}..{uZ.max():.2f} m')

# --- evolve ---------------------------------------------------------------
Br = anuga.Reflective_boundary(domain)
Bf = anuga.Flather_external_stage_zero_velocity_boundary(domain, lambda t: tide)
domain.set_boundary({'onshore': Br, 'bottom': Bf, 'ocean_east': Bf, 'top': Bf})

for t in domain.evolve(yieldstep=2*60., finaltime=a.finaltime):
    pass
runtime = time.time() - t_start
print(f'evolve done in {runtime:.0f} s')

# --- score ----------------------------------------------------------------
import utm
from scipy.spatial import cKDTree
from tsunami_observations import load_ttjs, in_extent

sp = anuga.SWW_plotter(domain.get_name() + '.sww', absolute=True)
sww_time = np.asarray(sp.time)
max_stage = np.asarray(sp.max_stage)
max_depth = np.asarray(sp.max_depth)

# DART 21418
de, dn, _, _ = utm.from_latlon(38.711, 148.694, force_zone_number=54)
cell = int(np.argmin((sp.xc - de)**2 + (sp.yc - dn)**2))
series = np.asarray(sp.stage)[:, cell]
dart_peak = float(series.max())
dart_tmin = float(sww_time[series.argmax()]/60.)

# Aida K / kappa over the close-up survey box
survey = in_extent(load_ttjs(types='IR', reliability='AB'), closeup_extent)
wet = max_depth >= md
tree = cKDTree(np.column_stack([sp.xc[wet], sp.yc[wet]]))
wstage = max_stage[wet]
d, j = tree.query(np.column_stack([survey['east'], survey['north']]),
                  distance_upper_bound=1500.0)
hit = np.isfinite(d)
h_model = np.full(survey['east'].shape, np.nan)
h_model[hit] = wstage[j[hit]]
h_obs = survey['height'] + tide

paired = np.isfinite(h_model) & (h_model > 0) & (h_obs > 0)
n_dry = int((~np.isfinite(h_model)).sum())
ratio = h_obs[paired]/h_model[paired]
log_K = np.mean(np.log(ratio))
K = float(np.exp(log_K))
kappa = float(np.exp(np.sqrt(np.mean((np.log(ratio) - log_K)**2))))
bias = float(np.mean(h_model[paired] - h_obs[paired]))
rms = float(np.sqrt(np.mean((h_model[paired] - h_obs[paired])**2)))

res = dict(tag=a.tag, push=bool(a.push), algorithm=a.flow_algorithm, dem=a.elevation_source, source=a.source, slip=a.slip, split_frac=a.split_frac, M0=a.M0, Mw=round(Mw,3), sig_u=a.sig_u, v0=a.v0, sig_v=a.sig_v,
           friction=a.friction, length=a.length, width=a.width, depth=a.depth,
           mean_slip=round(float(slips.mean()),2), peak_slip=round(float(slips.max()),2),
           uZ_max=round(float(uZ.max()),2), uZ_min=round(float(uZ.min()),2),
           dart_peak=round(dart_peak,3), dart_t_min=round(dart_tmin,1),
           K=round(K,3), kappa=round(kappa,3), bias=round(bias,3), rms=round(rms,3),
           n_paired=int(paired.sum()), n_dry=n_dry, triangles=len(domain),
           runtime_s=round(runtime,1))

print('\n' + '-'*72)
print(f"DART peak {dart_peak:.2f} m at {dart_tmin:.0f} min   (observed 1.87 m at 33 min)")
print(f"K {K:.2f}   kappa {kappa:.2f}   bias {bias:+.2f} m   RMS {rms:.2f}   dry {n_dry}/{survey['east'].size}")
print('-'*72)

with open(a.results, 'a') as f:
    f.write(json.dumps(res) + '\n')

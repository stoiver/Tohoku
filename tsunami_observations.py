"""
Observed 2011 Tohoku tsunami inundation and run-up heights, for validating the
modelled maximum stage against the field survey.

Primary source: the 2011 Tohoku Earthquake Tsunami Joint Survey (TTJS) Group
  - 5907 surveyed points along ~2000 km of coast, 31.5 N to 43.7 N
  - Inundation height (type 'I'), run-up height (type 'R'), harbour mark ('P')
  - Heights are tide-corrected, relative to the astronomical tide at the time
    the wave arrived, so they compare directly with ANUGA `stage` for a run
    with tide = 0.0 (see `tide` in the run scripts before comparing)
  - Release 20121229, http://www.coastal.jp/ttjt/ -- no key required
  - Local copy: observations/ttjs_survey_20121229.csv (converted to UTF-8),
    with the upstream readme in observations/ttjs_readme.txt

Secondary source: NOAA NCEI/WDS Global Historical Tsunami Database, run-up
table for event 5413 (2011 Honshu).  Coarser (tide gauges, eyewitness reports
and survey summaries worldwide) but it carries arrival times, which TTJS does
not.  See `download_ncei_runups()`.

Cite TTJS as required by their readme, e.g.
    "Data are from the 2011 Tohoku Earthquake Tsunami Joint Survey Group,
     release 20121229, http://www.coastal.jp/ttjt/"

Usage
-----
    from tsunami_observations import load_ttjs, in_extent

    obs = load_ttjs()                       # all points, UTM 54N added
    obs = load_ttjs(types='I')              # inundation heights only
    obs = load_ttjs(reliability='AB')       # drop the low-confidence points

    sendai = in_extent(obs, [484_000, 531_500, 4_209_000, 4_252_500])
    plt.scatter(sendai['east'], sendai['north'], c=sendai['height'])
"""

import os

import numpy as np

TTJS_CSV = os.path.join('observations', 'ttjs_survey_20121229.csv')

TTJS_URL = ('https://www.coastal.jp/ttjt/index.php?plugin=attach'
            '&refer=%E7%8F%BE%E5%9C%B0%E8%AA%BF%E6%9F%BB%E7%B5%90%E6%9E%9C'
            '&openfile=ttjt_survey_29-Dec-2012_tidecorrected_web.csv')

NCEI_RUNUP_URL = ('https://www.ngdc.noaa.gov/hazel/hazard-service/api/v1/'
                  'tsunamis/runups?tsunamiEventId=5413')


def load_ttjs(filename=TTJS_CSV, types='IR', reliability=None, zone=54):
    """
    Load the TTJS survey points as a dict of arrays.

    types        which measurement types to keep, as a string of the codes
                 'I' (inundation height), 'R' (run-up height), 'P' (harbour
                 mark), 'W' (too weak to measure).  None keeps everything.
    reliability  keep only these reliability grades, e.g. 'AB' for the clear
                 marks and reliable witness reports.  None keeps everything.

    Returns keys: lon, lat, east, north (UTM `zone`), height (m, tide-corrected
    by TTJS), height_msl (m above mean sea level), runup_distance (m inland
    from the shoreline, NaN where not surveyed), type, reliability, location.

    Heights are the 'height corrected by ttjt' column -- the one the survey
    group recommends -- and are NaN for the handful of points they could not
    correct.  Those rows are dropped.
    """
    import csv

    import utm

    with open(filename, newline='', encoding='utf-8') as fh:
        rows = list(csv.DictReader(fh))

    def col(row, key):
        try:
            return float(row[key])
        except (TypeError, ValueError):
            return np.nan

    keep = []
    for row in rows:
        t = row['type'].strip()
        r = row['reliability'].strip()
        if types is not None and t not in types:
            continue
        if reliability is not None and r not in reliability:
            continue
        if not np.isfinite(col(row, 'height corrected by ttjt [m]')):
            continue
        # Three points carry no position (lon and/or lat are NaN); they cannot
        # be projected or plotted, so they go too.
        if not (np.isfinite(col(row, 'lon [deg]')) and
                np.isfinite(col(row, 'lat [deg]'))):
            continue
        keep.append(row)

    lon = np.array([col(r, 'lon [deg]') for r in keep])
    lat = np.array([col(r, 'lat [deg]') for r in keep])

    east = np.empty_like(lon)
    north = np.empty_like(lat)
    for i, (la, lo) in enumerate(zip(lat, lon)):
        east[i], north[i], _, _ = utm.from_latlon(la, lo, force_zone_number=zone)

    return dict(
        lon=lon,
        lat=lat,
        east=east,
        north=north,
        height=np.array([col(r, 'height corrected by ttjt [m]') for r in keep]),
        height_msl=np.array([col(r, 'height from msl [m]') for r in keep]),
        runup_distance=np.array([col(r, 'runup distance [m]') for r in keep]),
        type=np.array([r['type'].strip() for r in keep]),
        reliability=np.array([r['reliability'].strip() for r in keep]),
        location=np.array([r['location'].strip() for r in keep]),
    )


def in_extent(obs, extent):
    """Subset `obs` to a UTM extent [xmin, xmax, ymin, ymax], as used by the
    notebook map plots."""
    xmin, xmax, ymin, ymax = extent
    m = ((obs['east'] >= xmin) & (obs['east'] <= xmax) &
         (obs['north'] >= ymin) & (obs['north'] <= ymax))
    return {k: v[m] for k, v in obs.items()}


def near(obs, gauge, radius=3000.0):
    """Subset `obs` to points within `radius` metres of a Gauge (anything with
    .east/.north), nearest first."""
    d = np.hypot(obs['east'] - gauge.east, obs['north'] - gauge.north)
    order = np.argsort(d)
    order = order[d[order] <= radius]
    out = {k: v[order] for k, v in obs.items()}
    out['distance'] = d[order]
    return out


def download_ttjs(filename=TTJS_CSV):
    """Re-download the TTJS survey CSV and store it as UTF-8.  The upstream
    file is Shift-JIS with a count line above the header; both are handled
    here so the stored copy is a plain UTF-8 CSV."""
    import csv
    import io

    import requests

    response = requests.get(TTJS_URL, timeout=120)
    response.raise_for_status()

    text = response.content.decode('shift_jis', errors='replace')
    rows = list(csv.reader(io.StringIO(text)))

    directory = os.path.dirname(filename)
    if directory:
        os.makedirs(directory, exist_ok=True)

    with open(filename, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        for row in rows[1:]:                     # row 0 is just the count
            writer.writerow([c.strip() for c in row])

    print(f'Wrote {filename}: {len(rows) - 2} survey points')
    return filename


def download_ncei_runups(filename=os.path.join('observations',
                                               'ncei_runups_2011honshu.json')):
    """Download the NCEI run-up records for the 2011 Honshu event (id 5413).
    Paged: the service returns 6427 records in pages of 100 by default."""
    import json

    import requests

    items = []
    page = 1
    while True:
        response = requests.get(f'{NCEI_RUNUP_URL}&page={page}', timeout=120)
        response.raise_for_status()
        payload = response.json()
        items.extend(payload['items'])
        if page >= payload['totalPages']:
            break
        page += 1

    directory = os.path.dirname(filename)
    if directory:
        os.makedirs(directory, exist_ok=True)
    with open(filename, 'w', encoding='utf-8') as fh:
        json.dump(items, fh)

    print(f'Wrote {filename}: {len(items)} runup records')
    return filename

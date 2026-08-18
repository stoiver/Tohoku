# Observed tsunami heights, 2011 Tohoku

Field data for validating the modelled maximum stage.  Loaded by
`tsunami_observations.py` in the repository root.

## `ttjs_survey_20121229.csv`

The 2011 Tohoku Earthquake Tsunami Joint Survey (TTJS) Group nationwide survey:
5907 points from 31.5 N to 43.7 N, measured by 299 researchers from 64
institutes.  Downloaded from <http://www.coastal.jp/ttjt/>, release 20121229,
and converted from Shift-JIS to UTF-8 with the upstream count line above the
header removed; columns are otherwise untouched, including the Japanese
location/target/group text.  `ttjs_readme.txt` is the upstream readme.

Columns that matter:

| Column | Meaning |
|--------|---------|
| `lon [deg]`, `lat [deg]` | position (WGS84) |
| `measured height [m]` | raw measurement |
| `height corrected by ttjt [m]` | **use this** — tide-corrected by the survey group, relative to the astronomical tide when the wave arrived |
| `height from msl [m]` | same, referred to mean sea level |
| `runup distance [m]` | distance inland from the shoreline |
| `type` | `I` inundation height, `R` run-up height, `P` harbour mark, `W` too weak to measure |
| `reliability` | `A` clear mark → `D` marginal |

Three points have no coordinates and a further ~120 have no tide-corrected
height; `load_ttjs()` drops both.

Citation, as required by the survey group's terms:

> Data are from the 2011 Tohoku Earthquake Tsunami Joint Survey Group,
> release 20121229, <http://www.coastal.jp/ttjt/>

Reference: Mori, N., T. Takahashi, T. Yasuda and H. Yanagisawa (2011), Survey
of 2011 Tohoku earthquake tsunami inundation and run-up, *Geophysical Research
Letters*, 38, L00G14, doi:10.1029/2011GL049210.

## NCEI run-ups (not stored)

`tsunami_observations.download_ncei_runups()` fetches the 6427 run-up records
NOAA NCEI holds for this event (tsunami event id 5413) as JSON.  Coarser than
TTJS along the Tohoku coast, but global in extent and it carries arrival times,
which TTJS does not.

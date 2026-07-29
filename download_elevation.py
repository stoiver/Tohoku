"""
Download openly available global bathymetry/topography data for use as ANUGA
elevation input.

Primary source: NOAA NCEI Global DEM Mosaic (GEBCO-based)
  - No API key required
  - Global coverage at 15 arc-second (~450 m) resolution
  - Land topography (SRTM-derived) + ocean bathymetry (GEBCO)
  - URL: https://gis.ngdc.noaa.gov/arcgis/rest/services/DEM_mosaics/DEM_global_mosaic/ImageServer

Alternative source: OpenTopography (free API key required, see opentopography.org)
  - Datasets: SRTMGL1 (1 arcsec, land only), SRTMGL3 (3 arcsec, land only),
    SRTM15Plus (15 arcsec, land + ocean), GEBCO (15 arcsec, land + ocean)
  - URL: https://portal.opentopography.org/API/globaldem

Usage
-----
    # Automatic bounds from project.py:
    download_noaa_dem('Tohoku_dem.tif')

    # Explicit bounds:
    download_noaa_dem('Tohoku_dem.tif', bounds=(140.0, 151.0, 35.0, 41.5))

    # In setup_simulation.py / notebooks — use the tif directly:
    domain.set_quantity('elevation', filename='Tohoku_dem.tif',
                        location='centroids', verbose=True)
"""

import os
import numpy as np
import requests
import rasterio
from pyproj import Transformer


def get_domain_bounds_latlon(padding_deg=0.5):
    """
    Return (lon_min, lon_max, lat_min, lat_max) for the domain defined in
    project.py, padded by `padding_deg` on each side.

    Requires project.py (bounding_polygon in UTM Zone 54N / EPSG:32654).
    """
    import project

    transformer = Transformer.from_crs('EPSG:32654', 'EPSG:4326', always_xy=True)
    poly = np.array(project.bounding_polygon)
    lons, lats = transformer.transform(poly[:, 0], poly[:, 1])

    return (
        float(np.min(lons)) - padding_deg,
        float(np.max(lons)) + padding_deg,
        float(np.min(lats)) - padding_deg,
        float(np.max(lats)) + padding_deg,
    )


def download_noaa_dem(output_filename, bounds=None, resolution_arcsec=15,
                      padding_deg=0.5, verbose=True):
    """
    Download global bathymetry/topography from the NOAA NCEI Global DEM Mosaic
    (GEBCO-based) as a GeoTIFF.

    No API key is required. The returned file is in EPSG:4326 and can be passed
    directly to domain.set_quantity('elevation', filename=...).

    Parameters
    ----------
    output_filename : str
        Path for the output GeoTIFF.
    bounds : tuple or None
        (lon_min, lon_max, lat_min, lat_max) in degrees WGS84. If None, derived
        automatically from project.py (requires project.py to be importable).
    resolution_arcsec : float
        Output resolution in arc-seconds. Default 15 matches GEBCO native
        resolution (~450 m at mid-latitudes). Use 30 or 60 for faster downloads.
    padding_deg : float
        Padding in degrees to add around `bounds` (or the project domain).
        Prevents edge effects during ANUGA interpolation.
    verbose : bool

    Returns
    -------
    str
        Path of the written GeoTIFF (same as output_filename).
    """
    if bounds is None:
        lon_min, lon_max, lat_min, lat_max = get_domain_bounds_latlon(padding_deg)
    else:
        lon_min, lon_max, lat_min, lat_max = bounds
        lon_min -= padding_deg
        lon_max += padding_deg
        lat_min -= padding_deg
        lat_max += padding_deg

    pixels_per_degree = 3600.0 / resolution_arcsec
    width  = int(round((lon_max - lon_min) * pixels_per_degree))
    height = int(round((lat_max - lat_min) * pixels_per_degree))

    if verbose:
        print(f'Downloading NOAA NCEI Global DEM (GEBCO-based)')
        print(f'  Lon: {lon_min:.3f} to {lon_max:.3f}')
        print(f'  Lat: {lat_min:.3f} to {lat_max:.3f}')
        print(f'  Resolution: {resolution_arcsec} arcsec  ({width} x {height} pixels)')

    url = ('https://gis.ngdc.noaa.gov/arcgis/rest/services'
           '/DEM_mosaics/DEM_global_mosaic/ImageServer/exportImage')
    params = {
        'bbox':          f'{lon_min},{lat_min},{lon_max},{lat_max}',
        'bboxSR':        4326,
        'size':          f'{width},{height}',
        'imageSR':       4326,
        'format':        'tiff',
        'pixelType':     'F32',
        'interpolation': 'RSP_NearestNeighbor',
        'f':             'image',
    }

    response = requests.get(url, params=params, stream=True, timeout=300)
    response.raise_for_status()

    content_type = response.headers.get('Content-Type', '')
    if 'json' in content_type or 'xml' in content_type:
        raise RuntimeError(f'Server returned an error instead of a TIFF:\n{response.text[:500]}')

    os.makedirs(os.path.dirname(os.path.abspath(output_filename)), exist_ok=True)
    with open(output_filename, 'wb') as fid:
        for chunk in response.iter_content(chunk_size=65536):
            fid.write(chunk)

    if verbose:
        _print_tif_summary(output_filename)

    return output_filename


def download_opentopography_dem(output_filename, api_key, bounds=None,
                                demtype='SRTM15Plus', padding_deg=0.5,
                                verbose=True):
    """
    Download elevation/bathymetry from OpenTopography using a free API key.

    API keys are free at https://opentopography.org/

    Parameters
    ----------
    output_filename : str
        Path for the output GeoTIFF.
    api_key : str
        Your OpenTopography API key.
    bounds : tuple or None
        (lon_min, lon_max, lat_min, lat_max). If None, derived from project.py.
    demtype : str
        Dataset identifier. Options with ocean bathymetry:
          'SRTM15Plus' — 15 arc-second, global land + ocean (Tozer et al. 2019)
          'GEBCO'      — 15 arc-second, GEBCO 2014
        Land-only (no bathymetry):
          'SRTMGL1'    — 1 arc-second (~30 m), SRTM v3
          'SRTMGL3'    — 3 arc-second (~90 m), SRTM v3
    padding_deg : float
    verbose : bool

    Returns
    -------
    str
        Path of the written GeoTIFF.
    """
    if bounds is None:
        lon_min, lon_max, lat_min, lat_max = get_domain_bounds_latlon(padding_deg)
    else:
        lon_min, lon_max, lat_min, lat_max = bounds
        lon_min -= padding_deg
        lon_max += padding_deg
        lat_min -= padding_deg
        lat_max += padding_deg

    if verbose:
        print(f'Downloading {demtype} from OpenTopography')
        print(f'  Lon: {lon_min:.3f} to {lon_max:.3f}')
        print(f'  Lat: {lat_min:.3f} to {lat_max:.3f}')

    params = {
        'demtype':      demtype,
        'south':        lat_min,
        'north':        lat_max,
        'west':         lon_min,
        'east':         lon_max,
        'outputFormat': 'GTiff',
        'API_Key':      api_key,
    }

    response = requests.get(
        'https://portal.opentopography.org/API/globaldem',
        params=params, stream=True, timeout=300,
    )
    response.raise_for_status()

    content_type = response.headers.get('Content-Type', '')
    if 'xml' in content_type or 'json' in content_type:
        raise RuntimeError(f'OpenTopography error:\n{response.text[:500]}')

    os.makedirs(os.path.dirname(os.path.abspath(output_filename)), exist_ok=True)
    with open(output_filename, 'wb') as fid:
        for chunk in response.iter_content(chunk_size=65536):
            fid.write(chunk)

    if verbose:
        _print_tif_summary(output_filename)

    return output_filename


def _print_tif_summary(filename):
    with rasterio.open(filename) as r:
        z = r.read(1).astype(float)
        if r.nodata is not None:
            z[z == r.nodata] = np.nan
        size_mb = os.path.getsize(filename) / 1e6
        print(f'  Saved {size_mb:.1f} MB → {filename}')
        print(f'  CRS: {r.crs}   Size: {r.width} x {r.height} px')
        print(f'  Bounds: lon [{r.bounds.left:.3f}, {r.bounds.right:.3f}]'
              f'  lat [{r.bounds.bottom:.3f}, {r.bounds.top:.3f}]')
        print(f'  Elevation range: {np.nanmin(z):.1f} m to {np.nanmax(z):.1f} m')


if __name__ == '__main__':
    # Download GEBCO-quality DEM covering the Tohoku simulation domain.
    # Output can be used directly: domain.set_quantity('elevation', filename='Tohoku_dem.tif')
    download_noaa_dem('Tohoku_dem.tif', resolution_arcsec=15, padding_deg=0.5)

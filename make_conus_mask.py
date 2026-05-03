"""
make_conus_mask.py  —  Generate usreg_half.txt

Creates a 116×50 binary mask on the 0.5° grid covering the contiguous US.
Grid cell centers:
    lon: ix*0.5 - 124.75  for ix=1..116  → -124.25°W to -67.25°W
    lat: iy*0.5 + 24.75   for iy=1..50   →  25.25°N  to  49.75°N

Output format (matches what us_dly_waves.f expects):
    50 rows, north→south (iy=50 first, iy=1 last)
    Each row: 116 space-separated integers (1=inside CONUS, 0=outside)

Run once; output file is consumed by us_dly_waves.py for IDW gridding.
"""

import numpy as np
import geopandas as gpd
from shapely.geometry import Point
import os

NX = 116
NY = 50
OUT_FILE = "usreg_half.txt"


def grid_centers():
    """Return (lons, lats) 1-D arrays for ix=1..116, iy=1..50."""
    lons = np.array([ix * 0.5 - 124.75 for ix in range(1, NX + 1)])
    lats = np.array([iy * 0.5 + 24.75  for iy in range(1, NY + 1)])
    return lons, lats


def load_conus_boundary_census():
    """
    Fallback: download Census Bureau CONUS shapefile on first call and cache it.
    Returns a shapely geometry for the 48 states.
    """
    import urllib.request, zipfile, pathlib, tempfile

    cache = pathlib.Path("cb_2018_us_nation_5m")
    if not cache.exists():
        url = ("https://www2.census.gov/geo/tiger/GENZ2018/shp/"
               "cb_2018_us_nation_5m.zip")
        print(f"Downloading CONUS boundary from Census Bureau…")
        with urllib.request.urlopen(url) as resp:
            data = resp.read()
        with zipfile.ZipFile(__import__("io").BytesIO(data)) as zf:
            zf.extractall(cache)
    gdf = gpd.read_file(list(cache.glob("*.shp"))[0])
    geom = gdf.geometry.iloc[0]
    # Filter to mainland: centroids inside CONUS lon/lat box
    from shapely.geometry import MultiPolygon
    if geom.geom_type == "MultiPolygon":
        parts = [p for p in geom.geoms
                 if -130 < p.centroid.x < -60 and 24 < p.centroid.y < 50]
        geom = MultiPolygon(parts) if len(parts) > 1 else parts[0]
    return geom


def build_mask(conus_geom):
    """Return (NY, NX) int array: 1 inside CONUS, 0 outside."""
    lons, lats = grid_centers()
    mask = np.zeros((NY, NX), dtype=int)
    for iy in range(NY):
        for ix in range(NX):
            pt = Point(lons[ix], lats[iy])
            if conus_geom.contains(pt):
                mask[iy, ix] = 1
    return mask


def write_mask(mask):
    """Write usreg_half.txt: 50 rows north→south, each 116 ints."""
    with open(OUT_FILE, "w") as f:
        for iy in range(NY - 1, -1, -1):   # iy=49 (north) down to iy=0 (south)
            row = " ".join(str(v) for v in mask[iy])
            f.write(row + "\n")
    print(f"Wrote {OUT_FILE}  ({NY} rows × {NX} cols)")


def report(mask):
    n_cells = int(mask.sum())
    total = NY * NX
    print(f"CONUS cells: {n_cells} / {total}  ({100*n_cells/total:.1f}%)")
    # Rough sanity: CONUS is ~7.6M km²; 0.5° cell at 37°N ≈ ~2700 km²
    # expect ~2500–3200 cells inside
    if not (2000 < n_cells < 3500):
        print("WARNING: cell count outside expected range 2000–3500 — check boundary.")


if __name__ == "__main__":
    print("Loading CONUS boundary…")
    conus = load_conus_boundary_census()
    print("  source: Census Bureau cb_2018_us_nation_5m (cached in cb_2018_us_nation_5m/)")

    print("Testing grid cells…")
    mask = build_mask(conus)
    report(mask)
    write_mask(mask)

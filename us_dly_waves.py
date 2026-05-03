#!/usr/bin/env python3
"""
Heat wave / cold wave frequency analysis for USHCN daily temperature data.
Translated from us_dly_waves.f (J.R. Christy, UAH/NSSTC, April 2026).

Usage:
    python us_dly_waves.py --var tmax --percentile 95 --min-run 3
    python us_dly_waves.py --var tmin --percentile 5  --min-run 4 --data-dir /path/to/data
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ── Constants ────────────────────────────────────────────────────────────────

NYB = 1899          # first year in dataset
NYE = 2025          # last year in dataset
N_YEARS = NYE - NYB + 1   # 127

NX = 116                # grid columns (lon)
NY = 50                 # grid rows (lat)
IDW_RADIUS_KM = 115.0   # Fortran uses 115 km search radius

# COOP state numeric codes 1-48 → 2-letter abbreviations (Fortran data ast).
# Index 48 (state "49") = CONUS average (set from regional grid, unavailable here).
STATE_ABBR = [
    'AL', 'AZ', 'AR', 'CA', 'CO', 'CN', 'DE', 'FL',   # 1-8
    'GA', 'ID', 'IL', 'IN', 'IA', 'KS', 'KE', 'LA',   # 9-16
    'ME', 'MD', 'MA', 'MI', 'MN', 'MS', 'MO', 'MT',   # 17-24
    'NB', 'NV', 'NH', 'NJ', 'NM', 'NY', 'NC', 'ND',   # 25-32
    'OH', 'OK', 'OR', 'PE', 'RI', 'SC', 'SD', 'TN',   # 33-40
    'TX', 'UT', 'VE', 'VA', 'WA', 'WV', 'WI', 'WY',   # 41-48
    'US',                                               # 49 (CONUS)
]

# Days per month (index 0 = Dec, 1 = Jan, …, 12 = Dec).
_MONTH_DAYS = [31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description='USHCN heat/cold wave frequency counter'
    )
    p.add_argument('--var', choices=['tmax', 'tmin'], required=True,
                   help='tmax → heat waves (≥ threshold);  tmin → cold waves (≤ threshold)')
    p.add_argument('--percentile', type=int, required=True,
                   help='Percentile threshold, e.g. 95 for heat waves or 5 for cold waves')
    p.add_argument('--min-run', type=int, required=True,
                   help='Minimum consecutive days to qualify as a wave event')
    p.add_argument('--data-dir', type=Path, default=Path('.'),
                   help='Folder containing data files (default: current directory)')
    return p.parse_args()


# ── Calendar helpers ─────────────────────────────────────────────────────────

def _is_leap(year: int) -> bool:
    # Fortran uses the simple divisible-by-4 rule (no century correction).
    return year % 4 == 0


def _days_in_month(year: int, month: int) -> int:
    if month == 2 and _is_leap(year):
        return 29
    return _MONTH_DAYS[month]


# ── Data loading ──────────────────────────────────────────────────────────────

def load_stations(data_dir: Path) -> dict:
    """Return {coop_id (int) -> (lat, lon)} from ushcn-v2-stations.txt."""
    result = {}
    with open(data_dir / 'ushcn-v2-stations.txt') as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    result[int(parts[0])] = (float(parts[1]), float(parts[2]))
                except ValueError:
                    pass
    return result


def load_temperature(data_dir: Path, var: str) -> pd.DataFrame:
    """Read tmax or tmin file into a DataFrame.

    Columns: stn_id, year, month, d1 … d31 (NaN where original value is -999).
    """
    fname = ('ushcn_jrc_tmax_260421.txt' if var == 'tmax'
             else 'ushcn_jrc_tmin_260421.txt')
    day_cols = [f'd{i}' for i in range(1, 32)]
    df = pd.read_csv(
        data_dir / fname,
        sep=r'\s+',
        header=None,
        names=['stn_id', 'year', 'month'] + day_cols,
        dtype={c: int for c in ['stn_id', 'year', 'month']},
    )
    # Replace sentinel missing value with NaN
    df[day_cols] = df[day_cols].where(df[day_cols] != -999, other=np.nan)
    return df


def station_order(df: pd.DataFrame) -> list:
    """Return station IDs in order of first appearance in the temperature file."""
    seen: set = set()
    order: list = []
    for sid in df['stn_id'].values:
        if sid not in seen:
            seen.add(sid)
            order.append(int(sid))
    return order


# ── Per-station seasonal array ────────────────────────────────────────────────

def build_jval(stn_df: pd.DataFrame, var: str) -> np.ndarray:
    """Build flat seasonal day array jval[year_idx, day_idx] for one station.

    Shape: (N_YEARS, jtot) where jtot = 122 (tmin) or 153 (tmax).
    NaN for missing days.

    tmin season: Dec (month 0) – Jan – Feb – Mar  (months [0,1,2,3])
                 December of year Y is remapped to month 0 of year Y+1.
    tmax season: May – Jun – Jul – Aug – Sep       (months [5,6,7,8,9])
    """
    if var == 'tmin':
        season_months = [0, 1, 2, 3]
        jtot = 122
    else:
        season_months = [5, 6, 7, 8, 9]
        jtot = 153

    day_cols = [f'd{i}' for i in range(1, 32)]

    # Build (year, month) → np.array(31,) lookup from the station's records.
    lookup: dict = {}
    years_np = stn_df['year'].values
    months_np = stn_df['month'].values
    data_np = stn_df[day_cols].values.astype(float)

    for k in range(len(years_np)):
        yr, mo = int(years_np[k]), int(months_np[k])
        lookup[(yr, mo)] = data_np[k]
        if var == 'tmin' and mo == 12:
            # December of year Y becomes month 0 of year Y+1 (Fortran lines 238-243).
            lookup[(yr + 1, 0)] = data_np[k]

    jval = np.full((N_YEARS, jtot), np.nan, dtype=np.float32)

    for yi in range(N_YEARS):
        year = NYB + yi
        pos = 0
        for mo in season_months:
            nd = _days_in_month(year, mo)
            arr = lookup.get((year, mo))
            end = min(pos + nd, jtot)
            if arr is not None:
                jval[yi, pos:end] = arr[:end - pos]
            pos += nd
            if pos >= jtot:
                break

    return jval


# ── Threshold computation ─────────────────────────────────────────────────────

def compute_thresholds(jval: np.ndarray, nperc: int) -> np.ndarray:
    """Compute per-day-of-season percentile thresholds (Fortran lines 275-308).

    For each calendar day idy, collect all non-missing values from a ±3-day
    window across all years, sort them, and take the nperc-th percentile.
    Returns ntval[jtot] (float32; NaN where fewer than 200 values exist).
    """
    _, jtot = jval.shape
    ntval = np.full(jtot, np.nan, dtype=np.float32)

    for idy in range(jtot):
        jd1 = max(0, idy - 3)
        jd2 = min(jtot - 1, idy + 3)
        window = jval[:, jd1:jd2 + 1].ravel()
        valid = window[~np.isnan(window)]
        icnt = len(valid)
        if icnt > 200:
            sorted_vals = np.sort(valid)
            # Fortran (1-indexed): xsort(int(icnt*nperc*0.01 + 1))
            # Python  (0-indexed): sorted_vals[int(icnt*nperc*0.01)]
            idx = min(int(icnt * nperc * 0.01), icnt - 1)
            ntval[idy] = sorted_vals[idx]

    return ntval


# ── Run counting ──────────────────────────────────────────────────────────────

def count_waves(
    jval: np.ndarray,
    ntval: np.ndarray,
    var: str,
    nrun: int,
    stn_id: int,
) -> np.ndarray:
    """Count wave days per year (Fortran lines 313-384).

    Returns nval[N_YEARS] (int32): total wave days per year.
    Years with < 70 % valid data are set to -99.
    Per-event detail lines are printed to stdout for years > 1984.
    """
    n_years, jtot = jval.shape
    nval = np.full(n_years, -99, dtype=np.int32)

    # Boolean flag array: True where day meets wave criterion.
    # NaN threshold or NaN data → False (not a wave day).
    with np.errstate(invalid='ignore'):
        if var == 'tmax':
            kval = (~np.isnan(ntval) & ~np.isnan(jval)
                    & (jval >= ntval)).astype(np.int8)
        else:
            kval = (~np.isnan(ntval) & ~np.isnan(jval)
                    & (jval <= ntval)).astype(np.int8)

    valid_per_year = np.sum(~np.isnan(jval), axis=1)
    min_valid = int(jtot * 0.7)

    for yi in range(n_years):
        if valid_per_year[yi] < min_valid:
            continue

        year = NYB + yi
        k = kval[yi]

        # Identify run starts (0→1 transition) and ends (1→0 transition).
        # Pad with 0 at both ends so edge runs are detected correctly.
        padded = np.empty(jtot + 2, dtype=np.int8)
        padded[0] = 0
        padded[1:-1] = k
        padded[-1] = 0
        diff = np.diff(padded.astype(np.int16))
        starts = np.where(diff == 1)[0]    # 0-indexed start of each run in k
        ends   = np.where(diff == -1)[0]   # 0-indexed: first position after run

        total = 0
        for s, e in zip(starts, ends):
            run_len = int(e - s)           # matches Fortran icnt
            if run_len >= nrun:
                total += run_len
                if year > 1984:
                    for jdy in range(s, e):
                        v = int(jval[yi, jdy]) if not np.isnan(jval[yi, jdy]) else -999
                        t = int(ntval[jdy])    if not np.isnan(ntval[jdy])    else -999
                        print(f'{stn_id:7d}{year:5d}{jdy + 1:4d}{run_len:4d}{v:4d}{t:4d}')

        nval[yi] = total

    return nval


# ── Aggregation ───────────────────────────────────────────────────────────────

def state_averages(nval_all: np.ndarray, stn_order: list) -> np.ndarray:
    """Average wave days per state (Fortran lines 390-409).

    Returns xsval[48, N_YEARS] (float32).  States are indexed 0-47 (= COOP
    state numbers 1-48).  Years with no valid stations get -99.
    """
    xsval = np.full((48, N_YEARS), -99.0, dtype=np.float32)

    # Group station array indices by COOP state number (coop_id // 10000).
    state_idx: dict = {}
    for i, sid in enumerate(stn_order):
        st = sid // 10000
        if 1 <= st <= 48:
            state_idx.setdefault(st, []).append(i)

    for st, indices in state_idx.items():
        arr = nval_all[indices, :]          # shape (n_stns_in_state, N_YEARS)
        mask = arr > -90                    # valid entries
        counts = mask.sum(axis=0).astype(float)
        totals = np.where(mask, arr, 0).sum(axis=0).astype(float)
        valid_yrs = counts > 0
        xsval[st - 1, valid_yrs] = totals[valid_yrs] / counts[valid_yrs]

    return xsval


# ── IDW CONUS gridding ────────────────────────────────────────────────────────

def load_mask(data_dir: Path):
    """Read usreg_half.txt → (NY, NX) int8 array, or None if file is absent."""
    path = data_dir / 'usreg_half.txt'
    if not path.exists():
        return None
    mask = np.zeros((NY, NX), dtype=np.int8)
    with open(path) as f:
        for j, line in enumerate(f):
            if j >= NY:
                break
            mask[j] = np.fromstring(line, sep=' ', dtype=np.int8)
    return mask


def _haversine_km(lat1, lon1, lat2, lon2):
    """Vectorized haversine distance in km (degrees input, broadcasts freely)."""
    R = 6371.0
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = (np.sin(dlat / 2) ** 2
         + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon / 2) ** 2)
    return R * 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))


def build_idw_weights(mask: np.ndarray, stations: dict, stn_order: list) -> np.ndarray:
    """Precompute IDW weight matrix (n_conus_cells, n_stns).

    For each CONUS grid cell, weight[cell, stn] = (radius/dist)^2 if dist ≤ radius,
    else 0.  Stations without lat/lon metadata get weight 0 everywhere.

    Grid cell coordinates (j=0 is northernmost row in file):
        lat = 49.75 - j*0.5
        lon = i*0.5 - 124.25
    """
    stn_lats = np.array([stations.get(s, (np.nan, np.nan))[0] for s in stn_order])
    stn_lons = np.array([stations.get(s, (np.nan, np.nan))[1] for s in stn_order])

    cell_lats, cell_lons = [], []
    for j in range(NY):
        for i in range(NX):
            if mask[j, i] > 0:
                cell_lats.append(49.75 - j * 0.5)
                cell_lons.append(i * 0.5 - 124.25)

    clat = np.array(cell_lats)[:, None]   # (n_cells, 1)
    clon = np.array(cell_lons)[:, None]

    dist = _haversine_km(clat, clon, stn_lats[None, :], stn_lons[None, :])
    # (n_cells, n_stns)

    weights = np.where(
        dist <= IDW_RADIUS_KM,
        (IDW_RADIUS_KM / np.maximum(dist, 0.1)) ** 2,
        0.0,
    )
    return weights.astype(np.float32)


def conus_idw_average(nval_all: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """IDW CONUS average wave days per year.

    nval_all : (n_stns, N_YEARS), -99 sentinel for missing years
    weights  : (n_cells, n_stns) precomputed from build_idw_weights()
    Returns  : conus_avg[N_YEARS] float32, -99 where no data
    """
    conus_avg = np.full(N_YEARS, -99.0, dtype=np.float32)
    weights = weights.astype(np.float64)

    for yi in range(N_YEARS):
        vals = nval_all[:, yi].astype(np.float64)
        valid = vals > -90
        vals[~valid] = 0.0

        w = weights * valid[None, :]          # zero-out missing stations
        w_sum = w.sum(axis=1)                 # (n_cells,)
        cell_vals = np.where(
            w_sum > 0,
            (w * vals[None, :]).sum(axis=1) / w_sum,
            np.nan,
        )

        good = ~np.isnan(cell_vals)
        if good.any():
            conus_avg[yi] = float(cell_vals[good].mean())

    return conus_avg


# ── Output printers ───────────────────────────────────────────────────────────

def print_station_table(nval_all: np.ndarray, stn_order: list) -> None:
    """Print per-station heat/cold wave day counts (Fortran format 300/301)."""
    n = len(stn_order)
    for start in range(0, n, 32):
        end = min(start + 32, n)
        chunk = stn_order[start:end]
        print('     ' + ''.join(f'{sid:7d}' for sid in chunk))
        for yi, year in enumerate(range(NYB, NYE + 1)):
            print(f'{year:5d}'
                  + ''.join(f'{int(nval_all[start + k, yi]):7d}'
                             for k in range(len(chunk))))
        print()


def print_state_table(
    xsval: np.ndarray,
    conus_avg: np.ndarray,
    nperc: int,
    nrun: int,
) -> None:
    """Print per-state averages (Fortran format 200/201/101/102)."""
    abbr = STATE_ABBR   # 49 entries, index 0-47 = states 1-48, index 48 = US

    print()
    print(f' State average {nperc} percentile {nrun} days')

    # Block 1: states 1-24 (indices 0-23)
    hdr1 = '     ' + ''.join(f'   {abbr[i]:2s} ' for i in range(24))
    print(hdr1)
    for yi, year in enumerate(range(NYB, NYE + 1)):
        print(f'{year:5d}' + ''.join(f'{xsval[i, yi]:6.2f}' for i in range(24)))
    print(hdr1)

    print()

    # Block 2: states 25-48 + CONUS IDW average (indices 24-47, then 48='US')
    hdr2 = ''.join(f'   {abbr[i]:2s} ' for i in range(24, 49))
    print(hdr2)
    for yi in range(N_YEARS):
        vals = list(xsval[24:48, yi])
        vals.append(float(conus_avg[yi]))
        print(''.join(f'{v:6.2f}' for v in vals))
    print(hdr2)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()
    data_dir = args.data_dir
    var      = args.var
    nperc    = args.percentile
    nrun     = args.min_run

    print('Loading station metadata …', file=sys.stderr)
    stations = load_stations(data_dir)
    print(f'  {len(stations)} stations in ushcn-v2-stations.txt', file=sys.stderr)

    print(f'Loading {var} temperature data …', file=sys.stderr)
    df_temp = load_temperature(data_dir, var)

    stn_order = station_order(df_temp)
    n_stns = len(stn_order)
    print(f'  {n_stns} stations in temperature file', file=sys.stderr)

    # Pre-group temperature records by station for fast per-station access.
    grouped = {int(sid): grp.reset_index(drop=True)
               for sid, grp in df_temp.groupby('stn_id')}

    nval_all = np.full((n_stns, N_YEARS), -99, dtype=np.int32)

    print('Processing stations …', file=sys.stderr)
    for i, stn_id in enumerate(stn_order):
        if (i + 1) % 200 == 0 or i == 0:
            print(f'  {i + 1}/{n_stns}', file=sys.stderr)

        stn_df = grouped.get(stn_id)
        if stn_df is None:
            continue

        jval  = build_jval(stn_df, var)
        ntval = compute_thresholds(jval, nperc)
        nval  = count_waves(jval, ntval, var, nrun, stn_id)
        nval_all[i] = nval

    print(f'  {n_stns}/{n_stns} done', file=sys.stderr)

    # ── IDW CONUS average ────────────────────────────────────────────────────
    mask = load_mask(data_dir)
    if mask is not None:
        print('Building IDW weight matrix …', file=sys.stderr)
        weights = build_idw_weights(mask, stations, stn_order)
        print('Computing IDW CONUS average …', file=sys.stderr)
        conus_avg = conus_idw_average(nval_all, weights)
    else:
        print('  usreg_half.txt not found — CONUS average set to -99', file=sys.stderr)
        print('  Run make_conus_mask.py to generate it.', file=sys.stderr)
        conus_avg = np.full(N_YEARS, -99.0, dtype=np.float32)

    # ── Output ──────────────────────────────────────────────────────────────
    print_station_table(nval_all, stn_order)

    xsval = state_averages(nval_all, stn_order)
    print_state_table(xsval, conus_avg, nperc, nrun)


if __name__ == '__main__':
    main()

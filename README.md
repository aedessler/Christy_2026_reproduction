# Christy (2026) — USHCN Daily Temperature Extremes Analysis

Replication and extension of:

> Christy, J. R. (2026). Declines in hot and cold daily temperature extremes in the
> conterminous US, 1899–2025. *Theoretical and Applied Climatology*, 157, 309.
> https://doi.org/10.1007/s00704-026-06200-3

---

## Overview

This directory contains data and Python scripts to reproduce the daily temperature
extreme analyses from Christy (2026). The study examines hot and cold extremes for
1,218 USHCN stations across the contiguous US (CONUS) over 1899–2025.

The original Fortran code was downloaded from [UAH/NSSTC](https://www.nsstc.uah.edu/data/ushcn_jrc/)
and converted into the Python scripts here.

Note that the original Fortran performs IDW regridding of station data to a 0.5° grid
before computing regional averages, but that requires `usreg_half.txt` which was not
included in the data archive. The Python version skips that step and uses a simple
station average instead, which closely reproduces Figures 3 and 4 from the paper.

The Python port also corrects an off-by-one bug in the Fortran's per-event output: the
original printed one extra (non-wave) day at the end of each reported event.

---

## Files

### Data

| File | Size | Description |
|------|------|-------------|
| `ushcn_jrc_tmax_260421.txt` | 128 MB | Daily T_Max (May–Sep, 1899–2025), 1,218 stations |
| `ushcn_jrc_tmin_260421.txt` | 102 MB | Daily T_Min (Dec–Mar, 1899–2025), 1,218 stations |
| `ushcn-v2-stations.txt` | 108 KB | Station metadata: COOP ID, latitude, longitude |
| `ushcn_stn_list_250617.txt` | 47 KB | Station list with state abbreviations |
| `ushcn_jrc_readme_260421.txt` | 7 KB | Data format documentation |
| `US_MidW_61stn_dly_ppt_1893_2024_250704.xlsx` | 12 MB | Midwest precipitation dataset (61 stations) |

> The two large temperature files (`ushcn_jrc_tmax_260421.txt` and
> `ushcn_jrc_tmin_260421.txt`) are excluded from the git repository due to GitHub's
> 100 MB file size limit. Download them from the
> [UAH/NSSTC data archive](https://www.nsstc.uah.edu/data/ushcn_jrc/).

**Data format** (temperature files): space-delimited, one record per station/year/month.

```
STNID  YEAR  MM  D1  D2  D3 … D31
11084  1899   5  85  83  87 …  88
```

Values are integer °F; `-999` = missing. Station IDs are 5-digit COOP IDs (leading zero
omitted for states 01–09).

### Scripts

| File | Description |
|------|-------------|
| `us_dly_waves.py` | Heat/cold wave frequency analysis (Python port of Fortran original) |
| `figures.py` | Reproduces Fig. 3 and Fig. 4 from the paper |

### Figures

| File | Description |
|------|-------------|
| `Fig3.jpg` | Reproduced Fig. 3 (all-time record fractions by year) |
| `Fig4.jpg` | Reproduced Fig. 4 (daily records per station per year) |

---

## Requirements

```
numpy
pandas
matplotlib
```

---

## Usage

### Heat/cold wave analysis (`us_dly_waves.py`)

Counts the number of heat wave (or cold wave) days per station per year. A "wave" is
defined as a consecutive run of days at or above the N-th percentile threshold (or at or
below for cold), where the threshold is computed from a ±3-day window across all years.

```bash
# Heat waves: T_Max ≥ 95th percentile for ≥ 3 consecutive days
python us_dly_waves.py --var tmax --percentile 95 --min-run 3 > output_tmax.txt

# Cold waves: T_Min ≤ 5th percentile for ≥ 4 consecutive days
python us_dly_waves.py --var tmin --percentile 5 --min-run 4 > output_tmin.txt
```

**Arguments:**

| Flag | Description |
|------|-------------|
| `--var` | `tmax` (heat waves) or `tmin` (cold waves) |
| `--percentile` | Integer threshold percentile (e.g. 95 or 5) |
| `--min-run` | Minimum consecutive days to qualify as a wave |
| `--data-dir` | Path to data files (default: current directory) |

**Output** (to stdout):
1. Per-event detail lines for years ≥ 1985: `stn_id year day run_len value threshold`
2. Station-by-station annual wave day counts (all years 1899–2025)
3. State averages for all 48 continental US states

> **Note:** Regional (CONUS grid) averages are skipped because `usreg_half.txt`
> (a 116×50 half-degree US region grid) is not included in the data archive.

---

### Figure reproduction (`figures.py`)

Reproduces Fig. 3 and Fig. 4 from Christy (2026) and writes `Fig3.jpg` and `Fig4.jpg`.

```bash
python figures.py
```

**Fig. 3 — All-time record fractions**
For each station, the year(s) in which its all-time hottest (T_Max) or coldest (T_Min)
temperature occurred receive a proportional share of 1.0. Values are averaged across
all stations to give the CONUS fraction per year.

**Fig. 4 — Daily records per station per year**
For each of the 153 (T_Max) or 122 (T_Min) seasonal days, the year with the highest
(or lowest) value receives one record count. Counts are averaged across stations.
An 11-year centered running average is overlaid.

> The paper uses IDW spatial interpolation to a 0.5° grid before computing CONUS
> averages (requires `usreg_half.txt`). Here a simple station average is used instead,
> which produces nearly identical results.

---

## Key results (from the paper, verified by these scripts)

| Metric | Value | Year |
|--------|-------|------|
| Largest CONUS hot-record fraction (Fig. 3) | ~22% | 1936 |
| Largest CONUS cold-record fraction (Fig. 3) | ~22% | 1899 |
| Most daily T_Max records per station (Fig. 4) | 6.7 | 1936 |
| Most daily T_Min records per station (Fig. 4) | 3.7 | 1899 |

Scripts reproduce these within ~5% (difference attributable to station vs. grid averaging).

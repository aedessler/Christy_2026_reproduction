#!/usr/bin/env python3
"""
Reproduces Fig. 3 and Fig. 4 from:
  Christy (2026) "Declines in hot and cold daily temperature extremes
  in the conterminous US, 1899–2025." Theor. Appl. Climatol. 157:309.

Fig. 3 – Proportion of CONUS stations with all-time record hot/cold temperature by year.
         Method: for each station find its all-time record value; distribute 1.0
         equally among all years in which that record was achieved; average across
         stations.

Fig. 4 – Per-station average number of daily TMax / TMin records per year, with
         11-year centered running average.
         Method: for each station and each seasonal day find the year with the
         hottest (coldest) value; count records per year; average across stations.

Outputs: Fig3.jpg and Fig4.jpg
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

DATA_DIR = Path(__file__).parent
NYB, NYE = 1899, 2025
YEARS     = np.arange(NYB, NYE + 1)   # 1899 … 2025
N_YEARS   = len(YEARS)                 # 127

RED  = '#CC2222'   # hot / TMax
BLUE = '#4488CC'   # cold / TMin

# ── Calendar helpers ──────────────────────────────────────────────────────────
# Index 0 = Dec (remapped), 1 = Jan, …, 12 = Dec (standard).
_MDAYS = [31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

def _ndays(year: int, month: int) -> int:
    return 29 if (month == 2 and year % 4 == 0) else _MDAYS[month]


# ── Data loading ──────────────────────────────────────────────────────────────

def load_temp(var: str) -> pd.DataFrame:
    """Read tmax or tmin file; replace -999 with NaN."""
    fname = ('ushcn_jrc_tmax_260421.txt' if var == 'tmax'
             else 'ushcn_jrc_tmin_260421.txt')
    dcols = [f'd{i}' for i in range(1, 32)]
    df = pd.read_csv(
        DATA_DIR / fname, sep=r'\s+', header=None,
        names=['stn', 'yr', 'mo'] + dcols,
        dtype={'stn': int, 'yr': int, 'mo': int},
    )
    df[dcols] = df[dcols].where(df[dcols] != -999, other=np.nan)
    return df


# ── Per-station seasonal day array ───────────────────────────────────────────

def jval_array(stn_df: pd.DataFrame, var: str) -> np.ndarray:
    """Return (N_YEARS, jtot) seasonal day array for one station; NaN = missing.

    tmin season : Dec(month 0, remapped) – Jan – Feb – Mar  [jtot = 122]
    tmax season : May – Jun – Jul – Aug – Sep               [jtot = 153]
    """
    seas = [0, 1, 2, 3] if var == 'tmin' else [5, 6, 7, 8, 9]
    jtot = 122             if var == 'tmin' else 153
    dcols = [f'd{i}' for i in range(1, 32)]

    raw  = stn_df[['yr', 'mo'] + dcols].values
    lk: dict = {}
    for row in raw:
        yr, mo = int(row[0]), int(row[1])
        vals = row[2:].astype(float)
        lk[(yr, mo)] = vals
        if var == 'tmin' and mo == 12:      # Dec of year Y → month 0 of year Y+1
            lk[(yr + 1, 0)] = vals

    out = np.full((N_YEARS, jtot), np.nan, dtype=np.float32)
    for yi in range(N_YEARS):
        yr  = NYB + yi
        pos = 0
        for mo in seas:
            nd  = _ndays(yr, mo)
            arr = lk.get((yr, mo))
            end = min(pos + nd, jtot)
            if arr is not None:
                out[yi, pos:end] = arr[:end - pos]
            pos += nd
            if pos >= jtot:
                break
    return out


# ── Metrics (single pass over stations) ──────────────────────────────────────

def compute_metrics(df: pd.DataFrame, var: str):
    """Compute Fig. 3 and Fig. 4 metrics in one pass.

    Returns:
        fracs   – Fig. 3: CONUS-average fraction with all-time record [N_YEARS]
        counts  – Fig. 4: per-station average daily record count [N_YEARS]
    """
    fracs  = np.zeros(N_YEARS)
    counts = np.zeros(N_YEARS)
    n_stn  = 0

    for _, grp in df.groupby('stn'):
        jv = jval_array(grp, var)           # (N_YEARS, jtot)

        alltime = np.nanmax(jv) if var == 'tmax' else np.nanmin(jv)
        if np.isnan(alltime):
            continue

        # ── Fig. 3: all-time single-event record ──────────────────────────
        # Find years where any seasonal day equalled the all-time record.
        has_alltime = np.any(jv == alltime, axis=1)    # (N_YEARS,)
        n_ties = int(has_alltime.sum())
        if n_ties:
            fracs[has_alltime] += 1.0 / n_ties

        # ── Fig. 4: per-day-of-season records ─────────────────────────────
        # For each of the jtot seasonal days, find the record value across years.
        rec_per_day = np.nanmax(jv, axis=0) if var == 'tmax' else np.nanmin(jv, axis=0)
        with np.errstate(invalid='ignore'):
            is_rec = (jv == rec_per_day[np.newaxis, :]) & ~np.isnan(jv)
        # Distribute count equally among tied years for each day.
        n_tied_per_day = is_rec.sum(axis=0).astype(float)
        n_tied_per_day[n_tied_per_day == 0] = 1.0      # avoid div-by-zero for NaN days
        weighted = is_rec / n_tied_per_day[np.newaxis, :]
        counts += weighted.sum(axis=1)                 # sum over days → (N_YEARS,)

        n_stn += 1

    if n_stn:
        fracs  /= n_stn
        counts /= n_stn

    return fracs, counts


def centered_mean(arr: np.ndarray, k: int) -> np.ndarray:
    """k-point centered running mean; NaN for the (k//2) years at each edge."""
    h   = k // 2
    out = np.full(len(arr), np.nan, dtype=float)
    for i in range(h, len(arr) - h):
        out[i] = np.nanmean(arr[i - h : i + h + 1])
    return out


# ── Figure 3 ──────────────────────────────────────────────────────────────────

def make_fig3(frac_hot: np.ndarray, frac_cold: np.ndarray) -> plt.Figure:
    """Proportion of CONUS with all-time record hot (above) / cold (below)."""
    fig, ax = plt.subplots(figsize=(10, 5.5))

    ax.bar(YEARS,  frac_hot,  width=0.85, color=RED,  zorder=3)
    ax.bar(YEARS, -frac_cold, width=0.85, color=BLUE, zorder=3)
    ax.axhline(0, color='black', linewidth=0.8, zorder=2)

    # Axes limits and ticks
    ax.set_xlim(1893, 2027)
    ax.set_xticks(np.arange(1895, 2026, 10))
    ax.tick_params(axis='x', labelsize=9)

    ax.set_ylim(-0.26, 0.26)
    ax.set_yticks(np.round(np.arange(-0.25, 0.251, 0.05), 2))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f'{abs(y):.2f}'))
    ax.set_ylabel('Fraction of CONUS with All Time Record Temperature', fontsize=9)
    ax.tick_params(axis='y', labelsize=8.5)

    # Horizontal gridlines
    ax.yaxis.grid(True, color='#CCCCCC', linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[['top', 'right']].set_visible(False)

    # Annotation boxes
    ax.text(1950, 0.125, 'Record Hot',  color=RED,  fontsize=10.5, fontweight='bold',
            ha='left', va='center',
            bbox=dict(facecolor='white', edgecolor=RED,  boxstyle='square,pad=0.3', lw=1.5))
    ax.text(1950, -0.135, 'Record Cold', color=BLUE, fontsize=10.5, fontweight='bold',
            ha='left', va='center',
            bbox=dict(facecolor='white', edgecolor=BLUE, boxstyle='square,pad=0.3', lw=1.5))

    # Legend
    ax.legend(
        handles=[Patch(fc=RED, label='CONUS_TMax'), Patch(fc=BLUE, label='CONUS_TMin')],
        loc='lower center', ncol=2, fontsize=9, frameon=True, framealpha=1.0,
        edgecolor='#AAAAAA', bbox_to_anchor=(0.5, 0.01),
    )

    # Caption below figure
    fig.text(
        0.08, 0.005,
        'Fig. 3  Proportion of the CONUS experiencing the All-Time extreme hottest '
        'and coldest day by year for 1899–2025',
        fontsize=8, style='italic', va='bottom',
    )

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    return fig


# ── Figure 4 ──────────────────────────────────────────────────────────────────

def make_fig4(cnt_hot: np.ndarray, cnt_cold: np.ndarray) -> plt.Figure:
    """Per-station daily TMax/TMin record count per year + 11-yr running mean."""
    rm_hot  = centered_mean(cnt_hot,  11)
    rm_cold = centered_mean(cnt_cold, 11)

    fig, ax = plt.subplots(figsize=(10, 5.5))

    ax.bar(YEARS,  cnt_hot,  width=0.85, color=RED,  zorder=3, label='CONUS Mx')
    ax.bar(YEARS, -cnt_cold, width=0.85, color=BLUE, zorder=3, label='CONUS Mi')
    ax.plot(YEARS,  rm_hot,  color='black',   linewidth=1.6, zorder=4,
            label='CONUS Mx 11-yr Avg')
    ax.plot(YEARS, -rm_cold, color='#001F7A', linewidth=1.6, zorder=4,
            label='CONUS Mi 11-yr Avg')
    ax.axhline(0, color='black', linewidth=0.8, zorder=2)

    # Axes limits and ticks
    ax.set_xlim(1893, 2027)
    ax.set_xticks(np.arange(1895, 2026, 10))
    ax.tick_params(axis='x', labelsize=9)

    # Asymmetric y-axis: +8 above, -4 below; labels show absolute values
    ax.set_ylim(-4.4, 8.8)
    ax.set_yticks([-4, -2, 0, 2, 4, 6, 8])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: str(int(abs(y)))))
    ax.set_ylabel('Per station Average Number of records per year', fontsize=9)
    ax.tick_params(axis='y', labelsize=8.5)

    # Horizontal gridlines
    ax.yaxis.grid(True, color='#CCCCCC', linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.spines[['top', 'right']].set_visible(False)

    # Legend (4 elements in a row)
    ax.legend(
        handles=[
            Patch(fc=RED,      label='CONUS Mx'),
            Patch(fc=BLUE,     label='CONUS Mi'),
            Line2D([0], [0], color='black',   linewidth=1.6, label='CONUS Mx 11-yr Avg'),
            Line2D([0], [0], color='#001F7A', linewidth=1.6, label='CONUS Mi 11-yr Avg'),
        ],
        loc='lower center', ncol=4, fontsize=9, frameon=True, framealpha=1.0,
        edgecolor='#AAAAAA', bbox_to_anchor=(0.5, 0.01),
    )

    # Caption below figure
    fig.text(
        0.08, 0.005,
        'Fig. 4  Average number of daily Tₘₐₓ (red) and Tₘᴵₙ (blue) records achieved '
        'in each year. Lines represent the 11-year centered average',
        fontsize=8, style='italic', va='bottom',
    )

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    return fig


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    print('Loading tmax data …',  file=sys.stderr)
    df_tmax = load_temp('tmax')
    print('Loading tmin data …',  file=sys.stderr)
    df_tmin = load_temp('tmin')

    print('Computing tmax metrics (Fig. 3 + Fig. 4) …', file=sys.stderr)
    frac_hot,  cnt_hot  = compute_metrics(df_tmax, 'tmax')
    print('Computing tmin metrics (Fig. 3 + Fig. 4) …', file=sys.stderr)
    frac_cold, cnt_cold = compute_metrics(df_tmin, 'tmin')

    # Quick sanity check against paper values
    yr1936 = 1936 - NYB
    yr1899 = 1899 - NYB
    print(f'  Fig. 3 check – 1936 hot frac: {frac_hot[yr1936]:.3f}  '
          f'(paper ~0.22)',                  file=sys.stderr)
    print(f'  Fig. 3 check – 1899 cold frac: {frac_cold[yr1899]:.3f} '
          f'(paper ~0.22)',                  file=sys.stderr)
    print(f'  Fig. 4 check – 1936 hot recs/stn: {cnt_hot[yr1936]:.2f}  '
          f'(paper ~6.7)',                   file=sys.stderr)
    print(f'  Fig. 4 check – 1899 cold recs/stn: {cnt_cold[yr1899]:.2f} '
          f'(paper ~3.7)',                   file=sys.stderr)

    out3 = DATA_DIR / 'Fig3.jpg'
    out4 = DATA_DIR / 'Fig4.jpg'

    print(f'Saving {out3} …', file=sys.stderr)
    make_fig3(frac_hot, frac_cold).savefig(
        out3, dpi=150, bbox_inches='tight',
        pil_kwargs={'quality': 92})
    plt.close('all')

    print(f'Saving {out4} …', file=sys.stderr)
    make_fig4(cnt_hot, cnt_cold).savefig(
        out4, dpi=150, bbox_inches='tight',
        pil_kwargs={'quality': 92})
    plt.close('all')

    print('Done.', file=sys.stderr)


if __name__ == '__main__':
    main()

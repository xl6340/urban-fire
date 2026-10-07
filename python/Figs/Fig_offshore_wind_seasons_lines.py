#!/usr/bin/env python3
"""
Offshore wind + fire seasonality: SoCal and NorCal on shared panels.

Three scatter–line panels (x = meteorological season):
  (a) Offshore wind days (fire-independent climatology)
      - SoCal: R1D-SAWRI > 1 (1990–2018)
      - NorCal: Diablo-like proxy from gridMET (NE wind + vs≥4 m/s over
        ≥10% of NorCal cells; 1990–2018; no RH — background wind days)
  (b) Fire counts (WUI / Wildland), all sizes
  (c) Fraction of fires ignited on offshore-wind days (no p90 split)

SoCal and NorCal are plotted together (different markers/linestyles),
not summed into one series.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

BASE = Path(__file__).resolve().parents[2]
SOCAL_FLAGS = BASE / "dataPrc/CalFire_SoCal_SAW_extended_flags.csv"
NORCAL_FLAGS = BASE / "dataPrc/CalFire_NorCal_Diablo_flags.csv"
SAWRI_PATH = BASE / "dataSrc/climate/SAWRI/R1DSAWRI_010148_123118_red.txt"
VS_DIR = BASE / "dataSrc/climate/gridMET/wind speed"
TH_DIR = BASE / "dataSrc/climate/gridMET/wind direction"
WINDDAY_CACHE = BASE / "dataPrc/offshore_wind_days_by_season.csv"
OUT_PNG = BASE / "Fig/Fig_offshore_wind_seasons_lines.png"
OUT_PDF = BASE / "Fig/Fig_offshore_wind_seasons_lines.pdf"
OUT_CSV = BASE / "dataPrc/CalFire_offshore_wind_seasons_lines.csv"

C_WUI = np.array([216, 118, 89]) / 255
C_WILD = np.array([41, 157, 143]) / 255
C_SOCAL = np.array([180, 90, 60]) / 255
C_NORCAL = np.array([40, 110, 130]) / 255

SEASON_MAP = {
    12: "Winter",
    1: "Winter",
    2: "Winter",
    3: "Spring",
    4: "Spring",
    5: "Spring",
    6: "Summer",
    7: "Summer",
    8: "Summer",
    9: "Autumn",
    10: "Autumn",
    11: "Autumn",
}
SEASONS = ["Spring", "Summer", "Autumn", "Winter"]
SEASON_LAB = {
    "Spring": "Spring\n(MAM)",
    "Summer": "Summer\n(JJA)",
    "Autumn": "Autumn\n(SON)",
    "Winter": "Winter\n(DJF)",
}
YEAR0, YEAR1 = 1990, 2018  # overlap with SAWRI for panel (a)


def load_sawri(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df.columns = ["dt", "SAWRI"]
    parsed = pd.to_datetime(df["dt"], format="%m/%d/%y")
    yy = parsed.dt.year % 100
    year = np.where(yy >= 48, 1900 + yy, 2000 + yy)
    dt = pd.to_datetime({"year": year, "month": parsed.dt.month, "day": parsed.dt.day})
    s = pd.Series(df["SAWRI"].astype(float).to_numpy(), index=dt).sort_index()
    return s[~s.index.duplicated(keep="first")]


def socal_saw_days() -> pd.Series:
    sawri = load_sawri(SAWRI_PATH)
    d = sawri[(sawri.index.year >= YEAR0) & (sawri.index.year <= YEAR1) & (sawri > 1)]
    season = d.index.month.map(SEASON_MAP)
    n_years = YEAR1 - YEAR0 + 1
    counts = season.value_counts()
    return pd.Series({s: counts.get(s, 0) / n_years for s in SEASONS}, name="SoCal")


def norcal_diablo_days(force: bool = False) -> pd.Series:
    """Mean Diablo-like days/year by season (gridMET areal proxy, 1990–2018)."""
    if WINDDAY_CACHE.exists() and not force:
        cache = pd.read_csv(WINDDAY_CACHE)
        if "NorCal_mean_days" in cache.columns:
            return cache.set_index("season").loc[SEASONS, "NorCal_mean_days"]

    rows = []
    for year in range(YEAR0, YEAR1 + 1):
        vs_p = VS_DIR / f"vs_{year}.nc"
        th_p = TH_DIR / f"th_{year}.nc"
        if not vs_p.exists() or not th_p.exists():
            print(f"  skip {year}: missing gridMET")
            continue
        vs = xr.open_dataset(vs_p)["wind_speed"].sel(
            lat=slice(42, 36.5), lon=slice(-124.5, -120)
        )
        th = xr.open_dataset(th_p)["wind_from_direction"].sel(
            lat=slice(42, 36.5), lon=slice(-124.5, -120)
        )
        vs = vs.coarsen(lat=4, lon=4, boundary="trim").mean().load()
        th = th.coarsen(lat=4, lon=4, boundary="trim").mean().load()
        ne = (th >= 315) | (th <= 90)
        frac = ((ne) & (vs >= 4.0)).mean(dim=("lat", "lon"))
        is_day = (frac >= 0.10).to_pandas()
        is_day.index = pd.to_datetime(is_day.index)
        seas = is_day.index.month.map(SEASON_MAP)
        for s in SEASONS:
            rows.append({"year": year, "season": s, "n_days": int(is_day[seas == s].sum())})
        print(f"  NorCal wind days {year}: {int(is_day.sum())}")

    df = pd.DataFrame(rows)
    mean_days = df.groupby("season")["n_days"].mean().reindex(SEASONS)
    out = pd.DataFrame(
        {
            "season": SEASONS,
            "SoCal_mean_days": socal_saw_days().reindex(SEASONS).to_numpy(),
            "NorCal_mean_days": mean_days.to_numpy(),
        }
    )
    out.to_csv(WINDDAY_CACHE, index=False)
    print("Saved", WINDDAY_CACHE)
    return mean_days


def fire_season_stats(path: Path, keep_col: str, wind_col: str, region: str) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["IDate"])
    df = df[df[keep_col]].copy()
    df = df[df["IDate"].notna() & df["FireType"].isin(["WUI", "Wildland"])].copy()
    df["season"] = df["IDate"].dt.month.map(SEASON_MAP)
    df["wind"] = df[wind_col].astype(bool)
    rows = []
    for season in SEASONS:
        for ft in ["WUI", "Wildland"]:
            g = df[(df["season"] == season) & (df["FireType"] == ft)]
            n = len(g)
            n_wind = int(g["wind"].sum()) if n else 0
            rows.append(
                {
                    "region": region,
                    "season": season,
                    "FireType": ft,
                    "n_fires": n,
                    "n_wind_fires": n_wind,
                    "frac_wind": (n_wind / n) if n else np.nan,
                }
            )
    return pd.DataFrame(rows)


def _style(ax):
    ax.tick_params(direction="out", labelsize=8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_xticks(range(len(SEASONS)))
    ax.set_xticklabels([SEASON_LAB[s] for s in SEASONS], fontsize=8)


def main():
    print("Panel (a): offshore wind days …")
    socal_days = socal_saw_days()
    print("  computing/loading NorCal Diablo-like days …")
    norcal_days = norcal_diablo_days()
    # refresh cache SoCal column if needed
    if WINDDAY_CACHE.exists():
        cache = pd.read_csv(WINDDAY_CACHE)
        cache["SoCal_mean_days"] = socal_days.reindex(SEASONS).to_numpy()
        cache["NorCal_mean_days"] = norcal_days.reindex(SEASONS).to_numpy()
        cache.to_csv(WINDDAY_CACHE, index=False)

    socal = fire_season_stats(SOCAL_FLAGS, "in_analysis", "SAW", "Southern California")
    norcal = fire_season_stats(NORCAL_FLAGS, "has_wx", "is_Diablo", "Northern California")
    fires = pd.concat([socal, norcal], ignore_index=True)
    fires.to_csv(OUT_CSV, index=False)
    print("Saved", OUT_CSV)
    print(fires.to_string(index=False))

    x = np.arange(len(SEASONS))
    fig, axes = plt.subplots(3, 1, figsize=(6.2, 7.2), sharex=True, facecolor="w")
    fig.subplots_adjust(hspace=0.22, left=0.14, right=0.98, top=0.92, bottom=0.08)

    # ── (a) wind days ─────────────────────────────────────────────────────
    ax = axes[0]
    ax.plot(
        x,
        [socal_days[s] for s in SEASONS],
        "-o",
        color=C_SOCAL,
        lw=1.8,
        ms=7,
        label="SoCal Santa Ana (SAWRI>1)",
    )
    ax.plot(
        x,
        [norcal_days[s] for s in SEASONS],
        "--s",
        color=C_NORCAL,
        lw=1.8,
        ms=7,
        label="NorCal Diablo-like (gridMET)",
    )
    ax.set_ylabel("Mean offshore wind\ndays per year", fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_ylim(bottom=0)
    _style(ax)
    ax.text(0.01, 0.96, "a", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    # ── (b) fire counts ───────────────────────────────────────────────────
    ax = axes[1]
    for region, ls, marker in [
        ("Southern California", "-", "o"),
        ("Northern California", "--", "s"),
    ]:
        for ft, color in [("WUI", C_WUI), ("Wildland", C_WILD)]:
            sub = fires[(fires.region == region) & (fires.FireType == ft)].set_index("season")
            y = [sub.loc[s, "n_fires"] for s in SEASONS]
            ax.plot(
                x,
                y,
                ls=ls,
                marker=marker,
                color=color,
                lw=1.6,
                ms=6,
                label=f"{'SoCal' if region.startswith('S') else 'NorCal'} {ft}",
            )
    ax.set_ylabel("Number of fires", fontsize=9)
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="upper right")
    ax.set_ylim(bottom=0)
    _style(ax)
    ax.text(0.01, 0.96, "b", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    # ── (c) fraction on offshore days ─────────────────────────────────────
    ax = axes[2]
    for region, ls, marker in [
        ("Southern California", "-", "o"),
        ("Northern California", "--", "s"),
    ]:
        for ft, color in [("WUI", C_WUI), ("Wildland", C_WILD)]:
            sub = fires[(fires.region == region) & (fires.FireType == ft)].set_index("season")
            y = [100 * sub.loc[s, "frac_wind"] for s in SEASONS]
            ax.plot(
                x,
                y,
                ls=ls,
                marker=marker,
                color=color,
                lw=1.6,
                ms=6,
                label=f"{'SoCal' if region.startswith('S') else 'NorCal'} {ft}",
            )
    ax.set_ylabel("Fires on offshore\nwind days (%)", fontsize=9)
    ax.legend(frameon=False, fontsize=7.5, ncol=2, loc="upper left")
    ax.set_ylim(0, 100)
    _style(ax)
    ax.text(0.01, 0.96, "c", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="w")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="w")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)


if __name__ == "__main__":
    main()

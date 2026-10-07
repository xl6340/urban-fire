#!/usr/bin/env python3
"""
Monthly climatology:

  (A) Statewide CalFire ignition counts, all sizes, 1990–2025
  (B) Statewide burned area, all sizes, 1990–2025
  (C) NLDN CG lightning (1999–2025) + Santa Ana (SAWRI>1, 1990–2018)
      (CalFire_fullSize.shp; Unknown coded as Human)
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, MultipleLocator, NullLocator
import numpy as np
import pandas as pd
import xarray as xr

BASE = Path(__file__).resolve().parents[2]
import sys

sys.path.insert(0, str(BASE / "code"))
from diablo_raws_catalog import PROXY as DIABLO_PROXY
from diablo_raws_catalog import YEAR0 as DIABLO_YEAR0
from diablo_raws_catalog import monthly_mean_days

SAWRI_PATH = BASE / "dataSrc/climate/SAWRI/R1DSAWRI_010148_123118_red.txt"
VS_DIR = BASE / "dataSrc/climate/gridMET/wind speed"
TH_DIR = BASE / "dataSrc/climate/gridMET/wind direction"
WINDDAY_CACHE = BASE / "dataPrc/offshore_wind_days_by_month.csv"
LTNG_DIR = Path(__file__).resolve().parents[2] / "dataPrc" / "lightning" / "california"  # place NLDN monthly files here, or set LTNG_DIR
CA_SHP = BASE / "dataSrc/boundary/ca_state/CA_State.shp"
LTNG_CACHE = BASE / "dataPrc/NLDN_CG_density_monthly.csv"
OUT_PNG = BASE / "Fig/Fig_offshore_wind_monthly.png"
OUT_PDF = BASE / "Fig/Fig_offshore_wind_monthly.pdf"
OUT_CSV = BASE / "dataPrc/CalFire_offshore_wind_monthly.csv"
OUT_WUI_PNG = BASE / "Fig/Fig_offshore_wind_monthly_WUI.png"
OUT_WUI_PDF = BASE / "Fig/Fig_offshore_wind_monthly_WUI.pdf"
OUT_WUI_CSV = BASE / "dataPrc/CalFire_offshore_wind_monthly_WUI.csv"

C_HUMAN = np.array([253, 216, 93]) / 255   # #fdd85d (Fig. 4 Human)
C_NAT = np.array([153, 214, 234]) / 255    # #99d6ea (Fig. 4 Natural)
C_FRAC = np.array([90, 90, 110]) / 255
C_LTNG = np.array([50, 130, 160]) / 255    # deeper cyan (NLDN; darker than Natural)
C_SAW = np.array([220, 160, 30]) / 255      # deeper gold yellow (Santa Ana)
C_DIABLO = np.array([40, 110, 130]) / 255

MONTHS = np.arange(1, 13)
MONTH_LAB = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
             "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
YEAR0, YEAR1 = 1990, 2025
N_YEARS = YEAR1 - YEAR0 + 1
LTNG_YEAR0, LTNG_YEAR1 = 1999, 2025
SAW_YEAR0, SAW_YEAR1 = 1990, 2018  # R1D-SAWRI ends 31 Dec 2018
N_SAW_YEARS = SAW_YEAR1 - SAW_YEAR0 + 1


def load_sawri(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df.columns = ["dt", "SAWRI"]
    parsed = pd.to_datetime(df["dt"], format="%m/%d/%y")
    yy = parsed.dt.year % 100
    year = np.where(yy >= 48, 1900 + yy, 2000 + yy)
    dt = pd.to_datetime({"year": year, "month": parsed.dt.month, "day": parsed.dt.day})
    s = pd.Series(df["SAWRI"].astype(float).to_numpy(), index=dt).sort_index()
    return s[~s.index.duplicated(keep="first")]


def socal_saw_days_monthly() -> pd.DataFrame:
    """Mean ± SD Santa Ana days/year by month (across years 1990–2018)."""
    sawri = load_sawri(SAWRI_PATH)
    sawri = sawri[(sawri.index.year >= SAW_YEAR0) & (sawri.index.year <= SAW_YEAR1)]
    years = np.arange(SAW_YEAR0, SAW_YEAR1 + 1)
    # year × month count of SAWRI > 1 days
    mat = np.zeros((len(years), 12), dtype=float)
    for i, y in enumerate(years):
        sub = sawri[(sawri.index.year == y) & (sawri > 1)]
        for m, c in sub.index.month.value_counts().items():
            mat[i, int(m) - 1] = float(c)
    return pd.DataFrame(
        {
            "month": MONTHS,
            "mean": mat.mean(axis=0),
            "std": mat.std(axis=0, ddof=1),
        }
    ).set_index("month")


def norcal_diablo_days_monthly(force: bool = False) -> pd.Series:
    """Mean Diablo days/year by month from hourly RAWS (2011–2018).

    Smith et al. (2018) NoBA: any of 6 RAWS, wind 315–135°, vs > 11.17 m/s,
    RH < 30%, ≥3 consecutive hours. IEM HADS archive starts in 2011.
    """
    if WINDDAY_CACHE.exists() and not force:
        cache = pd.read_csv(WINDDAY_CACHE)
        if (
            "NorCal_mean_days" in cache.columns
            and "diablo_proxy" in cache.columns
            and len(cache) == 12
            and cache["diablo_proxy"].astype(str).eq(DIABLO_PROXY).all()
        ):
            return cache.set_index("month")["NorCal_mean_days"]
    return monthly_mean_days(DIABLO_YEAR0, YEAR1)


def fire_counts_monthly(firetype: str | None = None) -> pd.DataFrame:
    """
    Statewide monthly counts and burned area from the full-size CalFire inventory.

    Means are across years (1990–2025). ``std_*`` / ``std_area_*`` are the
    year-to-year standard deviation of that month's annual total (ddof=1).
    """
    g = gpd.read_file(BASE / "dataPrc/firePrmt/CalFire_fullSize.shp", ignore_geometry=True)
    g["IDate"] = pd.to_datetime(g["IDate"], errors="coerce")
    g = g[
        g["IDate"].notna()
        & (g["year"] >= YEAR0)
        & (g["year"] <= YEAR1)
        & g["Ignition"].isin(["Human", "Natural"])
    ].copy()
    if firetype is not None:
        g = g[g["FireType"] == firetype].copy()
    g["month"] = g["IDate"].dt.month
    g["size"] = pd.to_numeric(g["size"], errors="coerce").fillna(0.0)

    years = np.arange(YEAR0, YEAR1 + 1)
    # year × month totals
    n_h = np.zeros((len(years), 12), dtype=float)
    n_n = np.zeros((len(years), 12), dtype=float)
    a_h = np.zeros((len(years), 12), dtype=float)
    a_n = np.zeros((len(years), 12), dtype=float)
    y_index = {int(y): i for i, y in enumerate(years)}
    for (y, m, ign), sub in g.groupby(["year", "month", "Ignition"]):
        if int(y) not in y_index:
            continue
        i, j = y_index[int(y)], int(m) - 1
        if ign == "Human":
            n_h[i, j] = len(sub)
            a_h[i, j] = float(sub["size"].sum())
        elif ign == "Natural":
            n_n[i, j] = len(sub)
            a_n[i, j] = float(sub["size"].sum())

    rows = []
    for j, m in enumerate(MONTHS):
        nh_tot = float(n_h[:, j].sum())
        nn_tot = float(n_n[:, j].sum())
        ah_tot = float(a_h[:, j].sum())
        an_tot = float(a_n[:, j].sum())
        n_known = nh_tot + nn_tot
        a_known = ah_tot + an_tot
        # annual human fraction for this month (years with any fires / any area)
        n_known_y = n_h[:, j] + n_n[:, j]
        frac_n_y = np.divide(
            n_h[:, j],
            n_known_y,
            out=np.full(len(years), np.nan),
            where=n_known_y > 0,
        )
        a_known_y = a_h[:, j] + a_n[:, j]
        frac_a_y = np.divide(
            a_h[:, j],
            a_known_y,
            out=np.full(len(years), np.nan),
            where=a_known_y > 0,
        )
        rows.append(
            {
                "month": int(m),
                "n_human": int(nh_tot),
                "n_natural": int(nn_tot),
                "mean_human": float(n_h[:, j].mean()),
                "mean_natural": float(n_n[:, j].mean()),
                "std_human": float(n_h[:, j].std(ddof=1)),
                "std_natural": float(n_n[:, j].std(ddof=1)),
                "frac_human": (nh_tot / n_known) if n_known else np.nan,
                "std_frac_human": float(np.nanstd(frac_n_y, ddof=1)),
                "area_human": ah_tot,
                "area_natural": an_tot,
                "mean_area_human": float(a_h[:, j].mean()),
                "mean_area_natural": float(a_n[:, j].mean()),
                "std_area_human": float(a_h[:, j].std(ddof=1)),
                "std_area_natural": float(a_n[:, j].std(ddof=1)),
                "frac_area_human": (ah_tot / a_known) if a_known else np.nan,
                "std_frac_area_human": float(np.nanstd(frac_a_y, ddof=1)),
            }
        )
    return pd.DataFrame(rows)


def nldn_density_monthly(force: bool = False) -> pd.DataFrame:
    """Mean ± SD NLDN CG flash density (flashes km⁻²) by month over CA land."""
    if LTNG_CACHE.exists() and not force:
        cache = pd.read_csv(LTNG_CACHE)
        if (
            "density" in cache.columns
            and "std_density" in cache.columns
            and len(cache) == 12
            and int(cache["year0"].iloc[0]) == LTNG_YEAR0
            and int(cache["year1"].iloc[0]) == LTNG_YEAR1
        ):
            return cache.set_index("month")

    ca = gpd.read_file(CA_SHP).to_crs(4326)
    ca_area_km2 = float(ca["ALAND"].iloc[0]) / 1e6
    ca_geom = ca.union_all() if hasattr(ca, "union_all") else ca.geometry.unary_union

    year_month = []  # list of (12,) density arrays
    for year in range(LTNG_YEAR0, LTNG_YEAR1 + 1):
        p = LTNG_DIR / f"lightning_CA_{year}.csv"
        if not p.exists():
            print(f"  skip lightning {year}: missing file")
            continue
        df = pd.read_csv(p, usecols=["ZDAY", "CENTERLON", "CENTERLAT", "TOTAL_COUNT"])
        print(f"  NLDN {year}: {len(df):,} grid-days")
        cells = df[["CENTERLON", "CENTERLAT"]].drop_duplicates()
        pts = gpd.GeoDataFrame(
            cells,
            geometry=gpd.points_from_xy(cells["CENTERLON"], cells["CENTERLAT"]),
            crs=4326,
        )
        pts["in_ca"] = pts.geometry.within(ca_geom)
        inside = pts.loc[pts["in_ca"], ["CENTERLON", "CENTERLAT"]]
        lit = df.merge(inside, on=["CENTERLON", "CENTERLAT"], how="inner")
        lit["month"] = pd.to_datetime(lit["ZDAY"].astype(str), format="%Y%m%d").dt.month
        totals = lit.groupby("month")["TOTAL_COUNT"].sum().reindex(MONTHS, fill_value=0)
        year_month.append((totals / ca_area_km2).to_numpy(dtype=float))

    mat = np.vstack(year_month)  # n_years × 12
    n_years = mat.shape[0]
    mean = mat.mean(axis=0)
    std = mat.std(axis=0, ddof=1)
    out = pd.DataFrame(
        {
            "month": MONTHS,
            "mean_flashes": (mean * ca_area_km2),
            "density": mean,
            "std_density": std,
            "ca_area_km2": ca_area_km2,
            "n_years": n_years,
            "year0": LTNG_YEAR0,
            "year1": LTNG_YEAR1,
        }
    )
    out.to_csv(LTNG_CACHE, index=False)
    print("Saved", LTNG_CACHE)
    return out.set_index("month")


def _style(ax, xlabel=True):
    ax.tick_params(direction="out", labelsize=8, labelbottom=xlabel)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.set_xticks(MONTHS)
    ax.set_xticklabels(MONTH_LAB if xlabel else [], fontsize=8)
    ax.set_xlim(0.5, 12.5)


def _panel_label(ax, letter):
    """(A)/(B) just outside, above the top-left corner of the axes frame."""
    ax.annotate(
        letter,
        xy=(0, 1),
        xycoords="axes fraction",
        xytext=(0, 10),
        textcoords="offset points",
        fontsize=10,
        ha="left",
        va="bottom",
        annotation_clip=False,
        clip_on=False,
    )


def _shade_between_months(ax, y_a, y_b, months=MONTHS, skip=(7, 8, 9), half=0.28):
    """Shade between two series; join contiguous shaded months (no cross-gap diagonals).
    Default: no shade for Jul–Sep.
    """
    color = 0.55 * C_HUMAN + 0.45 * C_NAT
    months = np.asarray(months, dtype=float)
    y_a = np.asarray(y_a, dtype=float)
    y_b = np.asarray(y_b, dtype=float)
    shade = ~np.isin(months.astype(int), list(skip))

    # Contiguous True runs
    n = len(months)
    i = 0
    while i < n:
        if not shade[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and shade[j + 1]:
            j += 1
        if j > i:
            # multi-month: continuous band along the two curves
            ax.fill_between(
                months[i : j + 1], y_a[i : j + 1], y_b[i : j + 1],
                color=color, alpha=0.35, lw=0, zorder=1,
            )
        else:
            # single month: thin vertical band (same-month comparison)
            m = months[i]
            lo, hi = sorted((y_a[i], y_b[i]))
            if np.isfinite(lo) and np.isfinite(hi) and lo != hi:
                ax.fill_between(
                    [m - half, m + half], [lo, lo], [hi, hi],
                    color=color, alpha=0.35, lw=0, zorder=1,
                )
        i = j + 1


def _fire_panel(
    ax,
    fires,
    letter,
    ylabel,
    y_human,
    y_natural,
    frac,
    frac_ylabel,
    xlabel=True,
    frac_ylim=(60, 105),
    legend_anchor=(0.02, 0.68),
    frac_legend="Human ignition proportion",
    show_legend=True,
    show_frac=True,
    shade_off_season=False,
    ylim=None,
    ylog=False,
    labelleft=True,
    labelright=True,
    std_human=None,
    std_natural=None,
    std_frac=None,
):
    y_h = fires[y_human].to_numpy(dtype=float)
    y_n = fires[y_natural].to_numpy(dtype=float)

    if (not ylog) and std_human is not None and std_human in fires.columns:
        s_h = fires[std_human].to_numpy(dtype=float)
        ax.fill_between(
            MONTHS, np.maximum(y_h - s_h, 0), y_h + s_h,
            color=C_HUMAN, alpha=0.22, lw=0, zorder=1,
        )
    if (not ylog) and std_natural is not None and std_natural in fires.columns:
        s_n = fires[std_natural].to_numpy(dtype=float)
        ax.fill_between(
            MONTHS, np.maximum(y_n - s_n, 0), y_n + s_n,
            color=C_NAT, alpha=0.22, lw=0, zorder=1,
        )

    if ylog:
        # true log: zeros drawn at axis floor (hollow markers); keep x aligned
        pos = np.concatenate([y_h[y_h > 0], y_n[y_n > 0]])
        # axis limits from panel's own positive values (ignore numerical dust <1e-2)
        pos_lim = pos[pos >= 1e-2] if np.any(pos >= 1e-2) else pos
        y0 = float(pos_lim.min()) * 0.5
        y1 = float(pos_lim.max()) * 2.0
        if ylim is not None:
            y0, y1 = float(ylim[0]), float(ylim[1])
        eps = y0
        # Connect through floor so zero months stay on the correct tick
        y_h_line = np.maximum(y_h, eps)
        y_n_line = np.maximum(y_n, eps)
        if shade_off_season:
            # Same-month band between human & natural; no shade for Jul–Sep
            _shade_between_months(ax, y_n_line, y_h_line)
            # Natural ignited to axis floor for Jun–Oct
            ax.fill_between(
                MONTHS[5:10], y_n_line[5:10], eps,
                color=C_NAT, alpha=0.28, lw=0, zorder=1,
            )
        (h1,) = ax.plot(
            MONTHS, y_h_line, "-", color=C_HUMAN, lw=1.8,
            label="Human ignited", zorder=3,
        )
        (h2,) = ax.plot(
            MONTHS, y_n_line, "-", color=C_NAT, lw=1.8,
            label="Natural ignited", zorder=3,
        )
        # Filled markers for positive months; hollow at floor for zeros
        pos_h = y_h >= eps
        pos_n = y_n >= eps
        ax.plot(
            MONTHS[pos_h], y_h[pos_h], "o", color=C_HUMAN, ms=6.5,
            zorder=4,
        )
        ax.plot(
            MONTHS[pos_n], y_n[pos_n], "s", color=C_NAT, ms=6.5,
            zorder=4,
        )
        if np.any(~pos_h):
            ax.plot(
                MONTHS[~pos_h], np.full((~pos_h).sum(), y0), "o", ms=6.5,
                mfc="white", mec=C_HUMAN, mew=1.2, zorder=4,
            )
        if np.any(~pos_n):
            ax.plot(
                MONTHS[~pos_n], np.full((~pos_n).sum(), y0), "s", ms=6.5,
                mfc="white", mec=C_NAT, mew=1.2, zorder=4,
            )
        ax.set_yscale("log")
        ax.set_ylim(y0, y1)
        ax.yaxis.set_major_locator(LogLocator(base=10))
        ax.yaxis.set_minor_locator(NullLocator())
    else:
        if shade_off_season:
            _shade_between_months(ax, y_n, y_h)
        (h1,) = ax.plot(
            MONTHS, y_h, "-o", color=C_HUMAN, lw=1.8, ms=6.5, label="Human ignited", zorder=3
        )
        (h2,) = ax.plot(
            MONTHS, y_n, "-s", color=C_NAT, lw=1.8, ms=6.5, label="Natural ignited", zorder=3
        )
        ax.set_ylim(bottom=0)
        if ylim is not None:
            ax.set_ylim(*ylim)

    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9)
    _style(ax, xlabel=xlabel)
    ax.tick_params(axis="y", labelleft=labelleft)

    handles = [h1, h2]
    if show_frac:
        ax.spines["right"].set_visible(True)
        ax2 = ax.twinx()
        frac_pct = 100 * fires[frac].to_numpy(dtype=float)
        if std_frac is not None and std_frac in fires.columns:
            s_f = 100 * fires[std_frac].to_numpy(dtype=float)
            ax2.fill_between(
                MONTHS,
                np.clip(frac_pct - s_f, 0, 100),
                np.clip(frac_pct + s_f, 0, 100),
                color=C_FRAC, alpha=0.15, lw=0, zorder=1,
            )
        (h3,) = ax2.plot(
            MONTHS,
            frac_pct,
            ":^",
            color=C_FRAC,
            lw=1.6,
            ms=5.5,
            label=frac_legend,
            zorder=3,
        )
        if frac_ylabel:
            ax2.set_ylabel(frac_ylabel, fontsize=9, color=C_FRAC)
        ax2.set_ylim(*frac_ylim)
        ax2.yaxis.set_major_locator(MultipleLocator(20))
        ax2.yaxis.set_minor_locator(NullLocator())
        ax2.tick_params(axis="y", labelsize=8, colors=C_FRAC, labelright=labelright)
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_color(C_FRAC)
        handles.append(h3)
    else:
        ax.spines["right"].set_visible(False)

    if show_legend:
        ax.legend(
            handles=handles,
            frameon=False,
            fontsize=8,
            loc="upper left",
            bbox_to_anchor=legend_anchor,
        )
    _panel_label(ax, letter)
    return (None if not show_frac else ax2), tuple(handles)


def plot_wui_wildland_monthly():
    """WUI vs wildland burned area (log), plus lightning / Santa Ana (stacked A–C)."""

    wui = fire_counts_monthly("Urban-edge")
    wild = fire_counts_monthly("Wildland")
    wui = wui.assign(FireType="WUI")
    wild = wild.assign(FireType="Wildland")
    pd.concat([wui, wild], ignore_index=True).to_csv(OUT_WUI_CSV, index=False)

    print("WUI\n", wui.to_string(index=False))
    print("Wildland\n", wild.to_string(index=False))

    print("Santa Ana days (SAWRI>1, 1990–2018) …")
    saw = socal_saw_days_monthly()
    print("NLDN CG density (1999–2025) …")
    ltng = nldn_density_monthly()

    fig, axes = plt.subplots(3, 1, figsize=(5.6, 6.6), sharex=True, facecolor="w")
    fig.subplots_adjust(hspace=0.36, left=0.15, right=0.86, top=0.96, bottom=0.08)

    _fire_panel(
        axes[0],
        wui,
        "(A) WUI fire",
        "Burned area\n(km$^2$ yr$^{-1}$)",
        "mean_area_human",
        "mean_area_natural",
        "frac_area_human",
        "Human proportion (%)",
        xlabel=True,
        show_legend=True,
        show_frac=False,
        ylog=True,
        shade_off_season=True,
        legend_anchor=(0.02, 0.98),
    )
    _fire_panel(
        axes[1],
        wild,
        "(B) Wildland fire",
        "Burned area\n(km$^2$ yr$^{-1}$)",
        "mean_area_human",
        "mean_area_natural",
        "frac_area_human",
        "Human proportion (%)",
        xlabel=True,
        show_legend=False,
        show_frac=False,
        ylog=True,
        shade_off_season=True,
    )

    ax = axes[2]
    ltng_mean = ltng["density"].to_numpy(dtype=float)
    # Shade lightning line to x-axis for Jun–Oct
    ax.fill_between(
        MONTHS[5:10], ltng_mean[5:10], 0,
        color=C_LTNG, alpha=0.22, lw=0, zorder=1,
    )
    (h_ltng,) = ax.plot(
        MONTHS,
        ltng_mean,
        "-D",
        color=C_LTNG,
        lw=1.8,
        ms=6,
        label="NLDN CG density",
        zorder=3,
    )
    ax.set_ylabel("Lightning density\n(flashes km$^{-2}$)", fontsize=9, color="k")
    ax.set_ylim(0, ltng_mean.max() * 1.08)
    ax.yaxis.set_major_locator(MultipleLocator(0.02))
    ax.yaxis.set_minor_locator(NullLocator())
    _style(ax, xlabel=True)
    ax.tick_params(axis="y", colors="k", labelsize=8)
    ax.spines["left"].set_color("k")
    ax.spines["right"].set_visible(True)

    ax_saw = ax.twinx()
    saw_mean = saw["mean"].to_numpy(dtype=float)
    # Shade Santa Ana line to x-axis for Jan–Apr and Sep–Dec
    for i0, i1 in ((0, 4), (8, 12)):  # Jan–Apr, Sep–Dec
        ax_saw.fill_between(
            MONTHS[i0:i1], saw_mean[i0:i1], 0,
            color=C_SAW, alpha=0.22, lw=0, zorder=1,
        )
    (h_saw,) = ax_saw.plot(
        MONTHS,
        saw_mean,
        "-o",
        color=C_SAW,
        lw=1.8,
        ms=6.5,
        label="Santa Ana (SAWRI > 1)",
        zorder=3,
    )
    ax_saw.set_ylabel("Santa Ana days / year", fontsize=9, color="k")
    ax_saw.set_ylim(0, saw_mean.max() * 1.08)
    ax_saw.yaxis.set_major_locator(MultipleLocator(5))
    ax_saw.yaxis.set_minor_locator(NullLocator())
    ax_saw.tick_params(axis="y", labelsize=8, colors="k")
    ax_saw.spines["top"].set_visible(False)
    ax_saw.spines["right"].set_color("k")
    ax.legend(
        handles=[h_ltng, h_saw],
        frameon=False,
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(0.14, 0.86),
    )
    _panel_label(ax, "(C) Lightning-Santa Ana")

    for pth in (
        OUT_WUI_PNG,
        OUT_WUI_PDF,
        OUT_PNG,  # also overwrite the path the user has open
        OUT_PDF,
        BASE / "Fig/FigS9.png",
        BASE / "Fig/FigS9.pdf",
    ):
        fig.savefig(pth, dpi=300, bbox_inches="tight", facecolor="w")
    print("Saved", OUT_WUI_PNG)
    print("Saved", OUT_PNG)
    plt.close(fig)


def main():
    print("Fires (statewide, all sizes, 1990–2025) …")
    fires = fire_counts_monthly()

    print("Santa Ana days (SAWRI>1, 1990–2018) …")
    saw = socal_saw_days_monthly()

    print("Diablo days (Smith 2018 NoBA RAWS, 2011–2018) …")
    diablo = norcal_diablo_days_monthly()

    print("NLDN CG density (1999–2025) …")
    ltng = nldn_density_monthly()

    out = fires.copy()
    out["Santa_Ana_mean_days"] = saw["mean"].reindex(MONTHS).to_numpy()
    out["Santa_Ana_std_days"] = saw["std"].reindex(MONTHS).to_numpy()
    out["Diablo_mean_days"] = [diablo.loc[m] for m in MONTHS]
    out["NLDN_CG_density"] = ltng["density"].reindex(MONTHS).to_numpy()
    out["NLDN_CG_std_density"] = ltng["std_density"].reindex(MONTHS).to_numpy()
    out.to_csv(OUT_CSV, index=False)
    wind_cache = pd.DataFrame(
        {
            "month": MONTHS,
            "SoCal_mean_days": saw["mean"].reindex(MONTHS).to_numpy(),
            "SoCal_std_days": saw["std"].reindex(MONTHS).to_numpy(),
            "NorCal_mean_days": [diablo.loc[m] for m in MONTHS],
            "diablo_proxy": DIABLO_PROXY,
        }
    )
    wind_cache.to_csv(WINDDAY_CACHE, index=False)
    print("Saved", OUT_CSV)
    print(out.to_string(index=False))

    fig, axes = plt.subplots(3, 1, figsize=(6.6, 8.2), sharex=True, facecolor="w")
    fig.subplots_adjust(hspace=0.38, left=0.14, right=0.88, top=0.93, bottom=0.07)

    _fire_panel(
        axes[0],
        fires,
        "(A)",
        "Mean fires per year",
        "mean_human",
        "mean_natural",
        "frac_human",
        "Human ignition proportion (%)",
        xlabel=False,
        std_human="std_human",
        std_natural="std_natural",
    )
    _fire_panel(
        axes[1],
        fires,
        "(B)",
        "Mean burned area\n(km$^2$ yr$^{-1}$)",
        "mean_area_human",
        "mean_area_natural",
        "frac_area_human",
        "Human burned-area proportion (%)",
        xlabel=False,
        frac_ylim=(30, 105),
        legend_anchor=(0.02, 0.68),
        frac_legend="Human burned-area proportion",
        std_human="std_area_human",
        std_natural="std_area_natural",
        std_frac="std_frac_area_human",
    )

    ax = axes[2]
    ltng_mean = ltng["density"].reindex(MONTHS).to_numpy(dtype=float)
    ltng_std = ltng["std_density"].reindex(MONTHS).to_numpy(dtype=float)
    ax.fill_between(
        MONTHS, np.maximum(ltng_mean - ltng_std, 0), ltng_mean + ltng_std,
        color=C_LTNG, alpha=0.22, lw=0, zorder=1,
    )
    (h_ltng,) = ax.plot(
        MONTHS, ltng_mean, "-D", color=C_LTNG, lw=1.8, ms=6,
        label="NLDN CG density", zorder=3,
    )
    ax.set_ylabel("Lightning density\n(flashes km$^{-2}$)", fontsize=9, color=C_LTNG)
    ax.set_ylim(0, (ltng_mean + ltng_std).max() * 1.08)
    _style(ax, xlabel=True)
    ax.tick_params(axis="y", colors=C_LTNG)
    ax.spines["left"].set_color(C_LTNG)
    ax.spines["right"].set_visible(True)

    ax_saw = ax.twinx()
    saw_mean = saw["mean"].reindex(MONTHS).to_numpy(dtype=float)
    saw_std = saw["std"].reindex(MONTHS).to_numpy(dtype=float)
    ax_saw.fill_between(
        MONTHS, np.maximum(saw_mean - saw_std, 0), saw_mean + saw_std,
        color=C_SAW, alpha=0.22, lw=0, zorder=1,
    )
    (h_saw,) = ax_saw.plot(
        MONTHS, saw_mean, "-o", color=C_SAW, lw=1.8, ms=6.5,
        label="Santa Ana (SAWRI > 1)", zorder=3,
    )
    ax_saw.set_ylabel("Santa Ana days / year", fontsize=9, color=C_SAW)
    ax_saw.set_ylim(0, (saw_mean + saw_std).max() * 1.08)
    ax_saw.tick_params(axis="y", labelsize=8, colors=C_SAW)
    ax_saw.spines["top"].set_visible(False)
    ax_saw.spines["right"].set_color(C_SAW)
    ax.legend(
        handles=[h_ltng, h_saw],
        frameon=False,
        fontsize=8,
        loc="upper left",
        bbox_to_anchor=(0.14, 0.86),
    )
    _panel_label(ax, "(C)")

    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="w")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="w")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)


if __name__ == "__main__":
    main()
    plot_wui_wildland_monthly()

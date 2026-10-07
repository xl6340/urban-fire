#!/usr/bin/env python3
"""
Offshore-wind ignition fractions by meteorological season.

Two stacked panels (all four seasons in each):
  (a) Southern California — Santa Ana (hybrid SAW flag, in_analysis fires)
  (b) Northern California — Diablo-like (gridMET proxy, has_wx fires)

Solid = fire size < p90; hatched = ≥ p90.
p90 is the within-season, within-FireType 90th percentile of fire size.
Season order: Spring, Summer, Autumn, Winter (MAM–JJA–SON–DJF).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

BASE = Path(__file__).resolve().parents[2]
SOCAL_FLAGS = BASE / "dataPrc/CalFire_SoCal_SAW_extended_flags.csv"
NORCAL_FLAGS = BASE / "dataPrc/CalFire_NorCal_Diablo_flags.csv"
OUT_PNG = BASE / "Fig/Fig_offshore_wind_fractions_seasons.png"
OUT_PDF = BASE / "Fig/Fig_offshore_wind_fractions_seasons.pdf"
OUT_CSV = BASE / "dataPrc/CalFire_offshore_wind_fractions_seasons.csv"

C_WUI = np.array([216, 118, 89]) / 255
C_WILD = np.array([41, 157, 143]) / 255

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
SEASON_SUB = {
    "Spring": "MAM",
    "Summer": "JJA",
    "Autumn": "SON",
    "Winter": "DJF",
}


def _load(path, keep_col, wind_col):
    df = pd.read_csv(path, parse_dates=["IDate"])
    df = df[df[keep_col]].copy()
    df = df[df["IDate"].notna() & df["FireType"].isin(["WUI", "Wildland"])].copy()
    df["season"] = df["IDate"].dt.month.map(SEASON_MAP)
    df["wind"] = df[wind_col].astype(bool)
    return df


def season_summary(df):
    rows = []
    for season in SEASONS:
        for ft in ["WUI", "Wildland"]:
            g = df[(df["season"] == season) & (df["FireType"] == ft)]
            p90 = g["size"].quantile(0.9) if len(g) else np.nan
            large = g[g["size"] >= p90] if pd.notna(p90) else g.iloc[0:0]
            below = g[g["size"] < p90] if pd.notna(p90) else g.iloc[0:0]
            n_below = len(below)
            n_large = len(large)
            rows.append(
                {
                    "season": season,
                    "FireType": ft,
                    "n": len(g),
                    "p90_km2": p90,
                    "n_below": n_below,
                    "n_wind_below": int(below["wind"].sum()) if n_below else 0,
                    "frac_below": float(below["wind"].mean()) if n_below else np.nan,
                    "n_large": n_large,
                    "n_large_wind": int(large["wind"].sum()) if n_large else 0,
                    "frac_large": float(large["wind"].mean()) if n_large else np.nan,
                }
            )
    return pd.DataFrame(rows)


def _bar_label(ax, bar, n_s, n, pad):
    if not np.isfinite(n) or n <= 0 or not np.isfinite(bar.get_height()):
        return
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + pad,
        f"{bar.get_height():.1f}%\n({int(n_s)}/{int(n)})",
        ha="center",
        va="bottom",
        fontsize=6,
        color="0.2",
        linespacing=1.0,
    )


def _vals(sub, col):
    return sub.set_index("FireType").loc[["WUI", "Wildland"], col].to_numpy()


def _panel(ax, df, ylabel):
    w = 0.32
    pair = 0.84
    block = 2.35
    centers = np.arange(len(SEASONS)) * block
    x_wui = centers - pair / 2
    x_wild = centers + pair / 2

    below = np.zeros((len(SEASONS), 2))
    large = np.zeros((len(SEASONS), 2))
    n_below = np.zeros((len(SEASONS), 2))
    n_wind_below = np.zeros((len(SEASONS), 2))
    n_large = np.zeros((len(SEASONS), 2))
    n_large_wind = np.zeros((len(SEASONS), 2))
    ok_below = np.zeros((len(SEASONS), 2), dtype=bool)
    ok_large = np.zeros((len(SEASONS), 2), dtype=bool)

    for i, season in enumerate(SEASONS):
        sub = df[df["season"] == season]
        below[i] = _vals(sub, "frac_below")
        large[i] = _vals(sub, "frac_large")
        n_below[i] = _vals(sub, "n_below")
        n_wind_below[i] = _vals(sub, "n_wind_below")
        n_large[i] = _vals(sub, "n_large")
        n_large_wind[i] = _vals(sub, "n_large_wind")
        ok_below[i] = np.isfinite(below[i]) & (n_below[i] > 0)
        ok_large[i] = np.isfinite(large[i]) & (n_large[i] > 0)

    below_pct = np.where(ok_below, below * 100, 0.0)
    large_pct = np.where(ok_large, large * 100, 0.0)
    xs = [x_wui, x_wild]
    colors = [C_WUI, C_WILD]

    ymax = float(np.nanmax(np.r_[below_pct[ok_below], large_pct[ok_large], 0.0]))
    extra = 16.0 if ymax >= 40 else max(9.0, ymax * 0.40)
    y_top = max(22.0, ymax + extra)
    ax.set_ylim(0, y_top)
    pad = y_top * 0.016

    for j, (x, color) in enumerate(zip(xs, colors)):
        bars1 = ax.bar(
            x - w / 2,
            below_pct[:, j],
            width=w,
            color=color,
            edgecolor="0.25",
            lw=0.6,
            zorder=2,
        )
        bars2 = ax.bar(
            x + w / 2,
            large_pct[:, j],
            width=w,
            color=color,
            edgecolor="0.25",
            lw=0.6,
            hatch="///",
            alpha=0.75,
            zorder=2,
        )
        for i, b in enumerate(bars1):
            if ok_below[i, j]:
                _bar_label(ax, b, n_wind_below[i, j], n_below[i, j], pad)
        for i, b in enumerate(bars2):
            if ok_large[i, j]:
                _bar_label(ax, b, n_large_wind[i, j], n_large[i, j], pad)

    for i in range(len(SEASONS) - 1):
        ax.axvline(
            (centers[i] + centers[i + 1]) / 2,
            color="0.85",
            lw=0.6,
            zorder=1,
        )

    tick_x = np.column_stack([x_wui, x_wild]).ravel()
    ax.set_xticks(tick_x)
    ax.set_xticklabels(["WUI", "Wildland"] * len(SEASONS), fontsize=8)
    ax.tick_params(axis="x", length=3, pad=2)
    ax.set_xlim(centers[0] - 1.15, centers[-1] + 1.15)

    for x, season in zip(centers, SEASONS):
        ax.annotate(
            f"{season}\n({SEASON_SUB[season]})",
            xy=(x, 0),
            xycoords=("data", "axes fraction"),
            xytext=(0, -22),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=9,
            color="0.15",
        )

    ax.set_ylabel(ylabel, fontsize=10)
    ax.tick_params(axis="y", direction="out", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    socal = season_summary(_load(SOCAL_FLAGS, "in_analysis", "SAW"))
    norcal = season_summary(_load(NORCAL_FLAGS, "has_wx", "is_Diablo"))
    socal.insert(0, "region", "Southern California")
    norcal.insert(0, "region", "Northern California")
    out = pd.concat([socal, norcal], ignore_index=True)
    out.to_csv(OUT_CSV, index=False)

    for name, df in [("SoCal", socal), ("NorCal", norcal)]:
        print(f"{name} by season (p90 within season × FireType):")
        for _, r in df.iterrows():
            fb = "NA" if pd.isna(r.frac_below) else f"{100 * r.frac_below:.1f}%"
            fl = "NA" if pd.isna(r.frac_large) else f"{100 * r.frac_large:.1f}%"
            print(
                f"  {r.season:6s} {r.FireType:8s}  <p90 {fb:>6s} "
                f"({int(r.n_wind_below)}/{int(r.n_below)}) | "
                f"≥p90 {fl:>6s} ({int(r.n_large_wind)}/{int(r.n_large)})  "
                f"n={int(r.n)}"
            )

    fig, axes = plt.subplots(2, 1, figsize=(8.8, 6.6), facecolor="white")
    _panel(axes[0], socal, "Fraction of ignitions on\nSanta Ana days (%)")
    _panel(axes[1], norcal, "Fraction of ignitions on\nDiablo-like days (%)")

    for ax, letter, region in zip(
        axes, ["a", "b"], ["Southern California", "Northern California"]
    ):
        ax.text(
            0.01,
            0.98,
            letter,
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
            ha="left",
        )
        ax.set_title(region, fontsize=11, pad=4)

    handles = [
        Patch(facecolor=C_WUI, edgecolor="0.25", label="Fire size < p90"),
        Patch(
            facecolor=C_WUI,
            edgecolor="0.25",
            hatch="///",
            alpha=0.75,
            label="Fire size ≥ p90",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=2,
        frameon=False,
        fontsize=9,
        bbox_to_anchor=(0.5, 1.02),
    )

    fig.tight_layout(rect=[0, 0.02, 1, 0.96], h_pad=1.6)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)
    print("Saved", OUT_CSV)


if __name__ == "__main__":
    main()

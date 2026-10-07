#!/usr/bin/env python3
"""
Fig 5 style: Shrub & Grassland RAP shrub NPP and herb NPP separately by decade.

Layout: 3 decades (rows) × 2 metrics (columns: shrub | herb).
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

BASE = Path(__file__).resolve().parents[2]
FUEL_CSV = BASE / "dataPrc/WUI_vs_wildland_fuel_abundance_by_fire.csv"
OUT_PNG = BASE / "Fig/Fig5_RAP_shrub_herb_split.png"
OUT_PDF = BASE / "Fig/Fig5_RAP_shrub_herb_split.pdf"

COLORS = {
    "WUI": np.array([216, 118, 89]) / 255,
    "Wildland": np.array([41, 157, 143]) / 255,
}
DECADES = ["2000s", "2010s", "2020s"]
FIRE_TYPES = ["WUI", "Wildland"]

# (column key, title, bin_w, xlim) — values in g C m⁻² yr⁻¹ (RAP DN / 10)
METRICS = [
    ("rap_shr_npp", r"RAP shrub NPP (g C m$^{-2}$ yr$^{-1}$)", 8.0, (0, 320)),
    ("rap_herb_npp", r"RAP herb NPP (g C m$^{-2}$ yr$^{-1}$)", 12.0, (0, 450)),
]


def load_data() -> pd.DataFrame:
    g = gpd.read_file(BASE / "dataPrc/firePrmt/CalFire.shp")
    fuel = pd.read_csv(FUEL_CSV)[["fid", "rap_shr_npp", "rap_herb_npp"]].drop_duplicates(
        "fid"
    )
    df = g.drop(columns="geometry").merge(fuel, on="fid", how="left")
    # RAP DN (kg C ha⁻¹) → g C m⁻² yr⁻¹
    df["rap_shr_npp"] = df["rap_shr_npp"].astype(float) / 10.0
    df["rap_herb_npp"] = df["rap_herb_npp"].astype(float) / 10.0
    return df


def annotate_gap(ax, x0, x1, y, label):
    lo, hi = (x0, x1) if x0 <= x1 else (x1, x0)
    ax.plot([lo, hi], [y, y], "-", color="0.5", lw=1.0, zorder=4)
    dy = 0.004 * ax.get_ylim()[1] / 0.15
    ax.plot([lo, lo], [y - dy, y + dy], "-", color="0.5", lw=1.0, zorder=4)
    ax.plot([hi, hi], [y - dy, y + dy], "-", color="0.5", lw=1.0, zorder=4)
    ax.text(
        0.5 * (lo + hi),
        y + 3.5 * dy,
        label,
        ha="center",
        va="bottom",
        fontsize=8,
        color="k",
        zorder=5,
    )


def plot_panel(ax, sg, decade, col, bin_w, xlim, show_legend=False, show_decade=False):
    x_min, x_max = xlim
    edges = np.arange(x_min, x_max + bin_w, bin_w)
    x_fit = np.linspace(x_min, x_max, 500)
    y_lim = 0.18

    means, p95s, groups = {}, {}, {}
    for ft in FIRE_TYPES:
        vals = sg.loc[
            (sg["decade"] == decade) & (sg["FireType"] == ft), col
        ].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals)]
        groups[ft] = vals
        if len(vals) < 3:
            continue
        mu, sigma = float(np.mean(vals)), float(np.std(vals, ddof=1))
        if not np.isfinite(sigma) or sigma == 0:
            sigma = 1.0
        means[ft] = mu
        p95s[ft] = float(np.percentile(vals, 95))
        c = COLORS[ft]
        ax.hist(
            vals,
            bins=edges,
            weights=np.ones_like(vals) / len(vals),
            color=0.8 * c,
            alpha=0.25,
            edgecolor="none",
            zorder=1,
        )
        ax.plot(x_fit, sps.norm.pdf(x_fit, mu, sigma) * bin_w, color=c, lw=1.5, label=ft, zorder=3)
        ax.axvline(mu, color=c, ls=":", lw=1.4, zorder=2)
        ax.axvline(p95s[ft], color=c, ls="--", lw=1.1, zorder=2)

    if "WUI" in means and "Wildland" in means:
        y_ann = 0.12
        annotate_gap(
            ax,
            means["Wildland"],
            means["WUI"],
            y_ann,
            f"{means['WUI'] - means['Wildland']:.0f}",
        )
        annotate_gap(
            ax,
            p95s["Wildland"],
            p95s["WUI"],
            y_ann,
            f"{p95s['WUI'] - p95s['Wildland']:.0f}",
        )

    u, w = groups.get("WUI", np.array([])), groups.get("Wildland", np.array([]))
    p = np.nan
    if len(u) >= 3 and len(w) >= 3:
        _, p = sps.ttest_ind(u, w, equal_var=False, nan_policy="omit")
        p_str = r"$p$ < 0.01" if p < 0.01 else rf"$p$ = {p:.2f}"
        ax.text(-0.02, 0.50, p_str, transform=ax.transAxes, fontsize=9, va="center")

    ax.set_xlim(*xlim)
    ax.set_ylim(0, y_lim)
    ax.set_yticks([0, 0.05, 0.10, 0.15])
    if show_decade:
        ax.text(0.97, 0.82, decade, transform=ax.transAxes, ha="right", va="center", fontsize=10)
    ax.tick_params(direction="out", labelsize=8)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    if show_legend:
        ax.legend(loc="upper right", frameon=False, fontsize=8, handlelength=1.4)
    return p


def main():
    df = load_data()
    sg = df[df["lc"] == "ShrubGrass"].copy()

    fig, axes = plt.subplots(
        len(DECADES),
        len(METRICS),
        figsize=(8.2, 5.6),
        sharey=True,
        facecolor="w",
    )
    fig.subplots_adjust(hspace=0.15, wspace=0.18, left=0.10, right=0.98, top=0.90, bottom=0.10)
    fig.suptitle("Shrub & Grassland", fontsize=12, fontweight="normal", y=0.97)

    for j, (col, title, bin_w, xlim) in enumerate(METRICS):
        axes[0, j].set_title(title, fontsize=11, pad=6)
        for i, decade in enumerate(DECADES):
            ax = axes[i, j]
            plot_panel(
                ax,
                sg,
                decade,
                col,
                bin_w,
                xlim,
                show_legend=(i == 0 and j == 1),
                show_decade=True,
            )
            if i == len(DECADES) - 1:
                ax.set_xlabel(title.replace("RAP ", ""), fontsize=10)
            else:
                ax.tick_params(labelbottom=False)

    fig.supylabel("Probability", fontsize=11)
    fig.savefig(OUT_PNG, dpi=300, facecolor="w", bbox_inches="tight")
    fig.savefig(OUT_PDF, facecolor="w", bbox_inches="tight")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)

    for col, title, *_ in METRICS:
        print(f"\n=== {title} ===")
        for decade in DECADES:
            u = sg.loc[(sg.decade == decade) & (sg.FireType == "WUI"), col].dropna()
            w = sg.loc[(sg.decade == decade) & (sg.FireType == "Wildland"), col].dropna()
            _, p = sps.ttest_ind(u, w, equal_var=False, nan_policy="omit")
            print(
                f"  {decade}: Δmean={u.mean()-w.mean():+.0f}  "
                f"Δp95={u.quantile(0.95)-w.quantile(0.95):+.0f}  p={p:.3g}"
            )


if __name__ == "__main__":
    main()

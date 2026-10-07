#!/usr/bin/env python3
"""
Fig 5-style decade distributions for Shrub & Grassland fires:
  1) RAP shrub NPP
  2) RAP herb NPP
  3) ESA CCI AGB
Each saved as its own 3-panel (2000s / 2010s / 2020s) figure.
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

BASE = Path(__file__).resolve().parents[2]
FUEL_CSV = BASE / "dataPrc/WUI_vs_wildland_fuel_abundance_by_fire.csv"
OUT_DIR = BASE / "Fig"

COLORS = {
    "WUI": np.array([216, 118, 89]) / 255,
    "Wildland": np.array([41, 157, 143]) / 255,
}
DECADES = ["2000s", "2010s", "2020s"]
FIRE_TYPES = ["WUI", "Wildland"]

# key, xlabel, outfile stem, bin_w, xlim, y_lim, gap_fmt
SPECS = [
    (
        "rap_shr_npp",
        r"RAP shrub NPP (g C m$^{-2}$ yr$^{-1}$)",
        "Fig5_RAP_shrub_NPP",
        8.0,
        (0, 320),
        0.18,
        "{:.0f}",
    ),
    (
        "rap_herb_npp",
        r"RAP herb NPP (g C m$^{-2}$ yr$^{-1}$)",
        "Fig5_RAP_herb_NPP",
        12.0,
        (0, 450),
        0.15,
        "{:.0f}",
    ),
    (
        "cci_agb",
        r"ESA CCI AGB (Mg ha$^{-1}$)",
        "Fig5_ESA_CCI_AGB",
        3.0,
        (0, 90),
        0.20,
        "{:.1f}",
    ),
]


def load_data() -> pd.DataFrame:
    g = gpd.read_file(BASE / "dataPrc/firePrmt/CalFire.shp")
    fuel = pd.read_csv(FUEL_CSV)[
        ["fid", "rap_shr_npp", "rap_herb_npp", "cci_agb"]
    ].drop_duplicates("fid")
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
        fontsize=9,
        color="k",
        zorder=5,
    )


def draw_one(df: pd.DataFrame, col: str, xlabel: str, stem: str, bin_w, xlim, y_lim, gap_fmt):
    sg = df[df["lc"] == "ShrubGrass"].copy()
    x_min, x_max = xlim
    edges = np.arange(x_min, x_max + bin_w, bin_w)
    x_fit = np.linspace(x_min, x_max, 600)

    fig, axes = plt.subplots(
        len(DECADES), 1, figsize=(4.2, 5.4), sharex=True, sharey=True, facecolor="w"
    )
    fig.subplots_adjust(hspace=0.12, left=0.18, right=0.96, top=0.92, bottom=0.10)

    print(f"\n=== {xlabel} ===")
    for ax, decade in zip(axes, DECADES):
        means, p95s, groups = {}, {}, {}
        for ft in FIRE_TYPES:
            vals = sg.loc[
                (sg["decade"] == decade) & (sg["FireType"] == ft), col
            ].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            groups[ft] = vals
            if len(vals) < 3:
                continue
            mu = float(np.mean(vals))
            sigma = float(np.std(vals, ddof=1))
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
            ax.plot(
                x_fit,
                sps.norm.pdf(x_fit, mu, sigma) * bin_w,
                color=c,
                lw=1.5,
                label=ft,
                zorder=3,
            )
            ax.axvline(mu, color=c, ls=":", lw=1.5, zorder=2)
            ax.axvline(p95s[ft], color=c, ls="--", lw=1.2, zorder=2)

        if "WUI" in means and "Wildland" in means:
            y_ann = min(0.10, 0.65 * y_lim)
            annotate_gap(
                ax,
                means["Wildland"],
                means["WUI"],
                y_ann,
                gap_fmt.format(means["WUI"] - means["Wildland"]),
            )
            annotate_gap(
                ax,
                p95s["Wildland"],
                p95s["WUI"],
                y_ann,
                gap_fmt.format(p95s["WUI"] - p95s["Wildland"]),
            )

        u, w = groups.get("WUI", np.array([])), groups.get("Wildland", np.array([]))
        if len(u) >= 3 and len(w) >= 3:
            _, p = sps.ttest_ind(u, w, equal_var=False, nan_policy="omit")
            p_str = r"$p$ < 0.01" if p < 0.01 else rf"$p$ = {p:.2f}"
            ax.text(-0.02, 0.50, p_str, transform=ax.transAxes, fontsize=10, va="center")
            print(
                f"  {decade}: nU/nW={len(u)}/{len(w)}  "
                f"Δmean={means['WUI']-means['Wildland']:+.2f}  "
                f"Δp95={p95s['WUI']-p95s['Wildland']:+.2f}  p={p:.3g}"
            )
        else:
            print(f"  {decade}: insufficient sample (nU/nW={len(u)}/{len(w)})")

        ax.set_xlim(*xlim)
        ax.set_ylim(0, y_lim)
        ax.set_yticks([0, 0.05, 0.10, 0.15] if y_lim <= 0.18 else [0, 0.05, 0.10, 0.15, 0.20])
        ax.text(0.97, 0.80, decade, transform=ax.transAxes, ha="right", va="center", fontsize=10)
        ax.tick_params(direction="out")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

        if decade == DECADES[0]:
            ax.set_title("Shrub & Grassland", fontweight="normal", fontsize=12)
            ax.legend(loc="center right", frameon=False, fontsize=9, handlelength=1.5)

    axes[-1].set_xlabel(xlabel)
    fig.supylabel("Probability", fontsize=11)

    out_png = OUT_DIR / f"{stem}.png"
    out_pdf = OUT_DIR / f"{stem}.pdf"
    fig.savefig(out_png, dpi=300, facecolor="w", bbox_inches="tight")
    fig.savefig(out_pdf, facecolor="w", bbox_inches="tight")
    plt.close(fig)
    print("Saved", out_png)


def main():
    df = load_data()
    for spec in SPECS:
        draw_one(df, *spec)


if __name__ == "__main__":
    main()

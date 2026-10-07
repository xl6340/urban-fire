#!/usr/bin/env python3
"""
Fig 5 style: Shrub & Grassland RAP shrub/herb NPP distributions by decade.

WUI vs Wildland histograms + normal PDF, mean (dotted) and 95th percentile
(dashed) lines, with gap annotations — same layout as Fig5 NDVI panels.
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

BASE = Path(__file__).resolve().parents[2]
FUEL_CSV = BASE / "dataPrc/WUI_vs_wildland_fuel_abundance_by_fire.csv"
OUT_PNG = BASE / "Fig/Fig5_RAP_ShrubHerb_NPP.png"
OUT_PDF = BASE / "Fig/Fig5_RAP_ShrubHerb_NPP.pdf"

COLORS = {
    "WUI": np.array([216, 118, 89]) / 255,
    "Wildland": np.array([41, 157, 143]) / 255,
}
DECADES = ["2000s", "2010s", "2020s"]
FIRE_TYPES = ["WUI", "Wildland"]


def load_data() -> pd.DataFrame:
    g = gpd.read_file(BASE / "dataPrc/firePrmt/CalFire.shp")
    fuel = pd.read_csv(FUEL_CSV)[["fid", "rap_shr_npp", "rap_herb_npp"]].drop_duplicates(
        "fid"
    )
    df = g.drop(columns="geometry").merge(fuel, on="fid", how="left")
    # RAP DN (kg C ha⁻¹) → g C m⁻² yr⁻¹
    df["rap"] = (
        df["rap_shr_npp"].astype(float) + df["rap_herb_npp"].astype(float)
    ) / 10.0
    return df


def annotate_gap(ax, x0, x1, y, label):
    """Horizontal bar between x0 and x1 with value centered above."""
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


def main():
    df = load_data()
    sg = df[(df["lc"] == "ShrubGrass") & np.isfinite(df["rap"])].copy()

    # Binning on RAP scale in g C m⁻² yr⁻¹
    bin_w = 15.0
    x_min, x_max = 0.0, 550.0
    edges = np.arange(x_min, x_max + bin_w, bin_w)
    x_fit = np.linspace(x_min, x_max, 600)
    y_lim = 0.15

    fig, axes = plt.subplots(
        len(DECADES), 1, figsize=(4.2, 5.4), sharex=True, sharey=True, facecolor="w"
    )
    fig.subplots_adjust(hspace=0.12, left=0.18, right=0.96, top=0.92, bottom=0.10)

    for ax, decade in zip(axes, DECADES):
        means = {}
        p95s = {}
        groups = {}
        for ft in FIRE_TYPES:
            vals = sg.loc[
                (sg["decade"] == decade) & (sg["FireType"] == ft), "rap"
            ].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            groups[ft] = vals
            if len(vals) < 3:
                continue
            mu, sigma = float(np.mean(vals)), float(np.std(vals, ddof=1))
            means[ft] = mu
            p95s[ft] = float(np.percentile(vals, 95))
            c = COLORS[ft]
            ax.hist(
                vals,
                bins=edges,
                density=False,
                weights=np.ones_like(vals) / len(vals),
                color=0.8 * c,
                alpha=0.25,
                edgecolor="none",
                zorder=1,
            )
            y_fit = sps.norm.pdf(x_fit, mu, sigma) * bin_w
            ax.plot(x_fit, y_fit, color=c, lw=1.5, label=ft, zorder=3)
            ax.axvline(mu, color=c, ls=":", lw=1.5, zorder=2)
            ax.axvline(p95s[ft], color=c, ls="--", lw=1.2, zorder=2)

        # Gap annotations (WUI − Wildland), same as Fig5.m
        if "WUI" in means and "Wildland" in means:
            y_ann = 0.10
            mean_gap = means["WUI"] - means["Wildland"]
            annotate_gap(
                ax, means["Wildland"], means["WUI"], y_ann, f"{mean_gap:.0f}"
            )
            pct_gap = p95s["WUI"] - p95s["Wildland"]
            annotate_gap(ax, p95s["Wildland"], p95s["WUI"], y_ann, f"{pct_gap:.0f}")

        u, w = groups.get("WUI", np.array([])), groups.get("Wildland", np.array([]))
        if len(u) >= 3 and len(w) >= 3:
            _, p = sps.ttest_ind(u, w, equal_var=False, nan_policy="omit")
            p_str = r"$p$ < 0.01" if p < 0.01 else rf"$p$ = {p:.2f}"
        else:
            p_str = ""
        if p_str:
            ax.text(-0.02, 0.50, p_str, transform=ax.transAxes, fontsize=10, va="center")

        ax.set_xlim(20, 520)
        ax.set_ylim(0, y_lim)
        ax.set_yticks([0, 0.05, 0.10, 0.15])
        ax.text(
            0.97, 0.80, decade, transform=ax.transAxes, ha="right", va="center", fontsize=10
        )
        ax.tick_params(direction="out")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

        if decade == DECADES[0]:
            ax.set_title("Shrub & Grassland", fontweight="normal", fontsize=12)
            ax.legend(loc="center right", frameon=False, fontsize=9, handlelength=1.5)

    axes[-1].set_xlabel(r"RAP shrub/herb NPP (g C m$^{-2}$ yr$^{-1}$)")
    fig.supylabel("Probability", fontsize=11)
    fig.savefig(OUT_PNG, dpi=300, facecolor="w", bbox_inches="tight")
    fig.savefig(OUT_PDF, facecolor="w", bbox_inches="tight")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)

    # summary
    for decade in DECADES:
        u = sg.loc[(sg.decade == decade) & (sg.FireType == "WUI"), "rap"]
        w = sg.loc[(sg.decade == decade) & (sg.FireType == "Wildland"), "rap"]
        _, p = sps.ttest_ind(u, w, equal_var=False, nan_policy="omit")
        print(
            f"{decade}: nU/nW={u.notna().sum()}/{w.notna().sum()}  "
            f"Δmean={u.mean()-w.mean():+.0f}  Δp95={u.quantile(0.95)-w.quantile(0.95):+.0f}  "
            f"p={p:.3g}"
        )


if __name__ == "__main__":
    main()

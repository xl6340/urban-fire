#!/usr/bin/env python3
"""
Forest fires: WUI vs wildland distributions of 4-day Canadian FWI
moisture codes (FFMC, DMC, DC). Same sample as Fig. 3B Forest panel.
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

BASE = Path(__file__).resolve().parents[2]
OUT_PNG = BASE / "Fig/Fig3B_Forest_FWI_dist.png"
OUT_PDF = BASE / "Fig/Fig3B_Forest_FWI_dist.pdf"

COLORS = {
    "WUI": np.array([216, 118, 89]) / 255,
    "Wildland": np.array([41, 157, 143]) / 255,
}

# (col, xlabel, bin_w, xlim, gap_fmt, overlay)
# FFMC is too peaked/skewed for a Gaussian; DMC/DC use a normal overlay.
PANELS = [
    ("FFMC", "FFMC", 2.0, (80, 110), "{:.1f}", "kde"),
    ("DMC", "DMC", 40.0, (0, 900), "{:.0f}", "norm"),
    ("DC", "DC", 80.0, (0, 1800), "{:.0f}", "norm"),
]


def load_forest() -> pd.DataFrame:
    g = gpd.read_file(BASE / "dataPrc/firePrmt/CalFire.shp")
    fwi = pd.read_csv(BASE / "dataPrc/CalFire_FWI_4day.csv")
    df = g.drop(columns="geometry").merge(
        fwi[["fire_id", "FFMC_4d", "DMC_4d", "DC_4d"]],
        left_on="fid",
        right_on="fire_id",
        how="left",
    )
    df["FFMC"] = df["FFMC_4d"]
    df["DMC"] = df["DMC_4d"]
    df["DC"] = df["DC_4d"]
    return df[df["lc"] == "Forest"].copy()


def annotate_gap(ax, x0, x1, y, label):
    lo, hi = (x0, x1) if x0 <= x1 else (x1, x0)
    ax.plot([lo, hi], [y, y], "-", color="0.45", lw=1.0, zorder=4)
    dy = 0.012 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.plot([lo, lo], [y - dy, y + dy], "-", color="0.45", lw=1.0, zorder=4)
    ax.plot([hi, hi], [y - dy, y + dy], "-", color="0.45", lw=1.0, zorder=4)
    ax.text(0.5 * (lo + hi), y + 2.8 * dy, label, ha="center", va="bottom", fontsize=9, zorder=5)


def main():
    fo = load_forest()
    fig, axes = plt.subplots(3, 1, figsize=(4.2, 6.2), facecolor="w")
    fig.subplots_adjust(hspace=0.38, left=0.20, right=0.96, top=0.92, bottom=0.08)

    for i, (ax, (col, xlabel, bin_w, xlim, gap_fmt, overlay)) in enumerate(zip(axes, PANELS)):
        x_min, x_max = xlim
        edges = np.arange(x_min, x_max + bin_w * 0.5, bin_w)
        x_fit = np.linspace(x_min, x_max, 500)
        means = {}
        for ft in ("WUI", "Wildland"):
            vals = fo.loc[fo["FireType"] == ft, col].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            mu, sigma = float(np.mean(vals)), float(np.std(vals, ddof=1))
            means[ft] = mu
            c = COLORS[ft]
            ax.hist(
                vals,
                bins=edges,
                weights=np.ones_like(vals) / len(vals),
                color=0.8 * c,
                alpha=0.28,
                edgecolor="none",
                zorder=1,
            )
            if overlay == "kde":
                y_fit = sps.gaussian_kde(vals)(x_fit) * bin_w
            else:
                y_fit = sps.norm.pdf(x_fit, mu, max(sigma, 1e-9)) * bin_w
            ax.plot(x_fit, y_fit, color=c, lw=1.6, label=ft, zorder=3)
            ax.axvline(mu, color=c, ls=":", lw=1.4, zorder=2)

        u = fo.loc[fo.FireType == "WUI", col].dropna()
        w = fo.loc[fo.FireType == "Wildland", col].dropna()
        _, p = sps.ttest_ind(u, w, equal_var=False, nan_policy="omit")
        p_str = r"$p$ < 0.01" if p < 0.01 else rf"$p$ = {p:.2f}"
        ax.text(0.08, 0.52, p_str, transform=ax.transAxes, fontsize=10, va="center", ha="left")

        y_max = ax.get_ylim()[1]
        ax.set_ylim(0, y_max * 1.12)
        annotate_gap(
            ax,
            means["Wildland"],
            means["WUI"],
            0.72 * ax.get_ylim()[1],
            gap_fmt.format(means["WUI"] - means["Wildland"]),
        )

        ax.set_xlim(*xlim)
        ax.set_xlabel(xlabel, fontsize=10)
        ax.text(0.97, 0.80, xlabel, transform=ax.transAxes, ha="right", va="center", fontsize=10)
        ax.tick_params(direction="out", labelsize=8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        print(
            f"{xlabel}: nU/nW={len(u)}/{len(w)}  "
            f"mean {means['WUI']:.3g} vs {means['Wildland']:.3g}  "
            f"Δ={means['WUI']-means['Wildland']:+.3g}  p={p:.3g}"
        )

        if i == 0:
            ax.set_title("Forest", fontweight="normal", fontsize=12)
            ax.legend(loc="upper left", frameon=False, fontsize=9, handlelength=1.5)

    fig.supylabel("Probability", fontsize=11)
    fig.savefig(OUT_PNG, dpi=300, facecolor="w", bbox_inches="tight")
    fig.savefig(OUT_PDF, facecolor="w", bbox_inches="tight")
    print("Saved", OUT_PNG)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Shrub & Grassland fires: WUI vs wildland vegetation distributions
(NDVI, RAP shrub/herb NPP, ESA CCI AGB). Pooled across years.
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import NullLocator
from scipy import stats as sps

BASE = Path(__file__).resolve().parents[2]
FUEL_CSV = BASE / "dataPrc/WUI_vs_wildland_fuel_abundance_by_fire.csv"
OUT_PNG = BASE / "Fig/Fig_SG_WUI_vs_wildland_veg_dist.png"
OUT_PDF = BASE / "Fig/Fig_SG_WUI_vs_wildland_veg_dist.pdf"
OUT_S8_PNG = BASE / "Fig/FigS8.png"
OUT_S8_PDF = BASE / "Fig/FigS8.pdf"

COLORS = {
    "WUI": np.array([216, 118, 89]) / 255,
    "Wildland": np.array([41, 157, 143]) / 255,
}

# RAP npp-partitioned-v3 stores DN ≈ kg C ha⁻¹ (= 0.1 g C m⁻²); divide by 10 for g C m⁻².
# (col, xlabel, bin_w, xlim, gap_fmt, overlay)
PANELS = [
    ("ndvi", "NDVI", 0.04, (0.10, 0.80), "{:.2f}", "norm"),
    ("rap", r"RAP shrub/herb NPP (g C m$^{-2}$ yr$^{-1}$)", 15.0, (0, 520), "{:.0f}", "norm"),
    ("cci_agb", r"ESA CCI AGB (Mg ha$^{-1}$)", 3.0, (0, 90), "{:.1f}", "lognorm"),
]


def load_data() -> pd.DataFrame:
    g = gpd.read_file(BASE / "dataPrc/firePrmt/CalFire.shp")
    fuel = pd.read_csv(FUEL_CSV)[
        ["fid", "rap_shr_npp", "rap_herb_npp", "cci_agb"]
    ].drop_duplicates("fid")
    df = g.drop(columns="geometry").merge(fuel, on="fid", how="left")
    # Convert RAP DN (kg C ha⁻¹) → g C m⁻² yr⁻¹
    df["rap"] = (
        df["rap_shr_npp"].astype(float) + df["rap_herb_npp"].astype(float)
    ) / 10.0
    return df[df["lc"] == "ShrubGrass"].copy()


def annotate_gap(ax, x0, x1, y, label, log_x=False):
    lo, hi = (x0, x1) if x0 <= x1 else (x1, x0)
    ax.plot([lo, hi], [y, y], "-", color="0.45", lw=1.0, zorder=4)
    dy = 0.012 * (ax.get_ylim()[1] - ax.get_ylim()[0])
    ax.plot([lo, lo], [y - dy, y + dy], "-", color="0.45", lw=1.0, zorder=4)
    ax.plot([hi, hi], [y - dy, y + dy], "-", color="0.45", lw=1.0, zorder=4)
    xmid = float(np.sqrt(lo * hi)) if log_x else 0.5 * (lo + hi)
    ax.text(xmid, y + 2.8 * dy, label, ha="center", va="bottom", fontsize=9, zorder=5)


def main():
    sg = load_data()
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.2), facecolor="w")
    fig.subplots_adjust(wspace=0.28, left=0.07, right=0.98, top=0.82, bottom=0.22)
    if not hasattr(axes, "__len__"):
        axes = [axes]
    fig.suptitle("Shrub & Grassland fires", fontsize=12, fontweight="normal", y=0.98)

    for ax, (col, xlabel, bin_w, xlim, gap_fmt, overlay), letter in zip(
        axes, PANELS, ("(A)", "(B)", "(C)")
    ):
        x_min, x_max = xlim
        log_x = overlay == "lognorm"
        if log_x:
            x_min = max(x_min, 1.0)
            n_bins = 18
            edges = np.logspace(np.log10(x_min), np.log10(x_max), n_bins + 1)
            x_fit = np.logspace(np.log10(x_min), np.log10(x_max), 400)
            dlog10 = float(np.log10(edges[1]) - np.log10(edges[0]))
        else:
            edges = np.arange(x_min, x_max + bin_w * 0.5, bin_w)
            x_fit = np.linspace(x_min, x_max, 500)
            dlog10 = np.nan
        means = {}
        for ft in ("WUI", "Wildland"):
            vals = sg.loc[sg["FireType"] == ft, col].to_numpy(dtype=float)
            vals = vals[np.isfinite(vals)]
            mu, sigma = float(np.mean(vals)), float(np.std(vals, ddof=1))
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
                loc = mu
            elif overlay == "lognorm":
                pos = vals[vals > 0]
                log_mu = float(np.mean(np.log(pos)))
                log_sd = float(np.std(np.log(pos), ddof=1))
                log10_mu = log_mu / np.log(10.0)
                log10_sd = max(log_sd / np.log(10.0), 1e-9)
                y_fit = sps.norm.pdf(np.log10(x_fit), log10_mu, log10_sd) * dlog10
                loc = float(np.exp(log_mu))  # geometric mean
            else:
                y_fit = sps.norm.pdf(x_fit, mu, max(sigma, 1e-9)) * bin_w
                loc = mu
            means[ft] = loc
            ax.plot(x_fit, y_fit, color=c, lw=1.6, label=ft, zorder=3)
            ax.axvline(loc, color=c, ls=":", lw=1.4, zorder=2)

        u = sg.loc[sg.FireType == "WUI", col].dropna()
        w = sg.loc[sg.FireType == "Wildland", col].dropna()
        if overlay == "lognorm":
            _, p = sps.ttest_ind(np.log(u[u > 0]), np.log(w[w > 0]), equal_var=False)
        else:
            _, p = sps.ttest_ind(u, w, equal_var=False, nan_policy="omit")
        p_str = r"$p$ < 0.01" if p < 0.01 else rf"$p$ = {p:.2f}"
        if log_x:
            ax.text(0.04, 0.36, p_str, transform=ax.transAxes, fontsize=10, va="center", ha="left")
        else:
            ax.text(0.97, 0.50, p_str, transform=ax.transAxes, fontsize=10, va="center", ha="right")

        y_max = ax.get_ylim()[1]
        ax.set_ylim(0, y_max * 1.08)
        annotate_gap(
            ax,
            means["Wildland"],
            means["WUI"],
            0.72 * ax.get_ylim()[1],
            gap_fmt.format(means["WUI"] - means["Wildland"]),
            log_x=log_x,
        )

        ax.set_xlim(x_min, x_max)
        if log_x:
            ax.set_xscale("log")
            ax.set_xticks([1, 3, 10, 30, 90])
            ax.set_xticklabels(["1", "3", "10", "30", "90"])
            ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xlabel(xlabel, fontsize=10)
        ax.set_ylabel("Probability" if ax is axes[0] else "", fontsize=10)
        ax.tick_params(direction="out", labelsize=8)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.text(
            -0.10,
            1.06,
            letter,
            transform=ax.transAxes,
            fontsize=12,
            va="bottom",
            ha="left",
            clip_on=False,
        )
        kind = "geomean" if overlay == "lognorm" else "mean"
        print(
            f"{xlabel}: nU/nW={len(u)}/{len(w)}  "
            f"{kind} {means['WUI']:.3g} vs {means['Wildland']:.3g}  "
            f"Δ={means['WUI']-means['Wildland']:+.3g}  p={p:.3g}"
        )

    axes[-1].legend(loc="upper left", frameon=False, fontsize=9, handlelength=1.4)
    for out in (OUT_PNG, OUT_S8_PNG):
        fig.savefig(out, dpi=300, facecolor="w", bbox_inches="tight")
    for out in (OUT_PDF, OUT_S8_PDF):
        fig.savefig(out, facecolor="w", bbox_inches="tight")
    print("Saved", OUT_PNG)
    print("Saved", OUT_S8_PNG)


if __name__ == "__main__":
    main()

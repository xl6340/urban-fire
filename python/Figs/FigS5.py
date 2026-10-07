#!/usr/bin/env python3
"""
Fig. S5: WUI vs wildland environmental distributions, forest-dominated fires.

3 x 4 panels: original climate / vegetation / terrain variables plus 4-day
mean Canadian FWI moisture codes (FFMC, DMC, DC). Style matches
code/FigSI_histogram.m (PDF histogram + Gaussian overlay + mean lines).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

BASE = Path(__file__).resolve().parents[2]
OUT_PNG = BASE / "Fig/FigS5.png"
OUT_PDF = BASE / "Fig/FigS5.pdf"
OUT_PNG_LEGACY = BASE / "Fig/FigSI_histogram_forest.png"

C_WILD = np.array([66, 157, 143]) / 255
C_WUI = np.array([231, 111, 81]) / 255

# (column, xlabel, xlim, vpd_hpa_to_kpa)
PANELS = [
    ("vpdmax", r"$VPD_{max}$ (kPa)", (0, 8), True),
    ("tmean", r"$T_{mean}$ (°C)", (0, 40), False),
    ("tmin", r"$T_{min}$ (°C)", (-5, 30), False),
    ("tmax", r"$T_{max}$ (°C)", (0, 50), False),
    ("ppt", "Precipitation (mm)", (0, 3000), False),
    ("vs", "Wind speed (m/s)", (0, 10), False),
    ("ndviM", "NDVI", (0.2, 1.0), False),
    ("slope", "Slope (°)", (0, 35), False),
    ("elevation", "Elevation (m)", (-500, 3000), False),
    ("FFMC", "FFMC", (70, 101), False),
    ("DMC", "DMC", (0, 900), False),
    ("DC", "DC", (0, 1800), False),
]
LETTERS = list("ABCDEFGHIJKL")


def load_forest() -> tuple[pd.DataFrame, pd.DataFrame]:
    urban = pd.read_csv(BASE / "dataFig/variable/Forest_Urban-edge_withFWI.csv")
    wild = pd.read_csv(BASE / "dataFig/variable/Forest_Wildland_withFWI.csv")
    return urban, wild


def p_string(p: float) -> str:
    if p < 0.001:
        return r"$\mathit{p}$<0.001"
    if p < 0.01:
        return rf"$\mathit{{p}}$={p:.3f}"
    return rf"$\mathit{{p}}$={p:.2f}"


def values(df: pd.DataFrame, col: str, to_kpa: bool) -> np.ndarray:
    v = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
    v = v[np.isfinite(v)]
    if to_kpa:
        v = v / 10.0
    return v


def main() -> None:
    urban, wild = load_forest()
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(3, 4, figsize=(13.2, 9.0), facecolor="white")
    fig.subplots_adjust(wspace=0.18, hspace=0.42, left=0.07, right=0.98, top=0.94, bottom=0.07)
    fig.supylabel("Probability Density", fontsize=12, x=0.015)

    legend_handles = None
    for ax, letter, (col, xlabel, xlim, to_kpa) in zip(axes.ravel(), LETTERS, PANELS):
        val_u = values(urban, col, to_kpa)
        val_w = values(wild, col, to_kpa)
        mu_u, sig_u = float(np.mean(val_u)), float(np.std(val_u, ddof=1))
        mu_w, sig_w = float(np.mean(val_w)), float(np.std(val_w, ddof=1))
        diff = mu_u - mu_w
        _, p = sps.ttest_ind(val_u, val_w, equal_var=False)

        x = np.linspace(xlim[0], xlim[1], 200)
        n_bins = 28
        h_w = ax.hist(
            val_w,
            bins=n_bins,
            range=xlim,
            density=True,
            color=C_WILD,
            edgecolor=(0.9, 0.9, 0.9),
            alpha=0.4,
            zorder=1,
        )
        ax.plot(x, sps.norm.pdf(x, mu_w, max(sig_w, 1e-9)), color=C_WILD, lw=2.0, zorder=3)
        ax.axvline(mu_w, color=C_WILD, ls=":", lw=2.0, zorder=2)

        h_u = ax.hist(
            val_u,
            bins=n_bins,
            range=xlim,
            density=True,
            color=C_WUI,
            edgecolor=(0.9, 0.9, 0.9),
            alpha=0.4,
            zorder=1,
        )
        ax.plot(x, sps.norm.pdf(x, mu_u, max(sig_u, 1e-9)), color=C_WUI, lw=2.0, zorder=3)
        ax.axvline(mu_u, color=C_WUI, ls=":", lw=2.0, zorder=2)

        if legend_handles is None:
            legend_handles = (h_w[2][0], h_u[2][0])

        ax.set_xlim(*xlim)
        ax.set_xlabel(xlabel, fontsize=10)
        ax.tick_params(axis="x", direction="out", labelsize=9, length=3.5)
        ax.tick_params(axis="y", length=0, labelleft=False)
        ax.yaxis.set_ticks([])
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_linewidth(0.8)
        ax.set_title(
            f"({letter})  Diff:{diff:+.1f}  ({p_string(p)})",
            loc="left",
            fontsize=10,
            fontweight="normal",
            pad=6,
        )
        print(
            f"({letter}) {col}: nU/nW={len(val_u)}/{len(val_w)}  "
            f"mean {mu_u:.3g} vs {mu_w:.3g}  Δ={diff:+.3g}  p={p:.3g}"
        )

    axes[0, 0].legend(
        legend_handles,
        ["Wildland", "Urban-edge"],
        loc="upper right",
        frameon=False,
        fontsize=9,
        handlelength=1.2,
        borderaxespad=0.2,
    )

    fig.savefig(OUT_PNG, dpi=300, facecolor="white")
    fig.savefig(OUT_PDF, facecolor="white")
    fig.savefig(OUT_PNG_LEGACY, dpi=300, facecolor="white")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)


if __name__ == "__main__":
    main()

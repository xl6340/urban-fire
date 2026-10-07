#!/usr/bin/env python3
"""Redraw Fig. S6 (forest fires by elevation quantile)."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
DAT = BASE / "dataFig/elevation"
OUT_PNG = BASE / "Fig/FigS6.png"
OUT_PDF = BASE / "Fig/FigS6.pdf"

C_WUI = np.array([216, 118, 89]) / 255
C_WILD = np.array([41, 157, 143]) / 255
C_HIST = "#7785ac"
C_Q = "#7209b7"

SPECS = [
    ("beta", r"$\beta$ value", (-1.6, -0.8), np.arange(-1.6, -0.79, 0.2)),
    ("VPDmax", r"$VPD_{max}$ (kPa)", (2, 5), np.arange(2, 6, 1)),
    ("Tmean", r"$T_{mean}$ ($^\circ$C)", (15, 25), np.arange(15, 26, 5)),
    ("Tmin", r"$T_{min}$ ($^\circ$C)", (8, 16), np.arange(8, 17, 4)),
    ("Tmax", r"$T_{max}$ ($^\circ$C)", (20, 35), np.arange(20, 36, 5)),
]


def style_ax(ax):
    ax.set_xlim(0, 3000)
    ax.tick_params(direction="out", length=4, width=0.8, labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_linewidth(0.8)


def errorbar_pair(ax, x, df):
    ax.errorbar(
        x,
        df["Mean_Urban-edge"],
        yerr=df["Err_Urban-edge"],
        fmt="-o",
        color=C_WUI,
        mfc=C_WUI,
        markersize=5,
        linewidth=1.2,
        capsize=3,
        capthick=0.8,
        elinewidth=0.8,
        label="WUI",
        zorder=3,
    )
    ax.errorbar(
        x,
        df["Mean_Wildland"],
        yerr=df["Err_Wildland"],
        fmt="-^",
        color=C_WILD,
        mfc=C_WILD,
        markersize=5,
        linewidth=1.2,
        capsize=3,
        capthick=0.8,
        elinewidth=0.8,
        label="Wildland",
        zorder=3,
    )


def _scale_hist_patches(patches, scale=1e4):
    for p in patches:
        p.set_height(p.get_height() * scale)


def main():
    elev = np.loadtxt(DAT / "elevations.csv")
    x = np.loadtxt(DAT / "eleCenters.csv")

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 10,
            "axes.labelsize": 11,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(3, 2, figsize=(6.2, 6.4), facecolor="white")
    fig.subplots_adjust(wspace=0.38, hspace=0.32, left=0.14, right=0.97, top=0.96, bottom=0.08)
    axes = axes.ravel()

    # Panel 1: pooled elevation PDF of forest fire perimeters
    ax = axes[0]
    _, _, patches = ax.hist(
        elev,
        bins=40,
        density=True,
        color=C_HIST,
        edgecolor=(0.9, 0.9, 0.9),
        alpha=0.45,
        linewidth=0.4,
        zorder=2,
    )
    _scale_hist_patches(patches)
    for xc in x:
        ax.axvline(xc, linestyle="-.", color=C_Q, linewidth=1.3, zorder=1)
    ax.set_ylabel(r"Density ($\times 10^{-4}$)")
    ax.set_ylim(0, 8)
    ax.set_yticks(np.arange(0, 9, 2))
    style_ax(ax)

    for i, (name, ylabel, ylim, yticks) in enumerate(SPECS, start=1):
        ax = axes[i]
        df = pd.read_csv(DAT / f"{name}.csv")
        errorbar_pair(ax, x, df)
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        ax.set_yticks(yticks)
        style_ax(ax)
        if i == 1:
            ax.legend(
                loc="upper right",
                frameon=False,
                fontsize=9,
                handlelength=1.6,
                borderpad=0.1,
                labelspacing=0.25,
            )

    for ax in axes[-2:]:
        ax.set_xlabel("Elevation (m)", fontsize=11)

    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=300)
    fig.savefig(OUT_PDF)
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_PDF}")


if __name__ == "__main__":
    main()

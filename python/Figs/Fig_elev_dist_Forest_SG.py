#!/usr/bin/env python3
"""Elevation PDFs: Forest and Shrub/grassland, WUI vs wildland (Fig. S5 style)."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

BASE = Path(__file__).resolve().parents[2]
VAR = BASE / "dataFig/variable"
OUT_PNG = BASE / "Fig/Fig_elev_dist_Forest_SG.png"
OUT_PDF = BASE / "Fig/Fig_elev_dist_Forest_SG.pdf"

C_WILD = np.array([66, 157, 143]) / 255
C_WUI = np.array([231, 111, 81]) / 255

XLIM = (-500, 3000)
N_BINS = 28

PANELS = [
    ("Forest", "Forest_Urban-edge_withFWI.csv", "Forest_Wildland_withFWI.csv"),
    ("Shrub/grassland", "ShrubGrass_Urban-edge_withFWI.csv", "ShrubGrass_Wildland_withFWI.csv"),
]


def p_string(p: float) -> str:
    if p < 0.001:
        return r"$\mathit{p}$ < 0.001"
    if p < 0.01:
        return rf"$\mathit{{p}}$ = {p:.3f}"
    return rf"$\mathit{{p}}$ = {p:.2f}"


def values(path: Path) -> np.ndarray:
    v = pd.to_numeric(pd.read_csv(path)["elevation"], errors="coerce").to_numpy(float)
    return v[np.isfinite(v)]


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.15), facecolor="white")
    fig.subplots_adjust(wspace=0.28, left=0.10, right=0.98, top=0.82, bottom=0.20)

    legend_handles = None
    letters = ("A", "B")
    for ax, letter, (title, wui_csv, wild_csv) in zip(axes, letters, PANELS):
        val_u = values(VAR / wui_csv)
        val_w = values(VAR / wild_csv)
        mu_u, sig_u = float(np.mean(val_u)), float(np.std(val_u, ddof=1))
        mu_w, sig_w = float(np.mean(val_w)), float(np.std(val_w, ddof=1))
        diff = mu_u - mu_w
        _, p = sps.ttest_ind(val_u, val_w, equal_var=False)

        x = np.linspace(*XLIM, 200)
        h_w = ax.hist(
            val_w,
            bins=N_BINS,
            range=XLIM,
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
            bins=N_BINS,
            range=XLIM,
            density=True,
            color=C_WUI,
            edgecolor=(0.9, 0.9, 0.9),
            alpha=0.4,
            zorder=1,
        )
        ax.plot(x, sps.norm.pdf(x, mu_u, max(sig_u, 1e-9)), color=C_WUI, lw=2.0, zorder=3)
        ax.axvline(mu_u, color=C_WUI, ls=":", lw=2.0, zorder=2)

        if legend_handles is None:
            legend_handles = (h_u[2][0], h_w[2][0])

        ax.set_xlim(*XLIM)
        ax.set_xlabel("Elevation (m)", fontsize=10)
        ax.set_ylabel("Density", fontsize=10)
        ax.tick_params(axis="both", direction="out", labelsize=9, length=3.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_linewidth(0.8)
            ax.spines[sp].set_visible(True)
        ax.set_title(
            f"({letter}) {title}   Diff:{diff:+.1f}  ({p_string(p)})",
            loc="left",
            fontsize=10,
            fontweight="normal",
            pad=6,
        )
        print(
            f"({letter}) {title}: nU/nW={len(val_u)}/{len(val_w)}  "
            f"mean {mu_u:.1f} vs {mu_w:.1f}  Δ={diff:+.1f}  p={p:.3g}"
        )

    axes[0].legend(
        legend_handles,
        ["WUI", "Wildland"],
        loc="upper right",
        frameon=False,
        fontsize=9,
        handlelength=1.2,
        borderaxespad=0.2,
    )

    fig.savefig(OUT_PNG, dpi=300, facecolor="white")
    fig.savefig(OUT_PDF, facecolor="white")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)


if __name__ == "__main__":
    main()

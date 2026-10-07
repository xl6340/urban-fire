#!/usr/bin/env python3
"""Fig 2B style: beta vs landscape VPDmax using fire-season (May–Oct) means."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

BASE = Path(__file__).resolve().parents[2]
OUT_PNG = BASE / "Fig/Fig2B_fireseason.png"
OUT_PDF = BASE / "Fig/Fig2B_fireseason.pdf"

COLOR_URBAN = np.array([216, 118, 89]) / 255
COLOR_WILD = np.array([41, 157, 143]) / 255


def main():
    beta = pd.read_excel(BASE / "dataFig/vpd/betaVPD.xlsx")
    fs = pd.read_csv(BASE / "dataPrc/vpdmax_landscape_allyear_vs_fireseason.csv")
    fs = fs[fs["season"] == "fire_season_MayOct"]
    df = beta.merge(
        fs[["Period", "vpd_wui_kPa", "vpd_wild_kPa"]],
        on="Period",
        how="inner",
        suffixes=("_pub", "_fs"),
    )

    vpd_wui = df["vpd_wui_kPa_fs"].to_numpy()
    vpd_wild = df["vpd_wild_kPa_fs"].to_numpy()
    beta_u = df["beta_urban"].to_numpy()
    err_u = df["beta_err.1"].to_numpy()
    beta_w = df["beta_wild"].to_numpy()
    err_w = df["beta_err.2"].to_numpy()
    periods = df["Period"].astype(str).tolist()

    fig, ax = plt.subplots(figsize=(4.2, 3.6), dpi=150, facecolor="white")

    ax.errorbar(
        vpd_wui, beta_u, yerr=err_u, fmt="o", color=COLOR_URBAN, mfc=COLOR_URBAN,
        markersize=6, linewidth=1, capsize=0, linestyle="none", zorder=3,
    )
    m = np.isfinite(vpd_wui) & np.isfinite(beta_u)
    r_u, p_u = pearsonr(vpd_wui[m], beta_u[m])
    coef = np.polyfit(vpd_wui[m], beta_u[m], 1)
    xfit = np.linspace(vpd_wui[m].min(), vpd_wui[m].max(), 100)
    ax.plot(xfit, np.polyval(coef, xfit), "-", color=COLOR_URBAN, lw=1.5, zorder=2)
    for x, y, lab in zip(vpd_wui, beta_u, periods):
        if np.isfinite(x) and np.isfinite(y):
            ax.text(x - 0.02, y - 0.03, lab, color="0.55", fontsize=7.5, ha="center", va="bottom")

    ax.errorbar(
        vpd_wild, beta_w, yerr=err_w, fmt="^", color=COLOR_WILD, mfc=COLOR_WILD,
        markersize=6, linewidth=1, capsize=0, linestyle="none", zorder=3,
    )
    m = np.isfinite(vpd_wild) & np.isfinite(beta_w)
    r_w, p_w = pearsonr(vpd_wild[m], beta_w[m])
    coef = np.polyfit(vpd_wild[m], beta_w[m], 1)
    xfit = np.linspace(vpd_wild[m].min(), vpd_wild[m].max(), 100)
    ax.plot(xfit, np.polyval(coef, xfit), "-", color=COLOR_WILD, lw=1.5, zorder=2)
    for x, y, lab in zip(vpd_wild, beta_w, periods):
        if np.isfinite(x) and np.isfinite(y):
            ax.text(x + 0.02, y - 0.01, lab, color="0.55", fontsize=7.5, ha="center", va="top")

    ax.text(0.05, 0.92, f"WUI: r={r_u:.2f}, p={p_u:.3f}", transform=ax.transAxes,
            color=COLOR_URBAN, fontsize=11, va="top")
    ax.text(0.05, 0.84, f"Wildland: r={r_w:.2f}, p={p_w:.3f}", transform=ax.transAxes,
            color=COLOR_WILD, fontsize=11, va="top")

    ax.set_xlabel(r"$VPD_{max}$ (kPa)", fontsize=12)
    ax.set_ylabel(r"$\beta$ value", fontsize=12)
    xmin = min(np.nanmin(vpd_wui), np.nanmin(vpd_wild)) - 0.05
    xmax = max(np.nanmax(vpd_wui), np.nanmax(vpd_wild)) + 0.05
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-1.7, -1.05)
    ax.set_yticks(np.arange(-1.7, -1.0, 0.2))
    ax.tick_params(direction="in", length=4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.text(-0.12, 1.02, "(B)", transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
    print("Saved", OUT_PNG)
    print(f"WUI r={r_u:.3f} p={p_u:.4f}; Wildland r={r_w:.3f} p={p_w:.4f}")


if __name__ == "__main__":
    main()

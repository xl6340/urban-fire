#!/usr/bin/env python3
"""
Fig. 4 composite: panels A–C in the layout of code/Figs/Fig4.m.

Panel A uses CalFire fires ≥ 1 km² (dataFig/ba/CalFire_ge1km2_annual.csv).
Trend lines and slope labels are drawn only when p < 0.01, as in Fig4.m.
Panels B and C use the original beta and facility tables.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(__file__).resolve().parents[2]
ANNUAL = BASE / "dataFig/ba/CalFire_ge1km2_annual.csv"
BETA = {
    "WUI": BASE / "dataFig/beta/2Igns-Urban-edge.csv",
    "Wildland": BASE / "dataFig/beta/2Igns-Wildland.csv",
}
FACILITY = BASE / "dataFig/facility.xlsx"
OUT_PNG = BASE / "Fig/Fig4.png"
OUT_PDF = BASE / "Fig/Fig4.pdf"

C_HUMAN = "#fdd85d"
C_NAT = "#99d6ea"
C_BAR = "#f2e9e4"


def fit_trend(x, y):
    res = stats.linregress(np.asarray(x, float), np.asarray(y, float))
    xfit = np.linspace(float(np.min(x)), float(np.max(x)), 100)
    return res.slope, res.pvalue, xfit, res.intercept + res.slope * xfit


def style_ax(ax, box=True, rotate_x=False):
    ax.tick_params(direction="in", labelsize=9, length=3.0, width=0.6, pad=1.5)
    if rotate_x:
        ax.tick_params(axis="x", rotation=30)
    for sp in ax.spines.values():
        sp.set_linewidth(0.6)
        sp.set_color("k")
    if not box:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)


def draw_A(axes, annual):
    x = annual["year"].to_numpy()
    specs = [
        (0, "WUI", "WUI_n_Human", "WUI_n_Natural", "WUI_area_Human", "WUI_area_Natural", "WUI fire"),
        (1, "Wildland", "Wildland_n_Human", "Wildland_n_Natural",
         "Wildland_area_Human", "Wildland_area_Natural", "Wildland fire"),
    ]
    series = [("Human", C_HUMAN), ("Natural", C_NAT)]
    for col, ft, n_h, n_n, a_h, a_n, title in specs:
        ax_n, ax_a = axes[0, col], axes[1, col]
        counts = {"Human": annual[n_h].to_numpy(), "Natural": annual[n_n].to_numpy()}
        areas = {
            "Human": annual[a_h].to_numpy() / 1000.0,
            "Natural": annual[a_n].to_numpy() / 1000.0,
        }
        for name, color in series:
            ax_n.plot(
                x, counts[name], "-o", color=color, markerfacecolor=color,
                markeredgecolor=color,                 ms=4, lw=0.4, label=f"{name} ignited",
                zorder=3, clip_on=True,
            )
            slope, p, xfit, yfit = fit_trend(x, counts[name])
            if p < 0.01:
                ax_n.plot(xfit, yfit, "-", color="0.7", lw=1.0, zorder=1)
                ax_n.text(
                    0.18, 0.40 if col == 0 else 0.58, f"s = {slope:.2f}",
                    transform=ax_n.transAxes, fontsize=9, color="k",
                )
            print(f"{ft:9s} count {name:8s}  s={slope:7.3f}  p={p:.3g}")

            ax_a.plot(
                x, areas[name], "-o", color=color, markerfacecolor=color,
                markeredgecolor=color,                 ms=4, lw=0.4, label=f"{name} ignited",
                zorder=3,
            )
            slope, p, xfit, yfit = fit_trend(x, areas[name])
            if p < 0.01:
                ax_a.plot(xfit, yfit, "-", color="0.7", lw=1.0, zorder=1)
                ax_a.text(
                    0.18, 0.40 if col == 0 else 0.58, f"s = {slope * 1000:.1f}",
                    transform=ax_a.transAxes, fontsize=9, color="k",
                )
            print(f"{ft:9s} area  {name:8s}  s={slope * 1000:7.2f}  p={p:.3g}")

        ax_n.set_title(title, fontsize=9, fontweight="normal", pad=2)
        ax_n.set_ylim(0, 400)
        ax_n.set_yticks([0, 100, 200, 300, 400])
        ax_a.set_ylim(0, 12)
        ax_a.set_yticks([0, 4, 8, 12])
        for ax in (ax_n, ax_a):
            ax.set_xlim(1990, 2025)
            ax.set_xticks(range(1990, 2026, 5))
            style_ax(ax, box=True, rotate_x=True)
            for lab in ax.get_xticklabels():
                lab.set_ha("right")
                lab.set_rotation_mode("anchor")
        ax_n.tick_params(labelbottom=False)
        if col == 1:
            ax_n.tick_params(labelleft=False)
            ax_a.tick_params(labelleft=False)
        else:
            ax_n.set_ylabel("Fire number (#)", fontsize=9)
            ax_a.set_ylabel(r"Burned area ($10^3$ km$^2$)", fontsize=9)

    axes[1, 1].legend(
        frameon=False, fontsize=9, loc="upper right",
        borderpad=0.1, handlelength=1.6, handletextpad=0.4, labelspacing=0.15,
    )


def draw_B(ax):
    order = [("WUI", 2), ("Wildland", 1)]
    for name, y in order:
        tab = pd.read_csv(BETA[name])
        betas = tab["beta"].to_numpy()
        errs = tab["betaErr"].to_numpy()
        ax.plot(betas, [y, y], "-", color="0.7", lw=1.0, zorder=1)
        for beta, err, color in zip(betas, errs, (C_HUMAN, C_NAT)):
            ax.errorbar(
                beta, y, xerr=err, fmt="o", color=color, markerfacecolor=color,
                markeredgecolor=color, ms=4, lw=1.5, capsize=6,
                markeredgewidth=0.6, zorder=3,
            )
    ax.set_yticks([1, 2])
    ax.set_yticklabels(["Wildland", "WUI"])
    ax.set_ylim(0.5, 2.5)
    ax.set_xlim(-1.8, -1.0)
    ax.set_xticks(np.round(np.arange(-1.8, -0.99, 0.2), 1))
    ax.set_xlabel(r"$\beta$ value", fontsize=11)
    style_ax(ax, box=True)


def draw_C(ax):
    fac = pd.read_excel(FACILITY)
    ax.barh(
        np.arange(len(fac)), fac["density"].to_numpy(),
        color=C_BAR, height=0.8, zorder=2, edgecolor=C_BAR,
    )
    ax.set_yticks(np.arange(len(fac)))
    ax.set_yticklabels(fac["landscape"])
    ax.set_ylim(len(fac) - 0.5, -0.5)
    ax.set_xlim(0, 4)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.set_xlabel(r"Facility density (#/100-km$^2$)", fontsize=9)
    style_ax(ax, box=False)


def main():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })
    annual = pd.read_csv(ANNUAL)

    # Original assembled page is 562 x 270 pt. Panel A is the left 2x2;
    # B and C share the right column and line up with A's two rows.
    fig = plt.figure(figsize=(562 / 72, 270 / 72), facecolor="w")
    outer = fig.add_gridspec(
        1, 2, width_ratios=[2.15, 1.0], wspace=0.16,
        left=0.09, right=0.985, top=0.93, bottom=0.16,
    )
    gs_a = outer[0].subgridspec(2, 2, wspace=0.22, hspace=0.32)
    gs_r = outer[1].subgridspec(2, 1, hspace=0.55)
    axes = np.empty((2, 2), dtype=object)
    for r in range(2):
        for c in range(2):
            axes[r, c] = fig.add_subplot(gs_a[r, c])
    ax_b = fig.add_subplot(gs_r[0])
    ax_c = fig.add_subplot(gs_r[1])

    draw_A(axes, annual)
    draw_B(ax_b)
    draw_C(ax_c)

    axes[0, 0].text(
        -0.34, 1.02, "(A)", transform=axes[0, 0].transAxes,
        fontsize=11, fontweight="bold", va="bottom", ha="left", clip_on=False,
    )
    ax_b.text(
        -0.22, 1.02, "(B)", transform=ax_b.transAxes,
        fontsize=11, fontweight="bold", va="bottom", ha="left", clip_on=False,
    )
    ax_c.text(
        -0.22, 1.06, "(C)", transform=ax_c.transAxes,
        fontsize=11, fontweight="bold", va="bottom", ha="left", clip_on=False,
    )

    fig.savefig(OUT_PNG, dpi=300, facecolor="w")
    fig.savefig(OUT_PDF, facecolor="w")
    print("Saved", OUT_PNG)


if __name__ == "__main__":
    main()

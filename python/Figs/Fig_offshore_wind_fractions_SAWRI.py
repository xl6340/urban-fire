#!/usr/bin/env python3
"""
Combined figure: both panels restricted to R1D-SAWRI era (≤2018).
(a) SoCal Santa Ana (SAWRI>1)
(b) NorCal Diablo (Smith et al. 2018 NoBA RAWS)
Legend: Fire size < p90 / ≥ p90.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

BASE = Path(__file__).resolve().parents[2]
SOCAL_FLAGS = BASE / "dataPrc/CalFire_SoCal_SAW_extended_flags.csv"
NORCAL_FLAGS = BASE / "dataPrc/CalFire_NorCal_Diablo_flags.csv"
OUT_PNG = BASE / "Fig/Fig_offshore_wind_fractions_SAWRI.png"
OUT_PDF = BASE / "Fig/Fig_offshore_wind_fractions_SAWRI.pdf"

C_WUI = np.array([216, 118, 89]) / 255
C_WILD = np.array([41, 157, 143]) / 255
ERA_END = "2018-12-31"


def era_summary(flags_path, wind_col, end=ERA_END, require_has_wx=False):
    f = pd.read_csv(flags_path, parse_dates=["IDate"])
    pre = f[f["IDate"] <= end].copy()
    if require_has_wx:
        pre = pre[pre["has_wx"]].copy()
    rows = []
    for ft in ["WUI", "Wildland"]:
        g = pre[pre["FireType"] == ft]
        p90 = g["size"].quantile(0.9)
        large = g[g["size"] >= p90]
        below = g[g["size"] < p90]
        rows.append(
            {
                "FireType": ft,
                "n_below": len(below),
                "n_wind_below": int(below[wind_col].sum()),
                "frac_below": below[wind_col].mean(),
                "n_large": len(large),
                "n_large_wind": int(large[wind_col].sum()),
                "frac_large": large[wind_col].mean(),
            }
        )
    return pd.DataFrame(rows)


def _panel(ax, df, ylabel):
    x = np.arange(2)
    w = 0.34
    df = df.set_index("FireType").loc[["WUI", "Wildland"]].reset_index()
    colors = [C_WUI, C_WILD]

    bars1 = ax.bar(
        x - w / 2,
        df["frac_below"] * 100,
        width=w,
        color=colors,
        edgecolor="0.25",
        lw=0.6,
        zorder=2,
    )
    bars2 = ax.bar(
        x + w / 2,
        df["frac_large"] * 100,
        width=w,
        color=colors,
        edgecolor="0.25",
        lw=0.6,
        hatch="///",
        alpha=0.75,
        zorder=2,
    )

    ymax = max(df["frac_below"].max(), df["frac_large"].max()) * 100
    y_top = max(22, ymax * 1.38)
    ax.set_ylim(0, y_top)
    pad = y_top * 0.02

    for b, n_s, n in zip(bars1, df["n_wind_below"], df["n_below"]):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + pad,
            f"{b.get_height():.1f}%\n({int(n_s)}/{int(n)})",
            ha="center",
            va="bottom",
            fontsize=8,
            color="0.2",
        )
    for b, n_s, n in zip(bars2, df["n_large_wind"], df["n_large"]):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + pad,
            f"{b.get_height():.1f}%\n({int(n_s)}/{int(n)})",
            ha="center",
            va="bottom",
            fontsize=8,
            color="0.2",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(["WUI", "Wildland"], fontsize=11)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.tick_params(direction="out", labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    socal = era_summary(SOCAL_FLAGS, "SAW_official")
    norcal = era_summary(NORCAL_FLAGS, "is_Diablo", require_has_wx=True)

    for name, df in [("SoCal ≤2018", socal), ("NorCal ≤2018", norcal)]:
        print(f"{name}:")
        for _, r in df.iterrows():
            print(
                f"  {r.FireType}: <p90 {100*r.frac_below:.1f}% "
                f"({int(r.n_wind_below)}/{int(r.n_below)}) | "
                f"≥p90 {100*r.frac_large:.1f}% "
                f"({int(r.n_large_wind)}/{int(r.n_large)})"
            )

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), facecolor="white")
    _panel(axes[0], socal, "Fraction of ignitions on\nSanta Ana days (%)")
    _panel(axes[1], norcal, "Fraction of ignitions on\nDiablo-like days (%)")

    for ax, letter, region in zip(
        axes, ["a", "b"], ["Southern California", "Northern California"]
    ):
        ax.text(
            0.02,
            0.98,
            letter,
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="top",
            ha="left",
        )
        ax.set_title(region, fontsize=11, pad=4)

    handles = [
        Patch(facecolor=C_WUI, edgecolor="0.25", label="Fire size < p90"),
        Patch(
            facecolor=C_WUI,
            edgecolor="0.25",
            hatch="///",
            alpha=0.75,
            label="Fire size ≥ p90",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=2,
        frameon=False,
        fontsize=9,
        bbox_to_anchor=(0.5, 1.02),
    )

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)


if __name__ == "__main__":
    main()

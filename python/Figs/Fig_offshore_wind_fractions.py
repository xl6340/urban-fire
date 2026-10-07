#!/usr/bin/env python3
"""
Southern California Santa Ana ignition fractions.
Solid bars = fires < size p90; hatched = fires ≥ p90. No title/footnote.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

BASE = Path(__file__).resolve().parents[2]
SOCAL = BASE / "dataPrc/CalFire_SoCal_SAW_extended_summary.csv"
OUT_PNG = BASE / "Fig/Fig_offshore_wind_fractions.png"
OUT_PDF = BASE / "Fig/Fig_offshore_wind_fractions.pdf"

C_WUI = np.array([216, 118, 89]) / 255
C_WILD = np.array([41, 157, 143]) / 255


def add_below_p90(df, n_wind_col, n_large_wind_col):
    """Derive <p90 counts/fractions from all vs ≥p90 columns."""
    out = df.copy()
    out["n_below"] = out["n"] - out["n_large"]
    out["n_wind_below"] = out[n_wind_col] - out[n_large_wind_col]
    out["frac_below"] = out["n_wind_below"] / out["n_below"]
    return out


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
    socal = pd.read_csv(SOCAL)
    socal = add_below_p90(socal, "n_SAW", "n_large_SAW")
    socal["frac_large"] = socal["frac_large"]
    socal["n_large_wind"] = socal["n_large_SAW"]

    print("SoCal <p90 / ≥p90:")
    for _, r in socal.iterrows():
        print(
            f"  {r.FireType}: <p90 {100*r.frac_below:.1f}% ({int(r.n_wind_below)}/{int(r.n_below)}) | "
            f"≥p90 {100*r.frac_large:.1f}% ({int(r.n_large_wind)}/{int(r.n_large)})"
        )

    fig, ax = plt.subplots(figsize=(3.8, 3.4), facecolor="white")
    _panel(ax, socal, "Fraction of ignitions on\nSanta Ana days (%)")

    handles = [
        Patch(facecolor="none", edgecolor="0.25", label="Fire size < p90"),
        Patch(
            facecolor="none",
            edgecolor="0.25",
            hatch="///",
            label="Fire size ≥ p90",
        ),
    ]
    ax.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(1.0, 0.88),
        frameon=False,
        fontsize=8.5,
    )
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)


if __name__ == "__main__":
    main()

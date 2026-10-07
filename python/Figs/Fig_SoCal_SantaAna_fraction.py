#!/usr/bin/env python3
"""
Fraction of SoCal CalFire ignitions on Santa Ana days (R1D-SAWRI > 1).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
FLAGS = BASE / "dataPrc/CalFire_SoCal_SAWRI_flags.csv"
OUT_PNG = BASE / "Fig/Fig_SoCal_SantaAna_fraction.png"
OUT_PDF = BASE / "Fig/Fig_SoCal_SantaAna_fraction.pdf"

C_WUI = np.array([216, 118, 89]) / 255
C_WILD = np.array([41, 157, 143]) / 255


def main():
    f = pd.read_csv(FLAGS, parse_dates=["IDate"])
    sub = f[f["in_SAWRI_period"]].copy()

    rows = []
    for ft, color in [("WUI", C_WUI), ("Wildland", C_WILD)]:
        g = sub[sub["FireType"] == ft]
        p90 = g["size"].quantile(0.9)
        large = g[g["size"] >= p90]
        rows.append(
            {
                "FireType": ft,
                "color": color,
                "frac_all": g["SAW_gt1"].mean(),
                "n_all": len(g),
                "n_saw_all": int(g["SAW_gt1"].sum()),
                "frac_large": large["SAW_gt1"].mean(),
                "n_large": len(large),
                "n_saw_large": int(large["SAW_gt1"].sum()),
            }
        )
    df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(5.2, 4.0), facecolor="white")
    x = np.arange(2)
    w = 0.34
    bars1 = ax.bar(
        x - w / 2,
        df["frac_all"] * 100,
        width=w,
        color=[df.loc[0, "color"], df.loc[1, "color"]],
        edgecolor="0.25",
        linewidth=0.6,
        label="All SoCal fires",
        alpha=0.95,
    )
    # hatch large-fire bars slightly differently via lighter face + hatch
    bars2 = ax.bar(
        x + w / 2,
        df["frac_large"] * 100,
        width=w,
        color=[df.loc[0, "color"], df.loc[1, "color"]],
        edgecolor="0.25",
        linewidth=0.6,
        hatch="///",
        label="Large fires (≥ size p90)",
        alpha=0.75,
    )

    # value labels
    for b, n_saw, n in zip(
        bars1,
        df["n_saw_all"],
        df["n_all"],
    ):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 1.2,
            f"{b.get_height():.1f}%\n({n_saw}/{n})",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color="0.2",
        )
    for b, n_saw, n in zip(
        bars2,
        df["n_saw_large"],
        df["n_large"],
    ):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 1.2,
            f"{b.get_height():.1f}%\n({n_saw}/{n})",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color="0.2",
        )

    ax.set_xticks(x)
    ax.set_xticklabels(["WUI", "Wildland"], fontsize=12)
    ax.set_ylabel("Fraction of ignitions on Santa Ana days (%)", fontsize=11)
    ax.set_ylim(0, 58)
    ax.set_title(
        "Southern California: Santa Ana ignition-day fraction\n"
        r"(R1D-SAWRI $> 1$; 1990–2018)",
        fontsize=12,
        pad=8,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out")
    ax.legend(frameon=False, fontsize=9, loc="upper right")

    # note under plot
    ax.text(
        0.0,
        -0.18,
        "Orange = WUI, teal = Wildland. Hatched bars = fires ≥ within-type size 90th percentile.",
        transform=ax.transAxes,
        fontsize=8,
        color="0.35",
        ha="left",
        va="top",
    )

    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)
    print(df[["FireType", "frac_all", "frac_large", "n_all", "n_large"]].to_string(index=False))


if __name__ == "__main__":
    main()

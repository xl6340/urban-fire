#!/usr/bin/env python3
"""
NorCal Diablo-like ignition-day fractions for CalFire WUI vs Wildland.

A fire is on a Diablo day if its IDate is in the Smith et al. (2018) NoBA
RAWS catalog (code/diablo_raws_catalog.py): any of 6 North-of-Bay stations,
wind 315–135°, vs > 11.17 m/s, RH < 30%, ≥3 consecutive hours.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[2]
import sys

sys.path.insert(0, str(BASE / "code"))
from diablo_raws_catalog import MIN_HOURS, RH_MAX, VS_MIN, load_diablo_catalog

OUT_FLAGS = BASE / "dataPrc/CalFire_NorCal_Diablo_flags.csv"
OUT_SUM = BASE / "dataPrc/CalFire_NorCal_Diablo_summary.csv"
OUT_PNG = BASE / "Fig/Fig_NorCal_Diablo_fraction.png"
OUT_PDF = BASE / "Fig/Fig_NorCal_Diablo_fraction.pdf"

C_WUI = np.array([216, 118, 89]) / 255
C_WILD = np.array([41, 157, 143]) / 255


def main():
    print("Building/loading Diablo-like day catalog …")
    catalog = load_diablo_catalog()
    diablo_dates = set(catalog.index[catalog])

    g = gpd.read_file(BASE / "dataPrc/firePrmt/CalFire.shp", ignore_geometry=True)
    g = g[g["NorCal"] == 1].copy()
    g["IDate"] = pd.to_datetime(g["IDate"]).dt.normalize()
    g["has_wx"] = g["IDate"].notna() & g["IDate"].isin(catalog.index)
    g["is_Diablo"] = g["has_wx"] & g["IDate"].isin(diablo_dates)

    size_col = "size_km2" if "size_km2" in g.columns else "size"
    flags = pd.DataFrame(
        {
            "FireType": g["FireType"].values,
            "fid": g["fid"].values,
            "IDate": g["IDate"].values,
            "year": g["year"].values,
            "size": g[size_col].values,
            "vs": g["vs"].values if "vs" in g.columns else np.nan,
            "vd": g["vd"].values if "vd" in g.columns else np.nan,
            "has_wx": g["has_wx"].values,
            "is_Diablo": g["is_Diablo"].values,
        }
    )
    flags.to_csv(OUT_FLAGS, index=False)
    print("Saved", OUT_FLAGS)

    sub = flags[flags["has_wx"]].copy()
    rows = []
    for ft in ["WUI", "Wildland"]:
        gg = sub[sub["FireType"] == ft]
        p90 = gg["size"].quantile(0.9)
        large = gg[gg["size"] >= p90]
        rows.append(
            {
                "FireType": ft,
                "n": len(gg),
                "n_Diablo": int(gg["is_Diablo"].sum()),
                "frac_Diablo": gg["is_Diablo"].mean(),
                "n_large": len(large),
                "n_large_Diablo": int(large["is_Diablo"].sum()),
                "frac_large_Diablo": large["is_Diablo"].mean(),
                "mean_size_Diablo": gg.loc[gg["is_Diablo"], "size"].mean(),
                "mean_size_non": gg.loc[~gg["is_Diablo"], "size"].mean(),
                "vs_min": VS_MIN,
                "rh_max": RH_MAX,
                "min_hours": MIN_HOURS,
            }
        )
        print(
            f"{ft}: Diablo-like {int(gg.is_Diablo.sum())}/{len(gg)} "
            f"= {100*gg.is_Diablo.mean():.1f}% | large {100*large.is_Diablo.mean():.1f}%"
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_SUM, index=False)
    print("Saved", OUT_SUM)

    fig, ax = plt.subplots(figsize=(5.2, 4.0), facecolor="white")
    x = np.arange(2)
    w = 0.34
    colors = [C_WUI, C_WILD]
    bars1 = ax.bar(
        x - w / 2,
        summary["frac_Diablo"] * 100,
        width=w,
        color=colors,
        edgecolor="0.25",
        linewidth=0.6,
        label="All NorCal fires",
    )
    bars2 = ax.bar(
        x + w / 2,
        summary["frac_large_Diablo"] * 100,
        width=w,
        color=colors,
        edgecolor="0.25",
        linewidth=0.6,
        hatch="///",
        alpha=0.75,
        label="Large fires (≥ size p90)",
    )
    for b, n_d, n in zip(bars1, summary["n_Diablo"], summary["n"]):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 0.8,
            f"{b.get_height():.1f}%\n({n_d}/{n})",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color="0.2",
        )
    for b, n_d, n in zip(bars2, summary["n_large_Diablo"], summary["n_large"]):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 0.8,
            f"{b.get_height():.1f}%\n({n_d}/{n})",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color="0.2",
        )

    ymax = max(summary["frac_Diablo"].max(), summary["frac_large_Diablo"].max()) * 100
    ax.set_ylim(0, max(25, ymax * 1.35))
    ax.set_xticks(x)
    ax.set_xticklabels(["WUI", "Wildland"], fontsize=12)
    ax.set_ylabel("Fraction of ignitions on Diablo-like days (%)", fontsize=11)
    ax.set_title(
        "Northern California: Diablo-like offshore ignition-day fraction\n"
        r"(Smith et al. 2018 NoBA RAWS)",
        fontsize=10,
        pad=8,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)


if __name__ == "__main__":
    main()

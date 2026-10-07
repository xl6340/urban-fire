#!/usr/bin/env python3
"""
Fig. 5A style: CalFire annual fire number and burned area, 1990–2025.
WUI vs Wildland × Human vs Natural. Fires with size < 1 km² excluded.

Layout and styling match code/Figs/Fig4.m (Fig 5A).
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(__file__).resolve().parents[2]
FIRE = BASE / "dataPrc/firePrmt/CalFire.shp"
OUT_PNG = BASE / "Fig/Fig5A_CalFire_ge1km2.png"
OUT_PDF = BASE / "Fig/Fig5A_CalFire_ge1km2.pdf"
OUT_CSV = BASE / "dataFig/ba/CalFire_ge1km2_annual.csv"

C_HUMAN = "#fdd85d"
C_NAT = "#99d6ea"
YEARS = np.arange(1990, 2025)
FIRE_TYPES = ["WUI", "Wildland"]
FT_LABEL = {"WUI": "WUI fire", "Wildland": "Wildland fire"}
SIZE_MIN = 1.0


def annual_table(g: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for year in YEARS:
        gy = g[g["year"] == year]
        rows.append(
            {
                "year": int(year),
                "Human": int((gy["Ignition"] == "Human").sum()),
                "Natural": int((gy["Ignition"] == "Natural").sum()),
                "area_Human": float(gy.loc[gy["Ignition"] == "Human", "size_km2"].sum()),
                "area_Natural": float(gy.loc[gy["Ignition"] == "Natural", "size_km2"].sum()),
            }
        )
    return pd.DataFrame(rows)


def fit_trend(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    res = stats.linregress(x, y)
    xfit = np.linspace(x.min(), x.max(), 100)
    yfit = res.intercept + res.slope * xfit
    return res.slope, res.pvalue, xfit, yfit


def main():
    g = gpd.read_file(FIRE, ignore_geometry=True)
    size = g["size_km2"] if "size_km2" in g.columns else g["size"]
    g = g[
        size.notna()
        & (size >= SIZE_MIN)
        & g["year"].between(YEARS[0], YEARS[-1])
        & g["Ignition"].isin(["Human", "Natural"])
        & g["FireType"].isin(FIRE_TYPES)
    ].copy()
    g["size_km2"] = size.loc[g.index]
    print(
        f"CalFire size ≥ {SIZE_MIN} km², {int(g['year'].min())}–{int(g['year'].max())}: "
        f"n={len(g)}  min={g['size_km2'].min():.3f}"
    )
    print(pd.crosstab(g["FireType"], g["Ignition"], margins=True))

    tables = {ft: annual_table(g[g["FireType"] == ft]) for ft in FIRE_TYPES}
    out = tables["WUI"].rename(
        columns={
            "Human": "WUI_n_Human",
            "Natural": "WUI_n_Natural",
            "area_Human": "WUI_area_Human",
            "area_Natural": "WUI_area_Natural",
        }
    ).merge(
        tables["Wildland"].rename(
            columns={
                "Human": "Wildland_n_Human",
                "Natural": "Wildland_n_Natural",
                "area_Human": "Wildland_area_Human",
                "area_Natural": "Wildland_area_Natural",
            }
        ),
        on="year",
    )
    out.to_csv(OUT_CSV, index=False)
    print("Saved", OUT_CSV)

    # code/Figs/Fig4.m uses MATLAB defaults: Box off (no top/right frame),
    # TickDir in, short ticks, axes LineWidth 0.5, MarkerSize 4, fontsize 9.
    # Fire-number ylim in that script is [0 400] for the unfiltered record.
    # Counts in the ≥1 km² sample stay below 100, so the axis is 0–120.
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "axes.unicode_minus": False,
        "pdf.fonttype": 42,
    })
    fig, axes = plt.subplots(2, 2, figsize=(500 / 96, 350 / 96), facecolor="w")
    fig.subplots_adjust(wspace=0.12, hspace=0.22, left=0.12, right=0.98, top=0.90, bottom=0.16)

    for col, ft in enumerate(FIRE_TYPES):
        t = tables[ft]
        x = t["year"].to_numpy()
        series = [
            ("Human", t["Human"].to_numpy(), t["area_Human"].to_numpy() / 1000.0, C_HUMAN),
            ("Natural", t["Natural"].to_numpy(), t["area_Natural"].to_numpy() / 1000.0, C_NAT),
        ]
        ax_n, ax_a = axes[0, col], axes[1, col]
        txt_y = 0.40 if col == 0 else 0.60

        for name, yn, ya, color in series:
            ax_n.plot(
                x, yn, "-o", color=color, markerfacecolor=color, markeredgecolor=color,
                ms=4, lw=0.4, label=f"{name} ignited",
            )
            slope, p, xfit, yfit = fit_trend(x, yn)
            if p < 0.01 and name == "Human":
                ax_n.plot(xfit, yfit, "-", color=(0.7, 0.7, 0.7), lw=0.8, zorder=1)
                ax_n.text(
                    0.20, txt_y, f"s = {slope:.2f}",
                    transform=ax_n.transAxes, fontsize=10, color="k",
                )
            print(f"{ft:9s} count {name:8s}  s={slope:.3f}  p={p:.3g}")

            ax_a.plot(
                x, ya, "-o", color=color, markerfacecolor=color, markeredgecolor=color,
                ms=4, lw=0.4, label=f"{name} ignited",
            )
            slope, p, xfit, yfit = fit_trend(x, ya)
            if p < 0.01 and name == "Human":
                ax_a.plot(xfit, yfit, "-", color=(0.7, 0.7, 0.7), lw=0.8, zorder=1)
                ax_a.text(
                    0.20, txt_y, f"s = {slope * 1000:.1f}",
                    transform=ax_a.transAxes, fontsize=10, color="k",
                )
            print(f"{ft:9s} area  {name:8s}  s={slope*1000:.2f} km²/yr  p={p:.3g}")

        ax_n.set_title(FT_LABEL[ft], fontsize=9, fontweight="normal", pad=2)
        ax_n.set_ylim(0, 120)
        ax_n.set_yticks([0, 40, 80, 120])
        ax_a.set_ylim(0, 12)
        ax_a.set_yticks([0, 4, 8, 12])
        for ax in (ax_n, ax_a):
            ax.set_xlim(1990, 2025)
            ax.set_xticks(range(1990, 2026, 5))
            ax.tick_params(direction="in", labelsize=9, length=2.0, width=0.4, pad=1.5)
            ax.minorticks_off()
            ax.tick_params(axis="x", rotation=30)
            for lab in ax.get_xticklabels():
                lab.set_ha("right")
                lab.set_rotation_mode("anchor")
            for sp in ax.spines.values():
                sp.set_linewidth(0.4)
                sp.set_color("k")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        ax_n.tick_params(labelbottom=False)
        if col == 1:
            ax_n.tick_params(labelleft=False)
            ax_a.tick_params(labelleft=False)
        else:
            ax_n.set_ylabel("Fire number (#)", fontsize=9)
            ax_a.set_ylabel(r"Burned area (10$^3$ km$^2$)", fontsize=9)

    axes[1, 1].legend(frameon=False, fontsize=9, loc="upper right")
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="w")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="w")
    print("Saved", OUT_PNG)


if __name__ == "__main__":
    main()

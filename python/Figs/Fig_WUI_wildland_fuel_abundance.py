#!/usr/bin/env python3
"""
Horizontal z-score bars: fuel/biomass abundance, WUI vs wildland.

Panels: Forest | Shrub & Grassland
Order (top→bottom): MODIS NDVI, ESA CCI tree woody AGB,
                    RAP shrub & herb NPP, LANDFIRE shrub & grass fraction
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

BASE = Path(__file__).resolve().parents[2]
BY_FIRE = BASE / "dataPrc/WUI_vs_wildland_fuel_abundance_by_fire.csv"
OUT_PNG = BASE / "Fig/Fig_WUI_wildland_fuel_abundance.png"
OUT_PDF = BASE / "Fig/Fig_WUI_wildland_fuel_abundance.pdf"
OUT_CSV = BASE / "dataPrc/WUI_vs_wildland_fuel_abundance_fig_stats.csv"

POS_C = "#c45c26"
NEG_C = "#2a6f7a"

# Fixed order (top → bottom); shrub + grass merged
VAR_SPECS = [
    ("ndvi", "MODIS NDVI"),
    ("cci_agb", "ESA CCI tree woody AGB"),
    ("rap_shr_herb_npp", "RAP shrub & herb NPP"),
    ("lf_gs", "LANDFIRE shrub & grass fraction"),
]


def zscore_wui_vs_wild(u: np.ndarray, w: np.ndarray):
    u = u[np.isfinite(u)]
    w = w[np.isfinite(w)]
    if len(u) < 3 or len(w) < 3:
        return np.nan, np.nan, np.nan, len(u), len(w)
    mu, sigma = float(np.mean(w)), float(np.std(w, ddof=1))
    if sigma == 0 or not np.isfinite(sigma):
        sigma = 1.0
    z = (u - mu) / sigma
    z_mean = float(np.mean(z))
    z_se = float(np.std(z, ddof=1) / np.sqrt(len(z)))
    _, p = sps.ttest_ind(u, w, equal_var=False, nan_policy="omit")
    return z_mean, z_se, float(p), len(u), len(w)


def stars(p):
    if not np.isfinite(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


def prepare_df(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["rap_shr_herb_npp"] = df["rap_shr_npp"].astype(float) + df["rap_herb_npp"].astype(float)
    if "lf_gs_fuel" in df.columns:
        df["lf_gs"] = df["lf_gs_fuel"].astype(float)
    else:
        df["lf_gs"] = df["lf_grass"].astype(float) + df["lf_shrub"].astype(float)
    return df


def compute_panel(df: pd.DataFrame, lc: str) -> pd.DataFrame:
    sub = df[df["lc"] == lc]
    u = sub[sub["FireType"] == "WUI"]
    w = sub[sub["FireType"] == "Wildland"]
    rows = []
    for key, label in VAR_SPECS:
        zm, se, p, nu, nw = zscore_wui_vs_wild(
            u[key].to_numpy(dtype=float), w[key].to_numpy(dtype=float)
        )
        rows.append(
            {
                "lc": lc,
                "variable": key,
                "label": label,
                "z": zm,
                "se": se,
                "p": p,
                "n_WUI": nu,
                "n_Wildland": nw,
            }
        )
    return pd.DataFrame(rows)


def plot_panel(ax, stats: pd.DataFrame, title: str, show_ylabel: bool = True):
    stats = stats.reset_index(drop=True)
    y = np.arange(len(stats))
    colors = [POS_C if z > 0 else NEG_C for z in stats["z"]]
    ax.barh(y, stats["z"], color=colors, height=0.65, edgecolor="none", zorder=2)
    ax.axvline(0, color="0.25", lw=0.9, zorder=1)
    ax.set_yticks(y)
    if show_ylabel:
        ax.set_yticklabels(stats["label"].tolist(), fontsize=9)
        ax.tick_params(axis="y", direction="out", left=True, length=3.5, labelleft=True)
    else:
        ax.set_yticklabels([])
        ax.tick_params(axis="y", left=False, labelleft=False, length=0)
    ax.invert_yaxis()
    for yi, z, p in zip(y, stats["z"], stats["p"]):
        if not np.isfinite(z):
            continue
        st = stars(p)
        ax.text(
            z + (0.04 if z >= 0 else -0.04),
            yi,
            f"{z:+.2f} {st}".rstrip(),
            va="center",
            ha="left" if z >= 0 else "right",
            fontsize=8,
            color="0.25",
        )
    ax.set_title(title, fontsize=12, fontweight="bold", pad=8)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xticks(np.arange(-0.9, 1.0, 0.3))
    ax.set_xlabel("WUI z-score vs wildland", fontsize=11)
    ax.tick_params(axis="x", direction="out")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_color("0.4")


def main():
    df = prepare_df(pd.read_csv(BY_FIRE))

    all_stats = pd.concat(
        [compute_panel(df, "Forest"), compute_panel(df, "ShrubGrass")],
        ignore_index=True,
    )
    all_stats.to_csv(OUT_CSV, index=False)
    print(all_stats.to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.4), sharex=True, facecolor="white")
    plot_panel(axes[0], all_stats[all_stats.lc == "Forest"], "Forest", show_ylabel=True)
    plot_panel(
        axes[1],
        all_stats[all_stats.lc == "ShrubGrass"],
        "Shrub & Grassland",
        show_ylabel=False,
    )
    fig.tight_layout(w_pad=1.5)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)


if __name__ == "__main__":
    main()

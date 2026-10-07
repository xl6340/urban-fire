#!/usr/bin/env python3
"""
Fig 3B: WUI vs wildland z-score dot plot (Forest | Shrub & Grassland).

Adds vegetation AGB (ESA CCI tree woody) and RAP (shrub + herb NPP)
alongside NDVI. Climate / terrain / ignition match published Fig 3B.
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sps

BASE = Path(__file__).resolve().parents[2]
FUEL_CSV = BASE / "dataPrc/WUI_vs_wildland_fuel_abundance_by_fire.csv"
OUT_PNG = BASE / "Fig/Fig3B.png"
OUT_PDF = BASE / "Fig/Fig3B.pdf"
OUT_CSV = BASE / "dataPrc/CalFire_zscore_with_AGB_RAP.csv"

CLIMATE = np.array([246, 198, 175]) / 255
VEG = np.array([181, 212, 190]) / 255
TERRAIN = np.array([175, 212, 227]) / 255
IGN = np.array([184, 185, 210]) / 255


def load_fires() -> pd.DataFrame:
    g = gpd.read_file(BASE / "dataPrc/firePrmt/CalFire.shp")
    fwi = pd.read_csv(BASE / "dataPrc/CalFire_FWI_4day.csv")
    df = g.drop(columns="geometry").merge(
        fwi[["fire_id", "FFMC_4d", "DMC_4d", "DC_4d"]],
        left_on="fid",
        right_on="fire_id",
        how="left",
    )
    df["FFMC"] = df["FFMC_4d"]
    df["DMC"] = df["DMC_4d"]
    df["DC"] = df["DC_4d"]
    df["ndviM"] = df["ndvi"]
    df["ignition"] = np.where(
        df["Ignition"] == "Human",
        1.0,
        np.where(df["Ignition"] == "Natural", 0.0, np.nan),
    )

    fuel = pd.read_csv(FUEL_CSV)[
        ["fid", "cci_agb", "rap_shr_npp", "rap_herb_npp"]
    ].drop_duplicates("fid")
    fuel["rap"] = fuel["rap_shr_npp"].astype(float) + fuel["rap_herb_npp"].astype(float)
    fuel = fuel.rename(columns={"cci_agb": "agb"})
    df = df.merge(fuel[["fid", "agb", "rap"]], on="fid", how="left")
    return df


def zscore_wui_vs_wild(u: np.ndarray, w: np.ndarray):
    u = u[np.isfinite(u)]
    w = w[np.isfinite(w)]
    if len(u) < 3 or len(w) < 3:
        return np.nan, np.nan, np.nan, len(u), len(w)
    mu, sigma = np.mean(w), np.std(w, ddof=1)
    if sigma == 0 or not np.isfinite(sigma):
        sigma = 1.0
    z = (u - mu) / sigma
    z_mean = float(np.mean(z))
    z_se = float(np.std(z, ddof=1) / np.sqrt(len(z)))
    _, p = sps.ttest_ind(u, w, equal_var=False, nan_policy="omit")
    return z_mean, z_se, float(p), len(u), len(w)


def compute_panel_stats(df: pd.DataFrame, lc: str) -> pd.DataFrame:
    sub = df[df["lc"] == lc]
    u = sub[sub["FireType"] == "WUI"]
    w = sub[sub["FireType"] == "Wildland"]

    specs = []

    def add(key, label, color, uvec, wvec):
        specs.append((key, label, color, np.asarray(uvec, float), np.asarray(wvec, float)))

    add("ignition", "%Human", IGN, u["ignition"], w["ignition"])
    add("elevation", "Elevation", TERRAIN, u["elevation"], w["elevation"])
    add("slope", "Slope", TERRAIN, u["slope"], w["slope"])
    add("ndviM", "NDVI", VEG, u["ndviM"], w["ndviM"])
    add("vs", "Wind", CLIMATE, u["vs"], w["vs"])
    add("ppt", "Precipitation", CLIMATE, u["ppt"], w["ppt"])
    add("tmax", r"$T_{max}$", CLIMATE, u["tmax"], w["tmax"])
    add("tmin", r"$T_{min}$", CLIMATE, u["tmin"], w["tmin"])
    add("tmean", r"$T_{mean}$", CLIMATE, u["tmean"], w["tmean"])
    add("vpdmax", r"$VPD_{max}$", CLIMATE, u["vpdmax"], w["vpdmax"])
    add("FFMC", "FFMC", CLIMATE, u["FFMC"], w["FFMC"])
    add("DMC", "DMC", CLIMATE, u["DMC"], w["DMC"])
    add("DC", "DC", CLIMATE, u["DC"], w["DC"])

    rows = []
    for key, label, color, uvec, wvec in specs:
        zm, se, p, nu, nw = zscore_wui_vs_wild(uvec, wvec)
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
                "color": color,
            }
        )
    return pd.DataFrame(rows)


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


def plot_panel(ax, stats: pd.DataFrame, title: str, show_ylabel: bool):
    stats = stats.sort_values("z", ascending=False).reset_index(drop=True)
    n = len(stats)
    for i in range(0, n, 2):
        ax.axvspan(i + 0.5, i + 1.5, color="0.96", zorder=0)
    ax.axhline(0, color="0.4", lw=0.8, zorder=1)

    for i, r in stats.iterrows():
        x = i + 1
        ci = 1.96 * r["se"] if np.isfinite(r["se"]) else np.nan
        if not np.isfinite(r["z"]):
            continue
        ax.plot([x, x], [r["z"] - ci, r["z"] + ci], color="0.35", lw=1.2, zorder=2)
        ax.plot([x - 0.15, x + 0.15], [r["z"] + ci, r["z"] + ci], color="0.35", lw=1.2)
        ax.plot([x - 0.15, x + 0.15], [r["z"] - ci, r["z"] - ci], color="0.35", lw=1.2)
        ax.scatter(
            x, r["z"], s=55, c=[r["color"]], edgecolors="0.3", linewidths=0.8, zorder=3
        )
        st = stars(r["p"])
        if st:
            y = r["z"] + ci + 0.06 if r["z"] >= 0 else r["z"] - ci - 0.06
            va = "bottom" if r["z"] >= 0 else "top"
            ax.text(x, y, st, ha="center", va=va, fontsize=9)

    ax.set_xticks(range(1, n + 1))
    ax.set_xticklabels(stats["label"].tolist(), rotation=45, ha="right", fontsize=9)
    ax.set_xlim(0.5, n + 0.5)
    ax.set_ylim(-1.1, 1.1)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.tick_params(direction="out")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    if show_ylabel:
        ax.set_yticks(np.arange(-0.9, 1.0, 0.3))
        ax.set_ylabel("Standardized deviations (Z-score)", fontsize=11)
    else:
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)


def main():
    df = load_fires()
    all_stats = pd.concat(
        [compute_panel_stats(df, "Forest"), compute_panel_stats(df, "ShrubGrass")],
        ignore_index=True,
    )
    all_stats.drop(columns=["color"]).to_csv(OUT_CSV, index=False)
    print("Saved", OUT_CSV)

    veg = all_stats[all_stats.variable.isin(["ndviM"])]
    print("\nVegetation z-scores:")
    print(veg[["lc", "label", "z", "p", "n_WUI", "n_Wildland"]].to_string(index=False))

    fig = plt.figure(figsize=(10.8, 4.4), facecolor="white")
    ax1 = fig.add_axes([0.08, 0.32, 0.41, 0.60])
    ax2 = fig.add_axes([0.53, 0.32, 0.41, 0.60])
    plot_panel(ax1, all_stats[all_stats.lc == "Forest"], "Forest", True)
    plot_panel(ax2, all_stats[all_stats.lc == "ShrubGrass"], "Shrub & Grassland", False)

    ax_leg = fig.add_axes([0.15, -0.01, 0.70, 0.09])
    ax_leg.set_xlim(0, 1)
    ax_leg.set_ylim(0, 1)
    ax_leg.axis("off")
    ax_leg.add_patch(plt.Rectangle((0, 0), 1, 1, fill=False, edgecolor="0.75", lw=0.8))
    cats = [
        ("Climate", CLIMATE),
        ("Vegetation", VEG),
        ("Terrain", TERRAIN),
        ("Ignition", IGN),
    ]
    xs = np.linspace(0.12, 0.78, len(cats))
    for x, (name, c) in zip(xs, cats):
        ax_leg.scatter([x], [0.5], s=55, c=[c], edgecolors="0.3", linewidths=0.8)
        ax_leg.text(x + 0.035, 0.5, name, va="center", fontsize=11)

    fig.savefig(OUT_PNG, dpi=300, facecolor="white", bbox_inches="tight")
    fig.savefig(OUT_PDF, facecolor="white", bbox_inches="tight")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Flag SoCal CalFire WUI / Wildland fires whose ignition day coincides with
Santa Ana wind activity using daily R1D-SAWRI (Guzman-Morales et al.).

SAWRI source: weclima UCSD Google Drive
  R1DSAWRI_010148_123118_red.txt  (1948-01-01 … 2018-12-31)

Thresholds (literature):
  SAWRI > 0  — liberal (any SAW activity in domain)
  SAWRI > 1  — Keeley et al. 2024 SAW-fire definition
"""

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
SAWRI_PATH = BASE / "dataSrc/climate/SAWRI/R1DSAWRI_010148_123118_red.txt"
OUT_DIR = BASE / "dataPrc"
OUT_FIRE = OUT_DIR / "CalFire_SoCal_SAWRI_flags.csv"
OUT_SUM = OUT_DIR / "CalFire_SoCal_SAWRI_summary.csv"


def load_sawri(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df.columns = ["dt", "SAWRI"]
    parsed = pd.to_datetime(df["dt"], format="%m/%d/%y")
    yy = parsed.dt.year % 100
    # File spans 01/01/48–12/31/18 → map 48–99 → 1948–1999, 00–18 → 2000–2018
    year = np.where(yy >= 48, 1900 + yy, 2000 + yy)
    dt = pd.to_datetime(
        {"year": year, "month": parsed.dt.month, "day": parsed.dt.day}
    )
    s = pd.Series(df["SAWRI"].astype(float).to_numpy(), index=dt, name="SAWRI")
    s = s[~s.index.duplicated(keep="first")].sort_index()
    assert s.index.min() == pd.Timestamp("1948-01-01"), s.index.min()
    assert s.index.max() == pd.Timestamp("2018-12-31"), s.index.max()
    return s


def load_fires(name: str) -> gpd.GeoDataFrame:
    g = gpd.read_file(BASE / f"dataPrc/firePrmt/fires/{name}.shp")
    g = g[g["SoCal"] == 1].copy()
    g["FireType"] = "WUI" if name == "Urban-edge" else "Wildland"
    g["IDate"] = pd.to_datetime(g["IDate"])
    return g


def flag_fires(g: gpd.GeoDataFrame, sawri: pd.Series) -> pd.DataFrame:
    df = pd.DataFrame(
        {
            "FireType": g["FireType"].values,
            "fid": g["fid"].values if "fid" in g.columns else g.index,
            "IDate": g["IDate"].values,
            "year": g["year"].values,
            "size": g["size"].values,
            "wind_mean": g["wind_mean"].values if "wind_mean" in g.columns else np.nan,
            "wind_dir": g["wind_dir"].values if "wind_dir" in g.columns else np.nan,
        }
    )
    df["IDate"] = pd.to_datetime(df["IDate"]).dt.normalize()
    df["in_SAWRI_period"] = (df["IDate"] >= sawri.index.min()) & (
        df["IDate"] <= sawri.index.max()
    )
    df["SAWRI"] = df["IDate"].map(sawri)
    df["SAW_gt0"] = df["SAWRI"] > 0
    df["SAW_gt1"] = df["SAWRI"] > 1  # Keeley et al. 2024
    return df


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    sub = df[df["in_SAWRI_period"]].copy()
    for ft, g in sub.groupby("FireType"):
        n = len(g)
        for thr_name, col in [("SAWRI>0", "SAW_gt0"), ("SAWRI>1", "SAW_gt1")]:
            n_saw = int(g[col].sum())
            sizes_saw = g.loc[g[col], "size"]
            sizes_nosaw = g.loc[~g[col], "size"]
            # "very large": top decile within SoCal fire type (period-restricted)
            p90 = g["size"].quantile(0.9)
            n_large = int((g["size"] >= p90).sum())
            n_large_saw = int(((g["size"] >= p90) & g[col]).sum())
            rows.append(
                {
                    "FireType": ft,
                    "threshold": thr_name,
                    "n_fires_SoCal_in_SAWRI_era": n,
                    "n_SAW_ignition_day": n_saw,
                    "frac_SAW": n_saw / n if n else np.nan,
                    "mean_size_SAW": float(sizes_saw.mean()) if len(sizes_saw) else np.nan,
                    "mean_size_nonSAW": float(sizes_nosaw.mean()) if len(sizes_nosaw) else np.nan,
                    "median_size_SAW": float(sizes_saw.median()) if len(sizes_saw) else np.nan,
                    "median_size_nonSAW": float(sizes_nosaw.median()) if len(sizes_nosaw) else np.nan,
                    "size_p90_km2": float(p90),
                    "n_large_ge_p90": n_large,
                    "n_large_on_SAW": n_large_saw,
                    "frac_large_that_are_SAW": n_large_saw / n_large if n_large else np.nan,
                    "frac_SAW_that_are_large": n_large_saw / n_saw if n_saw else np.nan,
                    "frac_nonSAW_that_are_large": (
                        (n_large - n_large_saw) / (n - n_saw) if (n - n_saw) else np.nan
                    ),
                    "mean_gridMET_vs_SAW": float(g.loc[g[col], "wind_mean"].mean()),
                    "mean_gridMET_vs_nonSAW": float(g.loc[~g[col], "wind_mean"].mean()),
                }
            )
    # also overall SoCal combined
    for thr_name, col in [("SAWRI>0", "SAW_gt0"), ("SAWRI>1", "SAW_gt1")]:
        n = len(sub)
        n_saw = int(sub[col].sum())
        rows.append(
            {
                "FireType": "All_SoCal",
                "threshold": thr_name,
                "n_fires_SoCal_in_SAWRI_era": n,
                "n_SAW_ignition_day": n_saw,
                "frac_SAW": n_saw / n if n else np.nan,
            }
        )
    return pd.DataFrame(rows)


def main():
    sawri = load_sawri(SAWRI_PATH)
    print(
        f"SAWRI days: {sawri.index.min().date()} → {sawri.index.max().date()}; "
        f"SAWRI>0: {(sawri > 0).sum():,}; SAWRI>1: {(sawri > 1).sum():,}"
    )

    parts = [flag_fires(load_fires(n), sawri) for n in ("Urban-edge", "Wildland")]
    fires = pd.concat(parts, ignore_index=True)
    fires.to_csv(OUT_FIRE, index=False)
    print("Saved", OUT_FIRE)

    # period coverage
    for ft, g in fires.groupby("FireType"):
        n = len(g)
        nin = int(g["in_SAWRI_period"].sum())
        print(f"{ft} SoCal: {n} fires; in SAWRI era (≤2018): {nin} ({100*nin/n:.1f}%)")

    summary = summarize(fires)
    summary.to_csv(OUT_SUM, index=False)
    print("Saved", OUT_SUM)

    print("\n=== SoCal ignition-day Santa Ana fractions (SAWRI era) ===")
    show = summary[summary["FireType"].isin(["WUI", "Wildland"])].copy()
    for _, r in show.iterrows():
        print(
            f"{r['FireType']:9s} {r['threshold']}: "
            f"{r['n_SAW_ignition_day']:.0f}/{r['n_fires_SoCal_in_SAWRI_era']:.0f} "
            f"= {100*r['frac_SAW']:.1f}% | "
            f"mean size SAW {r['mean_size_SAW']:.1f} vs non-SAW {r['mean_size_nonSAW']:.1f} km² | "
            f"large (≥p90) on SAW: {100*r['frac_large_that_are_SAW']:.1f}% | "
            f"gridMET vs SAW {r['mean_gridMET_vs_SAW']:.2f} vs non {r['mean_gridMET_vs_nonSAW']:.2f} m/s"
        )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Forest vs shrub/grassland composition *inside CalFire perimeters*,
stratified by the fire's dominant class (`lc`) and FireType.

Not landscape-wide WUI vs wildland. Each fire is a sample: among
pixels classified as forest or shrub/grassland, report the two fractions.

Default land cover: 2020 MCD12Q1 FAO-LCCS1 (LC_Prop1), same classes as
code/wui_vs_wildland_lc_fractions_2020.py.

  Forest          = 11–16, 21, 22
  Shrub/grassland = 31–32 (herbaceous) + 41–43 (shrublands)
"""

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.features import geometry_mask
from rasterio.windows import Window, from_bounds
from scipy.stats import mannwhitneyu

BASE = Path(__file__).resolve().parents[1]
FIRE = BASE / "dataPrc/firePrmt/CalFire.shp"
LCCS1 = BASE / "dataPrc/landcover/LCCS1_2020_3310.tif"
OUT = BASE / "dataPrc/CalFire_fire_perimeter_forest_SG_fractions.csv"

FOREST_LCCS = np.array([11, 12, 13, 14, 15, 16, 21, 22], dtype=np.uint8)
GRASS_LCCS = np.array([31, 32], dtype=np.uint8)
SHRUB_LCCS = np.array([41, 42, 43], dtype=np.uint8)
SG_LCCS = np.concatenate([GRASS_LCCS, SHRUB_LCCS])


def zonal_frac(g: gpd.GeoDataFrame, path: Path) -> pd.DataFrame:
    rows = []
    with rasterio.open(path) as src:
        H, W = src.height, src.width
        for _, r in g.iterrows():
            geom = r.geometry
            if geom is None or geom.is_empty:
                continue
            minx, miny, maxx, maxy = geom.bounds
            win = from_bounds(minx, miny, maxx, maxy, src.transform)
            win = win.intersection(Window(0, 0, W, H))
            if win.width <= 0 or win.height <= 0:
                continue
            win = win.round_offsets().round_lengths()
            if win.width <= 0 or win.height <= 0:
                continue
            arr = src.read(1, window=win)
            t = src.window_transform(win)
            mask = geometry_mask([geom], out_shape=arr.shape, transform=t, invert=True)
            vals = arr[mask]
            n_f = int(np.isin(vals, FOREST_LCCS).sum())
            n_g = int(np.isin(vals, GRASS_LCCS).sum())
            n_sh = int(np.isin(vals, SHRUB_LCCS).sum())
            n_s = n_g + n_sh
            n_v = n_f + n_s
            rows.append(
                {
                    "fid": int(r.fid),
                    "lc": r.lc,
                    "FireType": r.FireType,
                    "size_km2": float(r.size_km2),
                    "n_forest": n_f,
                    "n_grass": n_g,
                    "n_shrub": n_sh,
                    "n_sg": n_s,
                    "n_veg": n_v,
                    "frac_forest": n_f / n_v if n_v else np.nan,
                    "frac_grass": n_g / n_v if n_v else np.nan,
                    "frac_shrub": n_sh / n_v if n_v else np.nan,
                    "frac_sg": n_s / n_v if n_v else np.nan,
                }
            )
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame, name: str) -> pd.DataFrame:
    recs = []
    print(f"\n======== {name} ========")
    for lc in ["Forest", "ShrubGrass"]:
        print(f"\n-- fires dominated by {lc} --")
        w = df[(df.lc == lc) & (df.FireType == "WUI") & df.frac_sg.notna()]
        v = df[(df.lc == lc) & (df.FireType == "Wildland") & df.frac_sg.notna()]
        p = np.nan
        if len(w) >= 2 and len(v) >= 2:
            p = mannwhitneyu(w.frac_sg, v.frac_sg, alternative="two-sided").pvalue
        for ft, s in [("WUI", w), ("Wildland", v)]:
            n_v = int(s.n_forest.sum() + s.n_sg.sum())
            pooled_f = s.n_forest.sum() / n_v if n_v else np.nan
            pooled_s = s.n_sg.sum() / n_v if n_v else np.nan
            rec = {
                "scheme": name,
                "dominant_lc": lc,
                "FireType": ft,
                "n_fires": len(s),
                "mean_frac_forest": s.frac_forest.mean(),
                "mean_frac_grass": s.frac_grass.mean(),
                "mean_frac_shrub": s.frac_shrub.mean(),
                "mean_frac_sg": s.frac_sg.mean(),
                "median_frac_sg": s.frac_sg.median(),
                "pooled_frac_forest": pooled_f,
                "pooled_frac_grass": s.n_grass.sum() / n_v if n_v else np.nan,
                "pooled_frac_shrub": s.n_shrub.sum() / n_v if n_v else np.nan,
                "pooled_frac_sg": pooled_s,
                "p_mwu_frac_sg": p,
            }
            recs.append(rec)
            print(
                f"  {ft:9s} n={len(s):4d}  "
                f"mean forest={100 * s.frac_forest.mean():5.1f}%  "
                f"grass={100 * s.frac_grass.mean():5.1f}%  "
                f"shrub={100 * s.frac_shrub.mean():5.1f}%  "
                f"SG={100 * s.frac_sg.mean():5.1f}%  |  "
                f"area-pooled forest={100 * pooled_f:5.1f}%  SG={100 * pooled_s:5.1f}%"
            )
        dsg = w.frac_sg.mean() - v.frac_sg.mean()
        print(
            f"  Δ mean SG (WUI−wildland) = {100 * dsg:+.1f} pp  "
            f"Mann–Whitney p={p:.3g}"
        )
    return pd.DataFrame(recs)


def main():
    g = gpd.read_file(FIRE)
    g = g[g["lc"].isin(["Forest", "ShrubGrass"]) & g["FireType"].isin(["WUI", "Wildland"])].copy()
    g["fid"] = g["fid"].astype(int)
    print(f"Fires: {len(g)}  (Forest {int((g.lc=='Forest').sum())}, "
          f"ShrubGrass {int((g.lc=='ShrubGrass').sum())})")

    print("Zonal LCCS1 2020 inside perimeters …")
    lccs = zonal_frac(g, LCCS1)
    s1 = summarize(lccs, "LCCS1_2020")
    s1.to_csv(OUT, index=False)
    lccs.to_csv(BASE / "dataPrc/CalFire_fire_perimeter_LCCS1_by_fire.csv", index=False)
    print("\nSaved", OUT)


if __name__ == "__main__":
    main()

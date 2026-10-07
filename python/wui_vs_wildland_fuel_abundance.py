#!/usr/bin/env python3
"""Compare fuel/biomass abundance: WUI vs wildland fires (Line 229 framework).

Datasets (pre-fire year = ignition year - 1):
  1) ESA CCI forest AGB (Mg/ha) — exact map year == pre-fire year only
     (2007, 2010, 2015–2022); fires without a matching map year are excluded
  2) RAP herbaceous NPP / cover — ShrubGrass fires; years 1986–2021
  3) LANDFIRE LF2023 FVC + grass/shrub FBFM40 fraction — all fires (static map)

Mirrors NDVI analysis: zonal mean inside CalFire perimeters, stratified by
lc × FireType, then Wildland-referenced z-scores + Welch t-tests.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import ee
import geopandas as gpd
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
FIRE_SHP = ROOT / "dataPrc/firePrmt/CalFire.shp"
OUT_CSV = ROOT / "dataPrc/WUI_vs_wildland_fuel_abundance_by_fire.csv"
OUT_SUM = ROOT / "dataPrc/WUI_vs_wildland_fuel_abundance_summary.csv"
CACHE = ROOT / "dataPrc/_tmp_fuel_abundance_batches"
CCI_YEARS = [2007, 2010, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022]


def init_ee():
    ee.Initialize(project="globfire")


def exact_cci_year(pre_year: int) -> int | None:
    """Require an exact CCI map year equal to the pre-fire year; else exclude."""
    y = int(pre_year)
    return y if y in CCI_YEARS else None


def gdf_to_fc(gdf: gpd.GeoDataFrame) -> ee.FeatureCollection:
    gdf = gdf.to_crs(4326)
    feats = []
    for _, r in gdf.iterrows():
        geom = r.geometry
        if geom is None or geom.is_empty:
            continue
        # EE expects GeoJSON geometry dict
        gj = json.loads(gpd.GeoSeries([geom], crs=4326).to_json())["features"][0]["geometry"]
        feats.append(
            ee.Feature(
                ee.Geometry(gj),
                {
                    "fid": int(r["fid"]),
                    "year": int(r["year"]),
                    "pre_year": int(r["pre_year"]),
                    "FireType": str(r["FireType"]),
                    "lc": str(r["lc"]),
                    "size_km2": float(r["size_km2"]) if pd.notna(r["size_km2"]) else None,
                    "ndvi": float(r["ndvi"]) if pd.notna(r.get("ndvi")) else None,
                },
            )
        )
    return ee.FeatureCollection(feats)


def rap_images(pre_year: int):
    """Herbaceous abundance proxies for ShrubGrass fires."""
    y = int(np.clip(pre_year, 1986, 2021))
    npp = (
        ee.ImageCollection("projects/rap-data-365417/assets/npp-partitioned-v3")
        .filter(ee.Filter.eq("system:index", str(y)))
        .first()
    )
    cover = (
        ee.ImageCollection("projects/rap-data-365417/assets/vegetation-cover-v3")
        .filter(ee.Filter.eq("system:index", str(y)))
        .first()
    )
    mat = (
        ee.ImageCollection("projects/rap-data-365417/assets/gridmet-MAT")
        .filter(ee.Filter.eq("system:index", str(y)))
        .first()
        .select("MAT")
    )
    # NPP appears as ~gC/m^2 in this asset (CA medians hundreds).
    # Aboveground fraction from Hui & Jackson-style logistic used by RAP production apps.
    herb_npp = npp.select("afgNPP").add(npp.select("pfgNPP")).rename("rap_herb_npp")
    shr_npp = npp.select("shrNPP").rename("rap_shr_npp")
    ag_frac = mat.expression("1 / (1 + exp(-(0.0402 * MAT - 1.3964)))", {"MAT": mat})
    # gC/m2 -> lbs/acre dry matter: /0.45 C fraction * 8.92179
    herb_bio = herb_npp.multiply(ag_frac).divide(0.45).multiply(8.92179).rename("rap_herb_lbs_ac")
    herb_shr_cover = (
        cover.select("AFG")
        .add(cover.select("PFG"))
        .add(cover.select("SHR"))
        .rename("rap_herb_shr_cover")
    )
    return herb_npp.addBands([shr_npp, herb_bio, herb_shr_cover])


def cci_image(map_year: int):
    img = (
        ee.ImageCollection("projects/sat-io/open-datasets/ESA/ESA_CCI_AGB")
        .filter(ee.Filter.stringContains("system:index", str(int(map_year))))
        .first()
        .select("AGB")
        .rename("cci_agb")
    )
    return img, int(map_year)


def landfire_images():
    fvc = (
        ee.Image("projects/sat-io/open-datasets/landfire/FUEL/FVC/FVC_LC_240")
        .select("FVC")
        .rename("lf_fvc")
    )
    f40 = ee.Image("projects/sat-io/open-datasets/landfire/FUEL/FBFM40/F40_LC_240").select("F40")
    # Scott & Burgan: GR 101–109, GS 121–124, SH 141–149
    gs = (
        f40.gte(101)
        .And(f40.lte(109))
        .Or(f40.gte(121).And(f40.lte(124)))
        .Or(f40.gte(141).And(f40.lte(149)))
        .rename("lf_gs_fuel")
    )
    return fvc, gs


def reduce_batch(fc: ee.FeatureCollection, image: ee.Image, scale: int) -> list[dict]:
    reduced = image.reduceRegions(
        collection=fc,
        reducer=ee.Reducer.mean(),
        scale=scale,
        tileScale=4,
    )
    return reduced.getInfo()["features"]


def run_zonal():
    init_ee()
    CACHE.mkdir(parents=True, exist_ok=True)

    g = gpd.read_file(FIRE_SHP)
    g = g[g["lc"].isin(["Forest", "ShrubGrass"])].copy()
    g["fid"] = g["fid"].astype(int)
    g["year"] = g["year"].astype(int)
    g["pre_year"] = (g["year"] - 1).clip(lower=1986)
    g["geometry"] = g.geometry.simplify(30, preserve_topology=True)  # meters in EPSG:6414

    fvc, gs = landfire_images()
    rows = []

    # ---- LANDFIRE (static): one pass in year chunks for request size ----
    print("LANDFIRE FVC / grass-shrub fuel fraction …")
    for pre_year, sub in g.groupby("pre_year"):
        cache = CACHE / f"lf_{pre_year}.json"
        if cache.exists():
            batch = json.loads(cache.read_text())
        else:
            fc = gdf_to_fc(sub)
            img = fvc.addBands(gs)
            batch = reduce_batch(fc, img, scale=30)
            cache.write_text(json.dumps(batch))
            print(f"  LF pre_year={pre_year} n={len(sub)}")
            time.sleep(0.2)
        for f in batch:
            p = f["properties"]
            rows.append(
                {
                    "fid": p["fid"],
                    "lf_fvc": p.get("lf_fvc"),
                    "lf_gs_fuel": p.get("lf_gs_fuel"),
                }
            )
    lf = pd.DataFrame(rows).drop_duplicates("fid")

    # ---- RAP (ShrubGrass only; pre_year in 1986–2021) ----
    print("RAP herbaceous abundance …")
    rap_rows = []
    sg = g[(g["lc"] == "ShrubGrass") & (g["pre_year"] >= 1986) & (g["pre_year"] <= 2021)]
    for pre_year, sub in sg.groupby("pre_year"):
        cache = CACHE / f"rap_{pre_year}.json"
        if cache.exists():
            batch = json.loads(cache.read_text())
        else:
            fc = gdf_to_fc(sub)
            img = rap_images(int(pre_year))
            batch = reduce_batch(fc, img, scale=30)
            cache.write_text(json.dumps(batch))
            print(f"  RAP pre_year={pre_year} n={len(sub)}")
            time.sleep(0.2)
        for f in batch:
            p = f["properties"]
            rap_rows.append(
                {
                    "fid": p["fid"],
                    "rap_herb_npp": p.get("rap_herb_npp"),
                    "rap_shr_npp": p.get("rap_shr_npp"),
                    "rap_herb_lbs_ac": p.get("rap_herb_lbs_ac"),
                    "rap_herb_shr_cover": p.get("rap_herb_shr_cover"),
                    "rap_year": int(pre_year),
                }
            )
    rap = pd.DataFrame(rap_rows).drop_duplicates("fid")

    # ---- CCI AGB (exact pre-fire map year only; Forest + ShrubGrass) ----
    print("ESA CCI AGB (exact pre-fire year) …")
    cci_rows = []
    fo = g[g["lc"].isin(["Forest", "ShrubGrass"])].copy()
    fo["cci_year"] = fo["pre_year"].map(exact_cci_year)
    fo = fo[fo["cci_year"].notna()]
    for cci_year, sub in fo.groupby("cci_year"):
        cache = CACHE / f"cci_exact_{int(cci_year)}.json"
        if cache.exists():
            batch = json.loads(cache.read_text())
        else:
            fc = gdf_to_fc(sub)
            img, _ = cci_image(int(cci_year))
            img = img.updateMask(img.gt(0)).rename("cci_agb")
            batch = reduce_batch(fc, img, scale=100)
            cache.write_text(json.dumps(batch))
            print(f"  CCI year={cci_year} n={len(sub)}")
            time.sleep(0.2)
        for f in batch:
            p = f["properties"]
            cci_rows.append(
                {
                    "fid": p["fid"],
                    "cci_agb": p.get("cci_agb", p.get("mean")),
                    "cci_year": int(cci_year),
                }
            )
    cci = pd.DataFrame(cci_rows).drop_duplicates("fid")

    out = g.drop(columns="geometry").merge(lf, on="fid", how="left")
    out = out.merge(rap, on="fid", how="left")
    out = out.merge(cci, on="fid", how="left")
    out.to_csv(OUT_CSV, index=False)
    print("Wrote", OUT_CSV, "n=", len(out))
    return out


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    """Wildland-referenced z-score + Welch t-test, matching Fig3 NDVI logic."""
    # Note: raw LANDFIRE FVC class codes are thematic — do not treat mean FVC
    # as continuous abundance. Use lf_gs_fuel (grass/shrub FBFM40 fraction).
    metrics = [
        ("ShrubGrass", "rap_herb_lbs_ac", "RAP herbaceous biomass (lbs/ac)"),
        ("ShrubGrass", "rap_herb_npp", "RAP herbaceous NPP"),
        ("ShrubGrass", "rap_herb_shr_cover", "RAP AFG+PFG+SHR cover (%)"),
        ("ShrubGrass", "rap_shr_npp", "RAP shrub NPP"),
        ("ShrubGrass", "lf_gs_fuel", "LANDFIRE grass/shrub fuel frac (ShrubGrass)"),
        ("Forest", "cci_agb", "ESA CCI AGB (Mg/ha)"),
        ("Forest", "lf_gs_fuel", "LANDFIRE grass/shrub fuel frac (Forest)"),
        ("ShrubGrass", "ndvi", "MODIS NDVI (reference)"),
        ("Forest", "ndvi", "MODIS NDVI (reference)"),
    ]
    rows = []
    for lc, col, label in metrics:
        sub = df[df["lc"] == lc][["FireType", col]].dropna()
        wui = sub.loc[sub["FireType"] == "WUI", col].astype(float)
        wild = sub.loc[sub["FireType"] == "Wildland", col].astype(float)
        if len(wui) < 5 or len(wild) < 5:
            continue
        mu_w, sd_w = wild.mean(), wild.std(ddof=1)
        z = (wui.mean() - mu_w) / sd_w if sd_w > 0 else np.nan
        t, p = stats.ttest_ind(wui, wild, equal_var=False, nan_policy="omit")
        rows.append(
            {
                "lc": lc,
                "metric": col,
                "label": label,
                "n_WUI": int(wui.shape[0]),
                "n_Wildland": int(wild.shape[0]),
                "mean_WUI": float(wui.mean()),
                "mean_Wildland": float(wild.mean()),
                "diff_WUI_minus_Wildland": float(wui.mean() - wild.mean()),
                "z_WUI_vs_Wildland": float(z),
                "welch_t": float(t),
                "p_value": float(p),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_SUM, index=False)
    print("Wrote", OUT_SUM)
    return summary


def main():
    df = run_zonal()
    summary = summarize(df)
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.width", 160)
    print("\n=== Summary (positive z => higher in WUI, like NDVI Line 229) ===")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

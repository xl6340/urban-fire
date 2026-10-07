#!/usr/bin/env python3
"""
Half-decadal landscape VPDmax: all-year vs fire-season (May–Oct).

Downloads CA-wide PRISM daily-VPDmax temporal means from GEE (4 km),
then computes spatial means with local WUI / urban / desert masks.
"""

from __future__ import annotations

import io
import time
import urllib.request
from pathlib import Path

import ee
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio import features
from rasterio.transform import from_bounds
from rasterio.warp import Resampling, reproject

BASE = Path(__file__).resolve().parents[1]
CACHE = BASE / "dataPrc/_tmp_vpd_prism_means"
OUT = BASE / "dataPrc/vpdmax_landscape_allyear_vs_fireseason.csv"
OUT_CMP = BASE / "dataFig/vpd/betaVPD_fireseason_compare.csv"
LC_TYPE2 = BASE / "dataSrc/landcover/data/MCD12Q1.061_LC_Type2_doy2020001000000_aid0001.tif"

PERIODS = [
    (1990, 1994),
    (1995, 1999),
    (2000, 2004),
    (2005, 2009),
    (2010, 2014),
    (2015, 2019),
    (2020, 2024),
]
WUI_FOR_PERIOD = {1990: 1990, 1995: 1990, 2000: 2000, 2005: 2000, 2010: 2010, 2015: 2010, 2020: 2020}
SCALE = 4000  # m, PRISM-native-ish


def init_ee():
    ee.Initialize(project="globfire")


def ca_geometry():
    return (
        ee.FeatureCollection("TIGER/2018/States")
        .filter(ee.Filter.eq("NAME", "California"))
        .geometry()
    )


def prism_mean_image(start: int, end: int, fire_season_only: bool) -> ee.Image:
    col = (
        ee.ImageCollection("OREGONSTATE/PRISM/AN81d")
        .filterDate(f"{start}-01-01", f"{end + 1}-01-01")
        .select("vpdmax")
    )
    if fire_season_only:
        col = col.filter(ee.Filter.calendarRange(5, 10, "month"))
    return col.mean().rename("vpdmax").toFloat()


def download_geotiff(img: ee.Image, region: ee.Geometry, out_path: Path):
    if out_path.exists() and out_path.stat().st_size > 1000:
        return
    out_path.parent.mkdir(parents=True, exist_ok=True)
    url = img.clip(region).getDownloadURL(
        {
            "name": out_path.stem,
            "region": region,
            "scale": SCALE,
            "crs": "EPSG:4326",
            "format": "GEO_TIFF",
        }
    )
    print("  downloading", out_path.name)
    with urllib.request.urlopen(url, timeout=600) as r:
        data = r.read()
    out_path.write_bytes(data)


def rasterize_on(ds_ref: rasterio.io.DatasetReader, geoms) -> np.ndarray:
    if len(geoms) == 0:
        return np.zeros((ds_ref.height, ds_ref.width), dtype=bool)
    return features.geometry_mask(
        geoms,
        out_shape=(ds_ref.height, ds_ref.width),
        transform=ds_ref.transform,
        invert=True,
    )


def urban_on_grid(ds_ref: rasterio.io.DatasetReader) -> np.ndarray:
    """LC_Type2 == 13 warped to VPD grid."""
    dst = np.zeros((ds_ref.height, ds_ref.width), dtype=np.uint8)
    with rasterio.open(LC_TYPE2) as src:
        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=ds_ref.transform,
            dst_crs=ds_ref.crs,
            resampling=Resampling.nearest,
        )
    return dst == 13


def spatial_mean(arr: np.ndarray, mask: np.ndarray) -> float:
    m = mask & np.isfinite(arr)
    if not np.any(m):
        return float("nan")
    return float(np.mean(arr[m]))


def main():
    init_ee()
    CACHE.mkdir(parents=True, exist_ok=True)
    ca = ca_geometry()

    print("Loading vector masks …")
    wui = {
        y: gpd.read_file(BASE / f"dataPrc/WUI/WUI{y}.shp").to_crs(4326)
        for y in (1990, 2000, 2010, 2020)
    }
    epa = gpd.read_file(BASE / "dataPrc/boundary/EPALevel3.shp").to_crs(4326)
    deserts = epa[epa["ecoregion"].astype(str).str.lower().eq("deserts")]
    ca_gdf = gpd.read_file(
        # approximate from EE not needed; use unary of WUI+wild from state shapefile if available
        BASE / "dataPrc/boundary/EPALevel3.shp"
    ).to_crs(4326)
    # Better: dissolve CA from TIGER via local if exists
    states = None
    for p in [
        BASE / "dataPrc/boundary/CA_State.shp",
        BASE / "dataPrc/boundary/California.shp",
    ]:
        if p.exists():
            states = gpd.read_file(p).to_crs(4326)
            break
    if states is None:
        ca_poly = None
    else:
        ca_poly = states.union_all() if hasattr(states, "union_all") else states.unary_union

    rows = []
    urban_cache = {}

    for start, end in PERIODS:
        for season, fire_only in [("all_year", False), ("fire_season_MayOct", True)]:
            tag = f"vpdmax_{start}_{end}_{'fs' if fire_only else 'ay'}.tif"
            tif = CACHE / tag
            img = prism_mean_image(start, end, fire_only)
            download_geotiff(img, ca, tif)

            with rasterio.open(tif) as ds:
                vpd = ds.read(1).astype(float)
                nodata = ds.nodata
                if nodata is not None:
                    vpd[vpd == nodata] = np.nan
                vpd[~np.isfinite(vpd)] = np.nan
                # PRISM sometimes uses large negatives
                vpd[vpd < -100] = np.nan

                if id(ds.transform) not in urban_cache:
                    urban_cache["grid"] = urban_on_grid(ds)
                urban = urban_cache["grid"]
                # rebuild urban if grid shape differs
                if urban.shape != vpd.shape:
                    urban = urban_on_grid(ds)
                    urban_cache["grid"] = urban

                desert_m = rasterize_on(ds, list(deserts.geometry))
                if ca_poly is not None:
                    ca_m = rasterize_on(ds, [ca_poly])
                else:
                    ca_m = np.isfinite(vpd)

                wui_y = WUI_FOR_PERIOD[start]
                wui_m = rasterize_on(ds, list(wui[wui_y].geometry))

                wui_domain = ca_m & wui_m & ~urban & np.isfinite(vpd)
                wild_domain = ca_m & ~wui_m & ~urban & ~desert_m & np.isfinite(vpd)
                ca_exdesert = ca_m & ~desert_m & ~urban & np.isfinite(vpd)

                wui_hpa = spatial_mean(vpd, wui_domain)
                wild_hpa = spatial_mean(vpd, wild_domain)
                ca_hpa = spatial_mean(vpd, ca_exdesert)

            rows.append(
                {
                    "StartYear": start,
                    "EndYear": end,
                    "Period": f"{start}s",
                    "season": season,
                    "WUI_map_year": wui_y,
                    "vpd_wui_hPa": wui_hpa,
                    "vpd_wild_hPa": wild_hpa,
                    "vpd_ca_exdesert_hPa": ca_hpa,
                    "vpd_wui_kPa": wui_hpa / 10.0,
                    "vpd_wild_kPa": wild_hpa / 10.0,
                    "vpd_ca_exdesert_kPa": ca_hpa / 10.0,
                    "n_wui_px": int(wui_domain.sum()),
                    "n_wild_px": int(wild_domain.sum()),
                }
            )
            print(
                f"{start}-{end} {season:22s}  "
                f"WUI={wui_hpa/10:.3f}  Wild={wild_hpa/10:.3f}  CA={ca_hpa/10:.3f} kPa  "
                f"(px {int(wui_domain.sum())}/{int(wild_domain.sum())})"
            )

    df = pd.DataFrame(rows)
    df.to_csv(OUT, index=False)
    print("Wrote", OUT)

    beta = pd.read_excel(BASE / "dataFig/vpd/betaVPD.xlsx")
    beta = beta.dropna(subset=["vpd_wui_kPa"]).copy()
    beta["StartYear"] = beta["Period"].astype(str).str.replace("s", "", regex=False).astype(int)
    cmp = df.merge(
        beta[["StartYear", "vpd_wui_kPa", "vpd_wild_kPa", "vpd_region_kPa"]].rename(
            columns={
                "vpd_wui_kPa": "published_vpd_wui_kPa",
                "vpd_wild_kPa": "published_vpd_wild_kPa",
                "vpd_region_kPa": "published_vpd_region_kPa",
            }
        ),
        on="StartYear",
        how="left",
    )
    cmp.to_csv(OUT_CMP, index=False)
    print("Wrote", OUT_CMP)

    print("\nFire-season minus all-year (kPa):")
    ay = df[df.season == "all_year"].set_index("StartYear")
    fs = df[df.season == "fire_season_MayOct"].set_index("StartYear")
    delta = pd.DataFrame(
        {
            "dWUI": fs["vpd_wui_kPa"] - ay["vpd_wui_kPa"],
            "dWild": fs["vpd_wild_kPa"] - ay["vpd_wild_kPa"],
            "WUI_fs": fs["vpd_wui_kPa"],
            "Wild_fs": fs["vpd_wild_kPa"],
            "WUI_ay": ay["vpd_wui_kPa"],
            "Wild_ay": ay["vpd_wild_kPa"],
        }
    )
    print(delta.round(3).to_string())


if __name__ == "__main__":
    main()

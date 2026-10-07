#!/usr/bin/env python3
"""
Forest vs shrub/grassland proportions in California (2020) using
MCD12Q1 FAO-LCCS1 (LC_Prop1):

  Forest          = 11–16, 21, 22
  Shrub/grassland = 31–32 (herbaceous) + 41–43 (shrublands)

  WUI pixels   = inside WUI2020, exclude urban (MODIS LC_Type2 == 13 mask)
  Wildland     = outside WUI2020 & outside urban, exclude EPA Level-3 "Deserts"

Denominator = Forest + Shrub/grassland only
(barren/snow/water/unclassified excluded from both sides).
"""

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio import features
from rasterio.warp import Resampling, reproject

BASE = Path(__file__).resolve().parents[1]
LCCS1_PATH = BASE / "dataPrc/landcover/LCCS1_2020_3310.tif"
# Urban not in LCCS1 legend — mask with LC_Type2 class 13
LC_TYPE2_PATH = BASE / "dataSrc/landcover/data/MCD12Q1.061_LC_Type2_doy2020001000000_aid0001.tif"
OUT_CSV = BASE / "dataPrc/WUI_vs_wildland_forest_shrubgrass_fractions_2020_LCCS1.csv"

# Table 8 FAO-LCCS1
FOREST = {11, 12, 13, 14, 15, 16, 21, 22}
SHRUB_GRASS = {31, 32, 41, 42, 43}
GRASS = {31, 32}
SHRUB = {41, 42, 43}


def rasterize_mask(geoms, shape, transform):
    if len(geoms) == 0:
        return np.zeros(shape, dtype=bool)
    return features.geometry_mask(
        geoms, out_shape=shape, transform=transform, invert=True
    )


def urban_mask_on_grid(lc_type2_path, dst_crs, dst_transform, shape):
    dst = np.zeros(shape, dtype=np.uint8)
    with rasterio.open(lc_type2_path) as src:
        reproject(
            source=rasterio.band(src, 1),
            destination=dst,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=Resampling.nearest,
        )
    return dst == 13


def summarize(name, mask, lc):
    vals = lc[mask]
    n_forest = int(np.isin(vals, list(FOREST)).sum())
    n_grass = int(np.isin(vals, list(GRASS)).sum())
    n_shrub = int(np.isin(vals, list(SHRUB)).sum())
    n_sg = n_grass + n_shrub
    n_veg = n_forest + n_sg
    px_km2 = 0.25  # 500 m
    return {
        "domain": name,
        "scheme": "LCCS1_Forest_11-16_21-22__SG_31-32_41-43",
        "n_pixels_domain": int(vals.size),
        "n_forest": n_forest,
        "n_grass": n_grass,
        "n_shrub": n_shrub,
        "n_shrub_grass": n_sg,
        "n_forest_plus_shrubgrass": n_veg,
        "frac_forest": n_forest / n_veg if n_veg else np.nan,
        "frac_shrub_grass": n_sg / n_veg if n_veg else np.nan,
        "frac_grass_of_veg": n_grass / n_veg if n_veg else np.nan,
        "frac_shrub_of_veg": n_shrub / n_veg if n_veg else np.nan,
        "forest_km2": n_forest * px_km2,
        "shrub_grass_km2": n_sg * px_km2,
        "year_lc": 2020,
        "year_wui": 2020,
        "lc_source": "MCD12Q1.061 LC_Prop1 (FAO-LCCS1)",
    }


def main():
    with rasterio.open(LCCS1_PATH) as src:
        lc = src.read(1)
        transform = src.transform
        dst_crs = src.crs
        shape = lc.shape

    # Vectors are EPSG:6414 (NAD83(2011) Teale Albers); LCCS1 geotiff is tagged
    # EPSG:3310. Coordinates already align — do not reproject (6414→3310 breaks geoms).
    def as_raster_crs(gdf):
        return gdf.set_crs(dst_crs, allow_override=True)

    ca = as_raster_crs(gpd.read_file(BASE / "dataPrc/boundary/CA_State.shp"))
    wui = as_raster_crs(gpd.read_file(BASE / "dataPrc/WUI/WUI2020.shp"))
    epa = as_raster_crs(gpd.read_file(BASE / "dataPrc/boundary/EPALevel3.shp"))
    deserts = epa[epa["ecoregion"].astype(str).str.lower().eq("deserts")]

    ca_mask = rasterize_mask(list(ca.geometry), shape, transform)
    wui_mask = rasterize_mask(list(wui.geometry), shape, transform)
    desert_mask = rasterize_mask(list(deserts.geometry), shape, transform)
    print("Building urban mask from LC_Type2 class 13 …")
    urban_mask = urban_mask_on_grid(LC_TYPE2_PATH, dst_crs, transform, shape) & ca_mask

    # valid land inside CA (nodata 0 outside CA in this raster)
    in_ca = ca_mask & (lc != 0) & (lc != 255)

    wui_domain = in_ca & wui_mask & ~urban_mask
    wild_domain = in_ca & ~wui_mask & ~urban_mask & ~desert_mask

    rows = []
    for domain, mask in [
        ("WUI_excl_urban", wui_domain),
        ("Wildland_excl_deserts", wild_domain),
    ]:
        s = summarize(domain, mask, lc)
        rows.append(s)
        print(
            f"{domain}: Forest {100*s['frac_forest']:.1f}% | "
            f"Shrub/grass {100*s['frac_shrub_grass']:.1f}% "
            f"(grass {100*s['frac_grass_of_veg']:.1f}%, shrub {100*s['frac_shrub_of_veg']:.1f}%; "
            f"n_veg={s['n_forest_plus_shrubgrass']:,})"
        )

    # raw class counts among vegetated pixels in each domain
    print("\nLCCS1 class counts (WUI excl. urban, veg only):")
    v = lc[wui_domain]
    veg = np.isin(v, list(FOREST | SHRUB_GRASS))
    u, c = np.unique(v[veg], return_counts=True)
    print(dict(zip(u.tolist(), c.tolist())))
    print("LCCS1 class counts (Wildland excl. deserts, veg only):")
    v = lc[wild_domain]
    veg = np.isin(v, list(FOREST | SHRUB_GRASS))
    u, c = np.unique(v[veg], return_counts=True)
    print(dict(zip(u.tolist(), c.tolist())))

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print("\nSaved", OUT_CSV)


if __name__ == "__main__":
    main()

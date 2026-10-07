#!/usr/bin/env python3
"""
Compute Canadian FWI moisture codes (FFMC, DMC, DC) for CalFire fires.

For each fire with a valid IDate:
  1) Pull daily perimeter-mean gridMET (tmmx, rmin, vs, pr) for a spin-up window
  2) Run Canadian FWI (pyFWI) with default startup codes
  3) Mean FFMC/DMC/DC over IDate-3 .. IDate (4-day mean; same window as VPDmax)

Skips fires with missing IDate. Resumable via intermediate CSV.
"""

from __future__ import annotations

import json
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import ee
import geopandas as gpd
import numpy as np
import pandas as pd
from pyFWI import DC, DMC, FFMC

BASE = Path(__file__).resolve().parents[1]
CALFIRE = BASE / "dataPrc/firePrmt/CalFire.shp"
OUT_CSV = BASE / "dataPrc/CalFire_FWI_4day.csv"
OUT_GPKG = BASE / "dataPrc/CalFire_FWI_4day.gpkg"
CACHE_DIR = BASE / "dataPrc/fwi_cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

SPINUP_DAYS = 365  # days before IDate-3
SCALE = 4000  # gridMET ~4 km
MAX_WORKERS = 4
INIT_FFMC, INIT_DMC, INIT_DC = 85.0, 6.0, 15.0


def init_ee() -> None:
    ee.Initialize()


def load_fires() -> gpd.GeoDataFrame:
    gdf = gpd.read_file(CALFIRE)
    gdf = gdf[gdf["IDate"].notna()].copy()
    gdf["IDate"] = pd.to_datetime(gdf["IDate"])
    # stable id
    if "fid" in gdf.columns:
        gdf["fire_id"] = gdf["fid"].astype(int)
    else:
        gdf["fire_id"] = np.arange(len(gdf))
    gdf = gdf.to_crs(4326)
    return gdf.reset_index(drop=True)


def geom_to_ee(geom) -> ee.Geometry:
    # GeoJSON geometry from shapely
    return ee.Geometry(geom.__geo_interface__)


def fetch_gridmet_series(fire_id: int, geom, idate: pd.Timestamp) -> pd.DataFrame:
    """Perimeter-mean daily gridMET from (IDate-3-SPINUP) to IDate inclusive."""
    cache = CACHE_DIR / f"wx_{fire_id}.csv"
    if cache.exists():
        return pd.read_csv(cache, parse_dates=["date"])

    end = idate.normalize()
    start = end - pd.Timedelta(days=3 + SPINUP_DAYS)
    # EE end date is exclusive
    ee_start = start.strftime("%Y-%m-%d")
    ee_end = (end + pd.Timedelta(days=1)).strftime("%Y-%m-%d")

    g = geom_to_ee(geom)
    col = (
        ee.ImageCollection("IDAHO_EPSCOR/GRIDMET")
        .filterDate(ee_start, ee_end)
        .select(["tmmx", "rmin", "vs", "pr"])
    )

    def _extract(img):
        stats = img.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=g,
            scale=SCALE,
            maxPixels=1e9,
            bestEffort=True,
        )
        return ee.Feature(
            None,
            {
                "date": img.date().format("YYYY-MM-dd"),
                "tmmx": stats.get("tmmx"),
                "rmin": stats.get("rmin"),
                "vs": stats.get("vs"),
                "pr": stats.get("pr"),
            },
        )

    feat = ee.FeatureCollection(col.map(_extract))
    # retry getInfo
    last_err = None
    for attempt in range(4):
        try:
            data = feat.getInfo()["features"]
            break
        except Exception as e:
            last_err = e
            time.sleep(2 ** attempt)
    else:
        raise RuntimeError(f"EE extract failed for fire {fire_id}: {last_err}")

    rows = []
    for f in data:
        p = f["properties"]
        rows.append(
            {
                "date": p.get("date"),
                "tmmx": p.get("tmmx"),
                "rmin": p.get("rmin"),
                "vs": p.get("vs"),
                "pr": p.get("pr"),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError(f"Empty gridMET series for fire {fire_id}")
    df["date"] = pd.to_datetime(df["date"])
    df = df.sort_values("date").drop_duplicates("date").reset_index(drop=True)
    df.to_csv(cache, index=False)
    return df


def run_fwi(df: pd.DataFrame, lat: float) -> pd.DataFrame:
    """Run Canadian FWI on daily perimeter-mean meteorology."""
    out = df.copy()
    # gridMET units: tmmx Kelvin, rmin %, vs m/s, pr mm
    temp_c = out["tmmx"] - 273.15
    rh = out["rmin"].clip(0, 100)
    wind_kph = out["vs"] * 3.6
    rain = out["pr"].fillna(0).clip(lower=0)

    ffmc = np.full(len(out), np.nan)
    dmc = np.full(len(out), np.nan)
    dc = np.full(len(out), np.nan)

    ff_prev, dm_prev, dc_prev = INIT_FFMC, INIT_DMC, INIT_DC
    for i in range(len(out)):
        if not np.isfinite([temp_c.iloc[i], rh.iloc[i], wind_kph.iloc[i], rain.iloc[i]]).all():
            # carry previous codes if a day is missing
            ffmc[i], dmc[i], dc[i] = ff_prev, dm_prev, dc_prev
            continue
        month = int(out["date"].iloc[i].month)
        t = float(temp_c.iloc[i])
        r = float(rh.iloc[i])
        w = float(wind_kph.iloc[i])
        p = float(rain.iloc[i])
        ff_prev = float(FFMC(t, r, w, p, ff_prev))
        dm_prev = float(DMC(t, r, p, dm_prev, lat, month))
        dc_prev = float(DC(t, p, dc_prev, lat, month))
        ffmc[i], dmc[i], dc[i] = ff_prev, dm_prev, dc_prev

    out["FFMC"] = ffmc
    out["DMC"] = dmc
    out["DC"] = dc
    return out


def four_day_means(fwi_df: pd.DataFrame, idate: pd.Timestamp) -> dict:
    end = idate.normalize()
    start = end - pd.Timedelta(days=3)
    sub = fwi_df[(fwi_df["date"] >= start) & (fwi_df["date"] <= end)]
    if len(sub) == 0:
        return {"FFMC_4d": np.nan, "DMC_4d": np.nan, "DC_4d": np.nan, "n_days": 0}
    return {
        "FFMC_4d": float(sub["FFMC"].mean()),
        "DMC_4d": float(sub["DMC"].mean()),
        "DC_4d": float(sub["DC"].mean()),
        "n_days": int(len(sub)),
    }


def process_one(row: pd.Series) -> dict:
    fire_id = int(row["fire_id"])
    idate = pd.Timestamp(row["IDate"]).normalize()
    geom = row.geometry
    # latitude from representative point
    pt = geom.representative_point()
    lat = float(pt.y)

    result_cache = CACHE_DIR / f"fwi_{fire_id}.json"
    if result_cache.exists():
        return json.loads(result_cache.read_text())

    wx = fetch_gridmet_series(fire_id, geom, idate)
    fwi_df = run_fwi(wx, lat=lat)
    # save full series (optional, overwrite wx cache with codes)
    fwi_path = CACHE_DIR / f"fwi_series_{fire_id}.csv"
    fwi_df.to_csv(fwi_path, index=False)

    means = four_day_means(fwi_df, idate)
    out = {
        "fire_id": fire_id,
        "IDate": idate.strftime("%Y-%m-%d"),
        "lat": lat,
        **means,
    }
    result_cache.write_text(json.dumps(out))
    return out


def main(limit: int | None = None, workers: int = MAX_WORKERS) -> None:
    init_ee()
    fires = load_fires()
    print(f"CalFire with IDate: {len(fires)} (skipped missing IDate)")
    if limit is not None:
        fires = fires.iloc[:limit].copy()
        print(f"Limiting to first {limit} fires")

    # resume: skip completed
    pending = []
    done_rows = []
    for _, row in fires.iterrows():
        rc = CACHE_DIR / f"fwi_{int(row['fire_id'])}.json"
        if rc.exists():
            done_rows.append(json.loads(rc.read_text()))
        else:
            pending.append(row)

    print(f"Already done: {len(done_rows)}; pending: {len(pending)}")

    results = list(done_rows)
    if pending:
        with ThreadPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(process_one, row): int(row["fire_id"]) for row in pending}
            for i, fut in enumerate(as_completed(futs), 1):
                fid = futs[fut]
                try:
                    results.append(fut.result())
                    if i % 25 == 0 or i == len(futs):
                        print(f"  processed {i}/{len(futs)} (last fire_id={fid})")
                        pd.DataFrame(results).to_csv(OUT_CSV, index=False)
                except Exception as e:
                    print(f"  FAILED fire_id={fid}: {e}")

    out = pd.DataFrame(results).drop_duplicates("fire_id")
    out.to_csv(OUT_CSV, index=False)
    print(f"Saved {OUT_CSV} n={len(out)}")

    # join back to full CalFire (including missing IDate as NA FWI)
    full = gpd.read_file(CALFIRE)
    if "fid" in full.columns:
        full["fire_id"] = full["fid"].astype(int)
    else:
        full["fire_id"] = np.arange(len(full))
    merged = full.merge(out[["fire_id", "FFMC_4d", "DMC_4d", "DC_4d", "n_days"]], on="fire_id", how="left")
    # avoid fid conflict in gpkg
    if "fid" in merged.columns:
        merged = merged.rename(columns={"fid": "fire_fid"})
    merged.to_file(OUT_GPKG, driver="GPKG")
    print(f"Saved {OUT_GPKG}")
    print(out[["FFMC_4d", "DMC_4d", "DC_4d", "n_days"]].describe())


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--limit", type=int, default=None, help="Process only first N fires (test)")
    p.add_argument("--workers", type=int, default=MAX_WORKERS)
    args = p.parse_args()
    main(limit=args.limit, workers=args.workers)

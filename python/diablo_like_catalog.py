#!/usr/bin/env python3
"""
Daily Diablo-like catalog from gridMET, following Smith et al. (2018) and
Liu et al. (2021) as closely as daily data allow.

Papers (hourly RAWS / reanalysis):
  Smith: any NoBA station; wind 315–135°; vs > 11.17 m/s; RH < 30%; ≥3 h
  Liu:   any of two lee-side SFBA grids/stations (Knoxville Creek, Hawkeye);
         wind 315–90°; RH < 30%; vs > 8 (RAWS) / 7 (NARR) / 4 (ERA5) m/s; ≥6 h

Daily gridMET analog (this catalog):
  nearest gridMET cells to Knoxville Creek and Hawkeye;
  wind from 315–90°; vs ≥ 5 m/s; rmin < 30%; rmax < 60%;
  a day is Diablo-like if EITHER cell meets the criteria.

rmax is the daily stand-in for multi-hour dry air (Smith: both min and max RH
are depressed on Diablo days). vs ≥ 5 sits between Liu's ERA5 (4) and NARR (7)
thresholds because gridMET is a daily mean, not an hourly peak.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import xarray as xr

BASE = Path(__file__).resolve().parents[1]
VS_DIR = BASE / "dataSrc/climate/gridMET/wind speed"
TH_DIR = BASE / "dataSrc/climate/gridMET/wind direction"
RH_CACHE = BASE / "dataPrc/gridmet_diablo_rh"
DODS_RH = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET"
DAY_CSV = BASE / "dataPrc/diablo_like_days.csv"

# Liu et al. (2021) SFBA stations
STATIONS = {
    "Knoxville_Creek": (38.86, -122.42),
    "Hawkeye": (38.82, -122.51),
}
DIABLO_LAT = slice(39.6, 37.6)
DIABLO_LON = slice(-123.4, -121.7)
VS_MIN = 5.0
RMIN_MAX = 30.0
RMAX_MAX = 60.0
PROXY = "liu2grid_vs5_rmin30_rmax60_dir315-90"
YEAR0, YEAR1 = 1990, 2025


def _align_time(da: xr.DataArray, like: xr.DataArray) -> xr.DataArray:
    t_like = [d for d in like.dims if d not in ("lat", "lon")][0]
    t_da = [d for d in da.dims if d not in ("lat", "lon")][0]
    if t_da != t_like:
        da = da.rename({t_da: t_like})
    return da.reindex({t_like: like[t_like]}, method="nearest")


def _load_rh(var: str, year: int, like: xr.DataArray) -> xr.DataArray:
    RH_CACHE.mkdir(parents=True, exist_ok=True)
    p = RH_CACHE / f"{var}_{year}_nobay.nc"
    if p.exists():
        da = xr.open_dataarray(p).load()
    else:
        url = f"{DODS_RH}/{var}/{var}_{year}.nc"
        da = (
            xr.open_dataset(url)["relative_humidity"]
            .sel(lat=DIABLO_LAT, lon=DIABLO_LON)
            .load()
        )
        da.to_netcdf(p)
        print(f"    cached {p.name}")
    return _align_time(da, like.sel(lat=DIABLO_LAT, lon=DIABLO_LON))


def _station_frame(year: int) -> list[pd.DataFrame]:
    vs = xr.open_dataset(VS_DIR / f"vs_{year}.nc")["wind_speed"]
    th = xr.open_dataset(TH_DIR / f"th_{year}.nc")["wind_from_direction"]
    rmin = _load_rh("rmin", year, vs)
    rmax = _load_rh("rmax", year, vs)
    frames = []
    for lat, lon in STATIONS.values():
        df = pd.DataFrame(
            {
                "vs": vs.sel(lat=lat, lon=lon, method="nearest").to_pandas(),
                "th": th.sel(lat=lat, lon=lon, method="nearest").to_pandas(),
                "rmin": rmin.sel(lat=lat, lon=lon, method="nearest").to_pandas(),
                "rmax": rmax.sel(lat=lat, lon=lon, method="nearest").to_pandas(),
            }
        )
        df.index = pd.to_datetime(df.index).normalize()
        frames.append(df)
    return frames


def diablo_days_for_year(year: int) -> pd.Series:
    vs_p = VS_DIR / f"vs_{year}.nc"
    th_p = TH_DIR / f"th_{year}.nc"
    if not vs_p.exists() or not th_p.exists():
        return pd.Series(dtype=bool)
    bits = []
    for df in _station_frame(year):
        ne = (df["th"] >= 315) | (df["th"] <= 90)
        bits.append(
            ne & (df["vs"] >= VS_MIN) & (df["rmin"] < RMIN_MAX) & (df["rmax"] < RMAX_MAX)
        )
    is_day = bits[0] | bits[1]
    is_day.name = "is_Diablo"
    return is_day.astype(bool)


def load_diablo_catalog(force: bool = False) -> pd.Series:
    if DAY_CSV.exists() and not force:
        raw = pd.read_csv(DAY_CSV, parse_dates=["date"])
        if "proxy" in raw.columns and raw["proxy"].astype(str).eq(PROXY).all():
            s = raw.set_index("date")["is_Diablo"].astype(bool)
            s.index = pd.to_datetime(s.index).normalize()
            return s
    parts = []
    for year in range(YEAR0, YEAR1 + 1):
        d = diablo_days_for_year(year)
        if d.empty:
            print(f"  skip {year}: missing gridMET")
            continue
        parts.append(d)
        print(f"  Diablo-like {year}: {int(d.sum())}")
    cat = pd.concat(parts).sort_index()
    cat = cat[~cat.index.duplicated(keep="first")]
    pd.DataFrame(
        {"date": cat.index, "is_Diablo": cat.to_numpy(), "proxy": PROXY}
    ).to_csv(DAY_CSV, index=False)
    print("Saved", DAY_CSV)
    return cat


def monthly_mean_days(year0: int = 1990, year1: int = 2018) -> pd.Series:
    cat = load_diablo_catalog()
    cat = cat[(cat.index.year >= year0) & (cat.index.year <= year1)]
    n_years = year1 - year0 + 1
    counts = cat.groupby(cat.index.month).sum()
    return pd.Series({m: counts.get(m, 0) / n_years for m in range(1, 13)})

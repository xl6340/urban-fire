#!/usr/bin/env python3
"""
Hourly RAWS Diablo catalog following Smith et al. (2018) for Northern CA.

Smith NoBA stations (any one station qualifies the region-day):
  KNXC1  Knoxville Creek
  HWKC1  Hawkeye
  WISC1  County Line (Wilbur Springs)
  COWC1  Lyons Valley / Cow Mountain Ridge (same lat/lon as Smith Table 1)
  MASC1  Mendocino Pass
  EPKC1  Eagle Peak

Smith hourly rules:
  vs > 11.17 m/s (25 mph); wind 315–135°; RH < 30%; ≥3 consecutive hours.
A calendar day (America/Los_Angeles) is Diablo if any NoBA station has such a run.

Source: Iowa Environmental Mesonet HADS (public; first complete year 2011).
Smith's MesoWest record began in 1999; IEM does not go back that far.
HADS wind is mph; converted to m/s for the threshold.
"""

from __future__ import annotations

import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
RAW_DIR = BASE / "dataPrc/RAWS_HADS"
DAY_CSV = BASE / "dataPrc/diablo_raws_days.csv"
IEM = "https://mesonet.agron.iastate.edu/cgi-bin/request/hads.py"

STATIONS = {
    "KNXC1": "Knoxville_Creek",
    "HWKC1": "Hawkeye",
    "WISC1": "County_Line",
    "COWC1": "Lyons_Valley",
    "MASC1": "Mendocino_Pass",
    "EPKC1": "Eagle_Peak",
}
TZ = "America/Los_Angeles"
MPH_TO_MS = 0.44704
VS_MIN = 11.17  # Smith: vs > 11.17 m/s (25 mph)
RH_MAX = 30.0
DIR0, DIR1 = 315.0, 135.0
MIN_HOURS = 3
PROXY = "smith_noba_raws_vs11.17_rh30_dir315-135_h3"
YEAR0, YEAR1 = 2011, 2025  # first complete IEM year is 2011


def _iem_url(stid: str, year: int) -> str:
    q = urllib.parse.urlencode(
        {
            "network": "CA_DCP",
            "stations": stid,
            "year1": year,
            "month1": 1,
            "day1": 1,
            "year2": year,
            "month2": 12,
            "day2": 31,
            "delim": "comma",
            "what": "dl",
        }
    )
    return f"{IEM}?{q}"


def download_station_year(stid: str, year: int, force: bool = False) -> Path:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    path = RAW_DIR / f"{stid}_{year}.csv"
    if path.exists() and path.stat().st_size > 500 and not force:
        return path
    url = _iem_url(stid, year)
    print(f"  download {stid} {year} …")
    req = urllib.request.Request(url, headers={"User-Agent": "CalFire-Diablo/1.0"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        data = resp.read()
    if data[:1] in (b"{", b"<") or b"string_pattern_mismatch" in data[:200]:
        raise RuntimeError(f"IEM error for {stid} {year}: {data[:300]!r}")
    path.write_bytes(data)
    return path


def download_all(year0: int = YEAR0, year1: int = YEAR1, force: bool = False) -> None:
    for year in range(year0, year1 + 1):
        for stid in STATIONS:
            download_station_year(stid, year, force=force)
            time.sleep(0.3)


def _load_hourly(stid: str, year: int) -> pd.DataFrame:
    path = RAW_DIR / f"{stid}_{year}.csv"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path, low_memory=False)
    if df.empty or "utc_valid" not in df.columns:
        return pd.DataFrame()
    need = ["USIRGZ", "UDIRGZ", "XRIRGZ"]
    if any(c not in df.columns for c in need):
        return pd.DataFrame()
    out = pd.DataFrame(
        {
            "utc": pd.to_datetime(df["utc_valid"], utc=True),
            "vs_mph": pd.to_numeric(df["USIRGZ"], errors="coerce"),
            "wd": pd.to_numeric(df["UDIRGZ"], errors="coerce"),
            "rh": pd.to_numeric(df["XRIRGZ"], errors="coerce"),
        }
    )
    out = out.dropna(subset=["utc", "vs_mph", "wd", "rh"])
    if out.empty:
        return pd.DataFrame()
    out["vs"] = out["vs_mph"] * MPH_TO_MS
    out["station"] = stid
    local = out["utc"].dt.tz_convert(TZ)
    out["local"] = local.dt.tz_localize(None)
    out["hour"] = out["local"].dt.floor("h")
    out["date"] = out["hour"].dt.normalize()
    out = out.sort_values("local").drop_duplicates("hour", keep="last")
    return out.reset_index(drop=True)


def _qualifying(df: pd.DataFrame) -> pd.Series:
    ne = (df["wd"] >= DIR0) | (df["wd"] <= DIR1)
    return (df["vs"] > VS_MIN) & (df["rh"] < RH_MAX) & ne


def _event_hours(df: pd.DataFrame) -> pd.Series:
    """Hours that sit inside a ≥MIN_HOURS consecutive qualifying run."""
    if df.empty:
        return pd.Series(dtype=bool)
    s = df.set_index("hour").sort_index()
    # Insert missing hours so a data gap breaks the streak.
    full = pd.date_range(s.index.min(), s.index.max(), freq="h")
    ok = _qualifying(s).reindex(full, fill_value=False)
    grp = (ok != ok.shift(fill_value=False)).cumsum()
    runlen = ok.groupby(grp).transform("size")
    return (ok & (runlen >= MIN_HOURS)).rename("event")


def diablo_days_for_year(year: int) -> pd.Series:
    bits = []
    for stid in STATIONS:
        df = _load_hourly(stid, year)
        if df.empty:
            continue
        ev = _event_hours(df)
        if ev.empty:
            continue
        day = ev.groupby(ev.index.normalize()).any()
        bits.append(day.astype(bool))
    idx = pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")
    if not bits:
        return pd.Series(False, index=idx, name="is_Diablo")
    combined = pd.concat(bits, axis=1).any(axis=1)
    out = combined.reindex(idx, fill_value=False).astype(bool)
    out.name = "is_Diablo"
    return out


def load_diablo_catalog(force: bool = False, year0: int = YEAR0, year1: int = YEAR1) -> pd.Series:
    if DAY_CSV.exists() and not force:
        raw = pd.read_csv(DAY_CSV, parse_dates=["date"])
        if "proxy" in raw.columns and raw["proxy"].astype(str).eq(PROXY).all():
            s = raw.set_index("date")["is_Diablo"].astype(bool)
            s.index = pd.to_datetime(s.index).normalize()
            return s
    download_all(year0, year1, force=False)
    parts = []
    for year in range(year0, year1 + 1):
        d = diablo_days_for_year(year)
        print(f"  RAWS Diablo {year}: {int(d.sum())}")
        parts.append(d)
    cat = pd.concat(parts).sort_index()
    cat = cat[~cat.index.duplicated(keep="first")]
    pd.DataFrame(
        {"date": cat.index, "is_Diablo": cat.to_numpy(), "proxy": PROXY}
    ).to_csv(DAY_CSV, index=False)
    print("Saved", DAY_CSV)
    return cat


def monthly_mean_days(year0: int = YEAR0, year1: int = 2018) -> pd.Series:
    cat = load_diablo_catalog()
    cat = cat[(cat.index.year >= year0) & (cat.index.year <= year1)]
    n_years = year1 - year0 + 1
    counts = cat.groupby(cat.index.month).sum()
    return pd.Series({m: counts.get(m, 0) / n_years for m in range(1, 13)})


def _known_event_check() -> None:
    """Sanity: Tubbs 8 Oct 2017 and Camp 8 Nov 2018 should be True."""
    cat = load_diablo_catalog()
    for d in ["2017-10-08", "2017-10-09", "2018-11-08"]:
        ts = pd.Timestamp(d)
        print(f"  {d}: {bool(cat.loc[ts]) if ts in cat.index else 'missing'}")


if __name__ == "__main__":
    load_diablo_catalog(force=True)
    print("Monthly mean days 2011–2018:")
    m = monthly_mean_days(2011, 2018)
    print(" ".join(f"{m.loc[i]:4.2f}" for i in range(1, 13)))
    print("Known events:")
    _known_event_check()

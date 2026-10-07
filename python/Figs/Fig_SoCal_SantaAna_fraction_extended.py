#!/usr/bin/env python3
"""
Extend SoCal Santa Ana fire flags past 2018 with a gridMET proxy
(Keeley et al. 2024 style: Rolinski-like northeasterly + wind speed),
calibrated against R1D-SAWRI on 1990–2018 SoCal fires.

Santa Ana–affected definition (fire duration):
  A fire is SAW-affected if ANY day in [IDate, FDate] is a Santa Ana day.
  Missing FDate is filled with IDate (ignition-day only).
  Day-level SAW:
    - date ≤ 2018-12-31: official R1D-SAWRI > 1
    - date ≥ 2019-01-01: Diablo/SAW-like proxy
          wind_from_direction (fire-level vd) in [315, 360] ∪ [0, 90]
          AND daily vs ≥ vs_thr (calibrated F1 vs SAWRI>1 on ≤2018 ignition days)
          AND daily rmin < 30%
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import f1_score

BASE = Path(__file__).resolve().parents[2]
SAWRI_PATH = BASE / "dataSrc/climate/SAWRI/R1DSAWRI_010148_123118_red.txt"
WX_DIR = BASE / "dataPrc/fwi_cache"
OUT_FLAGS = BASE / "dataPrc/CalFire_SoCal_SAW_extended_flags.csv"
OUT_SUM = BASE / "dataPrc/CalFire_SoCal_SAW_extended_summary.csv"
OUT_PNG = BASE / "Fig/Fig_SoCal_SantaAna_fraction.png"
OUT_PDF = BASE / "Fig/Fig_SoCal_SantaAna_fraction.pdf"

C_WUI = np.array([216, 118, 89]) / 255
C_WILD = np.array([41, 157, 143]) / 255
RH_MAX = 30.0
SAWRI_END = pd.Timestamp("2018-12-31")


def load_sawri(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df.columns = ["dt", "SAWRI"]
    parsed = pd.to_datetime(df["dt"], format="%m/%d/%y")
    yy = parsed.dt.year % 100
    year = np.where(yy >= 48, 1900 + yy, 2000 + yy)
    dt = pd.to_datetime({"year": year, "month": parsed.dt.month, "day": parsed.dt.day})
    s = pd.Series(df["SAWRI"].astype(float).to_numpy(), index=dt, name="SAWRI")
    return s[~s.index.duplicated(keep="first")].sort_index()


def is_ne(vd) -> bool:
    if not np.isfinite(vd):
        return False
    return (vd >= 315.0) or (vd <= 90.0)


def load_wx(fid) -> pd.DataFrame | None:
    p = WX_DIR / f"wx_{int(fid)}.csv"
    if not p.exists():
        return None
    wx = pd.read_csv(p, parse_dates=["date"])
    wx["date"] = wx["date"].dt.normalize()
    return wx


def fire_window(idate, fdate):
    if pd.isna(idate):
        return None, None
    idate = pd.Timestamp(idate).normalize()
    if pd.isna(fdate):
        fdate = idate
    else:
        fdate = pd.Timestamp(fdate).normalize()
        if fdate < idate:
            fdate = idate
    return idate, fdate


def load_rmin_ign(fid, idate) -> float:
    wx = load_wx(fid)
    if wx is None or pd.isna(idate):
        return np.nan
    row = wx.loc[wx["date"] == pd.Timestamp(idate).normalize()]
    if row.empty:
        return np.nan
    return float(row["rmin"].iloc[0])


def calibrate_vs_thr(cal: pd.DataFrame) -> float:
    """Pick vs threshold maximizing F1 against ignition-day SAWRI>1 on ≤2018 fires."""
    y = cal["SAW_official_ign"].astype(bool).to_numpy()
    best_thr, best_f1 = 4.5, -1.0
    for thr in np.arange(2.0, 8.1, 0.25):
        pred = (
            cal["ne"].to_numpy()
            & (cal["vs_ign"].to_numpy() >= thr)
            & (cal["rmin_ign"].to_numpy() < RH_MAX)
        )
        if pred.sum() == 0:
            continue
        f1 = f1_score(y, pred, zero_division=0)
        if f1 > best_f1:
            best_f1, best_thr = f1, thr
    print(f"Calibrated vs_thr={best_thr:.2f} m/s (F1={best_f1:.3f} vs SAWRI>1 on ≤2018)")
    pred = cal["ne"] & (cal["vs_ign"] >= best_thr) & (cal["rmin_ign"] < RH_MAX)
    agree = (pred == cal["SAW_official_ign"]).mean()
    print(
        f"  agreement={100*agree:.1f}% | "
        f"proxy SAW rate={100*pred.mean():.1f}% | "
        f"official SAW rate={100*cal['SAW_official_ign'].mean():.1f}%"
    )
    return float(best_thr)


def saw_any_day(idate, fdate, vd, sawri, vs_thr, wx):
    """True if any day in [IDate, FDate] is Santa Ana (SAWRI>1 or post-2018 proxy)."""
    idate, fdate = fire_window(idate, fdate)
    if idate is None:
        return False, np.nan, "no_IDate"

    days = pd.date_range(idate, fdate, freq="D")
    for d in days:
        if d <= SAWRI_END:
            v = sawri.get(d, np.nan)
            if np.isfinite(v) and v > 1:
                return True, float(v), "SAWRI_anyday"

    post = days[days > SAWRI_END]
    if len(post) == 0:
        return False, float(sawri.get(idate, np.nan)) if idate <= SAWRI_END else np.nan, "SAWRI_window"
    if not is_ne(vd) or wx is None:
        return False, np.nan, "gridMET_proxy_window"
    sub = wx[wx["date"].isin(post)]
    if sub.empty:
        return False, np.nan, "gridMET_proxy_window"
    hit = (sub["vs"] >= vs_thr) & (sub["rmin"] < RH_MAX)
    if hit.any():
        return True, np.nan, "gridMET_proxy_anyday"
    return False, np.nan, "gridMET_proxy_window"


def main():
    sawri = load_sawri(SAWRI_PATH)
    g = gpd.read_file(BASE / "dataPrc/firePrmt/CalFire.shp")
    g = g[g["SoCal"] == 1].copy()
    g["IDate"] = pd.to_datetime(g["IDate"]).dt.normalize()
    g["FDate"] = pd.to_datetime(g["FDate"]).dt.normalize()
    size_col = "size_km2" if "size_km2" in g.columns else "size"

    print(f"Loading ignition-day rmin for {len(g)} SoCal fires …")
    g["rmin_ign"] = [load_rmin_ign(r.fid, r.IDate) for r in g.itertuples()]
    g["vs_ign"] = g["vs"]
    g["ne"] = g["vd"].map(is_ne)
    g["SAWRI_ign"] = g["IDate"].map(sawri)
    g["era_official"] = g["IDate"].notna() & (g["IDate"] <= SAWRI_END)
    g["SAW_official_ign"] = g["era_official"] & (g["SAWRI_ign"] > 1)

    cal = g[
        g["era_official"] & g["vs_ign"].notna() & g["vd"].notna() & g["rmin_ign"].notna()
    ].copy()
    vs_thr = calibrate_vs_thr(cal)

    print("Flagging duration-based Santa Ana …")
    saw_any, sawri_val, src = [], [], []
    wx_cache = {}
    for r in g.itertuples():
        fid = int(r.fid)
        if fid not in wx_cache:
            wx_cache[fid] = load_wx(fid)
        ok, sval, ssrc = saw_any_day(r.IDate, r.FDate, r.vd, sawri, vs_thr, wx_cache[fid])
        saw_any.append(ok)
        sawri_val.append(sval)
        src.append(ssrc)

    g["SAW"] = saw_any
    g["SAWRI"] = sawri_val
    g["SAW_source"] = src
    g["SAW_official"] = g["SAW_official_ign"]  # ignition-day reference
    g["SAW_proxy"] = (
        (~g["era_official"])
        & g["ne"]
        & (g["vs_ign"] >= vs_thr)
        & (g["rmin_ign"] < RH_MAX)
        & g["vs_ign"].notna()
        & g["rmin_ign"].notna()
    )
    g["in_analysis"] = g["era_official"] | (
        (~g["era_official"]) & g["vs_ign"].notna() & g["vd"].notna() & g["rmin_ign"].notna()
    )
    wins = [fire_window(i, f) for i, f in zip(g["IDate"], g["FDate"])]
    g["n_fire_days"] = [
        (1 + (b - a).days) if a is not None else np.nan for a, b in wins
    ]

    flags = pd.DataFrame(
        {
            "FireType": g["FireType"].values,
            "fid": g["fid"].values,
            "IDate": g["IDate"].values,
            "FDate": g["FDate"].values,
            "n_fire_days": g["n_fire_days"].values,
            "year": g["year"].values,
            "size": g[size_col].values,
            "vs": g["vs"].values,
            "vd": g["vd"].values,
            "rmin": g["rmin_ign"].values,
            "SAWRI": g["SAWRI"].values,
            "SAW_official": g["SAW_official"].values,
            "SAW_proxy": g["SAW_proxy"].values,
            "SAW": g["SAW"].values,
            "SAW_source": g["SAW_source"].values,
            "in_analysis": g["in_analysis"].values,
            "vs_thr": vs_thr,
        }
    )
    flags.to_csv(OUT_FLAGS, index=False)
    print("Saved", OUT_FLAGS)

    sub = flags[flags["in_analysis"]].copy()
    rows = []
    for ft, color in [("WUI", C_WUI), ("Wildland", C_WILD)]:
        gg = sub[sub["FireType"] == ft]
        p90 = gg["size"].quantile(0.9)
        large = gg[gg["size"] >= p90]
        pre = gg[gg["IDate"] <= SAWRI_END]
        post = gg[gg["IDate"] > SAWRI_END]
        rows.append(
            {
                "FireType": ft,
                "color": color,
                "n": len(gg),
                "n_SAW": int(gg["SAW"].sum()),
                "frac": gg["SAW"].mean(),
                "n_large": len(large),
                "n_large_SAW": int(large["SAW"].sum()),
                "frac_large": large["SAW"].mean() if len(large) else np.nan,
                "n_pre2019": len(pre),
                "frac_pre2019": pre["SAW"].mean() if len(pre) else np.nan,
                "n_post2018": len(post),
                "frac_post2018": post["SAW"].mean() if len(post) else np.nan,
                "n_post_SAW": int(post["SAW"].sum()) if len(post) else 0,
            }
        )
        print(
            f"{ft}: all {100*gg.SAW.mean():.1f}% ({int(gg.SAW.sum())}/{len(gg)}) | "
            f"large {100*large.SAW.mean():.1f}% | "
            f"≤2018 {100*pre.SAW.mean():.1f}% | ≥2019 {100*post.SAW.mean():.1f}% "
            f"({int(post.SAW.sum())}/{len(post)})"
        )
    summary = pd.DataFrame(rows)
    summary.drop(columns=["color"]).to_csv(OUT_SUM, index=False)
    print("Saved", OUT_SUM)

    # Figure
    fig, ax = plt.subplots(figsize=(5.2, 4.2), facecolor="white")
    x = np.arange(2)
    w = 0.34
    colors = [summary.loc[0, "color"], summary.loc[1, "color"]]
    bars1 = ax.bar(
        x - w / 2,
        summary["frac"] * 100,
        width=w,
        color=colors,
        edgecolor="0.25",
        lw=0.6,
        label="All SoCal fires",
    )
    bars2 = ax.bar(
        x + w / 2,
        summary["frac_large"] * 100,
        width=w,
        color=colors,
        edgecolor="0.25",
        lw=0.6,
        hatch="///",
        alpha=0.75,
        label="Large fires (≥ size p90)",
    )
    for b, n_s, n in zip(bars1, summary["n_SAW"], summary["n"]):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 1.0,
            f"{b.get_height():.1f}%\n({n_s}/{n})",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color="0.2",
        )
    for b, n_s, n in zip(bars2, summary["n_large_SAW"], summary["n_large"]):
        ax.text(
            b.get_x() + b.get_width() / 2,
            b.get_height() + 1.0,
            f"{b.get_height():.1f}%\n({n_s}/{n})",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color="0.2",
        )

    ymax = max(summary["frac"].max(), summary["frac_large"].max()) * 100
    ax.set_ylim(0, max(55, ymax * 1.3))
    ax.set_xticks(x)
    ax.set_xticklabels(["WUI", "Wildland"], fontsize=12)
    ax.set_ylabel("Santa Ana-affected fires (%)", fontsize=11)
    y0, y1 = int(pd.to_datetime(sub["IDate"]).dt.year.min()), int(
        pd.to_datetime(sub["IDate"]).dt.year.max()
    )
    ax.set_title(
        "Southern California: Santa Ana–affected fires\n"
        f"(any SAW day in IDate–FDate; SAWRI>1 through 2018; "
        f"gridMET proxy 2019–{y1})",
        fontsize=11,
        pad=8,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out")
    ax.legend(frameon=False, fontsize=9, loc="upper right")
    ax.text(
        0.0,
        -0.20,
        f"Proxy calibrated on ≤2018 ignition days (vs≥{vs_thr:.2f} m/s, RH<{RH_MAX:.0f}%, NE wind). "
        "Hatched = ≥ within-type size p90.",
        transform=ax.transAxes,
        fontsize=7.5,
        color="0.35",
        ha="left",
        va="top",
    )
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
    print("Saved", OUT_PNG)
    print("Saved", OUT_PDF)


if __name__ == "__main__":
    main()

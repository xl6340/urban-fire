#!/usr/bin/env python3
"""
Western U.S. MTBS WUI vs wildland size-distribution β.

California: Fig. 1C MTBS values (paper WUI / wildland flags).
Other states: WUI = perimeter intersects Radeloff 2020 WUI
              (WUIFLAG2020 = 1 intermix or 2 interface).
β = OLS slope of log10(PDF) vs log10(size) on log-spaced bins.
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyogrio import read_dataframe
from scipy import stats

BASE = Path(__file__).resolve().parents[1]
OUT_CSV = BASE / "dataPrc/MTBS_westcoast_WUI_wildland_beta.csv"
OUT_PNG = BASE / "Fig/Fig_westcoast_WUI_wildland_beta.png"
OUT_PDF = BASE / "Fig/Fig_westcoast_WUI_wildland_beta.pdf"

C_WUI = np.array([216, 118, 89]) / 255
C_WILD = np.array([41, 157, 143]) / 255
EDGES = 10.0 ** np.arange(-1.0, 5.05, 0.05)

# Pacific, then interior West (no pooled "West Coast" bar)
STATES = ["CA", "OR", "WA", "NV", "ID", "AZ", "UT", "CO", "NM", "MT", "WY"]
STATE_LABELS = {
    "CA": "CA",
    "OR": "OR",
    "WA": "WA",
    "NV": "NV",
    "ID": "ID",
    "AZ": "AZ",
    "UT": "UT",
    "CO": "CO",
    "NM": "NM",
    "MT": "MT",
    "WY": "WY",
}


def apply_fig1c_california(df):
    """Replace CA WUI/Wildland with Fig. 1C MTBS (4Srcs-*.csv), not WUI-2020 overlay."""
    wui = pd.read_csv(BASE / "dataFig/beta/4Srcs-Urban-edge.csv")
    wild = pd.read_csv(BASE / "dataFig/beta/4Srcs-Wildland.csv")
    wui = wui.set_index("Row").loc["MTBS"]
    wild = wild.set_index("Row").loc["MTBS"]
    out = df.copy()
    for ft, src in [("WUI", wui), ("Wildland", wild)]:
        m = (out["region"] == "CA") & (out["FireType"] == ft)
        out.loc[m, "n"] = int(src["count"])
        out.loc[m, "area"] = float(src["area"])
        out.loc[m, "beta"] = float(src["beta"])
        out.loc[m, "beta_se"] = float(src["betaErr"])
        out.loc[m, "r2"] = float(src["R2"])
        out.loc[m, ["median", "p90"]] = np.nan
    return out
    sizes = np.asarray(sizes, float)
    sizes = sizes[np.isfinite(sizes) & (sizes > 0)]
    out = dict(
        n=int(sizes.size),
        area=np.nan,
        median=np.nan,
        p90=np.nan,
        beta=np.nan,
        beta_se=np.nan,
        r2=np.nan,
    )
    if sizes.size < 10:
        return out
    counts, be = np.histogram(sizes, bins=EDGES)
    tot = counts.sum()
    if tot == 0:
        return out
    pdf = counts / (tot * np.diff(be))
    ctr = np.sqrt(be[:-1] * be[1:])
    m = pdf > 0
    res = stats.linregress(np.log10(ctr[m]), np.log10(pdf[m]))
    out.update(
        area=float(sizes.sum()),
        median=float(np.median(sizes)),
        p90=float(np.quantile(sizes, 0.9)),
        beta=float(res.slope),
        beta_se=float(res.stderr),
        r2=float(res.rvalue**2),
    )
    return out


def load_mtbs():
    like = " OR ".join(f"Event_ID LIKE '{st}%'" for st in STATES)
    fires = read_dataframe(
        BASE / "dataSrc/firePrmt/MTBS/mtbs_perims_DD.shp",
        columns=["Event_ID", "Incid_Type", "BurnBndAc", "Ig_Date"],
        where=like,
    )
    fires["state"] = fires["Event_ID"].astype(str).str[:2]
    fires["year"] = pd.to_datetime(fires["Ig_Date"], errors="coerce").dt.year
    fires["size"] = fires["BurnBndAc"] * 0.00404686
    fires = fires[
        fires.Incid_Type.eq("Wildfire")
        & fires.year.ge(1990)
        & fires.year.le(2024)
        & fires.state.isin(STATES)
    ].copy()
    return fires


def flag_wui(fires, wui, wui_crs):
    if fires.empty:
        return set()
    f = fires.to_crs(wui_crs)
    if not isinstance(wui, gpd.GeoDataFrame):
        wui = gpd.GeoDataFrame(wui, geometry="geometry", crs=wui_crs)
    elif str(wui.crs) != str(f.crs):
        wui = wui.to_crs(f.crs)
    joined = gpd.sjoin(
        f[["Event_ID", "geometry"]],
        wui[["geometry"]],
        predicate="intersects",
        how="inner",
    )
    return set(joined["Event_ID"].unique())


def load_wui_layer(path, fallback_crs=5070):
    wui = read_dataframe(path)
    if wui.crs is None:
        wui = wui.set_crs(fallback_crs)
        print(f"  assigned EPSG:{fallback_crs} to {path.name}")
    print(f"  {path.name}: n={len(wui)} crs={wui.crs}")
    return wui


def main():
    print("Loading MTBS western wildfires …")
    fires = load_mtbs()
    print("  n =", len(fires), fires.groupby("state").size().reindex(STATES).to_dict())

    print("Loading WUI 2020 layers …")
    wui_ca = load_wui_layer(BASE / "dataPrc/WUI/WUI2020.gpkg")
    wui_nw = load_wui_layer(BASE / "dataPrc/WUI/WUI2020_OR_WA.gpkg")
    wui_in = load_wui_layer(BASE / "dataPrc/WUI/WUI2020_interior_west.gpkg")

    ids = set()
    ids |= flag_wui(fires[fires.state.eq("CA")], wui_ca, wui_ca.crs)
    ids |= flag_wui(fires[fires.state.isin(["OR", "WA"])], wui_nw, wui_nw.crs)
    ids |= flag_wui(
        fires[fires.state.isin(["NV", "ID", "AZ", "UT", "CO", "NM", "MT", "WY"])],
        wui_in,
        wui_in.crs,
    )
    fires["FireType"] = np.where(fires["Event_ID"].isin(ids), "WUI", "Wildland")
    print(fires.groupby(["state", "FireType"]).size().unstack(fill_value=0).reindex(STATES))

    rows = []
    for st in STATES:
        sub = fires[fires.state.eq(st)]
        for ft in ["WUI", "Wildland"]:
            r = ols_beta(sub.loc[sub.FireType.eq(ft), "size"])
            r.update(region=st, FireType=ft)
            rows.append(r)
        r = ols_beta(sub["size"])
        r.update(region=st, FireType="All")
        rows.append(r)

    df = pd.DataFrame(rows)[
        ["region", "FireType", "n", "area", "median", "p90", "beta", "beta_se", "r2"]
    ]
    df = apply_fig1c_california(df)
    df.to_csv(OUT_CSV, index=False)
    print(df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print("Saved", OUT_CSV)

    plot_figure(df)


def plot_figure(df):
    plot_df = df[df.FireType.isin(["WUI", "Wildland"])].copy()
    wui = plot_df[plot_df.FireType.eq("WUI")].set_index("region").loc[STATES]
    wild = plot_df[plot_df.FireType.eq("Wildland")].set_index("region").loc[STATES]
    y = np.arange(len(STATES))[::-1]

    fig, ax = plt.subplots(figsize=(5.4, 6.0), facecolor="white")
    for yi, st in zip(y, STATES):
        ax.plot(
            [wild.loc[st, "beta"], wui.loc[st, "beta"]],
            [yi, yi],
            color="0.82",
            lw=1.1,
            zorder=1,
            solid_capstyle="round",
        )

    ax.errorbar(
        wild["beta"].to_numpy(),
        y,
        xerr=wild["beta_se"].to_numpy(),
        fmt="o",
        color=C_WILD,
        markerfacecolor=C_WILD,
        markersize=7,
        capsize=2.5,
        lw=1.2,
        label="Wildland",
        zorder=3,
    )
    ax.errorbar(
        wui["beta"].to_numpy(),
        y,
        xerr=wui["beta_se"].to_numpy(),
        fmt="o",
        color=C_WUI,
        markerfacecolor=C_WUI,
        markersize=7,
        capsize=2.5,
        lw=1.2,
        label="WUI",
        zorder=4,
    )

    ax.set_yticks(y)
    ax.set_yticklabels(
        [f"{STATE_LABELS[st]}  {int(wui.loc[st, 'n'])}/{int(wild.loc[st, 'n'])}" for st in STATES],
        fontsize=9,
    )
    ax.set_xlabel(r"Size-distribution $\beta$")
    ax.set_ylabel("State  (n WUI / n wildland)")
    xmin = np.nanmin(plot_df["beta"] - plot_df["beta_se"]) - 0.08
    xmax = np.nanmax(plot_df["beta"] + plot_df["beta_se"]) + 0.08
    ax.set_xlim(min(xmin, -2.05), max(xmax, -0.95))
    ax.set_ylim(-0.6, len(STATES) - 0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        [handles[1], handles[0]],
        [labels[1], labels[0]],
        loc="lower center",
        bbox_to_anchor=(0.5, 1.02),
        ncol=2,
        frameon=False,
        columnspacing=1.6,
        handletextpad=0.4,
        borderaxespad=0.0,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(OUT_PDF, bbox_inches="tight", facecolor="white")
    print("Saved", OUT_PNG)


if __name__ == "__main__":
    main()

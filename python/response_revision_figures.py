#!/usr/bin/env python3
"""
Response-revision analyses (COMMSENV-26-2928-T / CAfire-Scaling-CEE.docx)

One entry point for all NEW figures/analyses added in the response letter,
ordered by first appearance of each new figure in the document.

Usage
-----
  python response_revision_figures.py           # run all
  python response_revision_figures.py list      # print index
  python response_revision_figures.py S1 S2     # run selected keys
  python response_revision_figures.py R1a AIC   # etc.

Figure index (response order)
-----------------------------
  01  S1    Fig. S1     R1C10   WUI-ignite / WUI-spread / Wildland power law + WUI frac
  02  R1a   Fig. R1     R1C10   WUI-spread ≥10% within-WUI threshold sensitivity
  03  S6    Fig. S6     R2C8    Forest FWI moisture codes (FFMC, DMC, DC) histograms
  04  R1b   Fig. R1     R2C8    β vs regional FWI (FFMC / DMC / DC)
  05  S2    Fig. S2     R3C4    Trunc. exp. / lognormal / Weibull (4 datasets)
  05b AIC   (table)     R3C4    AIC + likelihood-ratio / Vuong (+ power law)
  06  S10   Fig. S10    R3C9    Santa Ana day fraction + size PDFs (SoCal)
  07  S11   Fig. S11    R3C13   SG vegetation: NDVI / RAP NPP / ESA CCI AGB
  08  S9    Fig. S9     R3C16   Monthly burned area (human/natural) + lightning + SAW
  09  R1c   Fig. R1     R3C22   Western U.S. MTBS WUI vs wildland β by state

Existing modular sources (imported below; do not duplicate large libraries here):
  code/Figs/Fig_CalFire_IgnWUI_powerlaw.py
  code/Figs/Fig3B_Forest_FWI_dist.py
  code/fit_size_distributions_alt.py
  code/size_dist_AIC_LR_comparison.py
  code/Figs/Fig_SoCal_SAW_powerlaw.py
  code/Figs/Fig_SG_WUI_vs_wildland_veg_dist.py
  code/Figs/Fig_offshore_wind_monthly.py
  code/westcoast_wui_wildland_beta.py

Analyses that previously lived only as one-off chat scripts are inlined here:
  R1a (spread ≥10%), R1b (β vs regional FWI).
"""

from __future__ import annotations

import importlib.util
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
CODE = BASE / "python"
FIGS = CODE / "Figs"
OUT_FIG = BASE / "Fig"
OUT_PRC = BASE / "dataPrc"

sys.path.insert(0, str(CODE))
sys.path.insert(0, str(FIGS))


def _load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def _copy_alias(src: Path, *aliases: Path):
    if not src.exists():
        print(f"  [warn] missing {src}")
        return
    for dst in aliases:
        shutil.copy2(src, dst)
        print(f"  alias → {dst.relative_to(BASE)}")


# =============================================================================
# 01  Fig. S1  |  R1C10  |  IgnWUI power-law + within-WUI burned-area fraction
# =============================================================================
def run_S1():
    print("\n=== 01 Fig. S1 (R1C10): IgnWUI power law ===")
    mod = _load_module(FIGS / "Fig_CalFire_IgnWUI_powerlaw.py", "fig_ignwui_s1")
    mod.main()
    _copy_alias(
        OUT_FIG / "CalFire_IgnWUI_powerlaw.png",
        OUT_FIG / "FigS1.png",
    )
    pdf = OUT_FIG / "CalFire_IgnWUI_powerlaw.pdf"
    if pdf.exists():
        _copy_alias(pdf, OUT_FIG / "FigS1.pdf")


# =============================================================================
# 02  Fig. R1a  |  R1C10  |  WUI-spread ≥10% within-WUI threshold (response-only)
# =============================================================================
def run_R1a():
    """Baseline vs WUI-spread ≥10% WUI area; dual-axis β + area-wt inset."""
    print("\n=== 02 Fig. R1a (R1C10): WUI-spread ≥10% threshold ===")
    import geopandas as gpd
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from scipy import stats

    FIRE_GPKG = OUT_PRC / "CalFire_IgnWUI_strict.gpkg"
    FRAC_CSV = OUT_PRC / "CalFire_IgnWUI_within_fire_WUI_frac.csv"
    OUT_PNG = OUT_FIG / "CalFire_IgnWUI_powerlaw_spread_frac_thresh.png"
    OUT_PDF = OUT_FIG / "CalFire_IgnWUI_powerlaw_spread_frac_thresh.pdf"
    OUT_CSV = OUT_PRC / "CalFire_IgnWUI_powerlaw_spread_frac_thresh.csv"
    FIG1_EDGES = 10.0 ** np.arange(-1.0, 5.0 + 1e-12, 0.05)

    TYPE_ORDER = ["WUIIgnite", "WUISpread", "Wildland"]
    LABELS = {"WUIIgnite": "WUI-ignite", "WUISpread": "WUI-spread", "Wildland": "Wildland"}
    XTICK = {"WUIIgnite": "Ignition", "WUISpread": "Spread", "Wildland": "Wildland"}
    COLORS = {"WUIIgnite": "#D8765A", "WUISpread": "#E0A35A", "Wildland": "#299A8F"}
    FRAC = "wui_burn_frac"

    fire = gpd.read_file(FIRE_GPKG)
    fire = fire[fire["n_AtlasIgn"] > 0].copy()
    fracs = pd.read_csv(FRAC_CSV)[["fire_fid", FRAC]].drop_duplicates("fire_fid")
    df = fire.drop(columns="geometry").merge(fracs, on="fire_fid", how="left")

    def fit_powerlaw_ols_fig1(sizes, size_min=1.0, edges=FIG1_EDGES):
        sizes = pd.Series(sizes).dropna()
        sizes = sizes[sizes >= size_min]
        counts, _ = np.histogram(sizes.to_numpy(), bins=edges)
        density = counts / (len(sizes) * np.diff(edges))
        centers = np.sqrt(edges[:-1] * edges[1:])
        mask = counts > 0
        log_x = np.log10(centers[mask])
        log_y = np.log10(density[mask])
        slope, intercept, r, p, se = stats.linregress(log_x, log_y)
        return dict(
            beta=float(slope),
            intercept=float(intercept),
            se=float(se),
            r2=float(r**2),
            n=int(len(sizes)),
            log_x=log_x,
            log_y=log_y,
        )

    def bootstrap_beta(sizes, n_iter=800, size_min=1.0, seed=1):
        rng = np.random.default_rng(seed)
        s = pd.Series(sizes).dropna()
        s = s[s >= size_min].to_numpy()
        betas = [
            fit_powerlaw_ols_fig1(rng.choice(s, size=len(s), replace=True))["beta"]
            for _ in range(n_iter)
        ]
        b = np.asarray(betas)
        return float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))

    def area_weighted_frac(sub):
        return float((sub[FRAC] * sub["size_km2"]).sum() / sub["size_km2"].sum())

    thresholds = [None, 0.10]
    titles = ["Baseline (all WUI-spread)", r"WUI-spread $\geq$ 10% WUI area"]
    rows = []

    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), sharey=True, facecolor="w")
    fig.subplots_adjust(wspace=0.16, left=0.09, right=0.98, top=0.82, bottom=0.15)
    handles = [Line2D([0], [0], color=COLORS[n], lw=2.2, label=LABELS[n]) for n in TYPE_ORDER]
    handles += [
        Line2D(
            [0],
            [0],
            marker="o",
            color="0.35",
            markerfacecolor="0.35",
            markersize=5.5,
            linestyle="None",
            label=r"$\beta$",
        ),
        Patch(facecolor="0.65", edgecolor="0.45", alpha=0.45, label="Within-WUI fraction"),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=5,
        frameon=False,
        fontsize=8.5,
        bbox_to_anchor=(0.55, 1.03),
        handlelength=1.5,
        columnspacing=1.1,
    )

    for ax, thr, title in zip(axes, thresholds, titles):
        d = df.copy()
        n_drop = 0
        if thr is not None:
            drop_mask = (d.Type == "WUISpread") & (d[FRAC] < thr)
            n_drop = int(drop_mask.sum())
            d = d.loc[~drop_mask].copy()

        results, aw = {}, {}
        for name in TYPE_ORDER:
            sub = d.loc[d.Type == name]
            res = fit_powerlaw_ols_fig1(sub["size_km2"])
            lo, hi = bootstrap_beta(sub["size_km2"])
            res["ci"] = [lo, hi]
            results[name] = res
            aw[name] = area_weighted_frac(sub)
            c = COLORS[name]
            ax.scatter(
                10 ** res["log_x"],
                10 ** res["log_y"],
                color=c,
                alpha=0.55,
                s=12,
                zorder=3,
                lw=0,
            )
            x_fit = np.linspace(res["log_x"].min(), res["log_x"].max(), 80)
            ax.plot(
                10**x_fit,
                10 ** (res["beta"] * x_fit + res["intercept"]),
                color=c,
                lw=2.0,
            )
            rows.append(
                {
                    "threshold": 0.0 if thr is None else thr,
                    "Type": name,
                    "n": res["n"],
                    "beta": res["beta"],
                    "ci_low": lo,
                    "ci_high": hi,
                    "r2": res["r2"],
                    "se": res["se"],
                    "aw_wui_frac": aw[name],
                    "n_spread_dropped": n_drop if name == "WUISpread" else 0,
                }
            )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"Fire size (km$^2$)")
        ax.set_title(title, fontsize=10, pad=4)
        ax.tick_params(direction="out")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)

        axins = ax.inset_axes([0.58, 0.58, 0.40, 0.38])
        axins_r = axins.twinx()
        xs = np.arange(3)
        bar_h = [aw[n] * 100 for n in TYPE_ORDER]

        axins_r.bar(
            xs,
            bar_h,
            width=0.55,
            color=[COLORS[n] for n in TYPE_ORDER],
            alpha=0.35,
            edgecolor=[COLORS[n] for n in TYPE_ORDER],
            linewidth=0.8,
            zorder=1,
        )
        for x, h in zip(xs, bar_h):
            axins_r.text(
                x,
                h + max(bar_h + [12]) * 0.04,
                f"{h:.1f}%",
                ha="center",
                va="bottom",
                fontsize=5.2,
                color="0.35",
            )

        for i, name in enumerate(TYPE_ORDER):
            b = results[name]["beta"]
            lo, hi = results[name]["ci"]
            axins.errorbar(
                i,
                b,
                yerr=[[b - lo], [hi - b]],
                fmt="o",
                color=COLORS[name],
                markersize=4.5,
                capsize=2.0,
                lw=1.1,
                markeredgecolor="white",
                markeredgewidth=0.5,
                zorder=3,
            )
            axins.text(
                i,
                lo - 0.035,
                f"n={results[name]['n']}",
                ha="center",
                va="top",
                fontsize=4.6,
                color="0.4",
            )

        axins.axhline(-1.5, color="0.7", ls=":", lw=0.8, zorder=0)
        axins.set_xticks(xs)
        axins.set_xticklabels(
            [XTICK[n] for n in TYPE_ORDER],
            fontsize=5.2,
            rotation=45,
            ha="right",
            rotation_mode="anchor",
        )
        axins.set_ylabel(r"$\beta$", fontsize=6.2, color="0.2")
        axins_r.set_ylabel("Within-WUI area (%)", fontsize=5.4, color="0.35")
        axins.tick_params(axis="y", labelsize=5.2, length=2, colors="0.2")
        axins_r.tick_params(axis="y", labelsize=5.0, length=2, colors="0.4")
        axins.set_xlim(-0.55, 2.55)

        los = [results[n]["ci"][0] for n in TYPE_ORDER]
        his = [results[n]["ci"][1] for n in TYPE_ORDER]
        axins.set_ylim(min(los) - 0.18, max(his) + 0.08)
        axins_r.set_ylim(0, max(bar_h + [15]) * 1.45)

        axins.set_facecolor("0.97")
        for sp in axins.spines.values():
            sp.set_linewidth(0.6)
            sp.set_color("0.55")
        for sp in axins_r.spines.values():
            sp.set_linewidth(0.6)
            sp.set_color("0.55")

    axes[0].set_ylabel("Probability density")
    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
    fig.savefig(OUT_PNG, dpi=300, facecolor="w", bbox_inches="tight")
    fig.savefig(OUT_PDF, facecolor="w", bbox_inches="tight")
    plt.close(fig)
    print("Saved", OUT_PNG)
    _copy_alias(OUT_PNG, OUT_FIG / "FigR1_IgnWUI_spread10pct.png")


# =============================================================================
# 03  Fig. S6  |  R2C8  |  Forest FWI moisture-code distributions
# =============================================================================
def run_S6():
    print("\n=== 03 Fig. S6 (R2C8): Forest FFMC / DMC / DC ===")
    mod = _load_module(FIGS / "Fig3B_Forest_FWI_dist.py", "fig_fwi_s6")
    mod.main()
    _copy_alias(
        OUT_FIG / "Fig3B_Forest_FWI_dist.png",
        OUT_FIG / "FigS6.png",
    )
    pdf = OUT_FIG / "Fig3B_Forest_FWI_dist.pdf"
    if pdf.exists():
        _copy_alias(pdf, OUT_FIG / "FigS6.pdf")


# =============================================================================
# 04  Fig. R1b  |  R2C8  |  β vs regional FWI (May–Oct half-decade means)
# =============================================================================
def run_R1b():
    print("\n=== 04 Fig. R1b (R2C8): β vs regional FWI ===")
    import matplotlib.pyplot as plt
    from scipy import stats

    daily = pd.read_csv(OUT_PRC / "CalFire_FWI_regional_daily.csv", parse_dates=["date"])
    daily["year"] = daily["date"].dt.year
    daily["month"] = daily["date"].dt.month
    fs = daily[daily["month"].between(5, 10)].copy()

    def half_decade_label(y):
        return f"{(y // 5) * 5}s"

    fs["half_decade"] = fs["year"].map(half_decade_label)
    keep = ["1990s", "1995s", "2000s", "2005s", "2010s", "2015s", "2020s"]
    fs = fs[fs["half_decade"].isin(keep)]

    agg = (
        fs.groupby("half_decade")
        .agg(
            FFMC_WUI=("FFMC_wui", "mean"),
            FFMC_Wildland=("FFMC_wild", "mean"),
            DMC_WUI=("DMC_wui", "mean"),
            DMC_Wildland=("DMC_wild", "mean"),
            DC_WUI=("DC_wui", "mean"),
            DC_Wildland=("DC_wild", "mean"),
            n_days=("date", "count"),
        )
        .reindex(keep)
        .reset_index()
    )
    agg.to_csv(OUT_PRC / "CalFire_FWI_regional_halfdecade_MayOct.csv", index=False)
    agg.to_csv(OUT_PRC / "CalFire_FWI_regional_halfdecade.csv", index=False)

    beta = pd.read_excel(BASE / "dataFig/vpd/betaVPD.xlsx", sheet_name="arith")
    beta = beta[beta["Period"].isin(keep)].copy()
    err_u, err_w = "beta_err.1", "beta_err.2"
    df = beta.merge(agg.rename(columns={"half_decade": "Period"}), on="Period", how="inner")

    colorUrban = np.array([216, 118, 89]) / 255
    colorWild = np.array([41, 157, 143]) / 255

    def plot_beta_vs(ax, x_u, x_w, y_u, y_w, yerr_u, yerr_w, periods, xlabel, panel):
        ax.errorbar(
            x_u,
            y_u,
            yerr=yerr_u,
            fmt="o",
            color=colorUrban,
            mfc=colorUrban,
            ms=7,
            lw=1,
            capsize=0,
            linestyle="none",
            label="WUI",
            zorder=3,
        )
        ax.errorbar(
            x_w,
            y_w,
            yerr=yerr_w,
            fmt="^",
            color=colorWild,
            mfc=colorWild,
            ms=7,
            lw=1,
            capsize=0,
            linestyle="none",
            label="Wildland",
            zorder=3,
        )

        def add_fit(x, y, color, label_prefix, ytext):
            mask = np.isfinite(x) & np.isfinite(y)
            xx, yy = np.asarray(x)[mask], np.asarray(y)[mask]
            slope, intercept, r, p, se = stats.linregress(xx, yy)
            xfit = np.linspace(xx.min(), xx.max(), 100)
            ax.plot(xfit, slope * xfit + intercept, "-", color=color, lw=1.6, zorder=2)
            pstr = f"{p:.3f}" if p >= 0.001 else "<0.001"
            ax.text(
                0.04,
                ytext,
                f"{label_prefix}: r={r:.2f}, p={pstr}",
                transform=ax.transAxes,
                color=color,
                fontsize=11,
                ha="left",
                va="top",
            )
            return r, p

        ru, pu = add_fit(x_u, y_u, colorUrban, "WUI", 0.96)
        rw, pw = add_fit(x_w, y_w, colorWild, "Wildland", 0.86)
        for i, per in enumerate(periods):
            ax.text(x_u[i], y_u[i] - 0.035, str(per), color="0.55", fontsize=7.5, ha="center", va="top")
            ax.text(x_w[i], y_w[i] - 0.02, str(per), color="0.55", fontsize=7.5, ha="center", va="top")
        ax.set_xlabel(xlabel, fontsize=12)
        ax.set_ylabel(r"$\beta$ value", fontsize=12)
        ax.set_title(panel, loc="left", fontsize=12)
        ax.tick_params(labelsize=10)
        return (ru, pu), (rw, pw)

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0), facecolor="white")
    specs = [
        ("FFMC_WUI", "FFMC_Wildland", r"Regional FFMC", "(A)"),
        ("DMC_WUI", "DMC_Wildland", r"Regional DMC", "(B)"),
        ("DC_WUI", "DC_Wildland", r"Regional DC", "(C)"),
    ]
    stats_rows = []
    for ax, (cu, cw, xlab, lab) in zip(axes, specs):
        (ru, pu), (rw, pw) = plot_beta_vs(
            ax,
            df[cu].values,
            df[cw].values,
            df["beta_urban"].values,
            df["beta_wild"].values,
            df[err_u].values,
            df[err_w].values,
            df["Period"].tolist(),
            xlab,
            lab,
        )
        ax.set_ylim(-1.75, -1.05)
        stats_rows += [
            {"index": xlab.replace("Regional ", ""), "group": "WUI", "r": ru, "p": pu, "season": "May-Oct", "n": len(df)},
            {"index": xlab.replace("Regional ", ""), "group": "Wildland", "r": rw, "p": pw, "season": "May-Oct", "n": len(df)},
        ]
    axes[0].legend(frameon=False, loc="lower right", fontsize=9)
    fig.tight_layout()
    for pth in [OUT_FIG / "Fig_beta_vs_regional_FWI.png", OUT_PRC / "CalFire_beta_vs_regional_FWI.png"]:
        fig.savefig(pth, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 4.3), facecolor="white")
    plot_beta_vs(
        ax,
        df["DC_WUI"].values,
        df["DC_Wildland"].values,
        df["beta_urban"].values,
        df["beta_wild"].values,
        df[err_u].values,
        df[err_w].values,
        df["Period"].tolist(),
        r"Regional DC",
        "(B)",
    )
    ax.set_ylim(-1.75, -1.05)
    ax.legend(frameon=False, loc="lower right", fontsize=10)
    fig.tight_layout()
    for pth in [OUT_FIG / "Fig_beta_vs_regional_DC.png", OUT_PRC / "CalFire_beta_vs_regional_DC.png"]:
        fig.savefig(pth, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    out_tbl = df[
        [
            "Period",
            "beta_urban",
            err_u,
            "beta_wild",
            err_w,
            "FFMC_WUI",
            "FFMC_Wildland",
            "DMC_WUI",
            "DMC_Wildland",
            "DC_WUI",
            "DC_Wildland",
        ]
    ].copy()
    out_tbl = out_tbl.rename(columns={err_u: "beta_err_WUI", err_w: "beta_err_Wildland"})
    out_tbl.to_csv(OUT_PRC / "CalFire_beta_vs_regional_FWI.csv", index=False)
    pd.DataFrame(stats_rows).to_csv(OUT_PRC / "CalFire_beta_vs_regional_FWI_stats.csv", index=False)
    print("Saved Fig_beta_vs_regional_FWI.png")
    _copy_alias(OUT_FIG / "Fig_beta_vs_regional_FWI.png", OUT_FIG / "FigR1_beta_vs_regional_FWI.png")


# =============================================================================
# 05  Fig. S2  |  R3C4  |  Trunc. exp. / lognormal / Weibull (4 sources)
# =============================================================================
def run_S2():
    print("\n=== 05 Fig. S2 (R3C4): alt. distributional forms ===")
    import fit_size_distributions_alt as alt

    data = alt.load_groups()
    alt.plot_trunc_cutoff_only(data)  # → FigS2_trunc.png (+ Fig_size_dist_trunc_cutoff)
    # also keep OLS power-law / LN / WB companion if needed
    # alt.plot_ols_vs_mle_calfire_mtbs(data)  # → FigS2.png (power law + LN + WB)
    _copy_alias(OUT_FIG / "FigS2_trunc.png", OUT_FIG / "FigS2.png")
    if (OUT_FIG / "FigS2_trunc.pdf").exists():
        _copy_alias(OUT_FIG / "FigS2_trunc.pdf", OUT_FIG / "FigS2.pdf")


# =============================================================================
# 05b  AIC / LR table  |  R3C4 companion (power law + three alt forms)
# =============================================================================
def run_AIC():
    print("\n=== 05b AIC + LR / Vuong (R3C4 companion) ===")
    import size_dist_AIC_LR_comparison as aic

    aic.main()


# =============================================================================
# 06  Fig. S10  |  R3C9  |  Santa Ana day fraction + size PDFs
# =============================================================================
def run_S10():
    print("\n=== 06 Fig. S10 (R3C9): Santa Ana winds ===")
    mod = _load_module(FIGS / "Fig_SoCal_SAW_powerlaw.py", "fig_saw_s10")
    mod.main()
    _copy_alias(OUT_FIG / "Fig_SoCal_SAW_powerlaw.png", OUT_FIG / "FigS10.png")
    if (OUT_FIG / "Fig_SoCal_SAW_powerlaw.pdf").exists():
        _copy_alias(OUT_FIG / "Fig_SoCal_SAW_powerlaw.pdf", OUT_FIG / "FigS10.pdf")


# =============================================================================
# 07  Fig. S11  |  R3C13  |  SG veg: NDVI / RAP NPP / ESA CCI AGB
# =============================================================================
def run_S11():
    print("\n=== 07 Fig. S11 (R3C13): SG vegetation abundance ===")
    mod = _load_module(FIGS / "Fig_SG_WUI_vs_wildland_veg_dist.py", "fig_veg_s11")
    mod.main()
    # module also writes FigS8; response labels this Fig. S11
    _copy_alias(
        OUT_FIG / "Fig_SG_WUI_vs_wildland_veg_dist.png",
        OUT_FIG / "FigS11.png",
        OUT_FIG / "FigS8.png",
    )
    if (OUT_FIG / "Fig_SG_WUI_vs_wildland_veg_dist.pdf").exists():
        _copy_alias(
            OUT_FIG / "Fig_SG_WUI_vs_wildland_veg_dist.pdf",
            OUT_FIG / "FigS11.pdf",
            OUT_FIG / "FigS8.pdf",
        )


# =============================================================================
# 08  Fig. S9  |  R3C16  |  Monthly burned area + lightning + Santa Ana days
# =============================================================================
def run_S9():
    print("\n=== 08 Fig. S9 (R3C16): monthly burned area / lightning / SAW ===")
    mod = _load_module(FIGS / "Fig_offshore_wind_monthly.py", "fig_monthly_s9")
    mod.plot_wui_wildland_monthly()
    _copy_alias(OUT_FIG / "Fig_offshore_wind_monthly_WUI.png", OUT_FIG / "FigS9.png")
    if (OUT_FIG / "Fig_offshore_wind_monthly_WUI.pdf").exists():
        _copy_alias(OUT_FIG / "Fig_offshore_wind_monthly_WUI.pdf", OUT_FIG / "FigS9.pdf")


# =============================================================================
# 09  Fig. R1c  |  R3C22  |  Western U.S. MTBS WUI vs wildland β
# =============================================================================
def run_R1c():
    print("\n=== 09 Fig. R1c (R3C22): western U.S. MTBS β ===")
    mod = _load_module(CODE / "westcoast_wui_wildland_beta.py", "fig_west_r1c")
    mod.main()  # fits + writes CSV + plot_figure
    _copy_alias(
        OUT_FIG / "Fig_westcoast_WUI_wildland_beta.png",
        OUT_FIG / "FigR1_westcoast_WUI_wildland_beta.png",
    )


# =============================================================================
# Registry
# =============================================================================
TASKS = {
    "S1": ("01 Fig. S1  | R1C10 | IgnWUI power law", run_S1),
    "R1a": ("02 Fig. R1a | R1C10 | WUI-spread ≥10%", run_R1a),
    "S6": ("03 Fig. S6  | R2C8  | Forest FWI codes", run_S6),
    "R1b": ("04 Fig. R1b | R2C8  | β vs regional FWI", run_R1b),
    "S2": ("05 Fig. S2  | R3C4  | Trunc.exp / LN / Weibull", run_S2),
    "AIC": ("05b table  | R3C4  | AIC + LR / Vuong", run_AIC),
    "S10": ("06 Fig. S10 | R3C9  | Santa Ana winds", run_S10),
    "S11": ("07 Fig. S11 | R3C13 | SG veg abundance", run_S11),
    "S9": ("08 Fig. S9  | R3C16 | Monthly BA / lightning / SAW", run_S9),
    "R1c": ("09 Fig. R1c | R3C22 | West-coast MTBS β", run_R1c),
}


def list_tasks():
    print("Response-revision figure index (run order):\n")
    for key, (desc, _) in TASKS.items():
        print(f"  {key:4s}  {desc}")
    print("\nExample: python response_revision_figures.py S1 S2 AIC")


def main(argv: list[str] | None = None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv == ["all"]:
        keys = list(TASKS.keys())
    elif argv[0] in ("list", "-h", "--help"):
        list_tasks()
        return
    else:
        keys = argv
        bad = [k for k in keys if k not in TASKS]
        if bad:
            print("Unknown keys:", bad)
            list_tasks()
            sys.exit(1)

    for k in keys:
        desc, fn = TASKS[k]
        print(f"\n{'=' * 72}\n{desc}\n{'=' * 72}")
        fn()
    print("\nDone.")


if __name__ == "__main__":
    main()

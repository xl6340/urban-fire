#!/usr/bin/env python3
"""
CalFire IgnWUI figure:
  Panel 1 — power-law PDF by Type, with beta±95% CI inset
  Panel 2 — within-fire burned-area fraction in WUI (box + points)

Power-law fits use the same log bin edges as Fig. 1:
  edges = 10.^(-1:0.05:5), geometric-mean bin centers, OLS on log10 PDF.
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

BASE = Path(__file__).resolve().parents[2]
FIRE_GPKG = BASE / "dataPrc/CalFire_IgnWUI.gpkg"
FRAC_CSV = BASE / "dataPrc/CalFire_IgnWUI_within_fire_WUI_frac.csv"
BETA_CSV = BASE / "dataPrc/CalFire_IgnWUI_powerlaw_byType_matched.csv"
OUT_PNG = BASE / "Fig/CalFire_IgnWUI_powerlaw.png"
OUT_PDF = BASE / "Fig/CalFire_IgnWUI_powerlaw.pdf"
OUT_PNG2 = BASE / "dataPrc/CalFire_IgnWUI_powerlaw.png"

TYPE_ORDER = ["WUIIgnite", "WUISpread", "Wildland"]
LABELS = {"WUIIgnite": "WUI-ignite", "WUISpread": "WUI-spread", "Wildland": "Wildland"}
COLORS = {"WUIIgnite": "#D8765A", "WUISpread": "#E0A35A", "Wildland": "#299A8F"}

# Fig. 1 MATLAB edges: 10.^(-1:0.05:5)
FIG1_EDGES = 10.0 ** np.arange(-1.0, 5.0 + 1e-12, 0.05)


def fit_powerlaw_ols_fig1(sizes, size_min=1.0, edges=FIG1_EDGES):
    """OLS power-law on log-binned PDF, matching Fig. 1 binning."""
    sizes = pd.Series(sizes).dropna()
    sizes = sizes[sizes >= size_min]
    if len(sizes) < 20:
        return None
    counts, _ = np.histogram(sizes.to_numpy(), bins=edges)
    bin_widths = np.diff(edges)
    density = counts / (len(sizes) * bin_widths)
    bin_centers = np.sqrt(edges[:-1] * edges[1:])  # geometric mean (Fig. 1)
    mask = counts > 0
    log_x = np.log10(bin_centers[mask])
    log_y = np.log10(density[mask])
    slope, intercept, r, p, se = stats.linregress(log_x, log_y)
    return {
        "beta": float(slope),
        "alpha": float(10**intercept),
        "intercept": float(intercept),
        "se": float(se),
        "r2": float(r**2),
        "p": float(p),
        "log_x": log_x,
        "log_y": log_y,
        "n": int(len(sizes)),
        "n_bins": int(mask.sum()),
        "size_min": float(size_min),
    }


def bootstrap_beta(sizes, n_iter=1000, size_min=1.0, seed=0):
    rng = np.random.default_rng(seed)
    s = pd.Series(sizes).dropna()
    s = s[s >= size_min].to_numpy()
    betas = []
    for _ in range(n_iter):
        sample = rng.choice(s, size=len(s), replace=True)
        res = fit_powerlaw_ols_fig1(sample, size_min=size_min)
        if res is not None:
            betas.append(res["beta"])
    betas = np.asarray(betas)
    return float(np.percentile(betas, 2.5)), float(np.percentile(betas, 97.5))


def main():
    fire = gpd.read_file(FIRE_GPKG)
    # Keep only fires with ≥1 Atlas/GFA ignition location inside the perimeter
    n_before = len(fire)
    fire = fire[fire["n_AtlasIgn"] > 0].copy()
    print(
        f"Restricting to fires with ignition location: {len(fire)}/{n_before} "
        f"(dropped {n_before - len(fire)} with n_AtlasIgn=0)"
    )
    print(fire["Type"].value_counts().reindex(TYPE_ORDER).to_string())

    fracs = pd.read_csv(FRAC_CSV)
    keep_cols = ["fire_fid", "Type", "wui_burn_frac", "size_km2"]
    if "n_AtlasIgn" in fracs.columns:
        fracs = fracs.loc[fracs["n_AtlasIgn"] > 0, keep_cols].copy()
    else:
        fracs = fracs[keep_cols].merge(
            fire[["fire_fid"]], on="fire_fid", how="inner"
        )
    if "Type" not in fracs.columns or fracs["Type"].isna().any():
        fracs = fire.drop(columns="geometry").merge(
            fracs[["fire_fid", "wui_burn_frac"]], on="fire_fid", how="left"
        )

    results = {}
    rows = []
    print("Fitting power laws with Fig. 1 bin edges + bootstrap CI (n=1000)...")
    for name in TYPE_ORDER:
        sizes = fire.loc[fire["Type"] == name, "size_km2"]
        res = fit_powerlaw_ols_fig1(sizes)
        if res is None:
            continue
        ci_low, ci_high = bootstrap_beta(sizes, n_iter=1000, seed=1)
        res["ci"] = [ci_low, ci_high]
        results[name] = res
        rows.append(
            {
                "Type": name,
                "n": res["n"],
                "beta": res["beta"],
                "se": res["se"],
                "ci_low": ci_low,
                "ci_high": ci_high,
                "r2": res["r2"],
                "p": res["p"],
                "n_bins": res["n_bins"],
                "size_min": res["size_min"],
                "dataset": "CalFire_matched_AtlasIgn",
                "binning": "Fig1_edges_10^(-1:0.05:5)",
                "note": "n_AtlasIgn>0 only",
            }
        )
        print(
            f"  {name}: n={res['n']}  beta={res['beta']:.4f}  "
            f"CI=[{ci_low:.4f}, {ci_high:.4f}]  bins={res['n_bins']}"
        )

    beta_tab = pd.DataFrame(rows)
    beta_tab.to_csv(BETA_CSV, index=False)
    print("Saved", BETA_CSV)

    fig = plt.figure(figsize=(11.8, 4.8), facecolor="white")
    ax1 = fig.add_axes([0.08, 0.15, 0.48, 0.78])
    ax2 = fig.add_axes([0.66, 0.15, 0.30, 0.78])

    # ── Panel 1: power-law ────────────────────────────────────────────────
    counts = {t: int((fire["Type"] == t).sum()) for t in TYPE_ORDER}
    for name, res in results.items():
        c = COLORS[name]
        ax1.scatter(
            10 ** res["log_x"], 10 ** res["log_y"],
            color=c, alpha=0.65, s=16, zorder=3, linewidths=0,
        )
        x_fit = np.linspace(res["log_x"].min(), res["log_x"].max(), 100)
        y_fit = res["beta"] * x_fit + res["intercept"]
        y_low = res["ci"][0] * x_fit + res["intercept"]
        y_high = res["ci"][1] * x_fit + res["intercept"]
        ax1.plot(
            10**x_fit, 10**y_fit, color=c, lw=2.2,
            label=f"{LABELS[name]} (n={counts[name]})",
        )
        ax1.fill_between(10**x_fit, 10**y_low, 10**y_high, color=c, alpha=0.12)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlabel(r"Fire size (km$^2$)")
    ax1.set_ylabel("Probability density")
    ax1.legend(fontsize=8.5, loc="lower left", frameon=False, handlelength=1.4)
    ax1.tick_params(direction="out")
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    ax1.text(
        -0.10, 1.02, "(A)", transform=ax1.transAxes,
        fontsize=13, fontweight="bold", va="bottom", ha="left",
    )

    # inset: beta comparison
    axins = ax1.inset_axes([0.62, 0.62, 0.34, 0.34])
    short = {"WUIIgnite": "Ignite", "WUISpread": "Spread", "Wildland": "Wild"}
    names = list(results.keys())
    xs = np.arange(len(names))
    for i, name in enumerate(names):
        b = results[name]["beta"]
        lo, hi = results[name]["ci"]
        axins.errorbar(
            i, b, yerr=[[b - lo], [hi - b]],
            fmt="o", color=COLORS[name], markersize=6.5, capsize=3.5, lw=1.4,
            markeredgecolor="white", markeredgewidth=0.6,
        )
    axins.set_xticks(xs)
    axins.set_xticklabels([short[n] for n in names], fontsize=7)
    axins.set_ylabel(r"$\beta$ value", fontsize=8)
    axins.tick_params(labelsize=7, length=3)
    axins.set_xlim(-0.55, len(names) - 0.45)
    axins.set_facecolor("0.97")
    for sp in axins.spines.values():
        sp.set_linewidth(0.7)
        sp.set_color("0.55")

    # ── Panel 2: box plot of within-fire WUI fraction ──────────────────────
    data_box = [
        fracs.loc[fracs["Type"] == t, "wui_burn_frac"].dropna().to_numpy()
        for t in TYPE_ORDER
    ]
    bp = ax2.boxplot(
        data_box,
        positions=np.arange(1, len(TYPE_ORDER) + 1),
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="0.12", lw=1.6),
        whiskerprops=dict(color="0.35", lw=1.0),
        capprops=dict(color="0.35", lw=1.0),
        boxprops=dict(lw=1.0),
    )
    for patch, name in zip(bp["boxes"], TYPE_ORDER):
        patch.set_facecolor(COLORS[name])
        patch.set_alpha(0.45)
        patch.set_edgecolor(COLORS[name])

    ax2.set_xticks(np.arange(1, len(TYPE_ORDER) + 1))
    ax2.set_xticklabels([LABELS[t] for t in TYPE_ORDER], fontsize=9)
    ax2.set_ylabel("Within-WUI burned area / total burned area")
    ax2.set_ylim(-0.02, 1.05)
    ax2.set_xlim(0.4, len(TYPE_ORDER) + 0.6)
    ax2.set_yticks(np.arange(0, 1.01, 0.25))
    ax2.yaxis.set_major_formatter(lambda v, _pos: f"{v:.0%}")
    ax2.tick_params(direction="out")
    for sp in ("top", "right"):
        ax2.spines[sp].set_visible(False)
    ax2.text(
        -0.14, 1.02, "(B)", transform=ax2.transAxes,
        fontsize=13, fontweight="bold", va="bottom", ha="left",
    )

    for out in (OUT_PNG, OUT_PDF, OUT_PNG2):
        fig.savefig(out, dpi=200 if out.suffix == ".png" else None, facecolor="white")
        print("Saved", out)

    print("\nWithin-fire WUI fraction summary:")
    for t in TYPE_ORDER:
        s = fracs.loc[fracs["Type"] == t, "wui_burn_frac"]
        print(
            f"  {t}: n={len(s)}  min={s.min():.1%}  p25={s.quantile(0.25):.1%}  "
            f"med={s.median():.1%}  p75={s.quantile(0.75):.1%}  max={s.max():.1%}"
        )


if __name__ == "__main__":
    main()

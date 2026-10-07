#!/usr/bin/env python3
"""
SoCal CalFire Santa Ana figure (Fig. S10):

  (A) Fraction of fires with ≥1 Santa Ana day during the burn window
      [IDate, FDate] — WUI vs Wildland
  (B) Mean fire size for four groups:
      WUI × Santa Ana / non–Santa Ana, Wildland × Santa Ana / non–Santa Ana

Uses any-SAW-day-during-fire flags (R1D-SAWRI era only, 1990–2018):
  dataPrc/CalFire_SoCal_SAW_extended_flags.csv
"""


from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

BASE = Path(__file__).resolve().parents[2]
FLAGS = BASE / "dataPrc/CalFire_SoCal_SAW_extended_flags.csv"
OUT_PNG = BASE / "Fig/Fig_SoCal_SAW_powerlaw.png"
OUT_PDF = BASE / "Fig/Fig_SoCal_SAW_powerlaw.pdf"
OUT_CSV = BASE / "dataPrc/CalFire_SoCal_SAW_mean_size.csv"
OUT_S10_PNG = BASE / "Fig/FigS10.png"
OUT_S10_PDF = BASE / "Fig/FigS10.pdf"

C_WUI = np.array([216, 118, 89]) / 255
C_WILD = np.array([41, 157, 143]) / 255
C_WUI_SAW = np.array([140, 70, 45]) / 255
C_WUI_NON = np.array([230, 170, 140]) / 255
C_WILD_SAW = np.array([25, 110, 100]) / 255
C_WILD_NON = np.array([140, 195, 185]) / 255


# Typography (pt): keep hierarchy tight so panel labels don't dominate
FS_PANEL = 9
FS_AXIS = 8.5
FS_TICK = 8
FS_ANNOT = 7
FS_LEGEND = 7.5


def group_stats(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ft in ["WUI", "Wildland"]:
        for saw, lab in [(True, "Santa Ana"), (False, "non–Santa Ana")]:
            s = pd.to_numeric(
                df.loc[(df["FireType"] == ft) & (df["SAW"] == saw), "size"],
                errors="coerce",
            )
            s = s[np.isfinite(s) & (s > 0)]
            rows.append(
                {
                    "FireType": ft,
                    "Santa_Ana": saw,
                    "label": lab,
                    "n": int(len(s)),
                    "mean_km2": float(s.mean()) if len(s) else np.nan,
                    "median_km2": float(s.median()) if len(s) else np.nan,
                    "se_km2": float(s.sem()) if len(s) else np.nan,
                    "std_km2": float(s.std(ddof=1)) if len(s) > 1 else np.nan,
                }
            )
    return pd.DataFrame(rows)


def draw_fraction_bar(ax, df: pd.DataFrame):
    """(A) % of fires with ≥1 Santa Ana day during burn duration."""
    rows = []
    for ft, color in [("WUI", C_WUI), ("Wildland", C_WILD)]:
        g = df[df["FireType"] == ft]
        n = len(g)
        n_saw = int(g["SAW"].sum())
        rows.append({"FireType": ft, "n": n, "n_SAW": n_saw, "frac": n_saw / n, "color": color})
    tab = pd.DataFrame(rows)
    vals = tab["frac"].to_numpy() * 100
    bars = ax.bar(
        np.arange(2), vals, width=0.58,
        color=tab["color"].tolist(), edgecolor="none", lw=0, zorder=2,
    )
    y_top = max(22, vals.max() * 1.35)
    ax.set_ylim(0, y_top)
    pad = y_top * 0.025
    for b, v, n_s, n in zip(bars, vals, tab["n_SAW"], tab["n"]):
        ax.text(
            b.get_x() + b.get_width() / 2, b.get_height() + pad,
            f"{v:.1f}%\n({int(n_s)}/{int(n)})",
            ha="center", va="bottom", fontsize=FS_ANNOT, color="0.2",
        )
    ax.set_xticks(np.arange(2))
    ax.set_xticklabels(["WUI", "Wildland"], fontsize=FS_TICK)
    ax.set_ylabel("Santa Ana-affected fires (%)", fontsize=FS_AXIS)
    ax.tick_params(direction="out", labelsize=FS_TICK)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return tab


def _fmt_p(p: float) -> str:
    if p < 0.001:
        return r"$p$ < 0.001"
    return rf"$p$ = {p:.2f}"


def draw_mean_size_bars(ax, stats: pd.DataFrame, df: pd.DataFrame):
    """(B) Mean fire size ± SE for four groups (grouped by FireType)."""
    # order within each FireType: Santa Ana, non–Santa Ana
    colors = {
        ("WUI", True): C_WUI_SAW,
        ("WUI", False): C_WUI_NON,
        ("Wildland", True): C_WILD_SAW,
        ("Wildland", False): C_WILD_NON,
    }
    x_centers = np.array([0.0, 1.15])
    width = 0.38
    offsets = [-width / 2, width / 2]

    # top of mean±SE + text block for each bar (for p bracket placement)
    label_tops = {}

    for i, ft in enumerate(["WUI", "Wildland"]):
        for j, saw in enumerate([True, False]):
            r = stats[(stats["FireType"] == ft) & (stats["Santa_Ana"] == saw)].iloc[0]
            x = x_centers[i] + offsets[j]
            ax.bar(
                x, r["mean_km2"], width=width * 0.95,
                color=colors[(ft, saw)], edgecolor="none", lw=0, zorder=2,
            )
            ax.errorbar(
                x, r["mean_km2"], yerr=r["se_km2"],
                fmt="none", ecolor="0.25", elinewidth=0.9, capsize=2.5, zorder=3,
            )
            y_lab = r["mean_km2"] + r["se_km2"] + 2.5
            ax.text(
                x, y_lab,
                f'{r["mean_km2"]:.1f}\n($n$={r["n"]})',
                ha="center", va="bottom", fontsize=FS_ANNOT, color="0.25",
            )
            # approximate text block height (~2 lines of FS_ANNOT)
            label_tops[(ft, saw)] = y_lab + 14.0

    # Mann–Whitney p above mean/n labels within each FireType
    for i, ft in enumerate(["WUI", "Wildland"]):
        a = pd.to_numeric(
            df.loc[(df["FireType"] == ft) & (df["SAW"]), "size"], errors="coerce"
        ).to_numpy()
        b = pd.to_numeric(
            df.loc[(df["FireType"] == ft) & (~df["SAW"]), "size"], errors="coerce"
        ).to_numpy()
        a = a[np.isfinite(a) & (a > 0)]
        b = b[np.isfinite(b) & (b > 0)]
        _, p = mannwhitneyu(a, b, alternative="two-sided")

        x0 = x_centers[i] + offsets[0]
        x1 = x_centers[i] + offsets[1]
        y = max(label_tops[(ft, True)], label_tops[(ft, False)])
        h = 3.5
        ax.plot([x0, x1], [y + h, y + h], color="0.3", lw=0.8, clip_on=False)
        ax.text(
            (x0 + x1) / 2, y + h + 1.0, _fmt_p(p),
            ha="center", va="bottom", fontsize=FS_ANNOT, color="0.25",
        )
        label_tops[(ft, "p")] = y + h + 10.0

    ax.set_xticks(x_centers)
    ax.set_xticklabels(["WUI", "Wildland"], fontsize=FS_TICK)
    ax.set_ylabel(r"Fire size (km$^2$)", fontsize=FS_AXIS)
    ymax = max(label_tops.values())
    ax.set_ylim(0, ymax * 1.08)
    ax.tick_params(direction="out", labelsize=FS_TICK)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # legend: Santa Ana vs non (hue from WUI pair; FireType shown on x-axis)
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor=C_WUI_SAW, edgecolor="none", label="Santa Ana"),
        Patch(facecolor=C_WUI_NON, edgecolor="none", label="non–Santa Ana"),
    ]
    ax.legend(
        handles=handles, frameon=False, fontsize=FS_LEGEND,
        loc="upper right", handlelength=1.2,
    )


def main():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": FS_TICK,
        "axes.labelsize": FS_AXIS,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND,
    })
    df = pd.read_csv(FLAGS, parse_dates=["IDate"])
    # Restrict to official R1D-SAWRI years (no post-2018 gridMET proxy)
    df = df[df["in_analysis"] & df["IDate"].notna() & (df["IDate"] <= "2018-12-31")].copy()
    df["size"] = pd.to_numeric(df["size"], errors="coerce")

    stats = group_stats(df)
    stats.to_csv(OUT_CSV, index=False)
    print(stats.to_string(index=False))
    print("Saved", OUT_CSV)

    fig, axes = plt.subplots(
        1, 2, figsize=(6.4, 3.1), facecolor="w",
        gridspec_kw={"width_ratios": [0.9, 1.35]},
    )
    fig.subplots_adjust(wspace=0.38, left=0.10, right=0.98, top=0.86, bottom=0.16)
    ax_a, ax_b = axes

    frac_tab = draw_fraction_bar(ax_a, df)
    draw_mean_size_bars(ax_b, stats, df)

    for ax, lab, xoff in ((ax_a, "(A)", -0.18), (ax_b, "(B)", -0.12)):
        ax.text(
            xoff, 1.03, lab, transform=ax.transAxes,
            fontsize=FS_PANEL, fontweight="normal", va="bottom", ha="left",
        )

    for pth in (OUT_PNG, OUT_PDF, OUT_S10_PNG, OUT_S10_PDF):
        fig.savefig(pth, dpi=300, bbox_inches="tight", facecolor="w")
    print("Saved", OUT_PNG)

    # brief print for manuscript
    for _, r in frac_tab.iterrows():
        print(f"{r['FireType']}: {100*r['frac']:.1f}% ({r['n_SAW']}/{r['n']})")
    for _, r in stats.iterrows():
        print(
            f"{r['FireType']:9s} {r['label']:14s}  "
            f"mean={r['mean_km2']:.2f} ± {r['se_km2']:.2f} km²  n={r['n']}"
        )
    plt.close(fig)


if __name__ == "__main__":
    main()

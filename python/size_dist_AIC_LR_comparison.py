#!/usr/bin/env python3
"""
Model comparison table for Reviewer 3 (Cumming / Reed–McKelvey):

Four forms matching FigS2_trunc + power law:
  1) Power law            p(x) ∝ x^β
  2) Trunc. exp. (cutoff) p(x) ∝ x^β exp(−x/xc)   [labeled Trunc. exp. in figure]
  3) Lognormal (quad)     ln f = ln a − β ln A − ψ [ln A]²
  4) Weibull (Reed–McKelvey Model I)

Reports:
  - R² from log-binned PDF fits (OLS for PL/LN/WB; cutoff curve vs bins for Trunc.exp.)
  - MLE log-likelihood, AIC, ΔAIC
  - Nested LR test: power law ⊂ trunc.exp. cutoff (χ², df=1)
  - Vuong LR for non-nested pairs vs power law

Bootstrap CIs / OLS SEs in the main text come from OLS on log-binned PDFs
(not MLE); AIC/LR here use MLE on the same xmin-truncated samples.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))

from fit_size_distributions_alt import (  # noqa: E402
    OUT_PRC,
    fit_lognormal_mle_orig,
    fit_lognormal_quad_ols,
    fit_pl_cutoff_mle,
    fit_powerlaw_mle_orig,
    fit_powerlaw_ols,
    fit_weibull_mle,
    fit_weibull_ols,
    load_groups,
    log_binned_pdf,
    pdf_pl_cutoff,
    pdf_weibull,
)

# Same completeness thresholds as FigS2_trunc / FigS2
XMIN_FIG = {"CalFire": 1.0, "MTBS": 4.0, "Atlas": 1.0, "FIRED": 1.0}
DATASETS = ["CalFire", "MTBS", "Atlas", "FIRED"]
TYPES = ["Urban-edge", "Wildland"]
MODELS = ["powerlaw", "truncexpon", "lognormal", "weibull"]
K_PARAMS = {"powerlaw": 1, "truncexpon": 2, "lognormal": 2, "weibull": 2}


def r2_log_pdf(x, y, yhat):
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(yhat) & (y > 0) & (yhat > 0)
    if int(m.sum()) < 3:
        return np.nan
    yt, yh = np.log10(y[m]), np.log10(yhat[m])
    ss_tot = np.sum((yt - yt.mean()) ** 2)
    if ss_tot <= 0:
        return np.nan
    return float(1.0 - np.sum((yt - yh) ** 2) / ss_tot)


def pointwise_ll_powerlaw(s, beta, xmin):
    # f(A)=α A^{-β}, α=(β-1) xmin^{β-1}; β is Clauset-style (>1)
    return np.log(beta - 1.0) - np.log(xmin) - beta * np.log(s / xmin)


def pointwise_ll_cutoff(s, beta, xc, xmin):
    from fit_size_distributions_alt import _pl_cutoff_logZ

    logz = _pl_cutoff_logZ(beta, xc, xmin)
    return beta * np.log(s) - s / xc - logz


def pointwise_ll_lognormal(s, beta, psi, xmin):
    from fit_size_distributions_alt import _trapz_logspace

    def log_u(A):
        lnA = np.log(A)
        return -beta * lnA - psi * lnA**2

    Z = _trapz_logspace(log_u, xmin)
    return log_u(s) - np.log(Z)


def pointwise_ll_weibull(s, alpha, beta, xmin):
    from fit_size_distributions_alt import ln_pdf_weibull_rm

    return ln_pdf_weibull_rm(s, alpha, beta, xmin)


def vuong_test(ll1, ll2):
    """Vuong (1989) non-nested LR; positive favors model 1."""
    m = np.asarray(ll1, float) - np.asarray(ll2, float)
    m = m[np.isfinite(m)]
    n = len(m)
    if n < 10:
        return np.nan, np.nan, np.nan
    lr = float(np.sum(m))
    omega = float(np.std(m, ddof=1))
    if omega <= 0:
        return lr, np.nan, np.nan
    v = lr / (np.sqrt(n) * omega)
    p = float(2.0 * (1.0 - stats.norm.cdf(abs(v))))
    return lr, float(v), p


def fit_group(sizes, xmin):
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s >= xmin)]
    n = len(s)
    if n < 20:
        return None

    x, y = log_binned_pdf(s)

    # --- OLS / display R² (figure convention) ---
    pl_ols = fit_powerlaw_ols(s, size_min=xmin)
    ln_ols = fit_lognormal_quad_ols(s, size_min=xmin)
    wb_ols = fit_weibull_ols(s, size_min=xmin)

    # --- MLE for AIC / LR ---
    pl_mle = fit_powerlaw_mle_orig(s, size_min=xmin)
    cut = fit_pl_cutoff_mle(s, size_min=xmin)
    ln_mle = fit_lognormal_mle_orig(s, size_min=xmin)
    wb_mle = fit_weibull_mle(s, size_min=xmin)

    out = {
        "n": n,
        "xmin": xmin,
        "powerlaw": {},
        "truncexpon": {},
        "lognormal": {},
        "weibull": {},
    }

    if pl_ols and pl_mle:
        out["powerlaw"] = {
            "beta_ols": pl_ols["beta"],
            "beta_ols_se": pl_ols["beta_se"],
            "r2": pl_ols["r2"],
            "beta_mle": -pl_mle["beta"],  # report as negative slope like figure
            "alpha_mle": pl_mle["beta"],  # Clauset α≡β_orig
            "ll": pl_mle["ll"],
            "aic": pl_mle["aic"],
            "k": 1,
            "ll_i": pointwise_ll_powerlaw(s, pl_mle["beta"], xmin),
        }

    if cut:
        yhat = pdf_pl_cutoff(x, cut["beta"], cut["xc"], xmin)
        out["truncexpon"] = {
            "beta_ols": cut["beta"],  # figure reports MLE β with bin R²
            "beta_ols_se": np.nan,
            "r2": r2_log_pdf(x, y, yhat),
            "beta_mle": cut["beta"],
            "xc": cut["xc"],
            "ll": cut["ll"],
            "aic": cut["aic"],
            "k": 2,
            "ll_i": pointwise_ll_cutoff(s, cut["beta"], cut["xc"], xmin),
        }

    if ln_ols and ln_mle:
        out["lognormal"] = {
            "beta_ols": -ln_ols["beta"],  # figure sign convention
            "beta_ols_se": ln_ols.get("beta_se", np.nan),
            "r2": ln_ols["r2"],
            "beta_mle": ln_mle["beta"],
            "psi": ln_mle["psi"],
            "ll": ln_mle["ll"],
            "aic": ln_mle["aic"],
            "k": 2,
            "ll_i": pointwise_ll_lognormal(s, ln_mle["beta"], ln_mle["psi"], xmin),
        }

    if wb_ols and wb_mle:
        # OLS R² from figure form; optional bin R² of MLE curve
        yhat_wb = pdf_weibull(x, wb_mle["alpha"], wb_mle["beta"], xmin)
        out["weibull"] = {
            "beta_ols": -wb_ols["beta"],
            "beta_ols_se": wb_ols.get("beta_se", np.nan),
            "r2": wb_ols["r2"],
            "r2_mle_curve": r2_log_pdf(x, y, yhat_wb),
            "beta_mle": -wb_mle["beta"],
            "alpha_mle": wb_mle["alpha"],
            "ll": wb_mle["ll"],
            "aic": wb_mle["aic"],
            "k": 2,
            "ll_i": pointwise_ll_weibull(s, wb_mle["alpha"], wb_mle["beta"], xmin),
        }

    return out


def main():
    data = load_groups()
    rows = []
    lr_rows = []

    for ds in DATASETS:
        xmin = XMIN_FIG[ds]
        for ft in TYPES:
            label = "WUI" if ft == "Urban-edge" else "Wildland"
            print(f"Fitting {ds} {label} (xmin={xmin})...")
            res = fit_group(data[ds][ft], xmin)
            if res is None:
                print("  skipped (n<20)")
                continue

            aics = {m: res[m]["aic"] for m in MODELS if res[m]}
            best = min(aics, key=aics.get)
            aic_min = aics[best]

            for m in MODELS:
                if not res[m]:
                    continue
                r = res[m]
                rows.append(
                    {
                        "dataset": ds,
                        "FireType": label,
                        "model": m,
                        "n": res["n"],
                        "xmin": xmin,
                        "beta_display": r.get("beta_ols", np.nan),
                        "R2": r.get("r2", np.nan),
                        "ll_MLE": r["ll"],
                        "AIC": r["aic"],
                        "delta_AIC": r["aic"] - aic_min,
                        "k_params": r["k"],
                        "xc": r.get("xc", np.nan),
                        "best_AIC": m == best,
                    }
                )

            # Nested LR: power law ⊂ trunc.exp. cutoff
            if res["powerlaw"] and res["truncexpon"]:
                ll0, ll1 = res["powerlaw"]["ll"], res["truncexpon"]["ll"]
                lr_stat = 2.0 * (ll1 - ll0)
                # boundary null (1/xc→0); report χ²(1) as conventional upper-bound p
                p_chi2 = float(stats.chi2.sf(max(lr_stat, 0.0), df=1))
                lr_rows.append(
                    {
                        "dataset": ds,
                        "FireType": label,
                        "comparison": "truncexpon_vs_powerlaw_nested",
                        "LR": lr_stat,
                        "df": 1,
                        "p_chi2": p_chi2,
                        "Vuong_V": np.nan,
                        "Vuong_p": np.nan,
                        "favor": "truncexpon" if lr_stat > 0 else "powerlaw",
                    }
                )

            # Vuong vs power law for non-nested / all pairs
            if res["powerlaw"]:
                for m in ["truncexpon", "lognormal", "weibull"]:
                    if not res[m]:
                        continue
                    lr_sum, v, p = vuong_test(res[m]["ll_i"], res["powerlaw"]["ll_i"])
                    lr_rows.append(
                        {
                            "dataset": ds,
                            "FireType": label,
                            "comparison": f"{m}_vs_powerlaw_Vuong",
                            "LR": lr_sum,
                            "df": np.nan,
                            "p_chi2": np.nan,
                            "Vuong_V": v,
                            "Vuong_p": p,
                            "favor": m
                            if (np.isfinite(v) and v > 0)
                            else ("powerlaw" if np.isfinite(v) else "inconclusive"),
                        }
                    )

    fit_df = pd.DataFrame(rows)
    lr_df = pd.DataFrame(lr_rows)

    out_fit = OUT_PRC / "size_dist_model_comparison_AIC_LR.csv"
    out_lr = OUT_PRC / "size_dist_model_LR_tests.csv"
    fit_df.to_csv(out_fit, index=False)
    lr_df.to_csv(out_lr, index=False)

    # Wide summary: R² | AIC | ΔAIC by model
    print("\n" + "=" * 88)
    print("R² (log-binned) | AIC (MLE) | ΔAIC  — four forms, same xmin as figure")
    print("=" * 88)
    for ds in DATASETS:
        for ft in ["WUI", "Wildland"]:
            sub = fit_df[(fit_df["dataset"] == ds) & (fit_df["FireType"] == ft)]
            if sub.empty:
                continue
            n = int(sub["n"].iloc[0])
            print(f"\n{ds}  {ft}  (n={n}, xmin={XMIN_FIG[ds]} km²)")
            print(f"  {'model':12s}  {'β_disp':>8s}  {'R²':>6s}  {'LL':>10s}  {'AIC':>10s}  {'ΔAIC':>7s}  best")
            for _, row in sub.iterrows():
                star = "*" if row["best_AIC"] else ""
                print(
                    f"  {row['model']:12s}  {row['beta_display']:8.3f}  {row['R2']:6.3f}  "
                    f"{row['ll_MLE']:10.1f}  {row['AIC']:10.1f}  {row['delta_AIC']:7.1f}  {star}"
                )

    print("\n" + "=" * 88)
    print("Likelihood-ratio / Vuong tests (vs power law)")
    print("=" * 88)
    print(lr_df.to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    print(f"\nSaved:\n  {out_fit}\n  {out_lr}")
    print(
        "\nNote for reviewer: figure β±SE / bootstrap CIs are from OLS on log-binned PDFs; "
        "AIC and LR use MLE on individual sizes (same xmin). "
        "'Trunc. exp.' in the figure is the power-law × exponential cutoff form."
    )


if __name__ == "__main__":
    main()

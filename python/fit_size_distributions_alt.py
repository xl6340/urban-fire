#!/usr/bin/env python3
"""
Compare fire-size distribution fits for main WUI vs wildland figures.

Models (x >= xmin), following Cumming (2001) / Reed & McKelvey (2002):
  1) Power law             — OLS on log-binned PDF (paper main method) + MLE
  2) Lognormal (quad)      — ln f(A) = ln a − β ln A − ψ [ln A]²
  3) Weibull (Reed–McKelvey) — ln f = ln α − β ln A − α/(1−β)(A^{1−β}−A₀^{1−β}), A₀=1 km²
  4) Truncated power law  — right-truncated Pareto on [xmin, xmax]

Outputs figures + tables under Fig/ and dataPrc/.
"""

from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import brentq, curve_fit, minimize, minimize_scalar

BASE = Path(__file__).resolve().parents[1]
OUT_FIG = BASE / "Fig"
OUT_PRC = BASE / "dataPrc"
OUT_FIG.mkdir(exist_ok=True)
OUT_PRC.mkdir(exist_ok=True)

COLOR = {
    "Urban-edge": np.array([216, 118, 89]) / 255,
    "WUI": np.array([216, 118, 89]) / 255,
    "Wildland": np.array([41, 157, 143]) / 255,
}
LS = {
    "powerlaw_ols": "-",
    "lognormal_ols": (0, (2.5, 1.5)),  # short dashes
    "truncpower_ols": "-.",
    "pl_cutoff": "-",
    "truncexpon_ols": "--",
    "weibull_ols": ":",
}


# Fig. 1 MATLAB edges: 10.^(-1:0.05:5); geometric-mean centers
FIG1_EDGES = 10.0 ** np.arange(-1.0, 5.0 + 1e-12, 0.05)


# ── empirical PDF (same bin edges as Fig. 1) ─────────────────────────────────
def log_binned_pdf(sizes: np.ndarray, n_bins: int | None = None, edges=None):
    """
    Log-binned PDF. Default uses Fig. 1 fixed edges ``10.^(-1:0.05:5)``.
    Pass ``n_bins`` for legacy adaptive bins from sample min–max.
    """
    sizes = np.asarray(sizes, float)
    sizes = sizes[np.isfinite(sizes) & (sizes > 0)]
    if len(sizes) < 10:
        return np.array([]), np.array([])
    if edges is None:
        if n_bins is not None:
            lo, hi = sizes.min(), sizes.max()
            if lo >= hi:
                return np.array([]), np.array([])
            edges = np.logspace(np.log10(lo), np.log10(hi), n_bins + 1)
        else:
            edges = FIG1_EDGES
    counts, _ = np.histogram(sizes, bins=edges)
    widths = np.diff(edges)
    centers = np.sqrt(edges[:-1] * edges[1:])
    dens = counts / (len(sizes) * widths)
    m = dens > 0
    return centers[m], dens[m]


# ── Power law: OLS on log-binned PDF (main-text β) ───────────────────────────
def fit_powerlaw_ols(sizes, n_bins=None, size_min=None):
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s)]
    if size_min is not None:
        s = s[s >= size_min]
    if len(s) < 20:
        return None
    x, y = log_binned_pdf(s, n_bins=n_bins)
    if len(x) < 5:
        return None
    log_x = np.log10(x)
    log_y = np.log10(y)
    slope, intercept, r, p, se = stats.linregress(log_x, log_y)
    # PDF ∝ x^β  (β negative); MLE-style α = -β
    # log10(pdf) = beta * log10(x) + intercept  =>  pdf = 10^intercept * x^beta
    return {
        "model": "powerlaw_ols",
        "beta": slope,
        "beta_se": se,
        "intercept": intercept,  # log10 scale
        "alpha": -slope,
        "r2": r**2,
        "n": len(s),
        "xmin": float(s.min()),
        "xmax": float(s.max()),
    }


def pdf_powerlaw_ols(x, beta, intercept):
    """OLS log-binned power-law PDF: pdf = 10^intercept * x^beta."""
    return (10.0 ** intercept) * np.power(x, beta)


def fit_powerlaw_mle(sizes, size_min=None):
    """Clauset-style continuous power-law MLE for fixed xmin."""
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s)]
    if size_min is not None:
        s = s[s >= size_min]
    else:
        size_min = float(s.min())
        s = s[s >= size_min]
    n = len(s)
    if n < 20:
        return None
    xmin = float(size_min)
    alpha = 1.0 + n / np.sum(np.log(s / xmin))
    # log-likelihood
    ll = n * np.log(alpha - 1) - n * np.log(xmin) - alpha * np.sum(np.log(s / xmin))
    aic = 2 * 1 - 2 * ll
    bic = np.log(n) * 1 - 2 * ll
    beta = -(alpha)  # log-log PDF slope ≈ -α
    # analytical SE for alpha
    se_alpha = (alpha - 1) / np.sqrt(n)
    return {
        "model": "powerlaw_mle",
        "alpha": alpha,
        "alpha_se": se_alpha,
        "beta": beta,
        "beta_se": se_alpha,
        "ll": ll,
        "aic": aic,
        "bic": bic,
        "n": n,
        "xmin": xmin,
        "k_params": 1,
    }


def pdf_powerlaw(x, alpha, xmin):
    return (alpha - 1) / xmin * (x / xmin) ** (-alpha)


def sf_powerlaw(x, alpha, xmin):
    return (x / xmin) ** (1 - alpha)


# ── Truncated power law (right-truncated Pareto) ─────────────────────────────
# p(x) = C x^{-α},  xmin ≤ x ≤ xmax
# C = (α-1) / (xmin^{1-α} - xmax^{1-α})   (α ≠ 1)
# On log–log: straight like power law, then hard cut at xmax.


def _tpl_norm_const(alpha, xmin, xmax):
    if abs(alpha - 1.0) < 1e-12:
        return 1.0 / np.log(xmax / xmin)
    return (alpha - 1.0) / (xmin ** (1.0 - alpha) - xmax ** (1.0 - alpha))


def pdf_truncpower(x, alpha, xmin, xmax):
    x = np.asarray(x, float)
    out = np.zeros_like(x, dtype=float)
    m = (x >= xmin) & (x <= xmax)
    C = _tpl_norm_const(alpha, xmin, xmax)
    out[m] = C * np.power(x[m], -alpha)
    return out


def sf_truncpower(x, alpha, xmin, xmax):
    """P(X > x) for truncated power law on [xmin, xmax]."""
    x = np.asarray(x, float)
    if abs(alpha - 1.0) < 1e-12:
        num = np.log(xmax / np.minimum(np.maximum(x, xmin), xmax))
        den = np.log(xmax / xmin)
        sf = num / den
    else:
        a = 1.0 - alpha
        num = np.power(np.minimum(np.maximum(x, xmin), xmax), a) - np.power(xmax, a)
        den = np.power(xmin, a) - np.power(xmax, a)
        sf = num / den
    return np.where(x < xmin, 1.0, np.where(x >= xmax, 0.0, sf))


def fit_truncpower_ols(sizes, n_bins=None, size_min=None):
    """
    Truncated power law via OLS on log-binned PDF (same as power-law OLS),
    with xmax = max(size). Returns normalized truncated-PDF parameters.
    """
    ols = fit_powerlaw_ols(sizes, n_bins=n_bins, size_min=size_min)
    if ols is None:
        return None
    xmin = float(size_min) if size_min is not None else ols["xmin"]
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s >= xmin)]
    xmax = float(np.max(s))
    alpha = float(ols["alpha"])
    if alpha <= 1.0:
        # still allow α<=1 for truncated support
        pass
    return {
        "model": "truncpower_ols",
        "alpha": alpha,
        "beta": -alpha,
        "beta_se": ols["beta_se"],
        "intercept": ols["intercept"],  # free OLS intercept (unnormalized)
        "r2": ols["r2"],
        "n": ols["n"],
        "xmin": xmin,
        "xmax": xmax,
    }


def fit_truncpower_mle(sizes, size_min=None):
    """MLE α for truncated power law on [xmin, xmax], xmax=max(size)."""
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s)]
    if size_min is not None:
        s = s[s >= size_min]
    else:
        size_min = float(s.min())
        s = s[s >= size_min]
    n = len(s)
    if n < 20:
        return None
    xmin = float(size_min)
    xmax = float(np.max(s))
    if xmax <= xmin:
        return None

    # Untruncated Clauset start, then maximize truncated log-likelihood
    alpha0 = 1.0 + n / np.sum(np.log(s / xmin))

    def neg_ll(alpha):
        alpha = float(alpha)
        if alpha <= 0:
            return 1e300
        C = _tpl_norm_const(alpha, xmin, xmax)
        if not np.isfinite(C) or C <= 0:
            return 1e300
        return -(n * np.log(C) - alpha * np.sum(np.log(s)))

    # allow α>1 typically; widen bounds slightly for truncated case
    res = minimize_scalar(neg_ll, bounds=(0.5, 5.0), method="bounded")
    alpha = float(res.x) if res.success else float(alpha0)
    ll = -neg_ll(alpha)
    aic = 2 * 1 - 2 * ll
    bic = np.log(n) * 1 - 2 * ll
    return {
        "model": "truncpower_mle",
        "alpha": alpha,
        "beta": -alpha,
        "ll": ll,
        "aic": aic,
        "bic": bic,
        "n": n,
        "xmin": xmin,
        "xmax": xmax,
        "k_params": 1,
    }


def pdf_truncpower_xc(x, beta, xc1, xc2):
    """
    User formula (β > 0):
      p(x) = (1-β) / (x_{c2}^{1-β} - x_{c1}^{1-β}) * x^{-β},  x_{c1} ≤ x ≤ x_{c2}
    Identical to pdf_truncpower(..., alpha=beta, xmin=xc1, xmax=xc2).
    """
    x = np.asarray(x, float)
    den = np.power(xc2, 1.0 - beta) - np.power(xc1, 1.0 - beta)
    C = (1.0 - beta) / den
    out = np.zeros_like(x, dtype=float)
    m = (x >= xc1) & (x <= xc2) & np.isfinite(C) & (C > 0)
    out[m] = C * np.power(x[m], -beta)
    return out


def fit_truncpower_xc(sizes, xc1):
    """
    Estimate β and upper cutoff x_{c2} in pdf_truncpower_xc, with x_{c1} fixed.

    x_{c2} cannot be < sample max (else some fires have p=0). Estimate it as a
    free parameter by Cumming / Hannon–Dahiya, then MLE β from this exact PDF.
    """
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s >= xc1)]
    if len(s) < 20:
        return None
    xc1 = float(xc1)
    xmax_obs = float(np.max(s))
    cum = fit_cumming_trunc_pareto(s, size_min=xc1)
    if cum is None:
        xc2 = xmax_obs
        beta0 = 1.0 + len(s) / np.sum(np.log(s / xc1))
    else:
        xc2 = float(max(cum["xmax"], xmax_obs))
        beta0 = float(cum["alpha"])

    n = len(s)
    logsum = float(np.sum(np.log(s)))

    def nll(beta):
        beta = float(beta)
        den = np.power(xc2, 1.0 - beta) - np.power(xc1, 1.0 - beta)
        C = (1.0 - beta) / den
        if not np.isfinite(C) or C <= 0:
            return 1e300
        return -(n * np.log(C) - beta * logsum)

    res = minimize_scalar(nll, bounds=(0.2, 4.5), method="bounded")
    beta = float(res.x) if res.success else float(beta0)
    return {
        "model": "truncpower_xc",
        "beta_pow": beta,          # β in x^{-β} (user formula, positive)
        "beta": -beta,             # plotted slope, negative
        "alpha": beta,
        "xc1": xc1,
        "xc2": xc2,
        "xmax_obs": xmax_obs,
        "ll": float(-nll(beta)),
        "n": n,
        "k_params": 2,
    }


def fit_cumming_trunc_pareto(sizes, size_min=None):
    """
    Cumming (2001) truncated exponential on log-size, eqs. [5]–[8].

    x = ln(z / t), t = xmin. Fit scale σ and upper bound b by the
    Hannon & Dahiya (1999) one-step update of b beyond the sample max,
    then transform back to a truncated Pareto on size:
        p(z) ∝ z^{-(1/σ + 1)},   t ≤ z ≤ t e^{b}.
    """
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s > 0)]
    if size_min is not None:
        s = s[s >= size_min]
    else:
        size_min = float(s.min())
        s = s[s >= size_min]
    n = len(s)
    if n < 20:
        return None
    t = float(size_min)
    x = np.log(s / t)  # natural log, as in Cumming
    xbar = float(np.mean(x))
    b = float(np.max(x))
    zmax_obs = float(np.max(s))
    if b <= 0 or xbar <= 0:
        return None

    def mean_trunc_exp(sigma, b_):
        # E[X] = σ − b / (e^{b/σ} − 1)  for Exp(scale=σ) truncated to [0, b]
        r = b_ / sigma
        if r > 80:
            return sigma
        return sigma - b_ / np.expm1(r)

    def solve_sigma(b_):
        # Solution exists iff 0 < x̄ < b/2
        if xbar >= 0.5 * b_ - 1e-12:
            return None

        def f(sig):
            return mean_trunc_exp(sig, b_) - xbar

        lo, hi = 1e-8, max(xbar * 50.0, 1.0)
        f_lo, f_hi = f(lo), f(hi)
        if not np.isfinite(f_lo) or not np.isfinite(f_hi):
            return None
        # f(lo) ≈ −x̄ < 0; f(∞) → b/2 − x̄ > 0
        if f_lo >= 0:
            return None
        if f_hi <= 0:
            for _ in range(20):
                hi *= 2.0
                f_hi = f(hi)
                if f_hi > 0:
                    break
            else:
                return None
        try:
            return float(brentq(f, lo, hi, maxiter=200))
        except ValueError:
            return None

    sigma = solve_sigma(b)
    if sigma is None:
        return None
    # One-step Hannon–Dahiya update of b (reproduces Cumming's 11.18 vs xmax=10.77)
    # b ← b + σ log(1 + (e^{b/σ} − 1)/n)
    b = b + sigma * np.log1p(np.expm1(b / sigma) / n)
    sigma2 = solve_sigma(b)
    if sigma2 is not None:
        sigma = sigma2
    zmax = t * np.exp(b)
    alpha = 1.0 + 1.0 / sigma  # p(z) ∝ z^{-α} on [t, zmax]
    ll = -n * (
        np.log(sigma) + xbar / sigma + np.log(-np.expm1(-b / sigma))
    )
    # log(1 − e^{-b/σ}) = log(-expm1(-b/σ))
    return {
        "model": "cumming_trunc_pareto",
        "sigma": float(sigma),
        "b": float(b),
        "alpha": float(alpha),
        "beta": float(-alpha),
        "xmin": t,
        "xmax": float(zmax),
        "xmax_obs": zmax_obs,
        "ll": float(ll),
        "n": n,
        "k_params": 2,
    }


# ── Power law with exponential cutoff: p(x) ∝ x^β exp(−x/xc), x ≥ xmin ───────
# Soft tail (Clauset-style), vs hard right-truncated Pareto above.


def _upper_inc_gamma(a, z):
    """Upper incomplete Γ(a, z) = ∫_z^∞ t^{a−1} e^{−t} dt, including a ≤ 0."""
    from scipy.special import exp1, gamma, gammaincc

    a = float(a)
    z = float(z)
    if not np.isfinite(z) or z <= 0:
        return np.nan
    if abs(a) < 1e-14:
        return float(exp1(z))
    if a > 0:
        g = float(gamma(a) * gammaincc(a, z))
        return g if np.isfinite(g) else np.nan
    val = _upper_inc_gamma(a + 1.0, z)
    if not np.isfinite(val):
        return np.nan
    return (val - np.exp(a * np.log(z) - z)) / a


def _pl_cutoff_logZ(beta, xc, xmin):
    """log ∫_{xmin}^{∞} x^β exp(−x/xc) dx = (β+1) log xc + log Γ(β+1, xmin/xc)."""
    from scipy.integrate import quad

    a = float(beta) + 1.0
    z0 = float(xmin) / float(xc)
    G = _upper_inc_gamma(a, z0)
    if np.isfinite(G) and G > 0:
        return a * np.log(float(xc)) + np.log(G)

    def f(x):
        return np.exp(beta * np.log(x) - x / xc)

    z, _ = quad(f, float(xmin), np.inf, limit=300, epsabs=1e-10, epsrel=1e-7)
    if not np.isfinite(z) or z <= 0:
        return np.nan
    return float(np.log(z))


def pdf_pl_cutoff(x, beta, xc, xmin):
    x = np.asarray(x, float)
    out = np.zeros_like(x, dtype=float)
    logz = _pl_cutoff_logZ(beta, xc, xmin)
    if not np.isfinite(logz):
        return out
    m = x >= xmin
    out[m] = np.exp(beta * np.log(x[m]) - x[m] / xc - logz)
    return out


def fit_pl_cutoff_mle(sizes, size_min=None):
    """MLE for p(x) ∝ x^β e^{−x/xc} on x ≥ xmin."""
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s > 0)]
    if size_min is not None:
        s = s[s >= size_min]
    else:
        size_min = float(s.min())
        s = s[s >= size_min]
    n = len(s)
    if n < 20:
        return None
    xmin = float(size_min)
    logx_sum = float(np.sum(np.log(s)))
    x_sum = float(np.sum(s))
    ols = fit_powerlaw_ols(s, size_min=xmin)
    beta0 = float(np.clip(ols["beta"] if ols is not None else -1.5, -4.0, -0.3))
    xc0 = float(max(np.quantile(s, 0.90), xmin * 10.0))

    def neg_ll(theta):
        beta, log_xc = float(theta[0]), float(theta[1])
        xc = float(np.exp(log_xc))
        if not np.isfinite(xc) or xc <= xmin * 0.5:
            return 1e300
        logz = _pl_cutoff_logZ(beta, xc, xmin)
        if not np.isfinite(logz):
            return 1e300
        return -( -n * logz + beta * logx_sum - x_sum / xc)

    bounds = [(-4.5, -0.15), (np.log(xmin * 2.0), np.log(float(s.max()) * 30.0))]
    best = None
    starts = (
        [beta0, np.log(xc0)],
        [beta0, np.log(max(float(np.median(s)) * 3.0, xmin * 10.0))],
        [-1.2, np.log(max(float(np.quantile(s, 0.95)), xmin * 20.0))],
    )
    for start in starts:
        res = minimize(neg_ll, start, method="L-BFGS-B", bounds=bounds)
        if res.success and (best is None or res.fun < best.fun):
            best = res
    if best is None:
        return None
    beta = float(best.x[0])
    xc = float(np.exp(best.x[1]))
    ll = -float(best.fun)
    return {
        "model": "pl_cutoff_mle",
        "beta": beta,
        "alpha": -beta,
        "xc": xc,
        "ll": ll,
        "aic": 2 * 2 - 2 * ll,
        "n": n,
        "xmin": xmin,
        "k_params": 2,
    }


# ── Truncated exponential on SIZE: UPPER truncation only ─────────────────────
# No lower truncation of small fires. Support starts at 0; cut large fires at xmax:
#   p(x) = λ exp(-λ x) / (1 - exp(-λ xmax)),   0 ≤ x ≤ xmax
# (If xmax is None: untruncated Exp from 0: p(x)=λ exp(-λ x), x≥0.)
# On log–log axes this curves down in the right tail (faster than power law).


def fit_truncexpon_mle(sizes, size_min=None):
    """MLE for upper-truncated exponential on size (support 0..xmax)."""
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s > 0)]
    if size_min is not None:
        s = s[s >= size_min]
    n = len(s)
    if n < 20:
        return None
    xmin = float(size_min) if size_min is not None else float(np.min(s))
    xmax = float(np.max(s))
    if xmax <= 0:
        return None
    # Untruncated MLE λ = 1/mean; refine for upper truncation on [0, xmax]
    lam = 1.0 / float(np.mean(s))
    for _ in range(80):
        t = lam * xmax
        if t > 50:
            corr = 0.0
        else:
            corr = xmax / (np.exp(t) - 1.0)
        # score: 1/λ - mean(x) - xmax/(e^{λ xmax}-1) = 0
        g = 1.0 / lam - float(np.mean(s)) - corr
        if t > 50:
            dg = -1.0 / lam**2
        else:
            e = np.exp(t)
            dg = -1.0 / lam**2 + (xmax**2 * e) / (e - 1.0) ** 2
        if abs(dg) < 1e-18:
            break
        lam_new = lam - g / dg
        if lam_new <= 0:
            lam_new = 0.5 * lam
        if abs(lam_new - lam) < 1e-12:
            lam = lam_new
            break
        lam = lam_new
    Z = 1.0 - np.exp(-lam * xmax)
    if Z <= 0:
        return None
    ll = n * np.log(lam) - lam * np.sum(s) - n * np.log(Z)
    aic = 2 * 1 - 2 * ll
    bic = np.log(n) * 1 - 2 * ll
    return {
        "model": "truncexpon_upper_mle",
        "lambda": lam,
        "sigma": 1.0 / lam,
        "xmax": xmax,
        "ll": ll,
        "aic": aic,
        "bic": bic,
        "n": n,
        "xmin": xmin,  # data filter only; model support starts at 0
        "k_params": 1,
    }


def pdf_truncexpon(x, lam, xmin=0.0, xmax=None):
    """
    Exponential on size from 0.
    xmax=None: p(x)=λ e^{-λx}, x≥0
    else:      p(x)=λ e^{-λx} / (1-e^{-λ xmax}), 0≤x≤xmax  (upper trunc. only)
    `xmin` kept for API compat; does not lower-truncate the density.
    """
    x = np.asarray(x, float)
    out = np.zeros_like(x, dtype=float)
    if xmax is None:
        m = x >= 0
        out[m] = lam * np.exp(-lam * x[m])
        return out
    m = (x >= 0) & (x <= xmax)
    Z = 1.0 - np.exp(-lam * xmax)
    if Z <= 0:
        return out
    out[m] = lam * np.exp(-lam * x[m]) / Z
    return out


def sf_truncexpon(x, lam, xmin=0.0, xmax=None):
    """Survival under upper-truncated (or untruncated) exponential on size."""
    x = np.asarray(x, float)
    if xmax is None:
        return np.where(x <= 0, 1.0, np.exp(-lam * x))
    Z = 1.0 - np.exp(-lam * xmax)
    xx = np.minimum(np.maximum(x, 0.0), xmax)
    sf = (np.exp(-lam * xx) - np.exp(-lam * xmax)) / Z
    return np.where(x < 0, 1.0, np.where(x >= xmax, 0.0, sf))


def fit_truncexpon_ols(sizes, n_bins=None, size_min=None):
    """
    OLS for upper-truncated exponential on size using log-binned PDF:
      log(p) ≈ const - λ x
    Fit log(pdf) vs x; λ = -slope. Model support [0, xmax], xmax=max(size).
    """
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s > 0)]
    if size_min is not None:
        s = s[s >= size_min]
    if len(s) < 20:
        return None
    xmin = float(size_min) if size_min is not None else float(np.min(s))
    xmax = float(np.max(s))
    x, y = log_binned_pdf(s, n_bins=n_bins)
    m = (x > 0) & (x <= xmax) & (y > 0)
    x, y = x[m], y[m]
    if len(x) < 5:
        return None
    log_y = np.log(y)
    slope, intercept, r, p, se = stats.linregress(x, log_y)
    lam = -slope
    if not np.isfinite(lam) or lam <= 0:
        return None
    return {
        "model": "truncexpon_upper_ols",
        "lambda": lam,
        "lambda_se": se,
        "sigma": 1.0 / lam,
        "intercept": intercept,
        "r2": r**2,
        "xmax": xmax,
        "n": len(s),
        "xmin": xmin,
    }


def pdf_truncexpon_ols(x, lam, xmin, xmax, intercept=None):
    """Plot upper-truncated Exp on size (normalized formula by default)."""
    if intercept is None:
        return pdf_truncexpon(x, lam, xmin=0.0, xmax=xmax)
    x = np.asarray(x, float)
    out = np.zeros_like(x, dtype=float)
    m = (x >= 0) & (x <= xmax)
    out[m] = np.exp(intercept - lam * x[m])
    return out


# ── Weibull / Reed–McKelvey Model I (paper form, log writing) ───────────────
# Paper: f(A)=α A^{-β} exp{ -α/(1-β) (A^{1-β} - A0^{1-β}) }, A≥A0
# Log:   ln f = ln(α) - β ln(A) - α/(1-β) (A^{1-β} - A0^{1-β})
# With A0=1 km²: ln f = ln(α) - β ln(A) - α/(1-β) (A^{1-β} - 1)
# (Handwritten [ln A]^{1-β} was the log-PDF rewrite intent; the paper term is A^{1-β}.)


def ln_pdf_weibull_rm(A, alpha, beta, A0):
    """Reed–McKelvey Model I log-density."""
    A = np.asarray(A, float)
    A0 = float(A0)
    alpha = float(alpha)
    beta = float(beta)
    if np.isclose(beta, 1.0):
        # b→1: f(A)=α A^{-1} exp{-α(ln A - ln A0)} = α A0^α A^{-(α+1)}
        return np.log(alpha) + alpha * np.log(A0) - (alpha + 1.0) * np.log(A)
    return (
        np.log(alpha)
        - beta * np.log(A)
        - (alpha / (1.0 - beta)) * (np.power(A, 1.0 - beta) - np.power(A0, 1.0 - beta))
    )


def pdf_weibull(x, alpha, beta, xmin):
    """Reed–McKelvey Model I PDF on A ≥ xmin (= A0)."""
    x = np.asarray(x, float)
    out = np.full_like(x, np.nan, dtype=float)
    m = np.isfinite(x) & (x >= xmin) & (alpha > 0)
    if np.any(m):
        out[m] = np.exp(ln_pdf_weibull_rm(x[m], alpha, beta, xmin))
    return out


def sf_weibull(x, alpha, beta, xmin):
    """Survival S(x)=exp{-α/(1-β)(x^{1-β}-A0^{1-β})}."""
    x = np.asarray(x, float)
    out = np.ones_like(x, dtype=float)
    m = np.isfinite(x) & (x >= xmin)
    if not np.any(m):
        return out
    if np.isclose(beta, 1.0):
        out[m] = np.power(xmin / x[m], alpha)
        return out
    cum = (alpha / (1.0 - beta)) * (
        np.power(x[m], 1.0 - beta) - np.power(xmin, 1.0 - beta)
    )
    out[m] = np.exp(-cum)
    return out


def fit_weibull_ols(sizes, n_bins=None, size_min=None):
    """
    OLS Weibull (Reed–McKelvey) on log-binned PDF:
      ln f(A) = ln(α) - β ln(A) - α/(1-β) (A^{1-β} - A0^{1-β}),  A0=size_min
    Tries β<1 and β>1; keeps higher R².
    """
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s)]
    if size_min is not None:
        s = s[s >= size_min]
    else:
        size_min = float(s.min())
        s = s[s >= size_min]
    if len(s) < 20:
        return None
    xmin = float(size_min)
    x, y = log_binned_pdf(s, n_bins=n_bins)
    msk = (x >= xmin) & (y > 0)
    x, y = x[msk], y[msk]
    if len(x) < 6:
        return None

    def model_ln(xv, alpha, beta):
        return ln_pdf_weibull_rm(xv, alpha, beta, xmin)

    b_abs = 1.2
    try:
        pl = fit_powerlaw_ols(s, size_min=xmin)
        if pl is not None:
            b_abs = float(np.clip(abs(pl["beta"]), 0.2, 2.5))
    except Exception:
        pass

    candidates = []
    starts = [
        ([1e-8, 0.05], [1e4, 0.95], min(b_abs, 0.8)),
        ([1e-8, 1.05], [1e4, 3.0], max(b_abs, 1.2)),
    ]
    for lo, hi, beta0 in starts:
        alpha0 = float(
            np.clip(np.exp(np.mean(np.log(y) + beta0 * np.log(x))), 1e-6, 1e3)
        )
        try:
            popt, pcov = curve_fit(
                model_ln,
                x,
                np.log(y),
                p0=[alpha0, beta0],
                bounds=(lo, hi),
                maxfev=30000,
            )
            alpha, beta = float(popt[0]), float(popt[1])
            yhat = model_ln(x, alpha, beta)
            lny = np.log(y)
            ss_res = np.sum((lny - yhat) ** 2)
            ss_tot = np.sum((lny - np.mean(lny)) ** 2)
            r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else -np.inf
            se_b = np.nan
            try:
                se_b = float(np.sqrt(max(float(pcov[1, 1]), 0.0)))
            except Exception:
                pass
            if np.isfinite(r2):
                candidates.append((r2, alpha, beta, se_b))
        except Exception:
            continue

    if not candidates:
        return None
    r2, alpha, beta, se_b = max(candidates, key=lambda t: t[0])
    return {
        "model": "weibull_ols",
        "alpha": alpha,
        "beta": beta,
        "beta_se": se_b,
        "shape": beta,
        "scale": alpha,
        "r2": float(r2),
        "n": len(s),
        "xmin": xmin,
    }


def fit_weibull_mle(sizes, size_min=None):
    """MLE for Reed–McKelvey Model I (α, β), A ≥ xmin."""
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s)]
    if size_min is not None:
        s = s[s >= size_min]
    else:
        size_min = float(s.min())
        s = s[s >= size_min]
    n = len(s)
    if n < 20:
        return None
    xmin = float(size_min)

    from scipy.optimize import minimize

    def nll(theta):
        alpha, beta = theta
        if alpha <= 0 or beta <= 0 or np.isclose(beta, 1.0):
            return 1e300
        ll = np.sum(ln_pdf_weibull_rm(s, alpha, beta, xmin))
        if not np.isfinite(ll):
            return 1e300
        return -ll

    candidates = []
    for lo, hi, x0 in [
        ([1e-8, 0.05], [1e4, 0.99], [0.1, 0.6]),
        ([1e-8, 1.01], [1e4, 3.0], [0.1, 1.3]),
    ]:
        try:
            res = minimize(nll, x0=x0, bounds=list(zip(lo, hi)), method="L-BFGS-B")
            if res.success and np.isfinite(res.fun):
                candidates.append(res)
        except Exception:
            continue
    if not candidates:
        return None
    res = min(candidates, key=lambda r: r.fun)
    alpha, beta = float(res.x[0]), float(res.x[1])
    ll = float(-res.fun)
    aic = 2 * 2 - 2 * ll
    bic = np.log(n) * 2 - 2 * ll
    return {
        "model": "weibull",
        "alpha": alpha,
        "beta": beta,
        "shape": beta,
        "scale": alpha,
        "ll": ll,
        "aic": aic,
        "bic": bic,
        "n": n,
        "xmin": xmin,
        "k_params": 2,
    }




# ── Lognormal via quadratic in ln(A) (log-log parabola) ───────────────────────
# ln(f(A)) = ln(a) - β ln(A) - ψ [ln(A)]²


# ── MLE on original-space densities (user equations) ─────────────────────────
# Power law:   f(A) = α A^{-β}
# Lognormal:   f(A) = α exp(-β ln A - ψ [ln A]²)
# Weibull:     f(A) = α A^{-β} exp(-α/(1-β) [ln A]^{1-β})
# MLE uses these kernels normalized on A ≥ xmin (Weibull needs A > 1 for real [ln A]^{1-β}).


def _trapz_logspace(log_unnorm, xmin, decades=8, n=6000):
    """∫_{xmin}^{xmin*10^decades} exp(log_unnorm(A)) dA via u=ln A."""
    umin = np.log(max(xmin, 1e-12))
    umax = umin + decades * np.log(10.0)
    u = np.linspace(umin, umax, n)
    A = np.exp(u)
    log_int = np.asarray(log_unnorm(A), float) + u  # +u for Jacobian
    m = np.isfinite(log_int)
    if not np.any(m):
        return np.nan
    c = float(np.max(log_int[m]))
    # keep exp(c) in log-space friendly range
    if c > 700:
        log_int = log_int - (c - 700)
        c = 700.0
    elif c < -700:
        return 0.0
    integ = float(np.trapz(np.exp(log_int[m] - c), u[m]))
    if not np.isfinite(integ) or integ <= 0:
        return np.nan
    return float(integ * np.exp(c))


def fit_powerlaw_mle_orig(sizes, size_min=None):
    """MLE for f(A)=α A^{-β} on [xmin,∞), β>1, α=(β-1) xmin^{β-1}."""
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s > 0)]
    if size_min is not None:
        s = s[s >= size_min]
    else:
        size_min = float(s.min())
        s = s[s >= size_min]
    n = len(s)
    if n < 20:
        return None
    xmin = float(size_min)
    beta = 1.0 + n / np.sum(np.log(s / xmin))  # Clauset α ≡ β here
    if not np.isfinite(beta) or beta <= 1.0:
        return None
    alpha = (beta - 1.0) * xmin ** (beta - 1.0)
    ll = n * np.log(beta - 1.0) - n * np.log(xmin) - beta * np.sum(np.log(s / xmin))
    return {
        "model": "powerlaw_mle_orig",
        "alpha": float(alpha),
        "beta": float(beta),
        "ll": float(ll),
        "aic": float(2 * 1 - 2 * ll),
        "n": n,
        "xmin": xmin,
        "k_params": 1,
    }


def pdf_powerlaw_orig(A, alpha, beta):
    A = np.asarray(A, float)
    return alpha * np.power(A, -beta)


def fit_lognormal_mle_orig(sizes, size_min=None):
    """
    MLE for f(A)=α exp(-β ln A - ψ [ln A]²), α from normalization on [xmin,∞).
    Free parameters: β, ψ>0.
    """
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s > 0)]
    if size_min is not None:
        s = s[s >= size_min]
    else:
        size_min = float(s.min())
        s = s[s >= size_min]
    n = len(s)
    if n < 20:
        return None
    xmin = float(size_min)
    ln = np.log(s)
    # method-of-moments start on ln A
    mu0 = float(np.mean(ln))
    v0 = float(max(np.var(ln), 1e-3))
    psi0 = 0.5 / v0
    beta0 = 1.0 - mu0 / v0  # matches lognormal β≈1-μ/σ² when ψ=1/(2σ²)

    def pack_ll(theta):
        beta, psi = theta
        if psi <= 0:
            return 1e300

        def log_u(A):
            lnA = np.log(A)
            return -beta * lnA - psi * lnA**2

        Z = _trapz_logspace(log_u, xmin)
        if not np.isfinite(Z) or Z <= 0:
            return 1e300
        ll = np.sum(log_u(s)) - n * np.log(Z)
        return -ll if np.isfinite(ll) else 1e300

    res = minimize(
        pack_ll,
        x0=[beta0, psi0],
        bounds=[(None, None), (1e-6, 50.0)],
        method="L-BFGS-B",
    )
    if not res.success and res.fun > 1e299:
        return None
    beta, psi = float(res.x[0]), float(res.x[1])

    def log_u(A):
        lnA = np.log(A)
        return -beta * lnA - psi * lnA**2

    Z = _trapz_logspace(log_u, xmin)
    alpha = 1.0 / Z
    ll = float(-res.fun)
    return {
        "model": "lognormal_mle_orig",
        "alpha": float(alpha),
        "beta": beta,
        "psi": psi,
        "ll": ll,
        "aic": float(2 * 2 - 2 * ll),
        "n": n,
        "xmin": xmin,
        "k_params": 2,
    }


def pdf_lognormal_orig(A, alpha, beta, psi):
    A = np.asarray(A, float)
    lnA = np.log(A)
    return alpha * np.exp(-beta * lnA - psi * lnA**2)


def fit_weibull_mle_orig(sizes, size_min=None):
    """
    MLE for f(A)=α A^{-β} exp(-α/(1-β) [ln A]^{1-β}) as written.
    Requires A > 1 so [ln A]^{1-β} is real → xmin_eff = max(xmin, 1+ε).
    Treats the written f as the density kernel and renormalizes on [xmin_eff, ∞)
    with free (β, γ) where γ=α/|1-β| controls curvature; α recovered after.
    """
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s > 0)]
    if size_min is not None:
        s = s[s >= size_min]
    xmin_in = float(size_min) if size_min is not None else float(s.min())
    xmin = max(xmin_in, 1.0 + 1e-9)  # ln A > 0
    s = s[s >= xmin]
    n = len(s)
    if n < 20:
        return None

    # Parameterize: β and γ>0 with
    #   ũ(A) = A^{-β} exp(-γ [ln A]^{1-β})   (β≠1)
    # then α from matching γ = α/|1-β| with sign of (1-β), and Z-normalize.
    def nll(theta):
        beta, gamma = theta
        if gamma <= 0 or beta <= 0 or np.isclose(beta, 1.0):
            return 1e300

        def log_u(A):
            lnA = np.log(A)
            return -beta * lnA - gamma * np.power(lnA, 1.0 - beta)

        Z = _trapz_logspace(log_u, xmin)
        if not np.isfinite(Z) or Z <= 0:
            return 1e300
        ll = np.sum(log_u(s)) - n * np.log(Z)
        return -ll if np.isfinite(ll) else 1e300

    candidates = []
    for x0, bounds in [
        ([0.7, 0.5], [(0.05, 0.99), (1e-4, 20.0)]),
        ([1.4, 0.5], [(1.01, 2.8), (1e-4, 20.0)]),
    ]:
        res = minimize(nll, x0=x0, bounds=bounds, method="L-BFGS-B")
        if np.isfinite(res.fun) and res.fun < 1e299:
            candidates.append(res)
    if not candidates:
        return None
    res = min(candidates, key=lambda r: r.fun)
    beta, gamma = float(res.x[0]), float(res.x[1])
    # recover α so that α/|1-β| = γ and prefactor matches after normalization
    alpha_kernel = gamma * abs(1.0 - beta)

    def log_u(A):
        lnA = np.log(A)
        return -beta * lnA - gamma * np.power(lnA, 1.0 - beta)

    Z = _trapz_logspace(log_u, xmin)
    ll = float(-res.fun)
    return {
        "model": "weibull_mle_orig",
        "alpha": float(alpha_kernel),
        "gamma": float(gamma),
        "beta": beta,
        "Z": float(Z),
        "ll": ll,
        "aic": float(2 * 2 - 2 * ll),
        "n": n,
        "xmin": xmin,
        "xmin_data": xmin_in,
        "k_params": 2,
    }


def pdf_weibull_orig(A, alpha, beta, Z=1.0, gamma=None):
    """Normalized Weibull PDF from user kernel."""
    A = np.asarray(A, float)
    lnA = np.log(np.maximum(A, 1.0 + 1e-12))
    if gamma is None:
        gamma = alpha / abs(1.0 - beta)
    u = np.power(A, -beta) * np.exp(-gamma * np.power(lnA, 1.0 - beta))
    return u / Z


def fit_lognormal_quad_ols(sizes, n_bins=None, size_min=None):
    """OLS: ln(f) = ln(a) - β ln(A) - ψ [ln(A)]² on log-binned PDF."""
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s > 0)]
    if size_min is not None:
        s = s[s >= size_min]
    if len(s) < 20:
        return None
    x, y = log_binned_pdf(s, n_bins=n_bins)
    if len(x) < 6:
        return None
    lnA = np.log(x)
    lnf = np.log(y)
    Xmat = np.column_stack([np.ones_like(lnA), lnA, lnA**2])
    coef, _, _, _ = np.linalg.lstsq(Xmat, lnf, rcond=None)
    ln_a, c1, c2 = coef
    beta = -c1
    psi = -c2
    yhat = Xmat @ coef
    ss_res = np.sum((lnf - yhat) ** 2)
    ss_tot = np.sum((lnf - np.mean(lnf)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    nobs = len(lnf)
    mse = ss_res / max(nobs - 3, 1)
    try:
        cov = mse * np.linalg.inv(Xmat.T @ Xmat)
        beta_se = float(np.sqrt(max(cov[1, 1], 0.0)))
    except np.linalg.LinAlgError:
        beta_se = np.nan
    mu = sigma = np.nan
    if psi > 0:
        sigma = float(np.sqrt(1.0 / (2.0 * psi)))
        mu = float((1.0 - beta) * sigma**2)
    return {
        "model": "lognormal_quad_ols",
        "ln_a": float(ln_a),
        "a": float(np.exp(ln_a)),
        "beta": float(beta),
        "psi": float(psi),
        "beta_se": beta_se,
        "mu": mu,
        "sigma": sigma,
        "r2": float(r2),
        "n": len(s),
        "xmin": float(size_min) if size_min is not None else float(s.min()),
    }


def pdf_lognormal_quad(A, ln_a, beta, psi):
    A = np.asarray(A, float)
    lnA = np.log(A)
    return np.exp(ln_a - beta * lnA - psi * lnA**2)


# ── bootstrap for key params ─────────────────────────────────────────────────
def bootstrap_param(sizes, fit_fn, key, n_iter=500, size_min=None, seed=0):
    rng = np.random.default_rng(seed)
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s)]
    if size_min is not None:
        s = s[s >= size_min]
    vals = []
    for _ in range(n_iter):
        sample = rng.choice(s, size=len(s), replace=True)
        res = fit_fn(sample, size_min=size_min)
        if res is not None and np.isfinite(res.get(key, np.nan)):
            vals.append(res[key])
    if len(vals) < 20:
        return np.nan, np.nan, np.nan
    return float(np.mean(vals)), *np.percentile(vals, [2.5, 97.5])


def load_groups():
    """Return dict dataset -> firetype -> sizes array (Fig1B sources)."""
    # Completeness thresholds used in main figures
    xmin_by_ds = {
        "CalFire": 1.0,
        "MTBS": 4.0,
        "Atlas": 0.21,
        "FIRED": 0.21,
    }
    out = {}
    for ds, xmin in xmin_by_ds.items():
        gdf = gpd.read_file(BASE / f"dataPrc/firePrmt/{ds}.shp")
        # firePrmt uses WUI; some IgnWUI layers use Urban-edge
        ft = gdf["FireType"].replace({"WUI": "Urban-edge"})
        size_col = "size_km2" if "size_km2" in gdf.columns else "size"
        out[ds] = {
            "Urban-edge": gdf.loc[ft == "Urban-edge", size_col].to_numpy(dtype=float),
            "Wildland": gdf.loc[ft == "Wildland", size_col].to_numpy(dtype=float),
            "xmin": xmin,
            "meta": gdf.assign(FireTypePlot=ft),
        }
    return out


def fit_all_for_group(sizes, xmin):
    """Return OLS fits for plotting + MLE for AIC."""
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s >= xmin)]
    ols = fit_powerlaw_ols(s, size_min=xmin)
    pl = fit_powerlaw_mle(s, size_min=xmin)
    tpl_ols = fit_truncpower_ols(s, size_min=xmin)
    tpl_mle = fit_truncpower_mle(s, size_min=xmin)
    te_ols = fit_truncexpon_ols(s, size_min=xmin)
    te_mle = fit_truncexpon_mle(s, size_min=xmin)
    ln_ols = fit_lognormal_quad_ols(s, size_min=xmin)
    wb_ols = fit_weibull_ols(s, size_min=xmin)
    wb_mle = fit_weibull_mle(s, size_min=xmin)
    return {
        "powerlaw_ols": ols,
        "powerlaw_mle": pl,
        "truncpower_ols": tpl_ols,
        "truncpower_mle": tpl_mle,
        "truncexpon_ols": te_ols,
        "truncexpon_mle": te_mle,
        "lognormal_ols": ln_ols,
        "weibull_ols": wb_ols,
        "weibull_mle": wb_mle,
    }


def extreme_probs(sizes, xmin, thresholds=(100.0, 1000.0)):
    """Empirical + model P(X > t) for each fit (OLS curves + MLE power law)."""
    s = np.asarray(sizes, float)
    s = s[np.isfinite(s) & (s >= xmin)]
    n = len(s)
    fits = fit_all_for_group(s, xmin)
    pl = fits["powerlaw_mle"]
    te = fits["truncexpon_ols"] or fits["truncexpon_mle"]
    wb = fits["weibull_ols"] or fits["weibull_mle"]
    rows = []
    for t in thresholds:
        emp = float(np.mean(s > t))
        row = {"threshold_km2": t, "empirical": emp, "n": n, "xmin": xmin}
        if pl:
            row["powerlaw_mle"] = float(sf_powerlaw(t, pl["alpha"], xmin)) if t >= xmin else 1.0
        if te:
            row["truncexpon"] = float(
                sf_truncexpon(t, te["lambda"], xmin=0.0, xmax=te.get("xmax"))
            )
        if wb:
            row["weibull"] = (
                float(sf_weibull(t, wb["alpha"], wb["beta"], xmin)) if t >= xmin else 1.0
            )
        rows.append(row)
    return rows, fits


def main():
    data = load_groups()
    summary_rows = []
    extreme_rows = []

    # ── Fit summary tables ───────────────────────────────────────────────────
    for ds, d in data.items():
        xmin = d["xmin"]
        for ft in ["Urban-edge", "Wildland"]:
            sizes = d[ft]
            fits = fit_all_for_group(sizes, xmin)
            ols, pl = fits["powerlaw_ols"], fits["powerlaw_mle"]
            te_ols, te_mle = fits["truncexpon_ols"], fits["truncexpon_mle"]
            wb_ols, wb_mle = fits["weibull_ols"], fits["weibull_mle"]
            # bootstrap key params (moderate iters for speed)
            b_mean, b_lo, b_hi = bootstrap_param(
                sizes, fit_powerlaw_ols, "beta", n_iter=400, size_min=xmin, seed=1
            )
            a_mean, a_lo, a_hi = bootstrap_param(
                sizes, fit_powerlaw_mle, "alpha", n_iter=400, size_min=xmin, seed=2
            )
            l_mean, l_lo, l_hi = bootstrap_param(
                sizes, fit_truncexpon_ols, "lambda", n_iter=400, size_min=xmin, seed=3
            )
            k_mean, k_lo, k_hi = bootstrap_param(
                sizes, fit_weibull_ols, "shape", n_iter=400, size_min=xmin, seed=4
            )

            row = {
                "dataset": ds,
                "FireType": ft,
                "n": int(np.sum(np.isfinite(sizes) & (sizes >= xmin))),
                "xmin": xmin,
                "beta_ols": ols["beta"] if ols else np.nan,
                "beta_ols_se": ols["beta_se"] if ols else np.nan,
                "beta_ols_ci_lo": b_lo,
                "beta_ols_ci_hi": b_hi,
                "alpha_mle": pl["alpha"] if pl else np.nan,
                "alpha_mle_ci_lo": a_lo,
                "alpha_mle_ci_hi": a_hi,
                "beta_mle": pl["beta"] if pl else np.nan,
                # size-exponential λ; smaller λ => heavier tails
                "lambda_ols": te_ols["lambda"] if te_ols else np.nan,
                "lambda_ols_ci_lo": l_lo,
                "lambda_ols_ci_hi": l_hi,
                "lambda_mle": te_mle["lambda"] if te_mle else np.nan,
                "xmax_te": te_ols["xmax"] if te_ols else np.nan,
                "weibull_k_ols": wb_ols["shape"] if wb_ols else np.nan,
                "weibull_k_ols_ci_lo": k_lo,
                "weibull_k_ols_ci_hi": k_hi,
                "weibull_scale_ols": wb_ols["scale"] if wb_ols else np.nan,
                "weibull_k_mle": wb_mle["shape"] if wb_mle else np.nan,
                "weibull_scale_mle": wb_mle["scale"] if wb_mle else np.nan,
                # used by param/decade panels
                "lambda_exp": te_ols["lambda"] if te_ols else np.nan,
                "lambda_ci_lo": l_lo,
                "lambda_ci_hi": l_hi,
                "weibull_k": wb_ols["shape"] if wb_ols else np.nan,
                "weibull_k_ci_lo": k_lo,
                "weibull_k_ci_hi": k_hi,
                "weibull_scale": wb_ols["scale"] if wb_ols else np.nan,
                "aic_powerlaw": pl["aic"] if pl else np.nan,
                "aic_truncexpon": te_mle["aic"] if te_mle else np.nan,
                "aic_weibull": wb_mle["aic"] if wb_mle else np.nan,
                "bic_powerlaw": pl["bic"] if pl else np.nan,
                "bic_truncexpon": te_mle["bic"] if te_mle else np.nan,
                "bic_weibull": wb_mle["bic"] if wb_mle else np.nan,
                "ll_powerlaw": pl["ll"] if pl else np.nan,
                "ll_truncexpon": te_mle["ll"] if te_mle else np.nan,
                "ll_weibull": wb_mle["ll"] if wb_mle else np.nan,
                "r2_powerlaw_ols": ols["r2"] if ols else np.nan,
                "r2_truncexpon_ols": te_ols["r2"] if te_ols else np.nan,
                "r2_weibull_ols": wb_ols["r2"] if wb_ols else np.nan,
            }
            # best by AIC
            aics = {
                "powerlaw_mle": row["aic_powerlaw"],
                "truncexpon": row["aic_truncexpon"],
                "weibull": row["aic_weibull"],
            }
            row["best_aic"] = min(aics, key=aics.get)
            summary_rows.append(row)

            erows, *_ = extreme_probs(sizes, xmin)
            for er in erows:
                er.update({"dataset": ds, "FireType": ft})
                extreme_rows.append(er)

            print(
                f"{ds:8s} {ft:12s} n={row['n']:4d}  "
                f"β_OLS={row['beta_ols']:.3f}  "
                f"λ_TE={row['lambda_ols']:.4f}  "
                f"k_OLS={row['weibull_k_ols']:.3f}  "
                f"R2[PL/TE/WB]={row['r2_powerlaw_ols']:.2f}/"
                f"{row['r2_truncexpon_ols']:.2f}/{row['r2_weibull_ols']:.2f}"
            )

    summary = pd.DataFrame(summary_rows)
    extremes = pd.DataFrame(extreme_rows)
    summary.to_csv(OUT_PRC / "size_dist_alt_fits.csv", index=False)
    extremes.to_csv(OUT_PRC / "size_dist_alt_extremes.csv", index=False)

    # ── Figure A: PDF overlays (4 sources), Fig1B layout ──────────────────────
    from matplotlib.lines import Line2D

    datasets = ["CalFire", "MTBS", "Atlas", "FIRED"]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6), facecolor="w", sharex=True, sharey=True)
    for ax, ds in zip(axes.ravel(), datasets):
        d = data[ds]
        xmin = d["xmin"]
        for ft in ["Urban-edge", "Wildland"]:
            s = d[ft]
            s = s[np.isfinite(s) & (s >= xmin)]
            if len(s) < 20:
                continue
            x, y = log_binned_pdf(s)
            c = COLOR[ft]
            ax.loglog(x, y, "o", color=c, ms=2.5, zorder=2)

            fits = fit_all_for_group(s, xmin)
            ols = fits["powerlaw_ols"]
            ln_ols = fits["lognormal_ols"]
            wb_ols = fits["weibull_ols"]
            xfit = np.logspace(np.log10(xmin), np.log10(s.max()), 400)
            if ols:
                ax.loglog(
                    xfit,
                    pdf_powerlaw_ols(xfit, ols["beta"], ols["intercept"]),
                    LS["powerlaw_ols"],
                    color=c,
                    lw=1.6,
                    alpha=0.95,
                    zorder=4,
                )
            if ln_ols:
                ax.loglog(
                    xfit,
                    pdf_lognormal_quad(
                        xfit, ln_ols["ln_a"], ln_ols["beta"], ln_ols["psi"]
                    ),
                    ls=LS["lognormal_ols"],
                    color=c,
                    lw=1.8,
                    alpha=0.95,
                    zorder=5,
                )
            if wb_ols:
                ax.loglog(
                    xfit,
                    pdf_weibull(xfit, wb_ols["alpha"], wb_ols["beta"], xmin),
                    LS["weibull_ols"],
                    color=c,
                    lw=1.5,
                    alpha=0.9,
                    zorder=3,
                )

        ax.set_xlim(1e-1, 1e5)
        ax.set_ylim(1e-7, 1e1)
        ax.set_xticks([1e-1, 1e1, 1e3, 1e5])
        ax.set_yticks([1e-7, 1e-4, 1e-1])
        ax.tick_params(which="both", direction="out")
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        ax.text(
            0.05,
            0.05,
            ds,
            transform=ax.transAxes,
            fontweight="bold",
            ha="left",
            va="bottom",
            fontsize=11,
        )

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color=COLOR["Urban-edge"],
            ls="-",
            lw=1.4,
            ms=5,
            label="WUI",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=COLOR["Wildland"],
            ls="-",
            lw=1.4,
            ms=5,
            label="Wildland",
        ),
        Line2D([0], [0], color="k", ls="-", lw=1.6, label="Power law"),
        Line2D([0], [0], color="k", ls=LS["lognormal_ols"], lw=1.8, label="Lognormal"),
        Line2D([0], [0], color="k", ls=":", lw=1.5, label="Weibull"),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=5,
        frameon=False,
        fontsize=9,
        bbox_to_anchor=(0.5, 1.02),
    )
    fig.supxlabel("Fire size (km$^2$)", fontsize=11)
    fig.supylabel("Probability density", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT_FIG / "Fig_size_dist_alt_PDF.png", dpi=300, bbox_inches="tight", facecolor="w")
    fig.savefig(OUT_PRC / "Fig_size_dist_alt_PDF.png", dpi=300, bbox_inches="tight", facecolor="w")
    print("Saved PDF overlay figure")

    # ── Figure B: comparable "tail heaviness" params by model ────────────────
    # Map each model to a signed metric where MORE POSITIVE / closer-to-zero
    # means heavier tails (like β in the paper):
    #   powerlaw: β_OLS
    #   truncexpon: -λ  (smaller λ → heavier tail)
    #   weibull: -k     (smaller k → heavier tail), or report k with reverse axis
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), facecolor="w")
    panels = [
        ("beta_ols", "beta_ols_ci_lo", "beta_ols_ci_hi", r"Power-law $\beta$ (OLS)", False),
        ("lambda_exp", "lambda_ci_lo", "lambda_ci_hi", r"Trunc. exp. $\lambda$ (size)", True),
        ("weibull_k", "weibull_k_ci_lo", "weibull_k_ci_hi", r"Weibull shape $k$", True),
    ]
    for ax, (key, lo, hi, title, invert) in zip(axes, panels):
        for i, ds in enumerate(["CalFire", "Atlas"]):
            sub = summary[summary["dataset"] == ds]
            y = i
            for ft, dy in [("Urban-edge", -0.12), ("Wildland", 0.12)]:
                r = sub[sub["FireType"] == ft].iloc[0]
                x = r[key]
                xerr = np.array([[x - r[lo]], [r[hi] - x]])
                ax.errorbar(
                    x,
                    y + dy,
                    xerr=xerr,
                    fmt="o",
                    color=COLOR[ft],
                    ms=6,
                    lw=1.4,
                    capsize=3,
                    label=ft if i == 0 else None,
                )
            # connector
            u = sub[sub["FireType"] == "Urban-edge"].iloc[0][key]
            w = sub[sub["FireType"] == "Wildland"].iloc[0][key]
            ax.plot([u, w], [y - 0.12, y + 0.12], "-", color="0.6", lw=0.8, zorder=0)
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["CalFire", "Atlas"])
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(title.split("(")[0].strip() if False else "")
        ax.axvline(np.nan, color="none")
        ax.grid(axis="x", color="0.9", lw=0.8)
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        if invert:
            ax.text(
                0.98,
                0.02,
                "← heavier tails",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=8,
                color="0.4",
            )
        else:
            ax.text(
                0.98,
                0.02,
                "heavier tails →",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=8,
                color="0.4",
            )
    axes[0].legend(frameon=False, fontsize=9, loc="best")
    # clean xlabels
    axes[0].set_xlabel(r"$\beta$")
    axes[1].set_xlabel(r"$\lambda$")
    axes[2].set_xlabel(r"$k$")
    fig.suptitle("WUI vs wildland under alternative size-distribution models", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_FIG / "Fig_size_dist_alt_params.png", dpi=300, bbox_inches="tight", facecolor="w")
    fig.savefig(OUT_PRC / "Fig_size_dist_alt_params.png", dpi=300, bbox_inches="tight", facecolor="w")
    print("Saved parameter comparison figure")

    # ── Figure C: P(X>100) and P(X>1000) empirical + models ──────────────────
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), facecolor="w")
    for ax, thr in zip(axes, [100.0, 1000.0]):
        sub = extremes[extremes["threshold_km2"] == thr]
        models = ["empirical", "powerlaw_mle", "truncexpon", "weibull"]
        model_labs = ["Empirical", "Power law", "Trunc. exp.", "Weibull"]
        x = np.arange(len(models))
        width = 0.18
        idx = 0
        for ds in ["CalFire", "Atlas"]:
            for ft in ["Urban-edge", "Wildland"]:
                r = sub[(sub["dataset"] == ds) & (sub["FireType"] == ft)].iloc[0]
                vals = [r[m] for m in models]
                offset = (idx - 1.5) * width
                ax.bar(
                    x + offset,
                    vals,
                    width=width,
                    color=COLOR[ft],
                    alpha=0.55 if ds == "Atlas" else 0.95,
                    edgecolor=COLOR[ft],
                    label=f"{ds} {('WUI' if ft=='Urban-edge' else ft)}",
                )
                idx += 1
        ax.set_xticks(x)
        ax.set_xticklabels(model_labs, fontsize=9)
        ax.set_ylabel(f"P(size > {thr:.0f} km$^2$)")
        ax.set_title(f"Threshold = {thr:.0f} km$^2$", fontsize=11)
        ax.set_yscale("log")
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
    axes[0].legend(frameon=False, fontsize=7, ncol=2, loc="upper right")
    fig.tight_layout()
    fig.savefig(OUT_FIG / "Fig_size_dist_alt_extremes.png", dpi=300, bbox_inches="tight", facecolor="w")
    fig.savefig(OUT_PRC / "Fig_size_dist_alt_extremes.png", dpi=300, bbox_inches="tight", facecolor="w")
    print("Saved extremes figure")

    # ── Figure D: CalFire decade trends under each model (main Fig1A style) ─
    cal = data["CalFire"]["meta"]
    cal = cal.assign(
        FireTypePlot=cal["FireType"].map({"WUI": "Urban-edge", "Wildland": "Wildland"})
    )
    decades = ["1990s", "2000s", "2010s", "2020s"]
    xmin = 1.0
    dec_rows = []
    for dec in decades:
        for ft in ["Urban-edge", "Wildland"]:
            s = cal.loc[(cal["decade"] == dec) & (cal["FireTypePlot"] == ft), "size_km2"].to_numpy()
            fits = fit_all_for_group(s, xmin)
            ols, pl = fits["powerlaw_ols"], fits["powerlaw_mle"]
            te, wb = fits["truncexpon_ols"], fits["weibull_ols"]
            te_mle, wb_mle = fits["truncexpon_mle"], fits["weibull_mle"]
            if ols is None or te is None or wb is None:
                continue
            _, blo, bhi = bootstrap_param(s, fit_powerlaw_ols, "beta", n_iter=300, size_min=xmin)
            _, llo, lhi = bootstrap_param(s, fit_truncexpon_ols, "lambda", n_iter=300, size_min=xmin)
            _, klo, khi = bootstrap_param(s, fit_weibull_ols, "shape", n_iter=300, size_min=xmin)
            dec_rows.append(
                {
                    "decade": dec,
                    "FireType": ft,
                    "n": int(np.sum(np.isfinite(s) & (s >= xmin))),
                    "beta_ols": ols["beta"],
                    "beta_lo": blo,
                    "beta_hi": bhi,
                    "lambda_exp": te["lambda"],
                    "lambda_lo": llo,
                    "lambda_hi": lhi,
                    "weibull_k": wb["shape"],
                    "k_lo": klo,
                    "k_hi": khi,
                    "aic_powerlaw": pl["aic"] if pl else np.nan,
                    "aic_truncexpon": te_mle["aic"] if te_mle else np.nan,
                    "aic_weibull": wb_mle["aic"] if wb_mle else np.nan,
                }
            )
            print(
                f"  {dec} {ft}: β_OLS={ols['beta']:.3f}, "
                f"λ_TE={te['lambda']:.4f}, k_OLS={wb['shape']:.3f}"
            )

    dec_df = pd.DataFrame(dec_rows)
    dec_df.to_csv(OUT_PRC / "size_dist_alt_by_decade_CalFire.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), facecolor="w")
    specs = [
        ("beta_ols", "beta_lo", "beta_hi", r"Power-law $\beta$ (OLS)"),
        ("lambda_exp", "lambda_lo", "lambda_hi", r"Trunc. exp. $\lambda$ (size)"),
        ("weibull_k", "k_lo", "k_hi", r"Weibull shape $k$"),
    ]
    x = np.arange(len(decades))
    for ax, (key, lo, hi, title) in zip(axes, specs):
        for ft, marker in [("Urban-edge", "o"), ("Wildland", "^")]:
            sub = dec_df[dec_df["FireType"] == ft]
            # align decade order
            sub = sub.set_index("decade").reindex(decades)
            y = sub[key].to_numpy()
            yerr = np.vstack([y - sub[lo].to_numpy(), sub[hi].to_numpy() - y])
            ax.errorbar(
                x,
                y,
                yerr=yerr,
                fmt=marker,
                color=COLOR[ft],
                ms=6,
                lw=1.3,
                capsize=3,
                label=("WUI" if ft == "Urban-edge" else ft),
            )
            ax.plot(x, y, "-", color=COLOR[ft], lw=1.0, alpha=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels(decades, rotation=20)
        ax.set_title(title, fontsize=11)
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        if key in ("weibull_k", "lambda_exp"):
            ax.text(
                0.02,
                0.02,
                "↓ heavier tails",
                transform=ax.transAxes,
                fontsize=8,
                color="0.4",
            )
        else:
            ax.text(
                0.02,
                0.02,
                "↑ heavier tails",
                transform=ax.transAxes,
                fontsize=8,
                color="0.4",
            )
    axes[0].legend(frameon=False, fontsize=9)
    axes[0].set_ylabel("Parameter value")
    fig.suptitle("CalFire decade trends under alternative fits", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_FIG / "Fig_size_dist_alt_decades.png", dpi=300, bbox_inches="tight", facecolor="w")
    fig.savefig(OUT_PRC / "Fig_size_dist_alt_decades.png", dpi=300, bbox_inches="tight", facecolor="w")
    print("Saved decade figure")

    # ── AIC Δ summary for reviewer ───────────────────────────────────────────
    print("\n=== AIC (lower better) ===")
    print(
        summary[
            [
                "dataset",
                "FireType",
                "aic_powerlaw",
                "aic_truncexpon",
                "aic_weibull",
                "best_aic",
            ]
        ].to_string(index=False)
    )
    print("\n=== Extreme P(X>1000) ===")
    print(
        extremes[extremes["threshold_km2"] == 1000][
            ["dataset", "FireType", "empirical", "powerlaw_mle", "truncexpon", "weibull"]
        ].to_string(index=False)
    )

    plot_mle_orig_figure(data)
    plot_ols_vs_mle_calfire_mtbs(data)


def plot_ols_vs_mle_calfire_mtbs(data=None):
    """2×2 OLS panels: CalFire | MTBS / Atlas | FIRED, with β and R²."""
    from matplotlib.lines import Line2D
    from matplotlib.ticker import NullLocator

    if data is None:
        data = load_groups()
    layout = [["CalFire", "MTBS"], ["Atlas", "FIRED"]]
    # Paper rewrite uses A0 = 1 km²; MTBS completeness threshold is ~4 km²
    xmin_fig = {"CalFire": 1.0, "MTBS": 4.0, "Atlas": 1.0, "FIRED": 1.0}

    def _fmt(beta, r2):
        return f"β={beta:.2f}, R$^2$={r2:.2f}"

    def _annotate(ax, lines):
        # columns: model | WUI | Wildland
        y0, dy = 0.02, 0.050
        x_lab, x_wui, x_wild = 0.02, 0.26, 0.62
        for i, (lab, su, sw) in enumerate(lines):
            yy = y0 + (len(lines) - 1 - i) * dy
            ax.text(
                x_lab, yy, f"{lab}", transform=ax.transAxes, ha="left", va="bottom",
                fontsize=7.5, color="0.2", zorder=11,
            )
            ax.text(
                x_wui, yy, su, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=7.5, color=COLOR["Urban-edge"], zorder=11,
            )
            ax.text(
                x_wild, yy, sw, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=7.5, color=COLOR["Wildland"], zorder=11,
            )

    def _style(ax, title):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-1, 1e5)
        ax.set_ylim(1e-7, 3e0)
        ax.set_xticks([1e-1, 1e1, 1e3, 1e5])
        ax.set_yticks([1e-7, 1e-5, 1e-3, 1e-1])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_minor_locator(NullLocator())
        ax.tick_params(which="major", direction="out", labelsize=8, width=0.7, length=3.5)
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        for sp in ["bottom", "left"]:
            ax.spines[sp].set_linewidth(0.7)
        ax.set_title(title, fontsize=8.5, fontweight="semibold", loc="left", pad=3)
        ax.set_xlabel("Fire size (km$^2$)", fontsize=8.5)
        ax.set_ylabel("Probability density", fontsize=8.5)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8), facecolor="w")
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.08, top=0.96, wspace=0.28, hspace=0.32)

    for i, row in enumerate(layout):
        for j, ds in enumerate(row):
            ax = axes[i, j]
            d = data[ds]
            xmin = xmin_fig[ds]
            params = {}
            for ft in ["Urban-edge", "Wildland"]:
                s = np.asarray(d[ft], float)
                s = s[np.isfinite(s) & (s >= xmin)]
                if len(s) < 20:
                    continue
                x, y = log_binned_pdf(s)
                c = COLOR[ft]
                ax.loglog(
                    x, y, "o", color=c, ms=3.2, zorder=2,
                    markerfacecolor="none", markeredgecolor=c, markeredgewidth=0.8,
                )
                fits = fit_all_for_group(s, xmin)
                ols, ln_ols, wb = fits["powerlaw_ols"], fits["lognormal_ols"], fits["weibull_ols"]
                xfit = np.logspace(np.log10(xmin), np.log10(s.max()), 400)
                p = {}
                if ols:
                    ax.loglog(
                        xfit, pdf_powerlaw_ols(xfit, ols["beta"], ols["intercept"]),
                        "-", color=c, lw=1.3, alpha=0.9, zorder=4,
                    )
                    p["pl_beta"], p["pl_r2"] = ols["beta"], ols["r2"]
                if ln_ols:
                    ax.loglog(
                        xfit,
                        pdf_lognormal_quad(xfit, ln_ols["ln_a"], ln_ols["beta"], ln_ols["psi"]),
                        ls=LS["lognormal_ols"], color=c, lw=1.35, alpha=0.9, zorder=5,
                    )
                    p["ln_beta"], p["ln_r2"] = -ln_ols["beta"], ln_ols["r2"]
                if wb:
                    ax.loglog(
                        xfit, pdf_weibull(xfit, wb["alpha"], wb["beta"], xmin),
                        ":", color=c, lw=1.2, alpha=0.85, zorder=3,
                    )
                    p["wb_beta"], p["wb_r2"] = -wb["beta"], wb["r2"]
                params[ft] = p
            _style(ax, ds)
            u, w = params["Urban-edge"], params["Wildland"]
            _annotate(
                ax,
                [
                    ("Power law", _fmt(u["pl_beta"], u["pl_r2"]), _fmt(w["pl_beta"], w["pl_r2"])),
                    ("Lognormal", _fmt(u["ln_beta"], u["ln_r2"]), _fmt(w["ln_beta"], w["ln_r2"])),
                    ("Weibull", _fmt(u["wb_beta"], u["wb_r2"]), _fmt(w["wb_beta"], w["wb_r2"])),
                ],
            )

    # Fitting forms → CalFire (upper right); WUI/Wildland → MTBS (upper right)
    axes[0, 0].legend(
        handles=[
            Line2D([0], [0], color="k", ls="-", lw=1.3, label=r"Power law: $y\sim x^{\beta}$"),
            Line2D(
                [0], [0], color="k", ls=LS["lognormal_ols"], lw=1.35,
                label=r"Lognormal: $y\sim x^{\beta}e^{-\psi(\ln x)^{2}}$",
            ),
            Line2D(
                [0], [0], color="k", ls=":", lw=1.2,
                label=r"Weibull: $y\sim x^{\beta}e^{-\frac{\alpha}{1+\beta}(x^{1+\beta}-x_{0}^{1+\beta})}$",
            ),
        ],
        loc="upper right", frameon=False, fontsize=8,
        bbox_to_anchor=(1.0, 1.08),
        handlelength=1.8, labelspacing=0.3, borderpad=0.2,
    )
    axes[0, 1].legend(
        handles=[
            Line2D(
                [0], [0], marker="o", color=COLOR["Urban-edge"], ls="-", lw=1.2, ms=4.5,
                markerfacecolor="none", markeredgecolor=COLOR["Urban-edge"], label="WUI",
            ),
            Line2D(
                [0], [0], marker="o", color=COLOR["Wildland"], ls="-", lw=1.2, ms=4.5,
                markerfacecolor="none", markeredgecolor=COLOR["Wildland"], label="Wildland",
            ),
        ],
        loc="upper right", frameon=False, fontsize=7.5,
        bbox_to_anchor=(1.0, 1.08),
        handlelength=1.8, labelspacing=0.25, borderpad=0.2,
    )
    for pth in [
        OUT_FIG / "FigS2.png",
        OUT_FIG / "FigS2.pdf",
        OUT_FIG / "Fig_size_dist_OLS_vs_MLE.png",
        OUT_PRC / "Fig_size_dist_OLS_vs_MLE.png",
        OUT_FIG / "Fig_size_dist_OLS_4src.png",
        OUT_FIG / "Fig_size_dist_OLS_4src.pdf",
        OUT_PRC / "Fig_size_dist_OLS_4src.png",
    ]:
        fig.savefig(pth, dpi=300, bbox_inches="tight", facecolor="w")
    print("Saved OLS 4-source figure")
    plt.close(fig)


def plot_powerlaw_vs_truncexpon(data=None):
    """
    2×2 figure matching Fig_size_dist_OLS_4src layout, but only:
      1) Power law (OLS on log-binned PDF)
      2) Trunc. exp. = power law × exp cutoff  p(x)∝ x^β e^{−x/xc} (MLE)
    """
    from matplotlib.lines import Line2D
    from matplotlib.ticker import NullLocator

    if data is None:
        data = load_groups()
    layout = [["CalFire", "MTBS"], ["Atlas", "FIRED"]]
    xmin_fig = {"CalFire": 1.0, "MTBS": 4.0, "Atlas": 1.0, "FIRED": 1.0}

    def _fmt(beta, r2):
        return f"β={beta:.2f}, R$^2$={r2:.2f}"

    def _r2_log(x, y, yhat):
        m = np.isfinite(x) & np.isfinite(y) & np.isfinite(yhat) & (y > 0) & (yhat > 0)
        if int(m.sum()) < 3:
            return np.nan
        yt, yh = np.log10(y[m]), np.log10(yhat[m])
        ss_tot = np.sum((yt - yt.mean()) ** 2)
        if ss_tot <= 0:
            return np.nan
        return float(1.0 - np.sum((yt - yh) ** 2) / ss_tot)

    def _annotate(ax, lines):
        y0, dy = 0.02, 0.055
        x_lab, x_wui, x_wild = 0.02, 0.28, 0.64
        for i, (lab, su, sw) in enumerate(lines):
            yy = y0 + (len(lines) - 1 - i) * dy
            ax.text(
                x_lab, yy, f"{lab}", transform=ax.transAxes, ha="left", va="bottom",
                fontsize=7.5, color="0.2", zorder=11,
            )
            ax.text(
                x_wui, yy, su, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=7.5, color=COLOR["Urban-edge"], zorder=11,
            )
            ax.text(
                x_wild, yy, sw, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=7.5, color=COLOR["Wildland"], zorder=11,
            )

    def _style(ax, title):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-1, 1e5)
        ax.set_ylim(1e-7, 3e0)
        ax.set_xticks([1e-1, 1e1, 1e3, 1e5])
        ax.set_yticks([1e-7, 1e-5, 1e-3, 1e-1])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_minor_locator(NullLocator())
        ax.tick_params(which="major", direction="out", labelsize=8, width=0.7, length=3.5)
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        for sp in ["bottom", "left"]:
            ax.spines[sp].set_linewidth(0.7)
        ax.set_title(title, fontsize=8.5, fontweight="semibold", loc="left", pad=3)
        ax.set_xlabel("Fire size (km$^2$)", fontsize=8.5)
        ax.set_ylabel("Probability density", fontsize=8.5)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8), facecolor="w")
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.08, top=0.90, wspace=0.28, hspace=0.32)

    for i, row in enumerate(layout):
        for j, ds in enumerate(row):
            ax = axes[i, j]
            d = data[ds]
            xmin = xmin_fig[ds]
            params = {}
            for ft in ["Urban-edge", "Wildland"]:
                s = np.asarray(d[ft], float)
                s = s[np.isfinite(s) & (s >= xmin)]
                if len(s) < 20:
                    continue
                x, y = log_binned_pdf(s)
                c = COLOR[ft]
                ax.loglog(
                    x, y, "o", color=c, ms=3.2, zorder=2,
                    markerfacecolor="none", markeredgecolor=c, markeredgewidth=0.8,
                )
                ols = fit_powerlaw_ols(s, size_min=xmin)
                cut = fit_pl_cutoff_mle(s, size_min=xmin)
                xfit = np.logspace(np.log10(xmin), np.log10(s.max()), 400)
                p = {}
                if ols:
                    ax.loglog(
                        xfit, pdf_powerlaw_ols(xfit, ols["beta"], ols["intercept"]),
                        "-", color=c, lw=1.3, alpha=0.9, zorder=4,
                    )
                    p["pl_beta"], p["pl_r2"] = ols["beta"], ols["r2"]
                if cut:
                    c_light = 0.45 * c + 0.55  # lighter tint of same hue
                    x_co = np.logspace(np.log10(xmin), 5.0, 500)
                    ax.loglog(
                        x_co, pdf_pl_cutoff(x_co, cut["beta"], cut["xc"], cut["xmin"]),
                        "-", color=c_light, lw=1.3, alpha=0.9, zorder=5,
                    )
                    yhat = pdf_pl_cutoff(x, cut["beta"], cut["xc"], cut["xmin"])
                    p["co_beta"] = cut["beta"]
                    p["co_xc"] = cut["xc"]
                    p["co_r2"] = _r2_log(x, y, yhat)
                params[ft] = p
            _style(ax, ds)
            u, w = params.get("Urban-edge", {}), params.get("Wildland", {})
            _annotate(
                ax,
                [
                    ("Power law", _fmt(u["pl_beta"], u["pl_r2"]), _fmt(w["pl_beta"], w["pl_r2"])),
                    ("Trunc. exp.", _fmt(u["co_beta"], u["co_r2"]), _fmt(w["co_beta"], w["co_r2"])),
                ],
            )

    handles = [
        Line2D(
            [0], [0], marker="o", color=COLOR["Urban-edge"], ls="-", lw=1.2, ms=4.5,
            markerfacecolor="none", markeredgecolor=COLOR["Urban-edge"], label="WUI",
        ),
        Line2D(
            [0], [0], marker="o", color=COLOR["Wildland"], ls="-", lw=1.2, ms=4.5,
            markerfacecolor="none", markeredgecolor=COLOR["Wildland"], label="Wildland",
        ),
        Line2D([0], [0], color="k", ls="-", lw=1.3, label="Power law"),
        Line2D([0], [0], color="0.65", ls="-", lw=1.3, label="Trunc. exp."),
    ]
    fig.legend(
        handles=handles, loc="upper center", ncol=4, frameon=False, fontsize=9.5,
        bbox_to_anchor=(0.5, 0.995),
    )
    for pth in [
        OUT_FIG / "Fig_size_dist_PL_vs_TruncExp.png",
        OUT_FIG / "Fig_size_dist_PL_vs_TruncExp.pdf",
        OUT_PRC / "Fig_size_dist_PL_vs_TruncExp.png",
        OUT_FIG / "FigS2_PL_TruncExp.png",
        OUT_FIG / "FigS2_PL_TruncExp.pdf",
    ]:
        fig.savefig(pth, dpi=300, bbox_inches="tight", facecolor="w")
    print("Saved Power law vs Trunc. exp. figure")
    plt.close(fig)


def plot_trunc_cutoff_only(data=None):
    """
    Separate 2×2 figure: only truncated-style tails.
      1) right-truncated Pareto  p(x)=C x^{−α} on [xmin, xmax]
      2) power law × exp cutoff  p(x)∝ x^β e^{−x/xc}
    """
    from matplotlib.lines import Line2D
    from matplotlib.ticker import NullLocator

    if data is None:
        data = load_groups()
    layout = [["CalFire", "MTBS"], ["Atlas", "FIRED"]]
    xmin_fig = {"CalFire": 1.0, "MTBS": 4.0, "Atlas": 1.0, "FIRED": 1.0}

    def _fmt_tp(beta, r2):
        return f"β={beta:.2f}, R$^2$={r2:.2f}"

    def _r2_log(x, y, yhat):
        m = np.isfinite(x) & np.isfinite(y) & np.isfinite(yhat) & (y > 0) & (yhat > 0)
        if int(m.sum()) < 3:
            return np.nan
        yt, yh = np.log10(y[m]), np.log10(yhat[m])
        ss_tot = np.sum((yt - yt.mean()) ** 2)
        if ss_tot <= 0:
            return np.nan
        return float(1.0 - np.sum((yt - yh) ** 2) / ss_tot)

    def _annotate(ax, lines):
        y0, dy = 1.03, 0.048
        x_lab, x_wui, x_wild = 0.0, 0.32, 0.72
        kw = dict(
            transform=ax.transAxes, ha="left", va="bottom",
            fontsize=6, clip_on=False, zorder=11,
        )
        for i, (lab, su, sw) in enumerate(lines):
            yy = y0 + (len(lines) - 1 - i) * dy
            ax.text(x_lab, yy, f"{lab}", color="0.2", **kw)
            ax.text(x_wui, yy, su, color=COLOR["Urban-edge"], **kw)
            ax.text(x_wild, yy, sw, color=COLOR["Wildland"], **kw)

    def _style(ax):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e-1, 1e5)
        ax.set_ylim(1e-7, 3e0)
        ax.set_xticks([1e-1, 1e1, 1e3, 1e5])
        ax.set_yticks([1e-7, 1e-5, 1e-3, 1e-1])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_minor_locator(NullLocator())
        ax.tick_params(which="major", direction="out", labelsize=8, width=0.45, length=2.8, pad=2)
        for sp in ax.spines.values():
            sp.set_linewidth(0.5)
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        ax.set_xlabel("Fire size (km$^2$)", fontsize=8.5, labelpad=2)
        ax.set_ylabel("Probability density", fontsize=8.5, labelpad=2)

    fig, axes = plt.subplots(2, 2, figsize=(6.4, 5.4), facecolor="w")
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.08, top=0.88, wspace=0.30, hspace=0.54)

    for i, row in enumerate(layout):
        for j, ds in enumerate(row):
            ax = axes[i, j]
            d = data[ds]
            xmin = xmin_fig[ds]
            params = {}
            for ft in ["Urban-edge", "Wildland"]:
                s = np.asarray(d[ft], float)
                s = s[np.isfinite(s) & (s >= xmin)]
                if len(s) < 20:
                    continue
                x, y = log_binned_pdf(s)
                c = COLOR[ft]
                ax.loglog(x, y, "o", color=c, ms=3.0, zorder=2)
                fits = fit_all_for_group(s, xmin)
                ln_ols, wb = fits["lognormal_ols"], fits["weibull_ols"]
                cut = fit_pl_cutoff_mle(s, size_min=xmin)
                xfit = np.logspace(np.log10(xmin), np.log10(s.max()), 400)
                p = {}
                if ln_ols:
                    ax.loglog(
                        xfit,
                        pdf_lognormal_quad(xfit, ln_ols["ln_a"], ln_ols["beta"], ln_ols["psi"]),
                        ls=LS["lognormal_ols"], color=c, lw=1.8, alpha=0.95, zorder=4,
                    )
                    p["ln_beta"], p["ln_r2"] = -ln_ols["beta"], ln_ols["r2"]
                if wb:
                    ax.loglog(
                        xfit, pdf_weibull(xfit, wb["alpha"], wb["beta"], xmin),
                        ":", color=c, lw=1.5, alpha=0.9, zorder=3,
                    )
                    p["wb_beta"], p["wb_r2"] = -wb["beta"], wb["r2"]
                if cut:
                    x_co = np.logspace(np.log10(xmin), 5.0, 500)
                    ax.loglog(
                        x_co, pdf_pl_cutoff(x_co, cut["beta"], cut["xc"], cut["xmin"]),
                        ls=LS["pl_cutoff"], color=c, lw=1.7, alpha=0.95, zorder=5,
                    )
                    yhat = pdf_pl_cutoff(x, cut["beta"], cut["xc"], cut["xmin"])
                    p["co_beta"] = cut["beta"]
                    p["co_xc"] = cut["xc"]
                    p["co_r2"] = _r2_log(x, y, yhat)
                    yc = float(
                        pdf_pl_cutoff(
                            np.array([cut["xc"]]), cut["beta"], cut["xc"], cut["xmin"]
                        )[0]
                    )
                    if np.isfinite(yc) and yc > 0:
                        y_lo = max(yc / 5.0, 1.2e-7)
                        y_hi = min(yc * 5.0, 2.0)
                        ax.plot(
                            [cut["xc"], cut["xc"]],
                            [y_lo, y_hi],
                            color=c,
                            ls="-.",
                            lw=1.3,
                            zorder=7,
                            solid_capstyle="round",
                        )
                        ax.text(
                            cut["xc"] * 1.12,
                            y_hi,
                            r"$x_c$",
                            color=c,
                            fontsize=7,
                            ha="left",
                            va="bottom",
                            zorder=8,
                        )
                params[ft] = p
            _style(ax)
            ax.text(
                0.03, 0.035, ds, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=10, fontweight="bold", zorder=11,
            )
            u, w = params["Urban-edge"], params["Wildland"]
            _annotate(
                ax,
                [
                    ("Trunc. exp.", _fmt_tp(u["co_beta"], u["co_r2"]), _fmt_tp(w["co_beta"], w["co_r2"])),
                    ("Lognormal", _fmt_tp(u["ln_beta"], u["ln_r2"]), _fmt_tp(w["ln_beta"], w["ln_r2"])),
                    ("Weibull", _fmt_tp(u["wb_beta"], u["wb_r2"]), _fmt_tp(w["wb_beta"], w["wb_r2"])),
                ],
            )

    h_data = [
        Line2D([0], [0], marker="o", color=COLOR["Urban-edge"], ls="none", ms=4.5, label="WUI"),
        Line2D([0], [0], marker="o", color=COLOR["Wildland"], ls="none", ms=4.5, label="Wildland"),
    ]
    h_fit = [
        Line2D([0], [0], color="k", ls="-", lw=1.5, label="Trunc. exp."),
        Line2D([0], [0], color="k", ls=LS["lognormal_ols"], lw=1.6, label="Lognormal"),
        Line2D([0], [0], color="k", ls=":", lw=1.4, label="Weibull"),
    ]
    ax0 = axes[0, 0]
    kw = dict(frameon=False, fontsize=6.8, borderaxespad=0.15, labelspacing=0.28, handletextpad=0.4)
    leg_data = ax0.legend(
        handles=h_data, loc="upper left", bbox_to_anchor=(0.36, 0.90),
        handlelength=1.1, **kw,
    )
    ax0.add_artist(leg_data)
    ax0.legend(
        handles=h_fit, loc="upper left", bbox_to_anchor=(0.62, 0.90),
        handlelength=1.8, **kw,
    )
    for pth in [
        OUT_FIG / "FigS2_trunc.png",
        OUT_FIG / "FigS2_trunc.pdf",
        OUT_FIG / "Fig_size_dist_trunc_cutoff.png",
        OUT_PRC / "Fig_size_dist_trunc_cutoff.png",
    ]:
        fig.savefig(pth, dpi=300, bbox_inches="tight", facecolor="w")
    print("Saved truncated-form figure")
    plt.close(fig)


def plot_mle_orig_figure(data=None):
    """
    Fit power law / lognormal / Weibull by MLE on original-space PDFs,
    then draw curves on log–log axes (Fig1B 2×2 layout).
    """
    from matplotlib.lines import Line2D

    if data is None:
        data = load_groups()
    datasets = ["CalFire", "MTBS", "Atlas", "FIRED"]
    rows = []
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.6), facecolor="w", sharex=True, sharey=True)
    for ax, ds in zip(axes.ravel(), datasets):
        d = data[ds]
        xmin = d["xmin"]
        for ft in ["Urban-edge", "Wildland"]:
            s = np.asarray(d[ft], float)
            s = s[np.isfinite(s) & (s >= xmin)]
            if len(s) < 20:
                continue
            x, y = log_binned_pdf(s)
            c = COLOR[ft]
            ax.loglog(x, y, "o", color=c, ms=2.5, zorder=2)

            pl = fit_powerlaw_mle_orig(s, size_min=xmin)
            ln = fit_lognormal_mle_orig(s, size_min=xmin)
            wb = fit_weibull_mle_orig(s, size_min=xmin)
            xfit = np.logspace(np.log10(xmin), np.log10(s.max()), 400)

            if pl:
                ax.loglog(
                    xfit,
                    pdf_powerlaw_orig(xfit, pl["alpha"], pl["beta"]),
                    "-",
                    color=c,
                    lw=1.6,
                    alpha=0.95,
                    zorder=4,
                )
            if ln:
                ax.loglog(
                    xfit,
                    pdf_lognormal_orig(xfit, ln["alpha"], ln["beta"], ln["psi"]),
                    ls=LS["lognormal_ols"],
                    color=c,
                    lw=1.8,
                    alpha=0.95,
                    zorder=5,
                )
            if wb:
                xwb = xfit[xfit >= wb["xmin"]]
                if len(xwb):
                    ax.loglog(
                        xwb,
                        pdf_weibull_orig(
                            xwb,
                            wb["alpha"],
                            wb["beta"],
                            Z=wb["Z"],
                            gamma=wb.get("gamma"),
                        ),
                        ":",
                        color=c,
                        lw=1.5,
                        alpha=0.9,
                        zorder=3,
                    )
            rows.append(
                {
                    "dataset": ds,
                    "FireType": ft,
                    "n": len(s),
                    "xmin": xmin,
                    "pl_beta": pl["beta"] if pl else np.nan,
                    "pl_ll": pl["ll"] if pl else np.nan,
                    "pl_aic": pl["aic"] if pl else np.nan,
                    "ln_alpha": ln["alpha"] if ln else np.nan,
                    "ln_beta": ln["beta"] if ln else np.nan,
                    "ln_psi": ln["psi"] if ln else np.nan,
                    "ln_ll": ln["ll"] if ln else np.nan,
                    "ln_aic": ln["aic"] if ln else np.nan,
                    "wb_alpha": wb["alpha"] if wb else np.nan,
                    "wb_beta": wb["beta"] if wb else np.nan,
                    "wb_ll": wb["ll"] if wb else np.nan,
                    "wb_aic": wb["aic"] if wb else np.nan,
                    "wb_xmin": wb["xmin"] if wb else np.nan,
                }
            )
            pl_s = f"PL β={pl['beta']:.3f}" if pl else "PL fail"
            ln_s = f"LN ψ={ln['psi']:.3f}" if ln else "LN fail"
            wb_s = f"WB β={wb['beta']:.3f}" if wb else "WB fail"
            print(f"MLE {ds:8s} {ft:12s}  {pl_s}  {ln_s}  {wb_s}")

        ax.set_xlim(1e-1, 1e5)
        ax.set_ylim(1e-7, 1e1)
        ax.set_xticks([1e-1, 1e1, 1e3, 1e5])
        ax.set_yticks([1e-7, 1e-4, 1e-1])
        ax.tick_params(which="both", direction="out")
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
        ax.text(
            0.05,
            0.05,
            ds,
            transform=ax.transAxes,
            fontweight="bold",
            ha="left",
            va="bottom",
            fontsize=11,
        )

    handles = [
        Line2D(
            [0], [0], marker="o", color=COLOR["Urban-edge"], ls="-", lw=1.4, ms=5, label="WUI"
        ),
        Line2D(
            [0], [0], marker="o", color=COLOR["Wildland"], ls="-", lw=1.4, ms=5, label="Wildland"
        ),
        Line2D([0], [0], color="k", ls="-", lw=1.6, label="Power law (MLE)"),
        Line2D([0], [0], color="k", ls=LS["lognormal_ols"], lw=1.8, label="Lognormal (MLE)"),
        Line2D([0], [0], color="k", ls=":", lw=1.5, label="Weibull (MLE)"),
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        ncol=5,
        frameon=False,
        fontsize=8.5,
        bbox_to_anchor=(0.5, 1.02),
    )
    fig.supxlabel("Fire size (km$^2$)", fontsize=11)
    fig.supylabel("Probability density", fontsize=11)
    fig.tight_layout()
    for pth in [
        OUT_FIG / "Fig_size_dist_alt_PDF_MLE.png",
        OUT_PRC / "Fig_size_dist_alt_PDF_MLE.png",
    ]:
        fig.savefig(pth, dpi=300, bbox_inches="tight", facecolor="w")
    pd.DataFrame(rows).to_csv(OUT_PRC / "size_dist_MLE_orig_fits.csv", index=False)
    print("Saved MLE original-equation figure + table")
    plt.close(fig)


if __name__ == "__main__":
    main()

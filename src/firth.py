"""Firth's penalized-likelihood logistic regression.

With 47 non-home-discharge events and small adjusted models (age acceleration plus a
few covariates), ordinary maximum likelihood is small-sample biased and unstable.
Firth's method adds the Jeffreys invariant prior 0.5*log|I(beta)| to the
log-likelihood, which removes the first-order bias and always yields finite
estimates. This module implements the penalized fit, profile-penalized-likelihood
confidence intervals, and penalized likelihood-ratio p-values, matching R's
``logistf`` (Heinze & Schemper 2002). It has no external dependency beyond
NumPy/SciPy so it runs in the offline analysis environment.

We report Firth estimates for every association odds ratio in the paper. Continuous
predictors are standardized to per-SD effects before fitting; binary predictors are
left as 0/1 and reported per unit.
"""
from __future__ import annotations

import numpy as np
from scipy.special import expit
from scipy.stats import chi2
from scipy.optimize import brentq

_QCHI = chi2.ppf(0.95, 1)  # 3.8415, threshold for a 95% profile-likelihood interval


def _penalized_loglik(X: np.ndarray, y: np.ndarray, beta: np.ndarray) -> float:
    """Firth penalized log-likelihood: l(beta) + 0.5*log|X'WX|."""
    eta = X @ beta
    # log-likelihood via logaddexp for numerical stability
    ll = float(np.sum(y * eta - np.logaddexp(0.0, eta)))
    p = expit(eta)
    w = np.clip(p * (1.0 - p), 1e-12, None)
    xtwx = X.T @ (X * w[:, None])
    sign, logdet = np.linalg.slogdet(xtwx)
    if sign <= 0:
        return -np.inf
    return ll + 0.5 * logdet


def _firth_solve(
    X: np.ndarray,
    y: np.ndarray,
    fixed: dict[int, float] | None = None,
    max_iter: int = 1000,
    tol: float = 1e-8,
):
    """Maximize the penalized log-likelihood by modified-score IRLS.

    ``fixed`` pins the given coefficient indices to fixed values (used for the
    profile likelihood); the remaining coefficients are updated with the
    Firth-modified score step restricted to the free block, with step-halving.
    Returns (beta, penalized_loglik, cov) where cov is the ML covariance
    (X'WX)^-1 at the solution (used only to bracket the profile search).
    """
    n, k = X.shape
    beta = np.zeros(k)
    if fixed:
        for j, v in fixed.items():
            beta[j] = v
    free = np.array([j for j in range(k) if not fixed or j not in fixed])

    pll = _penalized_loglik(X, y, beta)
    for _ in range(max_iter):
        eta = X @ beta
        p = expit(eta)
        w = np.clip(p * (1.0 - p), 1e-10, None)
        xtwx = X.T @ (X * w[:, None])
        try:
            xtwx_inv = np.linalg.inv(xtwx)
        except np.linalg.LinAlgError:
            xtwx_inv = np.linalg.pinv(xtwx)
        # hat-matrix diagonal of the weighted design
        h = w * np.einsum("ij,jk,ik->i", X, xtwx_inv, X)
        # Firth-modified score
        u = X.T @ (y - p + h * (0.5 - p))
        # Newton step on the free coordinates only
        hff = xtwx[np.ix_(free, free)]
        try:
            step_free = np.linalg.solve(hff, u[free])
        except np.linalg.LinAlgError:
            step_free = np.linalg.lstsq(hff, u[free], rcond=None)[0]

        # step-halving to guarantee ascent
        alpha = 1.0
        improved = False
        for _ in range(40):
            cand = beta.copy()
            cand[free] = beta[free] + alpha * step_free
            pll_cand = _penalized_loglik(X, y, cand)
            if pll_cand >= pll - 1e-12:
                beta, new_pll, improved = cand, pll_cand, True
                break
            alpha *= 0.5
        if not improved:
            break
        if abs(new_pll - pll) < tol and np.max(np.abs(alpha * step_free)) < tol:
            pll = new_pll
            break
        pll = new_pll

    eta = X @ beta
    w = np.clip(expit(eta) * (1.0 - expit(eta)), 1e-10, None)
    xtwx = X.T @ (X * w[:, None])
    try:
        cov = np.linalg.inv(xtwx)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(xtwx)
    return beta, pll, cov


def firth_logit(X: np.ndarray, y: np.ndarray):
    """Fit Firth logistic regression on design ``X`` (no intercept column) and
    binary ``y``. Returns a dict with per-coefficient estimate, 95% profile
    interval, and penalized likelihood-ratio p-value. Index 0 is the intercept.
    """
    y = np.asarray(y, dtype=float)
    Xc = np.column_stack([np.ones(len(y)), np.asarray(X, dtype=float)])
    k = Xc.shape[1]

    beta, pll_full, cov = _firth_solve(Xc, y)
    se = np.sqrt(np.clip(np.diag(cov), 1e-12, None))

    est, lo, hi, pval = (np.full(k, np.nan) for _ in range(4))
    est = beta.copy()
    for j in range(k):
        # penalized likelihood-ratio test for beta_j = 0
        _, pll_0, _ = _firth_solve(Xc, y, fixed={j: 0.0})
        lr = 2.0 * (pll_full - pll_0)
        pval[j] = float(chi2.sf(max(lr, 0.0), 1))

        # profile-penalized-likelihood 95% interval
        def g(c: float) -> float:
            _, pll_c, _ = _firth_solve(Xc, y, fixed={j: c})
            return 2.0 * (pll_full - pll_c) - _QCHI

        span = max(10.0 * se[j], 5.0)
        try:
            lo[j] = brentq(g, beta[j] - span, beta[j], maxiter=200, xtol=1e-6)
        except ValueError:
            lo[j] = beta[j] - span
        try:
            hi[j] = brentq(g, beta[j], beta[j] + span, maxiter=200, xtol=1e-6)
        except ValueError:
            hi[j] = beta[j] + span

    return {
        "beta": est,
        "or": np.exp(est),
        "ci_lo": np.exp(lo),
        "ci_hi": np.exp(hi),
        "p": pval,
        "loglik": pll_full,
        "se": se,
    }

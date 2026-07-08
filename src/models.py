"""Nested logistic models with leak-free repeated cross-validation.

The central correctness property: when ``leak_free`` is true (the reporting
default), the iterative imputer AND the standardizer are fit on the TRAINING
fold only and applied to the held-out fold. The legacy pipeline fit both on the
full dataset before splitting, leaking held-out information into every metric;
``leak_free=False`` reproduces that behavior for provenance only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from .features import model_specs

IMPUTER_MAXITER = 10  # small, well-conditioned predictor block


def _make_model(cfg: dict) -> LogisticRegression:
    # penalty defaults to L2 ridge; C controls shrinkage strength
    return LogisticRegression(C=cfg.get("l2_C", 1.0), max_iter=5000)


def preprocess_fold(Xtr, Xte, leak_free: bool, seed: int):
    """Fit imputation + standardization on the TRAIN fold only, apply to both.

    This is the anti-leakage contract in one place: when ``leak_free`` is true the
    iterative imputer is fit on ``Xtr`` and the scaler on the imputed ``Xtr``, so
    nothing about ``Xte`` can influence the fitted transformers. Returns
    (Xtr_t, Xte_t, scaler). In legacy mode imputation has already happened
    globally upstream, so only in-fold scaling remains.
    """
    if leak_free:
        imp = IterativeImputer(max_iter=IMPUTER_MAXITER, random_state=seed).fit(Xtr)
        Xtr, Xte = imp.transform(Xtr), imp.transform(Xte)
    scaler = StandardScaler().fit(Xtr)
    return scaler.transform(Xtr), scaler.transform(Xte), scaler


def cv_oof_predictions(
    df: pd.DataFrame, predictors: list[str], cfg: dict
) -> dict:
    """Repeated stratified K-fold out-of-fold predictions for one model.

    Returns a dict with the repeat-averaged OOF probability vector ``p``, the
    per-repeat AUC distribution, and per-repeat Brier scores. All preprocessing
    is fit inside folds when ``cfg['leak_free']`` is true.
    """
    outcome = cfg["outcome"]
    y = df[outcome].to_numpy(dtype=int)
    X = df[predictors].astype(float).to_numpy()
    n = len(y)
    leak_free = cfg["leak_free"]
    n_repeats, n_folds = cfg["n_repeats"], cfg["n_folds"]

    # Legacy (leaky) path: impute the full matrix once, before any splitting.
    if not leak_free:
        X = IterativeImputer(max_iter=20, random_state=0).fit_transform(X)

    seeder = np.random.default_rng(cfg["seed"])
    repeat_seeds = seeder.integers(0, 2**31 - 1, size=n_repeats)

    oof = np.full((n_repeats, n), np.nan)
    per_repeat_auc = np.empty(n_repeats)
    per_repeat_brier = np.empty(n_repeats)

    for r in range(n_repeats):
        rs = int(repeat_seeds[r])
        skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=rs)
        oof_r = np.empty(n)
        for k, (tr, te) in enumerate(skf.split(X, y)):
            Xtr_t, Xte_t, _ = preprocess_fold(X[tr], X[te], leak_free, rs + k)
            mod = _make_model(cfg).fit(Xtr_t, y[tr])
            pte = mod.predict_proba(Xte_t)[:, 1]
            oof[r, te] = pte
            oof_r[te] = pte
        per_repeat_auc[r] = roc_auc_score(y, oof_r)
        per_repeat_brier[r] = brier_score_loss(y, oof_r)

    p = np.nanmean(oof, axis=0)
    return {
        "p": p,
        "per_repeat_auc": per_repeat_auc,
        "per_repeat_brier": per_repeat_brier,
        "y": y,
    }


def run_nested_models(df: pd.DataFrame, cfg: dict) -> dict:
    """Fit all nested models; return per-model OOF results keyed by model name."""
    specs = model_specs()
    return {name: cv_oof_predictions(df, preds, cfg) for name, preds in specs.items()}


def _bootstrap_auc_ci(y, p, n_boot, seed):
    """Patient-level bootstrap 95% CI for AUC (sampling uncertainty)."""
    rng = np.random.default_rng(seed)
    n = len(y)
    aucs = []
    while len(aucs) < n_boot:
        bi = rng.integers(0, n, n)
        if len(np.unique(y[bi])) < 2:
            continue
        aucs.append(roc_auc_score(y[bi], p[bi]))
    return float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))


def model_metrics_table(results: dict, cfg: dict) -> pd.DataFrame:
    """Discrimination summary per model.

    Point estimate = AUC of the repeat-averaged out-of-fold predictions. The 95%
    CI is a patient-level bootstrap (sampling uncertainty). ``AUC_cv_sd`` reports
    the Monte-Carlo spread across CV repeats (estimator stability), a separate
    quantity from the sampling CI.
    """
    rows = []
    for name, res in results.items():
        y, p = res["y"], res["p"]
        aucs = res["per_repeat_auc"]
        lo, hi = _bootstrap_auc_ci(y, p, cfg["n_boot"], cfg["seed"] + 5)
        rows.append({
            "model": name,
            "n": len(y),
            "events": int(y.sum()),
            "AUC": float(roc_auc_score(y, p)),
            "AUC_lo": lo,
            "AUC_hi": hi,
            "AUC_cv_mean": float(aucs.mean()),
            "AUC_cv_sd": float(aucs.std()),
            "Brier": float(res["per_repeat_brier"].mean()),
        })
    return pd.DataFrame(rows)


def fit_full_logit(df: pd.DataFrame, cfg: dict, model: str = "M2_multimuscle"):
    """Fit the final model on all data for coefficient reporting (statsmodels).

    Uses in-sample iterative imputation + standardization (this is a descriptive
    coefficient table, not a performance estimate). Returns per-SD ORs with CIs.
    """
    import statsmodels.api as sm

    predictors = model_specs()[model]
    X = df[predictors].astype(float).to_numpy()
    X = IterativeImputer(max_iter=20, random_state=0).fit_transform(X)
    X = StandardScaler().fit_transform(X)
    y = df[cfg["outcome"]].to_numpy(dtype=int)
    res = sm.Logit(y, sm.add_constant(X)).fit(disp=0)
    params = np.asarray(res.params)[1:]
    ci = np.asarray(res.conf_int())[1:]
    pv = np.asarray(res.pvalues)[1:]
    return pd.DataFrame({
        "feature": predictors,
        "beta": params,
        "OR_per_SD": np.exp(params),
        "ci_lo": np.exp(ci[:, 0]),
        "ci_hi": np.exp(ci[:, 1]),
        "p": pv,
    })

"""Etiologic association analysis (Firth penalized logistic regression).

This module produces every *association* (odds-ratio) result in the paper. It
replaces the earlier unpenalized maximum-likelihood estimates with Firth's
penalized likelihood (see ``firth.py``), which is the appropriate estimator at
events-per-variable ~4.6 and under the near-separation present in the sparse
signal-intensity predictors.

Design (pre-specified; see ``features.py``):

* PRIMARY etiologic estimand — the total effect of preoperative iliopsoas volume
  on non-home discharge — is estimated from a parsimonious model adjusting only
  for confounders on the study DAG (age, sex, ASA class, hypertension, diabetes).
  Iliopsoas volume is fully observed, so this model is fit on complete cases with
  no imputation.
* The saturated multi-muscle model is retained only as a collinearity/suppression
  DIAGNOSTIC (Figure 1) and as a sensitivity analysis, not as the primary estimate.
* Continuous predictors are standardized to per-SD effects; binary predictors are
  left as 0/1 and reported per unit.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer

from .features import (BINARY, PRIMARY_EXPOSURE, model_specs, primary_spec,
                       sensitivity_specs)
from .firth import _firth_solve, firth_logit


def _design(df: pd.DataFrame, features: list[str], outcome: str, *,
            impute: bool = False, seed: int = 0):
    """Build a Firth design matrix for ``features``.

    Continuous predictors are standardized to per-SD; binary predictors are kept
    as 0/1. With ``impute=False`` the model is fit on complete cases (rows with
    any missing feature dropped); with ``impute=True`` missing values are filled
    by a single iterative imputation (disclosed as such in the manuscript).
    Returns (X, y, n).
    """
    cols = list(features)
    sub = df[cols + [outcome]].copy()
    y_all = sub[outcome].to_numpy(dtype=int)
    X_all = sub[cols].astype(float).to_numpy()

    if impute:
        X_all = IterativeImputer(max_iter=20, random_state=seed).fit_transform(X_all)
        mask = np.ones(len(sub), dtype=bool)
    else:
        mask = ~np.isnan(X_all).any(axis=1) & ~pd.isna(sub[outcome]).to_numpy()

    X = X_all[mask]
    y = y_all[mask]
    # standardize continuous columns only
    X = X.copy()
    for j, f in enumerate(cols):
        if f not in BINARY:
            sd = X[:, j].std(ddof=0)
            if sd > 0:
                X[:, j] = (X[:, j] - X[:, j].mean()) / sd
    return X, y, int(mask.sum())


def _or_rows(features: list[str], fit: dict) -> list[dict]:
    """Extract per-predictor OR/CI/p (skipping the intercept at index 0)."""
    rows = []
    for i, f in enumerate(features, start=1):
        rows.append({
            "feature": f, "unit": "per unit" if f in BINARY else "per SD",
            "OR": float(fit["or"][i]), "ci_lo": float(fit["ci_lo"][i]),
            "ci_hi": float(fit["ci_hi"][i]), "p": float(fit["p"][i]),
        })
    return rows


def firth_coefficients(df: pd.DataFrame, cfg: dict,
                       model: str = "M2_multimuscle") -> pd.DataFrame:
    """Full-model Firth coefficient table (per-SD/per-unit ORs with 95% CI)."""
    features = model_specs()[model]
    X, y, n = _design(df, features, cfg["outcome"], impute=True, seed=cfg["seed"])
    fit = firth_logit(X, y)
    out = pd.DataFrame(_or_rows(features, fit))
    out.insert(1, "n", n)
    out.insert(1, "model", model)
    return out


def firth_collinearity(df: pd.DataFrame, cfg: dict,
                       model: str = "M2_multimuscle") -> pd.DataFrame:
    """Crude (univariable Firth) vs adjusted (full Firth) per-predictor ORs, with
    VIF and a crude-to-adjusted sign-flip flag. A sign flip paired with a null
    crude association marks a suppression/collinearity artifact.

    Schema matches the legacy ``collinearity.csv`` so the forest figures and the
    manuscript continue to trace to the same columns.
    """
    features = model_specs()[model]
    outcome = cfg["outcome"]

    # adjusted = full mutually-adjusted Firth model (single imputation, disclosed)
    Xadj, yadj, _ = _design(df, features, outcome, impute=True, seed=cfg["seed"])
    adj = firth_logit(Xadj, yadj)
    vif = np.diag(np.linalg.inv(np.corrcoef(Xadj.T)))

    rows = []
    for i, f in enumerate(features):
        # crude = univariable Firth on complete cases for this feature
        Xu, yu, nu = _design(df, [f], outcome, impute=False)
        uni = firth_logit(Xu, yu)
        uni_or = float(uni["or"][1])
        adj_or = float(adj["or"][i + 1])
        rows.append({
            "feature": f, "n": nu,
            "uni_OR_per_SD": uni_or, "uni_lo": float(uni["ci_lo"][1]),
            "uni_hi": float(uni["ci_hi"][1]), "uni_p": float(uni["p"][1]),
            "adj_OR_per_SD": adj_or, "adj_lo": float(adj["ci_lo"][i + 1]),
            "adj_hi": float(adj["ci_hi"][i + 1]), "adj_p": float(adj["p"][i + 1]),
            "VIF": float(vif[i]),
            "sign_flip": bool((uni_or - 1) * (adj_or - 1) < 0),
        })
    return pd.DataFrame(rows)


def primary_model(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Primary confounder-adjusted etiologic model (Firth, complete cases).

    discharge ~ iliopsoas volume + age + sex + ASA + hypertension + diabetes.
    """
    features = primary_spec()
    X, y, n = _design(df, features, cfg["outcome"], impute=False)
    fit = firth_logit(X, y)
    out = pd.DataFrame(_or_rows(features, fit))
    out.insert(1, "n", n)
    return out


def iliopsoas_sensitivity(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Iliopsoas-volume OR across adjustment sets (stability of the primary effect)."""
    rows = []
    for name, features in sensitivity_specs().items():
        impute = name == "saturated_multimuscle"  # only this set touches missing muscles
        X, y, n = _design(df, features, cfg["outcome"], impute=impute, seed=cfg["seed"])
        fit = firth_logit(X, y)
        j = features.index(PRIMARY_EXPOSURE) + 1
        rows.append({
            "adjustment_set": name, "n_predictors": len(features), "n": n,
            "imputed": impute,
            "iliopsoas_vol_OR": float(fit["or"][j]),
            "ci_lo": float(fit["ci_lo"][j]), "ci_hi": float(fit["ci_hi"][j]),
            "p": float(fit["p"][j]),
        })
    return pd.DataFrame(rows)


def missingness(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Per-predictor missing-data summary (STROBE item 12c)."""
    feats = sorted(set(model_specs()["M2_multimuscle"]) | set(primary_spec()))
    n = len(df)
    rows = []
    for f in feats:
        n_obs = int(df[f].notna().sum())
        rows.append({
            "feature": f, "n_observed": n_obs, "n_missing": n - n_obs,
            "pct_missing": round(100 * (n - n_obs) / n, 1),
        })
    return pd.DataFrame(rows).sort_values("pct_missing", ascending=False)

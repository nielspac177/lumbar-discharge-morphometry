"""Internal validation: DeLong, reclassification (NRI/IDI) with CIs, optimism
correction, calibration, and coefficient stability.

All estimates operate on leak-free cross-validated out-of-fold predictions
except optimism correction, which runs its own bootstrap refit loop. NRI and IDI
are reported WITH bootstrap confidence intervals and framed cautiously (see
docs/adr and VALIDATION_REPORT) because they are known to be biased toward the
larger model (Pencina 2011/2012; Kerr 2014; Hilden 2014).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats as sps
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from .features import model_specs


# ======================================================================================
# DeLong test (fast algorithm, Sun & Xu 2014) — exact paired AUC comparison
# ======================================================================================

def _midrank(x: np.ndarray) -> np.ndarray:
    J = np.argsort(x)
    Z = x[J]
    N = len(x)
    T = np.zeros(N)
    i = 0
    while i < N:
        j = i
        while j < N and Z[j] == Z[i]:
            j += 1
        T[i:j] = 0.5 * (i + j - 1) + 1
        i = j
    T2 = np.empty(N)
    T2[J] = T
    return T2


def _fast_delong(preds_sorted_T: np.ndarray, m: int):
    n = preds_sorted_T.shape[1] - m
    pos = preds_sorted_T[:, :m]
    neg = preds_sorted_T[:, m:]
    k = preds_sorted_T.shape[0]
    tx = np.empty([k, m]); ty = np.empty([k, n]); tz = np.empty([k, m + n])
    for r in range(k):
        tx[r, :] = _midrank(pos[r, :])
        ty[r, :] = _midrank(neg[r, :])
        tz[r, :] = _midrank(preds_sorted_T[r, :])
    aucs = tz[:, :m].sum(axis=1) / m / n - (m + 1.0) / 2.0 / n
    v01 = (tz[:, :m] - tx) / n
    v10 = 1.0 - (tz[:, m:] - ty) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    delongcov = sx / m + sy / n
    return aucs, delongcov


def delong_test(y: np.ndarray, p1: np.ndarray, p2: np.ndarray) -> dict:
    """Two-sided DeLong test for the difference in AUC (model2 - model1)."""
    y = np.asarray(y)
    order = np.argsort(-y)              # positives (y=1) first
    m = int(y.sum())
    preds = np.vstack((p1, p2))[:, order]
    aucs, cov = _fast_delong(preds, m)
    var = cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]
    diff = aucs[1] - aucs[0]
    z = diff / np.sqrt(var) if var > 0 else 0.0
    pval = 2 * sps.norm.sf(abs(z))
    return {
        "auc1": float(aucs[0]), "auc2": float(aucs[1]),
        "delta_auc": float(diff), "se": float(np.sqrt(var)) if var > 0 else 0.0,
        "z": float(z), "p": float(pval),
    }


# ======================================================================================
# Reclassification — continuous NRI and IDI, with bootstrap CIs
# ======================================================================================

def nri_idi_point(y: np.ndarray, p_old: np.ndarray, p_new: np.ndarray) -> tuple[float, float]:
    ev, ne = y == 1, y == 0
    nri = ((p_new[ev] > p_old[ev]).mean() - (p_new[ev] < p_old[ev]).mean()) \
        + ((p_new[ne] < p_old[ne]).mean() - (p_new[ne] > p_old[ne]).mean())
    idi = (p_new[ev].mean() - p_old[ev].mean()) - (p_new[ne].mean() - p_old[ne].mean())
    return float(nri), float(idi)


def nri_idi_ci(y, p_old, p_new, n_boot: int, seed: int) -> dict:
    """Continuous NRI + IDI point estimates with stratified bootstrap 95% CIs."""
    y = np.asarray(y)
    nri, idi = nri_idi_point(y, p_old, p_new)
    rng = np.random.default_rng(seed)
    idx_ev = np.where(y == 1)[0]
    idx_ne = np.where(y == 0)[0]
    b_nri, b_idi = [], []
    for _ in range(n_boot):
        bi = np.concatenate([
            rng.choice(idx_ev, len(idx_ev), replace=True),
            rng.choice(idx_ne, len(idx_ne), replace=True),
        ])
        n, i = nri_idi_point(y[bi], p_old[bi], p_new[bi])
        b_nri.append(n); b_idi.append(i)
    return {
        "NRI": nri, "NRI_lo": float(np.percentile(b_nri, 2.5)),
        "NRI_hi": float(np.percentile(b_nri, 97.5)),
        "NRI_p": float(2 * min((np.array(b_nri) <= 0).mean(), (np.array(b_nri) >= 0).mean())),
        "IDI": idi, "IDI_lo": float(np.percentile(b_idi, 2.5)),
        "IDI_hi": float(np.percentile(b_idi, 97.5)),
        "IDI_p": float(2 * min((np.array(b_idi) <= 0).mean(), (np.array(b_idi) >= 0).mean())),
    }


# ======================================================================================
# Calibration — Cox recalibration slope & intercept, plus binned plot data
# ======================================================================================

def calibration(y: np.ndarray, p: np.ndarray, n_bins: int = 5) -> dict:
    """Calibration slope/intercept (logistic recalibration) + binned observations."""
    import statsmodels.api as sm

    eps = 1e-6
    logit_p = np.log(np.clip(p, eps, 1 - eps) / (1 - np.clip(p, eps, 1 - eps)))
    # Slope: logit(y) ~ a + b*logit(p)
    res = sm.Logit(y, sm.add_constant(logit_p)).fit(disp=0)
    intercept, slope = float(res.params[0]), float(res.params[1])
    # Intercept-in-the-large: offset logit(p) with slope 1
    res0 = sm.Logit(y, np.ones_like(y), offset=logit_p).fit(disp=0)
    citl = float(res0.params[0])
    # Binned observed vs predicted (deciles/quintiles)
    q = pd.qcut(p, n_bins, labels=False, duplicates="drop")
    binned = pd.DataFrame({"y": y, "p": p, "bin": q}).groupby("bin").agg(
        mean_pred=("p", "mean"), obs_rate=("y", "mean"), n=("y", "size")
    ).reset_index()
    return {"slope": slope, "intercept": intercept, "citl": citl, "binned": binned}


# ======================================================================================
# Optimism correction (Harrell/Efron bootstrap) for AUC
# ======================================================================================

def _fit_prob(Xtr, ytr, Xall, cfg, seed):
    imp = IterativeImputer(max_iter=10, random_state=seed).fit(Xtr)
    sc = StandardScaler().fit(imp.transform(Xtr))
    mod = LogisticRegression(C=cfg.get("l2_C", 1.0),
                             max_iter=5000).fit(sc.transform(imp.transform(Xtr)), ytr)
    return mod.predict_proba(sc.transform(imp.transform(Xall)))[:, 1]


def optimism_corrected_auc(df, predictors, cfg, n_boot: int | None = None) -> dict:
    """Bootstrap optimism-corrected apparent AUC (Harrell).

    apparent = AUC of model fit on full data, scored on full data.
    optimism  = mean over bootstraps of (AUC on bootstrap sample) - (AUC on original).
    corrected = apparent - optimism.
    """
    n_boot = n_boot or cfg["n_boot"] // 4
    y = df[cfg["outcome"]].to_numpy(dtype=int)
    X = df[predictors].astype(float).to_numpy()
    n = len(y)
    apparent = roc_auc_score(y, _fit_prob(X, y, X, cfg, cfg["seed"]))
    rng = np.random.default_rng(cfg["seed"] + 99)
    optimisms = []
    for b in range(n_boot):
        bi = rng.integers(0, n, n)
        if len(np.unique(y[bi])) < 2:
            continue
        p_boot = _fit_prob(X[bi], y[bi], X[bi], cfg, cfg["seed"] + b)
        p_orig = _fit_prob(X[bi], y[bi], X, cfg, cfg["seed"] + b)
        optimisms.append(roc_auc_score(y[bi], p_boot) - roc_auc_score(y, p_orig))
    optimism = float(np.mean(optimisms))
    return {"apparent_auc": float(apparent), "optimism": optimism,
            "corrected_auc": float(apparent - optimism), "n_boot": len(optimisms)}


# ======================================================================================
# Coefficient stability — bootstrap the sign/magnitude of each per-SD OR
# ======================================================================================

def collinearity_diagnostics(df, cfg, model: str = "M2_multimuscle"):
    """Univariable vs adjusted per-SD OR, and VIF, for each predictor.

    A predictor whose sign flips between the univariable and adjusted models with a
    null univariable association is a suppression/collinearity artifact, not an
    independent effect. Essential for honestly interpreting the deep-paraspinal
    volume coefficient.
    """
    import statsmodels.api as sm

    predictors = model_specs()[model]
    y = df[cfg["outcome"]].to_numpy(dtype=int)
    X = df[predictors].astype(float).to_numpy()
    X = IterativeImputer(max_iter=20, random_state=0).fit_transform(X)
    Xs = StandardScaler().fit_transform(X)

    # adjusted ORs with CIs
    adj = sm.Logit(y, sm.add_constant(Xs)).fit(disp=0)
    adj_or = np.exp(np.asarray(adj.params)[1:])
    adj_ci = np.exp(np.asarray(adj.conf_int())[1:])
    adj_p = np.asarray(adj.pvalues)[1:]

    # VIF from the standardized design correlation matrix
    vif = np.diag(np.linalg.inv(np.corrcoef(Xs.T)))

    rows = []
    for i, f in enumerate(predictors):
        uni = sm.Logit(y, sm.add_constant(Xs[:, i])).fit(disp=0)
        uni_or = float(np.exp(uni.params[1]))
        uni_ci = np.exp(np.asarray(uni.conf_int())[1])
        uni_p = float(uni.pvalues[1])
        rows.append({
            "feature": f, "n": len(y),
            "uni_OR_per_SD": uni_or, "uni_lo": float(uni_ci[0]),
            "uni_hi": float(uni_ci[1]), "uni_p": uni_p,
            "adj_OR_per_SD": float(adj_or[i]), "adj_lo": float(adj_ci[i, 0]),
            "adj_hi": float(adj_ci[i, 1]), "adj_p": float(adj_p[i]),
            "VIF": float(vif[i]),
            "sign_flip": bool((uni_or - 1) * (adj_or[i] - 1) < 0),
        })
    return pd.DataFrame(rows)


def coefficient_stability(df, cfg, model: str = "M2_multimuscle", n_boot: int | None = None):
    """Bootstrap Firth per-SD/per-unit ORs; report median OR, CI, and sign-consistency.

    Uses the same Firth point estimator as the main association analysis so the
    stability check is consistent with the reported coefficients. Only the point
    estimate is needed per resample, so the (slow) profile intervals are skipped.
    """
    from .features import BINARY
    from .firth import _firth_solve

    n_boot = n_boot or cfg["n_boot"] // 4
    predictors = model_specs()[model]
    cont = [j for j, f in enumerate(predictors) if f not in BINARY]
    y = df[cfg["outcome"]].to_numpy(dtype=int)
    X0 = df[predictors].astype(float).to_numpy()
    n = len(y)
    rng = np.random.default_rng(cfg["seed"] + 7)
    boot_or = []
    for b in range(n_boot):
        bi = rng.integers(0, n, n)
        if len(np.unique(y[bi])) < 2:
            continue
        Xb = IterativeImputer(max_iter=10, random_state=cfg["seed"] + b).fit_transform(X0[bi])
        Xb = Xb.copy()
        for j in cont:  # standardize continuous predictors; leave binary as 0/1
            sd = Xb[:, j].std(ddof=0)
            if sd > 0:
                Xb[:, j] = (Xb[:, j] - Xb[:, j].mean()) / sd
        try:
            beta, _, _ = _firth_solve(np.column_stack([np.ones(n), Xb]), y[bi].astype(float))
            boot_or.append(np.exp(beta[1:]))
        except Exception:
            continue
    boot_or = np.vstack(boot_or)
    rows = []
    for i, f in enumerate(predictors):
        col = boot_or[:, i]
        rows.append({
            "feature": f,
            "OR_median": float(np.median(col)),
            "OR_lo": float(np.percentile(col, 2.5)),
            "OR_hi": float(np.percentile(col, 97.5)),
            "frac_protective": float((col < 1).mean()),
            "sign_consistency": float(max((col < 1).mean(), (col > 1).mean())),
        })
    return pd.DataFrame(rows)

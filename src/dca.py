"""Decision-curve analysis with bootstrap uncertainty bands.

Net benefit = TP/n - FP/n * pt/(1-pt), evaluated over a threshold sweep for each
model plus the treat-all and treat-none references. Predictions are the leak-free
cross-validated out-of-fold probabilities, so the curves reflect internally
validated clinical utility rather than apparent (optimistic) utility.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _net_benefit(y: np.ndarray, p: np.ndarray, pt: np.ndarray) -> np.ndarray:
    n = len(y)
    w = pt / (1 - pt)
    nb = np.empty_like(pt)
    for i, t in enumerate(pt):
        flag = p >= t
        tp = np.sum(flag & (y == 1)) / n
        fp = np.sum(flag & (y == 0)) / n
        nb[i] = tp - fp * w[i]
    return nb


def decision_curve(y, preds: dict[str, np.ndarray], cfg,
                   lo: float = 0.05, hi: float = 0.80, step: float = 0.01) -> pd.DataFrame:
    """Net-benefit curves for each model + treat-all/none, with bootstrap 95% bands."""
    y = np.asarray(y)
    pt = np.arange(lo, hi + 1e-9, step)
    prev = y.mean()
    out = pd.DataFrame({"threshold": pt})
    out["nb_treat_none"] = 0.0
    out["nb_treat_all"] = prev - (1 - prev) * (pt / (1 - pt))
    for name, p in preds.items():
        out[f"nb_{name}"] = _net_benefit(y, p, pt)

    # Bootstrap bands for each model curve
    rng = np.random.default_rng(cfg["seed"] + 21)
    n = len(y)
    for name, p in preds.items():
        boot = np.empty((cfg["n_boot"] // 2, len(pt)))
        b = 0
        while b < boot.shape[0]:
            bi = rng.integers(0, n, n)
            if len(np.unique(y[bi])) < 2:
                continue
            boot[b] = _net_benefit(y[bi], p[bi], pt)
            b += 1
        out[f"nb_{name}_lo"] = np.percentile(boot, 2.5, axis=0)
        out[f"nb_{name}_hi"] = np.percentile(boot, 97.5, axis=0)
    return out


def delta_net_benefit(y, p_ref, p_new, cfg, lo=0.10, hi=0.40, step=0.01) -> pd.DataFrame:
    """Incremental net benefit of p_new OVER a reference model p_ref, with paired
    bootstrap 95% CIs. This is the correct test of a marker's ADDED clinical value
    (vs the clinical model), not the comparison against treat-all/none."""
    y = np.asarray(y)
    pt = np.arange(lo, hi + 1e-9, step)
    delta = _net_benefit(y, p_new, pt) - _net_benefit(y, p_ref, pt)
    rng = np.random.default_rng(cfg["seed"] + 41)
    n = len(y)
    boot = np.empty((cfg["n_boot"], len(pt)))
    b = 0
    while b < boot.shape[0]:
        bi = rng.integers(0, n, n)
        if len(np.unique(y[bi])) < 2:
            continue
        boot[b] = _net_benefit(y[bi], p_new[bi], pt) - _net_benefit(y[bi], p_ref[bi], pt)
        b += 1
    return pd.DataFrame({
        "threshold": pt, "delta_nb": delta,
        "delta_lo": np.percentile(boot, 2.5, axis=0),
        "delta_hi": np.percentile(boot, 97.5, axis=0),
    })


def delta_nb_summary(delta: pd.DataFrame, cfg, y, p_ref, p_new,
                     lo=0.10, hi=0.40, step=0.01) -> dict:
    """Mean incremental net benefit across the clinical range + P(mean>0)."""
    mean_delta = float(delta["delta_nb"].mean())
    pt = np.arange(lo, hi + 1e-9, step)
    rng = np.random.default_rng(cfg["seed"] + 43)
    y = np.asarray(y); n = len(y); pos = tot = 0
    for _ in range(cfg["n_boot"]):
        bi = rng.integers(0, n, n)
        if len(np.unique(y[bi])) < 2:
            continue
        md = float(np.mean(_net_benefit(y[bi], p_new[bi], pt) - _net_benefit(y[bi], p_ref[bi], pt)))
        pos += md > 0; tot += 1
    return {"mean_delta_nb_10_40": mean_delta, "extra_tp_per_100": mean_delta * 100,
            "prob_delta_positive": pos / tot}


def dca_summary(dca: pd.DataFrame, model_col: str, ref_col: str = "nb_treat_all",
                clin_lo: float = 0.10, clin_hi: float = 0.40) -> dict:
    """Summarize the model curve: positive-benefit range, mean gain, patients/100."""
    t = dca["threshold"].to_numpy()
    m = dca[model_col].to_numpy()
    all_ = dca[ref_col].to_numpy()
    positive = (m > 0) & (m >= all_) & (m >= dca["nb_treat_none"].to_numpy())
    band = (t >= clin_lo) & (t <= clin_hi)
    mean_gain_vs_all = float(np.mean(m[band] - all_[band]))
    # Interpretation: net-benefit gain × 100 ≈ extra true positives per 100 at no
    # net cost in false positives (standardized net benefit).
    return {
        "positive_range_lo": float(t[positive].min()) if positive.any() else float("nan"),
        "positive_range_hi": float(t[positive].max()) if positive.any() else float("nan"),
        "mean_gain_vs_treat_all_10_40": mean_gain_vs_all,
        "extra_tp_per_100": mean_gain_vs_all * 100,
    }

"""Tests for the Firth estimator and the aging-clock age-acceleration residual.

These run on a small self-contained synthetic frame (no real data, no scikit-learn)
so they are fast and safe to run anywhere.
"""
import numpy as np
import pandas as pd

from src.aging_clock import CLOCK_SPECS, build_features, fit_clock
from src.firth import firth_logit


def _synthetic_cohort(n=160, seed=0):
    """Minimal frame with the raw columns build_features() expects."""
    rng = np.random.default_rng(seed)
    age = rng.uniform(30, 85, n)
    # tissues degrade with age plus noise
    vert_vol = rng.normal(60, 8, n)
    vert_int = rng.normal(300, 40, n)
    df = pd.DataFrame({
        "age_yrs": age,
        "sex": rng.choice(["M", "F"], n),
        "asa": rng.integers(1, 4, n).astype(float),
        "rehab": rng.binomial(1, 0.25, n),
        "iliopsoas__vol_norm_vert": rng.normal(1.0, 0.2, n) - 0.004 * (age - 55),
        "iliopsoas__int_mean_mean": rng.normal(250, 40, n) + 1.5 * (age - 55),
        "deep_back__vol_norm_vert": rng.normal(2.5, 0.4, n) - 0.006 * (age - 55),
        "deep_back__int_mean_mean": rng.normal(260, 45, n) + 2.0 * (age - 55),
        "vertebra__vol_LM_cm3": vert_vol,
        "vertebra__int_mean": vert_int,
        "disc__vol_LM_cm3": rng.normal(45, 7, n),
        "disc__int_mean": rng.normal(400, 60, n) - 2.5 * (age - 55),
        "cord__vol_LM_cm3": rng.normal(15, 3, n),
        "cord__int_mean": rng.normal(500, 70, n),
    })
    return df


def test_build_features_columns():
    df = _synthetic_cohort()
    feats = build_features(df)
    for col in CLOCK_SPECS["scanner_robust"]:
        assert col in feats.columns


def test_age_acceleration_orthogonal_to_age():
    """The age-acceleration residual must be uncorrelated with chronological age."""
    df = _synthetic_cohort()
    cfg = {"outcome": "rehab", "seed": 20260518}
    clock, meta = fit_clock(df, cfg)
    r = np.corrcoef(clock["aar"], clock["age_yrs"])[0, 1]
    assert abs(r) < 1e-6, f"AAR not orthogonal to age (r={r})"
    assert abs(meta["corr_aar_age"]) < 1e-6


def test_firth_finite_under_separation():
    """Firth stays finite where maximum likelihood diverges."""
    y = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 1])
    x = np.arange(10, dtype=float)[:, None]
    fit = firth_logit(x, y)
    assert np.all(np.isfinite(fit["beta"]))
    assert np.all(np.isfinite(fit["ci_lo"])) and np.all(np.isfinite(fit["ci_hi"]))

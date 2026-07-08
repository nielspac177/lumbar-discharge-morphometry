"""Tiny mocked frames — NO real patient data. Enough columns/rows for the pipeline
to run so CI can smoke-test without the gitignored analytic dataset."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.features import CLINICAL, MUSCLES_M2, imaging_features


@pytest.fixture
def mock_cohort() -> pd.DataFrame:
    rng = np.random.default_rng(0)
    n = 120
    sex = rng.choice(["M", "F"], n)
    surgery = rng.integers(0, 2, n)
    df = pd.DataFrame({
        "rehab": rng.integers(0, 2, n),
        "sex": sex,
        "surgery": surgery,
        # derived clinical predictors (load_cohort adds these from sex/surgery)
        "female": (sex == "F").astype(int),
        "fusion": surgery.astype(int),
        "age_yrs": rng.normal(64, 11, n),
        "asa": rng.integers(1, 4, n).astype(float),
        "num_level": rng.integers(1, 4, n).astype(float),
        "bmi": rng.normal(28, 5, n),
        "tot_or_min": rng.normal(180, 60, n),
    })
    for f in imaging_features(MUSCLES_M2):
        df[f] = rng.normal(1.0, 0.3, n)
    for m in MUSCLES_M2:
        df[f"{m}__vol_LM_cm3_mean"] = rng.normal(100, 25, n)
        df[f"{m}__quality_svLM_mean"] = rng.uniform(0.4, 0.9, n)
    # inject some missingness in gluteus (matches real coverage gap)
    miss = rng.random(n) < 0.3
    df.loc[miss, "gluteus_medius__vol_norm_vert"] = np.nan
    df.loc[miss, "gluteus_medius__int_mean_mean"] = np.nan
    return df


@pytest.fixture
def smoke_cfg() -> dict:
    return {"outcome": "rehab", "seed": 1, "n_repeats": 2, "n_folds": 5,
            "n_boot": 50, "leak_free": True, "l2_C": 1.0}

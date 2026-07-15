"""Load the frozen analytic dataset and derive the modeling frame.

The analytic dataset (``data/analytic_cohort.csv``) is built once from the raw
segmentation workbook by ``build_dataset.py`` (PHI-touching, off the default
reproduction path). This module is the entry point for every analysis: it reads
that frozen CSV and produces the tidy per-patient frame used by the models.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from .features import CLINICAL, MUSCLES_M2, imaging_features


def load_config(path: str = "config.yaml") -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def load_cohort(cfg: dict) -> pd.DataFrame:
    """Read the frozen analytic CSV and build the modeling frame.

    Restricts to patients with a known discharge outcome and non-missing
    iliopsoas (the always-segmented reference muscle), derives the clinical
    predictors, and coerces the imaging block to numeric. Missing imaging values
    are left as NaN here — imputation happens INSIDE CV folds (see models.py).
    """
    df = pd.read_csv(Path(cfg["data_path"]))
    outcome = cfg["outcome"]

    df = df[df[outcome].notna()].copy()
    df[outcome] = df[outcome].astype(int)

    # Clinical derivations
    df["female"] = (df["sex"] == "F").astype(int)
    df["fusion"] = (df["surgery"] == 1).astype(int)
    for c in ["age_yrs", "asa", "num_level", "htn", "diabetes"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Drop implausible ages (data-entry errors, e.g. age recorded as ~0)
    df = df[df["age_yrs"].between(18, 100)].copy()

    # Imaging features → numeric (NaNs preserved for in-fold imputation)
    for c in imaging_features(MUSCLES_M2):
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Analytic cohort = has iliopsoas (reference muscle always present)
    df = df[df["iliopsoas__vol_norm_vert"].notna()].copy().reset_index(drop=True)
    return df

"""Covariate specification — the single source of truth for the model design.

Everything about *which* covariates enter the model and *how* they were selected
lives here so the manuscript Methods, the covariate-selection figure (WP-B3
Panel C), and the fitted models cannot drift apart. Do not hard-code predictor
lists anywhere else — import them from this module.

Covariate-selection rule (pre-specified; TRIPOD item 7a):
------------------------------------------------------------------------------
With only 51 non-home-discharge events, data-driven variable selection
(stepwise, univariate screening) is indefensible: it inflates optimism and
destabilizes coefficients. We therefore PRE-SPECIFY a small covariate set on
clinical and anatomical grounds and never let the data choose predictors. The
only model-building step is the pre-planned NESTING of blocks:

    M0  clinical            age, female, ASA, #levels, fusion
    M1  M0 + iliopsoas       + iliopsoas normalized volume + mean signal intensity
    M2  M1 + multi-muscle    + deep-paraspinal & gluteus-medius volume + intensity

Events-per-variable (EPV) for the largest model (M2, 11 predictors) is
51 / 11 = 4.6, below the classic 10 rule. Mitigations: (a) no free selection;
(b) ridge (L2) penalization for shrinkage; (c) internal validation with
optimism correction and leak-free repeated CV; (d) a parsimonious
clinical + iliopsoas-only sensitivity model (M1) reported alongside M2.
"""
from __future__ import annotations

# --- Clinical block (pre-specified from discharge-prediction literature) --------------
CLINICAL: list[str] = ["age_yrs", "female", "asa", "num_level", "fusion"]

# --- Imaging blocks (pre-specified anatomically) --------------------------------------
# Volume is normalized to L4 vertebral body volume (body-size scaling); intensity is the
# bilateral mean T2 signal (a fatty-infiltration / muscle-quality proxy).
MUSCLES_M1: list[str] = ["iliopsoas"]
MUSCLES_M2: list[str] = ["iliopsoas", "deep_back", "gluteus_medius"]


def imaging_features(muscles: list[str]) -> list[str]:
    """Normalized-volume + mean-intensity feature names for the given muscle groups."""
    return (
        [f"{m}__vol_norm_vert" for m in muscles]
        + [f"{m}__int_mean_mean" for m in muscles]
    )


# --- Nested model design (the ONLY model-building step) -------------------------------
def model_specs() -> dict[str, list[str]]:
    """Return the ordered predictor list for each nested model."""
    return {
        "M0_clinical": CLINICAL,
        "M1_iliopsoas": CLINICAL + imaging_features(MUSCLES_M1),
        "M2_multimuscle": CLINICAL + imaging_features(MUSCLES_M2),
    }


# --- Human-readable labels for figures/tables -----------------------------------------
LABELS: dict[str, str] = {
    "age_yrs": "Age (per SD)",
    "female": "Female sex",
    "asa": "ASA class (per SD)",
    "num_level": "No. of operated levels (per SD)",
    "fusion": "Fusion vs decompression",
    "iliopsoas__vol_norm_vert": "Iliopsoas volume (per SD)",
    "deep_back__vol_norm_vert": "Deep paraspinal volume (per SD)",
    "gluteus_medius__vol_norm_vert": "Gluteus medius volume (per SD)",
    "iliopsoas__int_mean_mean": "Iliopsoas mean signal (per SD)",
    "deep_back__int_mean_mean": "Deep paraspinal mean signal (per SD)",
    "gluteus_medius__int_mean_mean": "Gluteus medius mean signal (per SD)",
}

GROUP: dict[str, str] = {f: ("Clinical" if f in CLINICAL else "Imaging") for f in LABELS}


def epv(n_events: int, model: str = "M2_multimuscle") -> float:
    """Events-per-variable for a given nested model."""
    return n_events / len(model_specs()[model])

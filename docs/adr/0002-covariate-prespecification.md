# ADR 0002 — Pre-specified covariates and events-per-variable

**Status:** accepted · **Date:** 2026-07-08

## Context
With 51 events, the events-per-variable (EPV) for the largest model (M2, 11
predictors) is **51/11 ≈ 4.6**, below the traditional EPV ≥ 10 guideline
(Peduzzi 1996). Data-driven variable selection (stepwise, univariate screening)
at this scale is known to inflate optimism and destabilize coefficient signs.

## Decision
1. **Pre-specify** all covariates on clinical/anatomical grounds; never let the
   data choose predictors. The predictor lists live in `src/features.py` as the
   single source of truth.
   - Clinical block: age, female sex, ASA class, number of operated levels,
     fusion vs decompression.
   - Imaging block: for each pre-selected muscle group, L4-normalized volume and
     bilateral-mean T2 signal intensity.
2. The **only** model-building step is the pre-planned nesting of blocks:
   `M0 clinical → M1 +iliopsoas → M2 +multi-muscle (deep paraspinal, gluteus medius)`.
3. Mitigate low EPV with (a) no free selection, (b) **ridge (L2) penalization**
   for shrinkage, (c) leak-free internal validation with **optimism correction**,
   and (d) reporting the parsimonious M1 (clinical + iliopsoas) alongside M2.

## Consequences
- Coefficients are interpretable and, empirically, **stable** under bootstrap
  (see `results/coef_stability.csv`; e.g. iliopsoas signal sign-consistency 0.99).
- Reviewers get an explicit, defensible answer to "how were covariates selected?"
  This decision is depicted in the covariate-selection methods figure.
- We accept modest overfitting (optimism ≈ 0.04 AUC) rather than the larger,
  selection-induced optimism a stepwise search would introduce.

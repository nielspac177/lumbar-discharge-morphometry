# ADR 0003 — Leak-free internal validation

**Status:** accepted · **Date:** 2026-07-08

## Context
The original exploratory scripts fit the iterative imputer and the standardizer
on the **full dataset before** cross-validation splitting (e.g.
`03_abstract1_model.py`, `10_new_analyses.py`, `13_jama_forest_table.py`). This
leaks held-out information into every fold and optimistically biases AUC, NRI,
IDI, and net benefit.

## Decision
- **All preprocessing is fit inside each training fold** and applied to the
  held-out fold (`models.cv_oof_predictions`, `leak_free=True`, the reporting
  default). A `leak_free=False` path is retained solely to reproduce the legacy
  numbers for provenance.
- Discrimination is estimated by **repeated stratified 10-fold CV (×100)**; the
  point estimate is the AUC of the repeat-averaged out-of-fold predictions with a
  **patient-level bootstrap 95% CI**. Monte-Carlo spread across repeats is
  reported separately as estimator stability.
- **Bootstrap optimism correction** (Harrell/Efron, `validation.optimism_corrected_auc`)
  is reported as an independent overfitting check.
- A unit test (`tests/test_no_leakage.py`) fails if preprocessing responds to a
  held-out-row outlier, guarding against regression to the leaky pattern.

## Consequences
- Reported performance is internally validated and defensible.
- Empirically the corrected numbers closely match the legacy ones here (the small
  pre-specified ridge model was fairly robust), but the correction is required to
  *demonstrate* that, not assume it.

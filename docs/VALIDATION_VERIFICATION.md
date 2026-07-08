# Validation Verification Memo — Discharge Morphometry Statistics

**Date:** 2026-07-08
**Scope:** Synthesis of three independent adversarial reviews of the discharge-disposition
morphometry analysis (`results/`, `src/`, `docs/manuscript.md`).
**Reviews synthesized:** (1) Independent recompute from CSV outputs, (2) code red-team for
leakage/correctness, (3) biostatistician referee on methods and framing.

---

## 1. Overall Verdict

**Numbers: CONFIRMED.** Every requested quantity was independently re-derived from the
committed result files (`results/roc_data.csv`, `results/dca.csv`) without importing project
code, and all match to full displayed precision.

**Code: LEAK-FREE and CORRECT.** In the reporting default (`leak_free=True`), imputation and
standardization are fit strictly on the train fold inside each CV fold; every reported
internal-validation metric is computed on out-of-fold predictions. The fast DeLong,
continuous NRI/IDI, and Harrell optimism correction are faithful, correctly-signed
implementations. **No hard correctness bug was found.**

**Manuscript framing: ISSUES FOUND (not code bugs).** The statistics are honest, but the
central *clinical-utility* claim is not yet supported by the analyses shown. The added value
of morphometry over the clinical model M0 is essentially undemonstrated on discrimination
(ΔAUC +0.0055, p=0.79) and reclassification (IDI NS), while the headline decision-curve
"gain" is measured against the wrong comparator (treat-all instead of M0).

**Bottom line:** Nothing blocks the code. The required work is manuscript revision, not
re-computation.

---

## 2. What Was Checked

| # | Check | Reviewer | Status |
|---|-------|----------|--------|
| 1 | AUC of each model (M0=0.7652, M1=0.7549, M2=0.7707) vs `model_metrics.csv` | Recompute | PASS (exact) |
| 2 | DeLong ΔAUC M0→M2 = +0.005475, SE=0.020806, z=−0.2631, p=0.7924 | Recompute | PASS (exact; bootstrap p≈0.788 cross-check) |
| 3 | Continuous NRI M0→M2 = +0.5279 | Recompute | PASS (exact) |
| 4 | IDI M0→M2 = +0.01991 | Recompute | PASS (exact) |
| 5 | Net benefit of M2 at threshold 0.20 = 0.10366 | Recompute | PASS (exact) |
| 6 | Cohort n=205, events=51, prevalence=0.2488 (51/205) | Recompute | PASS |
| 7 | Imputation + scaling fit on train fold only (`preprocess_fold`, `models.py:39-43`) | Red-team | PASS |
| 8 | No full-data preprocessing contaminates any reported metric | Red-team | PASS |
| 9 | Fast DeLong (midrank, covariance, variance of difference) = Sun & Xu 2014 | Red-team | PASS |
| 10 | NRI/IDI continuous definitions + outcome-stratified bootstrap CIs | Red-team | PASS |
| 11 | Harrell optimism: sign, fit/predict targets (`corrected = apparent − optimism`) | Red-team | PASS |
| 12 | Determinism: all RNG derived from master seed; no global `np.random`/datetime | Red-team | PASS |
| 13 | Leak-free unit test asserts the anti-leakage contract | Red-team | PASS |
| 14 | Negative-control leakage test | Red-team | CONCERN (coverage gap, not a defect) |
| 15 | EPV / overfitting adequacy (EPV 4.64; slope 0.81) | Referee | CONCERN |
| 16 | Utility claim vs flat AUC + NS IDI | Referee | FAIL (framing) |
| 17 | DCA interpretation (comparator + averaging) | Referee | FAIL (framing) |
| 18 | Calibration slope <1 propagated to DCA/NRI | Referee | CONCERN |
| 19 | Single-center / external validation | Referee | FAIL (framing) |
| 20 | Multiple comparisons | Referee | FAIL (framing) |
| 21 | Deep-paraspinal opposite-direction OR (1.79, p=0.026) | Referee | FAIL (framing) |
| 22 | Whether category-free NRI should be reported | Referee | CONCERN (framing) |

Checks 1–13 concern the integrity of the computation and pipeline — all pass. Checks 15–22
concern how the results are interpreted and presented in the manuscript.

---

## 3. Prioritized Action List

### REQUIRED FIXES — code bugs

**None.** Both the independent recompute and the code red-team returned `confirmed` with an
empty `required_fixes` list. The numbers are reproducible and the pipeline is leak-free. No
change to `src/` is required for correctness.

### RECOMMENDED — coverage hardening (optional, non-blocking)

- **R-C1.** Add a test that exercises the real `leak_free=False` legacy path (global
  imputation at `models.py:63-64`), which currently has no test. The existing negative
  control uses an inline leaky scaler and never touches the legacy code path.
- **R-C2.** Extend the positive leak-free test to assert on the imputer directly, not only
  the scaler statistics.

### RECOMMENDED — manuscript framing (these are the substantive items)

Priority order (highest first):

1. **Show DCA of M2 vs M0 directly**, not vs treat-all. Plot both net-benefit curves and the
   delta-net-benefit with bootstrap CIs across 10–40%. The headline
   `mean_gain_vs_treat_all_10_40 = 0.1135` ("~11 extra TP per 100", `dca_summary.json`) is a
   bar M0 already clears; it cannot be attributed to muscle features. If M2−M0 delta-NB
   overlaps zero, revise the utility claim.
2. **Reframe the utility conclusion as hypothesis-generating** pending external/temporal
   validation. Single-center n=205 with flat added discrimination cannot support a
   clinical-utility claim.
3. **Resolve the deep-paraspinal opposite-direction OR** (`deep_back__vol_norm_vert`
   OR 1.79, p=0.026): report its univariable OR and collinearity diagnostics (correlation
   matrix / VIF among muscle predictors). Do not use bootstrap sign-consistency as evidence
   the effect is real — a sign flip from null/protective univariable to 1.79 adjusted is a
   textbook suppression signature and should not get a biological narrative.
4. **Drop category-free NRI as evidence** (or move to supplement). Its "CI excludes zero"
   (+0.528, 0.214–0.828) is being used as support while ΔAUC is flat and IDI is NS — exactly
   the miscalibration-driven false positive the cited literature (Pencina 2012, Kerr 2014,
   Hilden 2014) warns about. The NRI/IDI divergence is itself a red flag.
5. **Add a multiple-comparisons statement** and stop presenting borderline unadjusted ORs as
   robust: deep-paraspinal p=.026, iliopsoas volume p=.037, age p=.042 are all fragile under
   any correction. Pre-specification (ADR 0002) mitigates selection but not family-wise error.
6. **Confirm in text** that DCA, NRI, IDI, and calibration are all computed on out-of-fold
   predictions (the code confirms they are), and either compute DCA/net benefit from
   recalibrated/shrunken probabilities given the 0.81 calibration slope, or add a sensitivity
   analysis.
7. **Reconcile the two internal-validation M2 AUCs** into one consistently reported number:
   CV mean 0.771 (`model_metrics.csv`) vs optimism-corrected 0.790
   (`optimism_corrected.csv`). Also report a formal sample-size/shrinkage adequacy check
   (pmsampsize or equivalent) and the uniform shrinkage factor for EPV 4.6.

---

## 4. Referee's Top Objections vs How the Paper Answers Them

| Objection | Paper's current answer | Adequate? |
|-----------|------------------------|-----------|
| Utility claim rests on DCA vs treat-all, the wrong comparator | Reports `mean_gain_vs_treat_all_10_40 = 0.1135`; never shows M2 vs M0 delta-NB | **No.** Must show M2-vs-M0 net benefit with CIs. Given ΔAUC p=0.79 and NS IDI on the same information, the added net benefit is very likely negligible. |
| Added value of morphometry is undemonstrated | ADR 0002/0003/0004 honestly report flat AUC and NS IDI | **Honest, but conclusion overreaches.** The discrimination story is credible; the utility story is not yet supported. |
| Deep-paraspinal OR points the wrong way (sarcopenia framing predicts protective) | VALIDATION_REPORT §6: "stable, sign-consistency 0.986, CI excludes 1" | **No.** Bootstrap resamples the same 205 patients; it shows within-sample stability, not that the effect is real vs a suppression artifact. |
| NRI used as supportive evidence despite flat AUC/NS IDI | Already demoted below DCA with bias caveats (ADR 0004) | **Partial.** Demotion is good practice, but "CI excludes zero" language still leans on NRI as support — contradictory. Drop it. |
| EPV 4.64 too low; residual overfitting (slope 0.81) | Pre-specified covariates, ridge shrinkage, bootstrap optimism ~0.044 (VALIDATION_REPORT §2) | **Partial.** Optimism is corrected, but no formal shrinkage-adequacy check, and two M2 AUCs (0.771 vs 0.790) are quoted inconsistently. |
| Single-center, no external validation | Repeated stratified 10-fold ×100 + bootstrap optimism | **Partial.** Corrects optimism, says nothing about transportability. Reframe as hypothesis-generating. |
| No multiple-comparisons handling | Pre-specification (ADR 0002) | **Partial.** Pre-specification addresses data-driven selection, not family-wise error across AUC/NRI/IDI/DCA/~11 ORs. |

---

## Provenance

Values above verified against: `results/model_metrics.csv`, `results/reclassification.csv`,
`results/calibration.csv`, `results/optimism_corrected.csv`, `results/dca_summary.json`,
`results/model_coefficients.csv`, `results/run_meta.json`. Code paths cited from `src/models.py`,
`src/validation.py`, `tests/test_no_leakage.py`.

---

## Resolution (2026-07-08) — all required fixes implemented

The adversarial review's required fixes were acted on before release:

1. **Correct comparator added.** `src/dca.delta_net_benefit` now reports the
   incremental net benefit of the multi-muscle model **over the clinical model**
   (`results/delta_net_benefit.csv`): mean +0.002 across 10–40% (95% CI spanning
   zero; bootstrap P[>0] = 0.59). The former "+0.114 vs treat-all" is retained only
   as context (`FigS_decision_curve_vs_treatall`).
2. **Deep-paraspinal artifact resolved.** `src/validation.collinearity_diagnostics`
   (`results/collinearity.csv`) shows the univariable deep-paraspinal OR is null
   (0.96, P=.82), flipping to 1.79 only on adjustment (VIF/collinearity). It is now
   labeled a suppression artifact and not interpreted; Figure 1 flags it (†).
3. **Reframed to an etiologic association study.** The manuscript no longer claims
   added predictive utility. Continuous NRI is demoted (not used as evidence); the
   utility conclusion is stated as negative/hypothesis-generating pending external
   validation.
4. **Iliopsoas signal hedged.** Its univariable association is null (VIF ≈ 14);
   reported as adjusted-only/hypothesis-generating. Iliopsoas **volume** (robust
   univariably and adjusted, VIF 2.2) is the primary finding.

No code bugs were found; the pipeline is leak-free and deterministic.

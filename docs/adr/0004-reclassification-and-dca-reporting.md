# ADR 0004 — NRI/IDI and decision-curve reporting stance

**Status:** accepted · **Date:** 2026-07-08

## Context
Reclassification indices (continuous NRI, IDI) are widely reported but are known
to be biased toward the larger model and can be positive by chance even when the
added marker is useless (Pencina 2011/2012; Kerr 2014; Hilden & Gerds 2014).
Decision-curve analysis / net benefit (Vickers & Elkin 2006) is a more robust
measure of clinical utility.

## Decision
- **Foreground decision-curve analysis (net benefit)** as the primary utility
  metric, computed on leak-free out-of-fold predictions with **bootstrap
  uncertainty bands** (`dca.decision_curve`).
- Report continuous **NRI and IDI with bootstrap 95% CIs** (`validation.nri_idi_ci`),
  explicitly note their known biases, and treat them as **secondary/supportive**,
  not as primary evidence. Where an interval crosses zero (e.g. IDI here), say so
  plainly.
- Compare nested models' AUC with the **exact DeLong test** (`validation.delong_test`,
  Sun & Xu 2014), not a bootstrap approximation.
- Report **calibration** (slope, intercept, calibration-in-the-large) alongside
  discrimination, per the calibration-hierarchy argument (Van Calster 2019).

## Consequences
- The paper's utility claim rests on net benefit + calibration, which survive
  methodological scrutiny, rather than on reclassification indices alone.
- Honest reporting of a non-significant IDI strengthens, not weakens, credibility.

## Amendment (2026-07-08, after adversarial review)
Added value must be judged **against the clinical model**, not against treat-all. The
incremental decision-curve analysis (`dca.delta_net_benefit`, M2 vs M0) shows a mean
incremental net benefit of ≈0 across the 10–40% range (95% CI spanning zero;
bootstrap P(>0) ≈ 0.59). Combined with the flat ΔAUC and non-significant IDI, muscle
morphometry does **not** improve prediction. The paper is consequently reframed as an
**etiologic association study**: reclassification/DCA are reported as secondary
evidence that the association does *not* translate into added predictive utility.
Continuous NRI, though positive, is not used as evidence — it is the known-biased
index and contradicts every robust metric.

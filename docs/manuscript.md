---
title: "Preoperative Iliopsoas Muscle Morphometry Is Associated With Non-Home Discharge After Lumbar Spine Surgery: A Segmentation-Based Cohort Study"
short_title: "Iliopsoas Morphometry and Non-Home Discharge"
authors: "Jimena Gonzalez-Salido, MD*; Niels Pacheco-Barrios, MD* (*co-first authors). Department of Neurosurgery, Beth Israel Deaconess Medical Center / Harvard Medical School, Boston, MA."
reporting: "STROBE (cohort); prediction sub-analysis per TRIPOD"
bibliography: refs.bib
---

# Structured Abstract

**Importance** Sarcopenia and muscle quality are biologically plausible determinants
of recovery after spine surgery, but whether preoperative muscle morphometry is
associated with non-home discharge — and whether it adds anything to routine
clinical risk factors — is unclear.

**Objective** To examine the association between preoperative, MRI-segmentation–
derived muscle morphometry and non-home discharge after lumbar spine surgery, and,
secondarily, to test whether muscle morphometry improves prediction beyond clinical
factors.

**Design, Setting, and Participants** Retrospective single-center cohort of adults
undergoing lumbar spine surgery with a preoperative lumbar MRI. Iliopsoas, deep
paraspinal, and gluteus medius volume and T2 signal intensity were segmented at
L3–L5 (3D Slicer/TotalSegmentator). The analytic cohort comprised 205 patients.

**Main Outcomes and Measures** Non-home discharge (inpatient rehabilitation or
skilled nursing facility vs home). Associations were estimated with multivariable
logistic regression adjusting for age, sex, ASA class, number of operated levels,
and fusion, reporting per–standard-deviation (SD) odds ratios (ORs). The secondary
prediction analysis compared nested models with leak-free repeated cross-validation,
the DeLong test, reclassification indices, and decision-curve analysis.

**Results** Of 205 patients (median age 64.5 years; 48% female), 51 (24.9%) had a
non-home discharge. Lower preoperative iliopsoas volume was robustly associated with
higher odds of non-home discharge, both univariably (OR per SD, 0.57 [95% CI,
0.38–0.85]) and in the adjusted model (OR per SD, 0.52 [95% CI, 0.28–0.96]), with low
collinearity. Lower iliopsoas T2 signal intensity was associated in the adjusted
model (OR per SD, 0.11 [95% CI, 0.02–0.68]) but not univariably (OR, 0.82; P = .32)
and had high collinearity with the other signal measures, and is therefore
interpreted cautiously. Older age (OR, 1.66 [95% CI, 1.02–2.71]) and higher ASA class
(OR, 2.09 [95% CI, 1.31–3.35]) were also associated. The adjusted deep-paraspinal
volume coefficient was in the opposite direction (OR, 1.79) but had a null univariable
association (OR, 0.96; P = .82) and is attributable to collinearity with iliopsoas
volume rather than an independent effect. In the secondary prediction analysis,
adding muscle morphometry did **not** improve discrimination (multi-muscle vs
clinical area under the curve, 0.77 vs 0.77; ΔAUC +0.006, P = .79), reclassification
(integrated discrimination improvement +0.020; 95% CI, −0.017 to 0.054), or net
benefit over the clinical model (mean incremental net benefit across a 10%–40%
threshold range, +0.002; 95% CI, −0.018 to 0.020).

**Conclusions and Relevance** Lower preoperative iliopsoas volume was independently
associated with non-home discharge after lumbar spine surgery, consistent with a role
for sarcopenia, but muscle morphometry did not improve prediction beyond routine
clinical factors. These findings support muscle morphometry as a biologically
informative correlate rather than a
ready-to-use predictor, and illustrate that statistically significant associations
need not translate into added predictive utility. External validation is warranted.

---

# Introduction

Non-home discharge — transfer to inpatient rehabilitation or a skilled nursing
facility rather than home — is a frequent and consequential outcome after lumbar
spine surgery. Sarcopenia and paraspinal and psoas muscle quality are established,
biologically plausible determinants of adverse outcomes after surgery in general and
spine surgery in particular [@cruzjentoft2010ewgsop; @cruzjentoft2019ewgsop2;
@prado2008sarcopenia; @zakaria2015morphometrics; @bokshan2016sarcopenia;
@moskven2018sarcopenia]. Muscle degeneration is quantifiable from routinely acquired
MRI, both as volume and as fatty-infiltration–related signal intensity
[@fortin2014paraspinal; @takayama2016paraspinal]. Whether such measurements are
associated with the practical, resource-intensive outcome of non-home discharge —
and whether they add anything to the demographic and procedural variables that
existing tools already use [@oldmeadow2003rapt; @karhade2018discharge;
@ogink2019discharge; @stopa2019discharge] — has not been established.

We therefore examined, in a single-center cohort, the association between
preoperative multi-muscle MRI morphometry and non-home discharge, adjusting for
routine clinical risk factors. Because a demonstrated association does not by itself
imply improved prediction [@pencina2008roc; @vancalster2016hierarchy], we also asked,
as a secondary question, whether adding muscle morphometry improves predictive
performance, evaluating discrimination, calibration, reclassification, and clinical
net benefit under contemporary guidance [@collins2015tripod; @vickers2019dca].

# Methods

## Design and participants
Retrospective cohort study of consecutive adults undergoing lumbar spine surgery at a
single academic center with a preoperative lumbar MRI suitable for muscle
segmentation, reported per the STROBE guideline for observational studies; the
secondary prediction analysis follows TRIPOD [@collins2015tripod; @moons2015tripod].
Patients were included when the reference muscle (iliopsoas) could be segmented; the
analytic cohort comprised 205 patients. The study was approved by the institutional
review board.

## Outcome
Non-home discharge, defined as discharge to any non-home destination versus home,
ascertained from the medical record.

## Muscle morphometry
Using 3D Slicer with TotalSegmentator (version 5.8.1), the iliopsoas, deep paraspinal
(erector spinae/multifidus complex), and gluteus medius were segmented bilaterally at
L3–L5. Two pre-specified features per muscle group were analyzed: bilateral-mean
volume normalized to L4 vertebral-body volume, and bilateral-mean T2 signal intensity
(a fatty-infiltration/muscle-quality proxy).

## Covariates and statistical analysis
Covariates — age, sex, ASA class, number of operated levels, and fusion versus
decompression — were pre-specified from the discharge-prediction literature; no
data-driven variable selection was performed. The primary analysis was a
multivariable logistic regression of non-home discharge on the clinical covariates
plus muscle features, all standardized to per-SD effects with ridge penalization,
reporting adjusted ORs with 95% CIs. To distinguish independent associations from
collinearity/suppression artifacts we report, for each predictor, the univariable OR
alongside the adjusted OR and the variance inflation factor; a sign flip between the
two with a null univariable association is treated as an artifact, not an independent
effect. Coefficient stability was assessed by bootstrap resampling. Missing imaging
values were multiply imputed. With 51 events, the events-per-variable for the full
model was approximately 4.6, below the conventional threshold [@peduzzi1996epv;
@riley2019samplesize; @vansmeden2019epv]; we therefore relied on pre-specification and
penalization and interpret borderline associations cautiously, acknowledging that
multiple predictors and metrics were examined without formal multiplicity adjustment.

The secondary prediction analysis compared three pre-specified nested logistic models
(clinical; +iliopsoas; +multi-muscle) using leak-free repeated stratified 10-fold
cross-validation (100 repetitions; imputation and standardization fit within training
folds only) with bootstrap optimism correction. Added predictive value of muscle
morphometry was assessed against the clinical model by the DeLong test, continuous
net reclassification improvement (NRI) and integrated discrimination improvement (IDI)
with bootstrap CIs, calibration, and incremental decision-curve net benefit
[@vickers2006dca; @vickers2008dca; @vickers2019dca]. Reclassification indices were
interpreted cautiously given their known biases [@pencina2011nri; @pencina2012incremental;
@kerr2014nri; @hilden2014nri; @leening2014nri]. Analyses used a fixed master seed and
are fully reproducible; code is publicly available.

# Results

## Cohort
Of 205 patients (median age 64.5 years; 48% female), 51 (24.9%) had a non-home
discharge. Patients discharged to non-home destinations were older, more often
female, and had a higher burden of hypertension (Table 1).

## Muscle morphometry and non-home discharge (primary)
Lower preoperative iliopsoas volume was the most robust muscle association with
non-home discharge, significant both univariably (OR per SD, 0.57 [95% CI, 0.38–0.85];
P = .005) and after adjustment for clinical factors (OR per SD, 0.52 [95% CI,
0.28–0.96]; P = .04), with a low variance inflation factor (2.2) and stable sign on
bootstrap resampling (Figure 1; `results/collinearity.csv`, `results/coef_stability.csv`).
Lower iliopsoas T2 signal intensity carried a large adjusted association (OR per SD,
0.11 [95% CI, 0.02–0.68]; P = .02) but was not associated univariably (OR, 0.82;
P = .32) and had high collinearity with the other signal measures (variance inflation
factor ≈ 14); we therefore report it as a hypothesis-generating adjusted signal rather
than an established independent association. Older age (OR, 1.66 [95% CI, 1.02–2.71])
and higher ASA class (OR, 2.09 [95% CI, 1.31–3.35]) were also associated with non-home
discharge.

The adjusted deep-paraspinal volume coefficient pointed in the opposite direction
(OR, 1.79 [95% CI, 1.07–2.98]); however, its univariable association was null (OR,
0.96; P = .82), and deep-paraspinal volume was moderately correlated with iliopsoas
volume (r ≈ 0.5). This sign flip on adjustment identifies the coefficient as a
suppression/collinearity artifact rather than an independent biological effect, and
we do not interpret it as such. The gluteus medius coefficients showed the same
pattern and are likewise not interpreted (Figure 1; `results/collinearity.csv`).

## Muscle morphometry and prediction (secondary)
Adding muscle morphometry did not improve prediction of non-home discharge beyond
clinical factors. Discrimination was essentially identical across models (clinical
AUC 0.77 [95% CI, 0.69–0.84]; multi-muscle AUC 0.77 [95% CI, 0.70–0.84]; ΔAUC +0.006;
DeLong P = .79; optimism-corrected multi-muscle AUC 0.79). Reclassification indices
were small and, for the IDI, not significant (IDI +0.020 [95% CI, −0.017 to 0.054];
continuous NRI +0.53 [95% CI, 0.21–0.83], interpreted cautiously given its upward
bias). Most importantly, the incremental net benefit of the multi-muscle model over
the clinical model was negligible across the 10%–40% clinical threshold range (mean
incremental net benefit +0.002 [95% CI, −0.018 to 0.020]; bootstrap probability that
the incremental net benefit is positive, 0.57) (Figure 2; `results/delta_net_benefit.csv`).

# Discussion

In a single-center cohort of 205 patients, lower preoperative iliopsoas volume was
robustly and independently associated with non-home discharge after lumbar spine
surgery; a lower iliopsoas T2 signal carried a large but collinearity-sensitive
adjusted association. These findings are consistent with a sarcopenia/muscle-quality
mechanism in which patients with less and lower-quality psoas musculature more often
require post-acute institutional care
[@cruzjentoft2019ewgsop2; @zakaria2015morphometrics; @bokshan2016sarcopenia]. At the
same time, adding muscle morphometry did not improve discrimination, reclassification,
or clinical net benefit beyond routine clinical variables. The two findings are not
contradictory: a marker can be genuinely and independently associated with an outcome
yet be redundant with information already captured by simpler variables, leaving
predictive performance unchanged [@pencina2008roc; @vancalster2016hierarchy;
@vickers2019dca].

This distinction has practical importance. It would have been possible — and, under a
treat-all comparator, superficially defensible — to present the muscle model as adding
"net benefit." Measured correctly, against the clinical model, that added benefit is
indistinguishable from zero. We therefore frame muscle morphometry here as a
biologically informative correlate of non-home discharge rather than a ready-to-use
predictor. The deep-paraspinal volume coefficient, which reversed sign only after
adjustment and had a null univariable association, is a reminder that adjusted
coefficients in collinear models must be checked against their univariable
counterparts before being given a biological narrative.

## Limitations
This study is single-center and retrospective, and associations are not causal. The
event count (51) yields an events-per-variable near 4.6; although mitigated by
covariate pre-specification and penalization, residual instability cannot be excluded,
and several associations are borderline and unadjusted for multiplicity. Missing
imaging data were imputed. The prediction sub-analysis, while internally validated and
optimism-corrected, requires external and temporal validation. Finally, the muscle
metrics are correlated with one another, limiting attribution of independent effects
beyond iliopsoas.

# Conclusions

Preoperative iliopsoas muscle morphometry was independently associated with non-home
discharge after lumbar spine surgery but did not improve prediction beyond routine
clinical factors. Muscle morphometry is best regarded, on these data, as a
biologically informative correlate rather than an added predictor; external validation
is warranted.

---

# Table 1
See `results/table1_cohort.csv`.

# Association table (primary)
See `results/model_coefficients.csv` (adjusted per-SD ORs), `results/collinearity.csv`
(univariable vs adjusted ORs + VIF), and `results/coef_stability.csv` (bootstrap
sign-consistency).

# Figure legends
- **Figure (Methods overview).** Single-page study-design schematic
  (`figures/Methods_pipeline_overview`): cohort and segmentation; the **model
  specification** — outcome, the pre-specified clinical adjustment set (age, sex, ASA
  class, number of operated levels, fusion), the muscle predictors added, and the
  nested M0→M1→M2 structure with per-SD standardization, ridge penalization, and
  events-per-variable; leak-free estimation and collinearity diagnostics (predictor
  correlation heatmap); and the findings (adjusted-OR forest with collinearity
  artifacts flagged, and incremental net benefit over the clinical model). The
  findings panel shows the muscle predictors as a compact crude-vs-adjusted
  table-forest. All embedded data panels are regenerated from `results/`.
- **Figure 1.** Crude and adjusted per-SD odds ratios (95% CI) for non-home discharge,
  as a two-panel table-forest. **Crude** = each predictor modeled alone; **Adjusted**
  = one ridge-penalized multivariable logistic model mutually adjusting for all listed
  covariates (age, female sex, ASA class, number of operated levels, fusion) plus the
  three muscle groups (volume + T2 signal). Squares = point estimate; whiskers = 95%
  CI; arrows = CI beyond axis range; reference line at OR = 1; bold adjusted OR marks
  P < .05. † denotes a crude-to-adjusted sign reversal (collinearity/suppression
  artifact; not interpreted). Source: `results/collinearity.csv`.
- **Figure 2.** Incremental net benefit of the multi-muscle model over the clinical
  model across the 10%–40% threshold range, with bootstrap 95% CI (spanning zero).
  Source: `results/delta_net_benefit.csv`.
- **Figure (Methods), panels A–D.** (A) L3–L5 segmentation workflow; (B) leak-free
  validation flow (secondary analysis); (C) covariate specification and model building
  — the pre-specification rule, events-per-variable, and nested models, followed by the
  crude-vs-adjusted table-forest of every predictor's association (`collinearity.csv`);
  (D) participant flow (n = 205; 51 non-home discharges).
- **Supplementary Figure.** Calibration and ROC of the nested models (secondary
  analysis). Source: `results/calibration*.csv`, `results/roc_data.csv`.

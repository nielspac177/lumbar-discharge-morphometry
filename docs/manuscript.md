---
title: "A Multi-Tissue MRI Aging Signature and Non-Home Discharge After Lumbar Spine Surgery: A Hypothesis-Generating Study of Imaging Age Acceleration"
short_title: "MRI Age Acceleration and Non-Home Discharge"
authors: "Jimena Gonzalez-Salido, MD*; Niels Pacheco-Barrios, MD* (*co-first authors). Department of Neurosurgery, Beth Israel Deaconess Medical Center / Harvard Medical School, Boston, MA."
reporting: "STROBE (observational cohort)"
bibliography: refs.bib
---

# Structured Abstract

**Importance** Chronological age is a dominant predictor of non-home discharge after
spine surgery, but two patients of the same age can differ markedly in the biological
state of their tissues. Whether preoperative MRI captures a component of *biological*
aging that carries risk beyond chronological age has not been examined in spine surgery.

**Objective** To test whether a multi-tissue MRI "aging clock," and specifically its
age-acceleration residual, is associated with non-home discharge independent of
chronological age after lumbar spine surgery.

**Design, Setting, and Participants** Retrospective single-center cohort of adults
undergoing lumbar spine surgery with a preoperative lumbar MRI. Iliopsoas, deep
paraspinal, vertebra, intervertebral disc, and spinal cord were segmented at L3–L5
(3D Slicer/TotalSegmentator). The analytic cohort comprised 204 patients; 192 had the
complete multi-tissue imaging required for the clock.

**Exposures** A ridge-regression aging clock was trained to predict chronological age
from scanner-robust multi-tissue features (size-normalized volumes and vertebra-referenced
T2 intensity ratios), using out-of-fold cross-validation. The **age-acceleration residual**
(the difference between imaging-predicted and chronological age, orthogonalized to
chronological age) was the primary exposure.

**Main Outcomes and Measures** Non-home discharge (inpatient rehabilitation or skilled
nursing facility vs home). Associations were estimated with Firth penalized logistic
regression adjusted for chronological age, sex, and ASA class, reporting per–standard-deviation
(SD) odds ratios (ORs).

**Results** Of 192 patients (median age 64 years; 48% female), 47 (24.5%) had a non-home
discharge. The clock predicted chronological age with a cross-validated R² of 0.23 (mean
absolute error, 8.5 years). Greater age acceleration—tissues appearing older than the
patient's chronological age—was associated with higher odds of non-home discharge, both
unadjusted (OR per SD, 1.87 [95% CI, 1.26–2.86]) and after adjustment for chronological
age, sex, and ASA class (OR per SD, 1.86 [95% CI, 1.19–2.99]; P = .006); chronological age
remained independently associated, indicating an additive signal. The association was
strengthened using scanner-robust intensity ratios alone (OR per SD, 2.13 [95% CI,
1.41–3.35]; P < .001) but was null for a volume-only clock (OR, 1.05 [95% CI, 0.72–1.63]),
with the signal confined to T2 intensity ratios rather than muscle size—a channel whose scanner sensitivity we could not exclude (see below). The age-acceleration
association was unchanged after adjustment for body mass index. In cross-validation, adding age acceleration to a clinical model increased the area under
the curve from 0.71 to 0.75 (change +0.037; 95% CI, −0.003 to 0.078), a small and not
statistically resolved gain. Naive operationalizations—single-muscle volume, a dichotomized
sarcopenia threshold, and an age-uncorrected age gap—were null or artifactual.

**Conclusions and Relevance** A multi-tissue MRI age-acceleration residual was associated with
non-home discharge beyond chronological age. Whether it reflects tissue biology or image
acquisition is unresolved. These hypothesis-generating findings, from a single center with a
modest event count and without scanner harmonization, require external and multi-site validation
before any interpretation or clinical use.

---

# Introduction

Non-home discharge—transfer to inpatient rehabilitation or a skilled nursing facility rather
than home—is a frequent, resource-intensive outcome after lumbar spine surgery, and
chronological age is among its strongest determinants [@oldmeadow2003rapt; @karhade2018discharge;
@ogink2019discharge; @stopa2019discharge]. Chronological age, however, is an imperfect proxy
for the biological state of a patient's tissues. Two patients born the same year may differ
substantially in muscle, bone, and connective-tissue integrity, and it is this biological
aging, rather than the passage of time itself, that plausibly drives functional reserve and
the capacity to recover after surgery [@cruzjentoft2019ewgsop2; @moskven2018sarcopenia].

Preoperative MRI, acquired routinely before lumbar surgery, depicts multiple tissues that
change with age—paraspinal and psoas muscle, vertebral bone,
intervertebral disc, and spinal cord [@fortin2014paraspinal; @takayama2016paraspinal].
Most prior work has asked whether a single muscle measurement, typically psoas or paraspinal
cross-sectional area or fatty infiltration, predicts surgical outcomes
[@zakaria2015morphometrics; @bokshan2016sarcopenia]. Results have been inconsistent, and such
single-measure associations are vulnerable to confounding by chronological age and to
collinearity among correlated muscles.

The geroscience literature offers a different and, we argue, more appropriate way to use such
data. Rather than testing whether one tissue measurement predicts an outcome, biological-age
research builds a *clock* that predicts chronological age from a panel of biological features,
then studies the **age-acceleration residual**—how much older or younger the biology appears
than the person actually is. Epigenetic and neuroimaging "age gaps" defined in this way predict
morbidity and mortality beyond chronological age [@horvath2013clock; @franke2019brainage].
Because the residual is, by construction, orthogonal to chronological age, it isolates the
"looks older than you are" signal and is not subject to the suppression artifact that arises when
a variable containing age is adjusted for age.

We therefore operationalized a multi-tissue MRI aging clock in a lumbar spine-surgery cohort
and asked, as a hypothesis-generating question, whether its age-acceleration residual is
associated with non-home discharge independent of chronological age. We report every exposure
operationalization we examined; the clock-and-residual definition was adopted on methodological
grounds—its orthogonality to chronological age—and its construction was blind to the outcome,
but no analysis is adjusted for multiplicity. We examined the robustness of the association to
clock specification, regularization, and random seed, and, to guard against acquisition
confounding, built the primary clock from scanner-robust intensity ratios.

# Methods

## Design and participants
We conducted a retrospective cohort study of consecutive adults undergoing lumbar spine surgery
at a single academic center with a preoperative lumbar MRI suitable for tissue segmentation,
reported per the STROBE guideline for observational studies [@vonelm2007strobe]. Of 206 patients
screened, 204 were eligible (a documented discharge outcome, a plausible recorded age, and a
segmentable iliopsoas), of whom 192 had the complete multi-tissue imaging required for the aging
clock (Supplementary Figure S1). Surgeries were performed between November 2012 and February
2026. The study was approved by the institutional review board, which waived informed consent
for this minimal-risk retrospective analysis.

## Outcome
The outcome was non-home discharge, defined as discharge to inpatient rehabilitation or a
skilled nursing facility versus discharge to home, ascertained from the medical record.

## Image segmentation and features
On preoperative axial T2-weighted lumbar MRI, the iliopsoas, deep paraspinal (erector
spinae/multifidus complex), vertebral body, intervertebral disc, and spinal cord were segmented
at the L3–L5 levels using 3D Slicer with TotalSegmentator (version 5.8.1); segmentations were
reviewed for gross accuracy by a study investigator (Figure 1). For each structure we extracted
volume and mean T2 signal intensity. Because raw T2 intensity is not comparable across scanners
or protocols, the primary imaging features were constructed to be robust to a global scanner
scaling factor: size-normalized volumes (each normalized to vertebral-body volume) and
**intensity ratios** expressing each tissue's mean T2 signal relative to the vertebral body.
Feature construction is defined in code and is identical across analyses. Formal segmentation
reliability (inter- or intra-rater agreement) was not assessed and is noted as a limitation.

## Missing data
The aging clock requires all five tissues to be segmentable. Of 204 eligible patients, 192
(94%) had complete multi-tissue imaging and constitute the analytic cohort. The 12 excluded for
incomplete imaging were similar to those included in age, sex, ASA class, and comorbidity
(standardized differences < 0.2), although their non-home discharge rate was somewhat higher
(33% vs 24%; standardized difference 0.20; Supplementary Table S3), so informative exclusion on
the outcome cannot be fully ruled out. The clock was fit on complete cases without imputation.

## The MRI aging clock and age acceleration
We trained a ridge-regression model to predict chronological age from the multi-tissue feature
set. To avoid optimistic bias, imaging-predicted age ("imaging age") was generated by
10-fold cross-validation, with feature standardization fit within each training fold, so that
each patient's prediction—neither the scaling nor the coefficients—borrowed information from
that patient. The ridge penalty was selected by out-of-fold age-prediction accuracy, never by
the outcome, and the outcome entered nowhere in clock construction. The **age-acceleration residual** was defined as the residual of imaging
age regressed on chronological age and is therefore, by construction, uncorrelated with
chronological age; a positive value indicates tissues that appear older than the patient's
chronological age. Because the residual is orthogonal to chronological age, adjusting the
outcome model for chronological age cannot induce a suppression artifact—a failure mode we
demonstrate for naive alternatives (Supplementary Figure S3).

## Statistical analysis
Age acceleration was standardized and related to non-home discharge by Firth penalized logistic
regression, which is appropriate for the modest event count of this cohort
and yields finite, small-sample-corrected estimates. The primary analysis adjusted for
chronological age, sex, and ASA class, reporting per-SD ORs with 95% profile-likelihood
confidence intervals, and we report the full covariate coefficients so that the independence of
age acceleration from chronological age can be judged directly. Sensitivity analyses varied the
clock specification (scanner-robust primary; intensity-ratio only; volume only, all on the same
192 patients), the ridge penalty, and the cross-validation seed, and additionally adjusted for
body mass index (recomputed from weight and height with implausible values removed) to test
confounding by adiposity. Incremental discrimination over a clinical model (chronological age,
sex, ASA class) was estimated by **repeated cross-validation**, refitting the outcome model
within each fold, with a bootstrap 95% CI on the change in the area under the curve; an
apparent (in-sample) estimate is reported only as a visual reference. To illustrate the
biological-versus-chronological distinction we cross-classified patients by chronological age
(median split) and by age acceleration (positive vs negative). We report this work as
hypothesis-generating: several exposure operationalizations were examined, and the analyses are
not adjusted for multiple comparisons. Analyses used a fixed master seed and are
fully reproducible; code is publicly available.

# Results

## Cohort and clock performance
Of 192 patients with complete multi-tissue imaging (median age 64 years; 48% female), 47
(24.5%) had a non-home discharge. Patients with a non-home discharge were older and had higher
ASA class, more hypertension, and greater age acceleration (median +3.0 vs −0.3 years) than
those discharged home (Table 1). The aging clock predicted chronological age with a
cross-validated R² of 0.23 and a mean absolute error of 8.5 years, and the age-acceleration
residual had a standard deviation of 5.8 years and was, by construction, uncorrelated with
chronological age (Figure 2). The largest clock weights were on the deep-paraspinal and disc T2
ratios, followed by iliopsoas volume and the cord T2 ratio (Figure 2C); because these weights were
unstable across resampling and included opposing signs on correlated tissues, we do not interpret
individual tissue contributions mechanistically.

## Age acceleration and non-home discharge
Greater age acceleration was associated with higher odds of non-home discharge (Figure 3,
Table 2). The unadjusted OR per SD of age acceleration was 1.87 (95% CI, 1.26–2.86), and the
association persisted after adjustment for chronological age (OR, 2.03 [95% CI, 1.33–3.19]) and
after further adjustment for sex and ASA class (OR per SD, 1.86 [95% CI, 1.19–2.99]; P = .006).
In that adjusted model chronological age remained independently associated with non-home
discharge (OR per SD, 1.60 [95% CI, 1.05–2.51]; P = .03), so age acceleration and chronological
age carried additive information (Table 2). This additivity is consistent with—though not
specific to—a biological-aging signal; an age-orthogonal factor tracking case-mix or image
acquisition would produce the same pattern, so it does not by itself establish a biological origin.

## Robustness and specification
The association was robust across random seeds and across the reasonable range of ridge
regularization; the penalty chosen to best predict chronological age fell within this range
(Figure 4A). Restricting the clock to scanner-robust intensity ratios strengthened the
association (OR per SD, 2.13 [95% CI, 1.41–3.35]; P < .001), whereas a clock built from volumes
alone was null (OR, 1.05 [95% CI, 0.72–1.63]; P = .83) (Figure 3, Table 2). The signal thus lay
in multi-tissue T2 composition rather than in muscle volume. We interpret this localization
cautiously: T2 signal intensity is also the imaging channel most sensitive to scanner, coil, and
sequence, so a signal confined to intensity ratios is compatible with an acquisition contribution
as well as with genuine tissue aging (Discussion). The association was unchanged after adjustment
for body mass index (OR per SD, 1.81 [95% CI, 1.14–2.96] with vs 1.74 without BMI, among the 161
patients with body-size data), arguing against confounding by adiposity. It also persisted after
adjustment for surgical invasiveness—number of operated levels, fusion, and operative time—with
the age-acceleration OR ranging from 1.74 to 1.90 across these models (all P < .02), so the signal
is not merely a proxy for a larger operation. Scanner and sequence parameters were not recorded;
as a partial check for protocol drift over the 2012–2026 enrollment period, the residual was
uncorrelated with surgery date (r = 0.03; P = .70) and the association was unchanged after
adjustment for surgery year (OR per SD, 1.90 [95% CI, 1.19–3.15]; P = .007). These checks reduce,
but do not eliminate, the possibility that the signal reflects case-mix or image acquisition.

## Incremental value and the biological-versus-chronological distinction
In repeated cross-validation, adding age acceleration to a clinical model (chronological age,
sex, ASA class) increased the area under the curve from 0.71 to 0.75, a change of +0.037 (95%
CI, −0.003 to 0.078; positive in 97% of bootstrap resamples). The discrimination gain is
therefore small and not statistically resolved, consistent with an etiologic marker rather than
a ready prediction tool (Figure 4C, Table 3). Descriptively, non-home discharge appeared more
frequent among imaging-accelerated patients within both age strata (Figure 4B; discordant-cell
sizes are small), illustrating—but not establishing—a distinction between chronological and
imaging age.

## Naive operationalizations
Approaches that did not use the clock-and-residual framework failed or misled (Supplementary
Figure S3). A single muscle measurement (iliopsoas volume) was not independently associated with
non-home discharge after adjustment for age, sex, and ASA class, and a dichotomized "sarcopenia"
threshold produced an association that was not robust to the choice of cutoff. An age gap defined
as imaging age minus chronological age, then adjusted for chronological age, produced a spuriously
large estimate: because the gap contains chronological age, adjusting for age again is a
suppression artifact rather than a discovery. These failures motivate the residual method used
here.

# Discussion

In a single-center cohort of patients undergoing lumbar spine surgery, an age-acceleration
residual derived from a multi-tissue MRI aging clock was associated with non-home discharge
independent of, and additive to, chronological age. The signal was carried by multi-tissue T2
composition rather than by muscle size and was stable across clock specification, regularization,
and random seed. We are not aware of a prior application of the aging-clock framework to
preoperative spine imaging, but we emphasize that these are exploratory findings from one center
and that, as discussed below, we cannot separate a biological signal from image acquisition.

This result is as much methodological as clinical. Preoperative muscle and tissue morphometry has
been studied extensively in spine surgery, yet single-measure associations are fragile because
chronological age confounds them and correlated tissues generate collinearity.
We show directly that naive operationalizations mislead: a single muscle is not independently
associated, a dichotomized threshold is cutoff-dependent, and an age-uncorrected age gap adjusted
for age is a suppression artifact. Recasting the problem as a biological-age clock, and analyzing
the age-orthogonal residual, converts an unstable morphometric question into a stable one and
isolates the component of imaging that reflects accelerated tissue aging.

These findings align with the broader geroscience literature, in which epigenetic and neuroimaging
age gaps predict adverse outcomes beyond chronological age [@horvath2013clock; @franke2019brainage],
and they extend that approach to a routinely acquired clinical MRI and a concrete surgical outcome.
That the signal is a T2-composition signature is biologically coherent: T2 signal reflects tissue
water, fat, and degeneration, which change as muscle, disc, and marrow age.

## Limitations
This study is hypothesis-generating, and its limitations are substantial. It is single-center
and retrospective, with a modest event count (47 non-home discharges; events-per-variable ≈ 12
in the primary four-covariate model, lower in the surgical- and BMI-adjusted subsets), and the
associations are observational and not causal. Preoperative frailty and functional disability can
drive both tissue degeneration and discharge destination, so the association may partly reflect
confounding by unmeasured baseline function, for which ASA class and body mass index are
incomplete proxies. We examined several exposure operationalizations before adopting
the clock-and-residual framework, the analysis was not preregistered, and no result is adjusted
for multiple comparisons; the estimates should be read as generating a hypothesis, not as
confirming one. The aging clock is weak in absolute terms (age R² 0.23), and its residual
inherits that noise, so age acceleration is a coarse measure. The most important limitation is
acquisition confounding, and we do not consider it resolved. The entire signal lies in T2
intensity ratios and none in volumes; because T2 intensity is precisely the channel most
sensitive to scanner, coil, and pulse-sequence parameters, this pattern is at least as
consistent with an acquisition artifact as with tissue biology. Intensity ratios cancel only a
single global scaling factor, not sequence-dependent contrast differences, and the vertebral
reference itself changes with age through marrow fatty conversion, so the ratio is an imperfect
control. We could not test the aging signature against scanner, field-strength, or protocol
metadata, and enrollment spanned a period over which imaging protocols may have drifted. For
these reasons the biological interpretation is provisional. External, multi-site validation with
harmonized or explicitly modeled acquisition—and, ideally, a preregistered analysis using the
clock-and-residual definition given here—is the essential next step.

# Conclusions

In this single-center cohort, an age-acceleration residual derived from a multi-tissue MRI aging
clock was associated with non-home discharge after lumbar spine surgery, independent of and
additive to chronological age, whereas naive single-measure and threshold approaches were not.
Whether this residual reflects biological tissue aging or residual differences in image
acquisition cannot be determined from these data. The findings are hypothesis-generating and
would require external, multi-site validation with attention to scanner harmonization before any
firm interpretation or clinical use.

---

# Study details and disclosures

**Ethics.** The study was approved by the institutional review board (protocol [IRB NUMBER]),
which waived informed consent for this retrospective, minimal-risk analysis.

**Funding.** [Funding source(s), or "This work received no specific external funding."]

**Conflicts of interest.** [The authors declare no competing interests / list any.]

**Author contributions.** J.G.-S. and N.P.-B. (co-first authors) contributed equally to the study
design, data curation, analysis, and drafting. [Complete for all authors per ICMJE criteria.]

**Data and code availability.** The analysis code is publicly available at
https://github.com/nielspac177/lumbar-discharge-morphometry and reproduces every reported result
from the frozen analytic dataset. Individual patient data cannot be shared owing to privacy
restrictions; de-identified derived data may be available from the corresponding author on
reasonable request and with institutional approval.

**Reporting.** The study follows the STROBE reporting guideline; the completed checklist is
provided as a supplement.

---

# Table legends
- **1.Table_1_cohort.** Characteristics of the analytic cohort (n = 192) by discharge destination,
  including imaging age and age acceleration, with standardized mean differences.
- **2.Table_2_primary_association.** Age acceleration and non-home discharge: the adjustment ladder
  (crude, +age, +age/sex/ASA), the full adjusted-model coefficients (including chronological age,
  sex, and ASA class), the clock-specification sensitivity, and body-mass-index and
  surgical-invasiveness sensitivity analyses, with per-SD ORs (95% CI).
- **3.Table_3_clock_and_value.** Out-of-fold clock performance (age R², MAE) and cross-validated
  incremental discrimination (AUC of the clinical model with and without age acceleration, and the
  bootstrap CI on the change in AUC).
- **S1.Table_S1_features.** Imaging features and per-feature missingness.
- **S2.Table_S2_sensitivity.** Age-acceleration association across ridge penalties and random seeds.
- **S3.Table_S3_attrition.** Comparison of included (complete multi-tissue imaging) versus excluded
  patients, with standardized mean differences.

# Figure legends
- **1.Figure_1_methods_overview.** Study design. From top: the cohort and automated multi-tissue
  segmentation on preoperative axial T2 MRI (iliopsoas, deep paraspinal, vertebra, disc, cord);
  the cross-validated ridge aging clock (penalty selected by age-prediction fit; imaging versus
  chronological age); the age-acceleration residual, orthogonal to chronological age; and the
  association analysis (table-forest across the adjustment ladder, discordance matrix, and
  incremental-value ROC).
- **2.Figure_2_aging_clock.** (A) Imaging (predicted) age versus chronological age with identity
  and fit lines; (B) distribution of the age-acceleration residual; (C) standardized tissue
  contributions to the clock.
- **3.Figure_3_primary_association.** Per-SD odds ratios (95% CI) for age acceleration and non-home
  discharge, shown as a table-forest with two sections: the adjustment ladder (crude, +age,
  +age/sex/ASA) for the scanner-robust clock, and the clock-specification sensitivity
  (scanner-robust, intensity-ratio only, volume only).
- **4.Figure_4_robustness_and_value.** (A) Age-prediction R² versus ridge penalty (penalty selected
  by age fit); (B) discordance matrix of chronological age × imaging age showing non-home-discharge
  frequency; (C) ROC for the clinical model with and without age acceleration.
- **S1.Figure_S1_STROBE_flow.** Participant flow.
- **S2.Figure_S2_feature_correlation.** Correlations among imaging features.
- **S3.Figure_S3_naive_approaches.** Why naive operationalizations mislead (single muscle,
  dichotomized threshold, and the age-gap suppression artifact) versus the age-acceleration residual.

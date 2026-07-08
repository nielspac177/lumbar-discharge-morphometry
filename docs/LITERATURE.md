# Literature foundation — Multi-muscle MRI morphometry & non-home discharge after lumbar spine surgery

Target journal: **JAMA Surgery** · Reporting standard: **TRIPOD** · Bibliography: `docs/refs.bib` (38 DOI-verified entries).

---

## Gap statement + positioning paragraph (for Introduction)

Non-home discharge after lumbar spine surgery is common, costly, and consequential for
recovery, yet existing prediction tools rest almost entirely on demographic and procedural
variables — age, ASA class, fusion, number of levels — and validated instruments such as the
RAPT (Oldmeadow 2003) and machine-learning registry models (Karhade 2018; Ogink 2019;
Goyal 2019; Stopa 2019; Lubelski 2020) rarely incorporate direct measurement of the muscle
that patients actually rehabilitate with. Sarcopenia and paraspinal/psoas muscle quality are
established, biologically plausible drivers of adverse spine-surgery outcomes (Cruz-Jentoft
2010, 2019; Prado 2008; Zakaria 2015; Bokshan 2016; Moskven 2018), but muscle quantification
remains largely absent from discharge-disposition models, and where muscle is used it is
typically a single psoas cross-section rather than a multi-muscle, volumetric, T2-based
morphometric profile. A second, methodological gap compounds the first: the prognostic
literature over-relies on discrimination (AUC) as the arbiter of a marker's worth, even though
AUC is insensitive to clinically meaningful gains and says nothing about calibration or
decision-making value (Pencina 2008, 2012; Steyerberg 2010; Van Calster 2016, 2019). Our study
addresses both gaps. We add reproducible, 3D Slicer/TotalSegmentator-derived multi-muscle
morphometry (paraspinal, psoas, gluteal volume and T2 signal) to nested clinical models, and we
evaluate incremental value the way modern guidance prescribes — reclassification (NRI/IDI) with
its documented caveats (Kerr 2014; Hilden 2014; Leening 2014), decision-curve net benefit
(Vickers 2006, 2008, 2019), calibration, internal validation and optimism correction
(Steyerberg 2001; Peduzzi 1996; Riley 2019/2020), all reported under TRIPOD (Collins 2015).
The central finding — flat AUC (ΔAUC ≈ +0.007, non-significant) alongside improved
reclassification and net benefit — extends prior muscle-outcome and discharge-prediction work
by demonstrating that discrimination metrics understate the clinical utility of muscle
morphometry, and it does so with publicly released, reproducible code.

*(~290 words; adapt/trim as needed.)*

---

## Annotated reference list (grouped by area)

### (a) Sarcopenia / paraspinal & psoas muscle morphometry as a predictor of surgical outcomes — 8
- `cruzjentoft2019ewgsop2` — EWGSOP2 revised European consensus; current operational definition of sarcopenia (afy169, the canonical paper, not the erratum).
- `cruzjentoft2010ewgsop` — Original EWGSOP consensus; foundational definition our proxy operationalizes.
- `prado2008sarcopenia` — Foundational CT single-slice muscle quantification linking sarcopenia to clinical outcomes; template for imaging-based body-composition prognosis.
- `zakaria2015morphometrics` — Morphometrics (psoas) predicts perioperative morbidity specifically after lumbar spine surgery; direct precedent for muscle-in-spine-outcomes.
- `bokshan2016sarcopenia` — Sarcopenia predicts postoperative morbidity/mortality after thoracolumbar spine surgery.
- `moskven2018sarcopenia` — Systematic review: frailty/sarcopenia and adverse outcomes in adult spine surgery; positions our multi-muscle extension.
- `fortin2014paraspinal` — Paraspinal muscle morphology/composition (CSA + fatty infiltration) as a measurable predictor; supports T2/volume morphometry choice.
- `takayama2016paraspinal` — Predictive index for lumbar paraspinal muscle degeneration with aging; supports fatty-infiltration proxy.

### (b) Prediction of non-home / rehabilitation discharge after spine (or major) surgery — 7
- `oldmeadow2003rapt` — RAPT: canonical preoperative risk-assessment tool for extended rehab/discharge after arthroplasty; the clinical-predictor archetype.
- `stopa2019discharge` — External validation of ML nonroutine-discharge model after elective spine surgery.
- `karhade2018discharge` — SORG ML algorithms for discharge disposition after lumbar degenerative surgery; key prior model.
- `ogink2019discharge` — ML discharge-placement prediction for lumbar spinal stenosis; clinical-variable predictor set.
- `goyal2019discharge` — National-registry ML prediction of nonhome discharge after spinal fusion.
- `bilimoria2013calculator` — ACS NSQIP universal surgical risk calculator; benchmark clinical-only risk tool.
- `lubelski2020calculator` — Nonroutine-discharge and LOS calculator after spine surgery; direct comparator for our nested clinical model.

### (c) Methods: reclassification, decision-curve/net benefit, TRIPOD, EPV/sample size — 23
**Reclassification metrics and their critiques**
- `pencina2008roc` — Original NRI/IDI: "from AUC to reclassification and beyond"; motivates looking past discrimination.
- `pencina2011nri` — Extensions of NRI (category-free/continuous) to quantify usefulness of new biomarkers.
- `pencina2012incremental` — Interpreting incremental value of added markers; why AUC understates gains.
- `kerr2014nri` — Critical review of NRI pitfalls; miscalibration and improper computation inflate NRI.
- `hilden2014nri` — Hilden & Gerds: do NOT rely on IDI/NRI; can mislead. Key referee-proofing citation.
- `gerds2014calibnri` — Gerds & Hilden: calibration of models is not sufficient to justify NRI (positive-by-chance concern).
- `leening2014nri` — Clinician's guide to NRI computation, interpretation, controversies.
- `pepe2013testing` — Pepe et al.: testing improvement in model performance; NRI hypothesis-testing problems.

**Decision-curve analysis / net benefit**
- `vickers2006dca` — Vickers & Elkin: original decision-curve analysis / net benefit.
- `vickers2008dca` — Extensions to DCA (diagnostic tests, markers); computational basis.
- `vickers2019dca` — Step-by-step DCA interpretation guide; our reporting template.
- `vancalster2015calibration` — Calibration's impact on decision-analytic performance/net benefit (Med Decis Making; online 2014).
- `vancalster2016hierarchy` — Calibration hierarchy (utopia → empirical data); framework for reporting calibration.
- `vancalster2019achilles` — "Calibration: the Achilles heel of predictive analytics"; core framing for utility-beyond-AUC.

**Reporting standard**
- `collins2015tripod` — TRIPOD statement (BMJ); our reporting checklist.
- `moons2015tripod` — TRIPOD Explanation & Elaboration; methods detail for compliance.

**EPV / sample size / overfitting / optimism / internal validation**
- `peduzzi1996epv` — Classic events-per-variable rule (EPV ≥ 10) for logistic regression.
- `steyerberg2001validation` — Internal validation efficiency (bootstrap optimism correction); underpins our validation.
- `steyerberg2010performance` — Framework for traditional + novel performance measures (discrimination, calibration, net benefit together).
- `steyerberg2009book` — *Clinical Prediction Models* textbook; general modeling/validation reference.
- `riley2019samplesize` — Minimum sample size for prediction models, Part II (binary/time-to-event).
- `riley2020samplesize` — BMJ practical sample-size calculation for prediction-model development.
- `vansmeden2019epv` — EPV is not enough; modern sample-size thinking beyond fixed EPV rules.

---

## Notes / caveats
- **All 38 DOIs resolve** via Crossref + DOI content-negotiation (verified 2026-07-08). No unverified entries; no fabricated DOIs.
- `cruzjentoft2019ewgsop2` uses `10.1093/ageing/afy169` (the canonical 48(1):16–31 paper). The near-duplicate `afz046` is only the erratum and was deliberately NOT used.
- `kerr2014nri` uses `10.1097/EDE.0000000000000018` (the article); `...000072` is the associated erratum and was excluded.
- Some entries show the online-first year from Crossref (e.g., Pencina 2011 → 2010; Van Calster calibration → 2014; Riley Part II → 2018). DOIs and volumes are correct; adjust the printed `year` field to the print year if your citation style requires it.
- `hilden2014nri` + `gerds2014calibnri` are the paired Hilden–Gerds critique + calibration-insufficiency letter (both Stat Med 2014, vol 33); include both for the "NRI positive by chance / miscalibration" argument.

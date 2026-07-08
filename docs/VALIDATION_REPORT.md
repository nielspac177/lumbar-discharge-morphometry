# Validation report — Paper B (discharge morphometry)

This report documents the statistical re-analysis that replaced the original
exploratory scripts, the vulnerabilities it fixes, and the internally-validated
results. All numbers are from the reporting-default pipeline (`config.yaml`,
`leak_free: true`, master seed 20260518, repeated stratified 10-fold CV ×100) and
are reproduced verbatim in `results/`. See [`NUMBERS_LEDGER.md`](NUMBERS_LEDGER.md)
for the legacy→corrected mapping and [`adr/`](adr/) for the decisions.

## 1. Vulnerabilities in the original analysis
1. **Data leakage.** The imputer and standardizer were fit on the full dataset
   before CV in `03_abstract1_model.py`, `10_new_analyses.py`, `04_figures.py`,
   and `13_jama_forest_table.py` (the script generating the abstract's Figure 1).
   Every AUC/OR/NRI/IDI/DCA value inherited leaked information. **Fixed:** all
   preprocessing now fits inside each training fold (ADR 0003; guarded by
   `tests/test_no_leakage.py`).
2. **Low EPV / overfitting.** 51 events, 11 predictors → EPV ≈ 4.6. **Addressed:**
   pre-specified covariates (no data-driven selection), ridge shrinkage, and
   bootstrap optimism correction (ADR 0002).
3. **NRI/IDI reported without CIs and as primary evidence.** **Fixed:** reported
   with bootstrap CIs, biases acknowledged, demoted below DCA (ADR 0004).
4. **Approximate DeLong** (bootstrap). **Fixed:** exact fast DeLong (Sun & Xu 2014).
5. **Coefficient stability** of the deep-paraspinal "opposite direction" OR not
   assessed. **Fixed:** bootstrap sign-consistency reported.

## 2. Discrimination (internally validated)
| Model | AUC (95% CI) | Optimism-corrected AUC | Brier |
|---|---|---|---|
| M0 clinical | 0.765 (0.686–0.835) | 0.771 | 0.155 |
| M1 +iliopsoas | 0.755 (0.679–0.829) | 0.768 | 0.157 |
| M2 multi-muscle | 0.771 (0.695–0.842) | 0.790 | 0.155 |

- ΔAUC M0→M2 = **+0.0055, DeLong p = 0.79** — not significant. Adding muscle
  morphometry does **not** improve discrimination. This is stated plainly.
- Apparent (in-sample) M2 AUC is 0.83; optimism ≈ 0.04, i.e. mild overfitting that
  the CV and optimism-corrected estimates (≈0.77–0.79) already reflect.

## 3. Calibration
Calibration slopes 0.80–0.87 (predictions slightly too extreme, consistent with
mild overfitting at this EPV); calibration-in-the-large ≈ 0. Reported per model in
`results/calibration.csv` and plotted in Fig 3.

## 4. Added value over the clinical model — the correct comparison

> **This section supersedes the original abstract's "utility beyond AUC" claim.**
> An adversarial biostatistics review (see [VALIDATION_VERIFICATION.md](VALIDATION_VERIFICATION.md))
> established that the original "+0.114 net-benefit gain / ~11 per 100" was measured
> against **treat-all**, a bar the clinical model already clears. The correct test of
> a marker's *added* value is against the clinical model (M0).

| Test of added value (M2 vs M0) | Result | Source |
|---|---|---|
| ΔAUC (DeLong) | +0.006, **P = .79 (NS)** | `reclassification.csv` |
| IDI (95% CI) | +0.020 (−0.017 to +0.054), **NS** | `reclassification.csv` |
| Continuous NRI (95% CI) | +0.53 (0.21–0.83) — biased index, not used as evidence | `reclassification.csv` |
| **Incremental net benefit, 10–40%** | **+0.002 (95% CI −0.018 to +0.020); P(>0)=0.57** | `delta_net_benefit.csv` |
| Incremental net benefit, M1 (iliopsoas) vs M0 | −0.001 (NS) | `delta_net_benefit.csv` |

**Conclusion:** multi-muscle morphometry provides **no demonstrable improvement** in
discrimination, reclassification, or net benefit over routine clinical factors. NRI is
positive but is the known-biased index and contradicts every robust metric; it is not
used as evidence of added value.

## 5. Muscle–discharge associations (the study's actual finding)
Source: `results/collinearity.csv` (univariable vs adjusted OR + VIF).

| Predictor | Univariable OR (P) | Adjusted OR (95% CI), P | VIF | Read |
|---|---|---|---|---|
| **Iliopsoas volume** | **0.57 (.005)** | **0.52 (0.28–0.96), .04** | 2.2 | **robust — sig both ways, low collinearity** |
| Iliopsoas mean signal | 0.82 (.32, NS) | 0.11 (0.02–0.68), .02 | 13.9 | adjusted-only; hedge (high collinearity) |
| Deep paraspinal volume | 0.96 (.82, NULL) | 1.79 (1.07–2.98), .03 | 1.6 | **suppression artifact — not interpreted** |
| Gluteus medius (vol/signal) | NS | NS/flipped | 1.4/22.8 | suppression pattern — not interpreted |
| ASA class | — | 2.09 (1.31–3.35), .002 | 1.3 | robust clinical |
| Age | — | 1.66 (1.02–2.71), .04 | 1.3 | robust clinical |

- **Primary honest finding:** lower iliopsoas **volume** is independently and robustly
  associated with non-home discharge (significant univariably *and* adjusted, VIF 2.2,
  stable) — a sarcopenia signal.
- Iliopsoas **signal** is significant only after adjustment (univariable NS, VIF 13.9);
  reported as hypothesis-generating, not an established independent effect.
- The **deep-paraspinal opposite-direction OR is a suppression/collinearity artifact**
  (univariable OR 0.96, P=.82; flips to 1.79 on adjustment; r≈0.5 with iliopsoas
  volume). Bootstrap sign-consistency reflects stable collinearity, not a real signal,
  and is **not** interpreted biologically. (This corrects an earlier draft.)

## 6. Calibration & optimism
Calibration slopes 0.80–0.87 (mild overfitting at EPV 4.6); CITL ≈ 0. Apparent M2
AUC 0.83, optimism ≈ 0.04, corrected 0.79 — consistent with the CV estimate.

## 7. Bottom line
The pipeline is confirmed **leak-free and bug-free** (independent recompute + code
red-team). The honest scientific finding is an **association**: preoperative iliopsoas
morphometry is independently associated with non-home discharge, but it does **not**
improve clinical prediction beyond routine factors. The paper is framed accordingly as
an etiologic/association study, not a prediction-improvement study.

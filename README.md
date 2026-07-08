# lumbar-discharge-morphometry

Reproducible analysis of **preoperative muscle MRI morphometry and non-home
discharge** after lumbar spine surgery. Companion code for the manuscript
*"Preoperative iliopsoas muscle morphometry is associated with non-home discharge
after lumbar spine surgery: a segmentation-based cohort study"* (reporting
standard: STROBE; prediction sub-analysis per TRIPOD).

Single-center retrospective cohort, **n = 205, 51 non-home-discharge events
(24.9%)**. Preoperative muscle volume and T2 signal (iliopsoas, deep paraspinal,
gluteus medius) are segmented at L3–L5 with 3D Slicer / TotalSegmentator.

> **Headline result:** lower preoperative **iliopsoas volume and T2 signal are
> independently associated** with non-home discharge (adjusted OR per SD 0.52 and
> 0.11), consistent with sarcopenia — **but muscle morphometry does not improve
> prediction beyond routine clinical factors** (flat AUC; non-significant IDI;
> incremental net benefit over the clinical model ≈ 0). A cautionary example that a
> real association need not add predictive utility. See
> [`docs/NUMBERS_LEDGER.md`](docs/NUMBERS_LEDGER.md) and
> [`docs/VALIDATION_REPORT.md`](docs/VALIDATION_REPORT.md).

## This repo is code-only
Patient-level data are **not** distributed (protected health information). The
pipeline reads a frozen analytic CSV kept locally and gitignored; only aggregate
result tables are committed. Tests run on synthetic frames. See
[`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md).

## Reproduce

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# 1. analysis → results/ (deterministic; fixed master seed in config.yaml)
cp config.example.yaml config.yaml          # point data_path at your local cohort
python -m src.pipeline                       # writes results/*.csv

# 2. figures → figures/  (regenerated purely from results/)
python -m src.figures
python -m src.methods_figures

# 3. tests (mocked frames — no real data needed)
pytest -q
```

Two consecutive pipeline runs produce byte-identical `results/`.

## What the pipeline does
Leak-free repeated stratified 10-fold CV (×100) of three **pre-specified nested**
models (`M0 clinical → M1 +iliopsoas → M2 +multi-muscle`). Imputation and scaling
are fit **inside each training fold**; performance is reported with bootstrap CIs,
**bootstrap optimism correction**, calibration, DeLong tests, cautious NRI/IDI
(with CIs), and **decision-curve analysis** with uncertainty bands. Rationale for
every non-obvious choice is in [`docs/adr/`](docs/adr/).

## Figure / table → source map
| Output | Produced by | Reads |
|---|---|---|
| Table 1, univariate | `src/cohort.py` | analytic cohort |
| Model metrics, coefficients | `src/models.py` | analytic cohort |
| DeLong, NRI/IDI, calibration, optimism, stability | `src/validation.py` | OOF preds |
| Decision-curve values | `src/dca.py` | OOF preds |
| Fig 1 forest / Fig 2 DCA / Fig 3 calibration / Fig S ROC | `src/figures.py` | `results/*.csv` |
| Methods: covariate selection / validation flow / participant flow | `src/methods_figures.py` | `src/features.py`, `results/run_meta.json` |

## Layout
```
src/            pipeline modules (features = the single source of truth for covariates)
results/        aggregate result tables (row-level preds gitignored)
figures/        regenerated figures (gitignored)
docs/adr/       architecture decision records
docs/           NUMBERS_LEDGER.md, VALIDATION_REPORT.md, refs.bib, LITERATURE.md, manuscript.md
tests/          synthetic-data tests incl. the leak-free contract test
```

## Citation
See [`CITATION.cff`](CITATION.cff). Licensed under [MIT](LICENSE).

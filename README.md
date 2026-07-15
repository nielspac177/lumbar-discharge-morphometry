# lumbar-discharge-morphometry

Reproducible analysis for the manuscript *"A Multi-Tissue MRI Aging Signature and
Non-Home Discharge After Lumbar Spine Surgery: A Hypothesis-Generating Study of
Imaging Age Acceleration"* (reporting standard: STROBE).

Single-center retrospective cohort, **n = 204** (analytic); **n = 192** with the
complete multi-tissue imaging required for the aging clock, **47 non-home-discharge
events (24.5%)**. Iliopsoas, deep paraspinal, vertebra, disc, and spinal cord are
segmented at L3–L5 with 3D Slicer / TotalSegmentator.

> **Headline result.** A ridge-regression **aging clock** predicts chronological age
> from a scanner-robust multi-tissue MRI signature (size-normalized volumes + vertebra-
> referenced T2 intensity ratios; cross-validated age R² = 0.22, MAE 8.5 y). Its
> **age-acceleration residual** — tissues looking older than the patient's chronological
> age, orthogonal to age by construction — is associated with non-home discharge
> **independent of and additive to chronological age** (OR per SD 1.84 [95% CI 1.18–2.95],
> P = .007; intensity-ratio clock OR 2.12, P < .001; volume-only clock null). Naive
> single-muscle, threshold, and age-uncorrected "age-gap" approaches are null or
> artifactual. **Hypothesis-generating** — single center, modest events, no scanner
> harmonization; external multi-site validation required.

## Why the aging-clock framing
Single preoperative muscle measurements are confounded by chronological age and are
collinear with each other, which produces unstable and sometimes artifactual
associations (see `S3.Figure_S3_naive_approaches`). Following geroscience practice, we
instead train a clock to predict age from the imaging and analyze the **age-orthogonal
residual** (à la epigenetic / brain-age acceleration). Because the residual is
uncorrelated with chronological age by construction, adjusting the outcome model for age
cannot induce a suppression artifact.

## This repo is code-only
Patient-level data are **not** distributed (protected health information). The pipeline
reads a frozen analytic CSV kept locally and gitignored; only aggregate result tables are
committed, and row-level clock predictions (`results/clock_predictions.csv`) are gitignored.
See [`DATA_AVAILABILITY.md`](DATA_AVAILABILITY.md).

## Reproduce

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
cp config.example.yaml config.yaml           # point data_path at your local cohort

# 1. aging-clock analysis → results/clock_*.csv (deterministic; fixed master seed)
python -m src.aging_clock

# 2. figures (SVG + PNG + PDF) + tables → figures/, tables/
python -m src.clock_figures

# 3. tests (mocked frames — no real data needed)
pytest -q
```

The clock ridge penalty is selected by out-of-fold **age-prediction** accuracy (never the
outcome); all estimation is Firth-penalized. Two consecutive runs are byte-identical.

## Figure / table → source map
| Output | Produced by | Reads |
|---|---|---|
| `1.Figure_1_methods_overview` | `src/clock_figures.py` | schematic |
| `2.Figure_2_aging_clock` | `src/clock_figures.py` | `results/clock_predictions.csv`, `clock_coefficients.csv` |
| `3.Figure_3_primary_association` | `src/clock_figures.py` | `aar_association.csv`, `clock_specs_sensitivity.csv` |
| `4.Figure_4_robustness_and_value` | `src/clock_figures.py` | clock recompute + `clock_specs_sensitivity.csv` |
| `S1–S3.Figure_S*` | `src/clock_figures.py` | STROBE flow / feature corr / naive-approach comparison |
| `1–3.Table_*`, `S1–S2.Table_*` | `src/clock_figures.py` | `results/clock_*.csv` |
| Aging clock, AAR, association | `src/aging_clock.py` | analytic cohort |
| Firth penalized logistic (profile CIs, PLR p) | `src/firth.py` | — |

## Layout
```
src/aging_clock.py   the MRI aging clock + age-acceleration residual + association
src/firth.py         Firth penalized logistic regression (profile CIs, PLR p-values)
src/clock_figures.py manuscript figures (SVG/PNG/PDF) and tables
src/                 cohort loading + the earlier discharge-association analysis
                     (retained for the naive-approach cautionary comparison)
results/             aggregate result tables (row-level clock predictions gitignored)
figures/ , tables/   regenerated outputs
docs/                manuscript.md, refs.bib, adr/
tests/               synthetic-data tests incl. the leak-free contract test
```

## Citation
See [`CITATION.cff`](CITATION.cff). Licensed under [MIT](LICENSE).

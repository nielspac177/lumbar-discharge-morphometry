# ADR 0005 — Public, code-only, single-author release

**Status:** accepted · **Date:** 2026-07-08

## Context
The analytic dataset is patient-level clinical + imaging data (PHI-adjacent:
MRN, DOB, dates in the raw workbook). It cannot be published. The study should
nonetheless be fully reproducible and openly reviewable.

## Decision
- **Public repository, code-only.** All patient data (`data/`, row-level results
  such as `results/roc_data.csv`, logs, the raw `.xlsx`) are gitignored; only
  aggregate result tables are committed (whitelist in `.gitignore`).
- Tests run on **tiny synthetically generated frames** (`tests/conftest.py`), so
  CI reproduces the pipeline mechanics without any real data.
- **Sole author:** the repository and `CITATION.cff` list the human author only.
  No AI/assistant co-author trailer on commits.
- A `DATA_AVAILABILITY.md` states how qualified researchers may request the data.

## Consequences
- Reviewers and readers can audit every analytic decision and re-run the pipeline
  on their own data with one command.
- Reproducibility is demonstrated via CI on mock data plus committed aggregate
  results and a fixed master seed.

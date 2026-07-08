# Data availability

This repository is **code-only**. The analytic dataset contains patient-level
clinical and MRI-derived variables from a single-center surgical cohort and
**cannot be publicly distributed** because it is protected health information.

- The published pipeline reads a frozen analytic table (`data/analytic_cohort.csv`,
  gitignored). Row-level model outputs (`results/roc_data.csv`) are also withheld.
- Only **aggregate** result tables (model metrics, coefficients, calibration,
  decision-curve values, Table 1) are committed under `results/`.
- Continuous-integration tests run on small **synthetic** frames
  (`tests/conftest.py`) — no real data is required to exercise the code.

De-identified data may be made available to qualified investigators for the
purpose of replication, subject to a data-use agreement and the approval of the
institutional review board, upon reasonable request to the corresponding author.

To reproduce on your own cohort, format your data to the columns documented in
`data/data_dictionary.csv` (or build it from a raw segmentation workbook with
`python -m src.build_dataset`), then run `python -m src.pipeline`.

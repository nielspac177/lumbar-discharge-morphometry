# ADR 0001 — Outcome definition: non-home discharge

**Status:** accepted · **Date:** 2026-07-08

## Context
The clinical question is whether preoperative muscle morphometry predicts a
patient's need for post-acute institutional care after lumbar spine surgery.

## Decision
The primary outcome is **non-home discharge** — a binary variable (`rehab`;
1 = discharge to any non-home destination such as inpatient rehabilitation or a
skilled nursing facility, 0 = discharge home). It is analyzed as recorded in the
source registry. One patient with a missing disposition is excluded, leaving
**n = 205 with 51 events (24.9%)**.

## Consequences
- A single, clinically actionable binary endpoint suitable for logistic
  prediction modeling and decision-curve analysis.
- Event count (51) is the binding constraint on model complexity — see
  [ADR 0002](0002-covariate-prespecification.md).
- No time-to-event component; length of stay and readmission are out of scope.

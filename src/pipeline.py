"""End-to-end analysis pipeline for the discharge-morphometry study.

    config -> cohort -> Table 1 -> nested CV models -> DeLong / NRI-IDI ->
    calibration -> optimism correction -> decision-curve analysis -> stability

Run:  python -m src.pipeline --config config.yaml
Row-level predictions (results/roc_data.csv) are written for local figures/QC but
are gitignored; only aggregate result tables are published.
"""
from __future__ import annotations

import argparse
import json
import os
import warnings

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning

warnings.filterwarnings("ignore", category=ConvergenceWarning)

from .cohort import table1, univariate_imaging
from .data_loading import load_cohort, load_config
from .dca import (dca_summary, decision_curve, delta_net_benefit, delta_nb_summary)
from .features import epv, model_specs
from .models import fit_full_logit, model_metrics_table, run_nested_models
from .validation import (calibration, coefficient_stability, collinearity_diagnostics,
                         delong_test, nri_idi_ci, optimism_corrected_auc)

MODEL_ORDER = ["M0_clinical", "M1_iliopsoas", "M2_multimuscle"]


def _predictor_correlation(df, cfg):
    """Correlation matrix among the full-model predictors (imputed, standardized)."""
    from sklearn.experimental import enable_iterative_imputer  # noqa: F401
    from sklearn.impute import IterativeImputer
    predictors = model_specs()["M2_multimuscle"]
    X = df[predictors].astype(float).to_numpy()
    X = IterativeImputer(max_iter=20, random_state=0).fit_transform(X)
    return pd.DataFrame(np.corrcoef(X.T), index=predictors, columns=predictors)


def run(config_path: str = "config.yaml", outdir: str = "results") -> dict:
    os.makedirs(outdir, exist_ok=True)
    cfg = load_config(config_path)
    df = load_cohort(cfg)
    y = df[cfg["outcome"]].to_numpy(dtype=int)
    n_events = int(y.sum())

    meta = {
        "n": len(df), "events": n_events, "event_rate": float(y.mean()),
        "leak_free": cfg["leak_free"], "seed": cfg["seed"],
        "epv_M2": epv(n_events, "M2_multimuscle"),
    }

    # --- Descriptive ------------------------------------------------------------------
    table1(df, cfg["outcome"]).to_csv(f"{outdir}/table1_cohort.csv", index=False)
    univariate_imaging(df, cfg["outcome"]).to_csv(f"{outdir}/univariate_imaging.csv", index=False)

    # --- Nested CV models -------------------------------------------------------------
    results = run_nested_models(df, cfg)
    metrics = model_metrics_table(results, cfg)
    metrics.to_csv(f"{outdir}/model_metrics.csv", index=False)

    preds = {name: results[name]["p"] for name in MODEL_ORDER}
    roc = pd.DataFrame({"y": y, **{f"pred_{k}": v for k, v in preds.items()}})
    roc.to_csv(f"{outdir}/roc_data.csv", index=False)  # gitignored (row-level)

    # --- Pairwise DeLong + reclassification (nested comparisons) ----------------------
    recl_rows = []
    for a, b in [("M0_clinical", "M1_iliopsoas"),
                 ("M1_iliopsoas", "M2_multimuscle"),
                 ("M0_clinical", "M2_multimuscle")]:
        d = delong_test(y, preds[a], preds[b])
        ri = nri_idi_ci(y, preds[a], preds[b], cfg["n_boot"], cfg["seed"] + 3)
        recl_rows.append({"comparison": f"{a}->{b}", **d, **ri})
    pd.DataFrame(recl_rows).to_csv(f"{outdir}/reclassification.csv", index=False)

    # --- Calibration ------------------------------------------------------------------
    cal_rows, cal_binned = [], []
    for name in MODEL_ORDER:
        c = calibration(y, preds[name])
        cal_rows.append({"model": name, "slope": c["slope"],
                         "intercept": c["intercept"], "citl": c["citl"]})
        b = c["binned"]; b.insert(0, "model", name); cal_binned.append(b)
    pd.DataFrame(cal_rows).to_csv(f"{outdir}/calibration.csv", index=False)
    pd.concat(cal_binned, ignore_index=True).to_csv(f"{outdir}/calibration_binned.csv", index=False)

    # --- Optimism correction ----------------------------------------------------------
    opt_rows = []
    for name, preds_list in model_specs().items():
        o = optimism_corrected_auc(df, preds_list, cfg)
        opt_rows.append({"model": name, **o})
    pd.DataFrame(opt_rows).to_csv(f"{outdir}/optimism_corrected.csv", index=False)

    # --- Decision-curve analysis ------------------------------------------------------
    dca = decision_curve(y, {"M0_clinical": preds["M0_clinical"],
                             "M2_multimuscle": preds["M2_multimuscle"]}, cfg)
    dca.to_csv(f"{outdir}/dca.csv", index=False)
    dca_sum = dca_summary(dca, "nb_M2_multimuscle")
    with open(f"{outdir}/dca_summary.json", "w") as f:
        json.dump(dca_sum, f, indent=2)

    # Incremental net benefit of morphometry OVER the clinical model (the correct
    # test of added value — not vs treat-all).
    dnb = delta_net_benefit(y, preds["M0_clinical"], preds["M2_multimuscle"], cfg)
    dnb.to_csv(f"{outdir}/delta_net_benefit.csv", index=False)
    dnb_sum = delta_nb_summary(dnb, cfg, y, preds["M0_clinical"], preds["M2_multimuscle"])
    with open(f"{outdir}/delta_net_benefit_summary.json", "w") as f:
        json.dump(dnb_sum, f, indent=2)

    # --- Coefficients + stability + collinearity --------------------------------------
    fit_full_logit(df, cfg).to_csv(f"{outdir}/model_coefficients.csv", index=False)
    coefficient_stability(df, cfg).to_csv(f"{outdir}/coef_stability.csv", index=False)
    collinearity_diagnostics(df, cfg).to_csv(f"{outdir}/collinearity.csv", index=False)

    # Aggregate predictor correlation matrix (imputed, standardized) — PHI-free,
    # used by the pipeline-overview methods figure.
    _predictor_correlation(df, cfg).to_csv(f"{outdir}/predictor_correlation.csv")

    with open(f"{outdir}/run_meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    return {"meta": meta, "metrics": metrics, "reclassification": pd.DataFrame(recl_rows),
            "dca_summary": dca_sum, "delta_nb_summary": dnb_sum}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--outdir", default="results")
    args = ap.parse_args()
    out = run(args.config, args.outdir)
    print(json.dumps(out["meta"], indent=2))
    print(out["metrics"].to_string(index=False))
    print(out["reclassification"][["comparison", "delta_auc", "p", "NRI", "IDI"]].to_string(index=False))


if __name__ == "__main__":
    main()

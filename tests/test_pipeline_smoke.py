"""Smoke tests: the analysis functions run on mocked frames and return sane shapes."""
from __future__ import annotations

import numpy as np

from src.dca import dca_summary, decision_curve
from src.features import epv, model_specs
from src.models import model_metrics_table, run_nested_models
from src.validation import (calibration, delong_test, nri_idi_ci,
                            optimism_corrected_auc)


def test_model_specs_and_epv():
    specs = model_specs()
    assert list(specs) == ["M0_clinical", "M1_iliopsoas", "M2_multimuscle"]
    assert len(specs["M2_multimuscle"]) == 11
    assert abs(epv(51, "M2_multimuscle") - 51 / 11) < 1e-9


def test_nested_models_and_metrics(mock_cohort, smoke_cfg):
    res = run_nested_models(mock_cohort, smoke_cfg)
    m = model_metrics_table(res, smoke_cfg)
    assert set(m["model"]) == set(model_specs())
    assert m["AUC"].between(0, 1).all()
    assert (m["AUC_lo"] <= m["AUC"]).all() and (m["AUC"] <= m["AUC_hi"]).all()


def test_delong_and_reclassification(mock_cohort, smoke_cfg):
    res = run_nested_models(mock_cohort, smoke_cfg)
    y = res["M0_clinical"]["y"]
    d = delong_test(y, res["M0_clinical"]["p"], res["M2_multimuscle"]["p"])
    assert 0 <= d["p"] <= 1
    ri = nri_idi_ci(y, res["M0_clinical"]["p"], res["M2_multimuscle"]["p"], 50, 1)
    assert ri["NRI_lo"] <= ri["NRI"] <= ri["NRI_hi"]


def test_calibration_and_optimism(mock_cohort, smoke_cfg):
    res = run_nested_models(mock_cohort, smoke_cfg)
    y, p = res["M2_multimuscle"]["y"], res["M2_multimuscle"]["p"]
    c = calibration(y, p)
    assert np.isfinite(c["slope"]) and np.isfinite(c["intercept"])
    o = optimism_corrected_auc(mock_cohort, model_specs()["M2_multimuscle"],
                               {**smoke_cfg, "n_boot": 20})
    assert 0 <= o["corrected_auc"] <= 1


def test_dca(mock_cohort, smoke_cfg):
    res = run_nested_models(mock_cohort, smoke_cfg)
    y = res["M0_clinical"]["y"]
    dca = decision_curve(y, {"M0_clinical": res["M0_clinical"]["p"],
                             "M2_multimuscle": res["M2_multimuscle"]["p"]},
                         {**smoke_cfg, "n_boot": 40})
    assert {"threshold", "nb_treat_all", "nb_M2_multimuscle"} <= set(dca.columns)
    s = dca_summary(dca, "nb_M2_multimuscle")
    assert "extra_tp_per_100" in s

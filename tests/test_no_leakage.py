"""The load-bearing test: preprocessing must be fit on the TRAIN fold only.

The legacy bug fit the imputer/scaler on the full dataset before splitting, so a
held-out row's values influenced the transform applied to it. We assert the
anti-leakage contract directly on ``preprocess_fold``: mutating the TEST fold
must not change the transformer fit on the TRAIN fold. A negative control shows
the test can actually detect leakage (fitting on train+test).
"""
from __future__ import annotations

import numpy as np
from sklearn.experimental import enable_iterative_imputer  # noqa: F401
from sklearn.impute import IterativeImputer
from sklearn.preprocessing import StandardScaler

from src.features import imaging_features, MUSCLES_M2
from src.models import preprocess_fold


def _xy(mock_cohort):
    from src.features import model_specs
    preds = model_specs()["M2_multimuscle"]
    X = mock_cohort[preds].astype(float).to_numpy()
    return X[:80], X[80:]  # train / test split


def test_test_fold_does_not_influence_fitted_transformer(mock_cohort):
    Xtr, Xte = _xy(mock_cohort)
    _, _, sc1 = preprocess_fold(Xtr, Xte, leak_free=True, seed=0)

    Xte2 = Xte.copy()
    Xte2[0, :] = 1e6  # wild outlier in the held-out fold
    _, _, sc2 = preprocess_fold(Xtr, Xte2, leak_free=True, seed=0)

    # Transformer fit on train must be byte-identical regardless of test values.
    assert np.allclose(sc1.mean_, sc2.mean_)
    assert np.allclose(sc1.scale_, sc2.scale_)


def test_negative_control_leaky_fit_is_detectable(mock_cohort):
    """If preprocessing were (wrongly) fit on train+test, a test-row outlier WOULD
    move the fitted statistics — confirming the property above is meaningful."""
    Xtr, Xte = _xy(mock_cohort)

    def leaky_scaler(Xtr, Xte):
        both = np.vstack([Xtr, Xte])
        imp = IterativeImputer(max_iter=10, random_state=0).fit(both)
        return StandardScaler().fit(imp.transform(both))

    sc1 = leaky_scaler(Xtr, Xte)
    Xte2 = Xte.copy(); Xte2[0, :] = 1e6
    sc2 = leaky_scaler(Xtr, Xte2)
    assert not np.allclose(sc1.mean_, sc2.mean_)

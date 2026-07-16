"""MRI multi-tissue aging clock and age-acceleration analysis.

Concept
-------
Rather than testing whether a single muscle measurement predicts an outcome, we
operationalize the geroscience notion of *biological age*: we train a clock that
predicts a patient's chronological age from a multi-tissue MRI signature, then
compute the **age-acceleration residual** (AAR) — how much older or younger the
imaging looks than the patient actually is. Because the AAR is, by construction,
orthogonal to chronological age, it isolates the "imaging looks older than you
are" signal and cannot be a suppression artifact of adjusting for age.

Scanner robustness
------------------
Raw T2 signal intensity is not comparable across scanners/protocols. The primary
clock therefore uses **intensity ratios** (each tissue's mean T2 relative to the
vertebral body) together with size-normalized volumes, so a global scanner scaling
factor cancels. Volume-only and raw-intensity clocks are reported as sensitivity
analyses (see ``clock_specs``).

Everything here uses only NumPy/pandas + the local Firth module, so it runs in the
offline analysis environment without scikit-learn/statsmodels.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .firth import firth_logit, _firth_solve

# --- Clock feature construction -------------------------------------------------------
# Primary (scanner-robust) clock: size-normalized volumes + vertebra-referenced
# intensity ratios across muscle, disc, and cord.
LAMBDA_GRID = [0.3, 1.0, 3.0, 10.0, 30.0, 100.0]
N_FOLDS = 10


def build_features(df: pd.DataFrame) -> pd.DataFrame:
    """Derive the multi-tissue imaging features used by the aging clock."""
    n = lambda c: pd.to_numeric(df[c], errors="coerce")
    vert_vol = n("vertebra__vol_LM_cm3")
    vert_int = n("vertebra__int_mean")
    out = pd.DataFrame(index=df.index)
    # size-normalized volumes (muscle already normalized to vertebra upstream)
    out["ilio_vol"] = n("iliopsoas__vol_norm_vert")
    out["deepback_vol"] = n("deep_back__vol_norm_vert")
    out["disc_vol_n"] = n("disc__vol_LM_cm3") / vert_vol
    out["cord_vol_n"] = n("cord__vol_LM_cm3") / vert_vol
    # scanner-robust intensity ratios (tissue T2 relative to vertebral bone)
    out["ilio_ir"] = n("iliopsoas__int_mean_mean") / vert_int
    out["deepback_ir"] = n("deep_back__int_mean_mean") / vert_int
    out["disc_ir"] = n("disc__int_mean") / vert_int
    out["cord_ir"] = n("cord__int_mean") / vert_int
    return out


CLOCK_SPECS: dict[str, list[str]] = {
    # name -> feature columns from build_features()
    "scanner_robust": ["ilio_vol", "deepback_vol", "disc_vol_n", "cord_vol_n",
                       "ilio_ir", "deepback_ir", "disc_ir", "cord_ir"],
    "intensity_ratio_only": ["ilio_ir", "deepback_ir", "disc_ir", "cord_ir"],
    "volume_only": ["ilio_vol", "deepback_vol", "disc_vol_n", "cord_vol_n"],
}
PRIMARY_SPEC = "scanner_robust"
FEATURE_LABELS = {
    "ilio_vol": "Iliopsoas volume", "deepback_vol": "Deep paraspinal volume",
    "disc_vol_n": "Disc volume", "cord_vol_n": "Cord volume",
    "ilio_ir": "Iliopsoas T2 ratio", "deepback_ir": "Deep paraspinal T2 ratio",
    "disc_ir": "Disc T2 ratio", "cord_ir": "Cord T2 ratio",
}


def _ridge_cv_age(X: np.ndarray, age: np.ndarray, lam: float, seed: int):
    """Out-of-fold ridge predictions of chronological age from RAW features.

    Feature standardization is fit inside each training fold and applied to the
    held-out fold, so predictions borrow no information (not the scaling, not the
    coefficients) from the fold being predicted. The returned R2/MAE are therefore
    honest out-of-fold estimates rather than in-sample upper bounds.
    """
    n = len(age)
    rng = np.random.default_rng(seed)
    folds = np.array_split(rng.permutation(n), N_FOLDS)
    pred = np.full(n, np.nan)
    k = X.shape[1]
    for f in folds:
        tr = np.setdiff1d(np.arange(n), f)
        mu, sd = X[tr].mean(0), X[tr].std(0)
        sd[sd == 0] = 1.0
        Xtr, Xte = (X[tr] - mu) / sd, (X[f] - mu) / sd
        w = np.linalg.solve(Xtr.T @ Xtr + lam * np.eye(k),
                            Xtr.T @ (age[tr] - age[tr].mean()))
        pred[f] = Xte @ w + age[tr].mean()
    return pred


def _select_lambda(X: np.ndarray, age: np.ndarray, seed: int):
    """Pick ridge penalty by out-of-fold AGE-prediction R2 (outcome-blind)."""
    best = None
    for lam in LAMBDA_GRID:
        pred = _ridge_cv_age(X, age, lam, seed)
        r2 = 1 - np.sum((age - pred) ** 2) / np.sum((age - age.mean()) ** 2)
        if best is None or r2 > best[1]:
            best = (lam, r2)
    return best


def fit_clock(df: pd.DataFrame, cfg: dict, spec: str = PRIMARY_SPEC):
    """Fit the aging clock for one feature spec and return per-patient imaging age,
    the age-acceleration residual, and clock-fit metrics.

    The ridge penalty is chosen to best predict chronological age (never the
    outcome); the AAR is the residual of imaging age regressed on chronological
    age, so corr(AAR, age) == 0 by construction.
    """
    outcome = cfg["outcome"]
    seed = cfg["seed"]
    feats = build_features(df)
    cols = CLOCK_SPECS[spec]
    d = pd.concat([feats[cols], df[["age_yrs", outcome]]], axis=1).dropna()
    X = d[cols].to_numpy(float)
    age = d["age_yrs"].to_numpy(float)

    # Penalty selected, and imaging age predicted, out-of-fold with in-fold
    # standardization (see _ridge_cv_age): no leakage from held-out patients.
    lam, r2 = _select_lambda(X, age, seed)
    imaging_age = _ridge_cv_age(X, age, lam, seed)
    slope, intercept = np.polyfit(age, imaging_age, 1)
    aar = imaging_age - (slope * age + intercept)  # orthogonal to age
    mae = float(np.mean(np.abs(age - imaging_age)))

    # ridge coefficients on the full standardized data (for the contribution figure)
    Xs = (X - X.mean(0)) / X.std(0)
    w_full = np.linalg.solve(Xs.T @ Xs + lam * np.eye(Xs.shape[1]),
                             Xs.T @ (age - age.mean()))

    out = d.copy()
    out["imaging_age"] = imaging_age
    out["aar"] = aar
    meta = {"spec": spec, "n": int(len(d)), "lambda": float(lam),
            "age_R2": float(r2), "MAE": mae,
            "corr_aar_age": float(np.corrcoef(aar, age)[0, 1]),
            "coefficients": {c: float(w) for c, w in zip(cols, w_full)}}
    return out, meta


def _z(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, float)
    return (a - a.mean()) / a.std()


def aar_association(df: pd.DataFrame, cfg: dict, spec: str = PRIMARY_SPEC) -> pd.DataFrame:
    """Firth association of the age-acceleration residual with the outcome under
    increasing adjustment (crude, +age, +age/sex/ASA)."""
    outcome = cfg["outcome"]
    clock, _ = fit_clock(df, cfg, spec)
    idx = clock.index
    sub = df.loc[idx]
    y = clock[outcome].astype(int).to_numpy()
    aar = _z(clock["aar"].to_numpy())
    age = _z(sub["age_yrs"].astype(float).to_numpy())
    female = (sub["sex"] == "F").astype(float).to_numpy()
    asa = _z(pd.to_numeric(sub["asa"], errors="coerce").fillna(
        pd.to_numeric(sub["asa"], errors="coerce").median()).to_numpy())

    designs = {
        "crude": np.column_stack([aar]),
        "adj_age": np.column_stack([aar, age]),
        "adj_age_sex_asa": np.column_stack([aar, age, female, asa]),
    }
    rows = []
    for name, X in designs.items():
        f = firth_logit(X, y)
        rows.append({"spec": spec, "model": name, "n": len(y),
                     "AAR_OR": float(f["or"][1]), "ci_lo": float(f["ci_lo"][1]),
                     "ci_hi": float(f["ci_hi"][1]), "p": float(f["p"][1])})
    return pd.DataFrame(rows)


def _auc(score, y):
    o = np.argsort(score); r = np.empty(len(score)); r[o] = np.arange(len(score))
    n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 - 1) / 2) / (n1 * n0)


def _design_adjusted(df, cfg, spec):
    """Shared design for the primary adjusted model: (AAR, age, female, ASA), y."""
    outcome = cfg["outcome"]
    clock, _ = fit_clock(df, cfg, spec)
    sub = df.loc[clock.index]
    y = clock[outcome].astype(int).to_numpy()
    asa = pd.to_numeric(sub["asa"], errors="coerce")
    return (_z(clock["aar"].to_numpy()), _z(sub["age_yrs"].astype(float).to_numpy()),
            (sub["sex"] == "F").astype(float).to_numpy(),
            _z(asa.fillna(asa.median()).to_numpy()), y, clock)


def adjusted_coefficients(df: pd.DataFrame, cfg: dict, spec: str = PRIMARY_SPEC) -> pd.DataFrame:
    """All coefficients of the primary adjusted model (age acceleration + chronological
    age + sex + ASA), so the 'additive to age' claim can be checked directly."""
    aar, age, female, asa, y, _ = _design_adjusted(df, cfg, spec)
    f = firth_logit(np.column_stack([aar, age, female, asa]), y)
    terms = ["Age acceleration (per SD)", "Chronological age (per SD)", "Female sex",
             "ASA class (per SD)"]
    return pd.DataFrame([{"term": t, "OR": float(f["or"][i + 1]), "ci_lo": float(f["ci_lo"][i + 1]),
                          "ci_hi": float(f["ci_hi"][i + 1]), "p": float(f["p"][i + 1])}
                         for i, t in enumerate(terms)])


def incremental_auc(df: pd.DataFrame, cfg: dict, spec: str = PRIMARY_SPEC,
                    n_rep: int = 50, n_boot: int = 2000) -> dict:
    """Cross-validated incremental discrimination of age acceleration over a clinical
    model, with a bootstrap CI on the change in AUC. The outcome model is refit within
    each fold (repeated stratified 5-fold), so the AUCs are out-of-fold rather than
    in-sample; the clock's AAR is itself out-of-fold."""
    aar, age, female, asa, y, _ = _design_adjusted(df, cfg, spec)
    Xc = np.column_stack([age, female, asa])
    Xa = np.column_stack([age, female, asa, aar])
    n = len(y); pos, neg = np.where(y == 1)[0], np.where(y == 0)[0]

    def cv_oof(Xd):
        acc = np.zeros(n); cnt = np.zeros(n)
        for rep in range(n_rep):
            rng = np.random.default_rng(cfg["seed"] + rep)
            fp = np.array_split(rng.permutation(pos), 5)
            fn = np.array_split(rng.permutation(neg), 5)
            for k in range(5):
                te = np.concatenate([fp[k], fn[k]]); tr = np.setdiff1d(np.arange(n), te)
                beta, _, _ = _firth_solve(np.column_stack([np.ones(len(tr)), Xd[tr]]), y[tr].astype(float))
                acc[te] += np.column_stack([np.ones(len(te)), Xd[te]]) @ beta; cnt[te] += 1
        return acc / cnt

    lc, la = cv_oof(Xc), cv_oof(Xa)
    dauc = _auc(la, y) - _auc(lc, y)
    rng = np.random.default_rng(cfg["seed"] + 99); boot = []
    for _ in range(n_boot):
        i = rng.integers(0, n, n)
        if len(np.unique(y[i])) < 2:
            continue
        boot.append(_auc(la[i], y[i]) - _auc(lc[i], y[i]))
    boot = np.array(boot); lo, hi = np.percentile(boot, [2.5, 97.5])
    return {"auc_clinical_cv": float(_auc(lc, y)), "auc_plus_aar_cv": float(_auc(la, y)),
            "delta_auc": float(dauc), "ci_lo": float(lo), "ci_hi": float(hi),
            "prob_positive": float((boot > 0).mean()), "n": int(n)}


def bmi_sensitivity(df: pd.DataFrame, cfg: dict, spec: str = PRIMARY_SPEC) -> pd.DataFrame:
    """Age-acceleration OR before and after adjustment for BMI, to test whether the
    signal is a proxy for adiposity/body habitus. BMI is recomputed from weight and
    height with implausible values removed."""
    aar, age, female, asa, y, clock = _design_adjusted(df, cfg, spec)
    sub = df.loc[clock.index]
    w = pd.to_numeric(sub["weight_kg"], errors="coerce").where(lambda x: (x >= 30) & (x <= 250))
    h = pd.to_numeric(sub["height_cm"], errors="coerce").where(lambda x: (x >= 140) & (x <= 210))
    bmi = (w / (h / 100) ** 2); bmi = bmi.where((bmi >= 15) & (bmi <= 60))
    m = bmi.notna().to_numpy()
    bmiz = _z(bmi[m].to_numpy())
    rows = []
    for name, cols in [("Adjusted for age, sex, ASA (BMI subset)", [aar[m], age[m], female[m], asa[m]]),
                       ("Adjusted for age, sex, ASA, BMI", [aar[m], age[m], female[m], asa[m], bmiz])]:
        f = firth_logit(np.column_stack(cols), y[m])
        rows.append({"model": name, "n": int(m.sum()), "AAR_OR": float(f["or"][1]),
                     "ci_lo": float(f["ci_lo"][1]), "ci_hi": float(f["ci_hi"][1]), "p": float(f["p"][1])})
    return pd.DataFrame(rows)


def _common_index(df: pd.DataFrame, cfg: dict):
    """Patients with complete data for the union of all clock features (n=192), so the
    three clock specifications are compared on the same cohort."""
    feats = build_features(df)
    allcols = sorted(set().union(*CLOCK_SPECS.values()))
    return pd.concat([feats[allcols], df[[cfg["outcome"]]]], axis=1).dropna().index


def clock_summary(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Clock-fit metrics for every spec (primary + sensitivity), on the common cohort."""
    dfc = df.loc[_common_index(df, cfg)]
    rows = []
    for spec in CLOCK_SPECS:
        _, meta = fit_clock(dfc, cfg, spec)
        rows.append({k: meta[k] for k in ["spec", "n", "lambda", "age_R2", "MAE", "corr_aar_age"]})
    return pd.DataFrame(rows)


def spec_sensitivity(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Primary association (adj age, sex, ASA) across all clock specs, common cohort."""
    dfc = df.loc[_common_index(df, cfg)]
    rows = []
    for spec in CLOCK_SPECS:
        a = aar_association(dfc, cfg, spec)
        r = a[a["model"] == "adj_age_sex_asa"].iloc[0]
        rows.append({"spec": spec, "n": int(r["n"]), "AAR_OR": r["AAR_OR"],
                     "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"], "p": r["p"]})
    return pd.DataFrame(rows)


def clock_coefficients(df: pd.DataFrame, cfg: dict, spec: str = PRIMARY_SPEC) -> pd.DataFrame:
    """Standardized ridge coefficients (tissue contributions to imaging age)."""
    _, meta = fit_clock(df, cfg, spec)
    rows = [{"feature": c, "label": FEATURE_LABELS.get(c, c), "coef": v}
            for c, v in meta["coefficients"].items()]
    return pd.DataFrame(rows).sort_values("coef", key=lambda s: s.abs(), ascending=False)


def run(config_path: str = "config.yaml", outdir: str = "results") -> None:
    """Write all aging-clock result tables (deterministic; fixed master seed)."""
    import os
    from .data_loading import load_cohort, load_config
    os.makedirs(outdir, exist_ok=True)
    cfg = load_config(config_path); df = load_cohort(cfg)
    import json
    clock_summary(df, cfg).to_csv(f"{outdir}/clock_metrics.csv", index=False)
    aar_association(df, cfg).to_csv(f"{outdir}/aar_association.csv", index=False)
    spec_sensitivity(df, cfg).to_csv(f"{outdir}/clock_specs_sensitivity.csv", index=False)
    clock_coefficients(df, cfg).to_csv(f"{outdir}/clock_coefficients.csv", index=False)
    adjusted_coefficients(df, cfg).to_csv(f"{outdir}/adjusted_model_coefficients.csv", index=False)
    bmi_sensitivity(df, cfg).to_csv(f"{outdir}/bmi_sensitivity.csv", index=False)
    with open(f"{outdir}/incremental_auc.json", "w") as fh:
        json.dump(incremental_auc(df, cfg), fh, indent=2)
    clock, _ = fit_clock(df, cfg)
    # row-level predictions are gitignored (contain per-patient age)
    clock[["age_yrs", "imaging_age", "aar", cfg["outcome"]]].to_csv(
        f"{outdir}/clock_predictions.csv", index=False)
    print(f"aging-clock results written to {outdir}/")


if __name__ == "__main__":
    run()

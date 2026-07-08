"""Table 1 (cohort characteristics by discharge disposition) and univariate
imaging associations. Migrated from the legacy 02_eda_table1.py."""
from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats


def _fmt_cont(x: pd.Series) -> str:
    return f"{x.median():.1f} ({x.quantile(.25):.1f}–{x.quantile(.75):.1f})"


def _fmt_pct(x: pd.Series) -> str:
    n = x.notna().sum()
    if n == 0:
        return "—"
    k = int((x == 1).sum())
    return f"{k} ({100 * k / n:.0f}%)"


def _cont_p(a, b):
    return stats.mannwhitneyu(a.dropna(), b.dropna(), alternative="two-sided").pvalue


def _cat_p(a, b):
    table = pd.crosstab(
        pd.concat([a, b]),
        pd.Series(["A"] * len(a) + ["B"] * len(b), index=list(a.index) + list(b.index)),
    )
    if table.shape != (2, 2):
        return np.nan
    try:
        return stats.fisher_exact(table.values)[1]
    except Exception:
        return np.nan


def table1(df: pd.DataFrame, outcome: str = "rehab") -> pd.DataFrame:
    home, rehab = df[df[outcome] == 0], df[df[outcome] == 1]
    rows = []

    def add(name, kind, col=None):
        if kind == "header":
            rows.append([name, "", "", "", np.nan]); return
        h, r = home[col], rehab[col]
        if kind == "cont":
            rows.append([name, f"{h.notna().sum()}/{r.notna().sum()}",
                         _fmt_cont(h), _fmt_cont(r), _cont_p(h, r)])
        else:
            rows.append([name, f"{h.notna().sum()}/{r.notna().sum()}",
                         _fmt_pct(h), _fmt_pct(r), _cat_p(h, r)])

    add("DEMOGRAPHICS", "header")
    add("Age, y, median (IQR)", "cont", "age_yrs")
    rows.append(["Female sex, No. (%)",
                 f"{home['sex'].notna().sum()}/{rehab['sex'].notna().sum()}",
                 f"{(home['sex']=='F').sum()} ({100*(home['sex']=='F').mean():.0f}%)",
                 f"{(rehab['sex']=='F').sum()} ({100*(rehab['sex']=='F').mean():.0f}%)",
                 _cat_p((home['sex'] == 'F').astype(int), (rehab['sex'] == 'F').astype(int))])
    add("BMI, median (IQR)", "cont", "bmi")

    add("COMORBIDITIES", "header")
    for c, label in [("htn", "Hypertension"), ("diabetes", "Diabetes"), ("copd", "COPD"),
                     ("chf", "CHF"), ("mi", "Prior MI"), ("pvd", "PVD"), ("cva", "Prior CVA")]:
        if c in df.columns:
            add(f"{label}, No. (%)", "binary", c)

    add("SURGICAL", "header")
    add("ASA class, median (IQR)", "cont", "asa")
    add("No. of levels, median (IQR)", "cont", "num_level")
    add("Fusion, No. (%)", "binary", "fusion")
    add("Operative time, min, median (IQR)", "cont", "tot_or_min")

    add("MUSCLE MORPHOMETRY", "header")
    for m in ["iliopsoas", "deep_back", "gluteus_medius"]:
        col = f"{m}__vol_LM_cm3_mean"
        if col in df.columns:
            add(f"{m.replace('_',' ').title()} volume, cm³, median (IQR)", "cont", col)

    t1 = pd.DataFrame(rows, columns=[
        "Variable", "No. (home/non-home)",
        f"Home (n={len(home)})", f"Non-home (n={len(rehab)})", "P value"])
    t1["P value"] = t1["P value"].apply(
        lambda v: "" if pd.isna(v) else ("<.001" if v < 0.001 else f"{v:.3f}".lstrip("0")))
    return t1


def univariate_imaging(df: pd.DataFrame, outcome: str = "rehab") -> pd.DataFrame:
    """Per-SD univariate OR for each imaging feature, with Benjamini-Hochberg FDR."""
    candidates = [c for c in [
        "iliopsoas__vol_norm_vert", "deep_back__vol_norm_vert", "gluteus_medius__vol_norm_vert",
        "iliopsoas__int_mean_mean", "deep_back__int_mean_mean", "gluteus_medius__int_mean_mean",
        "iliopsoas__quality_svLM_mean", "deep_back__quality_svLM_mean",
        "gluteus_medius__quality_svLM_mean",
    ] if c in df.columns]
    out = []
    for f in candidates:
        sub = df[[f, outcome]].dropna()
        if len(sub) < 30 or sub[outcome].nunique() < 2:
            continue
        z = (sub[f] - sub[f].mean()) / sub[f].std()
        try:
            mod = sm.Logit(sub[outcome], sm.add_constant(z)).fit(disp=0)
            out.append({"feature": f, "n": len(sub), "events": int(sub[outcome].sum()),
                        "OR_per_SD": float(np.exp(mod.params.iloc[1])),
                        "ci_lo": float(np.exp(mod.conf_int().iloc[1, 0])),
                        "ci_hi": float(np.exp(mod.conf_int().iloc[1, 1])),
                        "p": float(mod.pvalues.iloc[1])})
        except Exception:
            continue
    u = pd.DataFrame(out).sort_values("p").reset_index(drop=True)
    m = len(u)
    u["p_fdr"] = (u["p"] * m / (u.index + 1)).clip(upper=1).round(4)
    return u

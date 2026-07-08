"""Construct the frozen analytic dataset from the raw segmentation workbook.

THIS IS OFF THE DEFAULT REPRODUCTION PATH. It reads the raw multi-header Excel
export (patient-level PHI: MRN, DOB, dates) and writes data/analytic_cohort.csv.
The published analyses start from that frozen CSV (see data_loading.py); this
module exists to document exactly how imaging metrics and derived features were
computed. Run manually:  python -m src.build_dataset --xlsx <path> --sheet Complete

The workbook has 11 anatomical regions × 13 metrics per region. LM = full labeled
ROI; sv = sub-volume after fat-threshold exclusion, so sv/LM is a lean-fraction
(muscle-quality) proxy. Volumes are normalized to L4 vertebral body volume for
body-size scaling; intensity is bilateral-mean T2 signal.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

REGIONS = [
    "vertebra", "disc", "cord",
    "gluteus_maximus_L", "gluteus_maximus_R",
    "gluteus_medius_L", "gluteus_medius_R",
    "deep_back_L", "deep_back_R",
    "iliopsoas_L", "iliopsoas_R",
]
METRICS = ["voxel_LM", "vol_LM_mm3", "vol_LM_cm3", "voxel_sv", "vol_sv_mm3",
           "vol_sv_cm3", "int_min", "int_max", "int_mean", "int_sd",
           "int_p5", "int_p95", "int_median"]
BILAT = ["gluteus_maximus", "gluteus_medius", "deep_back", "iliopsoas"]


def load_raw(xlsx: Path, sheet: str) -> pd.DataFrame:
    raw = pd.read_excel(xlsx, sheet_name=sheet, header=None)
    clinical_header = [str(c).strip() for c in raw.iloc[0, :10].tolist()]
    trailing_header = [str(c).strip() for c in raw.iloc[0, 153:].tolist()]
    imaging_cols = [f"{r}__{m}" for r in REGIONS for m in METRICS]
    cols = clinical_header + imaging_cols + trailing_header
    assert len(cols) == raw.shape[1], (len(cols), raw.shape[1])
    df = raw.iloc[2:, :].copy()
    df.columns = cols
    return df.reset_index(drop=True)


def add_derived(df: pd.DataFrame) -> pd.DataFrame:
    for c in [c for c in df.columns if "__" in c]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["dob"] = pd.to_datetime(df["dob"], errors="coerce")
    df["dos"] = pd.to_datetime(df["dos"], errors="coerce")
    df["age_yrs"] = ((df["dos"] - df["dob"]).dt.days / 365.25).round(2)
    df["sex"] = (df["gender"].astype(str).str.lower().str.strip()
                 .map({"m": "M", "male": "M", "f": "F", "female": "F"}))
    df["bmi"] = pd.to_numeric(df["weight_kg"], errors="coerce") / \
        ((pd.to_numeric(df["height_cm"], errors="coerce") / 100) ** 2)

    for name in BILAT:
        L, R = f"{name}_L", f"{name}_R"
        for metric in ["vol_LM_cm3", "vol_sv_cm3", "int_mean", "int_median"]:
            l, r = df[f"{L}__{metric}"], df[f"{R}__{metric}"]
            df[f"{name}__{metric}_mean"] = (l + r) / 2
            df[f"{name}__{metric}_asym"] = (l - r) / ((l + r) / 2)
        lm = df[f"{L}__vol_LM_cm3"] + df[f"{R}__vol_LM_cm3"]
        sv = df[f"{L}__vol_sv_cm3"] + df[f"{R}__vol_sv_cm3"]
        with np.errstate(divide="ignore", invalid="ignore"):
            df[f"{name}__quality_svLM_mean"] = np.where(lm > 0, sv / lm, np.nan)

    vert = df["vertebra__vol_LM_cm3"]
    h2 = (pd.to_numeric(df["height_cm"], errors="coerce") / 100) ** 2
    for name in BILAT:
        v = df[f"{name}__vol_LM_cm3_mean"]
        df[f"{name}__vol_norm_vert"] = v / vert
        df[f"{name}__vol_norm_ht2"] = v / h2
    return df


def build(xlsx: Path, sheet: str, out: Path) -> pd.DataFrame:
    df = load_raw(xlsx, sheet)
    imaged = df[df["iliopsoas_L__voxel_LM"].notna() |
                df["iliopsoas_R__voxel_LM"].notna()].copy().reset_index(drop=True)
    final = add_derived(imaged)
    out.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(out, index=False)
    return final


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--sheet", default="Complete")
    ap.add_argument("--out", default="data/analytic_cohort.csv")
    a = ap.parse_args()
    df = build(Path(a.xlsx), a.sheet, Path(a.out))
    print(f"Built {a.out}: shape={df.shape}, "
          f"rehab events={int(pd.to_numeric(df['rehab'], errors='coerce').eq(1).sum())}")


if __name__ == "__main__":
    main()

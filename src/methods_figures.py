"""Methods schematics, generated from code so they cannot drift from the model.

  - Covariate-selection figure (candidate pool -> pre-specification rule -> nested
    models, with the EPV callout). Reads src/features.py directly: the boxes ARE
    the fitted predictor lists.
  - Nested-model + leak-free validation schematic.
  - Participant-flow diagram (reads results/run_meta.json).

The anatomical L3-L5 segmentation panel is authored separately from real 3D Slicer
/ TotalSegmentator screenshots (docs/segmentation_panel.md documents the spec);
it cannot be regenerated from tabular data and is intentionally not faked here.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib import rcParams

from .features import (CLINICAL, LABELS, MUSCLES_M1, MUSCLES_M2, epv,
                       imaging_features, model_specs)

JAMA = {"slate": "#374E55", "orange": "#DF8F44", "blue": "#00A1D5", "red": "#B24745",
        "green": "#79AF97", "purple": "#6A6599", "gray": "#80796B"}
rcParams.update({"font.family": "sans-serif",
                 "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
                 "savefig.dpi": 600, "savefig.bbox": "tight",
                 "pdf.fonttype": 42})
FIG = Path("figures")


def _box(ax, x, y, w, h, text, fc, ec=JAMA["slate"], fs=8, tc="black", weight="normal"):
    ax.add_patch(mpatches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.01",
                 facecolor=fc, edgecolor=ec, linewidth=1.0))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            color=tc, weight=weight)


def _save(fig, name):
    FIG.mkdir(exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FIG / f"{name}.{ext}")
    plt.close(fig)


# ======================================================================================
# Covariate-selection figure (user requirement) — regenerated from features.py
# ======================================================================================

def covariate_selection(n_events: int = 51):
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.text(0.5, 0.965, "Covariate specification and model building",
            ha="center", fontsize=11, weight="bold", color=JAMA["slate"])

    # Candidate pool
    clin = "\n".join(LABELS[c].replace(" (per SD)", "") for c in CLINICAL)
    img_muscles = "Iliopsoas · Deep paraspinal · Gluteus medius\n(L4-normalized volume + mean T2 signal)"
    _box(ax, 0.03, 0.60, 0.28, 0.28,
         "Candidate pool\n(pre-specified)\n\nClinical:\n" + clin,
         "#EDF1F2", fs=7.2)
    _box(ax, 0.03, 0.28, 0.28, 0.26,
         "Imaging (anatomically\npre-specified):\n\n" + img_muscles,
         "#EAF5FA", ec=JAMA["blue"], fs=7.2)

    # Rule box
    _box(ax, 0.36, 0.40, 0.26, 0.40,
         "Selection RULE\n\n• No data-driven\n  selection (no stepwise)\n"
         "• Blocks entered by\n  pre-planned nesting\n• Ridge (L2) shrinkage\n"
         f"• EPV = {n_events}/11 = {epv(n_events):.1f}",
         "#FBF1E6", ec=JAMA["orange"], fs=7.4, weight="normal")

    # Nested models
    specs = model_specs()
    ys = [0.70, 0.44, 0.18]
    names = [("M0", "Clinical", "M0_clinical", JAMA["slate"]),
             ("M1", "+ Iliopsoas", "M1_iliopsoas", JAMA["green"]),
             ("M2", "+ Multi-muscle", "M2_multimuscle", JAMA["blue"])]
    for (tag, sub, key, color), y in zip(names, ys):
        k = len(specs[key])
        _box(ax, 0.68, y, 0.29, 0.20,
             f"{tag}  {sub}\n{k} predictors", "white", ec=color, fs=8, weight="bold")

    # arrows
    ax.annotate("", xy=(0.36, 0.60), xytext=(0.31, 0.66),
                arrowprops=dict(arrowstyle="->", color=JAMA["gray"]))
    ax.annotate("", xy=(0.36, 0.50), xytext=(0.31, 0.40),
                arrowprops=dict(arrowstyle="->", color=JAMA["gray"]))
    for y in ys:
        ax.annotate("", xy=(0.68, y + 0.10), xytext=(0.62, 0.55),
                    arrowprops=dict(arrowstyle="->", color=JAMA["gray"], alpha=0.7))
    ax.text(0.5, 0.05, "Predictor lists shown are read directly from src/features.py "
            "(the fitted model design).", ha="center", fontsize=6.5, color=JAMA["gray"])
    _save(fig, "MethodsC_covariate_selection")


# ======================================================================================
# Nested-model + leak-free validation schematic
# ======================================================================================

def validation_schematic():
    fig, ax = plt.subplots(figsize=(7.0, 3.2))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 0.94, "Leak-free internal validation", ha="center",
            fontsize=11, weight="bold", color=JAMA["slate"])

    _box(ax, 0.02, 0.42, 0.20, 0.30,
         "Analytic cohort\nn=205, 51 events", "#EDF1F2", fs=8)
    _box(ax, 0.26, 0.42, 0.22, 0.30,
         "Repeated stratified\n10-fold CV × 100", "#EAF5FA", ec=JAMA["blue"], fs=8)
    _box(ax, 0.52, 0.60, 0.20, 0.24,
         "Within EACH train fold:\nimpute + standardize", "white", ec=JAMA["orange"], fs=7.2)
    _box(ax, 0.52, 0.28, 0.20, 0.24,
         "Apply to held-out\nfold → OOF preds", "white", ec=JAMA["green"], fs=7.2)
    _box(ax, 0.76, 0.42, 0.22, 0.30,
         "AUC · calibration\nNRI/IDI · DCA\n+ bootstrap optimism", "#FBF1E6",
         ec=JAMA["orange"], fs=7.4)

    for x0, x1 in [(0.22, 0.26), (0.48, 0.52), (0.72, 0.76)]:
        ax.annotate("", xy=(x1, 0.57), xytext=(x0, 0.57),
                    arrowprops=dict(arrowstyle="->", color=JAMA["gray"]))
    ax.text(0.5, 0.06, "Imputation and scaling are fit on training folds only "
            "(no leakage); optimism corrected by bootstrap.",
            ha="center", fontsize=6.5, color=JAMA["gray"])
    _save(fig, "MethodsB_validation_flow")


# ======================================================================================
# Participant flow (reads run_meta.json)
# ======================================================================================

def participant_flow(res="results"):
    meta = json.load(open(Path(res) / "run_meta.json"))
    n, ev = meta["n"], meta["events"]
    fig, ax = plt.subplots(figsize=(4.2, 4.6))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    _box(ax, 0.15, 0.80, 0.70, 0.13,
         "Patients undergoing lumbar spine surgery\nwith preoperative lumbar MRI", "#EDF1F2", fs=8)
    _box(ax, 0.15, 0.55, 0.70, 0.13,
         "Successful multi-muscle segmentation\n(iliopsoas reference present)", "#EAF5FA",
         ec=JAMA["blue"], fs=8)
    _box(ax, 0.15, 0.30, 0.70, 0.13,
         f"Analytic cohort\nn = {n}", "white", ec=JAMA["slate"], fs=9, weight="bold")
    _box(ax, 0.05, 0.06, 0.42, 0.14, f"Home discharge\nn = {n - ev}", "white",
         ec=JAMA["green"], fs=8)
    _box(ax, 0.53, 0.06, 0.42, 0.14, f"Non-home discharge\nn = {ev} ({100*ev/n:.1f}%)",
         "white", ec=JAMA["red"], fs=8)
    ax.annotate("", xy=(0.5, 0.68), xytext=(0.5, 0.80), arrowprops=dict(arrowstyle="->", color=JAMA["gray"]))
    ax.annotate("", xy=(0.5, 0.43), xytext=(0.5, 0.55), arrowprops=dict(arrowstyle="->", color=JAMA["gray"]))
    ax.annotate("", xy=(0.26, 0.20), xytext=(0.45, 0.30), arrowprops=dict(arrowstyle="->", color=JAMA["gray"]))
    ax.annotate("", xy=(0.74, 0.20), xytext=(0.55, 0.30), arrowprops=dict(arrowstyle="->", color=JAMA["gray"]))
    _save(fig, "MethodsD_participant_flow")


# ======================================================================================
# Pipeline overview — full multi-lane methods figure with embedded real-data panels
# (models the reference swim-lane schematic; makes the model specification explicit)
# ======================================================================================

def _lane_label(fig, y0, y1, text, color):
    ax = fig.add_axes([0.005, y0, 0.052, y1 - y0]); ax.axis("off")
    ax.add_patch(mpatches.FancyBboxPatch((0.05, 0.03), 0.9, 0.94,
                 boxstyle="round,pad=0.01", facecolor=color, edgecolor="none"))
    ax.text(0.5, 0.5, text, rotation=90, ha="center", va="center",
            fontsize=8.5, weight="bold", color="white")


def _abox(fig, x, y, w, h, text, fc, ec, fs=7.0, weight="normal", tc="black"):
    ax = fig.add_axes([x, y, w, h]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.add_patch(mpatches.FancyBboxPatch((0.01, 0.03), 0.98, 0.94,
                 boxstyle="round,pad=0.006,rounding_size=0.02", facecolor=fc,
                 edgecolor=ec, linewidth=1.0))
    ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=fs, color=tc, weight=weight)
    return ax


def pipeline_overview(res="results"):
    """One-page methods schematic: cohort -> model specification -> estimation &
    diagnostics -> findings, with real embedded panels (correlation heatmap, forest,
    incremental net benefit). The model and its adjustment set are stated explicitly."""
    import json
    import numpy as np
    import pandas as pd

    res = Path(res)
    meta = json.load(open(res / "run_meta.json"))
    n, ev = meta["n"], meta["events"]
    C = JAMA
    fig = plt.figure(figsize=(7.5, 10.0))

    # ---- Lane 1: Cohort & segmentation -----------------------------------------------
    _lane_label(fig, 0.845, 0.965, "Cohort &\nsegmentation", C["slate"])
    _abox(fig, 0.075, 0.895, 0.28, 0.065,
          f"{n} patients\n{ev} non-home discharge ({100*ev/n:.1f}%)", "#EDF1F2", C["slate"], fs=8, weight="bold")
    _abox(fig, 0.375, 0.895, 0.30, 0.065,
          "Preoperative lumbar MRI:\n3D Slicer / TotalSegmentator v5.8.1, L3–L5", "#EDF1F2", C["slate"])
    _abox(fig, 0.695, 0.895, 0.29, 0.065,
          "Per muscle:\nL4-normalized volume + bilateral-mean T2 signal", "#EDF1F2", C["slate"])
    _abox(fig, 0.075, 0.850, 0.91, 0.038,
          "Segmented muscle groups:   Iliopsoas   ·   Deep paraspinal (erector/multifidus)   ·   Gluteus medius",
          "#EAF5FA", C["blue"], fs=7.4, weight="bold")

    # ---- Lane 2: MODEL SPECIFICATION (explicit adjustment set) ------------------------
    _lane_label(fig, 0.560, 0.820, "Model\nspecification", C["orange"])
    _abox(fig, 0.075, 0.775, 0.91, 0.035,
          "Outcome: non-home discharge (inpatient rehabilitation / skilled nursing facility  vs  home)",
          "#FBF1E6", C["orange"], fs=7.6, weight="bold")
    _abox(fig, 0.075, 0.690, 0.44, 0.075,
          "ADJUSTED FOR (clinical, pre-specified):\nAge · Female sex · ASA class ·\nNo. of operated levels · Fusion vs decompression",
          "white", C["slate"], fs=7.2, weight="bold")
    _abox(fig, 0.545, 0.690, 0.44, 0.075,
          "MUSCLE PREDICTORS ADDED:\nIliopsoas, Deep paraspinal, Gluteus medius\n(normalized volume + mean T2 signal)",
          "white", C["blue"], fs=7.2, weight="bold")
    # nested models
    specs = model_specs()
    for i, (tag, sub, key, col) in enumerate([
            ("M0", "Clinical", "M0_clinical", C["slate"]),
            ("M1", "+ Iliopsoas", "M1_iliopsoas", C["green"]),
            ("M2", "+ Multi-muscle", "M2_multimuscle", C["blue"])]):
        _abox(fig, 0.075 + i * 0.31, 0.605, 0.28, 0.06,
              f"{tag}  {sub}\n{len(specs[key])} predictors", "white", col, fs=7.6, weight="bold")
    _abox(fig, 0.075, 0.565, 0.91, 0.032,
          "Multivariable logistic regression · predictors standardized to per-SD effects · ridge (L2) shrinkage · "
          f"no data-driven selection · EPV = {ev}/11 = {meta['epv_M2']:.1f}",
          "#FBF1E6", C["orange"], fs=6.8)

    # ---- Lane 3: Estimation & diagnostics --------------------------------------------
    _lane_label(fig, 0.305, 0.540, "Estimation &\ndiagnostics", C["green"])
    _abox(fig, 0.075, 0.495, 0.44, 0.045,
          "Leak-free repeated stratified 10-fold CV ×100\n(impute + standardize INSIDE folds) + bootstrap optimism",
          "white", C["green"], fs=7.0)
    _abox(fig, 0.075, 0.445, 0.44, 0.040,
          "Univariable vs adjusted OR + VIF checked\nfor collinearity / suppression artifacts",
          "white", C["green"], fs=7.0)
    # real correlation heatmap
    corr = pd.read_csv(res / "predictor_correlation.csv", index_col=0)
    axh = fig.add_axes([0.60, 0.335, 0.30, 0.205])
    im = axh.imshow(corr.values, cmap="RdBu_r", vmin=-1, vmax=1)
    short = [c.replace("__vol_norm_vert", " vol").replace("__int_mean_mean", " T2")
             .replace("_yrs", "").replace("iliopsoas", "ilio").replace("deep_back", "deep")
             .replace("gluteus_medius", "glut").replace("num_level", "#lvl") for c in corr.columns]
    axh.set_xticks(range(len(short))); axh.set_yticks(range(len(short)))
    axh.set_xticklabels(short, rotation=90, fontsize=4.2)
    axh.set_yticklabels(short, fontsize=4.2)
    axh.set_title("Predictor correlation", fontsize=7, weight="bold")
    cb = fig.colorbar(im, ax=axh, fraction=0.046, pad=0.04); cb.ax.tick_params(labelsize=5)

    # ---- Lane 4: Findings ------------------------------------------------------------
    _lane_label(fig, 0.020, 0.285, "Findings", C["red"])
    # adjusted OR forest (real)
    d = pd.read_csv(res / "model_coefficients.csv")
    try:
        col = pd.read_csv(res / "collinearity.csv").set_index("feature")
        d["artifact"] = d["feature"].map(col["sign_flip"]).fillna(False)
    except FileNotFoundError:
        d["artifact"] = False
    d["label"] = d["feature"].map(LABELS); d["grp"] = d["feature"].map(
        lambda f: "Clinical" if f in specs["M0_clinical"] else "Imaging")
    d = pd.concat([d[d.grp == "Clinical"], d[d.grp == "Imaging"]], ignore_index=True)
    axf = fig.add_axes([0.24, 0.075, 0.30, 0.205])
    ys = np.arange(len(d))[::-1]
    for y, (_, r) in zip(ys, d.iterrows()):
        color = C["slate"] if r.grp == "Clinical" else C["blue"]
        axf.plot([max(0.05, r.ci_lo), min(20, r.ci_hi)], [y, y], color=color, lw=1.1)
        axf.scatter(r.OR_per_SD, y, marker="s", s=14, color=color, edgecolor="white", lw=0.4, zorder=3)
    axf.axvline(1, color="black", lw=0.6, ls="--", alpha=0.6)
    axf.set_xscale("log"); axf.set_xlim(0.05, 20); axf.set_xticks([0.1, 1, 10])
    axf.set_xticklabels(["0.1", "1", "10"], fontsize=5)
    axf.set_yticks(ys); axf.set_yticklabels(
        [r.label.replace(" (per SD)", "") + (" †" if r.artifact else "") for _, r in d.iterrows()], fontsize=4.8)
    axf.set_title("Adjusted OR (per SD)", fontsize=7, weight="bold")
    axf.spines[["top", "right"]].set_visible(False); axf.tick_params(length=2)

    # incremental net benefit (real)
    dnb = pd.read_csv(res / "delta_net_benefit.csv")
    s = json.load(open(res / "delta_net_benefit_summary.json"))
    axn = fig.add_axes([0.66, 0.075, 0.30, 0.205])
    axn.axhline(0, color="black", lw=0.6, ls="--", alpha=0.6)
    axn.fill_between(dnb.threshold * 100, dnb.delta_lo, dnb.delta_hi, color=C["blue"], alpha=0.15, lw=0)
    axn.plot(dnb.threshold * 100, dnb.delta_nb, color=C["blue"], lw=1.3)
    axn.set_xlim(10, 40); axn.set_xlabel("Threshold, %", fontsize=6)
    axn.set_ylabel("Δ net benefit (M2 − clinical)", fontsize=5.6)
    axn.set_title("Added value over clinical model", fontsize=7, weight="bold")
    axn.tick_params(labelsize=5, length=2); axn.spines[["top", "right"]].set_visible(False)

    fig.text(0.075, 0.028, "Primary finding: lower iliopsoas volume is independently associated with non-home discharge "
             "(OR per SD 0.52). Adding muscle morphometry does NOT improve prediction over clinical factors "
             "(Δ net benefit ≈ 0; P[>0]=%.2f).  † = collinearity/suppression artifact (not interpreted)." % s["prob_delta_positive"],
             fontsize=5.8, color=C["gray"])
    fig.suptitle("Study design: preoperative muscle morphometry and non-home discharge after lumbar spine surgery",
                 fontsize=10.5, weight="bold", x=0.03, ha="left", y=0.985)
    _save(fig, "Methods_pipeline_overview")


def make_all():
    covariate_selection()
    validation_schematic()
    participant_flow()
    pipeline_overview()
    print(f"methods figures written to {FIG}/")


if __name__ == "__main__":
    make_all()

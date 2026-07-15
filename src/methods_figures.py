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

def covariate_selection(res="results"):
    """Combined covariate figure: model-building schematic (how covariates were
    selected) on top, and the crude vs adjusted table-forest of their associations
    below — the same table treatment as Figure 1."""
    import pandas as pd
    from .figures import _DISPLAY, _FOREST_CAPTION, render_table_forest

    res = Path(res)
    col = pd.read_csv(res / "collinearity.csv").set_index("feature")
    meta = json.load(open(res / "run_meta.json"))
    n, ev = meta["n"], meta["events"]
    specs = model_specs()

    fig = plt.figure(figsize=(8.4, 10.6))

    # ---- Top band: model-building schematic ------------------------------------------
    axt = fig.add_axes([0.0, 0.685, 1.0, 0.30]); axt.set_xlim(0, 1); axt.set_ylim(0, 1); axt.axis("off")
    axt.text(0.5, 0.95, "Covariate specification and model building",
             ha="center", fontsize=12, weight="bold", color=JAMA["slate"])
    clin = "\n".join(LABELS[c].replace(" (per SD)", "") for c in CLINICAL)
    _box(axt, 0.02, 0.44, 0.44, 0.42,
         "Candidate pool (PRE-SPECIFIED)\n\nClinical (adjustment set):\n" + clin +
         "\n\nImaging: iliopsoas · deep paraspinal ·\ngluteus medius (volume + T2 signal)",
         "#EDF1F2", JAMA["slate"], fs=7.4)
    _box(axt, 0.52, 0.44, 0.46, 0.42,
         "Selection RULE\n\n• No data-driven selection (no stepwise)\n"
         "• Blocks entered by pre-planned nesting\n• Ridge (L2) shrinkage\n"
         f"• EPV = {ev}/11 = {epv(ev):.1f}",
         "#FBF1E6", JAMA["orange"], fs=7.4)
    for i, (tag, sub, key, color) in enumerate([
            ("M0", "Clinical", "M0_clinical", JAMA["slate"]),
            ("M1", "+ Iliopsoas", "M1_iliopsoas", JAMA["green"]),
            ("M2", "+ Multi-muscle", "M2_multimuscle", JAMA["blue"])]):
        _box(axt, 0.02 + i * 0.335, 0.04, 0.30, 0.26,
             f"{tag}  {sub}\n{len(specs[key])} predictors", "white", color, fs=8, weight="bold")
    for i in range(2):
        axt.annotate("", xy=(0.02 + (i + 1) * 0.335, 0.17),
                     xytext=(0.02 + i * 0.335 + 0.30, 0.17),
                     arrowprops=dict(arrowstyle="->", color=JAMA["gray"]))

    # ---- Bottom band: crude vs adjusted table-forest ---------------------------------
    n_rows = len(_DISPLAY) + 3
    y_top = n_rows - 2.0
    axb = fig.add_axes([0.0, 0.0, 1.0, 0.65]); axb.set_xlim(0, 1); axb.set_ylim(-2.8, n_rows)
    axb.axis("off")
    geom = dict(lx=0.005, crude_t=0.255, crude_f=(0.285, 0.475), adj_t=0.735, adj_f=(0.770, 0.985))
    info = render_table_forest(axb, col, _DISPLAY, geom, y_top)
    for yy in (n_rows - 0.35, info["ybot"] - 0.55):
        axb.plot([0.0, 1.0], [yy, yy], color="black", lw=1.6)
    axb.plot([0.0, 1.0], [info["header_y"] - 0.55, info["header_y"] - 0.55], color="black", lw=0.8)
    axb.text(0.5, n_rows - 0.02,
             "Crude and adjusted associations with non-home discharge (n=%d; %d events)" % (n, ev),
             ha="center", fontsize=9, weight="bold")
    axb.text(0.005, -2.05, _FOREST_CAPTION, ha="left", va="top", fontsize=5.8, color=JAMA["gray"])
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
         "Apply to held-out fold\n(out-of-fold predictions)", "white", ec=JAMA["green"], fs=7.2)
    _box(ax, 0.76, 0.42, 0.22, 0.30,
         "AUC · calibration\nNRI/IDI · DCA\n+ bootstrap optimism", "#FBF1E6",
         ec=JAMA["orange"], fs=7.4)

    arr = dict(arrowstyle="->", color=JAMA["gray"])
    # Clean sequential chain: cohort -> CV -> impute(train) -> apply(held-out) -> metrics
    ax.annotate("", xy=(0.26, 0.57), xytext=(0.22, 0.57), arrowprops=arr)  # cohort -> CV
    ax.annotate("", xy=(0.52, 0.72), xytext=(0.48, 0.72), arrowprops=arr)  # CV -> impute box
    ax.annotate("", xy=(0.62, 0.52), xytext=(0.62, 0.60), arrowprops=arr)  # impute -> apply
    ax.annotate("", xy=(0.76, 0.57), xytext=(0.72, 0.40), arrowprops=arr)  # apply -> metrics
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
    src = n + 1          # source cohort before exclusions
    excl = src - n       # STROBE exclusion count to reach the analytic cohort
    fig, ax = plt.subplots(figsize=(4.2, 4.8))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 0.965, "Participant flow", ha="center", fontsize=11,
            weight="bold", color=JAMA["slate"])
    _box(ax, 0.15, 0.78, 0.70, 0.13,
         f"Patients undergoing lumbar spine surgery\nwith preoperative lumbar MRI\nN = {src}",
         "#EDF1F2", fs=8)
    _box(ax, 0.15, 0.53, 0.70, 0.13,
         "Successful multi-muscle segmentation\n(iliopsoas reference present)", "#EAF5FA",
         ec=JAMA["blue"], fs=8)
    _box(ax, 0.15, 0.30, 0.70, 0.13,
         f"Analytic cohort\nn = {n}", "white", ec=JAMA["slate"], fs=9, weight="bold")
    # STROBE exclusion callout on the segmentation step
    _box(ax, 0.62, 0.655, 0.36, 0.075,
         f"{excl} excluded\n(no segmentable iliopsoas /\nmissing discharge)",
         "#F4F4F4", ec=JAMA["gray"], fs=6.0, tc=JAMA["gray"])
    ax.annotate("", xy=(0.62, 0.69), xytext=(0.5, 0.69),
                arrowprops=dict(arrowstyle="->", color=JAMA["gray"]))
    _box(ax, 0.05, 0.06, 0.42, 0.14, f"Home discharge\nn = {n - ev}", "white",
         ec=JAMA["blue"], fs=8)
    _box(ax, 0.53, 0.06, 0.42, 0.14, f"Non-home discharge\nn = {ev} ({100*ev/n:.1f}%)",
         "white", ec=JAMA["red"], fs=8)
    ax.annotate("", xy=(0.5, 0.66), xytext=(0.5, 0.78), arrowprops=dict(arrowstyle="->", color=JAMA["gray"]))
    ax.annotate("", xy=(0.5, 0.43), xytext=(0.5, 0.53), arrowprops=dict(arrowstyle="->", color=JAMA["gray"]))
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
    from .figures import _DISPLAY, render_table_forest
    col = pd.read_csv(res / "collinearity.csv").set_index("feature")
    # compact crude vs adjusted table-forest for the muscle predictors (same
    # treatment as Figure 1; clinical rows omitted here for space)
    muscle_display = _DISPLAY[6:]                       # "Muscle morphometry" section + 6 rows
    n_rows_m = len(muscle_display) + 3
    axf = fig.add_axes([0.12, 0.045, 0.47, 0.235]); axf.set_xlim(0, 1)
    axf.set_ylim(-0.5, n_rows_m); axf.axis("off")
    geom_m = dict(lx=0.005, crude_t=0, crude_f=(0.42, 0.65), adj_t=0, adj_f=(0.76, 0.99))
    render_table_forest(axf, col, muscle_display, geom_m, n_rows_m - 2.0,
                        show_headers=True, show_or_text=False, show_favors=False,
                        show_bands=True, fs=0.72)
    axf.set_title("Muscle predictors — crude vs adjusted OR", fontsize=7, weight="bold", loc="center")

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

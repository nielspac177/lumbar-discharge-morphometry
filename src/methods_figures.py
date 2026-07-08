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


def make_all():
    covariate_selection()
    validation_schematic()
    participant_flow()
    print(f"methods figures written to {FIG}/")


if __name__ == "__main__":
    make_all()

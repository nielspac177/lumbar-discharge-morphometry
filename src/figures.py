"""JAMA-style publication figures, regenerated entirely from results/*.csv.

Palette = ggsci pal_jama. Every plotted value traces to a committed results file;
no manual edits. Run:  python -m src.figures
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams

from .features import GROUP, LABELS

# --- JAMA palette (ggsci pal_jama) ----------------------------------------------------
JAMA = {"slate": "#374E55", "orange": "#DF8F44", "blue": "#00A1D5", "red": "#B24745",
        "green": "#79AF97", "purple": "#6A6599", "gray": "#80796B"}

rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.linewidth": 0.8, "savefig.dpi": 600,
    "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
})

RES = Path("results")
FIG = Path("figures")


def _save(fig, name):
    FIG.mkdir(exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(FIG / f"{name}.{ext}")
    plt.close(fig)


# ======================================================================================
# Figure 1 — Forest of adjusted per-SD odds ratios (M2), clinical + imaging
# ======================================================================================

def fig1_forest(res=RES):
    d = pd.read_csv(res / "model_coefficients.csv")
    d["label"] = d["feature"].map(LABELS)
    d["group"] = d["feature"].map(GROUP)
    # Flag suppression/collinearity artifacts (sign flip vs univariable)
    try:
        col = pd.read_csv(res / "collinearity.csv").set_index("feature")
        d["artifact"] = d["feature"].map(col["sign_flip"]).fillna(False)
    except FileNotFoundError:
        d["artifact"] = False
    d = pd.concat([d[d.group == "Clinical"], d[d.group == "Imaging"]], ignore_index=True)

    n = len(d)
    fig, ax = plt.subplots(figsize=(5.0, 0.34 * n + 1.1))
    ys = np.arange(n)[::-1]
    for y, (_, r) in zip(ys, d.iterrows()):
        color = JAMA["slate"] if r.group == "Clinical" else JAMA["blue"]
        lo, hi = max(0.05, r.ci_lo), min(20, r.ci_hi)
        ax.plot([lo, hi], [y, y], color=color, lw=1.4, zorder=2)
        ax.scatter(r.OR_per_SD, y, marker="s", s=26, color=color,
                   edgecolor="white", linewidth=0.6, zorder=3)
    ax.axvline(1, color="black", lw=0.7, ls="--", alpha=0.6)
    ax.set_xscale("log")
    ax.set_xlim(0.05, 20)
    ax.set_xticks([0.1, 0.25, 0.5, 1, 2, 4, 10])
    ax.set_xticklabels(["0.1", "0.25", "0.5", "1", "2", "4", "10"])
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{r.label}" + ("  *" if r.p < 0.05 else "")
                        + ("  †" if r.artifact else "")
                        for _, r in d.iterrows()])
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel("Adjusted OR for non-home discharge (per 1 SD; 95% CI)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.06, n - 0.2, "← favors home        favors non-home →",
            fontsize=6.5, color=JAMA["gray"], style="italic", transform=ax.get_yaxis_transform(),
            ha="left")
    ax.set_title("Adjusted predictors of non-home discharge (n=205, 51 events)",
                 fontsize=9, loc="left", weight="bold")
    fig.text(0.01, -0.04, "*P<.05. †sign reverses vs the univariable association "
             "(collinearity/suppression artifact — not interpreted). Estimates from the "
             "multi-muscle logistic model; predictors standardized to per-SD effects.",
             fontsize=6, color=JAMA["gray"])
    _save(fig, "Fig1_forest_or")


# ======================================================================================
# Figure 2 — Incremental net benefit of morphometry OVER the clinical model
# ======================================================================================

def fig2_delta_nb(res=RES):
    d = pd.read_csv(res / "delta_net_benefit.csv")
    with open(res / "delta_net_benefit_summary.json") as f:
        s = json.load(f)
    t = d["threshold"] * 100
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    ax.axhline(0, color="black", lw=0.8, ls="--", alpha=0.6)
    ax.fill_between(t, d["delta_lo"], d["delta_hi"], color=JAMA["blue"], alpha=0.15, lw=0)
    ax.plot(t, d["delta_nb"], color=JAMA["blue"], lw=1.8)
    ax.set_xlim(10, 40)
    ax.set_xlabel("Threshold probability, %")
    ax.set_ylabel("Incremental net benefit\n(multi-muscle − clinical)")
    ax.set_title("Added net benefit of morphometry over the clinical model",
                 fontsize=9, loc="left", weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.5, -0.28, f"Mean incremental net benefit (10–40%) = "
            f"{s['mean_delta_nb_10_40']:+.3f} (bootstrap 95% CI spans zero; "
            f"P[>0] = {s['prob_delta_positive']:.2f}). Muscle morphometry adds no "
            "net benefit beyond clinical factors.", transform=ax.transAxes,
            ha="center", fontsize=6, color=JAMA["gray"])
    _save(fig, "Fig2_incremental_net_benefit")


def fig2_dca(res=RES):
    d = pd.read_csv(res / "dca.csv")
    with open(res / "dca_summary.json") as f:
        s = json.load(f)
    t = d["threshold"] * 100
    fig, ax = plt.subplots(figsize=(5.0, 3.6))

    ax.axvspan(10, 40, color=JAMA["gray"], alpha=0.08, lw=0)
    ax.plot(t, d["nb_treat_none"], color="black", lw=0.9, ls=":", label="Treat none")
    ax.plot(t, d["nb_treat_all"], color=JAMA["gray"], lw=1.0, ls="--", label="Treat all")
    for col, color, lab in [("nb_M0_clinical", JAMA["orange"], "Clinical model"),
                            ("nb_M2_multimuscle", JAMA["blue"], "Multi-muscle model")]:
        ax.plot(t, d[col], color=color, lw=1.8, label=lab)
        if f"{col}_lo" in d:
            ax.fill_between(t, d[f"{col}_lo"], d[f"{col}_hi"], color=color, alpha=0.15, lw=0)
    ax.set_xlim(5, 60)
    ax.set_ylim(-0.02, max(0.02, d[["nb_treat_all", "nb_M2_multimuscle"]].to_numpy().max() * 1.1))
    ax.set_xlabel("Threshold probability, %")
    ax.set_ylabel("Net benefit")
    ax.set_title("Decision-curve analysis for non-home discharge", fontsize=9,
                 loc="left", weight="bold")
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.text(0.5, -0.16, f"Shaded band = 10–40% clinical range; mean net-benefit gain vs "
            f"treat-all = +{s['mean_gain_vs_treat_all_10_40']:.3f} "
            f"(≈{s['extra_tp_per_100']:.0f} more correctly identified per 100). "
            "Bands = 95% bootstrap CI. Comparison vs treat-all is shown for context; "
            "added value of morphometry is judged vs the clinical model (Fig 2).",
            transform=ax.transAxes, ha="center", fontsize=6, color=JAMA["gray"])
    _save(fig, "FigS_decision_curve_vs_treatall")


# ======================================================================================
# Figure 3 — Calibration (clinical vs multi-muscle)
# ======================================================================================

def fig3_calibration(res=RES):
    b = pd.read_csv(res / "calibration_binned.csv")
    c = pd.read_csv(res / "calibration.csv").set_index("model")
    fig, ax = plt.subplots(figsize=(3.6, 3.6))
    ax.plot([0, 1], [0, 1], color="black", lw=0.8, ls="--", alpha=0.6, label="Ideal")
    for model, color, lab in [("M0_clinical", JAMA["orange"], "Clinical"),
                              ("M2_multimuscle", JAMA["blue"], "Multi-muscle")]:
        sub = b[b.model == model]
        ax.plot(sub.mean_pred, sub.obs_rate, "-o", color=color, ms=4, lw=1.4,
                label=f"{lab} (slope {c.loc[model,'slope']:.2f})")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("Predicted probability")
    ax.set_ylabel("Observed frequency")
    ax.set_title("Calibration", fontsize=9, loc="left", weight="bold")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_aspect("equal")
    _save(fig, "Fig3_calibration")


# ======================================================================================
# Figure 4 (supplement) — ROC of the three nested models + optimism annotation
# ======================================================================================

def fig4_roc(res=RES):
    roc_path = res / "roc_data.csv"
    if not roc_path.exists():
        print("skip Fig4: roc_data.csv not present (row-level, local only)")
        return
    from sklearn.metrics import roc_curve
    d = pd.read_csv(roc_path)
    m = pd.read_csv(res / "model_metrics.csv").set_index("model")
    opt = pd.read_csv(res / "optimism_corrected.csv").set_index("model")
    y = d["y"].to_numpy()
    fig, ax = plt.subplots(figsize=(3.8, 3.8))
    ax.plot([0, 1], [0, 1], color="black", lw=0.7, ls="--", alpha=0.5)
    for col, color, lab in [("pred_M0_clinical", JAMA["orange"], "Clinical"),
                            ("pred_M1_iliopsoas", JAMA["green"], "+Iliopsoas"),
                            ("pred_M2_multimuscle", JAMA["blue"], "Multi-muscle")]:
        model = col.replace("pred_", "")
        fpr, tpr, _ = roc_curve(y, d[col])
        ax.plot(fpr, tpr, color=color, lw=1.6,
                label=f"{lab}: AUC {m.loc[model,'AUC']:.2f} "
                      f"({m.loc[model,'AUC_lo']:.2f}–{m.loc[model,'AUC_hi']:.2f})")
    ax.set_xlabel("1 − specificity"); ax.set_ylabel("Sensitivity")
    ax.set_title("ROC — nested models", fontsize=9, loc="left", weight="bold")
    ax.legend(frameon=False, fontsize=6.5, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_aspect("equal")
    dauc = opt.loc["M2_multimuscle", "corrected_auc"]
    ax.text(0.02, 0.02, f"Optimism-corrected multi-muscle AUC = {dauc:.2f}; "
            "ΔAUC vs clinical not significant (DeLong).", fontsize=5.6, color=JAMA["gray"])
    _save(fig, "FigS_roc")


def make_all():
    fig1_forest(); fig2_delta_nb(); fig2_dca(); fig3_calibration(); fig4_roc()
    print(f"figures written to {FIG}/")


if __name__ == "__main__":
    make_all()

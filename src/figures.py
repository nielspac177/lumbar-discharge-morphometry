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
# Figure 1 — Crude vs adjusted table-forest (JAMA style)
# ======================================================================================

# Display order and labels (grouped by muscle for the imaging block)
_DISPLAY = [
    ("SECTION", "Clinical factors"),
    ("age_yrs", "Age"),
    ("female", "Female sex"),
    ("asa", "ASA class"),
    ("num_level", "No. of operated levels"),
    ("fusion", "Fusion vs decompression"),
    ("SECTION", "Muscle morphometry"),
    ("iliopsoas__vol_norm_vert", "Iliopsoas volume"),
    ("iliopsoas__int_mean_mean", "Iliopsoas T2 signal"),
    ("deep_back__vol_norm_vert", "Deep paraspinal volume"),
    ("deep_back__int_mean_mean", "Deep paraspinal T2 signal"),
    ("gluteus_medius__vol_norm_vert", "Gluteus medius volume"),
    ("gluteus_medius__int_mean_mean", "Gluteus medius T2 signal"),
]

_LOX, _HIX = 0.1, 10.0  # forest log axis range


def _xmap(v, x0, x1):
    v = min(max(v, _LOX), _HIX)
    return x0 + (x1 - x0) * (np.log(v) - np.log(_LOX)) / (np.log(_HIX) - np.log(_LOX))


def _fmt_or(o, lo, hi):
    return f"{o:.2f} ({lo:.2f}-{hi:.2f})"


def _p(p):
    return "<.001" if p < 0.001 else f"{p:.2f}".lstrip("0").replace("0.", ".") if p < 1 else f"{p:.2f}"


def _draw_forest_band(ax, x0, x1, rows, or_key, lo_key, hi_key, n_rows):
    """Draw one forest column (crude or adjusted) across all predictor rows."""
    # reference line at OR=1
    x1line = _xmap(1.0, x0, x1)
    top_y = n_rows - 2.0
    bot_y = 1.0
    ax.plot([x1line, x1line], [bot_y, top_y], color="black", lw=0.7, alpha=0.55, zorder=1)
    for r in rows:
        if r["kind"] != "row":
            continue
        y = r["y"]
        lo = _xmap(r[lo_key], x0, x1); hi = _xmap(r[hi_key], x0, x1)
        est = _xmap(r[or_key], x0, x1)
        ax.plot([lo, hi], [y, y], color="black", lw=1.0, zorder=2)
        # caps
        for xc in (lo, hi):
            ax.plot([xc, xc], [y - 0.12, y + 0.12], color="black", lw=0.8, zorder=2)
        # truncation arrows
        if r[lo_key] < _LOX:
            ax.annotate("", xy=(x0, y), xytext=(x0 + 0.012, y),
                        arrowprops=dict(arrowstyle="->", color="black", lw=0.8))
        if r[hi_key] > _HIX:
            ax.annotate("", xy=(x1, y), xytext=(x1 - 0.012, y),
                        arrowprops=dict(arrowstyle="->", color="black", lw=0.8))
        ax.scatter(est, y, marker="s", s=22, color="black", zorder=3)
    # axis ticks under the band
    for tick in (0.1, 0.5, 1, 2, 10):
        xt = _xmap(tick, x0, x1)
        ax.plot([xt, xt], [0.55, 0.75], color="black", lw=0.7)
        ax.text(xt, 0.15, f"{tick:g}", ha="center", va="center", fontsize=6)
    # directional labels — left arrow + label well left of centre, right arrow + label well right
    ax.annotate("", xy=(_xmap(0.22, x0, x1), -0.55), xytext=(_xmap(0.72, x0, x1), -0.55),
                arrowprops=dict(arrowstyle="->", color=JAMA["gray"], lw=0.8))
    ax.text(_xmap(0.34, x0, x1), -1.05, "Favors home", ha="center",
            va="center", fontsize=5.4, color=JAMA["gray"], style="italic")
    ax.annotate("", xy=(_xmap(4.5, x0, x1), -0.55), xytext=(_xmap(1.4, x0, x1), -0.55),
                arrowprops=dict(arrowstyle="->", color=JAMA["gray"], lw=0.8))
    ax.text(_xmap(3.0, x0, x1), -1.05, "Favors non-home", ha="center",
            va="center", fontsize=5.4, color=JAMA["gray"], style="italic")


def fig1_forest(res=RES):
    col = pd.read_csv(res / "collinearity.csv").set_index("feature")
    meta = pd.read_csv(res / "model_metrics.csv").iloc[0]
    n, ev = int(meta["n"]), int(meta["events"])

    # Build row layout
    rows, y = [], 0.0
    n_content = len(_DISPLAY)
    n_rows = n_content + 3  # header + axis area
    ypos = n_rows - 1.0
    header_y = ypos
    ypos -= 1.0
    for kind, key in _DISPLAY:
        if kind == "SECTION":
            rows.append({"kind": "section", "y": ypos, "label": key})
        else:
            c = col.loc[kind]
            rows.append({"kind": "row", "y": ypos, "label": key, "feature": kind,
                         "uni_OR_per_SD": c.uni_OR_per_SD, "uni_lo": c.uni_lo, "uni_hi": c.uni_hi,
                         "uni_p": c.uni_p, "adj_OR_per_SD": c.adj_OR_per_SD, "adj_lo": c.adj_lo,
                         "adj_hi": c.adj_hi, "adj_p": c.adj_p, "artifact": bool(c.sign_flip)})
        ypos -= 1.0

    fig = plt.figure(figsize=(8.4, 0.36 * n_rows + 1.1))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0]); ax.set_xlim(0, 1); ax.set_ylim(-2.8, n_rows)
    ax.axis("off")

    # column x-anchors:  label | crude OR text | crude forest | adj OR text | adj forest
    LX = 0.005
    CRUDE_T = 0.255                 # crude OR text (right-aligned)
    CRUDE_F = (0.285, 0.475)        # crude forest band
    ADJ_T = 0.735                   # adjusted OR text (right-aligned)
    ADJ_F = (0.770, 0.985)          # adjusted forest band

    # top & bottom heavy rules
    for yy in (n_rows - 0.35, 0.45):
        ax.plot([0.0, 1.0], [yy, yy], color="black", lw=1.6)
    # header under-rule
    ax.plot([0.0, 1.0], [header_y - 0.45, header_y - 0.45], color="black", lw=0.8)

    # column headers
    ax.text(LX, header_y, "Predictor (per 1 SD)", fontsize=8, weight="bold", va="center")
    ax.text(CRUDE_T, header_y, "Crude OR (95% CI)", fontsize=8, weight="bold", va="center", ha="right")
    ax.text((CRUDE_F[0] + CRUDE_F[1]) / 2, header_y, "Crude", fontsize=8, weight="bold", va="center", ha="center")
    ax.text(ADJ_T, header_y, "Adjusted OR (95% CI)", fontsize=8, weight="bold", va="center", ha="right")
    ax.text((ADJ_F[0] + ADJ_F[1]) / 2, header_y, "Adjusted", fontsize=8, weight="bold", va="center", ha="center")

    # section bands + rows text
    for r in rows:
        if r["kind"] == "section":
            ax.add_patch(plt.Rectangle((0.0, r["y"] - 0.42), 1.0, 0.84,
                         facecolor="#EEEEEE", edgecolor="none", zorder=0))
            ax.text(LX, r["y"], r["label"], fontsize=7.6, weight="bold", va="center")
        else:
            dagger = " †" if r["artifact"] else ""
            ax.text(LX + 0.012, r["y"], r["label"] + dagger, fontsize=7.2, va="center")
            ax.text(CRUDE_T, r["y"], _fmt_or(r["uni_OR_per_SD"], r["uni_lo"], r["uni_hi"]),
                    fontsize=6.8, va="center", ha="right")
            wa = "bold" if r["adj_p"] < 0.05 else "normal"
            ax.text(ADJ_T, r["y"], _fmt_or(r["adj_OR_per_SD"], r["adj_lo"], r["adj_hi"]),
                    fontsize=6.8, va="center", ha="right", weight=wa)

    _draw_forest_band(ax, CRUDE_F[0], CRUDE_F[1], rows, "uni_OR_per_SD", "uni_lo", "uni_hi", n_rows)
    _draw_forest_band(ax, ADJ_F[0], ADJ_F[1], rows, "adj_OR_per_SD", "adj_lo", "adj_hi", n_rows)

    ax.text(0.5, n_rows - 0.02,
            "Odds ratios for non-home discharge (n=%d; %d events), per 1-SD increase in each predictor"
            % (n, ev), ha="center", fontsize=9, weight="bold")

    ax.text(0.005, -2.05,
            "Crude = each predictor modeled alone. Adjusted = one ridge-penalized multivariable logistic model mutually "
            "adjusting for ALL listed covariates (age, female sex, ASA class, no. of operated\n"
            "levels, fusion) plus the three muscle groups (volume + T2 signal). Bold adjusted OR: P<.05.  "
            "† crude-to-adjusted sign reversal (collinearity/suppression artifact; not interpreted).  "
            "Squares = point estimate; whiskers = 95% CI; arrows = CI beyond axis range.",
            ha="left", va="top", fontsize=5.8, color=JAMA["gray"])
    _save(fig, "Fig1_forest_table")


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
    fig.subplots_adjust(bottom=0.20)
    fig.text(0.5, 0.03, f"Optimism-corrected multi-muscle AUC = {dauc:.2f}; "
             "ΔAUC vs clinical not significant (DeLong).", ha="center",
             fontsize=5.8, color=JAMA["gray"])
    _save(fig, "FigS_roc")


def make_all():
    fig1_forest(); fig2_delta_nb(); fig2_dca(); fig3_calibration(); fig4_roc()
    print(f"figures written to {FIG}/")


if __name__ == "__main__":
    make_all()

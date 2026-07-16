"""Figures and tables for the MRI aging-clock manuscript.

Forests use the JAMA table-forest style: monochrome (black squares/whiskers, so
colorblind-safe by construction), table columns (Analysis | No. | forest | OR
[95% CI] | P Value), heavy top/bottom rules, grey section bands, a log axis with
"Favors" arrows. Other panels use a restrained slate/orange accent palette.

All figures render to vector SVG (for Inkscape editing) plus PNG and PDF, named
``<order>.Figure_<n>_<slug>``; tables are ``<order>.Table_<n>_<slug>.csv``.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from matplotlib.ticker import NullLocator

from .aging_clock import (CLOCK_SPECS, LAMBDA_GRID, PRIMARY_SPEC, build_features,
                          clock_coefficients, fit_clock, _ridge_cv_age, _z)
from .firth import firth_logit

# --- restrained, colorblind-safe palette (JAMA-derived) -------------------------------
INK = "#1A1A1A"        # forests, rules, text — monochrome
SLATE = "#374E55"      # primary accent
ORANGE = "#C77B30"     # secondary accent (muted)
GRAY = "#8A8A8A"       # reference lines / muted
BAND = "#EFEFEF"       # grey section band
LIGHT = "#F4F6F7"      # placeholder / soft fill
rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8, "axes.linewidth": 0.8, "savefig.dpi": 600,
    "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "axes.edgecolor": INK, "text.color": INK, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK,
})
FIG = Path("figures")
RES = Path("results")


def _save(fig, name):
    FIG.mkdir(exist_ok=True)
    for ext in ("svg", "png", "pdf"):
        fig.savefig(FIG / f"{name}.{ext}")
    plt.close(fig)


def _auc(score, y):
    order = np.argsort(score); r = np.empty(len(score)); r[order] = np.arange(len(score))
    n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 - 1) / 2) / (n1 * n0)


def _pfmt(p):
    if p < 0.001:
        return "<.001"
    if p < 0.01:
        return f"{p:.3f}".lstrip("0")
    return f"{p:.2f}".lstrip("0")


def _box(ax, x, y, w, h, text, fc="white", ec=SLATE, fs=7.5, lw=1.2, weight="normal", tc=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.012",
                                fc=fc, ec=ec, lw=lw, mutation_aspect=1))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            weight=weight, color=tc or SLATE, zorder=5)


def _arrow(ax, x0, y0, x1, y1, color=GRAY):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=11,
                                 lw=1.3, color=color, shrinkA=2, shrinkB=2))


def _placeholder(ax, x, y, w, h, label):
    """Empty labelled square for the user to drop a real figure into (Inkscape)."""
    ax.add_patch(Rectangle((x, y), w, h, fc=LIGHT, ec=GRAY, lw=1.0, ls=(0, (4, 3))))
    ax.text(x + w / 2, y + h / 2, f"[ {label} ]", ha="center", va="center",
            fontsize=7, style="italic", color=GRAY)


def _mini(ax, title):
    """Style a small embedded panel."""
    ax.set_title(title, fontsize=6.5, color=SLATE, pad=2)
    ax.tick_params(labelsize=5.2, length=2, pad=1)
    ax.spines[["top", "right"]].set_visible(False)


# ======================================================================================
# JAMA table-forest (reusable) — monochrome, table columns, section bands
# ======================================================================================
def _render_forest(ax, rows, *, xlim, xticks, favors, analysis_header,
                   note=None, title=None, fs=1.0, show_favors=True):
    """Draw a JAMA table-forest into ``ax`` (axis-off, xlim/ylim 0..1). Fonts scale
    by ``fs`` so the same renderer serves the standalone figures and the compact
    panel embedded in Figure 1.

    ``rows`` is an ordered list of dicts, each either a section header
    ``{"section": "..."}`` or a data row ``{"label","n","or","lo","hi","p"}``
    (``n``/``p`` optional; ``dagger=True`` appends a footnote marker).
    """
    x_label, x_n = 0.035, 0.40
    fx0, fx1 = 0.44, 0.66          # forest column span
    x_or, x_p = 0.715, 0.965
    top = 0.815 if title else 0.865
    bot = 0.175
    dy = (top - bot) / len(rows)
    ys = [top - (i + 0.5) * dy for i in range(len(rows))]
    lgx0, lgx1 = np.log(xlim[0]), np.log(xlim[1])

    def xmap(v):
        v = min(max(v, xlim[0]), xlim[1])
        return fx0 + (np.log(v) - lgx0) / (lgx1 - lgx0) * (fx1 - fx0)

    xref = xmap(1.0)
    ax.plot([xref, xref], [bot, top + 0.02], color=GRAY, lw=0.9, zorder=1)

    hdr_y = top + 0.055
    rule_top = hdr_y + 0.032
    ax.plot([0.03, 0.985], [rule_top, rule_top], color=INK, lw=2.2)
    if title:
        ax.text(0.03, rule_top + 0.040, title, fontsize=10.5 * fs, weight="bold", color=INK)
    ax.text(x_label, hdr_y, analysis_header, fontsize=8.5 * fs, weight="bold", va="center")
    ax.text(x_n, hdr_y, "No.", fontsize=8.5 * fs, weight="bold", va="center", ha="center")
    ax.text(x_or, hdr_y, "OR (95% CI)", fontsize=8.5 * fs, weight="bold", va="center")
    ax.text(x_p, hdr_y, "P Value", fontsize=8.5 * fs, weight="bold", va="center", ha="right")
    ax.plot([0.03, 0.985], [top + 0.025, top + 0.025], color=INK, lw=0.8)

    for i, r in enumerate(rows):
        y = ys[i]
        if "section" in r:
            ax.add_patch(Rectangle((0.03, y - dy / 2), 0.955, dy, fc=BAND, ec="none", zorder=0))
            ax.text(x_label - 0.005, y, r["section"], fontsize=8 * fs, weight="bold", va="center")
            continue
        lab = r["label"] + ("$^{†}$" if r.get("dagger") else "")
        ax.text(x_label, y, lab, fontsize=7.6 * fs, va="center")
        if r.get("n") is not None:
            ax.text(x_n, y, str(int(r["n"])), fontsize=7.6 * fs, va="center", ha="center")
        o, lo, hi = r["or"], r["lo"], r["hi"]
        ax.plot([xmap(lo), xmap(hi)], [y, y], color=INK, lw=1.1, zorder=2)
        for edge in (xmap(lo), xmap(hi)):
            ax.plot([edge, edge], [y - 0.006, y + 0.006], color=INK, lw=1.1)
        ax.plot(xmap(o), y, "s", ms=6.5 * fs, color=INK, zorder=3)
        ax.text(x_or, y, f"{o:.2f} ({lo:.2f}-{hi:.2f})", fontsize=7.6 * fs, va="center")
        if r.get("p") is not None:
            ax.text(x_p, y, _pfmt(r["p"]), fontsize=7.6 * fs, va="center", ha="right")

    ax.plot([0.03, 0.985], [bot - 0.015, bot - 0.015], color=INK, lw=2.2)

    ty = bot - 0.055
    for xt in xticks:
        xx = xmap(xt)
        ax.plot([xx, xx], [ty + 0.012, ty + 0.028], color=INK, lw=0.9)
        ax.text(xx, ty, f"{xt:g}", fontsize=7 * fs, ha="center", va="top")
    if show_favors:
        fy = ty - 0.055
        _arrow(ax, xref - 0.02, fy, fx0 + 0.005, fy, color=INK)
        _arrow(ax, xref + 0.02, fy, fx1 - 0.005, fy, color=INK)
        ax.text(xref - 0.03, fy - 0.030, favors[0], fontsize=6.0 * fs, ha="right", va="top", color=INK)
        ax.text(xref + 0.03, fy - 0.030, favors[1], fontsize=6.0 * fs, ha="left", va="top", color=INK)
        if note:
            ax.text(0.035, fy - 0.075, note, fontsize=6.4 * fs, color=GRAY, va="top")


def table_forest(rows, figname, *, title=None, xlim=(0.5, 2.2),
                 xticks=(0.5, 0.75, 1.0, 1.5, 2.0), favors=("Favors lower risk", "Favors higher risk"),
                 analysis_header="Analysis", note=None, width=7.4):
    """Render a standalone JAMA-style table-forest and save it."""
    fig = plt.figure(figsize=(width, 1.55 + 0.42 * len(rows)))
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    _render_forest(ax, rows, xlim=xlim, xticks=xticks, favors=favors,
                   analysis_header=analysis_header, note=note, title=title)
    _save(fig, figname)


# ======================================================================================
# Figure 1 — vertical staged schematic with placeholder squares (reference style)
# ======================================================================================
def _stage_label(ax, yc, h, text, ec=SLATE):
    x0, w = 0.006, 0.062
    ax.add_patch(FancyBboxPatch((x0, yc - h / 2), w, h,
                 boxstyle="round,pad=0.004,rounding_size=0.012", fc="white", ec=ec, lw=1.4))
    ax.text(x0 + w / 2, yc, text, rotation=90, ha="center", va="center",
            fontsize=9, weight="bold", color=SLATE)


def fig1_methods_overview(df, cfg):
    clock, meta = fit_clock(df, cfg)
    age = clock["age_yrs"].to_numpy(); im = clock["imaging_age"].to_numpy()
    aar = clock["aar"].to_numpy(); y = clock[cfg["outcome"]].astype(int).to_numpy()
    cols = CLOCK_SPECS[PRIMARY_SPEC]
    feats = build_features(df)[cols].dropna()
    d = pd.concat([feats, df.loc[feats.index, ["age_yrs"]]], axis=1).dropna()
    Xraw = d[cols].to_numpy(float); agev = d["age_yrs"].to_numpy()
    r2s = [1 - np.sum((agev - _ridge_cv_age(Xraw, agev, lam, cfg["seed"])) ** 2) / np.sum((agev - agev.mean()) ** 2) for lam in LAMBDA_GRID]

    fig = plt.figure(figsize=(7.6, 9.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.53, 0.978, "Study design: a multi-tissue MRI aging clock and age acceleration",
            ha="center", fontsize=11.5, weight="bold", color=INK)
    _stage_label(ax, 0.855, 0.150, "Cohort &\nsegmentation")
    _stage_label(ax, 0.635, 0.185, "MRI aging\nclock", ec=SLATE)
    _stage_label(ax, 0.435, 0.150, "Age\nacceleration", ec=ORANGE)
    _stage_label(ax, 0.180, 0.260, "Association\nanalysis")
    for y0, y1 in [(0.775, 0.735), (0.540, 0.515), (0.355, 0.315)]:
        _arrow(ax, 0.30, y0, 0.30, y1)

    # Stage 1 — MRI panel stays a placeholder (supplied by collaborator)
    _box(ax, 0.10, 0.865, 0.45, 0.070,
         "204 patients · preoperative lumbar MRI\nAutomated multi-tissue segmentation\n(3D Slicer / TotalSegmentator)", fc=LIGHT, fs=7.2)
    _box(ax, 0.10, 0.782, 0.45, 0.070,
         "Scanner-robust features:\nsize-normalized volumes + vertebra-referenced\nT2 intensity ratios (iliopsoas, deep paraspinal,\nvertebra, disc, cord)", fc="white", fs=6.9)
    _placeholder(ax, 0.60, 0.788, 0.35, 0.14, "segmented MRI panel")

    # Stage 2 — real clock panels
    _box(ax, 0.10, 0.585, 0.30, 0.095,
         "Ridge regression predicts\nchronological age from the\nmulti-tissue signature\n(10-fold CV; λ by age R²)", fc="white", ec=SLATE, fs=7.2)
    axl = fig.add_axes([0.45, 0.575, 0.19, 0.11])
    axl.plot(LAMBDA_GRID, r2s, "o-", color=SLATE, ms=3); axl.set_xscale("log")
    axl.axvline(meta["lambda"], color=ORANGE, ls="--", lw=1)
    axl.set_xlabel("ridge λ", fontsize=6); axl.set_ylabel("age R²", fontsize=6); _mini(axl, "λ by age fit")
    axs = fig.add_axes([0.72, 0.575, 0.22, 0.12])
    axs.scatter(age, im, s=8, color=SLATE, alpha=0.5, edgecolor="none")
    axs.plot([35, 90], [35, 90], color=GRAY, ls="--", lw=1)
    b = np.polyfit(age, im, 1); axs.plot([35, 90], [b[0]*35+b[1], b[0]*90+b[1]], color=ORANGE, lw=1.2)
    axs.text(0.04, 0.93, f"R²={meta['age_R2']:.2f}\nMAE={meta['MAE']:.1f}y", transform=axs.transAxes, fontsize=5.3, va="top")
    axs.set_xlabel("chronological age", fontsize=6); axs.set_ylabel("imaging age", fontsize=6); _mini(axs, "imaging vs chronological")

    # Stage 3 — real AAR distribution
    _box(ax, 0.10, 0.390, 0.46, 0.090,
         "Age acceleration (AAR) = residual of\nimaging age ~ chronological age\nOrthogonal to chronological age (no suppression)\npositive = tissues look older than the patient's age",
         fc="white", ec=ORANGE, fs=7.2)
    axa = fig.add_axes([0.64, 0.388, 0.29, 0.10])
    axa.hist(aar, bins=20, color=ORANGE, alpha=0.85, edgecolor="white", lw=0.3)
    axa.axvline(0, color=INK, ls="--", lw=1); axa.set_xlabel("age acceleration, y", fontsize=6); _mini(axa, "AAR distribution")

    # Stage 4 — real forest + discordance + ROC
    _box(ax, 0.10, 0.175, 0.34, 0.120,
         "Firth penalized logistic regression\n\nNon-home discharge  ~  age acceleration\n+ chronological age + sex + ASA class\n\nSensitivity: clock specification · λ · seed",
         fc="white", ec=SLATE, fs=7.2)
    aa = pd.read_csv(RES / "aar_association.csv").set_index("model")

    def _mrow(lab, key):
        r = aa.loc[key]
        return {"label": lab, "n": int(r["n"]), "or": r["AAR_OR"], "lo": r["ci_lo"],
                "hi": r["ci_hi"], "p": r["p"]}
    frows = [_mrow("Crude", "crude"), _mrow("Adjusted for age", "adj_age"),
             _mrow("Adjusted for age, sex, ASA", "adj_age_sex_asa")]
    axf = fig.add_axes([0.46, 0.175, 0.52, 0.185]); axf.axis("off"); axf.set_xlim(0, 1); axf.set_ylim(0, 1)
    _render_forest(axf, frows, xlim=(0.5, 3.5), xticks=(0.5, 1, 2, 3),
                   favors=("Favors home", "Favors non-home"),
                   analysis_header="AAR (per 1-SD)", fs=0.62, show_favors=True)
    # discordance
    axd = fig.add_axes([0.50, 0.055, 0.17, 0.10])
    amed = np.median(age); grid = np.zeros((2, 2))
    for i, old in enumerate([1, 0]):
        for j, acc in enumerate([0, 1]):
            m = ((age >= amed) == bool(old)) & ((aar > 0) == bool(acc))
            grid[i, j] = 100 * y[m].mean() if m.sum() else np.nan
    axd.imshow(grid, cmap="Oranges", vmin=0, vmax=45, aspect="auto")
    for i in range(2):
        for j in range(2):
            axd.text(j, i, f"{grid[i,j]:.0f}%", ha="center", va="center", fontsize=6.5,
                     color="white" if grid[i, j] > 25 else INK, weight="bold")
    axd.set_xticks([0, 1]); axd.set_xticklabels(["norm", "accel"], fontsize=5.2)
    axd.set_yticks([0, 1]); axd.set_yticklabels(["old", "young"], fontsize=5.2)
    axd.set_title("discordance (% non-home)", fontsize=6.0, color=INK, pad=2)
    # ROC
    axr = fig.add_axes([0.74, 0.055, 0.20, 0.11])
    sub = df.loc[clock.index]
    age_z = _z(sub["age_yrs"].astype(float).to_numpy()); fem = (sub["sex"] == "F").astype(float).to_numpy()
    asa_z = _z(pd.to_numeric(sub["asa"], errors="coerce").fillna(2).to_numpy()); aarz = _z(aar)
    for cm, color, ls, lab in [([age_z, fem, asa_z], SLATE, "-", "clinical"),
                               ([age_z, fem, asa_z, aarz], ORANGE, "--", "+ AAR")]:
        f = firth_logit(np.column_stack(cm), y); lp = np.column_stack([np.ones(len(y))] + cm) @ f["beta"]
        tpr = []; fpr = []
        for t in np.r_[np.inf, np.sort(lp)[::-1], -np.inf]:
            pred = lp >= t; tpr.append((pred & (y == 1)).sum() / (y == 1).sum()); fpr.append((pred & (y == 0)).sum() / (y == 0).sum())
        axr.plot(fpr, tpr, color=color, ls=ls, lw=1.4, label=f"{lab} {_auc(lp, y):.2f}")
    axr.plot([0, 1], [0, 1], color=GRAY, ls=":", lw=0.7); axr.set_aspect("equal")
    axr.set_xlabel("1 − spec.", fontsize=6); axr.set_ylabel("sens.", fontsize=6)
    axr.legend(frameon=False, fontsize=5.0, loc="lower right"); _mini(axr, "incremental value")
    _save(fig, "1.Figure_1_methods_overview")


# ======================================================================================
# Figure 2 — aging clock performance (restrained palette)
# ======================================================================================
def fig2_aging_clock(df, cfg):
    clock, meta = fit_clock(df, cfg)
    coef = clock_coefficients(df, cfg)
    age = clock["age_yrs"].to_numpy(); im = clock["imaging_age"].to_numpy(); aar = clock["aar"].to_numpy()
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7))
    ax = axs[0]
    ax.scatter(age, im, s=12, color=SLATE, alpha=0.55, edgecolor="none")
    lo, hi = 35, 90; ax.plot([lo, hi], [lo, hi], color=GRAY, ls="--", lw=1, label="Identity")
    b = np.polyfit(age, im, 1); ax.plot([lo, hi], [b[0]*lo+b[1], b[0]*hi+b[1]], color=ORANGE, lw=1.4, label="Fit")
    ax.set_xlabel("Chronological age, y"); ax.set_ylabel("Imaging (predicted) age, y")
    ax.set_title("A  Aging clock", loc="left", fontsize=9, weight="bold")
    ax.text(0.05, 0.92, f"R² = {meta['age_R2']:.2f}\nMAE = {meta['MAE']:.1f} y", transform=ax.transAxes, fontsize=7, va="top")
    ax.legend(frameon=False, fontsize=6.5, loc="lower right"); ax.spines[["top", "right"]].set_visible(False)
    ax = axs[1]
    ax.hist(aar, bins=22, color=SLATE, alpha=0.85, edgecolor="white", lw=0.4)
    ax.axvline(0, color=INK, lw=1, ls="--"); ax.set_xlabel("Age acceleration, y"); ax.set_ylabel("Patients")
    ax.set_title("B  Age-acceleration residual", loc="left", fontsize=9, weight="bold")
    ax.text(0.03, 0.92, "older-looking tissue", transform=ax.transAxes, fontsize=6.5, color=ORANGE, va="top")
    ax.spines[["top", "right"]].set_visible(False)
    ax = axs[2]
    c = coef.iloc[::-1]; ypos = np.arange(len(c))
    ax.barh(ypos, c["coef"], color=SLATE, alpha=0.85)
    ax.set_yticks(ypos); ax.set_yticklabels(c["label"], fontsize=6.5); ax.axvline(0, color=INK, lw=0.8)
    ax.set_xlabel("Ridge coefficient (higher = older)")
    ax.set_title("C  Tissue contributions", loc="left", fontsize=9, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); _save(fig, "2.Figure_2_aging_clock")


# ======================================================================================
# Figure 3 — primary association as a JAMA table-forest (two sections)
# ======================================================================================
def fig3_primary_association(df, cfg):
    a = pd.read_csv(RES / "aar_association.csv").set_index("model")
    s = pd.read_csv(RES / "clock_specs_sensitivity.csv").set_index("spec")
    rows = [{"section": "Age acceleration and non-home discharge (scanner-robust clock)"}]
    for m, lab in [("crude", "Crude"), ("adj_age", "Adjusted for chronological age"),
                   ("adj_age_sex_asa", "Adjusted for age, sex, ASA class")]:
        r = a.loc[m]; rows.append({"label": lab, "n": r["n"], "or": r["AAR_OR"], "lo": r["ci_lo"], "hi": r["ci_hi"], "p": r["p"]})
    rows.append({"section": "Clock specification (adjusted for age, sex, ASA class)"})
    for sp, lab in [("scanner_robust", "Scanner-robust (primary)"), ("intensity_ratio_only", "Intensity-ratio only"),
                    ("volume_only", "Volume only")]:
        r = s.loc[sp]; rows.append({"label": lab, "n": r["n"], "or": r["AAR_OR"], "lo": r["ci_lo"], "hi": r["ci_hi"], "p": r["p"]})
    table_forest(rows, "3.Figure_3_primary_association", xlim=(0.5, 3.5),
                 xticks=(0.5, 1.0, 2.0, 3.0), favors=("Favors home discharge", "Favors non-home discharge"),
                 analysis_header="Analysis (per 1-SD age acceleration)")


# ======================================================================================
# Figure 4 — robustness and incremental value (lambda selection, discordance, ROC)
# ======================================================================================
def fig4_robustness_and_value(df, cfg):
    feats = build_features(df); outcome = cfg["outcome"]
    cols = CLOCK_SPECS[PRIMARY_SPEC]
    d = pd.concat([feats[cols], df[["age_yrs", outcome]]], axis=1).dropna()
    X = d[cols].to_numpy(float); age = d["age_yrs"].to_numpy(float)
    r2s = [1 - np.sum((age - _ridge_cv_age(X, age, lam, cfg["seed"])) ** 2) / np.sum((age - age.mean()) ** 2) for lam in LAMBDA_GRID]
    clock, _ = fit_clock(df, cfg); idx = clock.index; sub = df.loc[idx]
    y = clock[outcome].astype(int).to_numpy(); aar_raw = clock["aar"].to_numpy(); agec = clock["age_yrs"].to_numpy()

    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7))
    ax = axs[0]
    ax.plot(LAMBDA_GRID, r2s, "o-", color=SLATE, ms=4); ax.set_xscale("log")
    best = LAMBDA_GRID[int(np.argmax(r2s))]; ax.axvline(best, color=ORANGE, ls="--", lw=1)
    ax.text(best, min(r2s), f" λ*={best:g}", fontsize=6.5, color=ORANGE)
    ax.set_xlabel("Ridge penalty λ"); ax.set_ylabel("Age-prediction R² (CV)")
    ax.set_title("A  λ chosen by age fit", loc="left", fontsize=9, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    # discordance
    ax = axs[1]
    amed = np.median(agec); grid = np.zeros((2, 2))
    for i, old in enumerate([1, 0]):
        for j, acc in enumerate([0, 1]):
            m = ((agec >= amed) == bool(old)) & ((aar_raw > 0) == bool(acc))
            grid[i, j] = 100 * y[m].mean() if m.sum() else np.nan
    im = ax.imshow(grid, cmap="Oranges", vmin=0, vmax=45, aspect="auto")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, f"{grid[i,j]:.0f}%", ha="center", va="center", fontsize=9,
                    color="white" if grid[i, j] > 25 else INK, weight="bold")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Normal", "Accelerated"], fontsize=7)
    ax.set_yticks([0, 1]); ax.set_yticklabels(["Older", "Younger"], fontsize=7)
    ax.set_xlabel("Imaging age", fontsize=8); ax.set_ylabel("Chronological age", fontsize=8)
    ax.set_title("B  Discordance (% non-home)", loc="left", fontsize=9, weight="bold")
    # ROC
    ax = axs[2]
    age_z = _z(sub["age_yrs"].astype(float).to_numpy()); fem = (sub["sex"] == "F").astype(float).to_numpy()
    asa_z = _z(pd.to_numeric(sub["asa"], errors="coerce").fillna(2).to_numpy()); aar = _z(aar_raw)
    for cols_m, color, ls, lab in [([age_z, fem, asa_z], SLATE, "-", "Clinical"),
                                    ([age_z, fem, asa_z, aar], ORANGE, "--", "Clinical + AAR")]:
        f = firth_logit(np.column_stack(cols_m), y); lp = np.column_stack([np.ones(len(y))] + cols_m) @ f["beta"]
        tpr = []; fpr = []
        for t in np.r_[np.inf, np.sort(lp)[::-1], -np.inf]:
            pred = lp >= t; tpr.append((pred & (y == 1)).sum() / (y == 1).sum()); fpr.append((pred & (y == 0)).sum() / (y == 0).sum())
        ax.plot(fpr, tpr, color=color, ls=ls, lw=1.7, label=f"{lab} (apparent AUC {_auc(lp, y):.2f})")
    ax.plot([0, 1], [0, 1], color=GRAY, ls=":", lw=0.8)
    ax.set_xlabel("1 − specificity"); ax.set_ylabel("Sensitivity")
    ax.set_title("C  Incremental value", loc="left", fontsize=9, weight="bold")
    ax.legend(frameon=False, fontsize=6, loc="lower right"); ax.set_aspect("equal"); ax.spines[["top", "right"]].set_visible(False)
    iauc = json.load(open(RES / "incremental_auc.json"))
    ax.text(0.5, -0.36, f"Cross-validated ΔAUC {iauc['delta_auc']:+.3f} (95% CI {iauc['ci_lo']:+.3f} to {iauc['ci_hi']:+.3f})",
            transform=ax.transAxes, ha="center", fontsize=5.6, color=GRAY)
    fig.tight_layout(); _save(fig, "4.Figure_4_robustness_and_value")


# ======================================================================================
# Supplemental figures
# ======================================================================================
def figS1_strobe(df, cfg):
    fig = plt.figure(figsize=(5.0, 5.6)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.5, 0.97, "Participant flow (STROBE)", ha="center", fontsize=10, weight="bold", color=INK)
    n_total = 206; n_out = int(df["rehab"].notna().sum())
    clock, _ = fit_clock(df, cfg); n_clock = len(clock); ev = int(clock["rehab"].sum())
    _box(ax, 0.15, 0.80, 0.70, 0.11, f"Adults undergoing lumbar spine surgery\nwith preoperative lumbar MRI (N = {n_total})", fs=8)
    _box(ax, 0.15, 0.55, 0.70, 0.11, f"Analytic cohort with discharge outcome,\nsegmentable iliopsoas, plausible age (n = {n_out})", fs=8)
    _box(ax, 0.15, 0.30, 0.70, 0.12, f"Complete multi-tissue imaging for the\naging clock (n = {n_clock})", fs=8, ec=SLATE)
    _box(ax, 0.18, 0.06, 0.30, 0.13, f"Home\n(n = {n_clock-ev})", fc=LIGHT, fs=8)
    _box(ax, 0.52, 0.06, 0.30, 0.13, f"Non-home\n(n = {ev}, {ev/n_clock*100:.0f}%)", fc="white", ec=ORANGE, fs=8)
    _arrow(ax, 0.5, 0.80, 0.5, 0.66); _arrow(ax, 0.5, 0.55, 0.5, 0.42); _arrow(ax, 0.5, 0.30, 0.33, 0.19); _arrow(ax, 0.5, 0.30, 0.67, 0.19)
    _box(ax, 0.60, 0.685, 0.34, 0.075, f"Excluded ({n_total-n_out})", fc="white", ec=GRAY, fs=6.5)
    _box(ax, 0.60, 0.435, 0.34, 0.075, f"Excluded: incomplete\nmulti-tissue imaging ({n_out-n_clock})", fc="white", ec=GRAY, fs=6.5)
    _save(fig, "S1.Figure_S1_STROBE_flow")


def figS2_feature_correlation(df, cfg):
    feats = build_features(df)[CLOCK_SPECS[PRIMARY_SPEC]].dropna()
    from .aging_clock import FEATURE_LABELS
    C = np.corrcoef(((feats - feats.mean()) / feats.std()).to_numpy().T)
    fig, ax = plt.subplots(figsize=(4.6, 4.0))
    im = ax.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1)
    labs = [FEATURE_LABELS[c] for c in feats.columns]
    ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, rotation=45, ha="right", fontsize=6)
    ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs, fontsize=6)
    for i in range(len(labs)):
        for j in range(len(labs)):
            ax.text(j, i, f"{C[i,j]:.1f}", ha="center", va="center", fontsize=5,
                    color="white" if abs(C[i, j]) > 0.6 else INK)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Pearson r")
    ax.set_title("Imaging-feature correlations", fontsize=9, weight="bold", color=INK)
    fig.tight_layout(); _save(fig, "S2.Figure_S2_feature_correlation")


def figS3_naive_approaches(df, cfg):
    """Cautionary table-forest: single-muscle, age-gap suppression artifact, residual."""
    outcome = cfg["outcome"]; sub = df.copy(); sub["female"] = (sub["sex"] == "F").astype(float)
    age = _z(pd.to_numeric(sub["age_yrs"], errors="coerce").to_numpy())
    asa = _z(pd.to_numeric(sub["asa"], errors="coerce").fillna(2).to_numpy())
    fem = sub["female"].to_numpy(); y = sub[outcome].astype(int).to_numpy()

    def adj(x):
        m = ~np.isnan(x); f = firth_logit(np.column_stack([_z(x[m]), age[m], fem[m], asa[m]]), y[m])
        return f["or"][1], f["ci_lo"][1], f["ci_hi"][1], f["p"][1], int(m.sum())
    ilio = pd.to_numeric(sub["iliopsoas__vol_norm_vert"], errors="coerce").to_numpy()
    o1, l1, h1, p1, n1 = adj(ilio)
    clock, _ = fit_clock(df, cfg); gap = (clock["imaging_age"] - clock["age_yrs"]).to_numpy()
    yg = clock[outcome].astype(int).to_numpy(); ag = _z(clock["age_yrs"].to_numpy())
    fg = firth_logit(np.column_stack([_z(gap), ag]), yg)
    aa = pd.read_csv(RES / "aar_association.csv"); rr = aa[aa["model"] == "adj_age_sex_asa"].iloc[0]
    rows = [
        {"label": "Single muscle (iliopsoas volume)", "n": n1, "or": o1, "lo": l1, "hi": h1, "p": p1},
        {"label": "Naive age-gap, adjusted for age", "dagger": True, "n": len(yg), "or": fg["or"][1], "lo": fg["ci_lo"][1], "hi": fg["ci_hi"][1], "p": fg["p"][1]},
        {"label": "Age acceleration (residual method)", "n": int(rr["n"]), "or": rr["AAR_OR"], "lo": rr["ci_lo"], "hi": rr["ci_hi"], "p": rr["p"]},
    ]
    table_forest(rows, "S3.Figure_S3_naive_approaches", xlim=(0.5, 5.0), xticks=(0.5, 1.0, 2.0, 4.0),
                 favors=("Favors home discharge", "Favors non-home discharge"),
                 title="Why naive operationalizations mislead",
                 note="† The age-gap adjusted for chronological age is a suppression artifact: the gap already contains age.")


# ======================================================================================
# Tables
# ======================================================================================
def _fmt(o, lo, hi, p):
    return f"{o:.2f} ({lo:.2f}-{hi:.2f})", f"{p:.3f}"


def _smd(a, b, binary):
    a, b = pd.to_numeric(a, errors="coerce").dropna(), pd.to_numeric(b, errors="coerce").dropna()
    if binary:
        p1, p2 = a.mean(), b.mean(); s = np.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
    else:
        s = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    return (a.mean() - b.mean()) / s if s > 0 else 0.0


def tables(df, cfg):
    RES.mkdir(exist_ok=True); TAB = FIG.parent / "tables"; TAB.mkdir(exist_ok=True)
    outcome = cfg["outcome"]
    df = df.copy(); df["female"] = (df["sex"] == "F").astype(float)
    clock, _ = fit_clock(df, cfg); ci = clock.index
    # Table 1 — clock cohort by discharge, with imaging age, AAR, and standardized differences
    t1df = df.loc[ci].copy(); t1df["imaging_age"] = clock["imaging_age"].values; t1df["aar"] = clock["aar"].values
    g0, g1 = t1df[t1df[outcome] == 0], t1df[t1df[outcome] == 1]
    rows = []
    for lab, col, binary in [("Age, y (median)", "age_yrs", False), ("Female, %", "female", True),
                             ("ASA class (median)", "asa", False), ("Hypertension, %", "htn", True),
                             ("Diabetes, %", "diabetes", True), ("Fusion, %", "fusion", True),
                             ("No. operated levels (median)", "num_level", False),
                             ("Imaging age, y (median)", "imaging_age", False),
                             ("Age acceleration, y (median)", "aar", False)]:
        a0, a1 = pd.to_numeric(g0[col], errors="coerce"), pd.to_numeric(g1[col], errors="coerce")
        home = f"{100*a0.mean():.0f}" if binary else f"{a0.median():.1f}"
        nonh = f"{100*a1.mean():.0f}" if binary else f"{a1.median():.1f}"
        rows.append({"Characteristic": lab, "Home": home, "Non-home": nonh, "SMD": f"{_smd(a0, a1, binary):+.2f}"})
    pd.DataFrame(rows).to_csv(TAB / "1.Table_1_cohort.csv", index=False)
    # Table S3 — attrition (included vs excluded for incomplete imaging)
    inc, exc = df.loc[ci], df.drop(index=ci)
    arows = []
    for lab, col, binary in [("Age, y (mean)", "age_yrs", False), ("Female, %", "female", True),
                             ("ASA class (mean)", "asa", False), ("Hypertension, %", "htn", True),
                             ("Diabetes, %", "diabetes", True), ("Non-home discharge, %", outcome, True)]:
        ai, ae = pd.to_numeric(inc[col], errors="coerce"), pd.to_numeric(exc[col], errors="coerce")
        vi = f"{100*ai.mean():.0f}" if binary else f"{ai.mean():.1f}"
        ve = f"{100*ae.mean():.0f}" if binary else f"{ae.mean():.1f}"
        arows.append({"Characteristic": lab, "Included": vi, "Excluded": ve, "SMD": f"{_smd(ai, ae, binary):+.2f}"})
    pd.DataFrame(arows).to_csv(TAB / "S3.Table_S3_attrition.csv", index=False)
    def _est(r): return f"{r['AAR_OR']:.2f} ({r['ci_lo']:.2f}-{r['ci_hi']:.2f})", _pfmt(r["p"])
    a = pd.read_csv(RES / "aar_association.csv"); s = pd.read_csv(RES / "clock_specs_sensitivity.csv")
    ac = pd.read_csv(RES / "adjusted_model_coefficients.csv"); bmi = pd.read_csv(RES / "bmi_sensitivity.csv")
    t2 = []
    for _, r in a.iterrows():
        est, p = _est(r)
        t2.append({"Section": "Age acceleration and non-home discharge (scanner-robust clock)",
                   "Analysis": {"crude": "Crude", "adj_age": "Adjusted for age",
                                "adj_age_sex_asa": "Adjusted for age, sex, ASA"}[r["model"]],
                   "OR (95% CI)": est, "P": p, "n": int(r["n"])})
    for _, r in ac.iterrows():
        t2.append({"Section": "Full adjusted model (each term)", "Analysis": r["term"],
                   "OR (95% CI)": f"{r['OR']:.2f} ({r['ci_lo']:.2f}-{r['ci_hi']:.2f})", "P": _pfmt(r["p"]), "n": 192})
    for _, r in s.iterrows():
        est, p = _est(r)
        t2.append({"Section": "Clock specification (adjusted for age, sex, ASA)",
                   "Analysis": {"scanner_robust": "Scanner-robust", "intensity_ratio_only": "Intensity-ratio only",
                                "volume_only": "Volume only"}[r["spec"]], "OR (95% CI)": est, "P": p, "n": int(r["n"])})
    for _, r in bmi.iterrows():
        est, p = _est(r)
        t2.append({"Section": "Sensitivity: BMI adjustment", "Analysis": r["model"],
                   "OR (95% CI)": est, "P": p, "n": int(r["n"])})
    for _, r in pd.read_csv(RES / "surgical_sensitivity.csv").iterrows():
        est, p = _est(r)
        t2.append({"Section": "Sensitivity: surgical invasiveness", "Analysis": r["model"],
                   "OR (95% CI)": est, "P": p, "n": int(r["n"])})
    pd.DataFrame(t2).to_csv(TAB / "2.Table_2_primary_association.csv", index=False)
    met = pd.read_csv(RES / "clock_metrics.csv"); iauc = json.load(open(RES / "incremental_auc.json"))
    t3 = [{"Metric": "Clock age R² (out-of-fold)", "Value": f"{met.iloc[0]['age_R2']:.2f}"},
          {"Metric": "Clock MAE, y (out-of-fold)", "Value": f"{met.iloc[0]['MAE']:.1f}"},
          {"Metric": "AUC clinical (age, sex, ASA), CV", "Value": f"{iauc['auc_clinical_cv']:.3f}"},
          {"Metric": "AUC clinical + age acceleration, CV", "Value": f"{iauc['auc_plus_aar_cv']:.3f}"},
          {"Metric": "ΔAUC (cross-validated)", "Value": f"{iauc['delta_auc']:+.3f} (95% CI {iauc['ci_lo']:+.3f} to {iauc['ci_hi']:+.3f})"}]
    pd.DataFrame(t3).to_csv(TAB / "3.Table_3_clock_and_value.csv", index=False)
    feats = build_features(df); from .aging_clock import FEATURE_LABELS
    n = len(df)
    st1 = [{"Feature": FEATURE_LABELS[c], "Spec": ", ".join(k for k, v in CLOCK_SPECS.items() if c in v),
            "n observed": int(feats[c].notna().sum()), "% missing": round(100*(n-feats[c].notna().sum())/n, 1)}
           for c in feats.columns]
    pd.DataFrame(st1).to_csv(TAB / "S1.Table_S1_features.csv", index=False)
    cols = CLOCK_SPECS[PRIMARY_SPEC]
    d = pd.concat([feats[cols], df[["age_yrs", outcome]]], axis=1).dropna()
    X = d[cols].to_numpy(float); agev = d["age_yrs"].to_numpy(float)
    yv = d[outcome].astype(int).to_numpy(); st2 = []
    for seed in [1, 7, 42, 2024]:
        for lam in [3.0, 10.0, 30.0]:
            pred = _ridge_cv_age(X, agev, lam, seed); b = np.polyfit(agev, pred, 1); aarv = pred - (b[0]*agev+b[1])
            f = firth_logit(np.column_stack([_z(aarv), _z(agev)]), yv)
            st2.append({"seed": seed, "lambda": lam, "AAR_OR": round(f["or"][1], 2), "P": round(f["p"][1], 3)})
    pd.DataFrame(st2).to_csv(TAB / "S2.Table_S2_sensitivity.csv", index=False)


def make_all():
    from .data_loading import load_cohort, load_config
    cfg = load_config("config.yaml"); df = load_cohort(cfg)
    # fig1_methods_overview is intentionally NOT called here: the committed Figure 1
    # is hand-finalized in Inkscape (it embeds the de-identified multi-tissue MRI
    # panel). Run fig1_methods_overview(df, cfg) explicitly to regenerate the code
    # baseline (this overwrites the hand-edited figure with the MRI placeholder).
    fig2_aging_clock(df, cfg)
    fig3_primary_association(df, cfg)
    fig4_robustness_and_value(df, cfg)
    # figS1_strobe is intentionally NOT called here: the committed STROBE flow
    # (S1) is hand-finalized in Inkscape. Run figS1_strobe(df, cfg) explicitly to
    # regenerate the code baseline (this overwrites the hand-edited S1 SVG).
    figS2_feature_correlation(df, cfg)
    figS3_naive_approaches(df, cfg)
    tables(df, cfg)
    print("figures + tables written")


if __name__ == "__main__":
    make_all()

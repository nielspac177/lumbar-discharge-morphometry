"""Figures and tables for the MRI aging-clock manuscript.

All figures render to vector SVG (for Inkscape editing) plus PNG and PDF, using the
JAMA (ggsci pal_jama) palette. Files follow the manuscript numbering convention
``<order>.Figure_<n>_<slug>`` / ``<order>.Table_<n>_<slug>``. Every plotted value
traces to a committed ``results/`` file or is computed here from the frozen cohort.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.ticker import NullLocator

from .aging_clock import (CLOCK_SPECS, LAMBDA_GRID, PRIMARY_SPEC, build_features,
                          clock_coefficients, fit_clock, _ridge_cv_age, _z)
from .firth import firth_logit

JAMA = {"slate": "#374E55", "orange": "#DF8F44", "blue": "#00A1D5", "red": "#B24745",
        "green": "#79AF97", "purple": "#6A6599", "gray": "#80796B"}
rcParams.update({
    # DejaVu Sans has full glyph coverage so raster previews never fail; svg.fonttype
    # 'none' keeps text editable in Inkscape, where it can be restyled to Helvetica/Arial.
    "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": 8, "axes.linewidth": 0.8, "savefig.dpi": 600,
    "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
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


def _box(ax, x, y, w, h, text, fc="white", ec=JAMA["slate"], fs=7.5, lw=1.2, weight="normal", tc=None):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.012",
                                fc=fc, ec=ec, lw=lw, mutation_aspect=1))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            weight=weight, color=tc or JAMA["slate"], zorder=5)


def _arrow(ax, x0, y0, x1, y1, color=None):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=11,
                                 lw=1.3, color=color or JAMA["gray"], shrinkA=2, shrinkB=2))


# ======================================================================================
# Figure 1 — Methods / pipeline overview (schematic)
# ======================================================================================
def fig1_methods_overview():
    fig = plt.figure(figsize=(7.2, 4.6))
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.5, 0.965, "Multi-tissue MRI aging clock and age acceleration",
            ha="center", fontsize=11, weight="bold", color=JAMA["slate"])
    # Stage A - segmentation
    _box(ax, 0.02, 0.60, 0.20, 0.28,
         "A. Preoperative lumbar MRI\n\nAutomated segmentation\n(3D Slicer / TotalSegmentator)\n\n• Iliopsoas\n• Deep paraspinal\n• Vertebra • Disc • Cord",
         fc="#F4F7F8", fs=7)
    # Stage B - clock
    _box(ax, 0.27, 0.60, 0.22, 0.28,
         "B. Multi-tissue signature\n\nScanner-robust features:\n- size-normalized volumes\n- T2 intensity ratios\n  (tissue / vertebra)\n\nyields ridge aging clock\n(10-fold CV; lambda by age R2)",
         fc="#F4F7F8", fs=7)
    _box(ax, 0.54, 0.66, 0.20, 0.16, "Imaging age\n(predicted from tissue)", fc="white",
         ec=JAMA["blue"], fs=7.5, weight="bold")
    # Stage C - AAR
    _box(ax, 0.78, 0.60, 0.20, 0.28,
         "C. Age acceleration\n\nresidual of\nimaging age ~ chronological age\n\nAAR independent of\nchronological age\n(no suppression)",
         fc="#FBF3EA", ec=JAMA["orange"], fs=7)
    _arrow(ax, 0.22, 0.74, 0.27, 0.74)
    _arrow(ax, 0.49, 0.74, 0.54, 0.74)
    _arrow(ax, 0.74, 0.74, 0.78, 0.74)
    # Stage D - association (bottom band)
    _box(ax, 0.06, 0.10, 0.40, 0.34,
         "D. Association analysis\n\nNon-home discharge  ~  age acceleration\nadjusted for chronological age, sex, ASA class\n\nFirth penalized logistic regression\n(small-sample robust)\n\nSensitivity: clock specification, lambda, seed;\nvolume-only and intensity-ratio-only clocks",
         fc="white", ec=JAMA["slate"], fs=7)
    # small illustrative AAR scatter inset
    axi = fig.add_axes([0.55, 0.09, 0.40, 0.34])
    rng = np.random.default_rng(0); a = rng.uniform(40, 85, 60); im = 0.42 * (a - 62) + 62 + rng.normal(0, 8, 60)
    axi.scatter(a, im, s=10, color=JAMA["blue"], alpha=0.6, edgecolor="none")
    axi.plot([40, 85], [40, 85], color=JAMA["gray"], ls="--", lw=1)
    for i in [5, 20, 40]:
        axi.plot([a[i], a[i]], [0.42 * (a[i] - 62) + 62, im[i]], color=JAMA["orange"], lw=0.8)
    axi.set_xlabel("Chronological age, y", fontsize=7); axi.set_ylabel("Imaging age, y", fontsize=7)
    axi.set_title("Age acceleration = vertical residual", fontsize=7.5, color=JAMA["slate"])
    axi.tick_params(labelsize=6); axi.spines[["top", "right"]].set_visible(False)
    _arrow(ax, 0.46, 0.27, 0.55, 0.27, color=JAMA["orange"])
    _save(fig, "1.Figure_1_methods_overview")


# ======================================================================================
# Figure 2 — Aging clock performance
# ======================================================================================
def fig2_aging_clock(df, cfg):
    clock, meta = fit_clock(df, cfg)
    coef = clock_coefficients(df, cfg)
    age = clock["age_yrs"].to_numpy(); im = clock["imaging_age"].to_numpy(); aar = clock["aar"].to_numpy()
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7))
    # A scatter
    ax = axs[0]
    ax.scatter(age, im, s=12, color=JAMA["blue"], alpha=0.6, edgecolor="none")
    lo, hi = 35, 90; ax.plot([lo, hi], [lo, hi], color=JAMA["gray"], ls="--", lw=1, label="Identity")
    b = np.polyfit(age, im, 1); ax.plot([lo, hi], [b[0] * lo + b[1], b[0] * hi + b[1]], color=JAMA["red"], lw=1.3, label="Fit")
    ax.set_xlabel("Chronological age, y"); ax.set_ylabel("Imaging (predicted) age, y")
    ax.set_title("A  Aging clock", loc="left", fontsize=9, weight="bold")
    ax.text(0.05, 0.92, f"R² = {meta['age_R2']:.2f}\nMAE = {meta['MAE']:.1f} y",
            transform=ax.transAxes, fontsize=7, va="top")
    ax.legend(frameon=False, fontsize=6.5, loc="lower right"); ax.spines[["top", "right"]].set_visible(False)
    # B AAR distribution
    ax = axs[1]
    ax.hist(aar, bins=22, color=JAMA["orange"], alpha=0.85, edgecolor="white", lw=0.4)
    ax.axvline(0, color=JAMA["slate"], lw=1, ls="--")
    ax.set_xlabel("Age acceleration, y"); ax.set_ylabel("Patients")
    ax.set_title("B  Age-acceleration residual", loc="left", fontsize=9, weight="bold")
    ax.text(0.03, 0.92, "older-looking tissue", transform=ax.transAxes, fontsize=6.5, color=JAMA["red"], va="top")
    ax.spines[["top", "right"]].set_visible(False)
    # C tissue contributions
    ax = axs[2]
    c = coef.iloc[::-1]; ypos = np.arange(len(c))
    cols = [JAMA["blue"] if v > 0 else JAMA["red"] for v in c["coef"]]
    ax.barh(ypos, c["coef"], color=cols, alpha=0.85)
    ax.set_yticks(ypos); ax.set_yticklabels(c["label"], fontsize=6.5)
    ax.axvline(0, color=JAMA["slate"], lw=0.8)
    ax.set_xlabel("Ridge coefficient (higher = older)")
    ax.set_title("C  Tissue contributions", loc="left", fontsize=9, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    _save(fig, "2.Figure_2_aging_clock")


# ======================================================================================
# Figure 3 — Primary association
# ======================================================================================
def _forest(ax, labels, ors, los, his, ps, title, ref=1.0, xlim=(0.5, 4.0)):
    y = np.arange(len(labels))[::-1]
    ax.axvline(ref, color=JAMA["gray"], lw=0.8, ls="--")
    for yi, o, lo, hi, p in zip(y, ors, los, his, ps):
        sig = p < 0.05
        ax.plot([lo, hi], [yi, yi], color=JAMA["slate"], lw=1.3)
        ax.plot(o, yi, "s", ms=6, color=JAMA["blue"] if sig else JAMA["gray"])
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7)
    ax.set_xscale("log"); ax.set_xlim(*xlim); ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticks([0.5, 1, 2, 4]); ax.set_xticklabels(["0.5", "1", "2", "4"])
    ax.set_xlabel("Odds ratio per SD (95% CI)")
    ax.set_title(title, loc="left", fontsize=9, weight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False); ax.tick_params(left=False)


def fig3_primary_association(df, cfg):
    assoc = pd.read_csv(RES / "aar_association.csv")
    sens = pd.read_csv(RES / "clock_specs_sensitivity.csv")
    clock, _ = fit_clock(df, cfg)
    fig = plt.figure(figsize=(7.4, 5.4))
    # A forest: adjustment ladder (primary spec)
    axA = fig.add_axes([0.30, 0.58, 0.66, 0.36])
    m = {"crude": "Crude", "adj_age": "+ chronological age", "adj_age_sex_asa": "+ age, sex, ASA"}
    a = assoc.set_index("model").loc[["crude", "adj_age", "adj_age_sex_asa"]]
    _forest(axA, [m[i] for i in a.index], a["AAR_OR"], a["ci_lo"], a["ci_hi"], a["p"],
            "A  Age acceleration and non-home discharge (adjustment)")
    for yi, (_, r) in zip(np.arange(len(a))[::-1], a.iterrows()):
        axA.text(4.15, yi, f"{r['AAR_OR']:.2f} ({r['ci_lo']:.2f}–{r['ci_hi']:.2f})  P={r['p']:.3f}",
                 fontsize=6.3, va="center")
    # B forest: clock specification robustness
    axB = fig.add_axes([0.30, 0.10, 0.40, 0.30])
    slabel = {"scanner_robust": "Scanner-robust\n(primary)", "intensity_ratio_only": "Intensity-ratio\nonly",
              "volume_only": "Volume only"}
    s = sens.set_index("spec").loc[["scanner_robust", "intensity_ratio_only", "volume_only"]]
    _forest(axB, [slabel[i] for i in s.index], s["AAR_OR"], s["ci_lo"], s["ci_hi"], s["p"],
            "B  Clock specification", xlim=(0.5, 4.0))
    # C discordance 2x2
    axC = fig.add_axes([0.80, 0.10, 0.17, 0.30])
    age = clock["age_yrs"].to_numpy(); aar = clock["aar"].to_numpy(); y = clock[cfg["outcome"]].to_numpy()
    amed = np.median(age)
    grid = np.zeros((2, 2))
    for i, old in enumerate([1, 0]):
        for j, acc in enumerate([0, 1]):
            mask = ((age >= amed) == bool(old)) & ((aar > 0) == bool(acc))
            grid[i, j] = 100 * y[mask].mean() if mask.sum() else np.nan
    im = axC.imshow(grid, cmap="Oranges", vmin=0, vmax=45, aspect="auto")
    for i in range(2):
        for j in range(2):
            axC.text(j, i, f"{grid[i,j]:.0f}%", ha="center", va="center", fontsize=8,
                     color="white" if grid[i, j] > 25 else JAMA["slate"], weight="bold")
    axC.set_xticks([0, 1]); axC.set_xticklabels(["Normal", "Accelerated"], fontsize=6.5)
    axC.set_yticks([0, 1]); axC.set_yticklabels(["Older", "Younger"], fontsize=6.5)
    axC.set_xlabel("Imaging age", fontsize=7); axC.set_ylabel("Chronological", fontsize=7)
    axC.set_title("C  Discordance\n(% non-home)", loc="left", fontsize=8, weight="bold")
    _save(fig, "3.Figure_3_primary_association")


# ======================================================================================
# Figure 4 — Robustness and incremental value
# ======================================================================================
def fig4_robustness_and_value(df, cfg):
    feats = build_features(df); outcome = cfg["outcome"]
    # A lambda vs age-R2 (primary spec)
    cols = CLOCK_SPECS[PRIMARY_SPEC]
    d = pd.concat([feats[cols], df[["age_yrs", outcome]]], axis=1).dropna()
    X = d[cols].to_numpy(float); age = d["age_yrs"].to_numpy(float); Xs = (X - X.mean(0)) / X.std(0)
    r2s = []
    for lam in LAMBDA_GRID:
        p = _ridge_cv_age(Xs, age, lam, cfg["seed"])
        r2s.append(1 - np.sum((age - p) ** 2) / np.sum((age - age.mean()) ** 2))
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7))
    ax = axs[0]
    ax.plot(LAMBDA_GRID, r2s, "o-", color=JAMA["blue"]); ax.set_xscale("log")
    best = LAMBDA_GRID[int(np.argmax(r2s))]
    ax.axvline(best, color=JAMA["orange"], ls="--", lw=1); ax.text(best, min(r2s), f" λ*={best:g}", fontsize=6.5, color=JAMA["orange"])
    ax.set_xlabel("Ridge penalty λ"); ax.set_ylabel("Age-prediction R² (CV)")
    ax.set_title("A  λ chosen by age fit", loc="left", fontsize=9, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    # B spec bars
    ax = axs[1]
    sens = pd.read_csv(RES / "clock_specs_sensitivity.csv").set_index("spec").loc[
        ["scanner_robust", "intensity_ratio_only", "volume_only"]]
    labs = ["Scanner-\nrobust", "Intensity-\nratio", "Volume\nonly"]
    cols_b = [JAMA["blue"] if p < 0.05 else JAMA["gray"] for p in sens["p"]]
    ax.bar(range(3), sens["AAR_OR"], color=cols_b, alpha=0.85)
    ax.errorbar(range(3), sens["AAR_OR"], yerr=[sens["AAR_OR"] - sens["ci_lo"], sens["ci_hi"] - sens["AAR_OR"]],
                fmt="none", ecolor=JAMA["slate"], lw=1, capsize=2)
    ax.axhline(1, color=JAMA["gray"], ls="--", lw=0.8)
    ax.set_xticks(range(3)); ax.set_xticklabels(labs, fontsize=6.5); ax.set_ylabel("AAR OR (adj)")
    ax.set_title("B  Clock specification", loc="left", fontsize=9, weight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    # C ROC clinical vs +AAR
    ax = axs[2]
    clock, _ = fit_clock(df, cfg); idx = clock.index; sub = df.loc[idx]
    y = clock[outcome].astype(int).to_numpy()
    age_z = _z(sub["age_yrs"].astype(float).to_numpy()); fem = (sub["sex"] == "F").astype(float).to_numpy()
    asa_z = _z(pd.to_numeric(sub["asa"], errors="coerce").fillna(2).to_numpy()); aar = _z(clock["aar"].to_numpy())
    for cols_m, color, lab in [([age_z, fem, asa_z], JAMA["orange"], "Clinical"),
                                ([age_z, fem, asa_z, aar], JAMA["blue"], "Clinical + AAR")]:
        f = firth_logit(np.column_stack(cols_m), y)
        lp = np.column_stack([np.ones(len(y))] + cols_m) @ f["beta"]
        thr = np.unique(lp); tpr = []; fpr = []
        for t in np.r_[-np.inf, np.sort(lp), np.inf][::-1]:
            pred = lp >= t; tpr.append((pred & (y == 1)).sum() / (y == 1).sum()); fpr.append((pred & (y == 0)).sum() / (y == 0).sum())
        ax.plot(fpr, tpr, color=color, lw=1.6, label=f"{lab} (AUC {_auc(lp, y):.2f})")
    ax.plot([0, 1], [0, 1], color=JAMA["gray"], ls=":", lw=0.8)
    ax.set_xlabel("1 − specificity"); ax.set_ylabel("Sensitivity")
    ax.set_title("C  Incremental value", loc="left", fontsize=9, weight="bold")
    ax.legend(frameon=False, fontsize=6, loc="lower right"); ax.set_aspect("equal"); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    _save(fig, "4.Figure_4_robustness_and_value")


# ======================================================================================
# Supplemental figures
# ======================================================================================
def figS1_strobe(df, cfg):
    fig = plt.figure(figsize=(5.0, 5.6)); ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.5, 0.97, "Participant flow (STROBE)", ha="center", fontsize=10, weight="bold", color=JAMA["slate"])
    n_total = 206; n_out = int(df["rehab"].notna().sum())
    clock, _ = fit_clock(df, cfg); n_clock = len(clock); ev = int(clock["rehab"].sum())
    _box(ax, 0.15, 0.80, 0.70, 0.11, f"Adults undergoing lumbar spine surgery\nwith preoperative lumbar MRI (N = {n_total})", fs=8)
    _box(ax, 0.15, 0.55, 0.70, 0.11, f"Analytic cohort with discharge outcome and\nsegmentable iliopsoas (n = {n_out})", fs=8)
    _box(ax, 0.15, 0.30, 0.70, 0.12, f"Complete multi-tissue imaging for the\naging clock (n = {n_clock})", fs=8, ec=JAMA["blue"])
    _box(ax, 0.18, 0.06, 0.30, 0.13, f"Home\n(n = {n_clock-ev})", fc="#F4F7F8", fs=8)
    _box(ax, 0.52, 0.06, 0.30, 0.13, f"Non-home\n(n = {ev}, {ev/n_clock*100:.0f}%)", fc="#FBF3EA", ec=JAMA["orange"], fs=8)
    _arrow(ax, 0.5, 0.80, 0.5, 0.66); _arrow(ax, 0.5, 0.55, 0.5, 0.42); _arrow(ax, 0.5, 0.30, 0.33, 0.19); _arrow(ax, 0.5, 0.30, 0.67, 0.19)
    _box(ax, 0.60, 0.685, 0.34, 0.075, f"Excluded: outcome missing\nor no iliopsoas ({n_total-n_out})", fc="white", ec=JAMA["gray"], fs=6.5)
    _box(ax, 0.60, 0.435, 0.34, 0.075, f"Excluded: incomplete\nmulti-tissue imaging ({n_out-n_clock})", fc="white", ec=JAMA["gray"], fs=6.5)
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
                    color="white" if abs(C[i, j]) > 0.6 else JAMA["slate"])
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Pearson r")
    ax.set_title("Imaging-feature correlations", fontsize=9, weight="bold", color=JAMA["slate"])
    fig.tight_layout(); _save(fig, "S2.Figure_S2_feature_correlation")


def figS3_naive_approaches(df, cfg):
    """Cautionary: single-muscle, threshold, and the age-gap suppression artifact."""
    outcome = cfg["outcome"]; sub = df.copy()
    sub["female"] = (sub["sex"] == "F").astype(float)
    age = _z(pd.to_numeric(sub["age_yrs"], errors="coerce").to_numpy())
    asa = _z(pd.to_numeric(sub["asa"], errors="coerce").fillna(2).to_numpy())
    fem = sub["female"].to_numpy(); y = sub[outcome].astype(int).to_numpy()
    def adj_or(x):
        m = ~np.isnan(x)
        f = firth_logit(np.column_stack([_z(x[m]), age[m], fem[m], asa[m]]), y[m])
        return f["or"][1], f["ci_lo"][1], f["ci_hi"][1], f["p"][1]
    ilio = pd.to_numeric(sub["iliopsoas__vol_norm_vert"], errors="coerce").to_numpy()
    # naive age-gap (raw imaging_age - age, adjusted for age -> suppression)
    clock, _ = fit_clock(df, cfg); idx = clock.index
    gap = (clock["imaging_age"] - clock["age_yrs"]).to_numpy()
    yg = clock[outcome].astype(int).to_numpy(); ag = _z(clock["age_yrs"].to_numpy())
    fgap = firth_logit(np.column_stack([_z(gap), ag]), yg)  # gap + age = suppression
    items = [("Single muscle\n(iliopsoas volume)", *adj_or(ilio)),
             ("Naive age-gap\n(adj. for age)†", fgap["or"][1], fgap["ci_lo"][1], fgap["ci_hi"][1], fgap["p"][1]),
             ("Age acceleration\n(residual method)", *adj_or_from_clock(df, cfg))]
    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    yp = np.arange(len(items))[::-1]
    for yi, (lab, o, lo, hi, p) in zip(yp, items):
        col = JAMA["red"] if "Naive" in lab else (JAMA["blue"] if p < 0.05 else JAMA["gray"])
        ax.plot([lo, hi], [yi, yi], color=JAMA["slate"], lw=1.2); ax.plot(o, yi, "s", ms=7, color=col)
        ax.text(hi + 0.1, yi, f"OR {o:.2f}, P={p:.3f}", fontsize=6.5, va="center")
    ax.axvline(1, color=JAMA["gray"], ls="--", lw=0.8)
    ax.set_yticks(yp); ax.set_yticklabels([i[0] for i in items], fontsize=7)
    ax.set_xscale("log"); ax.set_xlim(0.5, 6); ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticks([0.5, 1, 2, 4]); ax.set_xticklabels(["0.5", "1", "2", "4"])
    ax.set_xlabel("Odds ratio (95% CI)")
    ax.set_title("Why naive operationalizations mislead", fontsize=9, weight="bold", color=JAMA["slate"])
    ax.text(0, -0.28, "† Age-gap adjusted for chronological age is a suppression artifact (gap contains age).",
            transform=ax.transAxes, fontsize=5.8, color=JAMA["gray"])
    ax.spines[["top", "right", "left"]].set_visible(False); ax.tick_params(left=False)
    fig.tight_layout(); _save(fig, "S3.Figure_S3_naive_approaches")


def adj_or_from_clock(df, cfg):
    a = pd.read_csv(RES / "aar_association.csv"); r = a[a["model"] == "adj_age_sex_asa"].iloc[0]
    return r["AAR_OR"], r["ci_lo"], r["ci_hi"], r["p"]


# ======================================================================================
# Tables
# ======================================================================================
def _fmt(o, lo, hi, p):
    return f"{o:.2f} ({lo:.2f}–{hi:.2f})", f"{p:.3f}"


def tables(df, cfg):
    RES.mkdir(exist_ok=True); TAB = FIG.parent / "tables"; TAB.mkdir(exist_ok=True)
    outcome = cfg["outcome"]
    # Table 1 — cohort by discharge
    sub = df.copy(); sub["female"] = (sub["sex"] == "F").astype(int)
    grp = sub.groupby(outcome)
    def summ(col, kind="cont"):
        s = pd.to_numeric(sub[col], errors="coerce")
        if kind == "bin":
            return {g: f"{100*pd.to_numeric(d[col],errors='coerce').mean():.0f}%" for g, d in grp}
        return {g: f"{pd.to_numeric(d[col],errors='coerce').median():.1f}" for g, d in grp}
    rows = []
    for lab, col, kind in [("Age, y (median)", "age_yrs", "cont"), ("Female", "female", "bin"),
                           ("ASA class (median)", "asa", "cont"), ("Hypertension", "htn", "bin"),
                           ("Diabetes", "diabetes", "bin"), ("Fusion", "fusion", "bin"),
                           ("No. operated levels (median)", "num_level", "cont")]:
        d = summ(col, kind); rows.append({"Characteristic": lab, "Home": d.get(0, ""), "Non-home": d.get(1, "")})
    pd.DataFrame(rows).to_csv(TAB / "1.Table_1_cohort.csv", index=False)
    # Table 2 — primary association (ladder + specs)
    a = pd.read_csv(RES / "aar_association.csv"); s = pd.read_csv(RES / "clock_specs_sensitivity.csv")
    t2 = []
    for _, r in a.iterrows():
        est, p = _fmt(r["AAR_OR"], r["ci_lo"], r["ci_hi"], r["p"])
        t2.append({"Analysis": {"crude": "Crude", "adj_age": "Adjusted for age",
                                "adj_age_sex_asa": "Adjusted for age, sex, ASA"}[r["model"]],
                   "Clock": "Scanner-robust", "OR per SD (95% CI)": est, "P": p, "n": int(r["n"])})
    for _, r in s.iterrows():
        est, p = _fmt(r["AAR_OR"], r["ci_lo"], r["ci_hi"], r["p"])
        t2.append({"Analysis": "Adjusted for age, sex, ASA",
                   "Clock": {"scanner_robust": "Scanner-robust", "intensity_ratio_only": "Intensity-ratio only",
                             "volume_only": "Volume only"}[r["spec"]],
                   "OR per SD (95% CI)": est, "P": p, "n": int(r["n"])})
    pd.DataFrame(t2).to_csv(TAB / "2.Table_2_primary_association.csv", index=False)
    # Table 3 — clock performance + incremental value
    met = pd.read_csv(RES / "clock_metrics.csv")
    clock, _ = fit_clock(df, cfg); idx = clock.index; subc = df.loc[idx]
    y = clock[outcome].astype(int).to_numpy()
    age_z = _z(subc["age_yrs"].astype(float).to_numpy()); fem = (subc["sex"] == "F").astype(float).to_numpy()
    asa_z = _z(pd.to_numeric(subc["asa"], errors="coerce").fillna(2).to_numpy()); aar = _z(clock["aar"].to_numpy())
    def auc_model(cols):
        f = firth_logit(np.column_stack(cols), y); lp = np.column_stack([np.ones(len(y))] + cols) @ f["beta"]; return _auc(lp, y)
    auc_c = auc_model([age_z, fem, asa_z]); auc_a = auc_model([age_z, fem, asa_z, aar])
    t3 = [{"Metric": "Clock age R² (CV)", "Value": f"{met.iloc[0]['age_R2']:.2f}"},
          {"Metric": "Clock MAE, y", "Value": f"{met.iloc[0]['MAE']:.1f}"},
          {"Metric": "AUC clinical (age, sex, ASA)", "Value": f"{auc_c:.3f}"},
          {"Metric": "AUC clinical + age acceleration", "Value": f"{auc_a:.3f}"},
          {"Metric": "ΔAUC (in-sample)", "Value": f"{auc_a-auc_c:+.3f}"}]
    pd.DataFrame(t3).to_csv(TAB / "3.Table_3_clock_and_value.csv", index=False)
    # Supp Table 1 — features + missingness
    feats = build_features(df); from .aging_clock import FEATURE_LABELS
    n = len(df)
    st1 = [{"Feature": FEATURE_LABELS[c], "Spec": ", ".join(k for k, v in CLOCK_SPECS.items() if c in v),
            "n observed": int(feats[c].notna().sum()), "% missing": round(100*(n-feats[c].notna().sum())/n, 1)}
           for c in feats.columns]
    pd.DataFrame(st1).to_csv(TAB / "S1.Table_S1_features.csv", index=False)
    # Supp Table 2 — seed/lambda sensitivity
    cols = CLOCK_SPECS[PRIMARY_SPEC]
    d = pd.concat([feats[cols], df[["age_yrs", outcome]]], axis=1).dropna()
    X = d[cols].to_numpy(float); agev = d["age_yrs"].to_numpy(float); Xs = (X - X.mean(0)) / X.std(0)
    yv = d[outcome].astype(int).to_numpy()
    st2 = []
    for seed in [1, 7, 42, 2024]:
        for lam in [3.0, 10.0, 30.0]:
            pred = _ridge_cv_age(Xs, agev, lam, seed); b = np.polyfit(agev, pred, 1); aarv = pred - (b[0]*agev+b[1])
            f = firth_logit(np.column_stack([_z(aarv), _z(agev)]), yv)
            st2.append({"seed": seed, "lambda": lam, "AAR_OR": round(f["or"][1], 2), "P": round(f["p"][1], 3)})
    pd.DataFrame(st2).to_csv(TAB / "S2.Table_S2_sensitivity.csv", index=False)


def make_all():
    from .data_loading import load_cohort, load_config
    cfg = load_config("config.yaml"); df = load_cohort(cfg)
    fig1_methods_overview()
    fig2_aging_clock(df, cfg)
    fig3_primary_association(df, cfg)
    fig4_robustness_and_value(df, cfg)
    figS1_strobe(df, cfg)
    figS2_feature_correlation(df, cfg)
    figS3_naive_approaches(df, cfg)
    tables(df, cfg)
    print("figures + tables written")


if __name__ == "__main__":
    make_all()

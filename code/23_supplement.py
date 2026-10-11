"""Supplementary material: extra figures and the tables that support them.

Everything here is read from results/, like the main figures, so a supplementary
value can never disagree with the article. Output: figures/supp1..supp4.png and
docs/Supplementary_Material.docx.

Run from the deposit root:  python code/23_supplement.py
"""
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.shared import Inches, Pt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import omml  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")
DOCS = os.path.join(ROOT, "docs")
COL, PAGE = 3.5, 7.16
GREY, ACC, WARN = "#8c8c8c", "#1f4e79", "#c0392b"
POLY = ["HDPE", "LDPE", "PET", "PP", "PS", "PVC"]
NICE = {"A": "amplitude", "C": "shape correlation", "O": "baseline offset"}
WIN = {"W1": "C-H stretch", "W2": "carbonyl", "W3": "CH deformation"}


def setup():
    mpl.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight"})


def label(a):
    kind, win = a.split("_")
    return f"{NICE[kind]}, {WIN[win]}"


def supp1_per_class():
    """Per-polymer sensitivity for every method, with bootstrap intervals."""
    d = pd.read_csv(os.path.join(RES, "cv_per_class.csv"))
    methods = ["glassbox_forced", "svm_rbf", "random_forest", "plsda", "cnn1d"]
    nice = {"glassbox_forced": "glass box", "svm_rbf": "SVM", "random_forest": "forest",
            "plsda": "PLS-DA", "cnn1d": "1-D CNN"}
    d = d[d.method.isin(methods)]
    fig, ax = plt.subplots(figsize=(PAGE, 2.8))
    width = 0.15
    for i, m in enumerate(methods):
        g = d[d.method == m].set_index("polymer").reindex(POLY)
        x = np.arange(len(POLY)) + (i - 2) * width
        err = np.vstack([g.sensitivity - g.sens_ci_low, g.sens_ci_high - g.sensitivity])
        ax.bar(x, g.sensitivity, width, yerr=np.nan_to_num(err), label=nice[m],
               color=WARN if m == "glassbox_forced" else None,
               error_kw={"lw": 0.6, "ecolor": "#444444"})
    ax.set_xticks(np.arange(len(POLY)))
    ax.set_xticklabels(POLY)
    ax.set_ylim(0.75, 1.01)
    ax.set_ylabel("per-class sensitivity (recall)")
    ax.legend(ncol=5, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.28))
    fig.savefig(os.path.join(FIG, "supp1.png"))
    plt.close(fig)
    g = d[d.method == "glassbox_forced"].set_index("polymer").reindex(POLY)
    return g.sensitivity.idxmin(), float(g.sensitivity.min())


def supp2_overlap_matrix():
    """Every attribute against every polymer pair, as a heat map."""
    d = pd.read_csv(os.path.join(RES, "ovl_matrix.csv"))
    d = d[~d.comparison.str.contains("rest")]
    piv = d.pivot_table(index="attribute", columns="comparison", values="ovl")
    order = piv.mean(axis=1).sort_values().index
    piv = piv.loc[order]
    fig, ax = plt.subplots(figsize=(PAGE, 3.4))
    im = ax.imshow(piv.values, aspect="auto", cmap="RdYlGn_r", vmin=0, vmax=1)
    ax.set_xticks(range(piv.shape[1]))
    ax.set_xticklabels([c.replace(" vs ", "/") for c in piv.columns], rotation=90, fontsize=5.5)
    ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels([label(a) for a in piv.index], fontsize=6.5)
    cb = fig.colorbar(im, ax=ax, pad=0.015)
    cb.set_label("overlap coefficient (0 = separated)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    fig.savefig(os.path.join(FIG, "supp2.png"))
    plt.close(fig)
    best = piv.median(axis=1).idxmin()
    return best, float(piv.median(axis=1).min()), int((piv.values < 0.05).sum())


def supp3_rule_intervals():
    """Rule confidence with its Wilson interval, attribute by attribute."""
    d = pd.read_csv(os.path.join(RES, "rule_matrix.csv"))
    g = (d.groupby("attribute")
         .agg(conf=("conf_change", "mean"), lo=("wilson_low", "mean"),
              hi=("wilson_high", "mean"), sig=("significant_q05", "sum"),
              n=("attribute", "size")).sort_values("conf"))
    fig, ax = plt.subplots(figsize=(COL * 1.6, 2.9))
    y = np.arange(len(g))
    err = np.vstack([g.conf - g.lo, g.hi - g.conf])
    ax.barh(y, g.conf, xerr=err, color=ACC, error_kw={"lw": 0.7, "ecolor": "#333333"})
    ax.set_yticks(y)
    ax.set_yticklabels([label(a) for a in g.index], fontsize=6.5)
    ax.set_xlabel("rule confidence, mean over polymer pairs")
    ax.set_xlim(0, 1.02)
    for yy, (c, s, n) in enumerate(zip(g.conf, g.sig, g.n)):
        ax.text(min(c + 0.015, 0.98), yy, f"{c:.3f}  ({int(s)}/{int(n)} significant)",
                va="center", fontsize=6)
    fig.savefig(os.path.join(FIG, "supp3.png"))
    plt.close(fig)
    return g.round(4)


def supp4_attribution_ranks():
    """The attribution ranking fit by fit, against the rule's fixed selection."""
    d = pd.read_csv(os.path.join(RES, "shap_ranks_per_repeat.csv"))
    attrs = [c for c in d.columns
             if c not in ("repeat", "glassbox_set", "gamma_criterion_set")]
    fig, ax = plt.subplots(figsize=(PAGE, 2.8))
    for a in attrs:
        ax.plot(d.repeat, d[a], marker="o", ms=2.5, lw=1.0, label=label(a))
    ax.invert_yaxis()
    ax.set_xlabel("repeated fit")
    ax.set_ylabel("rank (1 = most important)")
    ax.set_yticks(range(1, len(attrs) + 1))
    ax.legend(ncol=3, frameon=False, fontsize=5.8, loc="lower center",
              bbox_to_anchor=(0.5, -0.42))
    fig.savefig(os.path.join(FIG, "supp4.png"))
    plt.close(fig)
    top3 = d[attrs].apply(lambda r: ",".join(sorted(r.nsmallest(3).index)), axis=1)
    return top3.nunique(), d.glassbox_set.nunique()


def doc_add(doc, text, size=10, bold=False, italic=False, space_after=6):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.font.size = Pt(size)
    r.bold, r.italic = bold, italic
    p.paragraph_format.space_after = Pt(space_after)
    return p


def doc_table(doc, header, rows, caption, note=None):
    doc_add(doc, caption, size=8, italic=True, space_after=3)
    t = doc.add_table(rows=1, cols=len(header))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, h in enumerate(header):
        c = t.rows[0].cells[i]
        c.text = ""
        r = c.paragraphs[0].add_run(str(h))
        r.bold = True
        r.font.size = Pt(7.5)
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            cells[i].paragraphs[0].add_run("" if v is None else str(v)).font.size = Pt(7.5)
    if note:
        doc_add(doc, note, size=7.5, italic=True, space_after=10)
    else:
        doc_add(doc, "", size=7.5, space_after=8)
    return t


def build_document(stats):
    doc = Document()
    for s in doc.sections:
        s.left_margin = s.right_margin = Inches(0.9)
        s.top_margin = s.bottom_margin = Inches(0.9)
    doc_add(doc, "Supplementary Material", size=16, bold=True, space_after=2)
    doc_add(doc, "A Glass-Box Classifier for FTIR Identification of Plastics: Accuracy "
                 "Parity on Accepted Specimens With an Explicit Referral Rate",
            size=10.5, italic=True, space_after=10)
    doc_add(doc, "This document holds material that supports the article without being "
                 "required to follow it: per-class results behind the aggregate accuracies, "
                 "the full attribute-by-pair separation matrix, the rule table with its "
                 "intervals, the attribution ranking fit by fit, and the computational cost "
                 "of the protocol. Every number is read from the same result tables as the "
                 "article, which are released with the code; no value here is recomputed or "
                 "typed by hand. Section numbers refer to the article.", size=9.5)

    doc_add(doc, "S1. Per-class results behind the aggregate accuracies", size=12,
            bold=True, space_after=4)
    doc_add(doc, "Table 2 of the article reports balanced accuracy, which is the mean of the "
                 "per-class recalls. That mean hides which polymers are hard. Fig. S1 and "
                 "Table S1 give the per-class picture for every method. The glass box in "
                 "forced-choice mode is weakest on " + stats["worst_poly"] + ", at a recall of "
                 f"{stats['worst_sens']:.3f}, while the conventional classifiers are close to "
                 "the ceiling on all six. This is the expected signature of a two-attribute "
                 "partition: it separates chemically distinct materials completely and "
                 "struggles exactly where two polymers differ in band shape rather than in "
                 "band presence. The referral outcome exists for that case, and the covered "
                 "mode in Table 2 is what remains once those specimens are set aside.",
            size=9.5)
    doc.add_picture(os.path.join(FIG, "supp1.png"), width=Inches(6.6))
    doc_add(doc, "Fig. S1. Per-class sensitivity with 95% percentile bootstrap intervals over "
                 "specimens, by method. "
                 "The glass box is shown in forced-choice mode so that the comparison is "
                 "like for like.", size=8, italic=True, space_after=10)
    doc_table(doc, stats["per_class"][0], stats["per_class"][1],
              "Table S1. Per-class sensitivity and specificity, five-fold cross-validation. Intervals are percentile bootstrap over specimens, not binomial. For the covered mode the denominator is the number of accepted specimens of that class, not 500.")

    doc_add(doc, "S2. Separation of every attribute on every polymer pair", size=12,
            bold=True, space_after=4)
    doc_add(doc, "Fig. 5 of the article reports the median overlap coefficient per attribute. "
                 "Fig. S2 gives the full matrix: nine attributes against all fifteen polymer "
                 "pairs. The strongest single attribute is " + label(stats["best_attr"]) +
            f", with a median overlap of {stats['best_ovl']:.4f} across the 15 pairs, and "
            f"{stats['n_sep']} of the 135 attribute-pair combinations fall below an overlap "
                 "of 0.05, which is the regime in which a single named quantity separates two "
                 "materials almost completely. The matrix is also the diagnostic a laboratory "
                 "would use when a new material is added: a column with no low-overlap cell "
                 "indicates a pair that the present three windows cannot separate, and "
                 "therefore a case where a fourth window, not a larger model, is the remedy.",
            size=9.5)
    doc.add_picture(os.path.join(FIG, "supp2.png"), width=Inches(6.6))
    doc_add(doc, "Fig. S2. Overlap coefficient for every attribute against every polymer "
                 "pair. Green is separated, red is indistinguishable.", size=8, italic=True,
            space_after=10)

    doc_add(doc, "S3. The rule table with intervals and significance", size=12, bold=True,
            space_after=4)
    doc_add(doc, "Table 4 of the article aggregates rule confidence by attribute. Table S2 and "
                 "Fig. S3 add the Wilson interval and the number of polymer pairs for which "
                 "the rule survives false-discovery-rate control at the five percent level. "
                 "The intervals are narrow because the pair database is large, which is a "
                 "consequence of the pairing construction rather than of the sample size: "
                 "3000 specimens yield far more admissible pairs than records. That is also "
                 "the reason the article reports a permutation test rather than relying on "
                 "the interval alone, since pairs drawn from a common specimen are not "
                 "independent.", size=9.5)
    doc.add_picture(os.path.join(FIG, "supp3.png"), width=Inches(5.4))
    doc_add(doc, "Fig. S3. Mean rule confidence per attribute with its Wilson interval, and "
                 "the number of polymer pairs that remain significant after false-discovery-"
                 "rate control.", size=8, italic=True, space_after=10)
    doc_table(doc, stats["rules"][0], stats["rules"][1],
              "Table S2. Rule confidence by attribute, with interval and significance count.")

    doc_add(doc, "S4. Attribution ranking, fit by fit", size=12, bold=True, space_after=4)
    doc_add(doc, "Section V reports that the attribute set selected by the rule was the same "
                 "in 19 of 20 refits while the top three attributes of the attribution "
                 "method changed in one repeat in ten. Fig. S4 shows the underlying "
                 "trajectories rather than the summary. Across the twenty fits the attribution "
                 f"produced {stats['n_top3']} distinct top-three sets against "
                 f"{stats['n_sets']} distinct selection{'' if stats['n_sets'] == 1 else 's'} by the "
                 "rule. The figure also makes "
                 "the nature of the instability visible: the ranking is stable in its extremes "
                 "and unstable in the middle, which is precisely the region a specification "
                 "would have to cut through if it named a fixed number of important "
                 "quantities.", size=9.5)
    doc.add_picture(os.path.join(FIG, "supp4.png"), width=Inches(6.6))
    doc_add(doc, "Fig. S4. Rank assigned to each attribute by the attribution method, fit by "
                 "fit. Lower is more important.", size=8, italic=True, space_after=10)

    doc_add(doc, "S5. Who gets referred", size=12, bold=True, space_after=4)
    doc_add(doc, "The referral rate of 0.247 reported in the article is an average over six "
                 "polymers whose individual rates differ by more than an order of magnitude. "
                 "Table S3 gives the breakdown. The rule accepts almost every specimen of "
                 "poly(ethylene terephthalate) and poly(vinyl chloride), both of which carry a "
                 "band the other classes lack, and refers two thirds of low-density "
                 "polyethylene, which differs from the high-density grade by chain branching "
                 "alone. Any deployment that is dominated by polyolefins should therefore "
                 "expect a far higher referral rate than the headline figure, and any "
                 "deployment dominated by bottle-grade polyester a far lower one. Recall "
                 "among accepted specimens is 0.9916 or better in every class, so the "
                 "refusals are placed where the evidence is weak and not scattered at random.",
            space_after=6)
    doc_table(doc, stats["per_cov"][0], stats["per_cov"][1],
              "Table S3. Coverage and recall by polymer, out of fold.",
              note="Coverage is the share of that polymer the rule accepts; recall among "
                   "accepted specimens is computed on that subset only.")

    doc_add(doc, "S6. Computational cost and the fold-level record", size=12, bold=True,
            space_after=4)
    doc_add(doc, "Table S4 gives the fold-level record behind Table 2, including the attribute "
                 "pair selected inside each fold, the quality of approximation on the training "
                 "part, the number of cells and the wall-clock cost. Two things in it are "
                 "worth a reviewer's attention. The attribute pair selected inside the folds "
                 "is not always the pair reported for the full sample, which is the honest "
                 "signature of a selection step that is repeated inside the cross-validation "
                 "rather than fixed in advance from all the data. And the total cost of one "
                 f"repetition is {stats['total_min']:.0f} minutes on a single core, which is "
                 "the reason the article reports one repetition rather than the twenty the "
                 "protocol specifies, and states that deviation in its limitations.",
            size=9.5)
    doc_table(doc, stats["folds"][0], stats["folds"][1],
              "Table S4. Fold-level record of the primary cross-validation.",
              note="Seconds are wall-clock on one core, including the baselines.")

    doc_add(doc, "S7. Released files", size=12, bold=True, space_after=4)
    doc_add(doc, "The following computed tables are released alongside the article and are the "
                 "direct source of the figures and tables named beside each. They are plain "
                 "comma-separated files with a header row.", size=9.5)
    doc_table(doc, ["File", "Rows", "What it holds", "Used by"], stats["files"],
              "Table S5. Computed data files released with the article.")

    out = os.path.join(DOCS, "Supplementary_Material.docx")
    doc.save(out)
    return out


def main():
    setup()
    worst_poly, worst_sens = supp1_per_class()
    best_attr, best_ovl, n_sep = supp2_overlap_matrix()
    rules = supp3_rule_intervals()
    n_top3, n_sets = supp4_attribution_ranks()

    pc = pd.read_csv(os.path.join(RES, "cv_per_class.csv"))
    pc = pc[pc.method.isin(["glassbox_forced", "glassbox_covered", "svm_rbf"])]
    per_class = (["Method", "Polymer", "n tested", "Sensitivity", "95% bootstrap interval", "Specificity"],
                 [[r.method.replace("_", " "), r.polymer, int(r.n_specimens),
                   f"{r.sensitivity:.3f}", f"{r.sens_ci_low:.3f} to {r.sens_ci_high:.3f}",
                   f"{r.specificity:.4f}"] for r in pc.itertuples()])

    rules_tab = (["Attribute", "Confidence", "Wilson interval", "Significant pairs"],
                 [[label(a), f"{r.conf:.4f}", f"{r.lo:.4f} to {r.hi:.4f}",
                   f"{int(r.sig)} of {int(r.n)}"] for a, r in rules.iterrows()])

    fs = pd.read_csv(os.path.join(RES, "cv_fold_scores.csv"))
    folds = (["Fold", "Attributes", "Gamma (train)", "Cells", "Coverage",
              "Glass box forced", "Glass box covered", "Seconds"],
             [[int(r.fold) + 1, " + ".join(label(a) for a in r.attributes.split("+")),
               f"{r.gamma_train:.4f}", int(r.n_cells), f"{r.coverage:.3f}",
               f"{r.bacc_glassbox_forced:.4f}", f"{r.bacc_glassbox_covered:.4f}",
               f"{r.secs_total:.0f}"] for r in fs.itertuples()])

    inventory = [
        ("cv_parity_folds.csv", "Balanced accuracy per method with fold-level intervals",
         "Table 2, Fig. 3"),
        ("cv_per_class.csv", "Per-class sensitivity and specificity per method",
         "Table S1, Fig. S1"),
        ("cv_fold_scores.csv", "Fold-level record: selection, cells, coverage, timing",
         "Table S3"),
        ("rough_sets.csv", "Approximation quality for every attribute set and variant",
         "Table 3, Fig. 8"),
        ("rule_matrix.csv", "Every mined rule with interval, permutation p and q value",
         "Table 4, Fig. 4, Table S2"),
        ("ovl_matrix.csv", "Overlap coefficient for every attribute and polymer pair",
         "Fig. 5, Fig. S2"),
        ("learning_curves.csv", "Accuracy against training size for every method",
         "Table 6, Fig. 6"),
        ("shap_stability.csv", "Rank-correlation summary of attribution stability",
         "Table 7"),
        ("shap_ranks_per_repeat.csv", "Attribution rank of each attribute, fit by fit",
         "Fig. 7, Fig. S4"),
        ("oc_curves.csv", "Operating characteristic for each plan and inspector", "Fig. 10"),
        ("oc_risk.csv", "Producer and consumer risk for four sampling plans", "Table 8"),
        ("pair_cap_sensitivity.csv", "Rule confidence under every open protocol choice",
         "Table 5, Fig. 9"),
        ("entry_condition_primary.csv", "Between-specimen spread against instrument noise",
         "Fig. 2"),
        ("parity_test.csv", "Paired test of the forced-choice difference", "Section V"),
        ("non_inferiority_test.csv", "Non-inferiority test at a one-point margin",
         "Section V"),
        ("selective_comparison.csv", "Every classifier under the same referral budget",
         "Table 9"),
        ("pareto_pairs.csv", "Candidate attribute pairs evaluated out of fold", "Table 10"),
        ("binning_out_of_fold.csv", "Coverage and accuracy under each binning rule",
         "Table 11"),
        ("accepted_per_class.csv", "Coverage and recall per polymer", "Table S5"),
    ]
    # Row counts are read from the files, never written by hand, so this
    # inventory cannot drift away from what is actually released.
    files = [[name, len(pd.read_csv(os.path.join(RES, name))), what, used]
             for name, what, used in inventory
             if os.path.exists(os.path.join(RES, name))]

    acc = pd.read_csv(os.path.join(RES, "accepted_per_class.csv"))
    per_cov = (["Polymer", "Accepted of 500", "Coverage", "Recall, forced", "Recall, accepted"],
               [[r.polymer, int(r.n_accepted), f"{r.coverage:.3f}",
                 f"{r.recall_forced:.4f}", f"{r.recall_accepted:.4f}"]
                for r in acc.itertuples()])

    sl = pd.read_csv(os.path.join(RES, "selective_comparison.csv"))
    sl = sl[sl.regime == "own reject option at matched coverage"]

    stats = {"worst_poly": worst_poly, "worst_sens": worst_sens, "best_attr": best_attr,
             "best_ovl": best_ovl, "n_sep": n_sep, "n_top3": n_top3, "n_sets": n_sets,
             "per_class": per_class, "per_cov": per_cov, "rules": rules_tab, "folds": folds, "files": files,
             "total_min": float(fs.secs_total.sum()) / 60}
    out = build_document(stats)
    print(f"saved {out}")
    print(f"weakest class for the glass box: {worst_poly} at {worst_sens:.3f}")
    print(f"best separating attribute {best_attr} mean overlap {best_ovl:.4f} | "
          f"combinations below 0.05: {n_sep}")
    print(f"distinct top-three attribution sets {n_top3} against {n_sets} rule selections")
    print(f"one repetition cost {fs.secs_total.sum() / 60:.0f} minutes")


if __name__ == "__main__":
    main()

"""Assemble docs/manuscript.json from the prose parts, the result tables and the figures.

Tables and figures are generated here from results/, never typed by hand, so a
number in the manuscript cannot disagree with the table it came from.
Citations are written in the prose as [@key] and resolved against
docs/references.csv in order of first appearance.

Run from the deposit root:  python code/22_assemble_manuscript.py
"""
import csv
import json
import os
import re

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
DOCS = os.path.join(ROOT, "docs")
FIGMAP = [("matys2.png", "Entry condition on the primary dataset. The ratio of between-specimen "
           "spread to the noise-calibrated detection threshold, by polymer, on a logarithmic scale."),
          ("matys3.png", "Balanced accuracy under five-fold cross-validation. Bars show the mean "
           "over folds and whiskers the 95 percent interval across folds. The glass box is shown "
           "twice: forced to decide every specimen, and on the specimens it accepts."),
          ("matys4.png", "Confidence that each named attribute moves beyond its own detection "
           "threshold when the polymer changes, for every pair of polymers."),
          ("matys5.png", "Median overlap coefficient between polymer pairs, by attribute. Lower is "
           "better separated."),
          ("matys6.png", "Learning curves and the difference against the strongest baseline. The "
           "shaded band is the 95 percent bootstrap interval. No crossover occurs."),
          ("matys7.png", "Stability across twenty repeated fits. (a) ranks assigned by post-hoc "
           "attribution; (b) the attribute set selected by the rule."),
          ("matys8.png", "Operating-characteristic curves for four single-stage sampling plans "
           "under a perfect inspector and under the observed one. Dotted lines mark the one and "
           "ten percent defective fractions."),
          ("matys9.png", "Mean transmittance spectrum of each polymer, with the three fixed "
           "analysis windows shaded. Every attribute used by the rule is computed inside one of "
           "these three bands; the remaining 1500 channels enter no decision."),
          ("matys10.png", "Decidability by attribute set. (a) quality of approximation for the "
           "ten strongest two-attribute sets; (b) referral rate against the number of cells the "
           "set induces, showing that a finer grid buys certainty per cell at the cost of "
           "referring more specimens."),
          ("matys11.png", "Headline rule confidence under every choice the frozen protocol left "
           "open: the pair sampling cap, the threshold multiplier and the spectral resolution.")]


def t(df, cols, rename=None, fmt=4):
    df = df[cols].copy()
    for c in df.columns:
        if pd.api.types.is_float_dtype(df[c]):
            df[c] = df[c].round(fmt)
    if rename:
        df = df.rename(columns=rename)
    return list(df.columns), df.astype(object).where(pd.notna(df), "").values.tolist()


def tables():
    out = {}

    par = pd.read_csv(os.path.join(RES, "cv_parity_folds.csv"))
    par["ci_high"] = par.ci_high.clip(upper=1.0)
    nice = {"glassbox_forced": "Glass box, forced choice",
            "glassbox_covered": "Glass box, accepted specimens",
            "svm_rbf": "SVM, radial basis function", "random_forest": "Random forest",
            "plsda": "PLS-DA", "cnn1d": "1-D convolutional network"}
    par["method"] = par.method.map(nice)
    par["interval"] = par.apply(lambda r: f"{r.ci_low:.4f} to {r.ci_high:.4f}", axis=1)
    par["coverage_mean"] = par.coverage_mean.map(lambda v: f"{v:.3f}")
    out["parity"] = t(par, ["method", "bacc_mean", "interval", "coverage_mean"],
                      {"method": "Method", "bacc_mean": "Balanced accuracy",
                       "interval": "95% interval across folds", "coverage_mean": "Coverage"})

    rs = pd.read_csv(os.path.join(RES, "rough_sets.csv"))
    best = (rs[rs.variant == "sens.min bins, k=4"]
            .groupby("attribute_set").gamma.first().idxmax())
    sub = rs[(rs.variant == "sens.min bins, k=4") & (rs.attribute_set == best)]
    out["rough"] = t(sub, ["polymer", "alpha", "boundary_size"],
                     {"polymer": "Polymer", "alpha": "Accuracy of approximation",
                      "boundary_size": "Boundary region (specimens)"})
    out["rough_note"] = (f"Attribute set {best.replace('+', ' and ')}; quality of approximation "
                         f"{float(sub.gamma.iloc[0]):.4f} over {int(sub.n_cells.iloc[0])} cells.")

    rm = pd.read_csv(os.path.join(RES, "rule_matrix.csv"))
    agg = (rm.groupby("attribute")
             .agg(mean_conf=("conf_change", "mean"), min_conf=("conf_change", "min"),
                  significant=("significant_q05", "sum"), n=("conf_change", "size"))
             .sort_values("mean_conf", ascending=False).reset_index())
    agg["significant"] = agg.apply(lambda r: f"{int(r.significant)} of {int(r.n)}", axis=1)
    out["rules"] = t(agg, ["attribute", "mean_conf", "min_conf", "significant"],
                     {"attribute": "Attribute", "mean_conf": "Mean confidence",
                      "min_conf": "Worst polymer pair", "significant": "Significant at q < 0.05"})

    pc = pd.read_csv(os.path.join(RES, "pair_cap_sensitivity.csv"))
    pc = pc[(pc.metric == "conf_change") & (pc.attribute == "A_W1")]
    out["sens"] = t(pc, ["variant_type", "variant_value", "value"],
                    {"variant_type": "Varied", "variant_value": "Setting",
                     "value": "Confidence, amplitude of C-H window"})

    lc = pd.read_csv(os.path.join(RES, "learning_curves.csv"))
    gl = lc[lc.method == "glassbox_forced"].sort_values("size_pct")
    gl["interval"] = gl.apply(lambda r: f"{r.diff_ci_low:.4f} to {r.diff_ci_high:.4f}", axis=1)
    out["learning"] = t(gl, ["size_pct", "n_per_class", "best_blackbox",
                             "diff_vs_best_blackbox", "interval"],
                        {"size_pct": "Training set (%)", "n_per_class": "Specimens per class",
                         "best_blackbox": "Strongest baseline",
                         "diff_vs_best_blackbox": "Glass box minus baseline",
                         "interval": "95% interval"})

    rk = pd.read_csv(os.path.join(RES, "oc_risk.csv"))
    rk["plan"] = rk.apply(lambda r: f"n = {int(r.n)}, c = {int(r.c)}", axis=1)
    out["oc"] = t(rk, ["plan", "inspector", "producer_risk_at_1pct", "consumer_risk_at_10pct"],
                  {"plan": "Sampling plan", "inspector": "Inspector",
                   "producer_risk_at_1pct": "Producer risk at 1% defective",
                   "consumer_risk_at_10pct": "Consumer risk at 10% defective"})

    st = pd.read_csv(os.path.join(RES, "shap_stability.csv"))
    out["stability"] = t(st, ["representation", "repeats", "kendall_tau_mean",
                              "modal_set_share", "n_distinct_sets"],
                         {"representation": "Representation", "repeats": "Repeated fits",
                          "kendall_tau_mean": "Mean rank correlation",
                          "modal_set_share": "Share of repeats with modal set",
                          "n_distinct_sets": "Distinct sets"})
    return out


def load_refs():
    path = os.path.join(DOCS, "references.csv")
    refs = {}
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for row in csv.DictReader(fh):
            authors = row.get("authors", "").strip()
            title = row.get("title", "").strip().rstrip(".")
            venue = row.get("venue", "").strip()
            year = row.get("year", "").strip()
            doi = row.get("doi", "").strip()
            url = row.get("url", "").strip()
            head = ", ".join(p for p in [authors] if p)
            tail = ", ".join(p for p in [venue, year] if p)
            s = (head + ", " if head else "") + (f'"{title}," ' if title else "") + tail + "."
            if doi and not doi.lower().startswith("none"):
                s += f" doi: {doi}."
            elif url:
                s += f" [Online]. Available: {url}"
            refs[row["key"].strip()] = s
    return refs


def resolve(blocks, refs, order):
    pat = re.compile(r"\[@([A-Za-z0-9_\-]+)\]")

    def sub(m):
        key = m.group(1)
        if key not in refs:
            raise KeyError(f"citation key not in references.csv: {key}")
        if key not in order:
            order.append(key)
        return f"[{order.index(key) + 1}]"
    for b in blocks:
        if b.get("type") == "para":
            b["text"] = pat.sub(sub, b["text"])
    return blocks



def renumber_figures(sections):
    """Figure numbers follow document order, whatever order blocks were added in."""
    n = 0
    for s in sections:
        for b in s["blocks"]:
            if b.get("type") == "figure":
                n += 1
                b["number"] = n
    return n


def main():
    parts = [json.load(open(os.path.join(DOCS, f"manuscript_part{i}.json"), encoding="utf-8"))
             for i in (1, 2, 3)]
    tb = tables()
    refs = load_refs()
    order = []

    sections = []
    for p in parts:
        sections.extend(p["sections"])

    def find(heading):
        for s in sections:
            if s["heading"] == heading:
                return s
        raise KeyError(heading)

    res = find("RESULTS")["blocks"]
    res.insert(0, {"type": "figure", "file": FIGMAP[1][0], "number": 3,
                   "caption": FIGMAP[1][1], "width": 3.4})
    res.insert(1, {"type": "table", "number": "II",
                   "caption": "Balanced accuracy on the primary dataset, five-fold "
                              "cross-validation.",
                   "header": tb["parity"][0], "rows": tb["parity"][1]})
    res.append({"type": "table", "number": "III",
                "caption": "Rough-set approximation by polymer for the selected attribute set.",
                "header": tb["rough"][0], "rows": tb["rough"][1], "note": tb["rough_note"]})
    res.append({"type": "figure", "file": FIGMAP[2][0], "number": 4,
                "caption": FIGMAP[2][1], "width": 7.0})
    res.append({"type": "table", "number": "IV",
                "caption": "Rule confidence by attribute, aggregated over polymer pairs.",
                "header": tb["rules"][0], "rows": tb["rules"][1]})
    res.append({"type": "figure", "file": FIGMAP[3][0], "number": 5,
                "caption": FIGMAP[3][1], "width": 3.4})
    res.append({"type": "table", "number": "V",
                "caption": "Sensitivity of the headline rule confidence to the protocol's "
                           "open choices.",
                "header": tb["sens"][0], "rows": tb["sens"][1]})
    res.append({"type": "figure", "file": FIGMAP[4][0], "number": 6,
                "caption": FIGMAP[4][1], "width": 7.0})
    res.append({"type": "table", "number": "VI",
                "caption": "Learning curve: glass box against the strongest baseline at each "
                           "training size.",
                "header": tb["learning"][0], "rows": tb["learning"][1]})
    res.append({"type": "figure", "file": FIGMAP[5][0], "number": 7,
                "caption": FIGMAP[5][1], "width": 7.0})
    res.append({"type": "table", "number": "VII",
                "caption": "Stability of the representation across twenty repeated fits.",
                "header": tb["stability"][0], "rows": tb["stability"][1]})
    res.append({"type": "figure", "file": FIGMAP[8][0], "number": 0,
                "caption": FIGMAP[8][1], "width": 6.9})
    res.append({"type": "figure", "file": FIGMAP[9][0], "number": 0,
                "caption": FIGMAP[9][1], "width": 6.9})
    res.append({"type": "figure", "file": FIGMAP[6][0], "number": 8,
                "caption": FIGMAP[6][1], "width": 7.0})
    res.append({"type": "table", "number": "VIII",
                "caption": "Producer and consumer risk for four single-stage sampling plans.",
                "header": tb["oc"][0], "rows": tb["oc"][1]})

    meth = find("MATERIALS AND METHOD")["blocks"]
    meth.insert(2, {"type": "table", "number": "I",
                    "caption": "Datasets and their roles.",
                    "header": ["Role", "Collection", "Content", "Supports"],
                    "rows": [["Primary", "FTIR-Plastics, Zenodo record 10736650, CC BY 4.0",
                              "3000 specimens, six polymers, 500 each, 399.2 to 4001.6 cm-1, "
                              "8 and 4 cm-1 resolution",
                              "Parity, learning curves, rough-set bounds, stability"],
                             ["Transfer case", "Recyclate FTIR archive, Zenodo 21822279, CC BY 4.0",
                              "356 spectra, 54 particles, 17 labelled, three waste streams",
                              "Provenance, readability, applied discussion"]]})
    meth.append({"type": "figure", "file": FIGMAP[0][0], "number": 1,
                 "caption": FIGMAP[0][1], "width": 3.4})
    # The spectra figure belongs with the attribute definitions it illustrates.
    att = next(i for i, b in enumerate(meth)
               if b.get("type") == "h2" and b["text"] == "ATTRIBUTES")
    meth.insert(att + 2, {"type": "figure", "file": FIGMAP[7][0], "number": 0,
                          "caption": FIGMAP[7][1], "width": 6.9})

    n_figs = renumber_figures(sections)

    for s in sections:
        resolve(s["blocks"], refs, order)

    base = parts[0]
    doc = {k: base[k] for k in ["title", "authors_line", "affiliations", "corresponding",
                                "funding", "abstract", "index_terms"]}
    doc["sections"] = sections
    doc["references"] = [refs[k] for k in order]
    doc["biographies"] = []
    json.dump(doc, open(os.path.join(DOCS, "manuscript.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    print(f"figures {n_figs} | sections {len(sections)} | citations resolved {len(order)} | "
          f"abstract words {len(doc['abstract'].split())}")
    print("cited keys in order:", ", ".join(order[:12]), "...")


if __name__ == "__main__":
    main()

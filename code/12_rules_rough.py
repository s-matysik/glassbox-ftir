"""Task 3: rule matrix, rough-set bounds, overlap coefficients, sensitivity.

These are structural descriptions of the primary dataset rather than
predictions, so they are computed on all 3000 specimens with one template and
one sens.min calibration; no held-out set is involved and none is implied.

The only object parameter that varies in FTIR-Plastics is polymer identity
(Technic, BR and RST are constant), so the body of every rule is a change of
polymer and the head is a move of an ACO attribute beyond sens.min. Each
specimen is its own group and carries one spectrum, so the replicate weight
1/(r_g * r_h) of concept Eq. (A12) equals one for every pair; the weighted and
unweighted support coincide here and that is stated rather than hidden.

Run from the deposit root:  python code/12_rules_rough.py
Writes results/rule_matrix.csv, results/rough_sets.csv,
results/ovl_matrix.csv, results/pair_cap_sensitivity.csv.
"""
from __future__ import annotations

import itertools
import os
import sys
import time
from dataclasses import replace

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ara_aco as aa  # noqa: E402
import plastics_data as pdata  # noqa: E402

CAP_DEFAULT = 2000
CAPS = [500, 2000, 5000]
KS = [3.0, 4.0, 5.0]
N_PERM = int(os.environ.get("N_PERM", 200))
N_MC = 100
EQUAL_BINS = 10
LOG = os.path.join(pdata.ROOT, "logs", "12_rules_rough.log")


def log(msg):
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(f"{time.strftime('%H:%M:%S')} {msg}\n")


def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    q = np.empty(m)
    prev = 1.0
    for rank, j in enumerate(order[::-1]):
        prev = min(prev, p[j] * m / (m - rank))
        q[j] = prev
    return q


def prepare(resolution):
    grid, Yraw, meta = pdata.load(resolution)
    y = meta["Polymer"].to_numpy().astype(str)
    Yp = aa.preprocess(Yraw)
    sigma = aa.instrument_noise(Yp)
    masks = aa.window_masks(grid, pdata.WINDOWS)
    template = aa.make_template(Yp, masks)
    F, names = aa.aco_matrix(Yp, masks, template)
    D = aa.delta(F, names)
    sens4 = aa.calibrate_sens_min(Yp, sigma, y, masks, template, names,
                                  k=4.0, n_mc=N_MC, seed=pdata.SEED)
    return dict(y=y, F=F, D=D, names=names, sens4=sens4,
                groups=np.arange(len(y)), classes=sorted(set(y)))


def sens_at(sens4, k):
    return replace(sens4, thresholds={n: k * v for n, v in sens4.sd_max.items()},
                   k=k, confidence=aa.implied_confidence(k, sens4.n_effective))


def rule_block(ctx, sens, cap, seed, n_perm, subset=None, body="any polymer change"):
    sel = np.arange(len(ctx["y"])) if subset is None else subset
    out = aa.pair_rules(ctx["F"][sel], ctx["names"], ctx["y"][sel], ctx["groups"][sel],
                        sens, max_pairs_per_stratum=cap, seed=seed, n_perm=n_perm)
    rows = []
    for r in out["rules"]:
        rows.append({"body": body, "attribute": r["attribute"],
                     "support_weighted": r["support"], "conf_change": r["conf_change"],
                     "conf_invariant": r["conf_invariant"],
                     "direction": "up" if r["direction"] > 0 else "down",
                     "conf_direction": r["conf_direction"],
                     "wilson_low": r["ci_low"], "wilson_high": r["ci_high"],
                     "n_pairs": r["n_pairs"], "p_perm": r.get("p_perm", np.nan)})
    return rows


def rough_block(ctx, sens, n_bins, label):
    rows = []
    names = ctx["names"]
    sets = [[a] for a in names] + [list(c) for c in itertools.combinations(names, 2)]
    for att in sets:
        m = aa.fit_rough(ctx["D"], names, ctx["y"], att, sens, n_bins=n_bins)
        for cls in m.classes:
            rows.append({"variant": label, "attribute_set": "+".join(att),
                         "n_attributes": len(att), "polymer": cls,
                         "alpha": m.alpha[cls], "boundary_size": m.boundary[cls],
                         "gamma": m.gamma, "n_cells": len(m.cell_class)})
    return rows


def main():
    t0 = time.time()
    ctx = prepare("c8")
    names, y, classes = ctx["names"], ctx["y"], ctx["classes"]
    log(f"prepared c8 in {time.time()-t0:.0f}s; sens.min k=4: "
        f"{ {k: round(v,5) for k,v in ctx['sens4'].thresholds.items()} }")

    # ---- rule matrix -----------------------------------------------------
    rows = rule_block(ctx, ctx["sens4"], CAP_DEFAULT, pdata.SEED, N_PERM)
    for a, b in itertools.combinations(classes, 2):
        sel = np.where((y == a) | (y == b))[0]
        rows += rule_block(ctx, ctx["sens4"], CAP_DEFAULT, pdata.SEED, N_PERM,
                           subset=sel, body=f"{a} -> {b}")
    rules = pd.DataFrame(rows)
    rules["q_fdr"] = bh(rules["p_perm"].fillna(1.0).to_numpy())
    rules["significant_q05"] = rules["q_fdr"] < 0.05
    pdata.write_csv(rules, "rule_matrix.csv")
    log(f"rule matrix done at {time.time()-t0:.0f}s")

    # ---- rough sets ------------------------------------------------------
    rough = pd.DataFrame(rough_block(ctx, ctx["sens4"], None, "sens.min bins, k=4"))
    pdata.write_csv(rough, "rough_sets.csv")
    log(f"rough sets done at {time.time()-t0:.0f}s")

    # ---- overlap coefficients -------------------------------------------
    ovl = []
    for j, nm in enumerate(names):
        col = ctx["D"][:, j]
        for a, b in itertools.combinations(classes, 2):
            ovl.append({"attribute": nm, "comparison": f"{a} vs {b}",
                        "ovl": aa.overlap_coefficient(col[y == a], col[y == b])})
        for c in classes:
            ovl.append({"attribute": nm, "comparison": f"{c} vs rest",
                        "ovl": aa.overlap_coefficient(col[y == c], col[y != c])})
    ovl = pd.DataFrame(ovl)
    pdata.write_csv(ovl, "ovl_matrix.csv")

    # ---- sensitivity -----------------------------------------------------
    sens_rows = []

    def emit(vtype, vvalue, metric, attribute, value):
        sens_rows.append({"variant_type": vtype, "variant_value": vvalue,
                          "metric": metric, "attribute": attribute, "value": value})

    def rule_summary(vtype, vvalue, sens, cap, ctx_local=None):
        c = ctx_local or ctx
        rr = pd.DataFrame(rule_block(c, sens, cap, pdata.SEED, N_PERM))
        rr["q_fdr"] = bh(rr["p_perm"].fillna(1.0).to_numpy())
        for _, r in rr.iterrows():
            emit(vtype, vvalue, "conf_change", r["attribute"], r["conf_change"])
            emit(vtype, vvalue, "wilson_low", r["attribute"], r["wilson_low"])
            emit(vtype, vvalue, "wilson_high", r["attribute"], r["wilson_high"])
        top = rr.loc[rr["conf_change"].idxmax()]
        emit(vtype, vvalue, "top_rule_conf", top["attribute"], top["conf_change"])
        emit(vtype, vvalue, "top_rule_wilson_low", top["attribute"], top["wilson_low"])
        emit(vtype, vvalue, "n_significant_q05", "", float((rr["q_fdr"] < 0.05).sum()))
        emit(vtype, vvalue, "mean_conf_change", "", float(rr["conf_change"].mean()))
        emit(vtype, vvalue, "n_pairs", "", float(rr["n_pairs"].iloc[0]))
        return rr

    def rough_summary(vtype, vvalue, sens, n_bins, ctx_local=None):
        c = ctx_local or ctx
        rb = pd.DataFrame(rough_block(c, sens, n_bins, f"{vtype}={vvalue}"))
        g = rb.groupby("attribute_set").agg(gamma=("gamma", "first"),
                                            n_attr=("n_attributes", "first"),
                                            alpha=("alpha", "mean"),
                                            cells=("n_cells", "first")).reset_index()
        for na, tag in ((1, "best_single"), (2, "best_pair")):
            sub = g[g.n_attr == na].sort_values("gamma", ascending=False).iloc[0]
            emit(vtype, vvalue, f"gamma_{tag}", sub["attribute_set"], sub["gamma"])
            emit(vtype, vvalue, f"alpha_mean_{tag}", sub["attribute_set"], sub["alpha"])
            emit(vtype, vvalue, f"n_cells_{tag}", sub["attribute_set"], sub["cells"])
        return g

    for cap in CAPS:
        rule_summary("pair_cap", str(cap), ctx["sens4"], cap)
    for k in KS:
        s = sens_at(ctx["sens4"], k)
        emit("k", str(k), "implied_confidence", "", s.confidence)
        rule_summary("k", str(k), s, CAP_DEFAULT)
        rough_summary("k", str(k), s, None)
    rough_summary("bins", "sens.min width", ctx["sens4"], None)
    rough_summary("bins", f"{EQUAL_BINS} equal width", ctx["sens4"], EQUAL_BINS)
    log(f"c8 sensitivity done at {time.time()-t0:.0f}s")

    ctx4 = prepare("c4")
    rule_summary("resolution", "c4 (4 cm-1)", ctx4["sens4"], CAP_DEFAULT, ctx_local=ctx4)
    rough_summary("resolution", "c4 (4 cm-1)", ctx4["sens4"], None, ctx_local=ctx4)
    rule_summary("resolution", "c8 (8 cm-1)", ctx["sens4"], CAP_DEFAULT)
    rough_summary("resolution", "c8 (8 cm-1)", ctx["sens4"], None)
    sens_df = pd.DataFrame(sens_rows)
    pdata.write_csv(sens_df, "pair_cap_sensitivity.csv")

    # ---- printed numbers -------------------------------------------------
    pooled = rules[rules.body == "any polymer change"].sort_values("conf_change", ascending=False)
    print(f"pair database: {int(pooled['n_pairs'].iloc[0])} pairs at cap {CAP_DEFAULT}, "
          f"replicate weight 1 for every pair")
    print("rule matrix, pooled polymer change (top rows):")
    print(pooled[["attribute", "conf_change", "wilson_low", "wilson_high",
                  "conf_invariant", "direction", "q_fdr"]].round(4).to_string(index=False))
    print(f"rules significant at q<0.05 over all {len(rules)} rules: "
          f"{int(rules['significant_q05'].sum())}")
    g = rough.groupby("attribute_set").agg(gamma=("gamma", "first"),
                                           n_attr=("n_attributes", "first"),
                                           alpha=("alpha", "mean")).reset_index()
    for na in (1, 2):
        b = g[g.n_attr == na].sort_values("gamma", ascending=False).iloc[0]
        print(f"best {na}-attribute rough-set model: {b['attribute_set']} "
              f"gamma {b['gamma']:.4f} mean alpha {b['alpha']:.4f}")
    print("lowest overlap coefficients against polymer identity:")
    print(ovl.sort_values("ovl").head(8).round(4).to_string(index=False))
    piv = sens_df[sens_df.metric == "top_rule_conf"]
    print("top-rule confidence under every sensitivity variant:")
    print(piv[["variant_type", "variant_value", "attribute", "value"]].round(4).to_string(index=False))
    print(sens_df[sens_df.metric.isin(["gamma_best_single", "gamma_best_pair"])]
          [["variant_type", "variant_value", "metric", "attribute", "value"]]
          .round(4).to_string(index=False))
    print(f"total wall {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()

"""Export the small JSON files the documentation site reads.

Everything is derived from results/ and from the primary spectra, so the site
cannot show a number that the manuscript does not. Output goes to docs/data/
and totals a few hundred kilobytes, which keeps the site archivable inside the
repository and inside the Zenodo record.

Run from the deposit root:  python code/20_export_web.py
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plastics_data as pdata  # noqa: E402

ROOT = pdata.ROOT
RES = os.path.join(ROOT, "results")
OUT = os.path.join(ROOT, "docs", "data")
POLY = ["HDPE", "LDPE", "PET", "PP", "PS", "PVC"]
N_SHOW = 12          # spectra kept per polymer for the explorer
N_POINTS = 600       # points per spectrum after decimation


def write(name, obj):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, separators=(",", ":"), allow_nan=False)
    return os.path.getsize(path)


def clean(v):
    if isinstance(v, float) and (np.isnan(v) or np.isinf(v)):
        return None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return None if not np.isfinite(v) else round(float(v), 6)
    return v


def frame(df):
    return {"columns": list(df.columns),
            "rows": [[clean(v) for v in row] for row in df.itertuples(index=False)]}


def spectra():
    grid, Y, meta = pdata.load("c8")
    step = max(1, len(grid) // N_POINTS)
    g = grid[::step]
    rng = np.random.default_rng(pdata.SEED)
    out = {"grid": [round(float(x), 2) for x in g], "windows": pdata.WINDOWS, "polymers": {}}
    for p in POLY:
        idx = np.where(meta.Polymer.to_numpy().astype(str) == p)[0]
        pick = rng.choice(idx, min(N_SHOW, len(idx)), replace=False)
        block = Y[:, ::step]
        out["polymers"][p] = {
            "mean": [round(float(v), 4) for v in block[idx].mean(axis=0)],
            "p05": [round(float(v), 4) for v in np.percentile(block[idx], 5, axis=0)],
            "p95": [round(float(v), 4) for v in np.percentile(block[idx], 95, axis=0)],
            "samples": [[round(float(v), 4) for v in block[i]] for i in pick]}
    return out


def main():
    sizes = {}

    sizes["spectra.json"] = write("spectra.json", spectra())

    par = pd.read_csv(os.path.join(RES, "cv_parity_folds.csv"))
    par["ci_high"] = par.ci_high.clip(upper=1.0)
    sizes["parity.json"] = write("parity.json", frame(par))

    rs = pd.read_csv(os.path.join(RES, "rough_sets.csv"))
    rs = rs[rs.variant == "sens.min bins, k=4"]
    gamma = (rs.groupby("attribute_set")
               .agg(gamma=("gamma", "first"), n_cells=("n_cells", "first"),
                    n_attributes=("n_attributes", "first"))
               .sort_values("gamma", ascending=False).reset_index())
    sizes["rough.json"] = write("rough.json", {
        "by_set": frame(gamma.head(20)),
        "by_polymer": frame(rs[["attribute_set", "polymer", "alpha", "boundary_size"]])})

    rm = pd.read_csv(os.path.join(RES, "rule_matrix.csv"))
    keep = ["body", "attribute", "conf_change", "conf_invariant", "wilson_low",
            "wilson_high", "n_pairs", "q_fdr", "significant_q05"]
    sizes["rules.json"] = write("rules.json", frame(rm[keep]))

    ov = pd.read_csv(os.path.join(RES, "ovl_matrix.csv"))
    sizes["ovl.json"] = write("ovl.json", frame(ov))

    lc = pd.read_csv(os.path.join(RES, "learning_curves.csv"))
    sizes["learning.json"] = write("learning.json", frame(lc))

    oc = pd.read_csv(os.path.join(RES, "oc_curves.csv"))
    sizes["oc.json"] = write("oc.json", frame(oc))

    st = pd.read_csv(os.path.join(RES, "shap_stability.csv"))
    rk = pd.read_csv(os.path.join(RES, "shap_ranks_per_repeat.csv"))
    sizes["stability.json"] = write("stability.json",
                                    {"summary": frame(st), "ranks": frame(rk)})

    ent = pd.read_csv(os.path.join(RES, "entry_condition_primary.csv"))
    sizes["entry.json"] = write("entry.json", frame(ent))

    sens = pd.read_csv(os.path.join(RES, "pair_cap_sensitivity.csv"))
    sizes["sensitivity.json"] = write("sensitivity.json", frame(sens))

    sel = pd.read_csv(os.path.join(RES, "selective_comparison.csv"))
    par = pd.read_csv(os.path.join(RES, "pareto_pairs.csv"))
    binn = pd.read_csv(os.path.join(RES, "binning_out_of_fold.csv"))
    acc = pd.read_csv(os.path.join(RES, "accepted_per_class.csv"))
    sizes["selective.json"] = write("selective.json",
                                    {"matched": frame(sel), "pairs": frame(par),
                                     "binning": frame(binn), "per_class": frame(acc)})

    total = sum(sizes.values())
    for k, v in sorted(sizes.items()):
        print(f"{k:20s} {v/1024:8.1f} kB")
    print(f"{'total':20s} {total/1024:8.1f} kB")


if __name__ == "__main__":
    main()

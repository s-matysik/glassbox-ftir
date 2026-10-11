"""Task 2: learning curves and the crossover point.

Training subsets of 1, 2, 5, 10, 25, 50, 75 and 100 percent of each training
fold, stratified by polymer, nested so that a smaller subset is contained in
every larger one. Template, sens.min, the selected attributes and the bin
edges are re-estimated on every subset; the baselines are re-tuned on every
subset. Test folds are untouched and identical to 10_cv_parity.py.

Crossover (PROTOCOL_v1.1 Section 6): the smallest training size at which the
95 percent interval of the difference in balanced accuracy between the glass
box and the best black box excludes zero in favour of the glass box. The
interval is a bootstrap over specimens, 1000 iterations, of the per-specimen
difference in seed-averaged correctness.

Run from the deposit root:  python code/11_learning_curves.py
Writes results/learning_curves.csv.
"""
from __future__ import annotations

# torch must be imported before scipy/sklearn in a worker process; see
# the note in 10_cv_parity.py and shard.py.
import torch  # noqa: E402,F401  (import order is load bearing)

import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ara_aco as aa  # noqa: E402
import baselines as bl  # noqa: E402
import glassbox as gb  # noqa: E402
import plastics_data as pdata  # noqa: E402
import shard  # noqa: E402

SIZES = [float(v) for v in os.environ.get(
    "LC_SIZES", "0.01,0.02,0.05,0.10,0.25,0.50,0.75,1.00").split(",")]
N_SEEDS = int(os.environ.get("LC_SEEDS", 5))
N_FOLDS = 5
N_JOBS = int(os.environ.get("N_JOBS", 5))
N_BOOT = 1000
K_SENS = 4.0
N_MC = 100
BASELINES = os.environ.get("LC_BASELINES", "svm_rbf,random_forest,plsda,cnn1d").split(",")
LOG = os.path.join(pdata.ROOT, "logs", "11_learning_curves.log")


def log(msg: str) -> None:
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(f"{time.strftime('%H:%M:%S')} {msg}\n")


def nested_subsets(y_tr, rng):
    """One permutation per class, drawn once; every subset size is a prefix of
    it, so subsets are nested and no generator is created inside the loop."""
    order = {}
    for c in np.unique(y_tr):
        idx = np.where(y_tr == c)[0]
        order[c] = idx[rng.permutation(len(idx))]
    out = {}
    for frac in SIZES:
        take = []
        for c, idx in order.items():
            n = max(int(round(frac * len(idx))), 2)
            take.append(idx[:min(n, len(idx))])
        out[frac] = np.sort(np.concatenate(take))
    return out


def run_fold(seed, fold, tr, te, Yraw, Yp, sigma, y, masks, names):
    torch.set_num_threads(1)
    rng = np.random.default_rng(pdata.SEED + 977 * seed + fold)
    subsets = nested_subsets(y[tr], rng)
    rows = []
    for frac, sub in subsets.items():
        t0 = time.time()
        idx = tr[sub]
        model = gb.fit(Yp[idx], sigma[idx], y[idx], masks, names, k=K_SENS,
                       n_mc=N_MC, seed=pdata.SEED + seed)
        dec, forced = model.predict(Yp[te])
        rec = {"seed": seed, "fold": fold, "size_frac": frac,
               "n_train": len(idx),
               "n_per_class": int(np.bincount(np.unique(y[idx], return_inverse=True)[1]).min()),
               "test_idx": te, "glass_decision": dec, "glass_forced": forced,
               "attributes": "+".join(model.attributes), "gamma": model.model.gamma}
        for name in BASELINES:
            pred, params = bl.METHODS[name](Yraw[idx], y[idx], Yraw[te], pdata.SEED + seed)
            rec[name] = np.asarray(pred).astype(str)
            rec[f"params_{name}"] = str(params)
        rec["secs"] = time.time() - t0
        rows.append(rec)
        log(f"seed {seed} fold {fold} size {frac:.2f} n={len(idx)} "
            f"attrs={rec['attributes']} {rec['secs']:.0f}s")
    return rows


def build():
    grid, Yraw, meta = pdata.load("c8")
    y = meta["Polymer"].to_numpy().astype(str)
    Yp = aa.preprocess(Yraw)
    sigma = aa.instrument_noise(Yp)
    masks = aa.window_masks(grid, pdata.WINDOWS)
    names = [f"{p}_{w}" for w in pdata.WINDOWS for p in ("A", "C", "O")]
    jobs = []
    for s in range(N_SEEDS):
        skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=pdata.SEED + s)
        for f, (tr, te) in enumerate(skf.split(np.zeros(len(y)), y)):
            jobs.append((s, f, tr, te))
    return Yraw, Yp, sigma, y, masks, names, jobs


def worker():
    idx, count = shard.worker_id()
    Yraw, Yp, sigma, y, masks, names, jobs = build()
    out = []
    for i, (s, f, tr, te) in shard.my_jobs(jobs, idx, count):
        out.append((i, run_fold(s, f, tr, te, Yraw, Yp, sigma, y, masks, names)))
    shard.dump("lc", idx, out)


def boot_diff(corr_a, corr_b, y, classes, rng, n_boot=N_BOOT):
    """Bootstrap over specimens of the balanced-accuracy difference a - b."""
    idx_all = np.arange(len(y))
    out = np.empty(n_boot)
    for b in range(n_boot):
        draw = rng.choice(idx_all, len(idx_all), replace=True)
        yd = y[draw]
        va, vb = [], []
        for c in classes:
            m = yd == c
            if m.any():
                va.append(corr_a[draw][m].mean())
                vb.append(corr_b[draw][m].mean())
        out[b] = np.mean(va) - np.mean(vb)
    return float(out.mean()), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


def main():
    Yraw, Yp, sigma, y, masks, names, jobs = build()
    classes = sorted(set(y))
    n = len(y)
    log(f"starting {len(jobs)} fold-runs x {len(SIZES)} sizes with {N_JOBS} workers")
    t0 = time.time()
    res = shard.run(os.path.abspath(__file__), "lc", N_JOBS)
    wall = time.time() - t0
    log(f"done in {wall/60:.1f} min")

    methods = ["glassbox_forced"] + BASELINES
    pred = {f: {m: np.empty((N_SEEDS, n), dtype=object) for m in methods} for f in SIZES}
    decision = {f: np.empty((N_SEEDS, n), dtype=object) for f in SIZES}
    meta_rows = []
    for fold_rows in res:
        for r in fold_rows:
            f, s, te = r["size_frac"], r["seed"], r["test_idx"]
            pred[f]["glassbox_forced"][s, te] = r["glass_forced"]
            decision[f][s, te] = r["glass_decision"]
            for m in BASELINES:
                pred[f][m][s, te] = r[m]
            meta_rows.append({"size_frac": f, "seed": s, "fold": r["fold"],
                              "n_train": r["n_train"], "n_per_class": r["n_per_class"],
                              "attributes": r["attributes"], "gamma": r["gamma"],
                              "secs": r["secs"]})
    md = pd.DataFrame(meta_rows)

    rng = np.random.default_rng(pdata.SEED + 3)
    rows = []
    for f in SIZES:
        corr = {m: (pred[f][m] == y[None, :]).mean(axis=0) for m in methods}
        bacc = {}
        for m in methods:
            bacc[m] = float(np.mean([gb.balanced_accuracy(y, pred[f][m][s], classes)
                                     for s in range(N_SEEDS)]))
        best_bb = max(BASELINES, key=lambda m: bacc[m])
        g = md[md.size_frac == f]
        cov = float(np.mean([(decision[f][s] != "refer").mean() for s in range(N_SEEDS)]))
        for m in methods:
            mean_b, lo, hi = boot_diff(corr[m], np.zeros(n), y, classes, rng)
            row = {"size_frac": f, "size_pct": 100 * f,
                   "n_train_total": int(g["n_train"].mean()),
                   "n_per_class": int(g["n_per_class"].mean()),
                   "method": m, "balanced_accuracy": bacc[m],
                   "bacc_ci_low": lo, "bacc_ci_high": hi,
                   "coverage": cov if m == "glassbox_forced" else 1.0,
                   "best_blackbox": best_bb,
                   "diff_vs_best_blackbox": np.nan,
                   "diff_ci_low": np.nan, "diff_ci_high": np.nan,
                   "glassbox_favoured": False}
            if m == "glassbox_forced":
                d, dlo, dhi = boot_diff(corr[m], corr[best_bb], y, classes, rng)
                row.update({"diff_vs_best_blackbox": d, "diff_ci_low": dlo,
                            "diff_ci_high": dhi, "glassbox_favoured": bool(dlo > 0)})
            rows.append(row)
    curves = pd.DataFrame(rows)
    pdata.write_csv(curves, "learning_curves.csv")

    gl = curves[curves.method == "glassbox_forced"].sort_values("size_frac")
    cross = gl[gl.glassbox_favoured]
    print(f"{N_SEEDS} seeds x {N_FOLDS} folds x {len(SIZES)} sizes, wall {wall/60:.1f} min")
    print(curves.pivot(index="size_pct", columns="method",
                       values="balanced_accuracy").round(4).to_string())
    print("glass box minus best black box, 95% bootstrap interval:")
    print(gl[["size_pct", "n_per_class", "best_blackbox", "diff_vs_best_blackbox",
              "diff_ci_low", "diff_ci_high"]].round(4).to_string(index=False))
    if len(cross):
        r = cross.iloc[0]
        print(f"crossover at {r.size_pct:.0f} percent "
              f"({r.n_train_total} training specimens, {r.n_per_class} per class)")
    else:
        print("no crossover: the interval of the difference never excludes zero "
              "in favour of the glass box at any training size")


if __name__ == "__main__":
    if shard.is_worker():
        worker()
    else:
        main()

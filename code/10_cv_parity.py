"""Task 1 of the IEEE Access computational core: accuracy parity.

Stratified 5-fold cross-validation of the glass box against four black boxes
on the primary dataset (FTIR-Plastics C8, 3000 specimens, six polymers),
repeated over N_SEEDS seeds. Protocol: PROTOCOL_v1.1.md Sections 5 and 6.

Everything that the glass box estimates -- template, sens.min, the at most two
condition attributes and the bin edges -- is estimated inside the training
fold. The baselines are tuned inside the same training fold. No quantity is
ever read from a test fold.

Run from the deposit root:  python code/10_cv_parity.py

Writes results/cv_parity.csv, results/cv_per_class.csv,
results/cv_fold_scores.csv and logs/cv_predictions.npz (the per-specimen
out-of-fold predictions reused by 13_shap_oc.py).
"""
from __future__ import annotations

# torch must be imported before scipy/sklearn. In a worker process spawned by
# another Python process (see shard.py) the reverse order makes torch's
# shm.dll fail to resolve a dependency: OSError WinError 127. Importing it
# first is the only change this costs.
import torch  # noqa: E402,F401  (import order is load bearing)

import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ara_aco as aa  # noqa: E402
import baselines as bl  # noqa: E402
import glassbox as gb  # noqa: E402
import plastics_data as pdata  # noqa: E402
import shard  # noqa: E402

N_SEEDS = int(os.environ.get("N_SEEDS", 20))
N_FOLDS = 5
N_BOOT = 1000
N_JOBS = int(os.environ.get("N_JOBS", 5))
K_SENS = 4.0
N_MC = 100
RESOLUTION = os.environ.get("RESOLUTION", "c8")
SUFFIX = "" if RESOLUTION == "c8" else f"_{RESOLUTION}"

BASELINES = ["svm_rbf", "random_forest", "plsda", "cnn1d"]
LOG = os.path.join(pdata.ROOT, "logs", f"10_cv_parity{SUFFIX}.log")


def log(msg: str) -> None:
    os.makedirs(os.path.dirname(LOG), exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as fh:
        fh.write(f"{time.strftime('%H:%M:%S')} {msg}\n")


def run_fold(seed, fold, tr, te, Yraw, Yp, sigma, y, masks, names):
    torch.set_num_threads(1)
    t0 = time.time()
    model = gb.fit(Yp[tr], sigma[tr], y[tr], masks, names, k=K_SENS, n_mc=N_MC,
                   seed=pdata.SEED + seed)
    dec, forced = model.predict(Yp[te])
    out = {"seed": seed, "fold": fold, "test_idx": te,
           "glass_decision": dec, "glass_forced": forced,
           "attributes": "+".join(model.attributes),
           "gamma": model.model.gamma,
           "inner_score": model.inner_score,
           "n_cells": len(model.model.cell_class),
           "params": {}}
    for name in BASELINES:
        t1 = time.time()
        pred, params = bl.METHODS[name](Yraw[tr], y[tr], Yraw[te], pdata.SEED + seed)
        out[name] = np.asarray(pred).astype(str)
        out["params"][name] = str(params)
        out[f"secs_{name}"] = time.time() - t1
    out["secs_total"] = time.time() - t0
    log(f"seed {seed} fold {fold} done in {out['secs_total']:.0f}s "
        f"attrs={out['attributes']} gamma={out['gamma']:.3f}")
    return out


def bacc(y_true, y_pred, classes):
    return gb.balanced_accuracy(y_true, y_pred, classes)


def bootstrap_ci(correct_by_class, y, classes, rng, n_boot=N_BOOT, subset=None):
    """Bootstrap over specimens: resample the specimens with replacement and
    recompute the seed-averaged balanced accuracy on the resample. ``subset``
    restricts the resampled population, which is what the referral-respecting
    interval needs: only the specimens the glass box actually decided on."""
    idx_all = np.arange(len(y)) if subset is None else np.asarray(subset)
    stats = np.empty(n_boot)
    for b in range(n_boot):
        draw = rng.choice(idx_all, len(idx_all), replace=True)
        yd = y[draw]
        vals = []
        for c in classes:
            m = yd == c
            if m.any():
                vals.append(correct_by_class[draw][m].mean())
        stats[b] = np.mean(vals)
    return float(stats.mean()), float(np.percentile(stats, 2.5)), float(np.percentile(stats, 97.5))


def nadeau_bengio(diffs, n_train, n_test):
    """Corrected resampled t-test (Nadeau and Bengio 2003)."""
    from scipy import stats as st
    d = np.asarray(diffs, float)
    n = len(d)
    var = d.var(ddof=1)
    if var == 0:
        return float(d.mean()), np.inf if d.mean() > 0 else -np.inf, 0.0
    t = d.mean() / np.sqrt(var * (1.0 / n + n_test / n_train))
    p = 2 * (1 - st.t.cdf(abs(t), df=n - 1))
    return float(d.mean()), float(t), float(p)


def build(resolution):
    grid, Yraw, meta = pdata.load(resolution)
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
    Yraw, Yp, sigma, y, masks, names, jobs = build(RESOLUTION)
    out = []
    for i, (s, f, tr, te) in shard.my_jobs(jobs, idx, count):
        out.append((i, run_fold(s, f, tr, te, Yraw, Yp, sigma, y, masks, names)))
    shard.dump(f"cv{SUFFIX}", idx, out)


def main():
    Yraw, Yp, sigma, y, masks, names, jobs = build(RESOLUTION)
    classes = sorted(set(y))
    n = len(y)
    log(f"starting {len(jobs)} fold-runs on {RESOLUTION} with {N_JOBS} workers")
    t0 = time.time()
    res = shard.run(os.path.abspath(__file__), f"cv{SUFFIX}", N_JOBS)
    wall = time.time() - t0
    log(f"all folds done in {wall/60:.1f} min")

    methods = ["glassbox_forced"] + BASELINES
    pred = {m: np.empty((N_SEEDS, n), dtype=object) for m in methods}
    decision = np.empty((N_SEEDS, n), dtype=object)
    fold_rows = []
    for r in res:
        te = r["test_idx"]
        pred["glassbox_forced"][r["seed"], te] = r["glass_forced"]
        decision[r["seed"], te] = r["glass_decision"]
        for m in BASELINES:
            pred[m][r["seed"], te] = r[m]
        row = {"seed": r["seed"], "fold": r["fold"], "attributes": r["attributes"],
               "gamma_train": r["gamma"], "inner_bacc": r["inner_score"],
               "n_cells": r["n_cells"], "secs_total": r["secs_total"]}
        cov = r["glass_decision"] != "refer"
        row["coverage"] = float(cov.mean())
        row["bacc_glassbox_forced"] = bacc(y[te], r["glass_forced"], classes)
        row["bacc_glassbox_covered"] = bacc(y[te][cov], r["glass_forced"][cov], classes) if cov.any() else np.nan
        for m in BASELINES:
            row[f"bacc_{m}"] = bacc(y[te], r[m], classes)
            row[f"params_{m}"] = r["params"][m]
        fold_rows.append(row)
    folds = pd.DataFrame(fold_rows).sort_values(["seed", "fold"]).reset_index(drop=True)
    pdata.write_csv(folds, f"cv_fold_scores{SUFFIX}.csv")

    np.savez_compressed(os.path.join(pdata.ROOT, "logs", f"cv_predictions{SUFFIX}.npz"),
                        y=y, decision=decision.astype(str),
                        **{m: pred[m].astype(str) for m in methods})

    rng = np.random.default_rng(pdata.SEED)
    rows = []
    n_train = n * (N_FOLDS - 1) / N_FOLDS
    n_test = n / N_FOLDS
    best_base = max(BASELINES, key=lambda m: folds[f"bacc_{m}"].mean())
    for m in methods:
        corr = (pred[m] == y[None, :]).mean(axis=0)
        mean_b, lo, hi = bootstrap_ci(corr, y, classes, rng)
        mcc = float(np.mean([matthews_corrcoef(y, pred[m][s]) for s in range(N_SEEDS)]))
        acc = float(np.mean([(pred[m][s] == y).mean() for s in range(N_SEEDS)]))
        if m == "glassbox_forced":
            d_mean, t, p = 0.0, np.nan, np.nan
        else:
            d_mean, t, p = nadeau_bengio(folds["bacc_glassbox_forced"] - folds[f"bacc_{m}"],
                                         n_train, n_test)
        rows.append({"method": m, "mode": "forced choice",
                     "balanced_accuracy": folds[f"bacc_{m}"].mean(),
                     "bacc_boot_mean": mean_b, "bacc_ci_low": lo, "bacc_ci_high": hi,
                     "mcc": mcc, "accuracy": acc, "coverage": 1.0,
                     "diff_vs_glassbox": -d_mean,
                     "nb_t": (-t if np.isfinite(t) else t) if m != "glassbox_forced" else np.nan,
                     "nb_p": p, "is_best_baseline": m == best_base})
    cov_mask = decision != "refer"
    corr_cov = np.where(cov_mask, pred["glassbox_forced"] == y[None, :], np.nan)
    seen = cov_mask.any(axis=0)
    corr_cov_mean = np.zeros(n)
    corr_cov_mean[seen] = np.nanmean(corr_cov[:, seen], axis=0)
    mean_b, lo, hi = bootstrap_ci(corr_cov_mean, y, classes, rng,
                                  subset=np.where(seen)[0])
    rows.insert(1, {"method": "glassbox_three_way", "mode": "referral respecting",
                    "balanced_accuracy": folds["bacc_glassbox_covered"].mean(),
                    "bacc_boot_mean": mean_b, "bacc_ci_low": lo, "bacc_ci_high": hi,
                    "mcc": float(np.mean([matthews_corrcoef(y[cov_mask[s]], pred["glassbox_forced"][s][cov_mask[s]])
                                          for s in range(N_SEEDS)])),
                    "accuracy": float(np.mean([(pred["glassbox_forced"][s][cov_mask[s]] == y[cov_mask[s]]).mean()
                                               for s in range(N_SEEDS)])),
                    "coverage": folds["coverage"].mean(),
                    "diff_vs_glassbox": np.nan, "nb_t": np.nan, "nb_p": np.nan,
                    "is_best_baseline": False})
    rows.append({"method": "majority_class", "mode": "reference",
                 "balanced_accuracy": 1.0 / len(classes), "bacc_boot_mean": np.nan,
                 "bacc_ci_low": np.nan, "bacc_ci_high": np.nan, "mcc": 0.0,
                 "accuracy": float(max((y == c).mean() for c in classes)),
                 "coverage": 1.0, "diff_vs_glassbox": np.nan, "nb_t": np.nan,
                 "nb_p": np.nan, "is_best_baseline": False})
    parity = pd.DataFrame(rows)
    pdata.write_csv(parity, f"cv_parity{SUFFIX}.csv")

    per_class = []
    for m in methods:
        for c in classes:
            mask = y == c
            sens = float(np.mean([(pred[m][s][mask] == c).mean() for s in range(N_SEEDS)]))
            spec = float(np.mean([(pred[m][s][~mask] != c).mean() for s in range(N_SEEDS)]))
            corr = (pred[m][:, mask] == c).mean(axis=0)
            boot = np.empty(N_BOOT)
            idx = np.arange(mask.sum())
            for b in range(N_BOOT):
                boot[b] = corr[rng.choice(idx, len(idx), replace=True)].mean()
            per_class.append({"method": m, "polymer": c, "n_specimens": int(mask.sum()),
                              "sensitivity": sens, "sens_ci_low": float(np.percentile(boot, 2.5)),
                              "sens_ci_high": float(np.percentile(boot, 97.5)),
                              "specificity": spec})
    for c in classes:
        mask = (y == c)
        cov_c = cov_mask[:, mask]
        sens = float(np.mean([(pred["glassbox_forced"][s][mask][cov_c[s]] == c).mean()
                              if cov_c[s].any() else np.nan for s in range(N_SEEDS)]))
        per_class.append({"method": "glassbox_three_way", "polymer": c,
                          "n_specimens": int(mask.sum()), "sensitivity": sens,
                          "sens_ci_low": np.nan, "sens_ci_high": np.nan,
                          "specificity": float(np.mean([
                              (pred["glassbox_forced"][s][~mask][cov_mask[s][~mask]] != c).mean()
                              for s in range(N_SEEDS)]))})
    pdata.write_csv(pd.DataFrame(per_class), f"cv_per_class{SUFFIX}.csv")

    attr_counts = folds["attributes"].value_counts()
    print(f"dataset {RESOLUTION}: {n} specimens, {len(classes)} polymers, "
          f"{N_SEEDS} seeds x {N_FOLDS} folds, wall {wall/60:.1f} min")
    print(parity[["method", "balanced_accuracy", "bacc_ci_low", "bacc_ci_high",
                  "mcc", "coverage", "nb_p"]].to_string(index=False))
    print(f"best baseline: {best_base}")
    print(f"glass-box attribute sets chosen (top 5 of {attr_counts.size}):")
    print(attr_counts.head(5).to_string())
    print(f"mean gamma on training folds: {folds['gamma_train'].mean():.4f}")


if __name__ == "__main__":
    if shard.is_worker():
        worker()
    else:
        main()

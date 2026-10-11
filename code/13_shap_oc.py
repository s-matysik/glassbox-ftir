"""Task 4: attribution stability, and acceptance sampling under an imperfect test.

Part 1 contrasts how stable a post-hoc attribution is against how stable the
glass box's own selected attributes are, over repeated fits on resampled
training data. Both live in the same nine-attribute ACO space, so the
comparison is like for like rather than nine attributes against 1869 raw
channels.

Attribution is computed with TreeSHAP when the `shap` package is present and
otherwise by exact interventional Shapley values over the nine attributes
(512 subsets, marginalising the complement against a background sample). The
path actually used is recorded in the output table.

Part 2 derives operating-characteristic curves for single-stage acceptance
sampling, first for a perfect inspector and then for an inspector with the
sensitivity and specificity measured in the cross-validation.

Run from the deposit root:  python code/13_shap_oc.py
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plastics_data as pdata
from ara_aco import (aco_matrix, calibrate_sens_min, fit_rough, instrument_noise,
                     make_template, preprocess, window_masks)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
N_REPEATS = int(os.environ.get("SHAP_REPEATS", 20))
N_TEST = 60
N_BG = 30
SEED = 17
WINDOWS = {"W1": (2800, 3000), "W2": (1650, 1800), "W3": (1350, 1500)}


def exact_shapley(model, X_test, X_bg, n_features):
    """Interventional Shapley values, exact over all 2^n subsets."""
    subsets = [frozenset(c) for r in range(n_features + 1)
               for c in itertools.combinations(range(n_features), r)]
    index = {s: i for i, s in enumerate(subsets)}
    rows, keep = [], []
    for s in subsets:
        mask = np.zeros(n_features, bool)
        mask[list(s)] = True
        X = np.repeat(X_bg[None, :, :], len(X_test), axis=0).copy()
        X[:, :, mask] = X_test[:, None, mask]
        rows.append(X.reshape(-1, n_features))
        keep.append(len(X_test) * len(X_bg))
    big = np.vstack(rows)
    proba = model.predict_proba(big).max(axis=1)
    v = {}
    off = 0
    for s, n in zip(subsets, keep):
        v[s] = proba[off:off + n].reshape(len(X_test), len(X_bg)).mean(axis=1)
        off += n
    phi = np.zeros((len(X_test), n_features))
    fact = {k: np.math.factorial(k) for k in range(n_features + 1)}
    for j in range(n_features):
        for s in subsets:
            if j in s:
                continue
            k = len(s)
            w = fact[k] * fact[n_features - k - 1] / fact[n_features]
            phi[:, j] += w * (v[s | {j}] - v[s])
    return phi


def main():
    os.makedirs(RES, exist_ok=True)
    rng = np.random.default_rng(SEED)
    grid, Y, meta = pdata.load("c8")
    y = meta["Polymer"].to_numpy().astype(str)
    Y = preprocess(Y)
    masks = window_masks(grid, WINDOWS)
    sigma = instrument_noise(Y)

    try:
        import shap
        explainer_kind = f"TreeSHAP (shap {shap.__version__})"
    except ImportError:
        shap = None
        explainer_kind = "exact interventional Shapley, 512 subsets"

    classes = sorted(set(y))
    shap_ranks, glass_sets = [], []
    for rep in range(N_REPEATS):
        idx = rng.choice(len(Y), int(0.7 * len(Y)), replace=False)
        hold = np.setdiff1d(np.arange(len(Y)), idx)
        template = make_template(Y[idx], masks)
        F, names = aco_matrix(Y, masks, template)
        sens = calibrate_sens_min(Y[idx][:200], sigma[idx][:200], y[idx][:200],
                                  masks, template, names, k=4.0, n_mc=30,
                                  seed=int(rng.integers(1e9)))
        rf = RandomForestClassifier(n_estimators=300, random_state=int(rng.integers(1e9)),
                                    n_jobs=1).fit(F[idx], y[idx])
        test = rng.choice(hold, min(N_TEST, len(hold)), replace=False)
        bg = rng.choice(idx, N_BG, replace=False)
        if shap is not None:
            values = shap.TreeExplainer(rf).shap_values(F[test])
            arr = np.abs(np.array(values))
            # Reduce every axis except the feature axis. Multiclass TreeSHAP
            # returns (samples, features, classes) in shap 0.4x and later, and a
            # list of per-class arrays in older versions; averaging over the
            # wrong axis would rank classes instead of attributes.
            feat_axis = [i for i, s in enumerate(arr.shape) if s == len(names)]
            if len(feat_axis) != 1:
                raise RuntimeError(f"cannot locate feature axis in SHAP output {arr.shape}")
            imp = arr.mean(axis=tuple(i for i in range(arr.ndim) if i != feat_axis[0]))
        else:
            imp = np.abs(exact_shapley(rf, F[test], F[bg], len(names))).mean(axis=0)
        shap_ranks.append(stats.rankdata(-np.asarray(imp).ravel()))

        best, best_gamma = None, -1.0
        for a in names:
            for b in names:
                if a >= b:
                    continue
                m = fit_rough(F[idx], names, y[idx].astype(str), [a, b], sens)
                if m.gamma > best_gamma:
                    best, best_gamma = (a, b), m.gamma
        glass_sets.append(frozenset(best))

    Rk = np.array(shap_ranks)
    taus = [stats.kendalltau(Rk[i], Rk[j]).statistic
            for i in range(len(Rk)) for j in range(i + 1, len(Rk))]
    top3 = [frozenset(np.argsort(r)[:3]) for r in Rk]
    modal_top3 = max(set(top3), key=top3.count)
    modal_glass = max(set(glass_sets), key=glass_sets.count)

    stab = pd.DataFrame([
        {"representation": "black box, post-hoc attribution", "method": explainer_kind,
         "repeats": N_REPEATS, "kendall_tau_mean": float(np.mean(taus)),
         "kendall_tau_min": float(np.min(taus)),
         "modal_set_share": top3.count(modal_top3) / len(top3),
         "n_distinct_sets": len(set(top3)),
         "modal_set": ", ".join(sorted(np.array(pdata.ATTR_NAMES)[list(modal_top3)]))
         if hasattr(pdata, "ATTR_NAMES") else str(sorted(modal_top3))},
        {"representation": "glass box, selected attributes", "method": "rough-set quality",
         "repeats": N_REPEATS, "kendall_tau_mean": np.nan, "kendall_tau_min": np.nan,
         "modal_set_share": glass_sets.count(modal_glass) / len(glass_sets),
         "n_distinct_sets": len(set(glass_sets)),
         "modal_set": ", ".join(sorted(modal_glass))}])
    stab.to_csv(os.path.join(RES, "shap_stability.csv"), index=False, encoding="utf-8-sig")

    per_repeat = pd.DataFrame(Rk, columns=names)
    per_repeat.insert(0, "repeat", np.arange(1, len(Rk) + 1))
    per_repeat["glassbox_set"] = [", ".join(sorted(s)) for s in glass_sets]
    per_repeat.to_csv(os.path.join(RES, "shap_ranks_per_repeat.csv"), index=False,
                      encoding="utf-8-sig")

    cv = pd.read_csv(os.path.join(RES, "cv_per_class.csv"))
    se = float(cv["sensitivity"].mean()) if "sensitivity" in cv else np.nan
    sp = se
    rows = []
    for n, c in [(50, 0), (80, 1), (125, 2), (200, 3)]:
        for p in np.round(np.arange(0.005, 0.155, 0.005), 3):
            for label, (s_e, s_p) in {"perfect inspector": (1.0, 1.0),
                                      "observed inspector": (se, sp)}.items():
                p_eff = p * s_e + (1 - p) * (1 - s_p)
                pa = float(stats.binom.cdf(c, n, p_eff))
                rows.append({"n": n, "c": c, "lot_fraction_defective": p,
                             "inspector": label, "sensitivity": s_e, "specificity": s_p,
                             "p_effective": p_eff, "prob_accept": pa})
    oc = pd.DataFrame(rows)
    oc.to_csv(os.path.join(RES, "oc_curves.csv"), index=False, encoding="utf-8-sig")

    risk = []
    for (n, c), g in oc.groupby(["n", "c"]):
        for label in g.inspector.unique():
            gg = g[g.inspector == label]
            prod = 1 - float(gg.loc[gg.lot_fraction_defective == 0.01, "prob_accept"].iloc[0])
            cons = float(gg.loc[gg.lot_fraction_defective == 0.10, "prob_accept"].iloc[0])
            risk.append({"n": n, "c": c, "inspector": label,
                         "producer_risk_at_1pct": prod, "consumer_risk_at_10pct": cons})
    RK = pd.DataFrame(risk)
    RK.to_csv(os.path.join(RES, "oc_risk.csv"), index=False, encoding="utf-8-sig")

    print(f"attribution path: {explainer_kind} | repeats {N_REPEATS}")
    print(stab.round(4).to_string(index=False))
    print(f"mean inspector sensitivity used for OC curves: {se:.4f}")
    print(RK.round(4).to_string(index=False))


if __name__ == "__main__":
    main()

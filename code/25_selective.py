"""Selective classification: the comparisons a reviewer asks for.

Three questions the main experiment leaves open.

1. The attribute pair the classifier selects is not the pair that maximises the
   quality of approximation on the full sample. Which is better out of fold,
   and what does each cost in referrals? This builds a small Pareto front over
   candidate pairs, evaluated with the same folds as the main run.

2. Accuracy on accepted specimens is not comparable with accuracy on all
   specimens. Two fair comparisons are made: every baseline restricted to
   exactly the specimens the glass box accepts, and every baseline given its
   own reject option tuned to the same coverage.

3. The discretisation rule changes the in-sample quality of approximation
   sharply. Does it change the out-of-fold accuracy and coverage as much?

Inputs : the primary dataset, logs/cv_predictions.npz and
         results/cv_fold_scores.csv from 10_cv_parity.py
Outputs: results/pareto_pairs.csv, results/selective_comparison.csv,
         results/accepted_per_class.csv, results/binning_out_of_fold.csv

Run from the deposit root:  python code/25_selective.py
"""
from __future__ import annotations

import ast
import os
import sys

import numpy as np
import pandas as pd
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score, matthews_corrcoef
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ara_aco as aa  # noqa: E402
import plastics_data as pdata  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
LOGS = os.path.join(ROOT, "logs")
N_FOLDS = 5
K = 4.0
N_MC = 100
# The pairs worth comparing: the one the classifier selects, the one the
# quality-of-approximation criterion selects, and the next best by either.
PAIRS = [["A_W2", "C_W3"], ["O_W1", "O_W2"], ["O_W1", "A_W3"], ["O_W1", "O_W3"],
         ["A_W1", "A_W3"], ["C_W3", "O_W2"]]


def folds(y):
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=pdata.SEED)
    return list(skf.split(np.zeros(len(y)), y))


def glass_fold(Yp, sigma, y, masks, names, tr, te, attrs, n_bins=None):
    """Fit the rule on the training fold with a FIXED attribute pair."""
    template = aa.make_template(Yp[tr], masks)
    F_tr, _ = aa.aco_matrix(Yp[tr], masks, template)
    F_te, _ = aa.aco_matrix(Yp[te], masks, template)
    D_tr, D_te = aa.delta(F_tr, names), aa.delta(F_te, names)
    sens = aa.calibrate_sens_min(Yp[tr], sigma[tr], y[tr], masks, template, names,
                                 k=K, n_mc=N_MC, seed=pdata.SEED)
    model = aa.fit_rough(D_tr, names, y[tr], attrs, sens, n_bins=n_bins)
    dec, forced = aa.three_way(model, D_te, names)
    return dec, forced, model.gamma


def selective_curve(score, correct, coverage):
    """Accuracy on the most confident fraction `coverage`."""
    n_keep = max(1, int(round(coverage * len(score))))
    keep = np.argsort(-np.asarray(score))[:n_keep]
    return float(np.mean(np.asarray(correct)[keep])), keep


def main():
    grid, Y, meta = pdata.load("c8")
    y = meta["Polymer"].to_numpy().astype(str)
    Yp = aa.preprocess(Y)
    sigma = aa.instrument_noise(Yp)
    masks = aa.window_masks(grid, pdata.WINDOWS)
    _, names = aa.aco_matrix(Yp[:2], masks, aa.make_template(Yp[:50], masks))
    fold_list = folds(y)
    classes = sorted(set(y))

    # ---- 1. Pareto front over candidate attribute pairs -------------------
    rows = []
    for attrs in PAIRS:
        dec_all = np.empty(len(y), dtype=object)
        forced_all = np.empty(len(y), dtype=object)
        gammas = []
        for tr, te in fold_list:
            dec, forced, g = glass_fold(Yp, sigma, y, masks, names, tr, te, attrs)
            dec_all[te], forced_all[te] = dec, forced
            gammas.append(g)
        acc = dec_all != "refer"
        rows.append({
            "attributes": "+".join(attrs),
            "coverage_out_of_fold": round(float(acc.mean()), 4),
            "bacc_covered": round(balanced_accuracy_score(y[acc], forced_all[acc]), 4),
            "bacc_forced": round(balanced_accuracy_score(y, forced_all), 4),
            "mcc_forced": round(matthews_corrcoef(y, forced_all.astype(str)), 4),
            "gamma_train_mean": round(float(np.mean(gammas)), 4)})
        print(f"pair {'+'.join(attrs)}: coverage {rows[-1]['coverage_out_of_fold']} "
              f"covered {rows[-1]['bacc_covered']} forced {rows[-1]['bacc_forced']}")
    pareto = pd.DataFrame(rows).sort_values("bacc_covered", ascending=False)
    pareto.to_csv(os.path.join(RES, "pareto_pairs.csv"), index=False, encoding="utf-8-sig")

    # ---- 2. Baselines on the same specimens, and with a reject option -----
    z = np.load(os.path.join(LOGS, "cv_predictions.npz"), allow_pickle=True)
    decision = np.asarray(z["decision"])[0].astype(str)
    accepted = decision != "refer"
    coverage = float(accepted.mean())
    fscores = pd.read_csv(os.path.join(RES, "cv_fold_scores.csv"))

    sel = []
    for method in ["glassbox_forced", "svm_rbf", "random_forest", "plsda", "cnn1d"]:
        pred = np.asarray(z[method])[0].astype(str)
        sel.append({"method": method, "regime": "all specimens",
                    "coverage": 1.0,
                    "balanced_accuracy": round(balanced_accuracy_score(y, pred), 4),
                    "mcc": round(matthews_corrcoef(y, pred), 4)})
        sel.append({"method": method, "regime": "specimens the rule accepts",
                    "coverage": round(coverage, 4),
                    "balanced_accuracy": round(
                        balanced_accuracy_score(y[accepted], pred[accepted]), 4),
                    "mcc": round(matthews_corrcoef(y[accepted], pred[accepted]), 4)})

    # Each baseline with its own reject option, abstaining on the least
    # confident specimens until it reaches the rule's coverage.
    oof_margin = {m: np.zeros(len(y)) for m in ["svm_rbf", "random_forest", "plsda"]}
    for (tr, te), (_, frow) in zip(fold_list, fscores.iterrows()):
        sc = StandardScaler().fit(Yp[tr])
        Xtr, Xte = sc.transform(Yp[tr]), sc.transform(Yp[te])
        par = ast.literal_eval(frow["params_svm_rbf"])
        svm = SVC(kernel="rbf", **par).fit(Xtr, y[tr])
        dfun = svm.decision_function(Xte)
        top2 = np.sort(dfun, axis=1)[:, -2:]
        oof_margin["svm_rbf"][te] = top2[:, 1] - top2[:, 0]
        rf = RandomForestClassifier(
            n_estimators=500, random_state=pdata.SEED, n_jobs=1).fit(Xtr, y[tr])
        p = rf.predict_proba(Xte)
        ps = np.sort(p, axis=1)
        oof_margin["random_forest"][te] = ps[:, -1] - ps[:, -2]
        ncomp = ast.literal_eval(frow["params_plsda"]).get("n_components", 10)
        Ytr = np.array([[1.0 if c == v else 0.0 for c in classes] for v in y[tr]])
        pls = PLSRegression(n_components=ncomp).fit(Xtr, Ytr)
        q = pls.predict(Xte)
        qs = np.sort(q, axis=1)
        oof_margin["plsda"][te] = qs[:, -1] - qs[:, -2]
        print(f"fold done, svm params {par}")

    for method, margin in oof_margin.items():
        pred = np.asarray(z[method])[0].astype(str)
        correct = (pred == y).astype(float)
        acc_at, keep = selective_curve(margin, correct, coverage)
        sel.append({"method": method, "regime": "own reject option at matched coverage",
                    "coverage": round(coverage, 4),
                    "balanced_accuracy": round(
                        balanced_accuracy_score(y[keep], pred[keep]), 4),
                    "mcc": round(matthews_corrcoef(y[keep], pred[keep]), 4)})
    gb = np.asarray(z["glassbox_forced"])[0].astype(str)
    sel.append({"method": "glassbox", "regime": "own reject option at matched coverage",
                "coverage": round(coverage, 4),
                "balanced_accuracy": round(
                    balanced_accuracy_score(y[accepted], gb[accepted]), 4),
                "mcc": round(matthews_corrcoef(y[accepted], gb[accepted]), 4)})
    SEL = pd.DataFrame(sel)
    SEL.to_csv(os.path.join(RES, "selective_comparison.csv"), index=False,
               encoding="utf-8-sig")
    print(SEL.to_string(index=False))

    # ---- 3. Accepted population per class ---------------------------------
    per = []
    for c in classes:
        m = y == c
        per.append({"polymer": c, "n": int(m.sum()),
                    "n_accepted": int((m & accepted).sum()),
                    "coverage": round(float((m & accepted).mean() / m.mean()), 4),
                    "recall_forced": round(float((gb[m] == c).mean()), 4),
                    "recall_accepted": round(
                        float((gb[m & accepted] == c).mean()), 4)})
    PER = pd.DataFrame(per)
    PER.to_csv(os.path.join(RES, "accepted_per_class.csv"), index=False,
               encoding="utf-8-sig")
    print(PER.to_string(index=False))

    # ---- 4. Does the binning rule move the out-of-fold result? ------------
    binrows = []
    for label, nb in [("sens.min width", None), ("10 equal width", 10)]:
        for attrs in (["A_W2", "C_W3"], ["O_W1", "O_W2"]):
            dec_all = np.empty(len(y), dtype=object)
            forced_all = np.empty(len(y), dtype=object)
            for tr, te in fold_list:
                dec, forced, _ = glass_fold(Yp, sigma, y, masks, names, tr, te,
                                            attrs, n_bins=nb)
                dec_all[te], forced_all[te] = dec, forced
            acc = dec_all != "refer"
            binrows.append({"binning": label, "attributes": "+".join(attrs),
                            "coverage_out_of_fold": round(float(acc.mean()), 4),
                            "bacc_covered": round(
                                balanced_accuracy_score(y[acc], forced_all[acc]), 4),
                            "bacc_forced": round(
                                balanced_accuracy_score(y, forced_all), 4)})
            print(f"binning {label} {'+'.join(attrs)} -> {binrows[-1]}")
    pd.DataFrame(binrows).to_csv(os.path.join(RES, "binning_out_of_fold.csv"),
                                 index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    main()

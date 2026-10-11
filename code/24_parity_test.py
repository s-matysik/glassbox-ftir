"""The significance test behind the parity claim, computed from the fold record.

The frozen protocol specifies the corrected resampled test of Nadeau and Bengio,
which needs repeated runs. One repetition was completed, so the test actually
applicable is a paired t-test over the five folds. This script computes it from
results/cv_fold_scores.csv and writes results/parity_test.csv, so the numbers in
the manuscript come from a file rather than from a one-off calculation.

Run from the deposit root:  python code/24_parity_test.py
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
BASELINES = ["svm_rbf", "plsda", "random_forest", "cnn1d"]


def paired(a, b):
    d = np.asarray(b, float) - np.asarray(a, float)
    n = len(d)
    se = d.std(ddof=1) / np.sqrt(n)
    t = d.mean() / se
    p = 2 * (1 - stats.t.cdf(abs(t), n - 1))
    half = stats.t.ppf(0.975, n - 1) * se
    return d.mean(), d.mean() - half, d.mean() + half, t, p, n


def main():
    fs = pd.read_csv(os.path.join(RES, "cv_fold_scores.csv"))
    rows = []
    for b in BASELINES:
        diff, lo, hi, t, p, n = paired(fs.bacc_glassbox_forced, fs["bacc_" + b])
        rows.append({"comparison": f"{b} minus glass box, forced choice", "n_folds": n,
                     "mean_difference": round(diff, 4), "ci_low": round(lo, 4),
                     "ci_high": round(hi, 4), "t": round(t, 2), "p": round(p, 5),
                     "test": "paired t-test across folds, one repetition completed"})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "parity_test.csv"), index=False, encoding="utf-8-sig")
    print(out.to_string(index=False))
    print(f"strongest baseline: {out.sort_values('mean_difference').iloc[-1].comparison}")


if __name__ == "__main__":
    main()

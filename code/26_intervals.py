"""Rebuild the fold-level summary with intervals that respect the metric range.

The first version of this table carried a Student interval on the raw scale.
Balanced accuracy is bounded above by 1 and the fold means sit within 0.002 of
that bound, so five of the six upper limits fell outside [0, 1] -- 1.0025 for
the convolutional network. The figure and document generators then clipped them
at 1, which hides an interval the data does not support rather than computing a
better one.

The interval here is computed on the logit scale and transformed back, with the
standard error carried across by the delta method. It cannot leave (0, 1), it is
asymmetric in the way a proportion near its bound should be, and it needs no
additional training: the fold scores are read from the record written by
code/10_cv_parity.py.

Input : results/cv_fold_scores.csv
Output: results/cv_parity_folds.csv

Run from the deposit root:  python code/26_intervals.py
"""
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")

ORDER = ["glassbox_forced", "glassbox_covered", "svm_rbf", "random_forest",
         "plsda", "cnn1d"]


def logit_interval(values, conf=0.95):
    """Two-sided interval for the mean of bounded scores, on the logit scale.

    The mean and its standard error are computed on the raw scale, then mapped
    through the logit with the delta-method factor 1 / (p (1 - p)). A mean of
    exactly 0 or 1 has no logit, and in that degenerate case the point estimate
    is returned as the interval, which is the honest answer for scores that
    never varied.
    """
    v = np.asarray(values, float)
    v = v[~np.isnan(v)]
    n = len(v)
    mean = v.mean()
    if n < 2 or mean <= 0.0 or mean >= 1.0:
        return mean, mean, mean, n
    se = v.std(ddof=1) / np.sqrt(n)
    if se == 0.0:
        return mean, mean, mean, n
    eta = np.log(mean / (1.0 - mean))
    se_eta = se / (mean * (1.0 - mean))
    half = stats.t.ppf(0.5 + conf / 2.0, n - 1) * se_eta
    lo = 1.0 / (1.0 + np.exp(-(eta - half)))
    hi = 1.0 / (1.0 + np.exp(-(eta + half)))
    return mean, lo, hi, n


def main():
    fs = pd.read_csv(os.path.join(RES, "cv_fold_scores.csv"))
    rows = []
    for m in ORDER:
        col = "bacc_" + m
        if col not in fs.columns:
            continue
        mean, lo, hi, n = logit_interval(fs[col])
        cov = (fs["coverage"].mean() if m == "glassbox_covered"
               else 1.0 if m.endswith("forced") or m != "glassbox_covered" else np.nan)
        rows.append({"method": m, "n_folds": n, "bacc_mean": mean,
                     "ci_low": lo, "ci_high": hi, "coverage_mean": cov,
                     "interval": "logit-scale Student interval, delta method"})
    out = pd.DataFrame(rows)
    assert ((out.ci_low >= 0) & (out.ci_high <= 1)).all(), "interval outside [0, 1]"
    assert (out.ci_low <= out.bacc_mean).all() and (out.bacc_mean <= out.ci_high).all()
    out.to_csv(os.path.join(RES, "cv_parity_folds.csv"), index=False,
               encoding="utf-8-sig")
    print(out.round(6).to_string(index=False))


if __name__ == "__main__":
    main()

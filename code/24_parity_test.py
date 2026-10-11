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



def non_inferiority(fs, margin=0.01):
    """Is the covered mode non-inferior to the strongest baseline on all specimens?

    One-sided paired test of H0: baseline minus glass box >= margin, against the
    alternative that the shortfall is smaller than the margin. A significant
    result supports a non-inferiority claim at that margin; it is not a test of
    exact equality.
    """
    rows = []
    for b in BASELINES:
        for mode in ("glassbox_covered", "glassbox_forced"):
            d = fs["bacc_" + b].to_numpy() - fs["bacc_" + mode].to_numpy()
            n = len(d)
            se = d.std(ddof=1) / np.sqrt(n)
            t = (d.mean() - margin) / se
            p = stats.t.cdf(t, n - 1)          # one-sided, lower tail
            rows.append({"baseline": b, "glass_box_mode": mode, "margin": margin,
                         "mean_shortfall": round(d.mean(), 4), "t": round(t, 2),
                         "p_one_sided": round(float(p), 5),
                         "non_inferior_at_margin": bool(p < 0.05)})
    return pd.DataFrame(rows)


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
    ni = non_inferiority(fs, margin=0.01)
    ni.to_csv(os.path.join(RES, "non_inferiority_test.csv"), index=False, encoding="utf-8-sig")
    print(ni.to_string(index=False))
    print(out.to_string(index=False))
    print(f"strongest baseline: {out.sort_values('mean_difference').iloc[-1].comparison}")


if __name__ == "__main__":
    main()

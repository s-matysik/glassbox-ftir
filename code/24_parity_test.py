"""The significance test behind the parity claim, computed from the fold record.

The frozen protocol specifies the corrected resampled test of Nadeau and Bengio.
That test is usually described for repeated resampling, but its content is a
variance correction for overlapping training sets, and that correction applies
to a single stratified k-fold run as well. Both forms are therefore computed
and reported: the uncorrected paired t-test and the corrected one. This script computes it from
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


def paired(a, b, corrected=False, n_test=600, n_train=2400):
    """Paired test over the folds, optionally with the Nadeau-Bengio correction.

    The uncorrected standard error treats the k fold differences as independent.
    They are not: the training sets overlap, so the variance of the mean
    difference is underestimated and the interval is too narrow. Nadeau and
    Bengio replace the factor 1/k by 1/k + n_test/n_train, which for this design
    inflates the variance by 2.25 and divides the statistic by 1.5. The
    correction requires no additional runs; it is a statement about the
    dependence between folds that have already been computed.
    """
    d = np.asarray(b, float) - np.asarray(a, float)
    n = len(d)
    if corrected:
        se = np.sqrt((1.0 / n + n_test / n_train) * d.var(ddof=1))
    else:
        se = d.std(ddof=1) / np.sqrt(n)
    t = d.mean() / se
    p = 2 * (1 - stats.t.cdf(abs(t), n - 1))
    half = stats.t.ppf(0.975, n - 1) * se
    return d.mean(), d.mean() - half, d.mean() + half, t, p, n



def non_inferiority(fs, margin=0.01):
    """Is the covered mode non-inferior to the strongest baseline on all specimens?

    The standard error carries the Nadeau-Bengio correction for overlapping
    training sets, so this is the conservative form of the test.

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
            se = np.sqrt((1.0 / n + 600.0 / 2400.0) * d.var(ddof=1))
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
        for corrected, label in [(False, "paired t-test across folds, uncorrected"),
                                 (True, "corrected resampled t-test, Nadeau and Bengio")]:
            diff, lo, hi, t, p, n = paired(fs.bacc_glassbox_forced, fs["bacc_" + b],
                                           corrected=corrected)
            rows.append({"comparison": f"{b} minus glass box, forced choice", "n_folds": n,
                         "mean_difference": round(diff, 4), "ci_low": round(lo, 4),
                         "ci_high": round(hi, 4), "t": round(t, 2), "p": round(p, 5),
                         "test": label})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RES, "parity_test.csv"), index=False, encoding="utf-8-sig")
    ni = non_inferiority(fs, margin=0.01)
    ni.to_csv(os.path.join(RES, "non_inferiority_test.csv"), index=False, encoding="utf-8-sig")
    print(ni.to_string(index=False))
    print(out.to_string(index=False))
    print(out[out.test.str.startswith("corrected")].to_string(index=False))
    print(f"strongest baseline: {out.sort_values('mean_difference').iloc[-1].comparison}")


if __name__ == "__main__":
    main()

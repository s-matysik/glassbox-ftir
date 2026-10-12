"""Stability of the selected attribute set, compared at matched set size.

The article first compared a two-attribute rule against the top three of a
post-hoc attribution, which are sets of different size. Size matters for every
set-stability measure, so the comparison has to be made at a common size. This
script recomputes both at size two and size three from the per-repeat ranks
already written by code/13_shap_oc.py, so it refits nothing.

Two measures are reported. The modal share is the fraction of repeats returning
the most common set, which is what the earlier text quoted. The mean pairwise
Jaccard index over all pairs of repeats is less brittle, because one swapped
member does not look the same as a complete change of set.

Input : results/shap_ranks_per_repeat.csv
Output: results/stability_matched.csv

Run from the deposit root:  python code/27_stability_matched.py
"""
import os
from collections import Counter
from itertools import combinations

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")

ATTR = ["A_W1", "C_W1", "O_W1", "A_W2", "C_W2", "O_W2", "A_W3", "C_W3", "O_W3"]


def summarise(sets, label, k):
    keys = [", ".join(sorted(s)) for s in sets]
    counts = Counter(keys)
    modal, n_modal = counts.most_common(1)[0]
    pairs = [len(a & b) / len(a | b) for a, b in combinations(sets, 2)]
    return {"selector": label, "set_size": k, "repeats": len(sets),
            "modal_set": modal, "modal_share": round(n_modal / len(sets), 4),
            "n_distinct_sets": len(counts),
            "jaccard_mean": round(float(np.mean(pairs)), 4),
            "jaccard_min": round(float(np.min(pairs)), 4)}


def main():
    rk = pd.read_csv(os.path.join(RES, "shap_ranks_per_repeat.csv"))
    missing = [c for c in ATTR if c not in rk.columns]
    assert not missing, f"rank file is missing attributes: {missing}"

    rows = []
    for k in (2, 3):
        # rank 1 is the most important attribute, so the top k are the k smallest
        shap_sets = [set(rk.loc[i, ATTR].astype(float).nsmallest(k).index)
                     for i in rk.index]
        rows.append(summarise(shap_sets, "post-hoc attribution, top k", k))

    gb = [set(s.split(", ")) for s in rk.glassbox_set]
    sizes = {len(s) for s in gb}
    assert sizes == {2}, f"glass box sets are not all of size two: {sizes}"
    rows.append(summarise(gb, "glass box, selected attributes", 2))

    gc = [set(s.split(", ")) for s in rk.gamma_criterion_set]
    rows.append(summarise(gc, "glass box, alternative criterion", 2))

    out = pd.DataFrame(rows).sort_values(["set_size", "selector"])
    out.to_csv(os.path.join(RES, "stability_matched.csv"), index=False,
               encoding="utf-8-sig")
    print(out.to_string(index=False))

    at2 = out[(out.set_size == 2)].set_index("selector")
    print(f"\nat size two the attribution is "
          f"{'more' if at2.loc['post-hoc attribution, top k', 'jaccard_mean'] > at2.loc['glass box, selected attributes', 'jaccard_mean'] else 'less'} "
          f"stable than the rule's own selection")


if __name__ == "__main__":
    main()

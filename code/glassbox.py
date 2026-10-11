"""The glass-box classifier of PROTOCOL_v1.1, assembled from ``ara_aco``.

Every step is a call into ``ara_aco``; nothing is reimplemented here. The
module only sequences the steps and holds the fitted state:

1. ACO attributes in the three protocol windows against a template that is the
   global median of the training rows (``make_template``, ``aco_matrix``).
2. The attributes are carried on the change scale of protocol Section 5,
   ``delta()``: ``ln A`` for amplitude, the raw value for correlation and
   offset. ``sens.min`` is calibrated on that same scale, so the rough-set
   bins -- whose width is ``sens.min`` -- are built on the delta-scaled matrix
   rather than on the raw attribute matrix. Mixing the two scales would give
   amplitude bins whose width is a log-scale threshold applied to a linear
   quantity.
3. ``sens.min`` from the Savitzky-Golay noise estimate propagated by Monte
   Carlo, worst case over polymers, on the training rows only.
4. At most two condition attributes, chosen by an inner stratified
   cross-validation of the training rows. Template and ``sens.min`` are
   re-estimated inside every inner fold, so the selection never sees a record
   it is scored on.
5. Rough-set model and the three-way decision at test time.

``preprocess`` (three-point median filter) and ``instrument_noise`` act on one
spectrum at a time and use no information from other records, so they are
applied once to the whole matrix by the calling script.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np
from sklearn.model_selection import StratifiedKFold

import ara_aco as aa


def candidate_sets(names: list[str], max_attrs: int = 2) -> list[list[str]]:
    out = [[n] for n in names]
    for r in range(2, max_attrs + 1):
        out += [list(c) for c in itertools.combinations(names, r)]
    return out


def balanced_accuracy(y_true: np.ndarray, y_pred: np.ndarray, classes) -> float:
    vals = []
    for c in classes:
        m = y_true == c
        if m.sum():
            vals.append(float((y_pred[m] == c).mean()))
    return float(np.mean(vals)) if vals else np.nan


@dataclass
class GlassBox:
    masks: dict
    names: list
    template: dict
    sens: aa.SensMin
    model: aa.RoughModel
    attributes: list
    inner_score: float
    classes: list
    candidate_scores: list = None  # (attribute set, inner balanced accuracy)

    def transform(self, Yp: np.ndarray) -> np.ndarray:
        F, _ = aa.aco_matrix(Yp, self.masks, self.template)
        return aa.delta(F, self.names)

    def predict(self, Yp: np.ndarray):
        """Return ``(decisions, forced)``: the three-way decision, in which
        records outside the lower approximation are referred, and the
        forced-choice label for the same records."""
        return aa.three_way(self.model, self.transform(Yp), self.names)


def _fit_core(Yp: np.ndarray, sigma: np.ndarray, y: np.ndarray, masks, names,
              k: float, n_mc: int, seed: int):
    template = aa.make_template(Yp, masks)
    F, _ = aa.aco_matrix(Yp, masks, template)
    sens = aa.calibrate_sens_min(Yp, sigma, y, masks, template, names,
                                 k=k, n_mc=n_mc, seed=seed)
    return template, aa.delta(F, names), sens


def fit(Yp: np.ndarray, sigma: np.ndarray, y: np.ndarray, masks, names,
        k: float = 4.0, n_mc: int = 100, seed: int = 0, n_inner: int = 3,
        n_bins: int | None = None, max_attrs: int = 2,
        candidates: list[list[str]] | None = None) -> GlassBox:
    """Fit on training rows only. ``seed`` fixes both the Monte Carlo draws of
    the noise propagation and the inner split used for attribute selection."""
    y = np.asarray(y).astype(str)
    classes = sorted(set(y))
    cands = candidates if candidates is not None else candidate_sets(names, max_attrs)

    counts = np.array([(y == c).sum() for c in classes])
    n_split = int(min(n_inner, counts.min()))
    scores = np.zeros(len(cands))
    coverage = np.zeros(len(cands))
    if n_split >= 2:
        inner = StratifiedKFold(n_splits=n_split, shuffle=True, random_state=seed)
        for itr, ite in inner.split(np.zeros(len(y)), y):
            t_i, D_i, s_i = _fit_core(Yp[itr], sigma[itr], y[itr], masks, names,
                                      k=k, n_mc=n_mc, seed=seed)
            F_h, _ = aa.aco_matrix(Yp[ite], masks, t_i)
            D_h = aa.delta(F_h, names)
            for j, att in enumerate(cands):
                m = aa.fit_rough(D_i, names, y[itr], att, s_i, n_bins=n_bins)
                dec, forced = aa.three_way(m, D_h, names)
                scores[j] += balanced_accuracy(y[ite], forced, classes)
                coverage[j] += float((dec != "refer").mean())
        scores /= n_split
        coverage /= n_split

    template, D, sens = _fit_core(Yp, sigma, y, masks, names, k=k, n_mc=n_mc, seed=seed)
    order = sorted(range(len(cands)), key=lambda j: (-scores[j], -coverage[j],
                                                     len(cands[j]), cands[j]))
    best = cands[order[0]]
    model = aa.fit_rough(D, names, y, best, sens, n_bins=n_bins)
    return GlassBox(masks=masks, names=names, template=template, sens=sens,
                    model=model, attributes=best, inner_score=float(scores[order[0]]),
                    classes=classes,
                    candidate_scores=[("+".join(c), float(s_)) for c, s_ in zip(cands, scores)])

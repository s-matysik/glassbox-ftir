"""Verification of the pipeline on synthetic spectra with a known structure.

Two parameters are planted so that each drives exactly one window: P1 scales
the amplitude of window W1, P2 shifts the band shape of window W2. A correct
implementation must recover that mapping and must not attribute either
parameter to the other's window.

Run from the deposit root:  python code/test_ara_aco.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ara_aco import (aco_matrix, calibrate_sens_min, exclusive_attributes,
                     fit_rough, instrument_noise, make_template,
                     overlap_coefficient, pair_rules, preprocess, three_way,
                     window_masks)

GRID = np.linspace(400, 4000, 900)
WINDOWS = {"W1": (2800, 3000), "W2": (1650, 1800), "W3": (1350, 1500)}


def gauss(centre, width, height):
    return height * np.exp(-0.5 * ((GRID - centre) / width) ** 2)


def synth(n_per_cell=40, noise=0.02, seed=0):
    rng = np.random.default_rng(seed)
    rows, p1, p2 = [], [], []
    for a in (0, 1, 2):                      # P1: amplitude of W1
        for b in (0, 1):                     # P2: band position inside W2
            for _ in range(n_per_cell):
                s = (gauss(2900, 40, 0.5 + 0.9 * a)
                     + gauss(1700 + 45 * b, 25, 0.8)
                     + gauss(1420, 30, 0.6)
                     + 0.25)
                rows.append(s + rng.normal(0, noise, len(GRID)))
                p1.append(a)
                p2.append(b)
    return np.array(rows), np.array(p1), np.array(p2)


def main():
    Y, p1, p2 = synth()
    Y = preprocess(Y)
    masks = window_masks(GRID, WINDOWS)
    template = make_template(Y, masks)
    F, names = aco_matrix(Y, masks, template)
    sigma = instrument_noise(Y)
    groups = np.arange(len(Y))

    sens = calibrate_sens_min(Y, sigma, p1, masks, template, names, k=4.0, n_mc=40)
    assert all(v > 0 for v in sens.thresholds.values()), "non-positive threshold"
    assert 0.9 < sens.confidence < 1.0, f"implied confidence off: {sens.confidence}"

    rules = {"P1": pair_rules(F, names, p1, groups, sens, max_pairs_per_stratum=400,
                              seed=1, holding=p2)["rules"],
             "P2": pair_rules(F, names, p2, groups, sens, max_pairs_per_stratum=400,
                              seed=2, holding=p1)["rules"]}
    excl = exclusive_attributes(rules, theta_s=0.80, theta_i=0.70)

    top = {p: [r["attribute"] for r in sorted(rs, key=lambda r: -r["conf_change"])][:3]
           for p, rs in rules.items()}
    assert any(a.endswith("W1") for a in top["P1"]), f"P1 not mapped to W1: {top['P1']}"
    assert any(a.endswith("W2") for a in top["P2"]), f"P2 not mapped to W2: {top['P2']}"

    a_w1 = names.index("A_W1")
    assert overlap_coefficient(F[p1 == 0, a_w1], F[p1 == 2, a_w1]) < 0.2, "W1 should separate P1"
    c_w1 = names.index("C_W1")
    assert overlap_coefficient(F[p2 == 0, c_w1], F[p2 == 1, c_w1]) > 0.5, "W1 should ignore P2"

    model = fit_rough(F, names, p1.astype(str), ["A_W1"], sens)
    assert model.gamma > 0.9, f"rough-set quality too low: {model.gamma}"
    dec, forced = three_way(model, F, names)
    coverage = float((dec != "refer").mean())
    acc = float((forced == p1.astype(str)).mean())
    assert acc > 0.9, f"forced-choice accuracy too low: {acc}"

    print(f"sens.min: {({k: round(v, 4) for k, v in sens.thresholds.items()})}")
    print(f"implied confidence of k=4: {sens.confidence:.4f}")
    print(f"top attributes P1 {top['P1']} | P2 {top['P2']}")
    print(f"exclusive: {excl}")
    print(f"rough sets on A_W1: gamma {model.gamma:.3f} | coverage {coverage:.3f} | "
          f"forced accuracy {acc:.3f}")
    print("all assertions passed")


if __name__ == "__main__":
    main()

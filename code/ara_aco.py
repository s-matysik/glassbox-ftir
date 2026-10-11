"""ARA-ACO-rough-set pipeline for FTIR spectra.

Implements Steps 0 to 9 of the article concept (xAI_GOZ rev. 3, Section 6.8)
under the frozen protocol PROTOCOL_v1.1.md. Equation labels (A1) to (A17)
follow the concept.

Nothing in this module estimates anything from test data: every calibration
function takes a training index and returns an object that is then applied
unchanged to held-out records.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from itertools import combinations

import numpy as np
from scipy import stats
from scipy.signal import medfilt, savgol_filter

# --------------------------------------------------------------------------
# Step 0: preprocessing
# --------------------------------------------------------------------------

def preprocess(Y: np.ndarray, median_width: int = 3) -> np.ndarray:
    """Three-point median filter against isolated spikes. No baseline removal
    and no smoothing beyond this: the offset attribute carries the baseline."""
    if median_width <= 1:
        return np.asarray(Y, float)
    return np.apply_along_axis(medfilt, 1, np.asarray(Y, float), median_width)


def instrument_noise(Y: np.ndarray, window: int = 11, order: int = 3) -> np.ndarray:
    """Per-spectrum noise: the high-frequency residual left by a
    Savitzky-Golay filter. Substitutes for replicate scans on datasets that
    have none (protocol Section 4)."""
    Y = np.asarray(Y, float)
    return (Y - savgol_filter(Y, window_length=window, polyorder=order, axis=1)).std(axis=1, ddof=1)


# --------------------------------------------------------------------------
# Step 1: ACO decomposition
# --------------------------------------------------------------------------

def window_masks(grid: np.ndarray, windows: dict[str, tuple[float, float]]) -> dict[str, np.ndarray]:
    return {name: (grid >= lo) & (grid <= hi) for name, (lo, hi) in windows.items()}


def aco_matrix(Y: np.ndarray, masks: dict[str, np.ndarray],
               template: dict[str, np.ndarray]) -> tuple[np.ndarray, list[str]]:
    """Nine attributes: amplitude, correlation and offset per window (A1)-(A3)."""
    Y = np.asarray(Y, float)
    cols, names = [], []
    for name, m in masks.items():
        S = Y[:, m]
        A = S.max(axis=1) - S.min(axis=1)
        O = S.min(axis=1)
        rng = np.where(A > 0, A, 1.0)
        Sn = (S - O[:, None]) / rng[:, None]
        T = template[name]
        Tc = T - T.mean()
        Sc = Sn - Sn.mean(axis=1, keepdims=True)
        denom = np.sqrt((Sc ** 2).sum(axis=1) * (Tc ** 2).sum()) + 1e-12
        C = (Sc @ Tc) / denom
        cols += [A, C, O]
        names += [f"A_{name}", f"C_{name}", f"O_{name}"]
    return np.column_stack(cols), names


def make_template(Y: np.ndarray, masks: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Global template: min-max normalised median spectrum of the training rows."""
    med = np.median(np.asarray(Y, float), axis=0)
    out = {}
    for name, m in masks.items():
        seg = med[m]
        rng = seg.max() - seg.min()
        out[name] = (seg - seg.min()) / (rng if rng > 0 else 1.0)
    return out


def delta(values: np.ndarray, names: list[str]) -> np.ndarray:
    """Change scale per attribute (A4): log for amplitude, absolute otherwise."""
    out = values.copy()
    for j, n in enumerate(names):
        if n.startswith("A_"):
            out[:, j] = np.log(np.clip(values[:, j], 1e-9, None))
    return out


# --------------------------------------------------------------------------
# Step 2: sens.min calibration
# --------------------------------------------------------------------------

def k_from_t(n: int, alpha: float = 0.05) -> float:
    """Concept Eq. (A7)."""
    df = max(n - 1, 1)
    return float(stats.t.ppf(1 - alpha / 2, df) * np.sqrt(1 + 1 / n) * np.sqrt(2))


def implied_confidence(k: float, n: int) -> float:
    """Concept Eq. (A8): the two-sided confidence that k actually delivers."""
    df = max(n - 1, 1)
    return float(2 * stats.t.cdf(k / (np.sqrt(1 + 1 / n) * np.sqrt(2)), df) - 1)


@dataclass
class SensMin:
    thresholds: dict[str, float]
    sd_max: dict[str, float]
    k: float
    n_effective: int
    confidence: float
    source: str


def calibrate_sens_min(Y: np.ndarray, sigma: np.ndarray, condition: np.ndarray,
                       masks, template, names, k: float = 4.0, n_mc: int = 100,
                       seed: int = 11) -> SensMin:
    """Propagate per-spectrum instrument noise into every attribute by Monte
    Carlo, take the worst-case condition (A5)-(A6), then sens.min = k * SD_max."""
    rng = np.random.default_rng(seed)
    Y = np.asarray(Y, float)
    sd_per_row = np.empty((len(Y), len(names)))
    for i in range(len(Y)):
        draws = Y[i][None, :] + rng.normal(0, sigma[i], size=(n_mc, Y.shape[1]))
        F, _ = aco_matrix(draws, masks, template)
        sd_per_row[i] = delta(F, names).std(axis=0, ddof=1)
    sd_max = {}
    for j, nm in enumerate(names):
        per_cond = [np.median(sd_per_row[condition == c, j]) for c in np.unique(condition)]
        sd_max[nm] = float(np.max(per_cond))
    return SensMin(thresholds={nm: k * v for nm, v in sd_max.items()},
                   sd_max=sd_max, k=k, n_effective=n_mc,
                   confidence=implied_confidence(k, n_mc),
                   source="savitzky-golay residual, Monte Carlo propagation")


# --------------------------------------------------------------------------
# Steps 3 and 4: pair database and weighted association rules
# --------------------------------------------------------------------------

def ternary_labels(d: np.ndarray, threshold: float) -> np.ndarray:
    """Concept Eq. (A11)."""
    return np.where(d > threshold, 1, np.where(d < -threshold, -1, 0))


def wilson(p: float, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return (max(0.0, c - h), min(1.0, c + h))


def pair_rules(F: np.ndarray, names: list[str], parameter: np.ndarray,
               groups: np.ndarray, sens: SensMin, weights: np.ndarray | None = None,
               max_pairs_per_stratum: int = 2000, seed: int = 7,
               n_perm: int = 0, holding: np.ndarray | None = None) -> dict:
    """Rules of the form: parameter changes -> attribute moves beyond sens.min.

    Concept Eq. (A10): a pair enters the database only if the two records
    differ in this parameter, come from different groups, and agree on every
    other parameter. `holding` carries those other parameters as one label per
    record; pairs are formed inside a held-fixed level, never across levels.
    Each parameter-value pair is capped at `max_pairs_per_stratum` draws so the
    database stays tractable and no stratum dominates (protocol Section 7).
    """
    rng = np.random.default_rng(seed)
    D = delta(F, names)
    w = np.ones(len(F)) if weights is None else np.asarray(weights, float)
    hold = np.zeros(len(F), dtype=int) if holding is None else np.asarray(holding)
    vals = np.unique(parameter)

    counts = {nm: {1: 0.0, 0: 0.0, -1: 0.0} for nm in names}
    total = 0.0
    pairs_used = 0
    for h in np.unique(hold):
        hm = hold == h
        idx_by_val = {v: np.where(hm & (parameter == v))[0] for v in vals}
        for a, b in combinations(vals, 2):
            ia, ib = idx_by_val[a], idx_by_val[b]
            if len(ia) == 0 or len(ib) == 0:
                continue
            n_draw = min(max_pairs_per_stratum, len(ia) * len(ib))
            pa = rng.choice(ia, n_draw, replace=len(ia) < n_draw)
            pb = rng.choice(ib, n_draw, replace=len(ib) < n_draw)
            keep = groups[pa] != groups[pb]
            pa, pb = pa[keep], pb[keep]
            if len(pa) == 0:
                continue
            pw = w[pa] * w[pb]
            pairs_used += len(pa)
            for j, nm in enumerate(names):
                lab = ternary_labels(D[pb, j] - D[pa, j], sens.thresholds[nm])
                for level in (1, 0, -1):
                    counts[nm][level] += float(pw[lab == level].sum())
            total += float(pw.sum())

    rules = []
    for nm in names:
        c = counts[nm]
        changed = c[1] + c[-1]
        conf_change = changed / total if total else np.nan
        direction = 1 if c[1] >= c[-1] else -1
        conf_dir = max(c[1], c[-1]) / total if total else np.nan
        lo, hi = wilson(conf_change, pairs_used)
        rules.append({"attribute": nm, "conf_change": conf_change,
                      "conf_invariant": c[0] / total if total else np.nan,
                      "direction": direction, "conf_direction": conf_dir,
                      "support": changed / total if total else np.nan,
                      "ci_low": lo, "ci_high": hi, "n_pairs": pairs_used})

    if n_perm > 0:
        obs = np.array([r["conf_change"] for r in rules])
        null = np.empty((n_perm, len(names)))
        for t in range(n_perm):
            perm = rng.permutation(parameter)
            sub = pair_rules(F, names, perm, groups, sens, weights,
                             max_pairs_per_stratum=max(200, max_pairs_per_stratum // 10),
                             seed=int(rng.integers(1e9)), n_perm=0)
            null[t] = [r["conf_change"] for r in sub["rules"]]
        p = ((null >= obs[None, :]).sum(axis=0) + 1) / (n_perm + 1)
        order = np.argsort(p)
        m = len(p)
        q = np.empty(m)
        prev = 1.0
        for rank, j in enumerate(order[::-1]):
            prev = min(prev, p[j] * m / (m - rank))
            q[j] = prev
        for j, r in enumerate(rules):
            r["p_perm"] = float(p[j])
            r["q_fdr"] = float(q[j])
    return {"rules": rules, "n_pairs": pairs_used, "total_weight": total}


def exclusive_attributes(rules_by_parameter: dict[str, list[dict]],
                         theta_s: float = 0.90, theta_i: float = 0.80) -> dict[str, list[str]]:
    """Concept Eq. (A14): sensitive to its own parameter, invariant to the rest."""
    out = {}
    for p, rules in rules_by_parameter.items():
        own = {r["attribute"]: r for r in rules}
        keep = []
        for nm, r in own.items():
            if r["conf_change"] < theta_s:
                continue
            if all(rules_by_parameter[q][[x["attribute"] for x in rules_by_parameter[q]].index(nm)]["conf_invariant"] >= theta_i
                   for q in rules_by_parameter if q != p):
                keep.append(nm)
        out[p] = sorted(keep, key=lambda nm: -own[nm]["conf_change"])
    return out


# --------------------------------------------------------------------------
# Step 6: rough sets
# --------------------------------------------------------------------------

@dataclass
class RoughModel:
    attributes: list[str]
    edges: dict[str, np.ndarray]
    cell_class: dict[tuple, str]
    cell_pure: dict[tuple, bool]
    alpha: dict[str, float]
    gamma: float
    boundary: dict[str, int]
    classes: list[str] = field(default_factory=list)

    def cells(self, F: np.ndarray, names: list[str]) -> list[tuple]:
        idx = [names.index(a) for a in self.attributes]
        out = []
        for row in F[:, idx]:
            out.append(tuple(int(np.digitize(v, self.edges[a])) for v, a in zip(row, self.attributes)))
        return out


def fit_rough(F: np.ndarray, names: list[str], y: np.ndarray, attributes: list[str],
              sens: SensMin, n_bins: int | None = None) -> RoughModel:
    """Bins of width sens.min anchored at the training minimum (protocol
    Section 7); `n_bins` switches to equal-width bins for the sensitivity run."""
    idx = [names.index(a) for a in attributes]
    edges = {}
    for a, j in zip(attributes, idx):
        col = F[:, j]
        if n_bins:
            edges[a] = np.linspace(col.min(), col.max(), n_bins + 1)[1:-1]
        else:
            width = sens.thresholds[a]
            n = max(int(np.ceil((col.max() - col.min()) / width)), 1)
            edges[a] = col.min() + width * np.arange(1, n)
    model = RoughModel(attributes=attributes, edges=edges, cell_class={},
                       cell_pure={}, alpha={}, gamma=0.0, boundary={},
                       classes=sorted(set(map(str, y))))
    cells = model.cells(F, names)
    by_cell: dict[tuple, list[str]] = {}
    for c, lab in zip(cells, map(str, y)):
        by_cell.setdefault(c, []).append(lab)
    certain = 0
    for c, labs in by_cell.items():
        uniq = set(labs)
        model.cell_pure[c] = len(uniq) == 1
        vals, cnt = np.unique(labs, return_counts=True)
        model.cell_class[c] = str(vals[np.argmax(cnt)])
        if len(uniq) == 1:
            certain += len(labs)
    for cls in model.classes:
        lower = sum(len(l) for c, l in by_cell.items() if model.cell_pure[c] and model.cell_class[c] == cls)
        upper = sum(len(l) for c, l in by_cell.items() if cls in set(l))
        model.alpha[cls] = lower / upper if upper else np.nan
        model.boundary[cls] = upper - lower
    model.gamma = certain / len(y)
    return model


def three_way(model: RoughModel, F: np.ndarray, names: list[str]):
    """Step 7: decide inside the lower approximation, refer otherwise."""
    decisions, forced = [], []
    for c in model.cells(F, names):
        if c in model.cell_class and model.cell_pure[c]:
            decisions.append(model.cell_class[c])
            forced.append(model.cell_class[c])
        elif c in model.cell_class:
            decisions.append("refer")
            forced.append(model.cell_class[c])
        else:
            decisions.append("refer")
            forced.append(_nearest_cell_class(model, c))
    return np.array(decisions), np.array(forced)


def _nearest_cell_class(model: RoughModel, cell: tuple) -> str:
    best, best_d = None, np.inf
    for c in model.cell_class:
        d = sum(abs(a - b) for a, b in zip(c, cell))
        if d < best_d:
            best, best_d = c, d
    return model.cell_class[best] if best is not None else model.classes[0]


def overlap_coefficient(x: np.ndarray, y: np.ndarray, n_grid: int = 512) -> float:
    """Concept Eq. (A17): Gaussian kernel density, Silverman bandwidth."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3 or len(y) < 3:
        return np.nan
    lo = min(x.min(), y.min())
    hi = max(x.max(), y.max())
    if hi <= lo:
        return 1.0
    grid = np.linspace(lo, hi, n_grid)
    fx = stats.gaussian_kde(x, bw_method="silverman")(grid)
    fy = stats.gaussian_kde(y, bw_method="silverman")(grid)
    return float(np.trapezoid(np.minimum(fx, fy), grid))

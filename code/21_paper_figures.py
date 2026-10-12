"""Figures for the manuscript, rebuilt from the tables in results/.

Every figure reads a CSV; none recomputes a result, so a figure can never
disagree with the table it illustrates. Column width 3.5 in, page width 7.16 in,
300 dpi, Times New Roman where available, per IEEE Access graphics guidance.

Run from the deposit root:  python code/21_paper_figures.py
"""
import os
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")
COL, PAGE = 3.5, 7.16
GREY = "#8c8c8c"
ACC = "#1f4e79"
WARN = "#c0392b"
POLY_ORDER = ["HDPE", "LDPE", "PET", "PP", "PS", "PVC"]
NICE = {"A": "amplitude", "C": "shape correlation", "O": "baseline offset"}
WIN = {"W1": "C-H stretch", "W2": "carbonyl", "W3": "CH deformation"}


def setup():
    os.makedirs(FIG, exist_ok=True)
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "DejaVu Serif"],
        "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight"})


def label(a):
    kind, win = a.split("_")
    return f"{NICE[kind]}, {WIN[win]}"


def fig_entry():
    d = pd.read_csv(os.path.join(RES, "entry_condition_primary.csv"))
    d = d.set_index("stratum").loc[POLY_ORDER].reset_index()
    fig, ax = plt.subplots(figsize=(COL, 2.1))
    ax.bar(d["stratum"], d["spread_over_noise"], color=ACC, width=0.62)
    ax.axhline(1, color=WARN, lw=1, ls="--")
    ax.text(5.4, 2.0, "method inapplicable below this line", ha="right",
            fontsize=6, color=WARN)
    ax.set_ylabel("between-specimen spread /\ninstrument noise")
    ax.set_yscale("log")
    for i, v in enumerate(d["spread_over_noise"]):
        ax.text(i, v * 1.08, f"{v:.0f}", ha="center", fontsize=6.5)
    lo, hi = d["spread_over_noise"].min(), d["spread_over_noise"].max()
    ax.set_title(f"Entry condition is met by {lo:.0f} to {hi:.0f} times")
    fig.savefig(os.path.join(FIG, "matys2.png"))
    plt.close(fig)
    return lo, hi


def fig_rules():
    d = pd.read_csv(os.path.join(RES, "rule_matrix.csv"))
    piv = d.pivot_table(index="attribute", columns="body", values="conf_change")
    order = piv.mean(axis=1).sort_values(ascending=False).index
    piv = piv.loc[order]
    fig, ax = plt.subplots(figsize=(PAGE, 2.6))
    im = ax.imshow(piv.values, aspect="auto", cmap="viridis", vmin=0.2, vmax=1.0)
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels([label(a) for a in piv.index], fontsize=6.5)
    ax.set_xticks(range(len(piv.columns)))
    ax.set_xticklabels([c.replace(" -> ", "/").replace("any polymer change", "pooled")
                        for c in piv.columns], rotation=90, fontsize=5.8)
    ax.set_title("Confidence that an attribute moves beyond its detection "
                 "threshold when the polymer changes")
    cb = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.012)
    cb.set_label("confidence", fontsize=7)
    fig.savefig(os.path.join(FIG, "matys4.png"))
    plt.close(fig)
    return float(np.nanmedian(piv.values))


def fig_ovl():
    d = pd.read_csv(os.path.join(RES, "ovl_matrix.csv"))
    # The file also holds six one-against-the-rest rows. Every summary in the
    # article is taken over the 15 genuine polymer pairs only, so that the
    # denominator quoted in the text is the one the figure shows.
    d = d[~d.comparison.str.contains("rest")]
    med = d.groupby("attribute").ovl.median().sort_values()
    fig, ax = plt.subplots(figsize=(COL, 2.3))
    ax.barh([label(a) for a in med.index], med.values,
            color=[ACC if v < 0.3 else GREY for v in med.values])
    ax.set_xlabel("median overlap between polymer pairs")
    ax.invert_yaxis()
    ax.set_title("Separation achieved by each named attribute")
    ax.tick_params(axis="y", labelsize=6.5)
    fig.savefig(os.path.join(FIG, "matys5.png"))
    plt.close(fig)
    return med.index[0], float(med.iloc[0])


def fig_oc():
    d = pd.read_csv(os.path.join(RES, "oc_curves.csv"))
    plans = sorted(d.groupby(["n", "c"]).groups)
    fig, axes = plt.subplots(1, 2, figsize=(PAGE, 2.5), sharey=True)
    for ax, inspector in zip(axes, ["perfect inspector", "observed inspector"]):
        for (n, c), colour in zip(plans, [ACC, "#4a90c2", "#e0a458", WARN]):
            g = d[(d.n == n) & (d.c == c) & (d.inspector == inspector)]
            ax.plot(g.lot_fraction_defective * 100, g.prob_accept, color=colour,
                    lw=1.3, label=f"n = {n}, c = {c}")
        ax.set_xlabel("lot fraction defective (%)")
        ax.set_title(inspector)
        ax.axvline(1, color=GREY, lw=0.7, ls=":")
        ax.axvline(10, color=GREY, lw=0.7, ls=":")
    axes[0].set_ylabel("probability of acceptance")
    axes[0].legend(frameon=False, loc="upper right")
    for ax, letter in zip(axes, "ab"):
        ax.annotate(letter, xy=(-0.08, 1.08), xycoords="axes fraction",
                    fontsize=9, fontweight="bold", va="top")
    fig.savefig(os.path.join(FIG, "matys8.png"))
    plt.close(fig)
    return d


def fig_parity():
    d = pd.read_csv(os.path.join(RES, "cv_parity_folds.csv"))
    nice = {"glassbox_forced": "glass box\nforced choice",
            "glassbox_covered": "glass box\non accepted",
            "svm_rbf": "SVM", "random_forest": "random\nforest",
            "plsda": "PLS-DA", "cnn1d": "1-D CNN"}
    d["lab"] = d.method.map(nice)
    d = d.rename(columns={"bacc_mean": "balanced_accuracy",
                          "coverage_mean": "coverage"})
    d = d.sort_values("balanced_accuracy")
    colours = [WARN if "glass" in m else GREY for m in d.method]
    fig, ax = plt.subplots(figsize=(COL, 2.4))
    err = np.vstack([(d.balanced_accuracy - d.ci_low).clip(lower=0),
                     (d.ci_high - d.balanced_accuracy).clip(lower=0)])
    ax.barh(d.lab, d.balanced_accuracy, xerr=err, color=colours,
            error_kw={"lw": 0.8, "ecolor": "#333333"})
    ax.set_xlim(0.90, 1.005)
    ax.set_xlabel("balanced accuracy")
    for i, (v, cov) in enumerate(zip(d.balanced_accuracy, d.coverage)):
        tag = f"{v:.4f}" + (f"  ({cov:.0%} decided)" if cov < 0.999 else "")
        ax.text(0.903, i, tag, va="center", fontsize=6.3, color="white"
                if v > 0.93 else "black")
    ax.set_title("Parity holds on accepted specimens")
    ax.tick_params(axis="y", labelsize=6.5)
    fig.savefig(os.path.join(FIG, "matys3.png"))
    plt.close(fig)
    return d


def fig_learning():
    d = pd.read_csv(os.path.join(RES, "learning_curves.csv"))
    nice = {"glassbox_forced": "glass box, forced choice", "svm_rbf": "SVM",
            "random_forest": "random forest", "plsda": "PLS-DA", "cnn1d": "1-D CNN"}
    fig, axes = plt.subplots(1, 2, figsize=(PAGE, 2.6))
    style = {"glassbox_forced": (WARN, "-", "o"), "svm_rbf": ("#1f4e79", "--", "s"),
             "plsda": ("#4a90c2", "-.", "^"), "random_forest": ("#6b6b6b", ":", "D"),
             "cnn1d": ("#2e7d32", (0, (3, 1, 1, 1)), "v")}
    for m, g in d.groupby("method"):
        g = g.sort_values("size_pct")
        colour, ls, mk = style.get(m, (GREY, "-", "o"))
        axes[0].plot(g.size_pct, g.balanced_accuracy, marker=mk, ms=3, lw=1.3,
                     color=colour, ls=ls, label=nice.get(m, m))
        axes[0].fill_between(g.size_pct, g.bacc_ci_low, g.bacc_ci_high,
                             color=colour, alpha=0.12, lw=0)
    axes[0].set_xscale("log")
    axes[0].set_xlabel("training specimens (% of full set)")
    axes[0].set_ylabel("balanced accuracy")
    axes[0].set_title("Learning curves")
    axes[0].legend(frameon=False, loc="lower right", fontsize=6)

    gl = d[d.method == "glassbox_forced"].sort_values("size_pct")
    axes[1].plot(gl.size_pct, gl.diff_vs_best_blackbox, marker="o", ms=3, lw=1.3,
                 color=WARN)
    axes[1].fill_between(gl.size_pct, gl.diff_ci_low, gl.diff_ci_high,
                         color=WARN, alpha=0.15, lw=0)
    axes[1].axhline(0, color="#333333", lw=0.8)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("training specimens (% of full set)")
    axes[1].set_ylabel("glass box minus best black box")
    axes[1].set_title("No crossover: the gap narrows but never closes")
    for i, (x, y, n) in enumerate(zip(gl.size_pct, gl.diff_vs_best_blackbox,
                                      gl.n_per_class)):
        if i % 2 == 0:
            axes[1].annotate(f"{int(n)}/class", xy=(x, y), xytext=(0, -11),
                             textcoords="offset points", ha="center", fontsize=5.6,
                             color=GREY)
    for ax, letter in zip(axes, "ab"):
        ax.annotate(letter, xy=(-0.07, 1.08), xycoords="axes fraction",
                    fontsize=9, fontweight="bold", va="top")
    fig.savefig(os.path.join(FIG, "matys6.png"))
    plt.close(fig)
    return gl[["size_pct", "n_per_class", "diff_vs_best_blackbox",
               "diff_ci_high"]].round(4)


def fig_shap():
    path = os.path.join(RES, "shap_ranks_per_repeat.csv")
    if not os.path.exists(path):
        return None
    d = pd.read_csv(path)
    attrs = [c for c in d.columns
             if c not in ("repeat", "glassbox_set", "gamma_criterion_set")]
    R = d[attrs].values
    fig, axes = plt.subplots(1, 2, figsize=(PAGE, 2.6))
    for j, a in enumerate(attrs):
        axes[0].plot(d.repeat, R[:, j], marker="o", ms=2.2, lw=0.8, label=label(a))
    axes[0].set_ylabel("rank (1 = most important)")
    axes[0].set_xlabel("repeated fit")
    axes[0].invert_yaxis()
    axes[0].set_title("Post-hoc attribution ranking moves between fits")
    axes[0].legend(frameon=False, fontsize=5.2, ncol=2, loc="lower right")
    sets = d.glassbox_set.values
    uniq = sorted(set(sets))
    axes[1].scatter(d.repeat, [uniq.index(s) for s in sets], color=WARN, s=14)
    axes[1].set_yticks(range(len(uniq)))
    axes[1].set_yticklabels(uniq, fontsize=6.5)
    axes[1].set_xlabel("repeated fit")
    axes[1].set_title("Glass-box attribute set does not")
    axes[1].set_ylim(-0.6, max(len(uniq) - 0.4, 0.6))
    for ax, letter in zip(axes, "ab"):
        ax.annotate(letter, xy=(-0.06, 1.08), xycoords="axes fraction",
                    fontsize=9, fontweight="bold", va="top")
    fig.savefig(os.path.join(FIG, "matys7.png"))
    plt.close(fig)
    return len(uniq)



def fig_windows():
    """Mean spectrum per polymer with the three analysis windows shaded."""
    import plastics_data as pdata
    grid, Y, meta = pdata.load("c8")
    y = meta["Polymer"].to_numpy().astype(str)
    from ara_aco import window_masks
    masks = window_masks(grid, pdata.WINDOWS)
    fig, ax = plt.subplots(figsize=(PAGE, 2.9))
    cols = {"HDPE": "#1f4e79", "LDPE": "#4a90c2", "PET": "#2e7d32",
            "PP": "#e0a458", "PS": "#8e5bab", "PVC": "#c0392b"}
    # The carbonyl and deformation windows are close together, so their labels
    # are staggered vertically rather than overprinted.
    for i, (name, m) in enumerate(masks.items()):
        lo, hi = grid[m].min(), grid[m].max()
        ax.axvspan(lo, hi, color=GREY, alpha=0.16, lw=0)
        ax.text((lo + hi) / 2, 1.02 + 0.07 * (i % 2), WIN[name], ha="center", va="bottom",
                transform=ax.get_xaxis_transform(), fontsize=7, color="#444444")
    for p in POLY_ORDER:
        ax.plot(grid, Y[y == p].mean(axis=0), lw=1.0, color=cols[p], label=p)
    ax.set_xlim(grid.max(), grid.min())
    ax.set_xlabel("wavenumber (cm$^{-1}$)")
    ax.set_ylabel("mean transmittance")
    ax.legend(ncol=6, frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.32))
    fig.savefig(os.path.join(FIG, "matys9.png"))
    plt.close(fig)
    return {k: (round(float(grid[m].min()), 1), round(float(grid[m].max()), 1))
            for k, m in masks.items()}


def fig_decidability():
    """How the certain fraction and the referral rate move with the attribute set."""
    rs = pd.read_csv(os.path.join(RES, "rough_sets.csv"))
    v0 = "sens.min bins, k=4"
    d = (rs[rs.variant == v0].groupby(["attribute_set", "n_attributes"])
         .agg(gamma=("gamma", "first"), cells=("n_cells", "first")).reset_index()
         .sort_values("gamma", ascending=False).head(10))
    d["referral"] = 1 - d.gamma
    fig, axes = plt.subplots(1, 2, figsize=(PAGE, 2.9))
    lab = [" + ".join(label(a) for a in s.split("+")) for s in d.attribute_set]
    ypos = np.arange(len(d))[::-1]
    axes[0].barh(ypos, d.gamma, color=[ACC if i else WARN for i in range(len(d))][::-1])
    axes[0].set_yticks(ypos)
    axes[0].set_yticklabels(lab, fontsize=6)
    axes[0].set_xlim(0, 1)
    axes[0].set_xlabel("decided with certainty")
    axes[0].set_title("Quality of approximation by attribute set")
    for yv, g in zip(ypos, d.gamma):
        axes[0].text(g + 0.012, yv, f"{g:.3f}", va="center", fontsize=6)
    axes[1].scatter(d.cells, d.referral, s=18, color=ACC, zorder=3)
    for c, r, n in zip(d.cells, d.referral, d.n_attributes):
        axes[1].annotate(f"{int(n)} attr", xy=(c, r), xytext=(0, 5),
                         textcoords="offset points", ha="center", fontsize=6, color=GREY)
    axes[1].set_xlabel("cells in the attribute grid")
    axes[1].set_ylabel("referral rate")
    axes[1].set_title("Price of certainty: finer grids refer more")
    axes[1].grid(alpha=0.25, lw=0.5)
    for ax, letter in zip(axes, "ab"):
        ax.annotate(letter, xy=(-0.09, 1.08), xycoords="axes fraction",
                    fontsize=9, fontweight="bold", va="top")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "matys10.png"))
    plt.close(fig)
    return d[["attribute_set", "n_attributes", "gamma", "referral", "cells"]].round(4)


def fig_sensitivity():
    """Headline rule confidence under every choice the protocol left open."""
    s = pd.read_csv(os.path.join(RES, "pair_cap_sensitivity.csv"))
    conf = s[s.metric == "conf_change"]
    kinds = list(dict.fromkeys(conf.variant_type))
    fig, axes = plt.subplots(1, len(kinds), figsize=(PAGE, 2.5), sharey=True)
    axes = np.atleast_1d(axes)
    nice_kind = {"pair_cap": "pair sampling cap", "k": "threshold multiplier k",
                 "resolution": "spectral resolution"}
    for ax, kind in zip(axes, kinds):
        g = conf[conf.variant_type == kind]
        top = (g.groupby("attribute").value.mean().sort_values(ascending=False)
               .head(4).index.tolist())
        for i, a in enumerate(top):
            gg = g[g.attribute == a].sort_values("variant_value")
            ax.plot(range(len(gg)), gg.value, marker="o", ms=3, lw=1.2,
                    label=label(a), color=[ACC, "#4a90c2", GREY, WARN][i])
            ax.set_xticks(range(len(gg)))
            ax.set_xticklabels([str(v) for v in gg.variant_value], fontsize=6.5)
        ax.set_xlabel(nice_kind.get(kind, kind))
        ax.grid(alpha=0.25, lw=0.5)
    axes[0].set_ylabel("rule confidence")
    axes[-1].legend(frameon=False, fontsize=6, loc="lower right")
    spread = conf.groupby("attribute").value.agg(lambda v: v.max() - v.min())
    fig.suptitle("Rule confidence is flat across every open protocol choice",
                 fontsize=8.5, y=1.03)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "matys11.png"))
    plt.close(fig)
    return float(spread.max())


def main():
    setup()
    lo, hi = fig_entry()
    par = fig_parity()
    med_rule = fig_rules()
    best_attr, best_ovl = fig_ovl()
    fig_oc()
    wins = fig_windows()
    dec = fig_decidability()
    spread = fig_sensitivity()
    lc = fig_learning()
    n_sets = fig_shap()
    print(lc.to_string(index=False))
    print(f"entry condition {lo:.1f} to {hi:.1f}")
    print(par[["method", "balanced_accuracy", "coverage"]].round(4).to_string(index=False))
    print(f"median rule confidence {med_rule:.4f}")
    print(f"most separating attribute {best_attr} with median overlap {best_ovl:.4f}")
    print(f"distinct glass-box attribute sets across repeats: {n_sets}")
    print(f"windows {wins}")
    print(dec.to_string(index=False))
    print(f"largest confidence spread across open choices {spread:.4f}")
    print("figures written:", sorted(f for f in os.listdir(FIG) if f.startswith("matys")))


if __name__ == "__main__":
    main()

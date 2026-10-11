# Frozen analysis protocol, version 1.1

Tag: `protocol-v1.1` · Frozen: 10 October 2026 · Status: frozen before any modelling run

Amendment 1.1, made before any modelling result was produced: the pair
database of Step 3 is capped at 2000 draws per ordered parameter-value pair
within each held-fixed stratum, with a fixed seed. The full database for the
primary task would contain about 3.75 million pairs, which is intractable
across 5 folds and 20 seeds; the cap also stops a large stratum from
dominating the weighted confidence. The cap is a sampling of the pair
database, not a change to its definition: the single-parameter-difference and
different-group constraints of Eq. (A10) are enforced before sampling.
Sensitivity to the cap is reported at 500, 2000 and 5000 draws.

This protocol implements Sections 6.7 and 6.8 of the article concept
(xAI_GOZ rev. 3) with the amendments forced by the data audit of Section 0.
Nothing below may be changed once a modelling result has been produced. Any
change requires a new version tag and a restatement of which results were
produced under which version.

## 1. What the audit changed, and why

The audit of the full recyclate archive (Zenodo 21822279) answered the seven
blocking questions of concept Section 0 and invalidated the planned unit of
analysis.

| Question | Finding |
|---|---|
| Does the metadata describe only the primary dataset? | No. The 394-row, 161-column metadata table carries the matched reference-library entry alongside each spectrum. `Organization` on rows that correspond to archive spectra resolves to Cowger and Primpke (107), FLOPP (43), Aline Carvalho (25), Primpke 2018 (11), RRUFF (7). Polymer identity, morphology and colour therefore describe the library match, not an independently established property of the particle. |
| Is n spectra or particles? | The dataset's own `object_id` convention is `NAME(i/N)`, which yields 54 particles. Spectrum file names use `NAME_i(N)` and `NAME_i`; reading only the first form gives 211 particles. The metadata convention is authoritative: `CompostRastko_film` is one particle with 119 metadata rows and 92 spectrum files. |
| How many particles carry a polymer label? | 172 spectra carry a polymer label, distributed over 17 particles under the authoritative convention. |
| 137 + 27 = 164, not 183 | 164 is the count of spectra in the PE, PP and PS classes. The remainder are other polymers: PVC 5, ABS 2, SAN 1. "183" never referred to particles. |
| How are multi-spectrum particles aggregated? | Median ACO per particle, with majority vote over spectra as sensitivity. Distribution is extreme: 11 particles with one spectrum, largest with 92. |
| FTIR-Plastics: 6000 or 3000 spectra? | 3000 specimens, each measured twice: 3000 spectra at 8 cm-1 (1869 points) and the same 3000 at 4 cm-1 (3736 points). "6000" counts both resolutions. The files are in record version 10736650, not 11176740. |
| How was pass/fail assigned? | By the archive's own automated quality flags: `match_threshold` 0.8 and `signal_threshold` 4.0, exposed as `good_matches`, `good_signal`, `good_correlations`. Not by visual inspection or expert assignment. On archive rows this gives 192 pass against 14 fail. The 0.90 cut used in earlier exploratory work is not the dataset's threshold. |

Consequence, and the decision taken: 17 labelled particles cannot support
stratified group k-fold, learning curves at five subset sizes, or group
bootstrap confidence intervals, and labels that are the output of library
matching cannot support a claim of independent identity verification.
**FTIR-Plastics becomes the primary dataset for all methodological claims; the
recyclate archive becomes the applied transfer case.**

## 2. Datasets and roles

| Role | Dataset | Content | What it supports |
|---|---|---|---|
| Primary | FTIR-Plastics, Zenodo record version 10736650, CC BY 4.0 | 3000 specimens, six polymers (PET, HDPE, PVC, LDPE, PP, PS), 500 per polymer, transmittance over 399.2 to 4001.6 cm-1, 32 scans; C8 at 8 cm-1 primary, C4 at 4 cm-1 as resolution sensitivity | Accuracy parity, learning curves, rough-set bounds, SHAP stability, black-box comparison |
| Transfer case | Recyclate FTIR archive, Zenodo 21822279, CC BY 4.0 | 356 spectra in three waste streams; 54 particles, 17 with a polymer label; archive quality flags | Provenance flag (leave-one-stream-out), inspection readability, applied discussion |

Three counts in the transfer case must not be confused, because they measure
different things. The 54 particles are counted from the names of the 356
spectrum files under the `object_id` convention, and the same figure is
obtained independently as the number of unique `PARTICLE_ID` values across all
394 metadata rows. The 172 labelled spectra in 17 particles are those files
that additionally carry a polymer label through the metadata join. A third
number, 19, is the count of distinct particles among the 206 metadata rows
whose file key matches a file present in the archive; it measures the coverage
of the metadata join, not the number of particles in the dataset, and is never
used as a sample size. Leave-one-particle-out splits run over the 54 particles
for unlabelled tasks and over the 17 labelled particles for polymer identity.

Labels in the primary dataset are assigned by the specimen's known polymer
identity, independently of any spectral matching procedure. This is the
property the recyclate archive lacks and the reason for the swap.

## 3. Unit of analysis and grouping

Primary dataset: the unit is the specimen (`IDE`), one spectrum per specimen
per resolution. The archive does not document whether the 500 specimens per
polymer are 500 distinct physical items or repeated presentations of fewer
items. The available evidence against acquisition-block structure is reported
and was computed before the freeze: the mean correlation between spectra
adjacent in file order equals the mean correlation between randomly drawn
pairs within the same polymer (HDPE 0.7370 against 0.7235; LDPE 0.9009 against
0.8808; PET 0.9291 against 0.9273; PP 0.7928 against 0.7308; PS 0.6482 against
0.6540; PVC 0.6559 against 0.6552). Splits are therefore stratified by polymer
without further nesting, and this limitation is stated in the manuscript.

Transfer case: the unit is the particle (`object_id` before `(i/N)`). Spectra
of one particle never cross a train and test boundary. Particle-level results
are primary; spectrum-level results are reported as secondary for
comparability with the literature.

## 4. Noise model and the sensitivity threshold

Concept Step 2 calibrates `sens.min` from pooled within-unit replicate
standard deviation. The primary dataset has no replicate acquisitions under
fixed conditions, so the following substitution is fixed here:

1. **Primary noise estimate.** Per spectrum, the high-frequency residual left
   by a Savitzky-Golay filter (window 11, polynomial order 3). This estimates
   instrument noise under fixed conditions. Measured before the freeze: median
   0.1082 and 95th percentile 0.1593 transmittance units.
2. **Upper-bound sensitivity.** The same specimen measured at 4 and 8 cm-1
   gives a two-condition comparison; the resulting disagreement bounds the
   noise from above because it also contains the resolution effect. Reported
   as a sensitivity analysis, never as the primary floor.
3. `sens.min(f) = k · SD_max(f)` with `k = 4` as in concept Eq. A9, and the
   implied confidence recomputed and reported for the degrees of freedom
   actually available. Sensitivity over `k` in {3, 4, 5}.
4. On the transfer case the original rule applies unchanged: pooled
   within-particle standard deviation over replicate spectra, `r_min = 3`.

Entry condition, computed before the freeze on the primary dataset: between-
specimen standard deviation exceeds the instrument-noise estimate by a factor
of 25.3 to 57.2 depending on polymer. The condition is satisfied with margin.

## 5. Features

Three spectral windows, fixed a priori from polymer chemistry and confirmed to
lie inside the measured range: W1 = 2800 to 3000 cm-1 (C-H stretching),
W2 = 1650 to 1800 cm-1 (carbonyl), W3 = 1350 to 1500 cm-1 (CH2 and CH3
deformation). Per window: A = max − min, O = min, C = Pearson correlation of
the min-max normalised window against the template. Nine attributes per
spectrum.

Template: global median spectrum of the training split (primary) and
per-polymer median (sensitivity). Never per stream, never estimated on test
data.

Change scale: `Δ ln A` for amplitude, absolute `ΔC` and `ΔO`, as in concept
Table 7.

## 6. Evaluation protocol

| Element | Specification |
|---|---|
| Primary task | Polymer identity, six classes, on FTIR-Plastics C8 |
| Split | Stratified k-fold, k = 5, repeated with 20 seeds |
| Resolution sensitivity | The same protocol re-run on C4 |
| Transfer tasks | Provenance: leave-one-stream-out over the three waste streams. Readability: the archive's own quality flags, particle level, leave-one-particle-out |
| Per-subset re-estimation | Template, `sens.min`, selected window and bins are estimated on each training subset only |
| Baselines | SVM with RBF kernel, random forest, PLS-DA, and a one-dimensional convolutional network, all on identical splits with nested tuning inside each training fold |
| Learning-curve sizes | 1, 2, 5, 10, 25, 50, 75 and 100 percent of training specimens, stratified by polymer, with the absolute count per class reported alongside the percentage |
| Metrics | Balanced accuracy and Matthews correlation coefficient as primary; per-class sensitivity; accuracy only alongside the majority-class baseline; referral rate for the glass box |
| Uncertainty | Bootstrap over specimens, 1000 iterations, mean and 95 percent interval |
| Comparison test | Corrected resampled t-test (Nadeau and Bengio) across repeated cross-validation; McNemar only on a single held-out set |
| Crossover definition | The smallest training size at which the 95 percent interval of the difference in balanced accuracy between the glass box and the best black box excludes zero in favour of the glass box. Absence of a crossover is reported as a result |
| Order of analyses | Cross-validation at full size, then learning curves, then the transfer case, then the resolution sensitivity. All reported regardless of outcome |

## 7. Rule mining and rough sets

Unchanged from concept Section 6.8: pair database restricted to pairs
differing in exactly one parameter, replicate weighting `1/(r_g · r_h)`,
ternary labels against `sens.min`, weighted support and confidence, Wilson 95
percent intervals, permutation test with Benjamini-Hochberg control,
exclusivity thresholds `θ_s = 0.90` and `θ_i = 0.80`, descriptors `SR` and
`conf_all`, rough-set approximations with `α_B(X)`, `γ_B(P)`, `BN_B(X)`, bins
of width `sens.min` with ten equal-width bins as sensitivity, at most two
condition attributes, overlap coefficient by Gaussian kernel density with
Silverman bandwidth on a 512-point grid, and the three-way decision at test
time.

## 8. What is claimed and what is not

Claimed: that a rule-based representation built on named spectral attributes
reaches the accuracy of black-box classifiers on the primary task, states
formally when an individual decision is certain, and behaves as specified
under sparse training data.

Not claimed: action rules in the sense of Raś and Wieczorkowska; fitness for
purpose or contractual acceptance of any material; batch-level decisions,
which are addressed conceptually through acceptance sampling; and any claim
about the recyclate archive that depends on its polymer labels being
independent of library matching.

Added after the external review of 11 October 2026, and also not claimed: that
the rule attains a higher selective accuracy than a black box operating under
the same referral budget. It does not. Section V of the article reports the
measured comparison, and the claim is confined to the legibility of the
refusal, at a stated accuracy cost.

## 9. Deviations from this protocol, as executed

This section is appended to the frozen protocol rather than altering it. The
protocol text above is unchanged; what follows records every place the
execution departed from it, with the reason.

1. **Repetitions of cross-validation.** Section 5 specifies twenty repetitions
   of stratified five-fold cross-validation. One repetition was completed, at
   a cost of 43 minutes on a single core. Intervals are therefore computed
   across the five folds of one repetition, and the comparison with the
   strongest baseline is a paired t-test across those folds instead of the
   corrected resampled test, which requires repeated runs.
2. **Learning-curve sizes.** Section 6 specifies training fractions of 1, 2,
   5, 10, 25, 50, 75 and 100 percent. The run completed covers 1 to 50
   percent. The conclusion that no crossover occurs is therefore established
   only up to half of the training data.
3. **Attribute-selection criterion.** The protocol does not state which of two
   defensible criteria selects the condition attributes. The validated
   classifier uses inner cross-validated balanced accuracy; Table III of the
   article reports the pair that maximises the training-set quality of
   approximation. The two disagree. Table X settles the disagreement out of
   fold in favour of the cross-validated criterion, and both are reported.
4. **Noise estimate in the entry condition.** The ratio quoted in an early
   draft, 25.3 to 57.2, came from an ad hoc calculation. The released figure,
   25.5 to 56.9, is produced by `code/09_entry_condition.py`, which is now
   part of the deposit. The difference is the estimator, not the data.
5. **Scaling inside the baseline pipelines.** Feature scaling is estimated on
   the whole training fold before the inner split used for hyperparameter
   search. Outer estimates are unaffected, since they are computed on folds
   the scaler never saw, but the tuned hyperparameters are mildly optimistic.
6. **Analyses added after the protocol was frozen.** The selective comparison
   at matched coverage, the out-of-fold comparison of candidate attribute
   pairs, the out-of-fold effect of the binning rule and the per-polymer
   coverage breakdown were added in response to external review. They are
   confirmatory analyses of an already-fitted model: no model was refitted
   with knowledge of a test fold, and every one of them uses the folds of the
   original run.

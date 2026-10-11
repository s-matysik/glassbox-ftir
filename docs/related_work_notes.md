# Related-work notes — working material

Corpus: 56 verified entries, built 10 Oct 2026 from OpenAlex (API key auth, no mailto) and Crossref.
Keys below are the `key` column of `references.csv` and the BibTeX keys of `references.bib`. Every key is built from
the **first author of the retrieved record** plus its year, so the key never implies an author who is not on the paper.
Every DOI was resolved against Crossref `/works/{doi}` (or, for the four arXiv DOIs, against OpenAlex, since
`10.48550/*` is DataCite-registered and not indexed by Crossref). Retraction was checked three ways for every
DOI entry: the OpenAlex `is_retracted` flag, the presence of a Crossref `update-to` of type `retraction` on the
record itself, and a Crossref `filter=updates:<doi>` query that returns any notice updating that DOI.
**No entry in this corpus is retracted and none carries a retraction notice.**

Per-section counts: A 15, B 10, C 7, D 8, E 8, F 8.

---

## A. Machine learning for FTIR polymer identification (15 entries)

The retrieved literature is dominated by supervised classifiers on spectral vectors, and the families split
cleanly: discriminant/chemometric (PLS-DA and SIMCA in `dasilva2020fpa`), tree ensembles (random decision
forests in `hufnagl2021mufti`, random forest alongside decision trees in `liu2023blended`, extremely randomised
trees in `rosalesmartinez2025`), kernel and instance methods (SVM, k-NN in `villegascamacho2025norm` and `yan2022ensemble`,
subspace k-NN in `williams2025`), and deep models (1-D CNN in `liu2023blended`, `jiang2021cnnatr` and `smolen2025pnas`,
a 2-D CNN on spectral images in `liu2023blended`, a dense feed-forward network in `herb2025humans`, a GAN-augmented
MLP in `liu2026sparse`). No retrieved paper uses a transformer on FTIR spectra for polymer identification,
which is worth stating explicitly in the manuscript rather than leaving as an unexamined gap.
**The split protocol is the weak point of this entire baseline and this is the argument the manuscript should
make.** Only three of the fifteen describe a split in retrievable text, and all three are random at the spectrum
level, not at the specimen or batch level: Monte Carlo cross-validation in `hufnagl2021mufti`, 10-fold
cross-validation in `villegascamacho2025norm`, and stratified repeated cross-validation plus a hold-out test set in
`rosalesmartinez2025`. Two papers do better by evaluating across acquisition conditions rather than across a random split
— `smolen2025pnas` trains on pristine spectra from 45 manufactured samples and tests on a noisy spiked-filter set,
and `williams2025` trains on 1800 spectra (900 µFTIR, 900 µRaman) and then tests against the independent external
FLOPP-E and SLOPP-E libraries, reporting an accuracy loss on real-world data. For the remaining ten the split is
not stated in the metadata that could be retrieved and is recorded as unknown; `kedzierski2019rf` returns no
abstract at all through OpenAlex and is recorded as unknown on every attribute except that it is an ML paper.
Reported accuracies cluster absurdly high (96.4 and 97.4 percent in `liu2023blended`, 99.5 percent in
`tian2023smalldata`, a mean F1 of 0.99995 and perfect hold-out classification in `rosalesmartinez2025`), which is
precisely what a spectrum-level random split produces when replicate spectra of the same specimen land on both
sides of the split; the manuscript should use the split column, not the accuracy column, as the comparison axis.
Uncertainty for an individual decision is almost absent. `herb2025humans` is the clearest positive case — it
returns a conditional class distribution over 16 microplastic categories per spectrum and uses disagreement with
human annotators to expose human labelling errors. `smolen2025pnas` gives a partial case: the similarity embedding
with an SVC separates compositions that were never in training, which is open-set detection rather than a
calibrated per-decision confidence. `williams2025` tests replicate consistency and reports degradation on external
libraries but produces no per-spectrum confidence. `bendell2024unc` addresses identification uncertainty only
through the choice of spectral library, which is current practice. Everything else reports aggregate accuracy,
precision, recall and F1. **No retrieved paper produces a per-decision accept/reject/defer gate**, which is the
space the manuscript's three-way decision rule occupies.
Three entries bear on the small-sample and class-imbalance regime that a recyclate study lives in:
`tian2023smalldata` compares models across data-set sizes with synthetic augmentation, `yan2023` uses XAI
attribution and Gaussian mixtures to locate the informative fingerprint regions and drive augmentation, and
`liu2026sparse` is the only retrieved item framed directly on post-consumer recyclate spectra with sparse,
imbalanced classes.

## B. Explainable AI and the interpretable-by-design argument (10 entries)

`rudin2019` is the load-bearing citation: for high-stakes decisions, build a model that is interpretable rather
than explaining a black box afterwards. `lipton2018` is needed immediately after it, because it shows that
"interpretability" is not one property — the manuscript must say which kind its amplitude-correlation-offset
attributes deliver (simulatability of a named band, not post-hoc attribution). Cite `lundberg2017` and
`ribeiro2016` as the methods being criticised, so the critique has a defined target rather than attacking a
strawman. The critique itself rests on four independent results that should be cited together: `slack2020` shows
LIME and SHAP can be deliberately fooled by an adversarially scaffolded model, so an explanation can be made to
hide the real decision rule; `kumar2020` shows Shapley-value attributions depend on a reference distribution
whose choice is rarely justified and whose assumptions are usually violated; `alvarezmelis2018` and `ghorbani2019`
show that explanations move under input perturbations that leave the prediction unchanged, which is the
instability claim in its strongest form; and `adebayo2018` shows several saliency methods fail sanity checks
because they are insensitive to model parameters and to the labels the model was trained on. `barredoarrieta2020`
supplies the taxonomy and vocabulary (transparent model versus post-hoc explanation) in which the manuscript
should locate its own method. Note the asymmetry worth exploiting: `yan2023` in section A uses post-hoc XAI to
*find* the informative fingerprint regions and then throws the explanation away; the manuscript's design keeps
those regions as the model.

## C. Association rules and action rules (7 entries)

`agrawal1993` supplies the support/confidence definitions and must be cited where those two numbers first
appear. `han2000` is the citation for the mining implementation and its complexity if FP-growth or a descendant
is used. The action-rule line proper runs `ras2003dear` (the DEAR system, extended action rules), `tzacheva2005`
(action-rule mining with the cost and feasibility of the recommended attribute change) and `ras2006simpl` (the
single-step simplified strategy) — `ras2006simpl` is the formulation closest to the pair-based construction the
manuscript uses, because it builds rules from objects that agree on stable attributes and differ in one flexible
one, which is exactly what cross-sectional spectral data permits. `tzacheva2005` is the right citation when the
manuscript argues which attribute changes are actually actionable in a recycling process, since not every
spectral change corresponds to an intervention a sorter can make. `geng2006` is needed wherever mined rules are
ranked, because support and confidence alone do not order rules usefully. `benedict2025ar` keeps the section
from looking like a 2000s-only reference list and shows the line is still active.

## D. Rough sets and three-way decisions (8 entries)

`pawlak1982` is the founding citation for indiscernibility and the lower/upper approximation; `pawlak2007rudiments`
is the better source for the actual definitions used (boundary region, accuracy and quality of approximation),
and the manuscript should cite the pair rather than 1982 alone. The three-way line is `yao2010prob` for the
derivation of accept, reject and defer regions from loss functions in probabilistic rough sets, `yao2012outline`
for the general trisecting-and-acting framing independent of rough sets, and `yao2011superior` for the explicit
argument that three regions beat two when misclassification costs are asymmetric — which is the direct
justification for a deferment region in an acceptance decision, where accepting a contaminated batch and
rejecting a good one carry different costs. `yao2019trilevel` is the current formulation and should be cited so
the section is not dominated by 2009–2012 work. `ziarko1993` supplies variable-precision rough sets, i.e. the
mechanism by which the inclusion threshold — and hence the width of the deferment region — becomes a tunable
parameter rather than a fixed strict-inclusion rule. `abdar2021skin3wd` is the one applied precedent in the set:
three-way decisions on top of Bayesian deep learning producing a deferral class in skin-cancer classification,
and it is the closest published analogue to what the manuscript proposes for spectra.

## E. Acceptance sampling, quality control and EN standards for recyclates (8 entries)

The standards side has moved recently and the manuscript must cite current editions. `en15347_1_2024`
(EN 15347-1:2024, general characterisation of sorted plastics wastes) replaced EN 15347:2007 and defines what a
supplier must declare about a sorted stream — this is the document that frames the acceptance problem.
`en15343_2007` remains the reference for traceability, conformity assessment and recycled content.
`en15344_2025` covers PE recyclate characterisation. Only one national-body record evidences the 2025 edition
(SIS, SS-EN 15344:2025); a 2022 Spanish adoption and 2007 Swedish and Estonian adoptions were also returned, so at
least three editions exist, but the edition number and supersession chain are unverified. Confirm against the CEN
catalogue before submission. `en18065_2025` (classification of recycled
plastics and data-quality levels) is the most useful recent item for this manuscript because it ties recyclate
classification to data quality for digital trading and product-passport use, which hands section F a
standards-based version of the information-asymmetry argument. On the sampling side, `iso2859_1_2026` gives the
lot-by-lot AQL system that a three-way rule should be compared against; a 2026 edition is the latest one evidenced, while
much of the literature still cites ISO 2859-1:1999 with Amd 1:2011, so check which edition applies, and `schilling2017` supplies the
vocabulary — OC curve, producer's risk, consumer's risk — in which the deferment region should be described.
`eriksen2018quality` and `eriksen2019closing` are the empirical anchors: they show for household streams and for
PET, PE and PP specifically that contamination and material properties, not collection rate, decide whether a
closed loop is achievable, which is the practical reason per-batch characterisation is worth automating.

## F. Circular economy and recyclate quality as an economic problem (8 entries)

`akerlof1970` is the theoretical core: when quality is unobservable at the point of sale and information is
asymmetric, buyers price for the average, good material leaves the market, and the market thins. The manuscript
should state this as the economic reading of its measurement problem — a glass-box, per-batch, auditable
classification is a quality-signalling device, not only a sensor method. `nicolli2012failures` applies the
market-failure analysis specifically to recycling markets and discusses technological innovation as a remedy,
which is the bridge from Akerlof to this domain. `huysveld2022substitut` is the strongest empirical support
available here: it separates technical from market substitutability of recycled materials and quantifies the gap,
i.e. it shows directly that a technically adequate recyclate can still fail to be accepted. `gall2020pickers`
supplies field evidence that recyclate quality rather than supply volume is the binding constraint, and
`hazen2017perceived` shows that buyer-side perceived quality is a measurable construct separable from actual
quality — together they make the asymmetry concrete rather than theoretical. `iacovidou2019quality` gives the
typology for distinguishing technical quality from market-accepted quality, which is the vocabulary the section
needs. `schyns2020mech` provides the mechanical-recycling background — contamination and polymer incompatibility —
explaining *why* quality is genuinely uncertain at the point of sale rather than merely poorly communicated.
`kusenberg2022opps` shows the same feedstock-quality constraint one step downstream, in pyrolysis-based chemical
recycling; it is a pyrolysis paper and must not be cited as a mechanical-sorting reference.

---

## Cautions for the author

1. Section A's `kedzierski2019rf` has no abstract in OpenAlex; classifier family, dataset size, split and
   uncertainty reporting are all recorded as unknown. If it is cited as a baseline, read the paper before
   describing its protocol.
2. Ten of fifteen section-A entries have an unknown split. Do not write "most studies use random splits" —
   write that the split is not reported for most of them, which is itself the reviewable finding.
3. `iso2859_1_2026`: only ISO catalogue titles and URLs were retrieved, never page bodies, so the edition number, the publication
   stage and the status of the earlier ISO 2859-1:1999 are unknown. An earlier draft of this corpus cited the 1999 edition as
   current on no evidence; confirm the applicable edition on iso.org before citing.
4. `en18065_2025` rests on a single catalogue record; confirm the exact title and publication month against the
   CEN catalogue before submission. The CEN webstore is not reachable from the environment this corpus was
   built in.
5. Polymer-specific recyclate standards: the only editions actually evidenced by the searches behind this corpus
   are EN 15344:2025 (PE, one record), EN 15348:2024 (PET, one record) and EN 15342:2025 (PS, one record). Nothing
   was retrieved for the current edition of EN 15345 (PP) or EN 15346 (PVC) — treat those as unknown. Only
   EN 15344:2025 is carried as an entry; if the manuscript names any of the others, verify its edition first.
6. Four section-B entries (`kumar2020`, `alvarezmelis2018`, `lundberg2017`, `adebayo2018`) carry arXiv DOIs because
   the venue (ICML, NeurIPS) registers no DOI. They resolve and are indexed in OpenAlex. If the editor objects
   to preprint citations, replace them with the proceedings page references, which exist but have no DOI.
7. The corpus contains no transformer-based FTIR polymer classifier because none was retrieved. State the
   absence; do not fill it with a Raman or NIR transformer paper and describe it as FTIR.

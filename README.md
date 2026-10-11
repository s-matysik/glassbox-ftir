# Glass-box FTIR classification of plastics

Replication material for the IEEE Access article *A Glass-Box Classifier for
FTIR Identification of Plastics: Accuracy Parity on Accepted Specimens With an
Explicit Referral Rate*, together with the interactive documentation site.

Every number quoted in the manuscript and its supplement is produced by a named
script in `code/` from a named table in `results/`. Nothing is typed by hand, so
a number in the text cannot drift away from the computation behind it.

## Reproducing the article, in order

```
python -m pip install -r requirements.txt
python code/00_download_data.py        # FTIR-Plastics, Zenodo 10736650, CC BY 4.0
python code/09_entry_condition.py      # Fig. 2
python code/10_cv_parity.py            # Tables II and S1, Fig. 3   (about 43 min)
python code/11_learning_curves.py      # Table VI, Fig. 6
python code/12_rules_rough.py          # Tables III, IV, V; Figs. 4, 5, 8, 9
python code/13_shap_oc.py              # Tables VII, VIII; Figs. 7, 10
python code/24_parity_test.py          # the significance and non-inferiority tests
python code/21_paper_figures.py        # every figure, from the tables above
python code/23_supplement.py           # the supplement and its four figures
python code/20_export_web.py           # the JSON the site reads
python code/22_assemble_manuscript.py && python code/build_manuscript.py
```

Verification scripts: `code/test_ara_aco.py` checks the pipeline against
synthetic spectra with planted structure, `code/test_equations.py` checks that
every manuscript equation becomes a real Word equation, and `code/check_site.js`
loads the site in a simulated browser and visits every view.

## What is here

```
code/        the method, the analyses, the document and figure builders, the tests
results/     every derived table the manuscript and supplement quote
figures/     the article figures (matys*.png) and the supplement figures (supp*.png)
docs/        the manuscript sources, the reference list and the documentation site
PROTOCOL_v1.1.md   the analysis protocol, frozen before any model was fitted
```

The raw spectra are not redistributed here: `code/00_download_data.py` fetches
them from the Zenodo record, which is the citable source.

## The IEEE Access article

The chain that builds the document:

```
code/ara_aco.py            the method: ACO attributes, threshold calibration,
                           pair rules, rough sets, three-way decision
code/test_ara_aco.py       verification on synthetic spectra with planted structure
code/10_cv_parity.py       five-fold cross-validation against four baselines
code/11_learning_curves.py learning curves (env: LC_SIZES, LC_SEEDS, LC_BASELINES)
code/12_rules_rough.py     rule matrix, rough sets, overlap, sensitivities
code/13_shap_oc.py         attribution stability and acceptance sampling
code/20_export_web.py      the JSON the documentation site reads
code/21_paper_figures.py   figures/matys1..10.png, each from a table in results/
code/22_assemble_manuscript.py  merges docs/manuscript_part*.json, builds the
                           tables from results/ and resolves [@key] citations
code/build_manuscript.py   writes docs/IEEE_Access_manuscript.docx inside the
                           official Access template, reusing its styles
```

Prose lives in `docs/manuscript_part1..3.json`; tables and figures are
generated. Editing a number in the prose is therefore not possible: change the
analysis, re-run, and the manuscript follows.

**Windows note.** The sandboxed Windows host refuses to create the named pipes
that `multiprocessing` needs, which makes `n_jobs=-1` hang rather than fail.
Run the heavy scripts with `N_JOBS=1` and `JOBLIB_MULTIPROCESSING=0`.

## Documentation site

`docs/` is a self-contained site published through GitHub Pages from this
folder on the default branch, and live at
<https://s-matysik.github.io/glassbox-ftir/>. It has ten views, including a
spectrum explorer, the rough-set decision view, the rule matrix, the learning
curves and an acceptance-sampling view, with cross-filtering by polymer and
per-chart download of the image and of the data behind it. It carries no
external dependency, so it also works from a downloaded copy of the Zenodo
record. Serve it over HTTP rather than opening `index.html` directly, because
browsers block local fetch requests:

```
cd docs && python -m http.server
```

Regenerate its data with `python code/20_export_web.py` after any analysis
change.

## Honest limitations

These belong in any manuscript built on this material.

- The action rules in concepts 2 and 3 are mined over object pairs and are
  symmetric by construction: the same confidence appears in both directions.
  They are evidence of association under held-fixed covariates, not of
  direction. Pairs are also not independent, since one object enters many
  pairs, so the pair count is not a sample size.
- Concept 3 has a single reporting year, so no temporal validation is possible,
  and large strata are subsampled to 260 records before pairing.
- Concept 4 rests on 50 buildings at one site, of which 32 are offices. The
  entry condition and the channel analysis use 144 buildings; the
  identification contrast uses the smaller rated subset. Why most buildings
  carry a non-numeric rating is not characterised.
- In concepts 1 and 2, 138 of 164 spectra come from one recovery matrix, and
  the identification test is binary, polyethylene against polypropylene.
- The validation deposit referenced in concept 1 covers different polymer
  classes, so it validates the procedure and not the methyl-band result.
- The required-reliability expression is derived here rather than taken from
  the literature, and should be checked against prior work before it is
  claimed as novel.

## Data sources

| Folder | Zenodo record | Licence |
|---|---|---|
| `data_raw/ftir` | 10.5281/zenodo.21822279 | CC BY 4.0 |
| `data_raw/industrial_co2` | 10.5281/zenodo.20767361 | CC BY 4.0 |
| `data_raw/choice_experiment` | 10.5281/zenodo.21740286 | CC BY 4.0 |
| `data_raw/buildings` | 10.5281/zenodo.3887306 | CC BY 4.0 |

Each folder carries a `zenodo_record.json` with the title, DOI, licence,
publication date and file listing as retrieved.

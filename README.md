# Brassica-Subgenome-Dominance

Two axes of subgenome dominance across the *Brassica* U's triangle: the structural asymmetry
between subgenomes is inherited from the diploid progenitors, while the expression bias between
homoeologs is largely rebuilt inside the hybrid. This repository accompanies the paper and
contains every script that turns the processed counts into the figures and tables of the
manuscript.

<p align="center"><img src="figures/Fig1_hero.png" alt="Evolutionary relationships of the U's triangle species and their subgenome composition" width="100%"></p>

The three allotetraploids of the U's triangle arose from independent allopolyploidizations of the
same three diploid genomes, so each genome occurs in two hybrids with a different partner each
time. That is what lets a feature of a *genome* be told apart from a feature of a *pairing*.

## What the analysis concludes

- Three pairings resolve to one order: **B is dominant wherever it occurs** and **C is
  non-dominant wherever it occurs**, across five fold-change thresholds and both mapping-quality
  settings. A is intermediate and tissue-dependent.
- **Structural asymmetry is inherited.** Each subgenome carries its own gene-proximal LTR density
  and its own K~a~/K~s~ into both hybrids containing it, agreeing to within 4-9%. The difference
  between the two *napus* subgenomes (1.46-fold) is the difference their progenitors already had
  (1.45-fold).
- **Expression bias is largely rebuilt.** The parental difference explains only 8-21% of the bias
  inside the hybrid, parental bias is dampened (slope 0.23-0.44), and allopolyploidy erases more
  bias than it creates. Stigma remodels more than pollen in every species.
- **No structural proxy predicts which homoeolog is suppressed gene by gene.** Flanking LTR
  density, centromere distance and K~a~/K~s~ (after controlling for the expression-rate
  anticorrelation) all give |rho| <= 0.05, while the same measures separate subgenomes 1.5- to
  2.1-fold on average.

## Layout

- `SUPPLEMENTARY.md` supplementary figures and tables, with the tables inline
- `figures/`
  - `style.py` one shared matplotlib style; fixed colour per genome (A blue, B vermilion,
    C green) so a genome is the same colour in every figure, tissue carried by marker
  - `make_fig2_expression_landscape.py` PCA, scree and replicate correlation, recomputed from the
    counts
  - `make_fig3_scoreboard.py` the dominance scoreboard with exact binomial tests and Wilson
    intervals, at five thresholds
  - `make_fig4_inheritance.py` allotetraploid bias against progenitor bias; slope, R2 and the
    four-way split
  - `make_fig6_genelevel.py` the three paired structural tests, with the subgenome-level contrast
  - `Fig1_evolutionary_relationships.png` the schematic (also `Fig1_hero.png`, downscaled)
- `tables/`
  - `make_tables.py` builds all four main tables as TSV and markdown
- `pipeline/` the analysis scripts, numbered in run order; see `pipeline/README.md`

## Reproducing the figures and tables

```bash
cd figures && python3 make_fig2_expression_landscape.py   # and make_fig3, make_fig4, make_fig6
cd ../tables && python3 make_tables.py
```

Both read the analysis outputs by relative path (`../expression_rebuild/...`), so that directory
must be present alongside this repository, or the paths adjusted at the top of each script. Every figure writes a `.png` and a `.pdf`; PDF text is editable (fonttype 42).

## Dependencies

Python 3 with numpy, pandas, scipy, scikit-learn, statsmodels, matplotlib and Pillow; R with
DESeq2 for the bias and inheritance analyses. The evolutionary metrics additionally use SynMap
output from the CoGe platform, RepeatMasker with a common TE library for the LTR annotations, and
BWA-MEM and featureCounts upstream.

## Data

Processed count matrices and the conserved homoeologous gene set are produced by the pipeline in
`pipeline/`. Raw sequencing reads are deposited under accession *[to be completed]*.

## Repository and citation

```bash
git clone https://github.com/USask-BINFO/Brassica-Subgenome-Dominance.git
```

The supporting figures and tables are in [`SUPPLEMENTARY.md`](SUPPLEMENTARY.md). The manuscript
itself is not distributed here while it is unpublished; every value it reports is reproduced by
the scripts in `pipeline/` and tabulated in `SUPPLEMENTARY.md`.

Two further repositories hold analysis this one calls out to:

| what | where |
|:-|:-|
| SynMap parsing, Gaussian-mixture fitting of K~s~, fractionation ratio, centromere distance | https://git.cs.usask.ca/lij313/analyzecoge |
| sequence similarity and K~s~ for orthologous gene pairs | https://git.cs.usask.ca/lij313/evolv |

Citation details will be added once the paper is published.

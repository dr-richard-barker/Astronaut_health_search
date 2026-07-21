# Astronaut Health Search — spaceflight oncogenic biomarkers and their reversal by food-derived flavonoids

**A from-scratch R pipeline for "opposite-forcing" drug discovery against spaceflight oncogenic biomarker signatures, plus the supporting JAXA6 human-cell RNA-seq analysis.**

> Independent research. Not affiliated with or endorsed by NASA or JAXA.

---

## Overview

Microgravity and ionising radiation are independently linked to hallmarks of
early oncogenesis. This repository asks whether a reproducible **oncogenic
biomarker signature** emerges across human and rodent spaceflight transcriptomes,
and whether existing drug or nutraceutical perturbation signatures can **reverse**
it ("opposite forcing").

The main analysis (`astronaut-opposite-forcing/`) is a from-scratch R pipeline
that:

1. Re-discovers a spaceflight oncogenic biomarker signature across **~23 NASA
   OSDR RNA-seq studies** (11 human, 12 rodent) using strict DE thresholds
   (padj < 0.05, |log2FC| > 1) and a random-effects meta-analysis (`metafor` REML).
2. Screens **420 LINCS L1000 drug perturbation signatures** for transcriptomic
   reversal via a tau-analog connectivity score.
3. Annotates top reversal hits with target/indication context (Broad Drug
   Repurposing Hub, PrimeKG).
4. Maps the signature to disease phenotypes (PrimeKG, Open Targets).
5. Flags nutraceutical candidates among the top reversal compounds.

## Key results

| Finding | Value |
|--------|-------|
| Meta-analysis genes scored | 40,211 (3,668 significant) |
| Oncogenic biomarker genes (∩ 11,241-gene oncogene union) | **1,610** |
| Validation vs. independent Python meta-signature (S2) | Spearman ρ = 0.60, 76% sign agreement |
| Compounds with maximal reversal (tau = −100) | 21 |
| Lead nutraceutical hits | **apigenin, luteolin, licochalcone A** (food-derived flavonoids) |
| Biomarker genes directly reversed by ≥1 flavonoid | 242 / 1,540 (15.7%); 20 reversed by all three |
| Disease enrichment | sarcomatoid carcinoma, coronary artery disease, prostate cancer, hereditary breast/ovarian cancer, fatty liver disease |

Full methods, figures, and the complete results narrative are in
[`astronaut-opposite-forcing/manuscript/manuscript.md`](astronaut-opposite-forcing/manuscript/manuscript.md),
with numbered supplementary tables (T3–T13) in
[`astronaut-opposite-forcing/tables/`](astronaut-opposite-forcing/tables/).

## Repository structure

```
.
├── astronaut-opposite-forcing/   # MAIN analysis
│   ├── manuscript/manuscript.md   # Full manuscript (prepared for npj Microgravity)
│   └── tables/                    # Supplementary tables T3–T13
│       ├── T3_meta_signature.tsv
│       ├── T4_oncogene_intersection.csv       # the 1,610 biomarker genes
│       ├── T6_lincs_reversal_scores.csv
│       ├── T7_top_opposite_forcing.csv        # the 21 reversal compounds
│       ├── T11_disease_enrichment_opentargets.csv
│       ├── T12_nutraceutical_flags.csv
│       ├── T13_gene_level_reversal.csv
│       └── ... (T5, T8–T10, T12b/c, T13b)
├── data/                          # Supporting JAXA6 human-cell RNA-seq inputs
│   ├── JAXA6_FL_vs_Pre_diff_expression_all_comparisons.csv
│   ├── JAXA5_venn_for_metaanalysis.csv
│   ├── JAXA_Matrix_TGB_050_..._targets_.csv
│   ├── JAXA6_TGB_050_..._SEM_4iDEP.csv.zip
│   └── JAXA_RPM_human_cells_RNAseq/
├── figures/                       # iDEP / KEGG figures (JAXA6 analysis)
│   ├── iDEP_DRB_JAXA6_V1.svg
│   ├── KEGG_fl_vs_GC.svg
│   └── Jaxa6_iDEP_figures_August/
├── execution_trace/PLAN.md        # Gene-level flavonoid-reversal analysis plan
├── LICENSE                        # CC0-1.0
└── README.md
```

## Data sources

- **NASA Open Science Data Repository (OSDR / GeneLab):** https://osdr.nasa.gov
  (cite the individual OSD accessions listed in the manuscript).
- **LINCS L1000** drug perturbation signatures (via GEO).
- **JAXA** human-cell RNA-seq (JAXA5/JAXA6 series) — supporting `data/` and
  `figures/`.
- Annotation: MSigDB C6, KEGG, Broad Drug Repurposing Hub, PrimeKG, Open Targets.

## Reproduction

The opposite-forcing pipeline is implemented in R; see the Methods section of
`astronaut-opposite-forcing/manuscript/manuscript.md` for package versions,
thresholds, and the connectivity-score formulation. The JAXA6 supporting
analysis was produced with [iDEP](http://bioinformatics.sdstate.edu/idep/).

## Code availability

The R pipeline — meta-analysis, LINCS L1000 opposite-forcing screen, and
gene-level flavonoid reversal — is in
[`astronaut-opposite-forcing/scripts/`](astronaut-opposite-forcing/scripts/),
with full methods in
[`astronaut-opposite-forcing/manuscript/manuscript.md`](astronaut-opposite-forcing/manuscript/manuscript.md).
The independent Python cross-study meta-signature (S2) used for validation lives
in the companion repository
[**astronaut-oncogene-biomarkers**](https://github.com/dr-richard-barker/astronaut-oncogene-biomarkers).

## Requirements

> **TODO:** pin exact package versions (add an `renv.lock` or the output of
> `sessionInfo()`).

- **R** (≥ 4.x) — the core opposite-forcing pipeline.
- **Key R packages** named in the Methods: `DESeq2` (per-study differential
  expression), `metafor` (REML random-effects meta-analysis), `biomaRt`
  (rodent → human ortholog mapping).
- **Reference resources:** LINCS L1000 drug perturbation signatures, MSigDB C6,
  KEGG, Broad Drug Repurposing Hub, PrimeKG, Open Targets.
- **Reproducibility:** add a pinned environment file so reviewers can restore the
  exact package set.

## Related work

- [**astronaut-oncogene-biomarkers**](https://github.com/dr-richard-barker/astronaut-oncogene-biomarkers)
  — the companion Python cross-study meta-signature (S2) used here for
  independent validation.

## License

**CC0 1.0 Universal** (public domain dedication) — see `LICENSE`. OSDR/LINCS/JAXA
source data remain subject to their own terms; cite the original accessions.

## Citation

> **TODO:** add co-authors (ORCID), affiliation, and DOI once available.

Barker, R. et al. (2026). *Spaceflight oncogenic biomarkers and their
transcriptomic reversal by food-derived flavonoids: an R pipeline for
opposite-forcing drug discovery.* Manuscript in preparation (targeted at *npj
Microgravity*). Correspondence: admin@cosecloud.com

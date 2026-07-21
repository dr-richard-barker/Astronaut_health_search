# Plan: Gene-level reversal of oncogenic biomarkers by luteolin, apigenin, and licochalcone A

## Summary

Identify which of the 1,610 oncogenic biomarker genes (T4) are directly reversed at the gene level by the three flavonoid LINCS signatures (GSE128097_1, GSE128097_2, GSE137934_1), using the LINCS z-score matrix. A gene is "directly reversed" when its spaceflight meta-log2FC direction is opposite to its LINCS compound z-score direction, with |z| ≥ 1.5 (moderate effect). Then create a manuscript.md for the repo and update README.md with the gene-level reversal findings.

## Data confirmed (Phase 1)

- **T4_oncogene_intersection.csv**: 1,610 genes with columns `gene, class, source, mean_lfc, padj, n_studies, n_species, direction_concordance`. Direction: 502 up-regulated, 1,108 down-regulated in spaceflight.
- **LINCS human_geo_sigs.tsv** (1.6 GB): 4,269 signatures (rows) × 35,238 genes (columns), signed z-scores. Three flavonoid signatures confirmed present:
  - `GSE128097_1` — apigenin+luteolin co-treatment, tau=-100 (reversal)
  - `GSE128097_2` — apigenin+luteolin second sub-signature, tau=+8.7 (NOT a reversal at the WCS level, but included per user request for per-gene comparison)
  - `GSE137934_1` — licochalcone A, tau=-100 (reversal)
- **T4 → LINCS overlap**: 1,540 of 1,610 T4 genes are present in the LINCS gene columns (70 missing, likely non-protein-coding or unmapped IDs).
- **Z-score distributions** (measured):
  - GSE128097_1: range [-7.5, 8.5], sd=0.82; among T4 genes: 433 at |z|≥1, 244 at |z|≥1.5, 140 at |z|≥2
  - GSE128097_2: range [-7.0, 9.4], sd=0.97; among T4: 494 at |z|≥1, 263 at |z|≥1.5, 159 at |z|≥2
  - GSE137934_1: range [-7.6, 8.3], sd=0.68; among T4: 358 at |z|≥1, 202 at |z|≥1.5, 114 at |z|≥2
- **Reversal ratios** (opposite sign / total with |z|≥threshold):
  - At |z|≥1.5: 0.52–0.55 (slight reversal enrichment)
  - At |z|≥2: 0.52–0.56
  - At |z|≥3: 0.60–0.70 (strong enrichment, but only ~27–39 genes)

## Methodology decisions (evidence-based)

1. **Signature scope**: All 3 flavonoid signatures (GSE128097_1, GSE128097_2, GSE137934_1), per user request. Report per-signature results so GSE128097_2 (non-reversal at WCS level) can be compared gene-by-gene.

2. **Z-score threshold**: **|z| ≥ 1.5** as primary (moderate LINCS effect; ~130–140 reversed genes per signature; reversal ratio 0.52–0.55). Report |z|≥2 and |z|≥3 as sensitivity thresholds in the same table.

3. **Significance filter**: A gene is "directly reversed" if ALL of:
   - (a) present in T4 (1,610 oncogenic biomarkers)
   - (b) present in the LINCS signature (non-zero z-score)
   - (c) `sign(mean_lfc) * sign(z) < 0` (opposite direction)
   - (d) `|z| ≥ 1.5`
   The full ranked table (all 1,540 overlapping T4 genes with their z-scores and reversal status at each threshold) is also output so readers can apply their own cutoff.

4. **Reversal strength metric**: `reversal_score = -sign(mean_lfc) * z * |mean_lfc|` — positive when reversed, larger when both the spaceflight effect and the compound counter-effect are strong. Used for ranking.

## Implementation

### New script: `scripts/07_gene_level_reversal.R`

**Inputs**:
- `results/tables/T4_oncogene_intersection.csv` (1,610 genes)
- `/mnt/datalake/LINCS1000/RNAseq_transcriptomics_signatures/human_geo_sigs.tsv` (1.6 GB)

**Process**:
1. Load T4; load LINCS matrix (data.table::fread, ~30s on heavy machine).
2. Identify the 3 flavonoid signatures by grep on row names.
3. For each signature, extract the z-score vector over all 35,238 genes.
4. Intersect with T4 genes (expect ~1,540 overlap).
5. For each overlapping T4 gene, compute:
   - `spaceflight_dir` = sign(mean_lfc)
   - `compound_z` = LINCS z-score
   - `compound_dir` = sign(z)
   - `reversed` = (spaceflight_dir * compound_dir < 0)
   - `reversed_z1.5`, `reversed_z2`, `reversed_z3` = reversed AND |z| ≥ threshold
   - `reversal_score` = -sign(mean_lfc) * z * |mean_lfc|
6. Bind per-signature results into a long table + a wide table (genes × signatures).

**Outputs** (to `results/tables/`):
- **T13_gene_level_reversal.csv** — long format: one row per (gene × signature), columns: gene, class, source, mean_lfc, padj, n_studies, direction_concordance, signature, compound_label (apigenin+luteolin_1 / apigenin+luteolin_2 / licochalcone_A), compound_z, reversed, reversed_z1.5, reversed_z2, reversed_z3, reversal_score. ~4,620 rows (1,540 × 3).
- **T13b_gene_level_reversal_summary.csv** — per-signature counts: n_T4_overlap, n_reversed_z1.5, n_reversed_z2, n_reversed_z3, n_concordant_z1.5, reversal_ratio_z1.5, etc. 3 rows.
- **T13c_reversed_genes_union.csv** — the union of genes reversed at |z|≥1.5 by ANY of the 3 signatures, with per-gene: which signature(s) reversed it, spaceflight direction, max |z| across signatures, max reversal_score. This is the headline "genes directly reversed by flavonoids" list.

**Outputs** (to `results/figures/`):
- **F11_gene_level_reversal.svg/png** — two-panel figure:
  - Panel A: bar chart of reversal counts per signature at |z|≥1.5, |z|≥2, |z|≥3 (grouped bars, reversal vs concordant).
  - Panel B: scatter plot of spaceflight mean_lfc (x) vs compound z-score (y) for GSE128097_1, colored by reversal status (reversed=blue, concordant=red, non-significant=grey), with quadrant labels. This visualizes the "opposite forcing" at the gene level.
- **F12_reversed_genes_heatmap.svg/png** — heatmap of the top 50 reversed genes (by reversal_score) × 3 signatures, showing z-scores. Rows annotated with spaceflight direction. ComplexHeatmap.

### Manuscript: `manuscript/manuscript.md` (new)

Create a new manuscript.md for the astronaut-opposite-forcing repo in npj Microgravity format (matching the reference Python repo's style). Structure:
- Title, author (Richard Barker, admin@cosecloud.com), affiliation TODO
- Abstract (~250 words) covering: pipeline purpose, 1,610 oncogenic biomarkers, 21 reversal candidates at tau=-100, gene-level reversal by flavonoids, disease enrichment
- Introduction (spaceflight oncogenic risk, opposite-forcing concept, LINCS L1000 rationale)
- Results sections:
  1. Spaceflight oncogenic biomarker signature (1,610 genes, metafor REML, validation vs S2)
  2. Opposite-forcing drug screen (21 tau=-100 candidates, top hits)
  3. **Gene-level reversal by flavonoids** (NEW — the T13 results: N genes reversed at |z|≥1.5 by apigenin+luteolin and licochalcone A, top reversed genes, quadrant analysis, F11/F12)
  4. Disease enrichment (PrimeKG + Open Targets)
  5. Nutraceutical flags
- Discussion (interpretation, limitations, vegan recipe future work)
- Methods (concise: OSDR query, DESeq2, metafor, LINCS tau-analog, gene-level reversal definition, PrimeKG, Open Targets)
- References (key: Lamb 2006, Subramanian 2017, NASA OSDR, LINCS L1000, PrimeKG, Open Targets)

### README update

Add a new section "## Gene-level reversal by flavonoids" after the "Key results" section, summarizing:
- N genes directly reversed at |z|≥1.5 by each signature
- The union count (T13c)
- Top reversed genes (table)
- Reference to F11/F12 and T13/T13b/T13c

## Compute/resource estimate

- LINCS matrix load: 1.6 GB TSV → ~2 GB RAM (data.table), ~30s on heavy machine (32 GB RAM, already provisioned).
- Gene-level computation: vectorized sign/magnitude ops over 1,540 × 3 = 4,620 cells — trivial (<1s).
- Figure generation: ggplot2 + ComplexHeatmap, <30s.
- Manuscript writing: no compute.
- **Execution target**: heavy machine (machine_id="heavy", 32 GB RAM) — already provisioned, sufficient for the 1.6 GB load.

## Acceptance criteria

1. `scripts/07_gene_level_reversal.R` runs clean and produces T13, T13b, T13c.
2. F11 and F12 generated as SVG + PNG with Liberation Sans font.
3. `manuscript/manuscript.md` created with a Results section incorporating the gene-level reversal findings (with actual numbers from T13).
4. README.md updated with the gene-level reversal section.
5. All new deliverables copied to `/mnt/results/astronaut-opposite-forcing/`.
6. All 14+1=15 R files syntax-check OK.

## Assumptions

- The 70 T4 genes not found in LINCS (1,610 → 1,540 overlap) are excluded with a note in the manuscript methods.
- GSE128097 contains apigenin AND luteolin as co-treatment — gene-level effects cannot be attributed to one compound alone; the manuscript will state this clearly.
- GSE128097_2 is included per user request even though its WCS-level tau is +8.7 (not a reversal); the per-gene analysis may still show individual reversed genes.
- The reversal_score is a descriptive ranking metric, not a statistical test (LINCS z-scores don't carry per-gene p-values).

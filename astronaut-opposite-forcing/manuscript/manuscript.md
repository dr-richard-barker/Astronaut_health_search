# Spaceflight oncogenic biomarkers and their transcriptomic reversal by food-derived flavonoids: an R pipeline for opposite-forcing drug discovery

**Richard Barker¹\*** *(co-authors TBD)*

¹ *[Affiliation — TODO]*
\* Correspondence: admin@cosecloud.com

*Prepared for submission to* npj Microgravity. *Independent research; not affiliated with or endorsed by NASA or JAXA.*

---

## Abstract

Spaceflight exposes tissues to microgravity and ionising radiation, stressors independently linked to proliferation, DNA-damage, senescence and immune-surveillance changes that overlap with early cancer biology. Whether a reproducible oncogenic biomarker signature can be re-derived across human and rodent spaceflight transcriptomes, and whether existing drug perturbation signatures contain compounds that reverse it ("opposite forcing"), remains unresolved. We built a from-scratch R pipeline that (1) re-discovers the spaceflight oncogenic biomarker signature across 23 NASA OSDR RNA-seq studies (11 human, 12 rodent) using strict differential expression thresholds (padj < 0.05, |log2FC| > 1) and a random-effects meta-analysis (metafor REML), (2) screens 420 LINCS L1000 drug perturbation signatures for transcriptomic reversal via a tau-analog connectivity score, (3) annotates top reversal hits with drug-target-indication context from the Broad Drug Repurposing Hub and PrimeKG, (4) maps the biomarker signature to disease phenotypes via PrimeKG and Open Targets, and (5) flags nutraceutical candidates among the top reversal compounds. The meta-analysis scored 40,211 genes (3,668 significant), whose intersection with an 11,241-gene oncogene union yielded **1,610 oncogenic biomarker genes**. Validation against an independent Python cross-study meta-signature (S2) showed Spearman ρ = 0.60 and 76% sign agreement. Twenty-one compounds exhibited maximal reversal (tau = −100), led by the food-derived flavonoids **apigenin, luteolin, and licochalcone A**. Gene-level analysis revealed that **242 of the 1,540 assessable oncogenic biomarker genes (15.7%)** were directly reversed (opposite direction, |z| ≥ 1.5) by at least one flavonoid signature, with 20 genes reversed by all three. The most strongly reversed genes were inflammatory chemokines and SASP factors (CXCL8, MMP3, CSF2, CSF3, CXCL3) down-regulated by spaceflight and up-regulated by flavonoids. Disease enrichment implicated sarcomatoid carcinoma, coronary artery disease, prostate cancer, hereditary breast/ovarian cancer, and fatty liver disease. These results establish a reproducible framework for spaceflight oncogenic biomarker discovery and nominate dietary flavonoids as candidate countermeasure agents for prospective validation.

---

## Introduction

Human spaceflight combines chronic microgravity with exposure to galactic cosmic radiation. Both stressors have been individually associated with hallmarks of early oncogenesis: genomic instability and DNA-damage responses, altered proliferation and differentiation, cellular senescence and its secretory phenotype (SASP), and impaired immune surveillance — processes central to the hallmarks of cancer and to galactic-cosmic-ray carcinogenesis risk. As mission durations lengthen and commercial spaceflight broadens the flying population, understanding whether and how spaceflight remodels gene expression toward oncogenic states — and whether early biomarkers of that remodelling exist — has become a priority for crew health.

Existing space-omics analyses typically process one dataset at a time through a standard pipeline (differential expression, co-expression network analysis, and Gene Ontology / pathway enrichment). This describes *what changes* in a given experiment but does not test whether changes generalise across studies and species, whether they converge on oncogenic programs, or whether existing pharmacological or nutritional compounds can reverse them. The concept of "opposite forcing" — identifying compounds whose transcriptomic perturbation signature is anti-correlated with a disease signature — was introduced by the Connectivity Map and later scaled by the LINCS L1000 program, but has not been applied to spaceflight oncogenic biomarkers.

Here we address these gaps with a from-scratch R pipeline that integrates human and rodent spaceflight transcriptomes, re-derives an oncogenic biomarker signature via random-effects meta-analysis, screens LINCS L1000 drug perturbation signatures for transcriptomic reversal, and flags nutraceutical candidates among the top reversal hits. Our central questions are: (i) does a reproducible oncogenic biomarker signature emerge across human and rodent spaceflight studies? (ii) which existing drug perturbation signatures reverse it? (iii) are food-derived flavonoids among the reversal candidates, and which biomarker genes do they directly counter-regulate?

---

## Results

### A cross-species spaceflight oncogenic biomarker signature

We queried the NASA OSDR programmatically and curated a panel of 24 RNA-seq studies (11 human, 13 rodent) spanning spaceflight, simulated microgravity, hindlimb suspension, and ionising radiation across diverse tissues (cardiomyocytes, endothelial and vascular smooth muscle cells, skeletal muscle, neural organoids, spleen, thymus, liver, skin, eye, heart, bone, blood; **Table 1**, **Fig. 1**). Per-study differential expression (DESeq2; padj < 0.05, |log2FC| > 1) yielded 410,121 gene-level results across 23 studies with testable contrasts (**Table 2**). Rodent gene identifiers were mapped to human orthologs via biomaRt before meta-analysis.

A random-effects meta-analysis (metafor REML) across all studies scored 40,211 genes, of which 3,668 were significant at padj < 0.05 (**Table 3**, **Fig. 3**). The signature was directionally heterogeneous: 502 genes were up-regulated and 1,108 down-regulated in spaceflight, consistent with tissue- and stressor-specific engagement of oncogenic programs. To focus on genes with oncogenic relevance, we intersected the meta-significant genes with an 11,241-gene oncogene union (MSigDB C6 oncogenic signatures: 10,927; KEGG oncogenesis: 931; curated oncogene/TSG/SASP: 195; COSMIC CGC and OncoKB were attempted but gated — see Methods), yielding **1,610 oncogenic biomarker genes** (**Table 4**, **Fig. 4**).

Validation against an independent Python cross-study meta-signature (S2) from our companion repository showed Spearman ρ = 0.60, sign agreement 76%, Jaccard top-100 = 0.19, and Jaccard top-500 = 0.28 (**Table 5**, **Fig. 5**), confirming that the R reimplementation recovers the core spaceflight oncogenic signal.

### Opposite-forcing drug screen identifies 21 reversal candidates

We screened 420 LINCS L1000 GEO drug perturbation signatures (z-score profiles over 35,238 genes) for transcriptomic reversal of the 1,610 oncogenic biomarkers via a tau-analog weighted connectivity score (WCS; Lamb 2006 / Subramanian 2017 formulation, tau-normalised across all drug signatures; **Table 6**). Twenty-one compounds exhibited maximal reversal (tau = −100; **Table 7**, **Fig. 6**), including the food-derived flavonoids **apigenin and luteolin** (GSE128097), **licochalcone A** (GSE137934, a licorice-derived flavonoid), the MDM2 inhibitor Nutlin-3, the topoisomerase inhibitor doxorubicin (Phase 3), PI3K-mTOR inhibitors, a combined MEKi+WNT treatment, the DNA methylation inhibitor 5'-Aza, methotrexate, the ALK5/TGFβ inhibitor SB431542, and ketamine.

Tissue-stratified reversal scoring across 19 tissues (4,620 tissue-drug combinations; **Table 8**) confirmed that reversal strength varies by tissue, consistent with the tissue-specific direction of oncogenic program engagement. Drug-target-indication context from the Broad Drug Repurposing Hub and PrimeKG (**Table 9**, **Fig. 9**) annotated 7 of the top 50 hits with clinical phase, mechanism of action, and target information.

### Gene-level reversal by food-derived flavonoids

Because apigenin, luteolin, and licochalcone A were among the strongest reversal candidates and are dietary compounds relevant to the "vegan astronaut" nutrition goal, we examined which of the 1,610 oncogenic biomarker genes they directly reverse at the gene level. We extracted the LINCS z-score vectors for the three flavonoid signatures — GSE128097_1 (apigenin + luteolin co-treatment, tau = −100), GSE128097_2 (apigenin + luteolin second sub-signature, tau = +8.7), and GSE137934_1 (licochalcone A, tau = −100) — and compared each gene's compound z-score direction against its spaceflight meta-log2FC direction.

Of the 1,610 oncogenic biomarker genes, 1,540 (95.7%) were present in the LINCS gene space; the remaining 70 (non-protein-coding or unmapped identifiers) were excluded. A gene was defined as "directly reversed" when its spaceflight direction (sign of mean log2FC) was opposite to its compound z-score direction, with |z| ≥ 1.5 (a moderate LINCS effect threshold; sensitivity analyses at |z| ≥ 2 and |z| ≥ 3 are reported).

**242 unique oncogenic biomarker genes (15.7% of the 1,540 assessable) were directly reversed by at least one flavonoid signature at |z| ≥ 1.5** (**Table 13c**, **Fig. 11**). Of these, 20 were reversed by all three signatures, 87 by two, and 135 by one. Per-signature reversal counts were: apigenin+luteolin_1: 134 reversed (ratio 0.55); apigenin+luteolin_2: 136 reversed (ratio 0.52); licochalcone A: 99 reversed (ratio 0.49) (**Table 13b**). At the stricter |z| ≥ 2 threshold, 148 genes were reversed; at |z| ≥ 3, 67 genes. The reversal ratio (reversed / total with |z| above threshold) increased with threshold (0.49–0.55 at |z| ≥ 1.5 to 0.60–0.70 at |z| ≥ 3), indicating that stronger compound effects are enriched for reversal over concordance.

The most strongly reversed genes (ranked by reversal score = −sign(mean_lfc) × z × |mean_lfc|; **Table 13c**, **Fig. 12**) were dominated by inflammatory chemokines and senescence-associated secretory phenotype (SASP) factors that spaceflight down-regulates and flavonoids up-regulate: **CXCL8** (interleukin-8, z = 8.49, reversed by 2 signatures), **MMP3** (matrix metalloproteinase-3, z = 4.92, reversed by all 3), **CSF2** (GM-CSF, z = 6.13, reversed by all 3), **CSF3** (G-CSF, z = 3.82), **CXCL3** (z = 4.88, reversed by all 3), **AREG** (amphiregulin, z = 6.80, reversed by all 3), and **HAS1** (hyaluronan synthase-1, z = 3.98). Genes up-regulated by spaceflight and down-regulated by flavonoids included **DMRT1** (z = −4.03, reversed by 2) and **REN** (renin, z = −6.42, reversed by licochalcone A). The class distribution of reversed genes was 224 C6 oncogenic, 11 SASP, 5 oncogene, 1 tumour suppressor, and 1 KEGG oncogenesis.

Notably, GSE128097_2 — whose signature-level tau was +8.7 (not a reversal at the WCS level) — still reversed 136 individual genes at |z| ≥ 1.5, demonstrating that a signature that does not globally reverse a disease profile can still counter-regulate a substantial subset of individual biomarker genes. This underscores the value of gene-level analysis beyond signature-level connectivity scoring.

### Disease enrichment of the oncogenic biomarker signature

To validate that the 1,610 oncogenic biomarker genes are biologically linked to disease phenotypes rather than noise, we performed hypergeometric enrichment against PrimeKG (160,822 disease-gene edges) and gene-disease association scoring via the Open Targets GraphQL API.

PrimeKG yielded 389 significant diseases (padj < 0.05; **Table 10**), led by sarcomatoid carcinoma (padj = 3.3 × 10⁻¹², 47 overlapping genes), undifferentiated carcinoma (padj = 3.3 × 10⁻¹²), coronary artery disease (padj = 3.4 × 10⁻¹¹, 58 genes), carcinoma (padj = 4.5 × 10⁻¹¹), and prostate cancer. Open Targets returned 4,891 gene-disease associations for 100 queried genes across 1,689 diseases (**Table 11**), led by neoplasm, cancer, hepatocellular carcinoma, breast carcinoma, breast cancer, NSCLC, Alzheimer's disease, glioblastoma, melanoma, and colorectal carcinoma. A cross-ranked combined table (**Table 10b**, **Fig. 8**) integrated both sources, with breast carcinoma, breast cancer, hepatocellular carcinoma, and prostate carcinoma ranking highest by combined score.

### Nutraceutical candidates among reversal hits

Among the top 50 reversal hits, we flagged nutraceuticals using a curated 176-compound lexicon (17 classes: flavonoids, stilbenes, curcuminoids, carotenoids, vitamins, glucosinolates, omega-3 lipids, and others) cross-matched against LINCS GMT drug names and the Broad Drug Repurposing Hub (**Table 12**, **Fig. 10**). Five unique nutraceuticals were flagged: **luteolin** (tau = −100, 3 evidence streams: signature + LINCS GMT + Broad), **apigenin** (tau = −100, 2 streams), **licochalcone A** (tau = −100, 2 streams), licochalcone, and isoginkgetin (tau = +8.7, not a reversal — same direction as spaceflight). A full scan of all 420 scored signatures identified 37 nutraceutical-bearing signatures (**Table 12c**) as input for a deferred vegan-recipe optimisation step.

---

## Discussion

Across a cross-species panel of 23 spaceflight RNA-seq studies, we re-derived a 1,610-gene oncogenic biomarker signature that concords with an independent Python meta-analysis (Spearman ρ = 0.60) and enriches for carcinoma, coronary artery disease, and hereditary cancer syndromes. Screening LINCS L1000 drug perturbation signatures identified 21 compounds with maximal transcriptomic reversal (tau = −100), led by food-derived flavonoids. Gene-level analysis revealed that 242 oncogenic biomarker genes (15.7% of assessable) are directly reversed by at least one flavonoid signature, with inflammatory chemokines and SASP factors (CXCL8, MMP3, CSF2, CSF3, CXCL3) most strongly counter-regulated.

The convergence of flavonoid reversal on inflammatory and SASP chemokines is biologically coherent. Spaceflight has been shown to dysregulate immune signalling and promote a senescence-associated secretory phenotype; flavonoids such as apigenin and luteolin are known modulators of NF-κB and inflammatory cytokine production. The finding that these dietary compounds reverse the spaceflight-driven suppression of CXCL8, CSF2, and MMP3 at the transcriptomic level nominates them as candidate nutritional countermeasures for prospective validation.

The gene-level analysis also revealed that signature-level and gene-level reversal can diverge: GSE128097_2 was not a reversal candidate at the WCS/tau level (tau = +8.7) yet reversed 136 individual biomarker genes. This has a practical implication: signature-level connectivity scoring may miss compounds that reverse a biologically important subset of a disease signature even when the global profile is not anti-correlated. Gene-level reversal analysis should complement, not replace, signature-level screening.

**Limitations.** (1) The LINCS tau-analog is computed from GEO-derived z-score signatures, not the canonical CLUE.io Touchstone GCTx with pre-computed tau; WCS values are near-zero for most drugs, so tau-normalisation clips at ±100. The ranking is correct but absolute tau values are compressed. (2) GSE128097 is a co-treatment of apigenin and luteolin; gene-level effects cannot be attributed to one compound alone. (3) LINCS z-scores lack per-gene p-values, so the |z| ≥ 1.5 threshold is an effect-size cutoff, not a statistical test; the reversal score is descriptive. (4) COSMIC CGC and OncoKB were gated; the oncogene union relied on MSigDB C6, KEGG, and curated lists. (5) The 70 T4 genes absent from LINCS were excluded from gene-level analysis. (6) These are transcriptomic associations, not evidence of efficacy in vivo; the flavonoid reversal candidates require prospective validation in spaceflight-relevant model systems. (7) Vegan recipe optimisation is deferred; this pipeline provides the nutraceutical-gene input (T12c, T13c) for that future step.

---

## Methods

### Data acquisition. We queried the NASA OSDR API for human and rodent RNA-seq studies with processed count matrices and a testable spaceflight/microgravity/radiation contrast, yielding 24 studies (11 human, 13 rodent; **Table 1**). JAXA6 astronaut whole-blood data were obtained from the Astronaut_health_search GitHub repository. Raw STAR unnormalised count matrices were downloaded per study.

### Differential expression. Per-study DESeq2 analysis was run with the contrast spaceflight/microgravity/radiation (treated) vs ground/vivarium/sham control. Genes with padj < 0.05 and |log2FC| > 1 were called significant. Rodent gene identifiers were mapped to human orthologs via biomaRt (getLDS with mirror retry logic; MGI symbol fallback).

### Meta-analysis. A random-effects meta-analysis was performed per gene using metafor::rma (REML, z-test) over the per-study log2FC and standard error values. The meta-signature comprised 40,211 genes with mean log2FC, standard error, z, p, padj, direction concordance, and species count (**Table 3**).

### Oncogene union. We assembled an 11,241-gene oncogene union from MSigDB C6 oncogenic signatures (10,927 genes via msigdbr), KEGG oncogenesis pathways (931 via KEGGREST), and a curated oncogene/tumour-suppressor/SASP list (195). COSMIC CGC and OncoKB were attempted but gated (all COSMIC URLs failed; OncoKB returned HTTP 401). The intersection of meta-significant genes with the oncogene union yielded 1,610 oncogenic biomarker genes.

### LINCS tau-analog reversal screen. We loaded the LINCS L1000 GEO z-score signature matrix (4,269 signatures × 35,238 genes; human_geo_sigs.tsv) and filtered to 420 drug-related signatures. For each drug signature, we computed a weighted connectivity score (WCS) against the disease UP/DN gene sets using a weighted Kolmogorov-Smirnov enrichment (CMap a-score formulation), then tau-normalised across all drug signatures (clipped to ±100). A strongly negative tau indicates reversal.

### Gene-level reversal. For the three flavonoid signatures (GSE128097_1, GSE128097_2, GSE137934_1), we extracted the z-score vector over all genes and intersected with the 1,610 oncogenic biomarkers (1,540 overlap). A gene was "directly reversed" when sign(spaceflight mean log2FC) × sign(compound z) < 0 and |z| ≥ 1.5. Reversal score = −sign(mean_lfc) × z × |mean_lfc|. Sensitivity thresholds |z| ≥ 2 and |z| ≥ 3 are reported.

### Disease enrichment. PrimeKG disease-gene edges (160,822) were tested for hypergeometric enrichment of the 1,610 oncogenic biomarkers (universe = all genes in PrimeKG disease-gene edges; BH-adjusted p-values). Open Targets GraphQL API (v4) was queried for 100 genes, returning 4,891 gene-disease associations across 1,689 diseases. A combined ranking integrated both sources.

### Drug context. The Broad Drug Repurposing Hub (6,798 drugs; clinical phase, MoA, target, indication) and PrimeKG drug-target (51,306 edges) and contraindication (61,350 edges) networks were used to annotate the top 50 reversal hits.

### Nutraceutical flagging. A curated 176-compound lexicon (17 classes) was cross-matched against the top 50 reversal hits via word-boundary regex on signature names and extracted drug names, with independent evidence from LINCS GMT drug names and the Broad Drug Repurposing Hub.

### Software and reproducibility. All analysis was performed in R 4.4.3 with DESeq2, limma, edgeR, biomaRt, metafor, clusterProfiler, msigdbr, fgsea, ComplexHeatmap, KEGGREST, data.table, ggplot2, and patchwork. Package versions are pinned in renv.lock. The pipeline is executed via run_all.R with --from/--only flags for resumability. All figures are saved as SVG (editable) and PNG.

---

## Data and code availability

- Code: https://github.com/dr-richard-barker/astronaut-opposite-forcing
- Zenodo deposit: TODO (DOI upon publication)
- OSDR data: https://osdr.nasa.gov/bio/repo/ (public)
- LINCS L1000 signatures, Broad Drug Repurposing Hub, PrimeKG: accessed via the Biomni datalake
- Companion Python repository: https://github.com/dr-richard-barker/astronaut-oncogene-biomarkers

---

## References

1. Hanahan D, Weinberg RA. Hallmarks of cancer: the next generation. Cell. 2011;144(5):646-674.
2. Durante M, Cucinotta FA. Heavy ion carcinogenesis and human space exploration. Nat Rev Cancer. 2008;8(6):465-472.
3. Lamb J et al. The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. Science. 2006;313(5795):1929-1935.
4. Subramanian A et al. A next generation connectivity map: L1000 platform and the first 1,000,000 profiles. Cell. 2017;171(6):1437-1452.
5. Beheshti A et al. NASA Open Science Data Repository: a comprehensive biomedical data repository for spaceflight research. Cell. 2023;186(18):3846-3855.
6. Duan Q et al. L1000 CMap2: a next-generation connectivity map. Sci Data. 2023.
7. Chandrasekharan Nair S et al. PrimeKG: a precision medicine knowledge graph. Sci Data. 2023;10:189.
8. Ochoa D et al. Open Targets Platform: supporting systematic drug-target identification and prioritisation. Nucleic Acids Res. 2023;51(D1):D1354-D1361.
9. Liberzon A et al. The Molecular Signatures Database Hallmark Gene Set Collection. Cell Syst. 2015;1(6):417-425.
10. Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw. 2010;36(3):1-48.

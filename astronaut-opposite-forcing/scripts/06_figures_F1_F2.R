#!/usr/bin/env Rscript
# ==============================================================================
# 06_figures_F1_F2.R
# astronaut-opposite-forcing pipeline
#
# Generates the two remaining overview/QC figures:
#   F1_dataset_panel_overview  -- study x tissue x species x factor panel
#   F2_qc_pca                  -- per-study vst PCA colored by condition
#
# Reads:
#   results/tables/T1_dataset_panel.csv
#   data/processed/sample_metadata.csv
#   data/raw/<OSD-*>/GLDS-*_rna_seq_STAR_Unnormalized_Counts_GLbulkRNAseq.csv
#
# Writes:
#   results/figures/F1_dataset_panel_overview.svg / .png
#   results/figures/F2_qc_pca.svg / .png
#
# Richard Barker -- 2026
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(DESeq2)
})

repo <- "/mnt/shared-workspace/astronaut-opposite-forcing"
fig_dir <- file.path(repo, "results", "figures")
raw_dir <- file.path(repo, "data", "raw")
proc_dir <- file.path(repo, "data", "processed")
t1_path <- file.path(repo, "results", "tables", "T1_dataset_panel.csv")
meta_path <- file.path(proc_dir, "sample_metadata.csv")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Phylo color palette
phylo <- c("#000000", "#ECE9E2", "#FAF9F3", "#E9ED4C", "#FF9400",
           "#75A025", "#FD9BED", "#0279EE")

# ------------------------------------------------------------------------------
# F1 -- Dataset panel overview
# ------------------------------------------------------------------------------
t1 <- fread(t1_path)
meta <- fread(meta_path)

# Compute per-study sample counts from metadata (T1 has NA for n_control/n_treated)
study_counts <- meta[, .(n_samples = .N,
                         n_control = sum(condition == "control", na.rm = TRUE),
                         n_treated = sum(condition == "treated", na.rm = TRUE)),
                     by = .(study, species, tissue, factor, system)]
setnames(study_counts, "study", "accession")
t1 <- merge(t1[, .(accession, tissue, system, factor, species, note)],
            study_counts, by = c("accession", "tissue", "system", "factor",
                                 "species"), all.x = TRUE)
# Order studies by species then factor then tissue
t1[, species_factor := factor(species, levels = c("human", "mouse"))]
t1[, factor_f := factor(factor, levels = c("spaceflight", "microgravity",
                                           "simulated_microgravity",
                                           "microgravity_radiation",
                                           "space_radiation", "ionizing_radiation",
                                           "hindlimb_suspension"))]
t1[, accession := factor(accession, levels = unique(accession[order(species_factor, factor_f, tissue)]))]

# Tile plot: study x tissue, fill = factor, size = n_samples
p1 <- ggplot(t1, aes(x = tissue, y = accession)) +
  geom_tile(aes(fill = factor_f), color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(is.na(n_samples), "", n_samples)),
            size = 3.2, color = "black", fontface = "bold") +
  scale_fill_manual(values = c(
    "spaceflight" = "#0279EE",
    "microgravity" = "#75A025",
    "simulated_microgravity" = "#E9ED4C",
    "microgravity_radiation" = "#FF9400",
    "space_radiation" = "#FD9BED",
    "ionizing_radiation" = "#000000",
    "hindlimb_suspension" = "#FAF9F3"
  ), name = "Factor") +
  facet_grid(species ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Spaceflight oncogenic biomarker dataset panel",
       subtitle = sprintf("%d studies (%d human, %d rodent), %d samples total",
                          nrow(t1), sum(t1$species == "human"),
                          sum(t1$species == "mouse"), sum(t1$n_samples, na.rm = TRUE)),
       x = "Tissue / model", y = "OSDR accession") +
  theme_minimal(base_family = "Liberation Sans") +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        strip.background = element_rect(fill = "grey92", color = NA),
        legend.position = "right")

ggsave(file.path(fig_dir, "F1_dataset_panel_overview.svg"), p1,
       width = 11, height = 8, dpi = 150)
ggsave(file.path(fig_dir, "F1_dataset_panel_overview.png"), p1,
       width = 11, height = 8, dpi = 150)
cat("[F1] wrote F1_dataset_panel_overview.svg/png\n")

# ------------------------------------------------------------------------------
# F2 -- QC PCA from per-study vst-transformed counts
# ------------------------------------------------------------------------------
# For each study with a STAR unnormalized count matrix, run DESeq2 vst on the
# raw counts (limited to the samples in metadata), compute PCA, and collect
# PC1/PC2 for a combined panel. We do NOT merge across studies (different
# gene spaces, species, platforms) -- each study gets its own PCA subplot.
studies <- unique(meta$study)
pca_list <- list()
n_processed <- 0
for (acc in studies) {
  raw_subdir <- file.path(raw_dir, acc)
  if (!dir.exists(raw_subdir)) next
  # Find STAR unnormalized counts file
  star_files <- list.files(raw_subdir,
                           pattern = "STAR_Unnormalized_Counts.*\\.csv$",
                           full.names = TRUE)
  if (length(star_files) == 0) {
    # Fall back to RSEM unnormalized
    star_files <- list.files(raw_subdir,
                             pattern = "RSEM_Unnormalized_Counts.*\\.csv$",
                             full.names = TRUE)
  }
  if (length(star_files) == 0) next
  counts_file <- star_files[1]
  # Read counts; first column is gene ID
  counts <- tryCatch(fread(counts_file), error = function(e) NULL)
  if (is.null(counts)) next
  # First column is gene ID; rename
  gene_col <- colnames(counts)[1]
  setnames(counts, gene_col, "gene_id")
  # Match samples to metadata
  meta_sub <- meta[study == acc & !is.na(condition)]
  if (nrow(meta_sub) < 4) next
  # Find count columns matching sample names (GSM IDs)
  sample_cols <- intersect(colnames(counts), meta_sub$sample_name)
  if (length(sample_cols) < 4) next
  # Build count matrix (genes x samples)
  mat <- as.matrix(counts[, ..sample_cols])
  rownames(mat) <- counts$gene_id
  # Filter very low-count genes
  keep <- rowSums(mat) >= 10
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) < 100) next
  # Build DESeq2 object
  meta_sub <- meta_sub[sample_name %in% sample_cols]
  meta_sub <- meta_sub[match(sample_cols, sample_name)]
  dds <- tryCatch({
    DESeqDataSetFromMatrix(countData = mat,
                           colData = data.frame(condition = meta_sub$condition,
                                                 tissue = meta_sub$tissue,
                                                 stringsAsFactors = FALSE),
                           design = ~ condition)
  }, error = function(e) NULL)
  if (is.null(dds)) next
  # vst (variance-stabilizing transformation) -- blind = TRUE for QC
  vsd <- tryCatch(vst(dds, blind = TRUE), error = function(e) NULL)
  if (is.null(vsd)) next
  # PCA
  pca <- tryCatch({
    pc <- prcomp(t(assay(vsd)), rank. = 2, scale. = FALSE)
    vars <- pc$sdev^2 / sum(pc$sdev^2)
    data.table(
      study = acc,
      species = meta_sub$species[1],
      tissue = meta_sub$tissue[1],
      factor = meta_sub$factor[1],
      sample = rownames(pc$x),
      PC1 = pc$x[, 1],
      PC2 = pc$x[, 2],
      var1 = vars[1],
      var2 = vars[2],
      condition = meta_sub$condition
    )
  }, error = function(e) NULL)
  if (!is.null(pca)) {
    pca_list[[acc]] <- pca
    n_processed <- n_processed + 1
    cat(sprintf("[F2] %s: %d samples, %d genes -> PCA OK\n",
                acc, nrow(pca), nrow(mat)))
  }
}

if (length(pca_list) > 0) {
  pca_dt <- rbindlist(pca_list)
  # Faceted PCA: one panel per study
  pca_dt[, study_label := sprintf("%s (%s)", study, tissue)]
  # Order panels by species then factor
  study_order <- unique(pca_dt[, .(study, species, factor, tissue)])
  study_order[, species_f := factor(species, levels = c("human", "mouse"))]
  study_order[, factor_f := factor(factor, levels = c("spaceflight", "microgravity",
                       "simulated_microgravity", "microgravity_radiation",
                       "space_radiation", "ionizing_radiation", "hindlimb_suspension"))]
  study_order <- study_order[order(species_f, factor_f, tissue)]
  pca_dt[, study_label := factor(study_label, levels = sprintf("%s (%s)",
                                                               study_order$study,
                                                               study_order$tissue))]

  p2 <- ggplot(pca_dt, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 1.8, alpha = 0.85) +
    stat_ellipse(level = 0.68, linetype = "dashed", show.legend = FALSE) +
    facet_wrap(~ study_label, scales = "free", ncol = 5) +
    scale_color_manual(values = c("control" = "#0279EE",
                                  "treated" = "#FF9400"),
                       name = "Condition") +
    labs(title = "Per-study QC PCA (DESeq2 vst, blind = TRUE)",
         subtitle = sprintf("%d studies with sufficient samples; dashed = 1 SD ellipse",
                            n_processed),
         x = "PC1", y = "PC2") +
    theme_minimal(base_family = "Liberation Sans") +
    theme(plot.title = element_text(face = "bold"),
          strip.text = element_text(size = 7.5),
          strip.background = element_rect(fill = "grey92", color = NA),
          panel.grid.minor = element_blank(),
          legend.position = "bottom")

  ggsave(file.path(fig_dir, "F2_qc_pca.svg"), p2,
         width = 14, height = 10, dpi = 150)
  ggsave(file.path(fig_dir, "F2_qc_pca.png"), p2,
         width = 14, height = 10, dpi = 150)
  cat(sprintf("[F2] wrote F2_qc_pca.svg/png (%d studies, %d samples)\n",
              n_processed, nrow(pca_dt)))
} else {
  cat("[F2] WARNING: no studies produced PCA; writing placeholder\n")
  p2 <- ggplot() + annotate("text", x = 0.5, y = 0.5,
                            label = "No studies with sufficient samples for PCA",
                            size = 6) +
    theme_void() + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(fig_dir, "F2_qc_pca.svg"), p2, width = 10, height = 6, dpi = 150)
  ggsave(file.path(fig_dir, "F2_qc_pca.png"), p2, width = 10, height = 6, dpi = 150)
}

cat("\n[06_figures_F1_F2.R] DONE\n")

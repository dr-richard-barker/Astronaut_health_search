#!/usr/bin/env Rscript
# ==============================================================================
# 07_gene_level_reversal.R
# astronaut-opposite-forcing pipeline
#
# Identify which of the 1,610 oncogenic biomarker genes (T4) are directly
# reversed at the gene level by the three flavonoid LINCS signatures:
#   GSE128097_1 -- apigenin + luteolin co-treatment (tau = -100, reversal)
#   GSE128097_2 -- apigenin + luteolin second sub-signature (tau = +8.7)
#   GSE137934_1 -- licochalcone A (tau = -100, reversal)
#
# A gene is "directly reversed" when:
#   (a) present in T4 (1,610 oncogenic biomarkers)
#   (b) present in the LINCS signature (non-zero z-score)
#   (c) sign(spaceflight mean_lfc) * sign(compound z) < 0  (opposite direction)
#   (d) |compound z| >= 1.5  (moderate LINCS effect threshold)
#
# Outputs:
#   results/tables/T13_gene_level_reversal.csv       -- long: gene x signature
#   results/tables/T13b_gene_level_reversal_summary.csv -- per-signature counts
#   results/tables/T13c_reversed_genes_union.csv     -- union of reversed genes
#   results/figures/F11_gene_level_reversal.svg/png  -- bar + scatter
#   results/figures/F12_reversed_genes_heatmap.svg/png -- top-50 heatmap
#
# Richard Barker -- 2026
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

repo <- "/mnt/shared-workspace/astronaut-opposite-forcing"
fig_dir <- file.path(repo, "results", "figures")
tab_dir <- file.path(repo, "results", "tables")
t4_path <- file.path(tab_dir, "T4_oncogene_intersection.csv")
lincs_path <- "/mnt/datalake/LINCS1000/RNAseq_transcriptomics_signatures/human_geo_sigs.tsv"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# Phylo palette
phylo_blue <- "#0279EE"
phylo_orange <- "#FF9400"
phylo_grey <- "#B0B0B0"
phylo_green <- "#75A025"
phylo_pink <- "#FD9BED"

# ------------------------------------------------------------------------------
# 1. Load T4 oncogenic biomarkers
# ------------------------------------------------------------------------------
t4 <- fread(t4_path)
cat(sprintf("[07] T4 oncogenic biomarkers: %d genes\n", nrow(t4)))
cat(sprintf("     direction: up=%d, down=%d\n",
            sum(t4$mean_lfc > 0), sum(t4$mean_lfc < 0)))

# ------------------------------------------------------------------------------
# 2. Load LINCS z-score matrix (1.6 GB, ~30s on heavy machine)
# ------------------------------------------------------------------------------
cat("[07] loading LINCS z-score matrix (1.6 GB)...\n")
mat <- fread(lincs_path, header = TRUE, sep = "\t", data.table = TRUE)
sig_names <- mat[[1]]
gene_names <- colnames(mat)[-1]
cat(sprintf("[07] matrix: %d signatures x %d genes\n",
            length(sig_names), length(gene_names)))

# ------------------------------------------------------------------------------
# 3. Identify the 3 flavonoid signatures
# ------------------------------------------------------------------------------
flav_idx <- grep("GSE128097|GSE137934", sig_names)
flav_sigs <- sig_names[flav_idx]
cat(sprintf("[07] flavonoid signatures found: %d\n", length(flav_sigs)))
for (s in flav_sigs) cat(sprintf("     - %s\n", substr(s, 1, 70)))

# Map to clean compound labels
compound_label <- function(sig_name) {
  if (grepl("GSE128097_1", sig_name)) "apigenin+luteolin_1"
  else if (grepl("GSE128097_2", sig_name)) "apigenin+luteolin_2"
  else if (grepl("GSE137934", sig_name)) "licochalcone_A"
  else "unknown"
}

# ------------------------------------------------------------------------------
# 4. Compute gene-level reversal for each signature
# ------------------------------------------------------------------------------
# T4 genes present in LINCS
t4_genes <- t4$gene
t4_in_lincs <- intersect(t4_genes, gene_names)
cat(sprintf("[07] T4 genes in LINCS: %d / %d (%d missing)\n",
            length(t4_in_lincs), nrow(t4), nrow(t4) - length(t4_in_lincs)))
missing_genes <- setdiff(t4_genes, gene_names)

# T4 direction lookup
t4_dir <- setNames(sign(t4$mean_lfc), t4$gene)
t4_lfc <- setNames(t4$mean_lfc, t4$gene)
t4_padj <- setNames(t4$padj, t4$gene)
t4_class <- setNames(t4$class, t4$gene)
t4_source <- setNames(t4$source, t4$gene)
t4_nstudies <- setNames(t4$n_studies, t4$gene)
t4_concord <- setNames(t4$direction_concordance, t4$gene)

# Build long table: one row per (gene x signature)
all_rows <- list()
for (fi in flav_idx) {
  sig_name <- sig_names[fi]
  clabel <- compound_label(sig_name)
  z <- as.numeric(mat[fi, -1])
  names(z) <- gene_names
  # restrict to T4 genes present in LINCS
  z_t4 <- z[t4_in_lincs]
  # skip genes with zero or NA z (no measurable perturbation)
  has_z <- !is.na(z_t4) & z_t4 != 0
  genes_use <- t4_in_lincs[has_z]
  z_use <- z_t4[has_z]
  sf_dir <- t4_dir[genes_use]
  sf_lfc <- t4_lfc[genes_use]
  reversed <- (sf_dir * sign(z_use)) < 0
  reversal_score <- -sf_dir * z_use * abs(sf_lfc)
  rows <- data.table(
    gene = genes_use,
    class = t4_class[genes_use],
    source = t4_source[genes_use],
    mean_lfc = sf_lfc,
    padj = t4_padj[genes_use],
    n_studies = t4_nstudies[genes_use],
    direction_concordance = t4_concord[genes_use],
    signature = sig_name,
    compound_label = clabel,
    compound_z = z_use,
    spaceflight_dir = ifelse(sf_dir > 0, "up", "down"),
    compound_dir = ifelse(sign(z_use) > 0, "up", "down"),
    reversed = reversed,
    reversed_z1.5 = reversed & abs(z_use) >= 1.5,
    reversed_z2 = reversed & abs(z_use) >= 2,
    reversed_z3 = reversed & abs(z_use) >= 3,
    concordant_z1.5 = !reversed & abs(z_use) >= 1.5,
    reversal_score = reversal_score
  )
  all_rows[[length(all_rows) + 1]] <- rows
  cat(sprintf("[07] %s: %d T4 genes with non-zero z, reversed@|z|>=1.5: %d\n",
              clabel, nrow(rows), sum(rows$reversed_z1.5)))
}

t13 <- rbindlist(all_rows)
cat(sprintf("[07] T13 long table: %d rows (%d genes x %d signatures)\n",
            nrow(t13), uniqueN(t13$gene), length(flav_sigs)))

# Write T13
fwrite(t13, file.path(tab_dir, "T13_gene_level_reversal.csv"))
cat(sprintf("[07] wrote T13_gene_level_reversal.csv\n"))

# ------------------------------------------------------------------------------
# 5. T13b -- per-signature summary
# ------------------------------------------------------------------------------
t13b <- t13[, .(
  n_T4_overlap = .N,
  n_reversed_z1.5 = sum(reversed_z1.5),
  n_reversed_z2 = sum(reversed_z2),
  n_reversed_z3 = sum(reversed_z3),
  n_concordant_z1.5 = sum(concordant_z1.5),
  n_concordant_z2 = sum(!reversed & abs(compound_z) >= 2),
  n_concordant_z3 = sum(!reversed & abs(compound_z) >= 3),
  reversal_ratio_z1.5 = sum(reversed_z1.5) / sum(abs(compound_z) >= 1.5),
  reversal_ratio_z2 = sum(reversed_z2) / sum(abs(compound_z) >= 2),
  reversal_ratio_z3 = sum(reversed_z3) / sum(abs(compound_z) >= 3),
  mean_reversal_score = mean(reversal_score, na.rm = TRUE),
  median_reversal_score = median(reversal_score, na.rm = TRUE)
), by = compound_label]
setorder(t13b, -n_reversed_z1.5)
fwrite(t13b, file.path(tab_dir, "T13b_gene_level_reversal_summary.csv"))
cat(sprintf("[07] wrote T13b_gene_level_reversal_summary.csv\n"))
print(t13b)

# ------------------------------------------------------------------------------
# 6. T13c -- union of genes reversed at |z|>=1.5 by ANY signature
# ------------------------------------------------------------------------------
rev_genes <- t13[reversed_z1.5 == TRUE]
# Per gene: which signatures reversed it, max |z|, max reversal_score
t13c <- rev_genes[, .(
  n_signatures_reversed = .N,
  signatures_reversed = paste(compound_label, collapse = ";"),
  spaceflight_dir = unique(spaceflight_dir),
  mean_lfc = unique(mean_lfc),
  padj = unique(padj),
  class = unique(class),
  source = unique(source),
  max_abs_z = max(abs(compound_z)),
  max_reversal_score = max(reversal_score),
  z_apigenin_luteolin_1 = compound_z[compound_label == "apigenin+luteolin_1"],
  z_apigenin_luteolin_2 = compound_z[compound_label == "apigenin+luteolin_2"],
  z_licochalcone_A = compound_z[compound_label == "licochalcone_A"]
), by = gene]
setorder(t13c, -max_reversal_score)
fwrite(t13c, file.path(tab_dir, "T13c_reversed_genes_union.csv"))
cat(sprintf("[07] wrote T13c_reversed_genes_union.csv (%d unique genes reversed at |z|>=1.5)\n",
            nrow(t13c)))
cat(sprintf("[07] top 10 reversed genes:\n"))
print(head(t13c, 10)[, .(gene, spaceflight_dir, mean_lfc, max_abs_z, max_reversal_score, n_signatures_reversed)])

# ------------------------------------------------------------------------------
# 7. F11 -- gene-level reversal overview (bar + scatter)
# ------------------------------------------------------------------------------
# Panel A: grouped bar chart of reversal vs concordant counts at 3 thresholds
bar_dt <- rbindlist(list(
  t13[, .(compound_label, threshold = "|z|>=1.5", count = sum(reversed_z1.5), type = "reversed")],
  t13[, .(compound_label, threshold = "|z|>=1.5", count = sum(concordant_z1.5), type = "concordant")],
  t13[, .(compound_label, threshold = "|z|>=2", count = sum(reversed_z2), type = "reversed")],
  t13[, .(compound_label, threshold = "|z|>=2", count = sum(!reversed & abs(compound_z) >= 2), type = "concordant")],
  t13[, .(compound_label, threshold = "|z|>=3", count = sum(reversed_z3), type = "reversed")],
  t13[, .(compound_label, threshold = "|z|>=3", count = sum(!reversed & abs(compound_z) >= 3), type = "concordant")]
))
bar_dt[, threshold := factor(threshold, levels = c("|z|>=1.5", "|z|>=2", "|z|>=3"))]
bar_dt[, compound_label := factor(compound_label,
                                  levels = c("apigenin+luteolin_1", "apigenin+luteolin_2", "licochalcone_A"))]

pA <- ggplot(bar_dt, aes(x = threshold, y = count, fill = type)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = count), position = position_dodge(width = 0.7),
            vjust = -0.3, size = 3) +
  facet_wrap(~ compound_label, ncol = 3) +
  scale_fill_manual(values = c("reversed" = phylo_blue, "concordant" = phylo_orange),
                    name = "Direction") +
  labs(title = "Gene-level reversal of oncogenic biomarkers by flavonoids",
       subtitle = "Reversed = spaceflight direction opposite to compound z-score",
       x = "LINCS |z| threshold", y = "Number of T4 genes") +
  theme_minimal(base_family = "Liberation Sans") +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold", size = 10),
        strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

# Panel B: scatter spaceflight mean_lfc vs compound z (GSE128097_1)
scatter_dt <- t13[compound_label == "apigenin+luteolin_1"]
scatter_dt[, status := fcase(
  reversed_z1.5 == TRUE, "Reversed (|z|>=1.5)",
  concordant_z1.5 == TRUE, "Concordant (|z|>=1.5)",
  abs(compound_z) >= 1.5, "Concordant (|z|>=1.5)",
  reversed == TRUE, "Reversed (weak)",
  reversed == FALSE, "Concordant (weak)",
  default = "No effect"
)]
# Simplify: reversed strong, concordant strong, weak/no effect
scatter_dt[, status_simple := fcase(
  reversed_z1.5 == TRUE, "Reversed (|z|>=1.5)",
  concordant_z1.5 == TRUE, "Concordant (|z|>=1.5)",
  default = "Below threshold"
)]

pB <- ggplot(scatter_dt, aes(x = mean_lfc, y = compound_z, color = status_simple)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = c(-1.5, 1.5), linetype = "dashed", color = "grey70", linewidth = 0.3) +
  scale_color_manual(values = c("Reversed (|z|>=1.5)" = phylo_blue,
                                "Concordant (|z|>=1.5)" = phylo_orange,
                                "Below threshold" = phylo_grey),
                     name = "Status") +
  annotate("text", x = -4, y = 5, label = "Reversed quadrant\n(SF down, compound up)",
           hjust = 0, size = 3, color = phylo_blue, fontface = "bold") +
  annotate("text", x = 3, y = -5, label = "Reversed quadrant\n(SF up, compound down)",
           hjust = 0, size = 3, color = phylo_blue, fontface = "bold") +
  labs(title = "Apigenin + Luteolin (GSE128097_1)",
       subtitle = sprintf("%d reversed, %d concordant at |z|>=1.5",
                          sum(scatter_dt$reversed_z1.5), sum(scatter_dt$concordant_z1.5)),
       x = "Spaceflight mean log2FC (T4 oncogenic biomarkers)",
       y = "LINCS compound z-score") +
  theme_minimal(base_family = "Liberation Sans") +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

p11 <- pA / pB + plot_layout(heights = c(1, 1.2))
ggsave(file.path(fig_dir, "F11_gene_level_reversal.svg"), p11,
       width = 12, height = 10, dpi = 150)
ggsave(file.path(fig_dir, "F11_gene_level_reversal.png"), p11,
       width = 12, height = 10, dpi = 150)
cat("[07] wrote F11_gene_level_reversal.svg/png\n")

# ------------------------------------------------------------------------------
# 8. F12 -- heatmap of top 50 reversed genes x 3 signatures
# ------------------------------------------------------------------------------
# Build z-score matrix: genes (top 50 by max_reversal_score) x 3 signatures
top_genes <- head(t13c$gene, 50)
heat_dt <- t13[gene %in% top_genes]
# Cast to wide: gene x compound_label
heat_wide <- dcast(heat_dt, gene ~ compound_label,
                   value.var = "compound_z")
setorder(heat_wide, -gene)
# Order genes by max reversal score (already in t13c order)
gene_order <- t13c$gene[t13c$gene %in% top_genes]
heat_wide[, gene := factor(gene, levels = rev(gene_order))]
heat_wide <- heat_wide[order(gene)]

# Build annotation: spaceflight direction
anno <- t13c[gene %in% top_genes, .(gene, spaceflight_dir, mean_lfc)]
anno[, gene := factor(gene, levels = rev(gene_order))]
anno <- anno[order(gene)]

# Use ComplexHeatmap if available, else ggplot tile
has_ch <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
          requireNamespace("circlize", quietly = TRUE)

if (has_ch) {
  library(ComplexHeatmap)
  library(circlize)
  # Select the 3 compound columns by name. Column names contain '+', so we
  # must use .SDcols (character index) rather than NSE expressions.
  comp_cols <- c("apigenin+luteolin_1", "apigenin+luteolin_2", "licochalcone_A")
  mat_z <- as.matrix(heat_wide[, .SD, .SDcols = comp_cols])
  rownames(mat_z) <- as.character(heat_wide$gene)
  # Clip extreme values for color scale
  mat_z[mat_z > 5] <- 5
  mat_z[mat_z < -5] <- -5

  col_fun <- colorRamp2(c(-5, -1.5, 0, 1.5, 5),
                        c(phylo_blue, "#82B5E8", "white", "#FFC266", phylo_orange))

  # Row annotation: spaceflight direction
  row_anno <- rowAnnotation(
    `SF dir` = anno$spaceflight_dir,
    `SF lFC` = anno$mean_lfc,
    col = list(`SF dir` = c("up" = phylo_pink, "down" = phylo_green)),
    annotation_legend_param = list(
      `SF dir` = list(title = "Spaceflight dir"),
      `SF lFC` = list(title = "Spaceflight lFC")
    )
  )

  # Column labels
  col_labels <- c("apigenin+luteolin_1" = "Api+Lut\n(GSE128097_1)",
                  "apigenin+luteolin_2" = "Api+Lut\n(GSE128097_2)",
                  "licochalcone_A" = "Licochalcone A\n(GSE137934_1)")

  ht <- Heatmap(mat_z,
                name = "z-score",
                col = col_fun,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = 7, fontfamily = "Liberation Sans"),
                column_labels = col_labels[colnames(mat_z)],
                column_names_gp = gpar(fontsize = 8, fontfamily = "Liberation Sans"),
                column_names_rot = 0,
                left_annotation = row_anno,
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                column_title = "Top 50 reversed oncogenic biomarker genes",
                column_title_gp = gpar(fontsize = 11, fontface = "bold",
                                       fontfamily = "Liberation Sans"))

  svg_path <- file.path(fig_dir, "F12_reversed_genes_heatmap.svg")
  png_path <- file.path(fig_dir, "F12_reversed_genes_heatmap.png")
  svglite::svglite(svg_path, width = 8, height = 10)
  draw(ht)
  dev.off()
  ggsave_png <- function(path, ht, w, h) {
    png(path, width = w, height = h, units = "in", res = 150)
    draw(ht)
    dev.off()
  }
  ggsave_png(png_path, ht, 8, 10)
  cat("[07] wrote F12_reversed_genes_heatmap.svg/png (ComplexHeatmap)\n")
} else {
  # Fallback: ggplot tile
  heat_long <- melt(heat_dt[, .(gene, compound_label, compound_z)],
                    id.vars = c("gene", "compound_label"))
  heat_long[, gene := factor(gene, levels = rev(gene_order))]
  heat_long[, compound_label := factor(compound_label,
                                       levels = c("apigenin+luteolin_1",
                                                  "apigenin+luteolin_2",
                                                  "licochalcone_A"))]
  heat_long[compound_z > 5, compound_z := 5]
  heat_long[compound_z < -5, compound_z := -5]

  p12 <- ggplot(heat_long, aes(x = compound_label, y = gene, fill = compound_z)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(low = phylo_blue, mid = "white", high = phylo_orange,
                         midpoint = 0, limits = c(-5, 5),
                         name = "z-score") +
    labs(title = "Top 50 reversed oncogenic biomarker genes",
         x = NULL, y = NULL) +
    theme_minimal(base_family = "Liberation Sans") +
    theme(plot.title = element_text(face = "bold"),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(angle = 0, hjust = 0.5),
          panel.grid = element_blank())

  ggsave(file.path(fig_dir, "F12_reversed_genes_heatmap.svg"), p12,
         width = 8, height = 10, dpi = 150)
  ggsave(file.path(fig_dir, "F12_reversed_genes_heatmap.png"), p12,
         width = 8, height = 10, dpi = 150)
  cat("[07] wrote F12_reversed_genes_heatmap.svg/png (ggplot fallback)\n")
}

# ------------------------------------------------------------------------------
# 9. Console summary
# ------------------------------------------------------------------------------
cat("\n========== GENE-LEVEL REVERSAL SUMMARY ==========\n")
cat(sprintf("T4 oncogenic biomarkers: %d (up=%d, down=%d)\n",
            nrow(t4), sum(t4$mean_lfc > 0), sum(t4$mean_lfc < 0)))
cat(sprintf("T4 genes in LINCS: %d (%d missing)\n",
            length(t4_in_lincs), length(missing_genes)))
cat(sprintf("Flavonoid signatures: %d\n", length(flav_sigs)))
cat("\nPer-signature reversal counts:\n")
print(t13b)
cat(sprintf("\nUnion of reversed genes (|z|>=1.5): %d\n", nrow(t13c)))
cat(sprintf("Reversed by all 3 signatures: %d\n",
            sum(t13c$n_signatures_reversed == 3)))
cat(sprintf("Reversed by 2 signatures: %d\n",
            sum(t13c$n_signatures_reversed == 2)))
cat(sprintf("Reversed by 1 signature: %d\n",
            sum(t13c$n_signatures_reversed == 1)))
cat("\nTop 15 reversed genes (by max reversal score):\n")
print(head(t13c, 15)[, .(gene, spaceflight_dir, mean_lfc, max_abs_z,
                         max_reversal_score, n_signatures_reversed)])
cat("=================================================\n")

cat("\n[07_gene_level_reversal.R] DONE\n")

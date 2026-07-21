#' Script 03 — Drug screening: LINCS tau-analog reversal + Broad repurposing hub + PrimeKG drug context.
#'
#' The "opposite forcing" screen: find compounds whose perturbational transcriptomic
#' signature is anti-correlated with the spaceflight oncogenic biomarker signature.
#'
#' Pipeline:
#'   A. Load the R meta-signature (T3) -> define disease UP / DN gene sets.
#'   B. LINCS tau-analog: compute weighted connectivity score (WCS) of the biomarker
#'      signature against each drug perturbation signature in the datalake's
#'      human_geo_sigs.tsv (~380 drug-related signatures), then tau-normalize.
#'      Strongly NEGATIVE tau = reversal candidate ("opposite forcing").
#'   C. fgsea cross-check: rank biomarker signature, run fgsea against the 271-drug
#'      up/down GMT (single_drug_perturbations-v1.0.gmt). Reversal drugs show negative
#'      NES (biomarker-up enriched in drug-down and vice versa).
#'   D. Tissue-stratified reversal: per-tissue DEG lists -> per-tissue tau-analog.
#'   E. Drug context: annotate top reversal hits with Broad Drug Repurposing Hub
#'      (MoA/target/phase/indication) + PrimeKG (drug-target/contraindication).
#'
#' Outputs:
#'   results/tables/T6_lincs_reversal_scores.csv
#'   results/tables/T7_top_opposite_forcing.csv
#'   results/tables/T8_tissue_stratified_reversal.csv
#'   results/tables/T9_drug_context.csv
#'   results/figures/F6_lincs_reversal_ranking.svg/png
#'   results/figures/F7_tissue_stratified_reversal.svg/png
#'   results/figures/F9_drug_context.svg/png

source("R/lincs_reversal.R")
source("R/primekg.R")

PADJ <- 0.05
LFC  <- 1.0
TOP_N <- 50  # top reversal compounds to annotate

# Datalake paths (mounted read-only)
LINCS_SIGS <- "/mnt/datalake/LINCS1000/RNAseq_transcriptomics_signatures/human_geo_sigs.tsv"
LINCS_DRUG_GMT <- "/mnt/datalake/LINCS1000/RNAseq_transcriptomics_genesets/single_drug_perturbations-v1.0.gmt"
BROAD_PHASE <- "/mnt/datalake/broad_drug_repurposing_hub/broad_repurposing_hub_phase_moa_target_info.parquet"
BROAD_PHASE_CSV <- file.path("data", "processed", "broad_phase_moa_target.csv")
BROAD_MOL <- "/mnt/datalake/broad_drug_repurposing_hub/broad_repurposing_hub_molecule_with_smiles.parquet"
PRIMEKG_CSV <- "/mnt/datalake/primekg/primekg.csv"

main <- function() {
  root <- getwd()
  tab_dir <- file.path(root, "results", "tables")
  fig_dir <- file.path(root, "results", "figures")
  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- A. Load meta-signature, define disease UP / DN ----
  cat("=== Loading meta-signature ===\n")
  meta <- read.table(file.path(tab_dir, "T3_meta_signature.tsv"), header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
  # disease UP = genes significantly up-regulated in treated (spaceflight/radiation)
  # disease DN = genes significantly down-regulated
  sig <- meta[meta$padj < PADJ & abs(meta$mean_lfc) > LFC, ]
  disease_up <- sig$gene[sig$mean_lfc > 0]
  disease_dn <- sig$gene[sig$mean_lfc < 0]
  cat(sprintf("meta-signature: %d total, %d significant (padj<0.05, |lFC|>1)\n", nrow(meta), nrow(sig)))
  cat(sprintf("  disease UP: %d genes, disease DN: %d genes\n", length(disease_up), length(disease_dn)))

  # also build a ranked vector for fgsea (by mean_lfc, signed)
  ranked <- meta$mean_lfc
  names(ranked) <- meta$gene
  ranked <- sort(ranked, decreasing = TRUE)

  # ---- B. LINCS tau-analog reversal screen ----
  cat("\n=== LINCS tau-analog reversal screen ===\n")
  reversal <- run_reversal_screen(LINCS_SIGS, disease_up, disease_dn)
  write.csv(reversal, file.path(tab_dir, "T6_lincs_reversal_scores.csv"), row.names = FALSE)
  cat(sprintf("Wrote T6_lincs_reversal_scores.csv (%d drug signatures)\n", nrow(reversal)))

  # ---- C. fgsea cross-check ----
  cat("\n=== fgsea cross-check (271-drug GMT) ===\n")
  fgsea_res <- run_fgsea_crosscheck(ranked, LINCS_DRUG_GMT)
  if (!is.null(fgsea_res) && nrow(fgsea_res) > 0) {
    write.csv(fgsea_res, file.path(tab_dir, "T6b_fgsea_crosscheck.csv"), row.names = FALSE)
    cat(sprintf("Wrote T6b_fgsea_crosscheck.csv (%d drug sets)\n", nrow(fgsea_res)))
    # merge fgsea NES into the reversal table (by best-matching drug name)
    reversal <- merge_fgsea(reversal, fgsea_res)
  }

  # ---- D. Tissue-stratified reversal ----
  cat("\n=== Tissue-stratified reversal ===\n")
  tissue_rev <- run_tissue_stratified(meta, LINCS_SIGS, tab_dir)
  if (!is.null(tissue_rev) && nrow(tissue_rev) > 0) {
    write.csv(tissue_rev, file.path(tab_dir, "T8_tissue_stratified_reversal.csv"), row.names = FALSE)
    cat(sprintf("Wrote T8_tissue_stratified_reversal.csv (%d tissue-drug scores)\n", nrow(tissue_rev)))
  }

  # ---- E. Drug context annotation ----
  cat("\n=== Drug context annotation ===\n")
  # dedupe reversal by signature before taking top N (fgsea merge may have duplicated)
  reversal <- reversal[!duplicated(reversal$signature), ]
  top_drugs <- head(reversal[order(reversal$tau), ], TOP_N)
  # extract clean drug names (improved: pull compound after "treatment with" etc.)
  top_drugs$drug_name <- sapply(top_drugs$signature, extract_drug_name)
  # annotate by signature (unique key) to avoid merge duplication
  context <- annotate_drug_context(top_drugs$drug_name, top_drugs)
  context$signature <- top_drugs$signature  # carry the unique key
  write.csv(context, file.path(tab_dir, "T9_drug_context.csv"), row.names = FALSE)
  cat(sprintf("Wrote T9_drug_context.csv (%d drugs annotated)\n", nrow(context)))

  # combine top reversal + context -> T7 (merge on signature, the unique key)
  top_combined <- merge(top_drugs, context, by = "signature", all.x = TRUE, suffixes = c("",".ctx"))
  top_combined <- top_combined[order(top_combined$tau), ]
  write.csv(top_combined, file.path(tab_dir, "T7_top_opposite_forcing.csv"), row.names = FALSE)
  cat(sprintf("Wrote T7_top_opposite_forcing.csv (%d top compounds)\n", nrow(top_combined)))

  # ---- Figures ----
  cat("\n=== Generating figures ===\n")
  make_figures(reversal, top_combined, tissue_rev, context, fig_dir)

  cat("\nNext: Rscript scripts/04_enrichment.R\n")
}

# ---- C. fgsea cross-check ----

run_fgsea_crosscheck <- function(ranked, gmt_path) {
  if (!file.exists(gmt_path)) {
    cat("  GMT not found, skipping fgsea\n")
    return(NULL)
  }
  # parse GMT: each line = name \t description \t gene1 \t gene2 ...
  lines <- readLines(gmt_path)
  pathways <- lapply(lines, function(l) {
    parts <- strsplit(l, "\t")[[1]]
    if (length(parts) < 3) return(NULL)
    list(name = parts[1], genes = parts[-c(1,2)])
  })
  pathways <- Filter(function(x) !is.null(x), pathways)
  names(pathways) <- sapply(pathways, function(x) x$name)
  pathways <- lapply(pathways, function(x) x$genes)
  cat(sprintf("  loaded %d drug perturbation gene sets\n", length(pathways)))

  res <- fgsea::fgsea(pathways = pathways, stats = ranked, minSize = 5, maxSize = 500)
  res <- as.data.frame(res)
  res$leadingEdge <- sapply(res$leadingEdge, function(x) paste(x, collapse = "; "))
  res <- res[order(res$NES), ]
  res
}

#' Merge fgsea NES into the reversal table by matching drug names.
merge_fgsea <- function(reversal, fgsea_res) {
  # fgsea has -up and -dn sets per drug; compute a combined reversal NES
  # (average of -up NES and -dn NES, both expected negative for reversal)
  fgsea_res$base <- sub("-up$|-dn$", "", fgsea_res$pathway)
  fgsea_res$dir <- ifelse(grepl("-up$", fgsea_res$pathway), "up", "dn")
  # for reversal: biomarker-up should be enriched in drug-DN (negative NES for -up set? no)
  # Actually: fgsea ranks biomarker by lfc (up first). A drug-down set enriched at the TOP
  # of biomarker-up means the drug DOWN-regulates biomarker-UP genes = reversal => positive NES for -dn set.
  # Simpler: take the -dn set NES (drug-down genes vs biomarker ranking); positive NES = reversal.
  dn_sets <- fgsea_res[fgsea_res$dir == "dn", c("base","NES","padj")]
  colnames(dn_sets) <- c("drug_base","fgsea_NES_dn","fgsea_padj_dn")
  # match by lowercased drug name
  reversal$drug_base <- tolower(sub(" .*", "", reversal$signature))
  dn_sets$drug_base <- tolower(dn_sets$drug_base)
  merged <- merge(reversal, dn_sets, by = "drug_base", all.x = TRUE)
  merged$drug_base <- NULL
  merged
}

# ---- D. Tissue-stratified reversal ----

run_tissue_stratified <- function(meta, sigs_path, tab_dir) {
  # load per-study DE to get tissue-specific gene lists
  de_all <- read.csv(file.path(tab_dir, "T2_per_study_DE.csv"), stringsAsFactors = FALSE)
  # load sample metadata for study -> tissue mapping
  sm <- read.csv(file.path("data", "processed", "sample_metadata.csv"), stringsAsFactors = FALSE)
  study_tissue <- unique(sm[, c("study","tissue")])

  # for each study (tissue proxy), define UP/DN and run reversal
  tissues <- unique(de_all$study)
  cat(sprintf("  running reversal for %d tissues/studies\n", length(tissues)))

  # load the signature matrix once (it's 1.6 GB)
  cat("  loading LINCS signature matrix (shared across tissues)...\n")
  mat <- data.table::fread(sigs_path, header = TRUE, sep = "\t",
                           data.table = FALSE, nThread = parallel::detectCores())
  sig_names <- mat[[1]]; mat[[1]] <- NULL
  gene_names <- colnames(mat)
  is_drug <- grepl("drug|perturb|treatment|compound|inhibitor|agonist|exposure|dose|administer", sig_names, ignore.case = TRUE)
  sig_names_d <- sig_names[is_drug]
  mat_d <- mat[is_drug, , drop = FALSE]
  cat(sprintf("  matrix loaded: %d drug signatures x %d genes\n", nrow(mat_d), length(gene_names)))

  out <- list()
  for (tstudy in tissues) {
    de <- de_all[de_all$study == tstudy, ]
    tissue_name <- study_tissue$tissue[study_tissue$study == tstudy]
    if (is.na(tissue_name) || length(tissue_name) == 0) tissue_name <- tstudy
    up <- de$gene[de$padj < PADJ & de$log2FoldChange > LFC]
    dn <- de$gene[de$padj < PADJ & de$log2FoldChange < -LFC]
    if (length(up) < 10 || length(dn) < 10) {
      cat(sprintf("  [%s:%s] too few DEGs (up=%d, dn=%d), skipping\n", tstudy, tissue_name, length(up), length(dn)))
      next
    }
    cat(sprintf("  [%s:%s] up=%d, dn=%d\n", tstudy, tissue_name, length(up), length(dn)))
    # compute WCS per drug
    wcs <- sapply(seq_len(nrow(mat_d)), function(i) {
      z <- as.numeric(mat_d[i, ]); names(z) <- gene_names
      compute_wcs(z, up, dn)
    })
    valid <- !is.na(wcs)
    if (sum(valid) < 10) next
    tau <- tau_normalize(wcs[valid])
    df <- data.frame(study = tstudy, tissue = tissue_name,
                     signature = sig_names_d[valid], tau = tau,
                     stringsAsFactors = FALSE)
    out[[length(out)+1]] <- df
  }
  if (length(out) == 0) return(NULL)
  do.call(rbind, out)
}

# ---- E. Drug context annotation ----

annotate_drug_context <- function(drug_names, top_drugs) {
  # Broad Drug Repurposing Hub (read from CSV if available, else try parquet)
  cat("  loading Broad Drug Repurposing Hub...\n")
  broad <- NULL
  if (file.exists(BROAD_PHASE_CSV)) {
    broad <- tryCatch(read.csv(BROAD_PHASE_CSV, stringsAsFactors = FALSE),
                      error = function(e) { cat("  Broad CSV load failed:", conditionMessage(e), "\n"); NULL })
  }
  if (is.null(broad)) {
    broad <- tryCatch(
      data.table::fread(BROAD_PHASE, header = TRUE, data.table = FALSE),
      error = function(e) { cat("  Broad parquet load failed:", conditionMessage(e), "\n"); NULL }
    )
  }
  if (is.null(broad) || nrow(broad) == 0) {
    cat("  WARNING: Broad Repurposing Hub not loadable; using PrimeKG only\n")
  }
  # PrimeKG
  cat("  loading PrimeKG (filtered)...\n")
  pkg <- tryCatch(
    load_primekg_filtered(PRIMEKG_CSV, drug_names),
    error = function(e) { cat("  PrimeKG load failed:", conditionMessage(e), "\n"); NULL }
  )

  out <- data.frame(drug_name = drug_names, stringsAsFactors = FALSE)
  # match to Broad (case-insensitive partial)
  if (!is.null(broad) && nrow(broad) > 0) {
    broad$pert_iname_l <- tolower(broad$pert_iname)
    out$broad_match <- sapply(tolower(drug_names), function(d) {
      idx <- which(grepl(d, broad$pert_iname_l, fixed = TRUE) | grepl(broad$pert_iname_l, d, fixed = TRUE))
      if (length(idx) == 0) NA else broad$pert_iname[idx[1]]
    })
    out$moa <- sapply(seq_len(nrow(out)), function(i) {
      if (is.na(out$broad_match[i])) return(NA)
      idx <- which(broad$pert_iname == out$broad_match[i])
      if (length(idx) == 0) NA else broad$moa[idx[1]]
    })
    out$target <- sapply(seq_len(nrow(out)), function(i) {
      if (is.na(out$broad_match[i])) return(NA)
      idx <- which(broad$pert_iname == out$broad_match[i])
      if (length(idx) == 0) NA else broad$target[idx[1]]
    })
    out$clinical_phase <- sapply(seq_len(nrow(out)), function(i) {
      if (is.na(out$broad_match[i])) return(NA)
      idx <- which(broad$pert_iname == out$broad_match[i])
      if (length(idx) == 0) NA else broad$clinical_phase[idx[1]]
    })
    out$indication <- sapply(seq_len(nrow(out)), function(i) {
      if (is.na(out$broad_match[i])) return(NA)
      idx <- which(broad$pert_iname == out$broad_match[i])
      if (length(idx) == 0) NA else broad$indication[idx[1]]
    })
    out$disease_area <- sapply(seq_len(nrow(out)), function(i) {
      if (is.na(out$broad_match[i])) return(NA)
      idx <- which(broad$pert_iname == out$broad_match[i])
      if (length(idx) == 0) NA else broad$disease_area[idx[1]]
    })
  } else {
    out$broad_match <- out$moa <- out$target <- out$clinical_phase <- out$indication <- out$disease_area <- NA
  }

  # PrimeKG drug-target / contraindication
  if (!is.null(pkg)) {
    out$primekg_targets <- sapply(tolower(drug_names), function(d) {
      hits <- pkg$drug_target[grepl(d, pkg$drug_target$drug_name, ignore.case = TRUE) |
                              grepl(pkg$drug_target$drug_name, d, ignore.case = TRUE), ]
      if (nrow(hits) == 0) return(NA)
      paste(unique(c(hits$target_gene)), collapse = "; ")
    })
    out$primekg_contraindications <- sapply(tolower(drug_names), function(d) {
      hits <- pkg$drug_contra[grepl(d, pkg$drug_contra$drug_name, ignore.case = TRUE), ]
      if (nrow(hits) == 0) return(NA)
      paste(unique(hits$disease), collapse = "; ")
    })
  } else {
    out$primekg_targets <- out$primekg_contraindications <- NA
  }
  out
}

#' Load only the PrimeKG edges relevant to a set of drug names (filtered, memory-efficient).
load_primekg_filtered <- function(path, drug_names) {
  cat("  filtering PrimeKG for drug-target + contraindication edges...\n")
  # read only the columns we need
  df <- data.table::fread(path, header = TRUE, sep = ",",
                          select = c("relation","x_type","x_name","y_type","y_name"),
                          data.table = FALSE)
  drug_target <- df[(df$relation == "drug_protein" |
                     (df$x_type == "drug" & df$y_type == "gene/protein") |
                     (df$x_type == "gene/protein" & df$y_type == "drug")), ]
  if (nrow(drug_target) > 0) {
    drug_target$drug_name <- ifelse(drug_target$x_type == "drug", drug_target$x_name, drug_target$y_name)
    drug_target$target_gene <- ifelse(drug_target$x_type == "gene/protein", drug_target$x_name, drug_target$y_name)
    drug_target <- drug_target[, c("drug_name","target_gene")]
  }
  drug_contra <- df[df$relation == "contraindication", ]
  if (nrow(drug_contra) > 0) {
    drug_contra$drug_name <- drug_contra$x_name
    drug_contra$disease <- drug_contra$y_name
    drug_contra <- drug_contra[, c("drug_name","disease")]
  }
  cat(sprintf("  PrimeKG: %d drug-target, %d contraindication edges\n", nrow(drug_target), nrow(drug_contra)))
  list(drug_target = drug_target, drug_contra = drug_contra)
}

# ---- Figures ----

make_figures <- function(reversal, top_combined, tissue_rev, context, fig_dir) {
  th <- ggplot2::theme_minimal(base_family = "Liberation Sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 7))

  # F6: top 30 reversal compounds (most negative tau)
  top30 <- head(reversal[order(reversal$tau), ], 30)
  top30$short_name <- sapply(top30$signature, function(s) substr(s, 1, 50))
  # dedupe short names (truncation can create duplicates); append index if needed
  if (any(duplicated(top30$short_name))) {
    top30$short_name <- make.unique(top30$short_name)
  }
  top30$short_name <- factor(top30$short_name, levels = top30$short_name[order(top30$tau)])
  p6 <- ggplot2::ggplot(top30, ggplot2::aes(x = short_name, y = tau, fill = tau)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_gradient2(low = "#0279EE", mid = "grey90", high = "#c0392b", midpoint = 0) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Top 30 opposite-forcing compounds (LINCS tau-analog)",
                  x = NULL, y = "tau (negative = reversal candidate)", fill = "tau") +
    th + ggplot2::theme(legend.position = "right")
  save_fig(p6, file.path(fig_dir, "F6_lincs_reversal_ranking"))

  # F7: tissue-stratified reversal heatmap
  if (!is.null(tissue_rev) && nrow(tissue_rev) > 0) {
    # pick top 20 drugs by overall tau, show across tissues
    top_drug_names <- head(reversal$signature, 20)
    tr_sub <- tissue_rev[tissue_rev$signature %in% top_drug_names, ]
    if (nrow(tr_sub) > 0) {
      tr_sub$short_name <- sapply(tr_sub$signature, function(s) substr(s, 1, 40))
      # average tau per tissue x drug
      agg <- aggregate(tau ~ tissue + short_name, data = tr_sub, FUN = mean)
      p7 <- ggplot2::ggplot(agg, ggplot2::aes(x = tissue, y = short_name, fill = tau)) +
        ggplot2::geom_tile(color = "white") +
        ggplot2::scale_fill_gradient2(low = "#0279EE", mid = "grey95", high = "#c0392b", midpoint = 0) +
        ggplot2::labs(title = "Tissue-stratified reversal (tau-analog)",
                      x = "tissue", y = "drug", fill = "tau") +
        th + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8))
      save_fig(p7, file.path(fig_dir, "F7_tissue_stratified_reversal"))
    }
  }

  # F9: drug context summary (top compounds with MoA)
  if (!is.null(top_combined) && nrow(top_combined) > 0) {
    tc <- head(top_combined[order(top_combined$tau), ], 20)
    tc$short_name <- sapply(tc$signature, function(s) substr(s, 1, 40))
    if (any(duplicated(tc$short_name))) tc$short_name <- make.unique(tc$short_name)
    tc$short_name <- factor(tc$short_name, levels = tc$short_name[order(tc$tau)])
    tc$has_context <- !is.na(tc$moa)
    p9 <- ggplot2::ggplot(tc, ggplot2::aes(x = short_name, y = tau, fill = has_context)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_manual(values = c("TRUE" = "#75A025", "FALSE" = "grey70"),
                                 labels = c("TRUE" = "annotated", "FALSE" = "no context")) +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "Top 20 reversal compounds with drug-context availability",
                    x = NULL, y = "tau", fill = "context") + th
    save_fig(p9, file.path(fig_dir, "F9_drug_context"))
  }
}

save_fig <- function(plot, base_path) {
  svg_path <- paste0(base_path, ".svg")
  png_path <- paste0(base_path, ".png")
  svglite::svglite(svg_path, width = 10, height = 7)
  print(plot)
  dev.off()
  ggplot2::ggsave(png_path, plot, width = 10, height = 7, dpi = 150, bg = "white")
  cat(sprintf("  saved %s + %s\n", basename(svg_path), basename(png_path)))
}

main()

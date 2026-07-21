#' Script 04 — Disease enrichment: PrimeKG (local) + Open Targets (API).
#'
#' Maps the spaceflight oncogenic biomarker signature to known disease phenotypes,
#' providing statistical validation that the identified genes are biologically linked
#' to radiation- or microgravity-induced pathology (not just noise).
#'
#' Pipeline:
#'   A. Load the oncogenic biomarker intersection (T4) as the query gene set.
#'   B. PrimeKG: hypergeometric over-representation test of biomarker genes against
#'      each disease's associated gene set (~160k disease-gene edges).
#'   C. Open Targets: GraphQL API query for target-disease associations per biomarker
#'      gene; aggregate to a disease-level ranking.
#'   D. Cross-rank PrimeKG + Open Targets into a combined disease enrichment table.
#'
#' Outputs:
#'   results/tables/T10_disease_enrichment_primekg.csv
#'   results/tables/T11_disease_enrichment_opentargets.csv
#'   results/figures/F8_disease_enrichment.svg/png

source("R/primekg.R")
source("R/opentargets.R")

PADJ <- 0.05
PRIMEKG_CSV <- "/mnt/datalake/primekg/primekg.csv"
UNIVERSE_SIZE <- 20000  # approx total human protein-coding gene universe

main <- function() {
  root <- getwd()
  tab_dir <- file.path(root, "results", "tables")
  fig_dir <- file.path(root, "results", "figures")
  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- A. Load biomarker genes ----
  cat("=== Loading oncogenic biomarker genes ===\n")
  intersection <- read.csv(file.path(tab_dir, "T4_oncogene_intersection.csv"), stringsAsFactors = FALSE)
  biomarker_genes <- unique(intersection$gene)
  cat(sprintf("biomarker genes: %d (from T4 oncogene intersection)\n", length(biomarker_genes)))

  # also load the full significant meta-signature for a broader enrichment
  meta <- read.table(file.path(tab_dir, "T3_meta_signature.tsv"), header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
  sig_genes <- meta$gene[meta$padj < PADJ]
  cat(sprintf("full significant meta-signature genes: %d\n", length(sig_genes)))

  # ---- B. PrimeKG disease enrichment ----
  cat("\n=== PrimeKG disease enrichment ===\n")
  cat("loading PrimeKG (this takes ~60-90s)...\n")
  pkg <- tryCatch(load_primekg(PRIMEKG_CSV),
                  error = function(e) { cat("PrimeKG load failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(pkg)) {
    # enrichment on the oncogenic biomarker subset (primary)
    cat("\n[PrimeKG] enriching oncogenic biomarker subset...\n")
    enr_onc <- disease_enrichment(biomarker_genes, pkg, min_genes = 5, universe_size = UNIVERSE_SIZE)
    enr_onc$query_set <- "oncogenic_biomarkers"
    # enrichment on the full significant meta-signature (broader)
    cat("\n[PrimeKG] enriching full significant meta-signature...\n")
    enr_full <- disease_enrichment(sig_genes, pkg, min_genes = 5, universe_size = UNIVERSE_SIZE)
    enr_full$query_set <- "full_meta_signature"
    # combine
    enr_primekg <- rbind(enr_onc, enr_full)
    write.csv(enr_primekg, file.path(tab_dir, "T10_disease_enrichment_primekg.csv"), row.names = FALSE)
    cat(sprintf("Wrote T10_disease_enrichment_primekg.csv (%d disease associations)\n", nrow(enr_primekg)))
    cat(sprintf("  oncogenic biomarkers: %d significant diseases (padj<0.05)\n",
                sum(enr_onc$padj < 0.05, na.rm = TRUE)))
    cat(sprintf("  full meta-signature: %d significant diseases (padj<0.05)\n",
                sum(enr_full$padj < 0.05, na.rm = TRUE)))
    # show top 10
    cat("\n  Top 10 diseases (oncogenic biomarkers):\n")
    print(head(enr_onc[enr_onc$padj < 0.05, c("disease","n_overlap","n_disease_genes","padj")], 10))
  } else {
    enr_primekg <- NULL
  }

  # ---- C. Open Targets disease associations ----
  cat("\n=== Open Targets disease associations ===\n")
  # query the oncogenic biomarker genes (cap at 100 to limit API calls)
  query_genes <- biomarker_genes
  if (length(query_genes) > 100) {
    cat(sprintf("  capping Open Targets query to 100 genes (of %d)\n", length(query_genes)))
    query_genes <- head(query_genes, 100)
  }
  ot_assoc <- tryCatch(ot_batch(query_genes, sleep_sec = 0.3),
                       error = function(e) { cat("Open Targets failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(ot_assoc) && nrow(ot_assoc) > 0) {
    ot_ranking <- ot_disease_ranking(ot_assoc)
    write.csv(ot_ranking, file.path(tab_dir, "T11_disease_enrichment_opentargets.csv"), row.names = FALSE)
    cat(sprintf("Wrote T11_disease_enrichment_opentargets.csv (%d diseases)\n", nrow(ot_ranking)))
    cat("\n  Top 10 diseases (Open Targets):\n")
    print(head(ot_ranking[, c("disease_name","n_genes","mean_score","max_score")], 10))
  } else {
    ot_ranking <- NULL
    cat("Open Targets returned no associations\n")
  }

  # ---- D. Cross-rank PrimeKG + Open Targets ----
  cat("\n=== Cross-ranking PrimeKG + Open Targets ===\n")
  combined <- cross_rank_diseases(enr_primekg, ot_ranking)
  if (!is.null(combined) && nrow(combined) > 0) {
    write.csv(combined, file.path(tab_dir, "T10b_disease_enrichment_combined.csv"), row.names = FALSE)
    cat(sprintf("Wrote T10b_disease_enrichment_combined.csv (%d diseases)\n", nrow(combined)))
  }

  # ---- Figure ----
  cat("\n=== Generating figure ===\n")
  make_figure(enr_primekg, ot_ranking, combined, fig_dir)

  cat("\nNext: Rscript scripts/05_nutraceutical_flags.R\n")
}

#' Cross-rank diseases from PrimeKG and Open Targets.
#' @param enr_primekg data.frame from PrimeKG enrichment (or NULL)
#' @param ot_ranking data.frame from Open Targets (or NULL)
#' @return data.frame: disease, primekg_padj, primekg_n_overlap, ot_n_genes, ot_mean_score, combined_rank
cross_rank_diseases <- function(enr_primekg, ot_ranking) {
  # normalize disease names for matching (lowercase, strip whitespace)
  norm_name <- function(x) tolower(trimws(x))

  rows <- list()
  if (!is.null(enr_primekg)) {
    # use the oncogenic biomarker subset only
    enr_onc <- enr_primekg[enr_primekg$query_set == "oncogenic_biomarkers", ]
    enr_onc$disease_norm <- norm_name(enr_onc$disease)
    for (i in seq_len(nrow(enr_onc))) {
      rows[[enr_onc$disease_norm[i]]] <- list(
        disease = enr_onc$disease[i],
        primekg_padj = enr_onc$padj[i],
        primekg_n_overlap = enr_onc$n_overlap[i],
        primekg_n_disease_genes = enr_onc$n_disease_genes[i],
        ot_n_genes = NA_integer_, ot_mean_score = NA_real_, ot_max_score = NA_real_
      )
    }
  }
  if (!is.null(ot_ranking)) {
    ot_ranking$disease_norm <- norm_name(ot_ranking$disease_name)
    for (i in seq_len(nrow(ot_ranking))) {
      key <- ot_ranking$disease_norm[i]
      if (key %in% names(rows)) {
        rows[[key]]$ot_n_genes <- ot_ranking$n_genes[i]
        rows[[key]]$ot_mean_score <- ot_ranking$mean_score[i]
        rows[[key]]$ot_max_score <- ot_ranking$max_score[i]
      } else {
        rows[[key]] <- list(
          disease = ot_ranking$disease_name[i],
          primekg_padj = NA_real_, primekg_n_overlap = NA_integer_,
          primekg_n_disease_genes = NA_integer_,
          ot_n_genes = ot_ranking$n_genes[i],
          ot_mean_score = ot_ranking$mean_score[i],
          ot_max_score = ot_ranking$max_score[i]
        )
      }
    }
  }
  if (length(rows) == 0) return(NULL)
  out <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
  out$disease <- as.character(out$disease)
  # combined rank: rank by -log10(primekg_padj) + ot_n_genes (normalized)
  out$primekg_score <- ifelse(is.na(out$primekg_padj), 0, -log10(out$primekg_padj))
  out$ot_score <- ifelse(is.na(out$ot_n_genes), 0, out$ot_n_genes)
  # normalize each to [0,1] then sum
  if (max(out$primekg_score, na.rm = TRUE) > 0) {
    out$primekg_score <- out$primekg_score / max(out$primekg_score, na.rm = TRUE)
  }
  if (max(out$ot_score, na.rm = TRUE) > 0) {
    out$ot_score <- out$ot_score / max(out$ot_score, na.rm = TRUE)
  }
  out$combined_score <- out$primekg_score + out$ot_score
  out <- out[order(-out$combined_score), ]
  rownames(out) <- NULL
  out
}

# ---- Figure ----

make_figure <- function(enr_primekg, ot_ranking, combined, fig_dir) {
  th <- ggplot2::theme_minimal(base_family = "Liberation Sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"),
                   axis.text.y = ggplot2::element_text(size = 8))

  # F8: disease enrichment dot plot (top 20 by combined score)
  if (!is.null(combined) && nrow(combined) > 0) {
    top20 <- head(combined, 20)
    # truncate long disease names
    top20$disease_short <- sapply(top20$disease, function(x) {
      if (nchar(x) > 45) paste0(substr(x, 1, 42), "...") else x
    })
    top20$disease_short <- factor(top20$disease_short, levels = top20$disease_short[order(top20$combined_score)])
    p8 <- ggplot2::ggplot(top20, ggplot2::aes(x = combined_score, y = disease_short,
                                               size = primekg_n_overlap, color = ot_mean_score)) +
      ggplot2::geom_point() +
      ggplot2::scale_color_gradient(low = "#0279EE", high = "#c0392b", na.value = "grey80") +
      ggplot2::labs(title = "Disease enrichment of spaceflight oncogenic biomarkers\n(PrimeKG + Open Targets combined)",
                    x = "combined score (PrimeKG -log10 padj + Open Targets n_genes, normalized)",
                    y = NULL, size = "PrimeKG overlap", color = "OT mean score") +
      th + ggplot2::theme(legend.position = "right")
    save_fig(p8, file.path(fig_dir, "F8_disease_enrichment"))
  } else if (!is.null(enr_primekg)) {
    # fallback: PrimeKG only
    enr_onc <- enr_primekg[enr_primekg$query_set == "oncogenic_biomarkers" & enr_primekg$padj < 0.05, ]
    if (nrow(enr_onc) > 0) {
      top20 <- head(enr_onc[order(enr_onc$padj), ], 20)
      top20$disease_short <- sapply(top20$disease, function(x) {
        if (nchar(x) > 45) paste0(substr(x, 1, 42), "...") else x
      })
      top20$disease_short <- factor(top20$disease_short, levels = top20$disease_short[order(-log10(top20$padj))])
      p8 <- ggplot2::ggplot(top20, ggplot2::aes(x = -log10(padj), y = disease_short, size = n_overlap)) +
        ggplot2::geom_point(color = "#0279EE") +
        ggplot2::labs(title = "Disease enrichment of spaceflight oncogenic biomarkers (PrimeKG)",
                      x = "-log10 adjusted p", y = NULL, size = "overlap") + th
      save_fig(p8, file.path(fig_dir, "F8_disease_enrichment"))
    }
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

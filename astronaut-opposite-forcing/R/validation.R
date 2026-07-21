#' Validation: compare the new R meta-signature against the Python repo's S2.
#'
#' S2_cross_study_meta_signature.csv columns:
#'   gene, n_studies, stouffer_z, combined_p, mean_log2fc, direction_concordance,
#'   n_sig_up, n_sig_down, combined_padj
#'
#' The R meta-signature (meta_signature.tsv) columns:
#'   gene, n_studies, mean_lfc, se, z, p, padj, direction_concordance, n_species
#'
#' Concordance metrics:
#'   - Spearman rho of mean log2FC (R vs Python)
#'   - Sign agreement rate
#'   - Jaccard overlap of top-100 and top-500 genes (by |z| or padj)
#'   - Per-gene scatter

#' @param r_sig data.frame (R meta-signature)
#' @param s2_path path to S2_cross_study_meta_signature.csv
#' @return list(metrics = data.frame, joined = data.frame)
validate_vs_s2 <- function(r_sig, s2_path) {
  cat(sprintf("[validation] loading S2: %s\n", s2_path))
  s2 <- read.csv(s2_path, stringsAsFactors = FALSE)
  # join on gene symbol
  joined <- merge(r_sig[, c("gene","mean_lfc","z","padj")],
                  s2[, c("gene","mean_log2fc","stouffer_z","combined_padj")],
                  by = "gene", suffixes = c("_R","_S2"))
  cat(sprintf("[validation] joined genes: %d (R: %d, S2: %d)\n",
              nrow(joined), nrow(r_sig), nrow(s2)))

  # Spearman rho of log2FC
  rho <- cor(joined$mean_lfc, joined$mean_log2fc, method = "spearman", use = "pairwise.complete.obs")
  # Sign agreement
  sign_agree <- mean(sign(joined$mean_lfc) == sign(joined$mean_log2fc), na.rm = TRUE)
  # Jaccard top-100 / top-500 by |z|
  top_r_100  <- head(joined$gene[order(-abs(joined$z))], 100)
  top_s2_100 <- head(joined$gene[order(-abs(joined$stouffer_z))], 100)
  top_r_500  <- head(joined$gene[order(-abs(joined$z))], 500)
  top_s2_500 <- head(joined$gene[order(-abs(joined$stouffer_z))], 500)
  jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  j100 <- jacc(top_r_100, top_s2_100)
  j500 <- jacc(top_r_500, top_s2_500)

  metrics <- data.frame(
    metric = c("n_joined", "spearman_rho_lfc", "sign_agreement",
               "jaccard_top100", "jaccard_top500"),
    value = c(nrow(joined), round(rho, 4), round(sign_agree, 4),
              round(j100, 4), round(j500, 4))
  )
  cat("[validation] concordance metrics:\n")
  print(metrics)
  list(metrics = metrics, joined = joined)
}

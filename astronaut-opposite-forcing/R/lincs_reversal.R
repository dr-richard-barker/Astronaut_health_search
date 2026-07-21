#' LINCS L1000 "opposite forcing" reversal screen (tau-analog).
#'
#' The datalake does NOT contain the canonical CLUE.io Touchstone GCTx with
#' pre-computed tau-normalized connectivity scores. Instead it holds GEO-derived
#' drug perturbation signatures as a gene x signature z-score matrix
#' (human_geo_sigs.tsv, ~4270 signatures x ~20k genes). We compute a tau-analog
#' of the CMap weighted connectivity score (WCS) per Lamb (2006) / Subramanian
#' (2017), then tau-normalize across all drug signatures.
#'
#' A strongly NEGATIVE tau = the drug's perturbation profile is anti-correlated
#' with the disease signature = "opposite forcing" / reversal candidate.
#'
#' Method:
#'   For each drug signature d (a named vector of z-scores over genes):
#'     1. Rank d by absolute value (most extreme changes first).
#'     2. a_up  = weighted KS enrichment of disease-UP genes at the top of |d|-ranked d
#'        a_dn  = weighted KS enrichment of disease-DN genes at the top of |d|-ranked d
#'     3. WCS  = 0.5 * (a_up - a_dn)  if a_up and a_dn have the same sign (connectivity)
#'               else 0
#'     4. tau  = 100 * (WCS - mean(WCS_all)) / sd(WCS_all), clipped to [-100, 100]
#'
#' Negative tau => reversal candidate.

# ---- weighted KS enrichment (CMap a-score) ----
#' @param gene_ranks integer ranks of genes in the drug profile, ordered by |z| desc
#' @param gene_z signed z-scores of the drug profile, in the same gene order as gene_ranks
#' @param query_set character vector of disease UP (or DN) genes
#' @param universe character vector of all genes in the profile
#' @return numeric scalar a (weighted enrichment, signed)
weighted_ks_a <- function(gene_order, gene_z, query_set, universe) {
  # gene_order: genes ordered by |z| descending (the ranked list)
  in_q <- gene_order %in% query_set
  n_q <- sum(in_q)
  if (n_q == 0) return(0)
  # sign of the drug at each position
  s <- sign(gene_z[match(gene_order, names(gene_z))])
  # weighted cumulative enrichment: up-weight top-ranked genes
  # CMap uses a correlation-weighted KS; here we weight by rank position
  N <- length(gene_order)
  # V_j = sum of |z| over query genes up to position j, normalized
  weights <- abs(gene_z[match(gene_order, names(gene_z))])
  weights <- weights / sum(weights)
  cum_q <- cumsum(ifelse(in_q, weights, 0))
  cum_nonq <- cumsum(ifelse(!in_q, 1/(N - n_q), 0))
  a <- max(cum_q - cum_nonq)  # a^+ (enrichment at top)
  # signed by the dominant direction of the drug at query positions
  # CMap: a = a^+ if the query genes' drug z's are mostly positive, else -a^+
  q_sign <- mean(s[in_q])
  if (is.na(q_sign) || q_sign == 0) return(0)
  sign(q_sign) * a
}

#' Compute WCS for one drug signature against a disease signature.
#' @param drug_z named numeric vector (gene -> z-score) for the drug
#' @param disease_up character vector of disease UP genes
#' @param disease_dn character vector of disease DN genes
#' @return numeric WCS
compute_wcs <- function(drug_z, disease_up, disease_dn) {
  drug_z <- drug_z[!is.na(drug_z) & names(drug_z) != "" & is.finite(drug_z)]
  if (length(drug_z) < 50) return(NA_real_)
  universe <- names(drug_z)
  disease_up <- intersect(disease_up, universe)
  disease_dn <- intersect(disease_dn, universe)
  if (length(disease_up) < 5 || length(disease_dn) < 5) return(NA_real_)
  # rank by |z| descending
  ord <- order(abs(drug_z), decreasing = TRUE)
  gene_order <- names(drug_z)[ord]
  a_up <- weighted_ks_a(gene_order, drug_z, disease_up, universe)
  a_dn <- weighted_ks_a(gene_order, drug_z, disease_dn, universe)
  # connectivity condition: WCS nonzero only if a_up and a_dn have the same sign
  if (sign(a_up) == sign(a_dn) && a_up != 0 && a_dn != 0) {
    0.5 * (a_up - a_dn)
  } else {
    0
  }
}

#' tau-normalize a vector of WCS scores across all drugs.
#' @param wcs numeric vector
#' @return numeric vector of tau (clipped to [-100, 100])
tau_normalize <- function(wcs) {
  mu <- mean(wcs, na.rm = TRUE)
  sdv <- sd(wcs, na.rm = TRUE)
  if (sdv == 0 || is.na(sdv)) return(rep(0, length(wcs)))
  tau <- 100 * (wcs - mu) / sdv
  pmax(pmin(tau, 100), -100)
}

#' Run the full reversal screen.
#' @param sig_matrix path to human_geo_sigs.tsv (genes x signatures z-score matrix)
#' @param disease_up character vector of disease UP genes (e.g. up-regulated DEGs)
#' @param disease_dn character vector of disease DN genes
#' @param drug_filter_regex regex to select drug-related signatures (default matches
#'        drug|perturb|treatment|compound|inhibitor|agonist)
#' @return data.frame: signature, wcs, tau (sorted by tau ascending = strongest reversal first)
run_reversal_screen <- function(sig_matrix_path, disease_up, disease_dn,
                                drug_filter_regex = "drug|perturb|treatment|compound|inhibitor|agonist|exposure|dose|administer") {
  cat(sprintf("[lincs] loading signature matrix: %s\n", sig_matrix_path))
  # data.table::fread is fast for the 1.6 GB TSV
  mat <- data.table::fread(sig_matrix_path, header = TRUE, sep = "\t",
                           data.table = FALSE, nThread = parallel::detectCores())
  # first column = signature names
  sig_names <- mat[[1]]
  mat[[1]] <- NULL
  gene_names <- colnames(mat)
  cat(sprintf("[lincs] matrix: %d signatures x %d genes\n", length(sig_names), length(gene_names)))

  # filter to drug-related signatures
  is_drug <- grepl(drug_filter_regex, sig_names, ignore.case = TRUE)
  cat(sprintf("[lincs] drug-related signatures: %d / %d\n", sum(is_drug), length(sig_names)))
  sig_names_d <- sig_names[is_drug]
  mat_d <- mat[is_drug, , drop = FALSE]

  # build named z-vectors per signature
  cat("[lincs] computing WCS per drug signature...\n")
  wcs <- sapply(seq_len(nrow(mat_d)), function(i) {
    z <- as.numeric(mat_d[i, ])
    names(z) <- gene_names
    compute_wcs(z, disease_up, disease_dn)
  })
  names(wcs) <- sig_names_d
  valid <- !is.na(wcs)
  wcs <- wcs[valid]
  cat(sprintf("[lincs] valid WCS scores: %d\n", length(wcs)))

  tau <- tau_normalize(wcs)
  out <- data.frame(signature = names(wcs), wcs = wcs, tau = tau,
                    stringsAsFactors = FALSE)
  out <- out[order(out$tau), ]  # most negative tau (reversal) first
  rownames(out) <- NULL
  cat(sprintf("[lincs] top 5 reversal candidates (most negative tau):\n"))
  print(head(out, 5))
  out
}

#' Extract a clean drug name from a verbose GEO signature name.
#' @param sig_name e.g. "Transcriptome analysis of HEK 293 cells after treatment with Quinolinol B23 GSE130876_1"
#' @return short drug name (best-effort)
extract_drug_name <- function(sig_name) {
  # try patterns in priority order
  patterns <- c(
    "(?:treatment with|exposure to|treated with|administered|dose of)\\s+([A-Za-z0-9\\-]+)",
    "(?:inhibitor|agonist|antagonist)\\s+([A-Za-z0-9\\-]+)",
    "([A-Z][A-Za-z0-9\\-]{3,})\\s+(?:inhibitor|agonist|antagonist|treatment)"
  )
  for (pat in patterns) {
    m <- regmatches(sig_name, regexec(pat, sig_name, ignore.case = TRUE))
    if (length(m[[1]]) > 1 && nchar(m[[1]][2]) >= 3) return(m[[1]][2])
  }
  # fallback: extract the GSE ID (always present, unique)
  gse <- regmatches(sig_name, regexpr("GSE\\d+", sig_name))
  if (length(gse) > 0) return(gse[1])
  # last fallback: first capitalized word
  m2 <- regmatches(sig_name, regexpr("\\b[A-Z][A-Za-z0-9\\-]{3,}", sig_name))
  if (length(m2) > 0) return(m2[1])
  substr(sig_name, 1, 30)
}

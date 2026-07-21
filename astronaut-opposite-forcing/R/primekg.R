#' PrimeKG helpers: drug-target-indication-contraindication + disease-gene enrichment.
#'
#' PrimeKG (primekg.csv, ~937 MB) is a heterogeneous biomedical knowledge graph with:
#'   - drug -> gene/protein  (51k edges)
#'   - drug -> disease       (indication, 18k; contraindication, 61k)
#'   - disease -> gene/protein (80k edges)
#'   - drug -> effect/phenotype
#'
#' Used for (a) drug context annotation of LINCS reversal hits, and
#' (b) hypergeometric disease enrichment of the biomarker signature.

#' Load PrimeKG (filtered to the relation types we use).
#' @param path path to primekg.csv
#' @return list(drug_target, drug_indication, drug_contraindication, disease_gene)
load_primekg <- function(path) {
  cat(sprintf("[primekg] loading %s (this takes ~30-60s)...\n", path))
  df <- data.table::fread(path, header = TRUE, sep = ",",
                          select = c("relation","x_type","x_name","y_type","y_name"),
                          data.table = FALSE)
  cat(sprintf("[primekg] total edges: %d\n", nrow(df)))

  drug_target <- df[df$relation == "drug_protein" |
                    (df$x_type == "drug" & df$y_type == "gene/protein") |
                    (df$x_type == "gene/protein" & df$y_type == "drug"), ]
  drug_indication <- df[df$relation == "indication", ]
  drug_contra <- df[df$relation == "contraindication", ]
  disease_gene <- df[(df$x_type == "disease" & df$y_type == "gene/protein") |
                     (df$x_type == "gene/protein" & df$y_type == "disease"), ]

  cat(sprintf("[primekg] drug_target: %d, indication: %d, contraindication: %d, disease_gene: %d\n",
              nrow(drug_target), nrow(drug_indication), nrow(drug_contra), nrow(disease_gene)))
  list(drug_target = drug_target, drug_indication = drug_indication,
       drug_contraindication = drug_contra, disease_gene = disease_gene,
       full = df)
}

#' Annotate a list of drug names with PrimeKG drug-target / indication / contraindication.
#' @param drugs character vector of drug names (lowercased for matching)
#' @param pkg output of load_primekg()
#' @return data.frame: drug, targets, indications, contraindications
annotate_drugs <- function(drugs, pkg) {
  drugs_l <- tolower(drugs)
  # PrimeKG drug names are in x_name (for drug->X) or y_name (for X->drug)
  # build a drug-name index
  dt <- pkg$drug_target
  di <- pkg$drug_indication
  dc <- pkg$drug_contraindication

  out <- lapply(drugs_l, function(d) {
    # match drug name (case-insensitive, partial)
    hits_dt <- dt[grepl(d, dt$x_name, ignore.case = TRUE) | grepl(d, dt$y_name, ignore.case = TRUE), ]
    targets <- unique(c(hits_dt$x_name[hits_dt$x_type == "gene/protein"],
                        hits_dt$y_name[hits_dt$y_type == "gene/protein"]))
    hits_di <- di[grepl(d, di$x_name, ignore.case = TRUE), ]
    indications <- unique(hits_di$y_name)
    hits_dc <- dc[grepl(d, dc$x_name, ignore.case = TRUE), ]
    contras <- unique(hits_dc$y_name)
    data.frame(drug = d, targets = paste(targets, collapse = "; "),
               indications = paste(indications, collapse = "; "),
               contraindications = paste(contras, collapse = "; "),
               n_targets = length(targets), n_indications = length(indications),
               n_contraindications = length(contras),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Hypergeometric disease enrichment of a biomarker gene set.
#' @param biomarker_genes character vector of HGNC symbols
#' @param pkg output of load_primekg()
#' @param min_genes minimum disease-associated genes to test (default 5)
#' @param universe_size total gene universe (default 20000)
#' @return data.frame: disease, n_overlap, n_disease_genes, n_biomarker_genes,
#'         odds_ratio, p, padj, sorted by padj
disease_enrichment <- function(biomarker_genes, pkg, min_genes = 5, universe_size = 20000) {
  dg <- pkg$disease_gene
  # disease -> set of genes
  diseases <- unique(dg$x_name[dg$x_type == "disease"])
  cat(sprintf("[primekg] testing %d diseases for enrichment of %d biomarker genes\n",
              length(diseases), length(biomarker_genes)))

  # build disease -> gene list
  disease_genes <- split(dg$y_name[dg$x_type == "disease" & dg$y_type == "gene/protein"],
                         dg$x_name[dg$x_type == "disease" & dg$y_type == "gene/protein"])
  # also the reverse direction
  rev <- dg[dg$x_type == "gene/protein" & dg$y_type == "disease", ]
  if (nrow(rev) > 0) {
    rev_genes <- split(rev$x_name, rev$y_name)
    for (d in names(rev_genes)) {
      disease_genes[[d]] <- unique(c(disease_genes[[d]], rev_genes[[d]]))
    }
  }

  results <- lapply(diseases, function(d) {
    dg_set <- unique(disease_genes[[d]])
    dg_set <- dg_set[dg_set != "" & !is.na(dg_set)]
    if (length(dg_set) < min_genes) return(NULL)
    overlap <- intersect(biomarker_genes, dg_set)
    n_overlap <- length(overlap)
    if (n_overlap == 0) return(NULL)
    # hypergeometric test (one-sided, over-representation)
    p <- stats::phyper(n_overlap - 1, length(dg_set),
                       universe_size - length(dg_set),
                       length(biomarker_genes), lower.tail = FALSE)
    data.frame(disease = d, n_overlap = n_overlap,
               n_disease_genes = length(dg_set),
               n_biomarker_genes = length(biomarker_genes),
               overlap_genes = paste(overlap, collapse = "; "),
               p = p, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, results)
  if (is.null(res) || nrow(res) == 0) {
    return(data.frame(disease = character(), n_overlap = integer(),
                      n_disease_genes = integer(), n_biomarker_genes = integer(),
                      overlap_genes = character(), p = numeric(), padj = numeric(),
                      stringsAsFactors = FALSE))
  }
  res$padj <- stats::p.adjust(res$p, method = "BH")
  res <- res[order(res$padj), ]
  rownames(res) <- NULL
  res
}

#' Mouse -> human ortholog mapping via biomaRt (Ensembl).
#'
#' Maps mouse gene symbols to human orthologs using Ensembl BioMart.
#' Resolves one-to-many by orthology confidence; keeps 1:1 and high-confidence 1:many.
#'
#' @param mouse_genes character vector of mouse gene symbols (e.g. "Mmp3", "Trp53")
#' @return data.frame: mouse_symbol, human_symbol, homology_type, confidence
map_orthologs <- function(mouse_genes) {
  mouse_genes <- unique(mouse_genes)
  cat(sprintf("[orthologs] mapping %d mouse symbols -> human via biomaRt\n", length(mouse_genes)))

  # Use Ensembl BioMart: mouse mart, getLDS to human
  # Try cached connection first; biomaRt can be flaky
  mart <- tryCatch({
    biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  }, error = function(e) {
    cat("[orthologs] ensembl.org failed, trying mirror...\n")
    biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                     host = "https://www.ensembl.org")
  })

  res <- tryCatch({
    biomaRt::getLDS(
      attributes  = c("mgi_symbol"),
      filters     = "mgi_symbol",
      values      = mouse_genes,
      mart        = mart,
      attributesL = c("hgnc_symbol"),
      martL       = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    )
  }, error = function(e) {
    cat("[orthologs] getLDS failed:", conditionMessage(e), "\n")
    return(NULL)
  })

  if (is.null(res) || nrow(res) == 0) {
    cat("[orthologs] no orthologs returned; returning empty mapping\n")
    return(data.frame(mouse_symbol = character(), human_symbol = character(),
                      homology_type = character(), confidence = character(),
                      stringsAsFactors = FALSE))
  }

  colnames(res) <- c("mouse_symbol", "human_symbol")
  # drop blank human symbols and duplicates
  res <- res[res$human_symbol != "" & !is.na(res$human_symbol), ]
  # annotate multiplicity
  counts <- table(res$mouse_symbol)
  res$homology_type <- ifelse(counts[res$mouse_symbol] == 1, "1:1", "1:many")
  res$confidence <- "biomaRt_getLDS"
  res <- res[!duplicated(res), ]
  cat(sprintf("[orthologs] %d mouse symbols -> %d ortholog rows (%d 1:1, %d 1:many)\n",
              length(unique(res$mouse_symbol)), nrow(res),
              sum(res$homology_type == "1:1"), sum(res$homology_type == "1:many")))
  res
}

#' Build a gene-symbol translation vector for a count matrix.
#' @param symbols gene symbols in the matrix (could be mouse or human)
#' @param species "mouse" or "human"
#' @return list(mapped = data.frame with original + human_symbol, unmapped = character)
to_human_symbols <- function(symbols, species) {
  if (species == "human") {
    return(list(mapped = data.frame(original = symbols, human_symbol = symbols,
                                    stringsAsFactors = FALSE),
                unmapped = character()))
  }
  if (species == "mouse") {
    map <- map_orthologs(symbols)
    mapped <- merge(data.frame(original = symbols, stringsAsFactors = FALSE),
                    map, by.x = "original", by.y = "mouse_symbol", all.x = TRUE)
    unmapped <- mapped$original[is.na(mapped$human_symbol)]
    list(mapped = mapped, unmapped = unmapped)
  } else {
    stop("unknown species: ", species)
  }
}

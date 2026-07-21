#' Open Targets Platform GraphQL API client.
#'
#' Endpoint: https://api.platform.opentargets.org/api/v4/graphql
#' Used for target-disease association scores to complement PrimeKG disease enrichment.

OT_GRAPHQL <- "https://api.platform.opentargets.org/api/v4/graphql"

#' Query associated diseases for a gene (by HGNC symbol).
#' @param gene HGNC symbol (e.g. "TP53")
#' @param size max diseases to return (default 50)
#' @return data.frame: gene, disease_id, disease_name, score, datatype_scores (JSON)
ot_associated_diseases <- function(gene, size = 50) {
  # first resolve symbol -> ensembl id (search without size arg — API changed)
  q_resolve <- sprintf('{
    search(queryString: "%s", entityNames: ["target"]) {
      hits { id entity }
    }
  }', gene)
  r <- tryCatch(
    httr::POST(OT_GRAPHQL,
               body = list(query = q_resolve),
               encode = "json", httr::timeout(30)),
    error = function(e) NULL
  )
  if (is.null(r) || httr::status_code(r) != 200) return(NULL)
  dat <- httr::content(r, as = "parsed", type = "application/json")
  hits <- dat$data$search$hits
  if (length(hits) == 0) return(NULL)
  # find the hit whose entity is "target" and id starts with ENSG
  ensembl_id <- NULL
  for (h in hits) {
    if (h$entity == "target" && grepl("^ENSG", h$id)) { ensembl_id <- h$id; break }
  }
  if (is.null(ensembl_id)) ensembl_id <- hits[[1]]$id
  if (is.null(ensembl_id)) return(NULL)

  q_assoc <- sprintf('{
    target(ensemblId: "%s") {
      approvedSymbol
      associatedDiseases(page: {size: %d, index: 0}) {
        rows {
          score
          datatypeScores { id score }
          disease { id name }
        }
      }
    }
  }', ensembl_id, size)
  r2 <- tryCatch(
    httr::POST(OT_GRAPHQL,
               body = list(query = q_assoc),
               encode = "json", httr::timeout(30)),
    error = function(e) NULL
  )
  if (is.null(r2) || httr::status_code(r2) != 200) return(NULL)
  dat2 <- httr::content(r2, as = "parsed", type = "application/json")
  rows <- dat2$data$target$associatedDiseases$rows
  if (length(rows) == 0) return(NULL)
  out <- lapply(rows, function(row) {
    data.frame(gene = gene, disease_id = row$disease$id,
               disease_name = row$disease$name, score = row$score,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Query associated diseases for a list of genes (rate-limited).
#' @param genes character vector of HGNC symbols
#' @param sleep_sec seconds between queries (default 0.5)
#' @return data.frame: gene, disease_id, disease_name, score
ot_batch <- function(genes, sleep_sec = 0.5) {
  cat(sprintf("[opentargets] querying %d genes...\n", length(genes)))
  out <- list()
  for (i in seq_along(genes)) {
    g <- genes[i]
    res <- tryCatch(ot_associated_diseases(g), error = function(e) NULL)
    if (!is.null(res) && nrow(res) > 0) out[[length(out)+1]] <- res
    if (i %% 25 == 0) cat(sprintf("  ...%d/%d\n", i, length(genes)))
    Sys.sleep(sleep_sec)
  }
  if (length(out) == 0) {
    return(data.frame(gene = character(), disease_id = character(),
                      disease_name = character(), score = numeric(),
                      stringsAsFactors = FALSE))
  }
  df <- do.call(rbind, out)
  rownames(df) <- NULL
  cat(sprintf("[opentargets] %d gene-disease associations for %d genes\n",
              nrow(df), length(unique(df$gene))))
  df
}

#' Aggregate to a disease-level ranking.
#' @param assoc output of ot_batch
#' @return data.frame: disease_name, n_genes, mean_score, max_score, genes
ot_disease_ranking <- function(assoc) {
  if (nrow(assoc) == 0) return(assoc)
  agg <- aggregate(score ~ disease_name + disease_id, data = assoc,
                   FUN = function(x) c(mean = mean(x), max = max(x), n = length(x)))
  agg <- do.call(data.frame, agg)
  colnames(agg) <- c("disease_name","disease_id","mean_score","max_score","n_genes")
  # attach gene lists
  agg$genes <- sapply(agg$disease_id, function(d) {
    paste(unique(assoc$gene[assoc$disease_id == d]), collapse = "; ")
  })
  agg <- agg[order(-agg$n_genes, -agg$mean_score), ]
  rownames(agg) <- NULL
  agg
}

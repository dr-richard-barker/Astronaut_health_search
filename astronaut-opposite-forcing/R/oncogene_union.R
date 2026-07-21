#' Build the union oncogene / tumor-suppressor reference list.
#'
#' Sources (in priority order; each is tried and recorded as succeeded/failed):
#'   1. COSMIC Cancer Gene Census (open-access TSV from cancer.sanger.ac.uk)
#'   2. OncoKB open gene list (if downloadable without login)
#'   3. MSigDB C6 oncogenic signatures (189 sets, via msigdbr) -> gene union
#'   4. KEGG oncogenesis pathways (via clusterProfiler / KEGG REST)
#'   5. Curated oncogene / TSG / SASP list (extending the discovery repo's set)
#'
#' @return list(genes = data.frame(gene, class, source), log = data.frame(source, status, n_genes))

build_oncogene_union <- function() {
  log_rows <- list()
  gene_rows <- list()

  add_log <- function(source, status, n_genes) {
    log_rows[[length(log_rows) + 1]] <<- data.frame(source = source, status = status,
                                                     n_genes = n_genes, stringsAsFactors = FALSE)
  }
  add_genes <- function(genes, class, source) {
    if (length(genes) == 0) return(invisible(NULL))
    gene_rows[[length(gene_rows) + 1]] <<- data.frame(gene = genes, class = class,
                                                       source = source, stringsAsFactors = FALSE)
  }

  # ---- 1. COSMIC CGC open-access ----
  cosmic_res <- tryCatch({
    # COSMIC CGC open-access: the cancer_gene_census.csv requires a COSMIC login.
    # Try the public download endpoint; if it fails, fall back to the CGC TSV via
    # the Sanger downloads page (sometimes accessible without login).
    urls <- c(
      "https://cancer.sanger.ac.uk/cosmic/file_download/GRCh38/cancer_gene_census.csv",
      "https://cancer.sanger.ac.uk/cosmic-download/GRCh38/cosmic/v97/CancerGeneCensus.csv"
    )
    for (url in urls) {
      r <- tryCatch(httr::GET(url, httr::timeout(60),
                              httr::user_agent("astronaut-opposite-forcing/0.1 (R)")),
                    error = function(e) NULL)
      if (!is.null(r) && httr::status_code(r) == 200) {
        txt <- httr::content(r, as = "text", encoding = "UTF-8")
        df <- tryCatch(read.csv(text = txt, stringsAsFactors = FALSE), error = function(e) NULL)
        if (!is.null(df) && "Gene.Symbol" %in% colnames(df)) {
          genes <- df$Gene.Symbol[df$Gene.Symbol != "" & !is.na(df$Gene.Symbol)]
          cls <- ifelse(grepl("oncogene", df$Role.in.Cancer, ignore.case = TRUE), "oncogene",
                        ifelse(grepl("TSG|tumour suppressor", df$Role.in.Cancer, ignore.case = TRUE), "TSG",
                               "cancer_gene"))
          add_log("COSMIC_CGC", "OK", length(genes))
          add_genes(genes, cls, "COSMIC_CGC")
          return(TRUE)
        }
      }
    }
    add_log("COSMIC_CGC", "FAILED_all_urls", 0)
    FALSE
  }, error = function(e) {
    add_log("COSMIC_CGC", paste0("ERR:", conditionMessage(e)), 0)
    FALSE
  })

  # ---- 2. OncoKB open gene list ----
  onkokb_res <- tryCatch({
    # OncoKB public API (may require registration; try without auth first)
    r <- httr::GET("https://www.oncokb.org/api/v1/utils/allActionableGenes",
                   httr::timeout(60),
                   httr::add_headers(`Authorization` = "Bearer dummy"))
    if (httr::status_code(r) == 200) {
      genes <- unlist(httr::content(r, as = "parsed"))
      add_log("OncoKB", "OK", length(genes))
      add_genes(genes, "actionable_oncogene", "OncoKB")
      TRUE
    } else {
      add_log("OncoKB", paste0("HTTP_", httr::status_code(r)), 0)
      FALSE
    }
  }, error = function(e) {
    add_log("OncoKB", paste0("ERR:", conditionMessage(e)), 0)
    FALSE
  })

  # ---- 3. MSigDB C6 oncogenic signatures (always available via msigdbr) ----
  c6_res <- tryCatch({
    m <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C6")
    genes <- unique(m$gene_symbol)
    add_log("MSigDB_C6", "OK", length(genes))
    add_genes(genes, "C6_oncogenic", "MSigDB_C6")
    TRUE
  }, error = function(e) {
    add_log("MSigDB_C6", paste0("ERR:", conditionMessage(e)), 0)
    FALSE
  })

  # ---- 4. KEGG oncogenesis pathways ----
  kegg_res <- tryCatch({
    # Use KEGGREST::keggGet to fetch pathway entries, parse GENE field
    paths <- c("hsa05200","hsa05202","hsa05206","hsa05210","hsa05224","hsa05223","hsa05219")
    all_genes <- character()
    for (pid in paths) {
      g <- tryCatch(KEGGREST::keggGet(pid), error = function(e) NULL)
      if (!is.null(g) && "GENE" %in% names(g[[1]])) {
        gene_field <- g[[1]]$GENE
        # entries alternate: id, "symbol; description"
        defs <- gene_field[seq(2, length(gene_field), by = 2)]
        syms <- sub(";.*", "", defs)
        all_genes <- c(all_genes, syms)
      }
    }
    all_genes <- unique(all_genes[all_genes != ""])
    add_log("KEGG_oncogenesis", "OK", length(all_genes))
    add_genes(all_genes, "KEGG_oncogenesis", "KEGG")
    TRUE
  }, error = function(e) {
    add_log("KEGG_oncogenesis", paste0("ERR:", conditionMessage(e)), 0)
    FALSE
  })

  # ---- 5. Curated oncogene / TSG / SASP list (extending discovery repo) ----
  curated_oncogenes <- c(
    "KRAS","NRAS","HRAS","MYC","MYCN","MYCL","EGFR","ERBB2","MET","BRAF",
    "PIK3CA","AKT1","AKT2","CCND1","CCNE1","CDK4","CDK6","MDM2","MDM4",
    "BCL2","MCL1","JUN","FOS","MYB","ETS1","CTNNB1","NOTCH1","JAK2",
    "STAT3","FLT3","KIT","ABL1","YAP1","E2F1","SRC","RAF1","MTOR","GLI1",
    "WNT1","FGFR1","PDGFRA","REL","SKP2","AURKA","PLK1","BIRC5",
    "ERBB3","ERBB4","PIK3CB","PIK3CD","PIK3CG","AKT3","CCND2","CCND3",
    "CDK2","CDK1","CHEK1","CHEK2","WEE1","MAP2K1","MAP2K2","MAPK1","MAPK3",
    "ARAF","GNA11","GNAQ","GNAS","SMO","SUFU","IDH1","IDH2",
    "FGFR2","FGFR3","FGFR4","RET","NTRK1","NTRK2","NTRK3","ALK","ROS1",
    "AXL","MER","TYRO3","NPM1","TERT","PROM1","CD44","ALDH1A1","SOX2","NANOG"
  )
  curated_tsgs <- c(
    "TP53","RB1","PTEN","CDKN2A","CDKN1A","CDKN1B","CDKN2B","APC","VHL",
    "BRCA1","BRCA2","NF1","NF2","SMAD4","STK11","ATM","TSC1","TSC2",
    "WT1","MEN1","PTCH1","FBXW7","KEAP1","MLH1","MSH2","DCC","FH",
    "SDHB","CDH1","GADD45A","BAX","PML","RUNX3",
    "ARID1A","BAP1","BMPR1A","CDKN2C","FANCA","FANCC","FANCD2",
    "PALB2","RAD51","RAD51C","SLX4","XRCC2","RNF43","SMAD2","SMAD3",
    "ACVR2A","RHOA","KMT2C","KDM6A","STAG2","BCOR",
    "SDHA","SDHC","SDHD"
  )
  curated_sasp <- c(
    "IL6","IL1A","IL1B","CXCL8","IL8","CXCL1","CXCL2","CXCL10","CCL2",
    "CCL20","MMP1","MMP3","MMP10","MMP12","MMP19","TIMP1","TIMP2",
    "IGFBP3","IGFBP7","SERPINE1","IL15","CSF2","HGF","FGF2","VEGFA",
    "PLAU","PLAUR","ICAM1","TNFRSF1A","IL13","CXCL12",
    "GDF15","INHBA","CSF1","IGF1","IGF2","DLL4","NOTCH1","EFNB2","EPHB4",
    "TGFB1","TGFB2","TGFB3","IL18","IL33","CXCL13","CCL5","CCL19","CCL21"
  )
  add_log("Curated", "OK", length(unique(c(curated_oncogenes, curated_tsgs, curated_sasp))))
  add_genes(curated_oncogenes, "oncogene", "Curated")
  add_genes(curated_tsgs, "TSG", "Curated")
  add_genes(curated_sasp, "SASP", "Curated")

  # ---- combine ----
  all_genes <- do.call(rbind, gene_rows)
  log_df <- do.call(rbind, log_rows)
  if (is.null(all_genes) || nrow(all_genes) == 0) {
    return(list(genes = data.frame(gene = character(), class = character(), source = character(),
                                   stringsAsFactors = FALSE),
                log = log_df))
  }
  # dedupe keeping the most informative class (oncogene > TSG > SASP > C6 > KEGG > cancer_gene)
  class_rank <- c(oncogene = 1, TSG = 2, SASP = 3, C6_oncogenic = 4, KEGG_oncogenesis = 5,
                  actionable_oncogene = 1, cancer_gene = 6)
  all_genes$rank <- class_rank[all_genes$class]
  all_genes <- all_genes[order(all_genes$gene, all_genes$rank), ]
  all_genes <- all_genes[!duplicated(all_genes$gene), ]
  all_genes$rank <- NULL

  cat("\n[oncogene_union] sources:\n")
  print(log_df)
  cat(sprintf("[oncogene_union] total unique genes: %d\n", nrow(all_genes)))
  list(genes = all_genes, log = log_df)
}

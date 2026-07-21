#' Script 02 — Biomarker identification: DESeq2 + ortholog mapping + metafor meta + oncogene union.
#'
#' Pipeline:
#'   A. Load each study's STAR unnormalized counts + runsheet (or sample-name parsing).
#'   B. Assign condition (control vs treated) from runsheet Factor Value fields or sample names.
#'   C. For rodent studies: map Ensembl mouse gene IDs -> human HGNC symbols via biomaRt.
#'   D. Run DESeq2 per study (treated vs control), strict thresholds: padj < 0.05, |log2FC| > 1.0.
#'   E. metafor random-effects meta-analysis per gene (accounts for cross-species/tissue heterogeneity).
#'   F. Build oncogene union list (COSMIC CGC open-access + OncoKB + MSigDB C6 + KEGG + curated).
#'   G. Intersect meta-signature with oncogene union -> oncogenic biomarker subset.
#'   H. Validate R signature vs Python S2 (concordance).
#'
#' Outputs:
#'   results/tables/T2_per_study_DE.csv       (per-study DESeq2 results, concatenated)
#'   results/tables/T3_meta_signature.tsv     (R random-effects meta-signature)
#'   results/tables/T4_oncogene_intersection.csv
#'   results/tables/T5_validation_concordance.csv
#'   results/figures/F3_meta_signature.svg/png
#'   results/figures/F4_oncogene_intersection.svg/png
#'   results/figures/F5_validation_concordance.svg/png
#'   data/processed/sample_metadata.csv       (full sample metadata with conditions)

source("R/orthologs.R")
source("R/oncogene_union.R")
source("R/validation.R")

PADJ <- 0.05
LFC  <- 1.0
MIN_SAMPLES_PER_GROUP <- 2

# ---- A/B. Load counts + assign conditions ----

#' Find the STAR unnormalized counts file for a study.
find_counts_file <- function(study_dir) {
  # prefer STAR unnormalized counts (gene-level, raw counts for DESeq2)
  candidates <- list.files(study_dir, pattern = "STAR_Unnormalized_Counts.*\\.csv$", full.names = TRUE)
  if (length(candidates) == 0) {
    # fallback: any unnormalized counts
    candidates <- list.files(study_dir, pattern = "Unnormalized_Counts.*\\.csv$", full.names = TRUE)
  }
  if (length(candidates) == 0) {
    candidates <- list.files(study_dir, pattern = ".*_Unnormalized_Counts\\.csv$", full.names = TRUE)
  }
  if (length(candidates) == 0) return(NA_character_)
  candidates[1]  # take first match
}

#' Find the runsheet for a study.
find_runsheet <- function(study_dir) {
  candidates <- list.files(study_dir, pattern = "runsheet.*\\.csv$", full.names = TRUE)
  if (length(candidates) == 0) return(NA_character_)
  candidates[1]
}

#' Assign condition (control vs treated) from runsheet Factor Value fields.
#' @param runsheet data.frame
#' @return data.frame: sample_name, condition (control/treated)
assign_condition_from_runsheet <- function(runsheet) {
  # find Factor Value columns
  fcols <- grep("^Factor Value", colnames(runsheet), value = TRUE)
  if (length(fcols) == 0) return(NULL)
  # the sample name column
  samp_col <- grep("^Sample Name$", colnames(runsheet), value = TRUE)
  if (length(samp_col) == 0) samp_col <- "Sample Name"
  samples <- runsheet[[samp_col[1]]]

  # control tokens (checked first): ground control, vivarium, sham, non-irradiated,
  #   hindlimb loaded, normal loading, 1G on earth, 1G by centrifugation, 1G in clinostat
  control_tokens <- c("ground control","vivarium control","sham-irradiated","non-irradiated",
                      "hindlimb loaded control","normal loading control","baseline","bsl",
                      "1g on earth","1g by centrifugation","1g in 3d clinostat","1g in clinostat",
                      "control")
  # treated tokens: space flight, hindlimb unloading, ionizing radiation, uG, microgravity, etc.
  treated_tokens <- c("space flight","spaceflight","hindlimb unloading","hindlimb unloaded",
                      "ionizing radiation","fe-56","proton","o-16","si-28","mixed radiation",
                      "ug by 3d clinostat","ug in clinostat","ug","microgravity","flight",
                      "irradiated","unloaded")

  # IMPORTANT: some studies (e.g. OSD-811) have ALL samples as "Space Flight" in one
  # factor column, but split by "Altered Gravity" (uG vs 1G centrifugation). So we must
  # check ALL factor columns for each sample and prefer the most specific control/treated
  # designation. Strategy: for each sample, collect all factor values across columns,
  # then assign control if ANY column has a control token, else treated if ANY has a
  # treated token. This lets "1G by centrifugation" override "Space Flight".
  conditions <- rep(NA_character_, length(samples))
  for (i in seq_along(samples)) {
    all_vals <- tolower(unlist(runsheet[i, fcols, drop = TRUE]))
    has_control <- any(sapply(control_tokens, function(t) any(grepl(t, all_vals, fixed = TRUE))))
    has_treated <- any(sapply(treated_tokens, function(t) any(grepl(t, all_vals, fixed = TRUE))))
    if (has_control && !has_treated) {
      conditions[i] <- "control"
    } else if (has_treated && !has_control) {
      conditions[i] <- "treated"
    } else if (has_control && has_treated) {
      # both present: prefer control if a strong control token (ground/vivarium/sham/non-irr/1G)
      # is found, else treated. This handles OSD-811 (Space Flight + 1G centrifugation = control)
      strong_control <- c("ground control","vivarium control","sham-irradiated","non-irradiated",
                          "1g on earth","1g by centrifugation","1g in 3d clinostat","1g in clinostat",
                          "hindlimb loaded control","normal loading control")
      if (any(sapply(strong_control, function(t) any(grepl(t, all_vals, fixed = TRUE))))) {
        conditions[i] <- "control"
      } else {
        conditions[i] <- "treated"
      }
    }
  }
  data.frame(sample_name = samples, condition = conditions, stringsAsFactors = FALSE)
}

#' Assign condition from sample name (fallback when no runsheet).
#' Sample name patterns: FLT=flight/treated, GC/CC/BSL/VIV=control, uG=microgravity/treated
assign_condition_from_name <- function(sample_names) {
  conds <- rep(NA_character_, length(sample_names))
  sn <- tolower(sample_names)
  # treated: FLT (flight), uG (microgravity), Unloaded, irradiated
  conds[grepl("_flt_|_flt\\.|flt_c|_ug_|_ug\\.|unloaded|irradiated|_irr_|fe-56|proton|si-28|o-16", sn)] <- "treated"
  # control: GC, CC, BSL, VIV, loaded, sham, non-irr
  conds[grepl("_gc_|_gc\\.|gc_c|_cc_|_cc\\.|cc_c|_bsl_|_bsl\\.|bsl_c|_viv_|_viv\\.|viv_c|loaded|sham|non_irr|nonirr", sn)] <- "control"
  data.frame(sample_name = sample_names, condition = conds, stringsAsFactors = FALSE)
}

#' Load one study: counts matrix + condition assignment.
#' @return list(counts = matrix, conditions = data.frame, gene_ids = character, species = character)
load_study <- function(acc, species) {
  study_dir <- file.path("data", "raw", acc)
  counts_file <- find_counts_file(study_dir)
  if (is.na(counts_file)) {
    cat(sprintf("  [%s] no counts file found, skipping\n", acc))
    return(NULL)
  }
  counts <- read.csv(counts_file, check.names = FALSE, row.names = 1)
  gene_ids <- rownames(counts)
  samples <- colnames(counts)

  # assign conditions
  runsheet_file <- find_runsheet(study_dir)
  conds <- NULL
  if (!is.na(runsheet_file)) {
    rs <- read.csv(runsheet_file, check.names = FALSE, stringsAsFactors = FALSE)
    conds <- assign_condition_from_runsheet(rs)
    # match runsheet samples to counts columns
    if (!is.null(conds)) {
      m <- match(samples, conds$sample_name)
      conditions <- conds$condition[m]
      # if any NA, fall back to name-based for those
      if (any(is.na(conditions))) {
        name_conds <- assign_condition_from_name(samples)
        conditions[is.na(conditions)] <- name_conds$condition[is.na(conditions)]
      }
      conds <- data.frame(sample_name = samples, condition = conditions, stringsAsFactors = FALSE)
    }
  }
  if (is.null(conds)) {
    conds <- assign_condition_from_name(samples)
  }

  n_ctrl <- sum(conds$condition == "control", na.rm = TRUE)
  n_trt  <- sum(conds$condition == "treated", na.rm = TRUE)
  n_na   <- sum(is.na(conds$condition))
  cat(sprintf("  [%s] %d samples: %d control, %d treated, %d unassigned\n",
              acc, length(samples), n_ctrl, n_trt, n_na))

  if (n_ctrl < MIN_SAMPLES_PER_GROUP || n_trt < MIN_SAMPLES_PER_GROUP) {
    cat(sprintf("  [%s] WARNING: insufficient samples in one group (ctrl=%d, trt=%d)\n",
                acc, n_ctrl, n_trt))
  }

  list(counts = counts, conditions = conds, gene_ids = gene_ids,
       species = species, n_control = n_ctrl, n_treated = n_trt)
}

# ---- C. Ortholog mapping for rodent studies ----

#' Map Ensembl mouse gene IDs -> human HGNC symbols.
#' @param gene_ids character vector of ENSMUSG... IDs
#' @return data.frame(ensembl_mouse, hgnc_symbol)
map_mouse_ensembl_to_human <- function(gene_ids) {
  cat(sprintf("[orthologs] mapping %d mouse Ensembl IDs -> human HGNC\n", length(gene_ids)))
  cache_file <- file.path("data", "processed", "mouse_ensembl_to_human_cache.tsv")
  if (file.exists(cache_file)) {
    cache <- read.table(cache_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (all(gene_ids %in% cache$ensembl_mouse)) {
      cat("[orthologs] using cached mapping\n")
      return(cache[match(gene_ids, cache$ensembl_mouse), ])
    }
  }
  # use biomaRt: mouse mart with external_gene_name (mgi) + getLDS to human hgnc
  # try multiple mirrors (Ensembl can redirect to status page during outages)
  mart <- tryCatch(
    biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
    error = function(e) NULL
  )
  if (is.null(mart)) {
    for (host in c("https://www.ensembl.org","https://useast.ensembl.org","https://asia.ensembl.org")) {
      cat(sprintf("[orthologs] trying mirror: %s\n", host))
      mart <- tryCatch(biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host),
                       error = function(e) NULL)
      if (!is.null(mart)) break
      Sys.sleep(2)
    }
  }
  if (is.null(mart)) {
    cat("[orthologs] all Ensembl mirrors failed; using gene-ID fallback\n")
    return(data.frame(ensembl_mouse = gene_ids, hgnc_symbol = gene_ids, stringsAsFactors = FALSE))
  }
  # first: ensembl_mouse -> mgi_symbol (with mirror retry)
  res1 <- NULL
  for (attempt in 1:3) {
    res1 <- tryCatch({
      biomaRt::getBM(attributes = c("ensembl_gene_id","mgi_symbol"),
                     filters = "ensembl_gene_id",
                     values = gene_ids, mart = mart)
    }, error = function(e) {
      cat(sprintf("[orthologs] getBM attempt %d failed: %s\n", attempt, conditionMessage(e)))
      # retry with a mirror
      for (host in c("https://www.ensembl.org","https://useast.ensembl.org","https://asia.ensembl.org")) {
        cat(sprintf("[orthologs] retrying getBM via mirror: %s\n", host))
        mart <<- tryCatch(biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host),
                          error = function(e2) NULL)
        if (!is.null(mart)) {
          res <- tryCatch(biomaRt::getBM(attributes = c("ensembl_gene_id","mgi_symbol"),
                                          filters = "ensembl_gene_id", values = gene_ids, mart = mart),
                          error = function(e2) NULL)
          if (!is.null(res)) return(res)
        }
      }
      NULL
    })
    if (!is.null(res1) && nrow(res1) > 0) break
    Sys.sleep(3)
  }
  if (is.null(res1) || nrow(res1) == 0) {
    cat("[orthologs] all getBM attempts failed; using Ensembl-ID fallback (no symbol mapping)\n")
    # fallback: return Ensembl IDs as-is (will yield 0 DE genes but won't crash)
    out <- data.frame(ensembl_mouse = gene_ids, hgnc_symbol = NA_character_, stringsAsFactors = FALSE)
    return(out)
  }
  # drop blank mgi symbols
  res1 <- res1[res1$mgi_symbol != "" & !is.na(res1$mgi_symbol), ]
  # second: mgi_symbol -> hgnc_symbol via getLDS
  martH <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  res2 <- tryCatch({
    biomaRt::getLDS(attributes = "mgi_symbol", filters = "mgi_symbol",
                    values = res1$mgi_symbol, mart = mart,
                    attributesL = "hgnc_symbol", martL = martH)
  }, error = function(e) {
    cat("[orthologs] getLDS failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  if (is.null(res2) || nrow(res2) == 0) {
    # fallback: just use mgi symbols as human (many are same)
    out <- data.frame(ensembl_mouse = res1$ensembl_gene_id,
                      hgnc_symbol = toupper(res1$mgi_symbol), stringsAsFactors = FALSE)
  } else {
    colnames(res2) <- c("mgi_symbol","hgnc_symbol")
    merged <- merge(res1, res2, by = "mgi_symbol")
    merged <- merged[merged$hgnc_symbol != "" & !is.na(merged$hgnc_symbol), ]
    # dedupe (keep first)
    merged <- merged[!duplicated(merged$ensembl_gene_id), ]
    out <- data.frame(ensembl_mouse = merged$ensembl_gene_id, hgnc_symbol = merged$hgnc_symbol,
                      stringsAsFactors = FALSE)
  }
  dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
  write.table(out, cache_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("[orthologs] %d mouse Ensembl -> %d human HGNC (cached)\n", length(gene_ids), nrow(out)))
  out
}

#' Map Ensembl human gene IDs -> HGNC symbols (for human studies).
map_human_ensembl_to_hgnc <- function(gene_ids) {
  cat(sprintf("[hgnc] mapping %d human Ensembl IDs -> HGNC\n", length(gene_ids)))
  cache_file <- file.path("data", "processed", "human_ensembl_to_hgnc_cache.tsv")
  if (file.exists(cache_file)) {
    cache <- read.table(cache_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (all(gene_ids %in% cache$ensembl_human)) {
      cat("[hgnc] using cached mapping\n")
      return(cache[match(gene_ids, cache$ensembl_human), ])
    }
  }
  mart <- tryCatch(
    biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
    error = function(e) NULL
  )
  if (is.null(mart)) {
    for (host in c("https://www.ensembl.org","https://useast.ensembl.org","https://asia.ensembl.org")) {
      cat(sprintf("[hgnc] trying mirror: %s\n", host))
      mart <- tryCatch(biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host),
                       error = function(e) NULL)
      if (!is.null(mart)) break
      Sys.sleep(2)
    }
  }
  if (is.null(mart)) {
    cat("[hgnc] all Ensembl mirrors failed; using gene-ID fallback\n")
    return(data.frame(ensembl_human = gene_ids, hgnc_symbol = gene_ids, stringsAsFactors = FALSE))
  }
  res <- tryCatch({
    biomaRt::getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                   filters = "ensembl_gene_id", values = gene_ids, mart = mart)
  }, error = function(e) {
    cat("[hgnc] getBM failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  if (is.null(res) || nrow(res) == 0) {
    # fallback: strip version, use as-is
    return(data.frame(ensembl_human = gene_ids, hgnc_symbol = gene_ids, stringsAsFactors = FALSE))
  }
  res <- res[res$hgnc_symbol != "" & !is.na(res$hgnc_symbol), ]
  res <- res[!duplicated(res$ensembl_gene_id), ]
  colnames(res) <- c("ensembl_human","hgnc_symbol")
  dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
  write.table(res, cache_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("[hgnc] %d mappings (cached)\n", nrow(res)))
  res
}

# ---- D. DESeq2 per study ----

run_deseq2 <- function(counts, conditions, gene_symbols) {
  # keep only control + treated samples
  keep <- conditions$condition %in% c("control","treated")
  counts <- counts[, keep, drop = FALSE]
  conditions <- conditions[keep, , drop = FALSE]
  if (ncol(counts) < 4 || sum(conditions$condition=="control") < MIN_SAMPLES_PER_GROUP ||
      sum(conditions$condition=="treated") < MIN_SAMPLES_PER_GROUP) {
    return(NULL)
  }
  # drop genes with no symbol mapping
  valid <- !is.na(gene_symbols) & gene_symbols != ""
  if (sum(valid) < 100) {
    cat(sprintf("  DESeq2 skipped: only %d genes mapped to symbols\n", sum(valid)))
    return(NULL)
  }
  counts <- counts[valid, , drop = FALSE]
  gene_symbols <- gene_symbols[valid]
  # aggregate counts by gene symbol (sum duplicates) -> unique rows
  counts <- rowsum(as.matrix(counts), group = gene_symbols)
  # round to integers (DESeq2 requirement)
  counts <- round(counts)
  # filter low-count genes
  keep_genes <- rowSums(counts >= 10) >= max(2, ncol(counts) * 0.1)
  counts <- counts[keep_genes, , drop = FALSE]

  coldata <- data.frame(condition = factor(conditions$condition, levels = c("control","treated")),
                        row.names = conditions$sample_name)
  # ensure coldata rows match counts columns
  coldata <- coldata[colnames(counts), , drop = FALSE]

  tryCatch({
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    res <- DESeq2::results(dds, contrast = c("condition","treated","control"))
    res <- as.data.frame(res)
    res$gene <- rownames(res)
    res
  }, error = function(e) {
    cat(sprintf("  DESeq2 failed: %s\n", conditionMessage(e)))
    NULL
  })
}

# ---- E. metafor random-effects meta-analysis ----

run_meta_analysis <- function(de_list) {
  cat("\n[meta] running random-effects meta-analysis per gene...\n")
  # collect per-gene: log2FC + SE across studies
  # de_list is a list of data.frames (one per study), each with gene, log2FoldChange, lfcSE
  all_genes <- unique(unlist(lapply(de_list, function(d) d$gene)))
  cat(sprintf("[meta] %d unique genes across %d studies\n", length(all_genes), length(de_list)))

  # build a long table: gene, study, lfc, se
  long <- do.call(rbind, lapply(seq_along(de_list), function(i) {
    d <- de_list[[i]]
    data.frame(gene = d$gene, study = names(de_list)[i],
               lfc = d$log2FoldChange, se = d$lfcSE, stringsAsFactors = FALSE)
  }))
  long <- long[!is.na(long$lfc) & !is.na(long$se) & long$se > 0, ]

  # per-gene random-effects meta (metafor::rma)
  # this is the slow step; parallelize with parallel::mclapply
  genes_with_data <- unique(long$gene)
  cat(sprintf("[meta] running rma on %d genes (parallel)...\n", length(genes_with_data)))

  ncores <- max(1, parallel::detectCores() - 1)
  results <- parallel::mclapply(genes_with_data, function(g) {
    sub <- long[long$gene == g, ]
    if (nrow(sub) < 2) {
      # single study: use the single estimate
      return(data.frame(gene = g, n_studies = nrow(sub), mean_lfc = sub$lfc[1],
                        se = sub$se[1], z = sub$lfc[1]/sub$se[1],
                        p = 2*pnorm(-abs(sub$lfc[1]/sub$se[1])),
                        direction_concordance = NA_real_, n_species = NA_integer_))
    }
    tryCatch({
      fit <- metafor::rma(yi = sub$lfc, sei = sub$se, method = "REML", test = "z")
      # direction concordance: fraction of studies with same sign as pooled mean
      dir_conc <- mean(sign(sub$lfc) == sign(fit$b), na.rm = TRUE)
      data.frame(gene = g, n_studies = nrow(sub), mean_lfc = as.numeric(fit$b),
                 se = fit$se, z = as.numeric(fit$zval), p = as.numeric(fit$pval),
                 direction_concordance = dir_conc, n_species = NA_integer_)
    }, error = function(e) {
      # fallback: Stouffer's
      z <- sub$lfc / sub$se
      combined_z <- sum(z) / sqrt(length(z))
      data.frame(gene = g, n_studies = nrow(sub), mean_lfc = mean(sub$lfc),
                 se = sqrt(mean(sub$se^2)), z = combined_z,
                 p = 2*pnorm(-abs(combined_z)),
                 direction_concordance = mean(sign(sub$lfc) == sign(combined_z)),
                 n_species = NA_integer_)
    })
  }, mc.cores = ncores)

  meta <- do.call(rbind, results)
  meta$padj <- p.adjust(meta$p, method = "BH")
  # add n_species (how many species contributed)
  study_species <- sapply(names(de_list), function(s) {
    # determine species from study accession (rodent panel has species in T1)
    if (s %in% names(RODENT_ACC)) "mouse" else "human"
  })
  for (i in seq_len(nrow(meta))) {
    g <- meta$gene[i]
    contributing <- unique(long$study[long$gene == g])
    meta$n_species[i] <- length(unique(study_species[contributing]))
  }
  meta <- meta[order(meta$padj), ]
  rownames(meta) <- NULL
  cat(sprintf("[meta] done: %d genes, %d significant (padj<0.05)\n",
              nrow(meta), sum(meta$padj < 0.05)))
  meta
}

# ---- main ----

main <- function() {
  root <- getwd()
  raw_dir <- file.path(root, "data", "raw")
  proc_dir <- file.path(root, "data", "processed")
  tab_dir <- file.path(root, "results", "tables")
  fig_dir <- file.path(root, "results", "figures")
  dir.create(proc_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

  # load panel
  panel <- read.csv(file.path(tab_dir, "T1_dataset_panel.csv"), stringsAsFactors = FALSE)
  RODENT_ACC <<- panel$accession[panel$species == "mouse"]
  HUMAN_ACC  <<- panel$accession[panel$species == "human"]

  # ---- A/B/C/D. Load + DESeq2 per study ----
  cat("\n=== Per-study DESeq2 ===\n")
  de_list <- list()
  sample_meta <- list()
  for (i in seq_len(nrow(panel))) {
    acc <- panel$accession[i]
    species <- panel$species[i]
    cat(sprintf("\n[%s] (%s, %s)\n", acc, species, panel$tissue[i]))
    study <- load_study(acc, species)
    if (is.null(study)) next

    # gene ID mapping
    if (species == "mouse") {
      # mouse Ensembl -> human HGNC
      if (startsWith(study$gene_ids[1], "ENSMUSG")) {
        map <- map_mouse_ensembl_to_human(study$gene_ids)
        gene_symbols <- map$hgnc_symbol[match(study$gene_ids, map$ensembl_mouse)]
      } else {
        # already symbols (mouse) -> map to human via orthologs
        map <- map_orthologs(study$gene_ids)
        gene_symbols <- map$human_symbol[match(study$gene_ids, map$mouse_symbol)]
      }
    } else {
      # human Ensembl -> HGNC
      if (startsWith(study$gene_ids[1], "ENSG")) {
        map <- map_human_ensembl_to_hgnc(study$gene_ids)
        gene_symbols <- map$hgnc_symbol[match(study$gene_ids, map$ensembl_human)]
      } else {
        gene_symbols <- study$gene_ids  # already HGNC
      }
    }
    cat(sprintf("  [%s] %d/%d genes mapped to HGNC\n", acc, sum(!is.na(gene_symbols)), length(gene_symbols)))

    # DESeq2
    de <- run_deseq2(study$counts, study$conditions, gene_symbols)
    if (!is.null(de) && nrow(de) > 0) {
      de_list[[acc]] <- de
      cat(sprintf("  [%s] DESeq2: %d genes, %d significant (padj<0.05, |lFC|>1)\n",
                  acc, nrow(de), sum(de$padj < PADJ & abs(de$log2FoldChange) > LFC, na.rm = TRUE)))
    }

    # sample metadata
    sm <- study$conditions
    sm$study <- acc
    sm$tissue <- panel$tissue[i]
    sm$species <- species
    sm$factor <- panel$factor[i]
    sm$system <- panel$system[i]
    sample_meta[[length(sample_meta)+1]] <- sm
  }

  # save per-study DE
  de_all <- do.call(rbind, lapply(names(de_list), function(s) {
    d <- de_list[[s]]; d$study <- s; d
  }))
  write.csv(de_all, file.path(tab_dir, "T2_per_study_DE.csv"), row.names = FALSE)
  cat(sprintf("\nWrote T2_per_study_DE.csv (%d rows from %d studies)\n", nrow(de_all), length(de_list)))

  # save sample metadata
  sm_all <- do.call(rbind, sample_meta)
  write.csv(sm_all, file.path(proc_dir, "sample_metadata.csv"), row.names = FALSE)
  cat(sprintf("Wrote sample_metadata.csv (%d samples)\n", nrow(sm_all)))

  # ---- E. metafor random-effects meta-analysis ----
  meta <- run_meta_analysis(de_list)
  write.table(meta, file.path(tab_dir, "T3_meta_signature.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("Wrote T3_meta_signature.tsv (%d genes)\n", nrow(meta)))

  # ---- F. oncogene union ----
  cat("\n=== Building oncogene union ===\n")
  onc <- build_oncogene_union()
  write.csv(onc$log, file.path(tab_dir, "T4a_oncogene_sources.csv"), row.names = FALSE)

  # ---- G. intersect meta-signature with oncogene union ----
  sig_genes <- meta$gene[meta$padj < PADJ]
  cat(sprintf("\n[intersection] %d significant meta-signature genes (padj<0.05)\n", length(sig_genes)))
  intersection <- onc$genes[onc$genes$gene %in% sig_genes, ]
  # attach meta-signature stats
  intersection <- merge(intersection, meta[, c("gene","mean_lfc","padj","n_studies","n_species","direction_concordance")],
                        by = "gene", all.x = TRUE)
  intersection <- intersection[order(intersection$padj), ]
  write.csv(intersection, file.path(tab_dir, "T4_oncogene_intersection.csv"), row.names = FALSE)
  cat(sprintf("Wrote T4_oncogene_intersection.csv (%d genes)\n", nrow(intersection)))

  # ---- H. validation vs S2 ----
  cat("\n=== Validation vs Python S2 ===\n")
  s2_path <- file.path(root, "data", "external", "S2_cross_study_meta_signature.csv")
  if (!file.exists(s2_path)) {
    # try the discovery repo path
    s2_path <- "/workspace/astronaut-oncogene-biomarkers/results/tables/S2_cross_study_meta_signature.csv"
  }
  if (file.exists(s2_path)) {
    val <- validate_vs_s2(meta, s2_path)
    write.csv(val$metrics, file.path(tab_dir, "T5_validation_concordance.csv"), row.names = FALSE)
    write.csv(val$joined, file.path(tab_dir, "T5b_R_vs_S2_joined.csv"), row.names = FALSE)
    cat("Wrote T5_validation_concordance.csv + T5b_R_vs_S2_joined.csv\n")
  } else {
    cat("WARNING: S2 not found, skipping validation\n")
    val <- NULL
  }

  # ---- Figures ----
  cat("\n=== Generating figures ===\n")
  make_figures(meta, intersection, val, fig_dir)

  cat("\nNext: Rscript scripts/03_drug_screening.R\n")
}

# ---- Figures ----

make_figures <- function(meta, intersection, val, fig_dir) {
  # common theme
  th <- ggplot2::theme_minimal(base_family = "Liberation Sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 12, face = "bold"))

  # F3: meta-signature volcano + top-gene heatmap
  meta$sig <- ifelse(meta$padj < PADJ & abs(meta$mean_lfc) > LFC, "significant", "ns")
  p_volcano <- ggplot2::ggplot(meta, ggplot2::aes(x = mean_lfc, y = -log10(padj), color = sig)) +
    ggplot2::geom_point(alpha = 0.4, size = 0.8) +
    ggplot2::scale_color_manual(values = c("significant" = "#c0392b", "ns" = "grey70")) +
    ggplot2::geom_vline(xintercept = c(-LFC, LFC), linetype = "dashed", color = "grey50") +
    ggplot2::geom_hline(yintercept = -log10(PADJ), linetype = "dashed", color = "grey50") +
    ggplot2::labs(title = "R random-effects meta-signature (human + rodent)",
                  x = "mean log2 fold change (treated vs control)",
                  y = "-log10 adjusted p") +
    ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) + th
  save_fig(p_volcano, file.path(fig_dir, "F3_meta_signature_volcano"))

  # F4: oncogene intersection barplot by source
  if (nrow(intersection) > 0) {
    src_counts <- as.data.frame(table(intersection$source), stringsAsFactors = FALSE)
    colnames(src_counts) <- c("source","n_genes")
    src_counts <- src_counts[order(-src_counts$n_genes), ]
    p_onc <- ggplot2::ggplot(src_counts, ggplot2::aes(x = reorder(source, n_genes), y = n_genes)) +
      ggplot2::geom_col(fill = "#0279EE") +
      ggplot2::coord_flip() +
      ggplot2::labs(title = "Oncogenic biomarker intersection by source",
                    x = NULL, y = "n genes") + th
    save_fig(p_onc, file.path(fig_dir, "F4_oncogene_intersection"))
  }

  # F5: validation concordance scatter
  if (!is.null(val) && nrow(val$joined) > 0) {
    j <- val$joined
    p_val <- ggplot2::ggplot(j, ggplot2::aes(x = mean_log2fc, y = mean_lfc)) +
      ggplot2::geom_point(alpha = 0.3, size = 0.6, color = "#0279EE") +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
      ggplot2::geom_hline(yintercept = 0, color = "grey80") +
      ggplot2::geom_vline(xintercept = 0, color = "grey80") +
      ggplot2::labs(title = "R vs Python S2 concordance (mean log2FC)",
                    x = "Python S2 mean log2FC", y = "R meta mean log2FC") + th
    save_fig(p_val, file.path(fig_dir, "F5_validation_concordance"))
  }
}

save_fig <- function(plot, base_path) {
  svg_path <- paste0(base_path, ".svg")
  png_path <- paste0(base_path, ".png")
  svglite::svglite(svg_path, width = 8, height = 6)
  print(plot)
  dev.off()
  ggplot2::ggsave(png_path, plot, width = 8, height = 6, dpi = 150, bg = "white")
  cat(sprintf("  saved %s + %s\n", basename(svg_path), basename(png_path)))
}

main()

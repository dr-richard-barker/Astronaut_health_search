#' NASA OSDR (Open Science Data Repository) API client for R.
#'
#' Mirrors the Python osdr_client.py in the discovery repo. Fetches only
#' *processed* lightweight products (count matrices, DGE tables, ISA metadata);
#' skips raw FASTQ/BAM. Idempotent / resume-safe (skips already-downloaded files).
#'
#' Endpoints (verified 2026-07-19):
#'   Search : https://osdr.nasa.gov/osdr/data/search?term=&type=cgene&size=N&ffield=organism&fvalue=...
#'   Files  : https://osdr.nasa.gov/osdr/data/osd/files/{numeric_id}
#'            -> JSON: studies[OSD-xxx].study_files[] with file_name, category,
#'               subcategory, remote_url, file_size
#'   Download: https://osdr.nasa.gov + remote_url

BASE <- "https://osdr.nasa.gov"
SEARCH <- paste0(BASE, "/osdr/data/search")
FILES  <- paste0(BASE, "/osdr/data/osd/files/{num}")

# product types to keep (drop everything else, esp. raw sequence)
KEEP_NAME <- "(norm.*counts|unnormalized.*counts|counts.*csv|differential.*expression|dge|_DGE_|contrasts|ISA\\.zip|metadata|runsheet|samples?\\.txt)"
SKIP_NAME <- "\\.(fastq|fq|bam|bai|sam|cram|bigwig|bw)(\\.gz)?$"
MAX_MB <- 60  # safety cap per processed file

#' Extract numeric id from accession (e.g. "OSD-940" -> "940")
osd_num <- function(accession) {
  m <- regmatches(accession, regexpr("\\d+", accession))
  if (length(m) == 0) stop("no numeric id in ", sQuote(accession))
  m
}

#' List all study files for an OSD accession.
#' @return data.frame with columns: file_name, category, subcategory,
#'         remote_url, file_size, study
list_files <- function(accession, session = NULL) {
  url <- sub("{num}", osd_num(accession), FILES, fixed = TRUE)
  if (is.null(session)) {
    r <- httr::GET(url, httr::timeout(60),
                   httr::user_agent("astronaut-opposite-forcing/0.1 (R)"))
  } else {
    r <- httr::GET(url, session, httr::timeout(60))
  }
  httr::stop_for_status(r)
  data <- httr::content(r, as = "parsed", type = "application/json")
  out <- list()
  studies <- data[["studies"]] %||% list()
  for (sid in names(studies)) {
    sfiles <- studies[[sid]][["study_files"]] %||% list()
    for (f in sfiles) {
      f[["study"]] <- sid
      out[[length(out) + 1]] <- f
    }
  }
  if (length(out) == 0) return(data.frame())
  do.call(rbind, lapply(out, as.data.frame, stringsAsFactors = FALSE))
}

#' Filter a file listing to the lightweight processed products we want.
select_processed <- function(files) {
  if (nrow(files) == 0) return(files)
  size_mb <- (as.numeric(files$file_size) %||% 0) / 1e6
  keep <- grepl(KEEP_NAME, files$file_name, ignore.case = TRUE, perl = TRUE)
  skip <- grepl(SKIP_NAME, files$file_name, ignore.case = TRUE, perl = TRUE)
  files[keep & !skip & size_mb <= MAX_MB, , drop = FALSE]
}

#' Download one file (resume-safe: skip if already present and non-empty).
download <- function(f, dest_dir, session = NULL) {
  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
  remote <- f[["remote_url"]]
  url <- if (startsWith(remote, "/")) paste0(BASE, remote) else remote
  dest <- file.path(dest_dir, f[["file_name"]])
  if (file.exists(dest) && file.info(dest)$size > 0) return(dest)
  if (is.null(session)) {
    r <- httr::GET(url, httr::write_disk(paste0(dest, ".part"), overwrite = TRUE),
                   httr::timeout(300),
                   httr::user_agent("astronaut-opposite-forcing/0.1 (R)"))
  } else {
    r <- httr::GET(url, session, httr::write_disk(paste0(dest, ".part"), overwrite = TRUE),
                   httr::timeout(300))
  }
  httr::stop_for_status(r)
  # use file.copy + file.remove (file.rename fails on some FUSE mounts)
  file.copy(paste0(dest, ".part"), dest, overwrite = TRUE)
  file.remove(paste0(dest, ".part"))
  Sys.sleep(0.3)  # be polite
  dest
}

#' Search OSDR for studies by organism.
#' @param organism e.g. "Homo sapiens", "Mus musculus"
#' @param size max results (default 250)
#' @return data.frame: accession, study_title, assay_tech, organism, factor_type, factor_value
search_osdr <- function(organism, size = 250, session = NULL) {
  url <- paste0(SEARCH, "?term=&type=cgene&size=", size,
                "&ffield=organism&fvalue=", utils::URLencode(organism))
  if (is.null(session)) {
    r <- httr::GET(url, httr::timeout(60),
                   httr::user_agent("astronaut-opposite-forcing/0.1 (R)"))
  } else {
    r <- httr::GET(url, session, httr::timeout(60))
  }
  httr::stop_for_status(r)
  data <- httr::content(r, as = "parsed", type = "application/json")
  hits <- data[["hits"]][["hits"]]
  if (length(hits) == 0) return(data.frame())
  rows <- lapply(hits, function(h) {
    s <- h[["_source"]]
    data.frame(
      accession       = s[["Accession"]] %||% NA_character_,
      study_title     = s[["Study Title"]] %||% NA_character_,
      assay_tech      = s[["Study Assay Technology Type"]] %||% NA_character_,
      organism        = s[["organism"]] %||% NA_character_,
      factor_type     = paste(s[["Study Factor Type"]] %||% "", collapse = "; "),
      factor_value    = paste(s[["Factor Value"]] %||% "", collapse = "; "),
      material_type   = s[["Material Type"]] %||% NA_character_,
      study_id        = s[["Study Identifier"]] %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# null-coalescing helper
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

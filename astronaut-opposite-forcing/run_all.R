#!/usr/bin/env Rscript
# ==============================================================================
# run_all.R
# astronaut-opposite-forcing pipeline -- one-command reproduction driver
#
# Executes the full pipeline in order:
#   01_OSD_query.R              -- download OSDR human + rodent RNA-seq
#   02_biomarker_identification.R -- DESeq2 + metafor meta + oncogene union
#   03_drug_screening.R         -- LINCS tau-analog reversal + Broad + PrimeKG
#   04_enrichment.R             -- PrimeKG + Open Targets disease enrichment
#   05_nutraceutical_flags.R    -- flag nutraceuticals among top reversal hits
#   06_figures_F1_F2.R          -- dataset panel + QC PCA figures
#
# Usage:
#   Rscript run_all.R            # run all steps
#   Rscript run_all.R --from 03  # resume from step 03
#   Rscript run_all.R --only 04  # run only step 04
#
# Prerequisites:
#   - R 4.4+ with packages listed in renv.lock
#   - Internet access for OSDR API, biomaRt, Open Targets
#   - Datalake mounted at /mnt/datalake (LINCS, Broad, PrimeKG)
#
# Richard Barker -- 2026
# ==============================================================================

repo_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(sub("^--file=", "", file_arg[1]))
  else "."
}, error = function(e) ".")

scripts <- c(
  "01_OSD_query.R",
  "02_biomarker_identification.R",
  "03_drug_screening.R",
  "04_enrichment.R",
  "05_nutraceutical_flags.R",
  "06_figures_F1_F2.R",
  "07_gene_level_reversal.R"
)

# Parse CLI flags
cli <- commandArgs(trailingOnly = TRUE)
from_step <- 1
only_step <- NA
if ("--from" %in% cli) {
  idx <- which(cli == "--from")
  from_step <- as.integer(cli[idx + 1])
}
if ("--only" %in% cli) {
  idx <- which(cli == "--only")
  only_step <- as.integer(cli[idx + 1])
}

cat("============================================================\n")
cat(" astronaut-opposite-forcing pipeline\n")
cat(" R version:", R.version.string, "\n")
cat(" Repo:", repo_dir, "\n")
cat(" Steps:", paste(scripts, collapse = " -> "), "\n")
cat("============================================================\n\n")

run_script <- function(script) {
  path <- file.path(repo_dir, "scripts", script)
  if (!file.exists(path)) {
    stop("Script not found: ", path)
  }
  cat(sprintf("\n>>> [%s] %s\n", format(Sys.time(), "%H:%M:%S"), script))
  cat(strrep("-", 60), "\n")
  t0 <- Sys.time()
  exit_code <- system2("Rscript", args = path, wait = TRUE)
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  if (exit_code != 0) {
    cat(sprintf("\n[FAIL] %s exited with code %d after %.1f min\n",
                script, exit_code, as.numeric(elapsed)))
    stop("Pipeline halted at ", script)
  }
  cat(sprintf("\n[OK] %s completed in %.1f min\n",
              script, as.numeric(elapsed)))
  invisible(exit_code)
}

# Determine which steps to run
if (!is.na(only_step)) {
  steps_to_run <- only_step
} else {
  steps_to_run <- from_step:length(scripts)
}

for (i in steps_to_run) {
  if (i < 1 || i > length(scripts)) next
  run_script(scripts[i])
}

cat("\n============================================================\n")
cat(" Pipeline complete.\n")
cat(" Figures: results/figures/F*.svg + F*.png\n")
cat(" Tables:  results/tables/T*.csv + T*.tsv\n")
cat("============================================================\n")

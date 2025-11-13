#!/usr/bin/env Rscript

## STEP 0 – Initialise environment
# Start the script quietly so it can be run from the command line without
# printing library attachment messages.
suppressPackageStartupMessages({
  library(utils)
})

## STEP 1 – Define helper utilities
# Utility helpers are defined once and reused by the main routine.
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- args[startsWith(args, file_arg)]
  if (length(script_path) == 0) {
    return(getwd())
  }
  normalizePath(dirname(sub(file_arg, "", script_path[1])), winslash = "/")
}

check_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    stop(sprintf("Missing required packages: %s", paste(missing, collapse = ", ")))
  }
}

## STEP 2 – Run the orchestration routine
# The main routine aligns working directories, validates dependencies, and
# invokes each script in sequence.
main <- function() {
  ## STEP 2.1 – Align working directory with repository root
  repo_dir <- get_script_dir()
  setwd(repo_dir)

  ## STEP 2.2 – Check required R packages
  required_pkgs <- c("apollo", "data.table", "dplyr", "readr", "tidyr")
  check_packages(required_pkgs)

  ## STEP 2.3 – Ensure data directories exist and seed raw input
  data_dir <- file.path(repo_dir, "DATA")
  processed_dir <- file.path(data_dir, "processed")
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
  if (!dir.exists(processed_dir)) dir.create(processed_dir, recursive = TRUE)

  raw_source <- file.path(repo_dir, "apollo_modeChoiceData.csv")
  raw_target <- file.path(data_dir, "apollo_modeChoiceData.csv")
  if (!file.exists(raw_source)) {
    stop("Unable to locate 'apollo_modeChoiceData.csv' in the repository root.")
  }
  if (!file.exists(raw_target)) {
    file.copy(raw_source, raw_target, overwrite = FALSE)
  }

  ## STEP 2.4 – Run preprocessing script to build the long-format data
  message("==> Running data processing (R00_Data Processing.R)")
  source("R00_Data Processing.R", local = FALSE)

  ## STEP 2.5 – Estimate the independent MMNL model
  message("==> Running MMNL estimation (R_02_MMNL Independent.R)")
  autorun_old <- getOption("mmnl.skip.autorun")
  on.exit(options(mmnl.skip.autorun = autorun_old), add = TRUE)
  options(mmnl.skip.autorun = TRUE)
  source("R_02_MMNL Independent.R", local = FALSE)
  options(mmnl.skip.autorun = autorun_old)
  results <- mmnl_model(
    processed_path = file.path("DATA", "processed", "modechoice_long.csv"),
    output_dir = "output"
  )

  ## STEP 2.6 – Report completion status to the console
  message(sprintf("==> Estimation finished. Summary written to %s", results$summary_path))
}

## STEP 3 – Handle runtime errors gracefully
# Wrap the main routine in a guard that returns a non-zero exit code when the
# pipeline fails.
tryCatch(
  main(),
  error = function(cond) {
    message("Run failed: ", conditionMessage(cond))
    quit(status = 1L)
  }
)

#!/usr/bin/env Rscript

## STEP 0 – Initialise environment
suppressPackageStartupMessages({
  library(utils)
})

## STEP 1 – Helper utilities ---------------------------------------------------

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
main <- function() {
  ## STEP 2.1 – Align working directory with repository root
  repo_dir <- get_script_dir()
  setwd(repo_dir)

  ## STEP 2.2 – Check required R packages
  required_pkgs <- c("apollo")
  check_packages(required_pkgs)

  ## STEP 2.3 – Copy raw data into DATA/ if needed
  data_dir <- file.path(repo_dir, "DATA")
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

  raw_source <- file.path(repo_dir, "apollo_modeChoiceData.csv")
  raw_target <- file.path(data_dir, "apollo_modeChoiceData.csv")
  if (!file.exists(raw_source)) {
    stop("Unable to locate 'apollo_modeChoiceData.csv' in the repository root.")
  }
  if (!file.exists(raw_target)) {
    file.copy(raw_source, raw_target, overwrite = FALSE)
  }

  ## STEP 2.4 – Run baseline MNL (wide data)
  message("==> Running R_01_MNL Base Line.R")
  source("R_01_MNL Base Line.R", local = FALSE)

  ## STEP 2.5 – Run MMNL independent (wide data)
  message("==> Running R_02_MMNL Independent.R")
  source("R_02_MMNL Independent.R", local = FALSE)

  ## STEP 2.6 – Run MMNL dependent (wide data)
  message("==> Running R_03_MMNL Dependent.R")
  source("R_03_MMNL Dependent.R", local = FALSE)

  message("==> All models have finished running.")
}

## STEP 3 – Execute pipeline
main()


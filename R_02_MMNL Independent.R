##############################################################
# R02_MMNL_Independent.R
# ------------------------------------------------------------
# Project   : Mode Choice Modelling with Apollo
# Researcher: Karina
# Date      : 30 October 2025
#
# Description:
# This script estimates a Mixed Multinomial Logit (MMNL) model
# with independent random coefficients for travel time and cost.
# The code mirrors the baseline MNL structure and extends it
# with Sobol draws for two normally distributed coefficients.
# ------------------------------------------------------------
# Dependencies: apollo, readr, dplyr, tidyr (optional)
##############################################################


### ==========================================================
### STEP 0 – Load packages and initialise Apollo
### ==========================================================

## --- Hard-reset the error handler to the correct form
options(error = function(e) {
  # Print it first so we can see if there's a real error
  print(e)
  stop(e)  # Then rethrow the error
})


library(apollo)
library(readr)
library(dplyr)
library(tidyr)   # optional — only needed if reshape or clean data again

apollo_initialise()

apollo_control <- list(
  modelName  = "MMNL_independent",
  modelDescr = "Mixed Logit with independent normal coefficients",
  indivID    = "ID",
  panelData  = TRUE,
  mixing     = TRUE,
  nCores     = 1
)


### ==========================================================
### STEP 1 – Read WIDE data and locate project paths
### ==========================================================

get_script_path <- function(){
  p <- NULL

  # Prefer the active document path when running inside RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && isTRUE(rstudioapi::isAvailable())) {
    p_try <- try(rstudioapi::getActiveDocumentContext()$path, silent = TRUE)
    if (!inherits(p_try, "try-error") && is.character(p_try) && nzchar(p_try)) p <- p_try
  }

  # Fall back to --file argument when sourced via Rscript
  if (is.null(p) || !nzchar(p)) {
    args <- commandArgs(trailingOnly = FALSE)
    fileArg <- grep("^--file=", args, value = TRUE)
    if (length(fileArg) > 0) p <- sub("^--file=", "", fileArg[1])
  }

  # Final fallback: current working directory
  if (is.null(p) || !nzchar(p)) {
    p <- normalizePath(".", winslash = "/", mustWork = FALSE)
  }

  normalizePath(p, winslash = "/", mustWork = FALSE)
}

.script_path  <- get_script_path()
.script_dir   <- if (dir.exists(.script_path)) .script_path else dirname(.script_path)
.project_root <- normalizePath(file.path(.script_dir, ".."), winslash = "/", mustWork = FALSE)

cat(".script_dir   : ", .script_dir, "\n")
cat(".project_root : ", .project_root, "\n")

candidate_paths <- c(
  file.path(.script_dir,   "DATA", "processed", "modechoice_long.csv"),
  file.path(.project_root, "DATA", "processed", "modechoice_long.csv"),
  file.path(.script_dir,   "data", "processed", "modechoice_long.csv"),
  file.path(.project_root, "data", "processed", "modechoice_long.csv"),
  file.path(.script_dir,   "modechoice_long.csv"),
  file.path(.project_root, "modechoice_long.csv")
)
hit <- which(file.exists(candidate_paths))
if (length(hit) == 0) {
  cat("Checked these paths:\n"); print(candidate_paths)
  stop("modechoice_long.csv not found. Make sure DATA/processed is in the right place.")
}

data_path <- candidate_paths[hit[1]]
cat("Using data file:\n", data_path, "\n")

database <- readr::read_csv(data_path, show_col_types = FALSE)

cat("\nColumns available in dataset:\n"); print(names(database))
cat("\nMissing values summary:\n"); print(colSums(is.na(database)))


### ==========================================================
### STEP 2 – Build choice variables and check consistency
### ==========================================================

alt_map <- c(car = 1L, bus = 2L, air = 3L, rail = 4L)

database <- database %>%
  mutate(alt = tolower(trimws(alt)))

stopifnot("choice" %in% names(database))

choice_map <- database %>%
  group_by(ID, choice_id) %>%
  summarise(
    .groups    = "drop",
    idx        = which.max(choice),
    alt_sel    = tolower(trimws(alt[idx])),
    choice_num = unname(alt_map[alt_sel])
  ) %>%
  select(ID, choice_id, choice_num)

cat("Head of choice_map:\n"); print(head(choice_map))

# Harmonise key types before joining

database <- database %>%
  mutate(
    ID        = as.integer(ID),
    choice_id = as.integer(choice_id)
  )

choice_map <- choice_map %>%
  mutate(
    ID        = as.integer(ID),
    choice_id = as.integer(choice_id)
  )

database <- database %>%
  select(-dplyr::any_of("choice_num")) %>%
  left_join(choice_map, by = c("ID", "choice_id"), suffix = c("", ".cm"))

if (("choice_num.cm" %in% names(database)) && !("choice_num" %in% names(database))) {
  database <- database %>%
    dplyr::rename(choice_num = `choice_num.cm`)
}

if (!("choice_num" %in% names(database))) {
  cat("\n[DIAGNOSTIC] Missing choice_num after join.\n")
  missing_keys <- dplyr::anti_join(
    database %>% dplyr::distinct(ID, choice_id),
    choice_map %>% dplyr::distinct(ID, choice_id),
    by = c("ID", "choice_id")
  )
  print(head(missing_keys, 10))
  stop("choice_num missing after join. Check that each (ID, choice_id) has a chosen alternative.")
}

if (any(is.na(database$choice_num))) {
  bad <- database %>%
    filter(is.na(choice_num)) %>%
    distinct(ID, choice_id) %>%
    head(10)
  print(bad)
  stop("Found NA in choice_num. Please verify the original data.")
}

database$choice_num <- as.integer(database$choice_num)
stopifnot(!any(is.na(database$choice_num)))


### ==========================================================
### STEP 3 – Define starting values and fixed parameters
### ==========================================================

apollo_beta <- c(
  asc_bus    = -1.2,
  asc_air    = -0.2,
  asc_rail   =  0.3,
  mu_time    = -0.05,
  mu_cost    = -0.02,
  sigma_time =  0.03,
  sigma_cost =  0.02,
  b_access   =  0.01,
  b_service  =  0.20
)

apollo_fixed <- c()


### ==========================================================
### STEP 4 – Specify simulation draws for random coefficients
### ==========================================================

apollo_draws <- list(
  interDrawsType = "sobol",
  interNDraws    = 200,
  interNormDraws = c("draws_time", "draws_cost")
)


### ==========================================================
### STEP 5 – Define random coefficients and utility functions
### ==========================================================

apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  with(as.list(c(apollo_beta, apollo_inputs$draws)), {
    b_time <- mu_time + sigma_time * draws_time
    b_cost <- mu_cost + sigma_cost * draws_cost
    return(list(
      b_time = b_time,
      b_cost = b_cost
    ))
  })
}

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs), add = TRUE)

  V <- list()
  V[["car"]]  = 0 + b_time  * (time   * (alt == "car")) + b_cost  * (cost   * (alt == "car"))
  V[["bus"]]  = asc_bus + b_time  * (time   * (alt == "bus")) + b_cost  * (cost   * (alt == "bus")) + b_access* (access * (alt == "bus"))
  V[["air"]]  = asc_air + b_time  * (time   * (alt == "air")) + b_cost  * (cost   * (alt == "air")) + b_access* (access * (alt == "air")) + b_service*(service* (alt == "air"))
  V[["rail"]] = asc_rail + b_time  * (time   * (alt == "rail")) + b_cost  * (cost   * (alt == "rail")) + b_access* (access * (alt == "rail")) + b_service*(service* (alt == "rail"))

  avail_list <- list(car = 1, bus = 1, air = 1, rail = 1)

  mnl_settings <- list(
    alternatives = c(car = 1, bus = 2, air = 3, rail = 4),
    avail        = avail_list,
    choiceVar    = choice_num,
    V            = V
  )

  P <- list()
  P[["model"]] <- apollo_mnl(mnl_settings, functionality)
  # Panel product over repeated observations for each individual
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  # Mixed logit (average over draws)
  P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
  # Final preparation of probabilities
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  
  return(P)
  
}


### ==========================================================
### STEP 6 – Validate Apollo inputs
### ==========================================================

apollo_inputs <- apollo_validateInputs()

invisible(apollo_probabilities(apollo_beta, apollo_inputs, functionality = "validate"))
cat("Validation passed for MMNL model.\n")


### ==========================================================
### STEP 7 – Estimate MMNL model and save output
### ==========================================================

mmnl_model <- apollo_estimate(apollo_beta, apollo_fixed,
                              apollo_probabilities, apollo_inputs)

apollo_modelOutput(mmnl_model)
apollo_saveOutput(mmnl_model)

cat("\nMMNL estimation completed successfully.\n")
##############################################################

##############################################################
# R02_MMNL_Independent.R
# ------------------------------------------------------------
# Project   : Mode Choice Modelling with Apollo
# Researcher: Karina
# Date      : 25 November 2025
#
# Description:
# Mixed Multinomial Logit (MMNL) with *independent* random
# coefficients for travel time and cost, using the original
# WIDE Apollo example data:
#   apollo_modeChoiceData.csv
#
# Alternatives: car (1), bus (2), air (3), rail (4)
# Data columns used (wide format):
#   ID, av_car, av_bus, av_air, av_rail,
#   time_car, time_bus, time_air, time_rail,
#   cost_car, cost_bus, cost_air, cost_rail,
#   choice
#
# Random coefficients (independent normals):
#   b_time ~ Normal(mu_time, sigma_time^2)
#   b_cost ~ Normal(mu_cost, sigma_cost^2)
#
# Structure follows the Apollo v0.3.6 manual (Section 6)
##############################################################


### ==========================================================
### STEP 0 – Load packages and initialise Apollo
### ==========================================================

# Important:
# - Do NOT set options(error = ...) in this script.
# - Restart the R session after experimenting with the error handler.

library(apollo)
library(readr)
library(dplyr)

apollo_initialise()

apollo_control <- list(
  modelName  = "MMNL_independent",
  modelDescr = "Mixed Logit (time & cost, independent normals, WIDE data)",
  indivID    = "ID",
  panelData  = TRUE,
  mixing     = TRUE,
  nCores     = 4    
)


### ==========================================================
### STEP 1 – Locate data file and read WIDE dataset
### ==========================================================

get_script_path <- function(){
  p <- NULL
  
  # RStudio: use the active document path if available.
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    if (tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE)) {
      p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                    error = function(e) NULL)
    }
  }
  
  # If still NULL, try the --file argument (Rscript).
  if (is.null(p) || !nzchar(p)) {
    args    <- commandArgs(trailingOnly = FALSE)
    fileArg <- grep("^--file=", args, value = TRUE)
    if (length(fileArg) > 0){
      p <- sub("^--file=", "", fileArg[1])
    }
  }
  
  # Final fallback: the current working directory.
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
  file.path(.script_dir,   "apollo_modeChoiceData.csv"),
  file.path(.project_root, "apollo_modeChoiceData.csv"),
  file.path(.script_dir,   "DATA", "apollo_modeChoiceData.csv"),
  file.path(.project_root, "DATA", "apollo_modeChoiceData.csv")
)

hit <- which(file.exists(candidate_paths))
if (length(hit) == 0) {
  cat("Checked these paths for apollo_modeChoiceData.csv:\n")
  print(candidate_paths)
  stop("apollo_modeChoiceData.csv not found. Please place it in the repo root or DATA/ folder.")
}

data_path <- candidate_paths[hit[1]]
cat("Using data file:\n", data_path, "\n")

database <- readr::read_csv(data_path, show_col_types = FALSE)

cat("\nColumns in WIDE dataset:\n"); print(names(database))
cat("\nMissing values summary:\n"); print(colSums(is.na(database)))

stopifnot("choice" %in% names(database))


### ==========================================================
### STEP 2 – Define parameters & fixed params
### ==========================================================

# ASC_car = 0 (base)
apollo_beta <- c(
  asc_bus    = -1.0,
  asc_air    = -0.3,
  asc_rail   =  0.2,
  
  mu_time    = -0.05,   # mean time coefficient
  mu_cost    = -0.02,   # mean cost coefficient
  
  sigma_time =  0.02,   # std dev time
  sigma_cost =  0.01    # std dev cost
  # There is no rho in the independent model
)

apollo_fixed <- c()   # all wil be estimate


### ==========================================================
### STEP 3 – Simulation draws (independent normals)
### ==========================================================

apollo_draws <- list(
  interDrawsType = "sobol",
  interNDraws    = 200,
  interNormDraws = c("draws_time", "draws_cost")   # two N(0,1) independent
)


### ==========================================================
### STEP 4 – Random coefficients (independent)
### ==========================================================

# Standard random coefficient:
#   b_time = mu_time + sigma_time * z1
#   b_cost = mu_cost + sigma_cost * z2
# with z1,z2 ~ N(0,1) independent

apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  
  par   <- as.list(apollo_beta)
  draws <- apollo_inputs$draws
  
  mu_time    <- par$mu_time
  mu_cost    <- par$mu_cost
  sigma_time <- par$sigma_time
  sigma_cost <- par$sigma_cost
  
  z1 <- draws$draws_time   # N(0,1)
  z2 <- draws$draws_cost   # N(0,1)
  
  b_time <- mu_time + sigma_time * z1
  b_cost <- mu_cost + sigma_cost * z2
  
  return(list(
    b_time = b_time,
    b_cost = b_cost
  ))
}


### ==========================================================
### STEP 5 – Utility functions & probabilities
### ==========================================================

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){
  
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs), add = TRUE)
  
  # 5.1 Utility functions (wide)
  V <- list()
  
  V[["car"]]  <- 0 +
    b_time * time_car +
    b_cost * cost_car
  
  V[["bus"]]  <- asc_bus +
    b_time * time_bus +
    b_cost * cost_bus
  
  V[["air"]]  <- asc_air +
    b_time * time_air +
    b_cost * cost_air
  
  V[["rail"]] <- asc_rail +
    b_time * time_rail +
    b_cost * cost_rail
  
  # 5.2 Availability
  avail <- list(
    car  = av_car,
    bus  = av_bus,
    air  = av_air,
    rail = av_rail
  )
  
  # 5.3 MNL kernel
  mnl_settings <- list(
    alternatives = c(car = 1, bus = 2, air = 3, rail = 4),
    avail        = avail,
    choiceVar    = choice,
    V            = V,
    componentName = "MNL"
  )
  
  P <- list()
  P[["model"]] <- apollo_mnl(mnl_settings, functionality)
  
  # Panel product over repeated choice situations per ID
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  
  # Average over draws (Mixed Logit)
  P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
  
  # Final prep
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  
  return(P)
}


### ==========================================================
### STEP 6 – Validate inputs
### ==========================================================

apollo_inputs <- apollo_validateInputs()

invisible(apollo_probabilities(apollo_beta, apollo_inputs, functionality = "validate"))
cat("Validation passed for MMNL_independent.\n")


### ==========================================================
### STEP 7 – Estimate & save output
### ==========================================================

mmnl_indep_model <- apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

apollo_modelOutput(mmnl_indep_model)
apollo_saveOutput(mmnl_indep_model)

cat("\nMMNL_independent estimation completed.\n")
##############################################################

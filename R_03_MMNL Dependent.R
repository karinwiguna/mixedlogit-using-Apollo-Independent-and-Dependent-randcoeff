##############################################################
# R03_MMNL_Dependent.R
# ------------------------------------------------------------
# Project   : Mode Choice Modelling with Apollo
# Researcher: Karina
# Date      : 25 November 2025
#
# Description:
# Mixed Multinomial Logit (MMNL) with *correlated* random
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
# Random coefficients:
#   b_time ~ Normal(mu_time, sigma_time^2)
#   b_cost ~ Normal(mu_cost, sigma_cost^2)
#   corr(b_time, b_cost) = rho   (estimated)
#
# Implementation follows the logic of Section 6 (random
# heterogeneity) of the Apollo manual v0.3.6.
##############################################################


### ==========================================================
### STEP 0 – Load packages and initialise Apollo
### ==========================================================

# IMPORTANT:
# - Do NOT set options(error = ...) in this script.
# - Ensure no custom error handler is active.

library(apollo)
library(readr)
library(dplyr)

apollo_initialise()

apollo_control <- list(
  modelName  = "MMNL_dependent",
  modelDescr = "Mixed Logit (time & cost, correlated normals, WIDE data)",
  indivID    = "ID",      # panel individual ID
  panelData  = TRUE,      # repeated choices per ID
  mixing     = TRUE,      # we use random coefficients
  nCores     = 4          
)


### ==========================================================
### STEP 1 – Locate data file and read WIDE dataset
### ==========================================================

get_script_path <- function(){
  p <- NULL
  
  # RStudio: use the active document path if available
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    if (tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE)) {
      p <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                    error = function(e) NULL)
    }
  }
  
  # If still NULL, try the --file argument (Rscript)
  if (is.null(p) || !nzchar(p)) {
    args    <- commandArgs(trailingOnly = FALSE)
    fileArg <- grep("^--file=", args, value = TRUE)
    if (length(fileArg) > 0){
      p <- sub("^--file=", "", fileArg[1])
    }
  }
  
  # Final fallback: the current working directory
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

# Search for the wide-format file in several possible locations:
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

# === Load WIDE data (exactly as in Apollo example) ===
database <- readr::read_csv(data_path, show_col_types = FALSE)

cat("\nColumns in WIDE dataset:\n"); print(names(database))
cat("\nMissing values summary:\n"); print(colSums(is.na(database)))

# We assume choice is coded 1=car, 2=bus, 3=air, 4=rail, as in Apollo example
stopifnot("choice" %in% names(database))


### ==========================================================
### STEP 2 – Define parameters (starting values) & fixed params
### ==========================================================

# ASC_car is used as the base (0), so it is not included.
# Random coefficients: b_time, b_cost (correlated normals)
apollo_beta <- c(
  asc_bus    = -1.0,
  asc_air    = -0.3,
  asc_rail   =  0.2,
  
  mu_time    = -0.05,   # mean of time coefficient
  mu_cost    = -0.02,   # mean of cost coefficient
  
  sigma_time =  0.02,   # std dev of time coefficient
  sigma_cost =  0.01,   # std dev of cost coefficient
  
  rho_raw    =  0       # transformed to rho in (-1,1)
)

# All parameters are estimated (none are fixed).
apollo_fixed <- c()


### ==========================================================
### STEP 3 – Simulation draws for random coefficients
### ==========================================================

# Using inter-individual (panel) draws.
apollo_draws <- list(
  interDrawsType = "sobol",                        # Sobol draws (manual section 6)
  interNDraws    = 200,                            # This can be increased later (e.g., 500 or 1000)
  interNormDraws = c("draws_time", "draws_cost")   # two standard normal draws
)


### ==========================================================
### STEP 4 – Random coefficients with correlation
### ==========================================================

# b_time and b_cost ~ bivariate normal with correlation rho.
# Standard implementation (following Train and the Apollo manual):
#   z1 = draws_time ~ N(0,1)
#   z2 = draws_cost ~ N(0,1) independent
#   b_time = mu_time + sigma_time * z1
#   b_cost = mu_cost + sigma_cost * ( rho * z1 + sqrt(1 - rho^2) * z2 )

apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  
  # access parameters (mu, sigma, rho_raw) and the draws (draws_time, draws_cost)
  par   <- as.list(apollo_beta)
  draws <- apollo_inputs$draws
  
  mu_time    <- par$mu_time
  mu_cost    <- par$mu_cost
  sigma_time <- par$sigma_time
  sigma_cost <- par$sigma_cost
  
  # transform rho_raw --> rho in (-1,1) using the logistic function
  rho_raw <- par$rho_raw
  rho     <- 2 * plogis(rho_raw) - 1   # rho = 2*logit^-1(rho_raw) - 1
  
  z1 <- draws$draws_time   # N(0,1)
  z2 <- draws$draws_cost   # N(0,1), independent
  
  b_time <- mu_time + sigma_time * z1
  b_cost <- mu_cost + sigma_cost * ( rho * z1 + sqrt(1 - rho^2) * z2 )
  
  return(list(
    b_time = b_time,
    b_cost = b_cost
  ))
}


### ==========================================================
### STEP 5 – Define utility functions & likelihood
### ==========================================================

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){
  
  # Attach:
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs), add = TRUE)
  
  # --------------------------------------------------------
  # 5.1 Utility functions in WIDE format
  # --------------------------------------------------------
  V <- list()
  
  # car (base ASC = 0)
  V[["car"]]  <- 0 +
    b_time * time_car +
    b_cost * cost_car
  
  # bus
  V[["bus"]]  <- asc_bus +
    b_time * time_bus +
    b_cost * cost_bus
  
  # air
  V[["air"]]  <- asc_air +
    b_time * time_air +
    b_cost * cost_air
  
  # rail
  V[["rail"]] <- asc_rail +
    b_time * time_rail +
    b_cost * cost_rail
  
  # --------------------------------------------------------
  # 5.2 Availability
  # --------------------------------------------------------
  avail <- list(
    car  = av_car,
    bus  = av_bus,
    air  = av_air,
    rail = av_rail
  )
  
  # --------------------------------------------------------
  # 5.3 MNL kernel inside Mixed Logit
  # --------------------------------------------------------
  mnl_settings <- list(
    alternatives = c(car = 1, bus = 2, air = 3, rail = 4),
    avail        = avail,
    choiceVar    = choice,
    V            = V,
    componentName = "MNL"
  )
  
  P <- list()
  P[["model"]] <- apollo_mnl(mnl_settings, functionality)
  
  # Panel product over repeated choice situations for each ID
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  
  # Average over inter-individual draws (Mixed Logit)
  P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
  
  # Prepare final probabilities
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  
  return(P)
}

### ==========================================================
### STEP 6 – Validate inputs
### ==========================================================

apollo_inputs <- apollo_validateInputs()

# Optional quick validation run
invisible(apollo_probabilities(apollo_beta, apollo_inputs, functionality = "validate"))
cat("Validation passed for MMNL_dependent.\n")


### ==========================================================
### STEP 7 – Estimate & save output
### ==========================================================

mmnl_dep_model <- apollo_estimate(
  apollo_beta,
  apollo_fixed,
  apollo_probabilities,
  apollo_inputs
)

apollo_modelOutput(mmnl_dep_model)
apollo_saveOutput(mmnl_dep_model)

cat("\nMMNL_dependent estimation completed.\n")
##############################################################

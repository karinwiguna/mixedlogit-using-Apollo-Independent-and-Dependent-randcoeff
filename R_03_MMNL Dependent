##############################################################
# R03_MMNL_Dependent.R
# ------------------------------------------------------------
# Project   : Mode Choice Modelling with Apollo
# Researcher: Karina
# Date      : 30 October 2025
#
# Description:
# Mixed Multinomial Logit (MMNL) with *correlated* (dependent)
# random coefficients for travel time and cost, using the
# original wide-format Apollo example data:
#   DATA/apollo_modeChoiceData.csv
#
# Alternatives : car, bus, air, rail
# Choice var   : choice (1=car, 2=bus, 3=air, 4=rail)
##############################################################

### ==========================================================
### STEP 0 – Load packages and initialise Apollo
### ==========================================================

library(apollo)

apollo_initialise()

apollo_control <- list(
  modelName  = "MMNL_dependent",
  modelDescr = "Mixed Logit with correlated time & cost coefficients (wide data)",
  indivID    = "ID",       # panel over repeated RP/SP tasks
  panelData  = TRUE,
  mixing     = TRUE,
  nCores     = 1           
)

### ==========================================================
### STEP 1 – Locate project paths and read WIDE data
### ==========================================================

get_script_path <- function(){
  p <- NULL

  # 1) Kalau jalan di RStudio, pakai path dokumen aktif
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    if (tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE)) {
      p_try <- try(rstudioapi::getActiveDocumentContext()$path, silent = TRUE)
      if (!inherits(p_try, "try-error") && is.character(p_try) && nzchar(p_try)) {
        p <- p_try
      }
    }
  }

  # 2) Kalau via Rscript, pakai argumen --file
  if (is.null(p) || !nzchar(p)) {
    args    <- commandArgs(trailingOnly = FALSE)
    fileArg <- grep("^--file=", args, value = TRUE)
    if (length(fileArg) > 0) p <- sub("^--file=", "", fileArg[1])
  }

  # 3) Fallback: current working directory
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

# Cari file data dalam DATA/ atau root (wide format)
candidate_paths <- c(
  file.path(.script_dir,   "DATA", "apollo_modeChoiceData.csv"),
  file.path(.project_root, "DATA", "apollo_modeChoiceData.csv"),
  file.path(.script_dir,   "apollo_modeChoiceData.csv"),
  file.path(.project_root, "apollo_modeChoiceData.csv")
)

hit <- which(file.exists(candidate_paths))
if (length(hit) == 0) {
  cat("Checked these paths:\n"); print(candidate_paths)
  stop("apollo_modeChoiceData.csv not found in expected locations (DATA/ or script folder).")
}

data_path <- candidate_paths[hit[1]]
cat("Using data file:\n", data_path, "\n")

database <- read.csv(data_path, header = TRUE)

cat("\nColumns available in dataset:\n"); print(names(database))
cat("\nMissing values summary:\n"); print(colSums(is.na(database)))

# Pastikan choice ada dan bernilai 1..4
stopifnot("choice" %in% names(database))
database$choice <- as.integer(database$choice)
stopifnot(all(database$choice %in% 1:4))

### ==========================================================
### STEP 2 – Define starting values & fixed parameters
### ==========================================================

# Parameter tetap + mean & st.dev random coeff + korelasi
apollo_beta <- c(
  # Asc (relative to car)
  asc_bus      = -1.2,
  asc_air      = -0.2,
  asc_rail     =  0.3,

  # Mean random coefficients (time & cost)
  mu_time      = -0.05,
  mu_cost      = -0.02,

  # Std dev (heterogeneity)
  sigma_time   =  0.03,
  sigma_cost   =  0.02,

  # Korelasi "raw" antara time & cost (nanti di-transform ke [-1,1])
  rho_raw_tc   =  0.3,

  # Deterministic coefficients
  b_access     =  0.01,
  b_service    =  0.20
)

# Semua parameter diestimasi (tidak ada yang difix)
apollo_fixed <- c()

### ==========================================================
### STEP 3 – Simulation draws (for random coefficients)
### ==========================================================

apollo_draws <- list(
  interDrawsType = "sobol",                  # Sobol sequence
  interNDraws    = 200,                      # bisa dinaikkan 500/1000 nanti
  interNormDraws = c("draws_t", "draws_c")   # dua N(0,1) independen
)

### ==========================================================
### STEP 4 – Random coefficients with correlation
### ==========================================================

apollo_randCoeff <- function(apollo_beta, apollo_inputs){
  # apollo_inputs$draws berisi draws_t & draws_c (N(0,1) Sobol)
  with(as.list(c(apollo_beta, apollo_inputs$draws)), {

    # Transformasi rho_raw_tc -> rho_tc di [-1,1] pakai tanh
    rho_tc <- tanh(rho_raw_tc)

    # b_time ~ N(mu_time, sigma_time^2)
    b_time <- mu_time + sigma_time * draws_t

    # b_cost = mu_cost + sigma_cost * (rho * draws_t + sqrt(1-rho^2) * draws_c)
    # Ini implementasi korelasi standar untuk dua Normal
    b_cost <- mu_cost + sigma_cost * (rho_tc * draws_t + sqrt(1 - rho_tc^2) * draws_c)

    return(list(
      b_time = b_time,
      b_cost = b_cost
    ))
  })
}

### ==========================================================
### STEP 5 – Utility functions & probability structure
### ==========================================================

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){

  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs), add = TRUE)

  # Utility per alternatif (pakai variabel WIDE: time_car, cost_car, dst.)
  V <- list()

  V[["car"]]  <- 0 +
    b_time * time_car  +
    b_cost * cost_car

  V[["bus"]]  <- asc_bus +
    b_time * time_bus  +
    b_cost * cost_bus  +
    b_access  * access_bus

  V[["air"]]  <- asc_air +
    b_time * time_air  +
    b_cost * cost_air  +
    b_access  * access_air +
    b_service * service_air

  V[["rail"]] <- asc_rail +
    b_time * time_rail +
    b_cost * cost_rail +
    b_access  * access_rail +
    b_service * service_rail

  # Availability (pakai av_car, av_bus, ...)
  avail <- list(
    car  = av_car,
    bus  = av_bus,
    air  = av_air,
    rail = av_rail
  )

  mnl_settings <- list(
    alternatives = c(car = 1, bus = 2, air = 3, rail = 4),
    avail        = avail,
    choiceVar    = choice,   # sudah numeric 1..4
    V            = V
  )

  P <- list()
  P[["model"]] <- apollo_mnl(mnl_settings, functionality)

  # Panel product (karena tiap ID punya beberapa RP/SP tasks)
  P <- apollo_panelProd(P, apollo_inputs, functionality)

  # Rata-rata atas draws → mixed logit
  P <- apollo_avgInterDraws(P, apollo_inputs, functionality)

  # Final preparation
  P <- apollo_prepareProb(P, apollo_inputs, functionality)

  return(P)
}

### ==========================================================
### STEP 6 – Validate inputs
### ==========================================================

apollo_inputs <- apollo_validateInputs()
invisible(apollo_probabilities(apollo_beta, apollo_inputs, functionality = "validate"))
cat("Validation passed for MMNL dependent model.\n")

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

cat("\nMMNL dependent estimation completed.\n")
##############################################################

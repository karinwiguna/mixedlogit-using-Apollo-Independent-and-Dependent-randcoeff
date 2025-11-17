##############################################################
# R01_MNL_Baseline.R
# ------------------------------------------------------------
# Project   : Mode Choice Modelling with Apollo
# Researcher: Karina
# Date      : 30 October 2025
#
# Description:
# This script estimates a baseline Multinomial Logit (MNL) model
# using the Apollo package in R. The model uses open-source 
# mode choice data with 4 alternatives: car, bus, air, and rail.
# ------------------------------------------------------------
# Dependencies: apollo, readr, dplyr
##############################################################


### ==========================================================
### 0. INITIAL SETUP
### ==========================================================

# Load libraries
library(apollo)
library(readr)
library(dplyr)
library(tidyr)   # optional — only needed if reshape or clean data again

# Initialize Apollo
apollo_initialise()

# Define control settings
apollo_control = list(
  modelName  = "MNL_baseline",
  modelDescr = "Multinomial Logit baseline (Apollo processed long data)",
  indivID    = "ID",
  nCores     = 1
)


### ==========================================================
### 1. LOAD & CHECK DATA
### ==========================================================

get_script_path <- function(){
  p <- NULL
  if (requireNamespace("rstudioapi", quietly=TRUE)) {
    if (tryCatch(rstudioapi::isAvailable(), error=function(e) FALSE)) {
      p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error=function(e) NULL)
    }
  }
  if (is.null(p) || !nzchar(p)) {
    # fallback jika run via terminal/source
    args <- commandArgs(trailingOnly = FALSE)
    fileArg <- grep("^--file=", args, value = TRUE)
    if (length(fileArg) > 0) p <- sub("^--file=", "", fileArg[1])
  }
  if (is.null(p) || !nzchar(p)) p <- normalizePath(".", winslash = "/", mustWork = FALSE)
  return(normalizePath(p, winslash = "/", mustWork = FALSE))
}

.script_path  <- get_script_path()
.script_dir   <- if (dir.exists(.script_path)) .script_path else dirname(.script_path)
.project_root <- normalizePath(file.path(.script_dir, ".."), winslash = "/", mustWork = FALSE)

# Buat path data tanpa setwd:
cat(".script_dir   : ", .script_dir, "\n")
cat(".project_root : ", .project_root, "\n")

# Coba beberapa kemungkinan lokasi DATA:
candidates <- c(
  file.path(.script_dir,   "DATA", "processed", "modechoice_long.csv"), # DATA di dalam folder skrip (kasusmu)
  file.path(.project_root, "DATA", "processed", "modechoice_long.csv"), # DATA di parent
  file.path(.script_dir,   "data", "processed", "modechoice_long.csv"), # fallback huruf kecil
  file.path(.project_root, "data", "processed", "modechoice_long.csv")
)
hit <- which(file.exists(candidates))
if (length(hit) == 0) {
  cat("Checked these paths:\n"); print(candidates)
  stop("modechoice_long.csv not found. Make sure DATA/processed is next to the R01 script *or* one level up.")
}
data_path <- candidates[hit[1]]
cat("Using data file:\n", data_path, "\n")


# 1) Load data via portable path
database <- readr::read_csv(data_path, show_col_types = FALSE)

# 2) Cek cepat
cat("\nColumns available in dataset:\n"); print(names(database))
cat("\nMissing values summary:\n"); print(colSums(is.na(database)))

### ==========================================================
### 2. DEFINE MODEL PARAMETERS
### ==========================================================

## Alt mapping --> Code 1,2,3,4
alt_map <- c(car=1L, bus=2L, air=3L, rail=4L)

# Normalise alt to avoid 'car ' vs 'Car' issues
database <- database %>%
  mutate(alt = tolower(trimws(alt)))
cat("Names in database:\n"); print(names(database))
cat("Class of keys in database:\n"); print(sapply(database[,c("ID","choice_id")], class))

# Pastikan 'choice' benar2 ada dan numeric/binary
stopifnot("choice" %in% names(database))

# Take choices per (ID, choice_id) from baris yang choice==1, lalu lekatkan ke semua baris choice set tsb
choice_map <- database %>%
  group_by(ID, choice_id) %>%
  summarise(
    .groups = "drop",
    # baris pilihan (kalau tidak ada 1, tetap ambil idx max -> asumsi ada satu alt terpilih)
    idx = which.max(choice),
    alt_sel = tolower(trimws(alt[idx])),
    choice_num = unname(alt_map[alt_sel])
  ) %>%
  select(ID, choice_id, choice_num)

cat("Head of choice_map:\n"); print(head(choice_map))
cat("Class of keys in choice_map:\n"); print(sapply(choice_map[,c("ID","choice_id")], class))

# Samakan tipe key dulu
database <- database %>%
  mutate(
    ID = as.integer(ID),
    choice_id = as.integer(choice_id)
  )

choice_map <- choice_map %>%
  mutate(
    ID = as.integer(ID),
    choice_id = as.integer(choice_id)
  )

# Bersihkan sisa kolom choice_num lama (kalau ada), lalu JOIN dengan suffix jelas
database <- database %>%
  select(-dplyr::any_of("choice_num")) %>%
  left_join(choice_map, by = c("ID","choice_id"), suffix = c("", ".cm"))

# Jika kolom dari choice_map bernama choice_num.cm, normalkan kembali jadi 'choice_num'
if (("choice_num.cm" %in% names(database)) && !("choice_num" %in% names(database))) {
  database <- database %>%
    dplyr::rename(choice_num = `choice_num.cm`)
}

# Diagnostik jika tetap tidak ada choice_num
if (!("choice_num" %in% names(database))) {
  cat("\n[DIAGNOSTIC] Keys present in database but missing in choice_map (showing 10):\n")
  miss_keys <- dplyr::anti_join(
    database %>% dplyr::distinct(ID, choice_id),
    choice_map %>% dplyr::distinct(ID, choice_id),
    by = c("ID","choice_id")
  )
  print(head(miss_keys, 10))
  stop("choice_num missing after join. Kemungkinan: (1) 'choice_id' tidak ada/berbeda; (2) tipe data key beda; (3) nilai key tidak match.")
}

# Cek NA pada choice_num
if (any(is.na(database$choice_num))) {
  bad <- database %>% dplyr::filter(is.na(choice_num)) %>% dplyr::distinct(ID, choice_id) %>% head(10)
  print(bad)
  stop("Found NA in choice_num. Cek apakah 'alt' hanya {car,bus,air,rail} dan tiap (ID,choice_id) punya tepat satu pilihan (choice==1).")
}

# Pastikan integer
database$choice_num <- as.integer(database$choice_num)

# --- make sure choice_num exists & is clean
stopifnot("choice_num" %in% names(database))
database$choice_num <- as.integer(database$choice_num)
stopifnot(!any(is.na(database$choice_num)))

# Coefficients (starting values)
#Starting values (car jadi base → ASC_car = 0)
apollo_beta <- c(
  asc_bus   = -1.5,
  asc_air   = -0.3,
  asc_rail  =  0.3,
  l_time    = log(0.005),  # ~ -5.298; b_time = -exp(l_time) < 0
  l_cost    = log(0.010),  # ~ -4.605; b_cost = -exp(l_cost) < 0
  b_access  =  0.01,    # ada kolom access di data
  b_service =  0.20     # ada kolom service di data
)

# No fixed parameters (all will be estimated)
apollo_fixed <- c()


### ==========================================================
### 3. VALIDATE INPUTS
### ==========================================================

# Validate all Apollo inputs before estimation
apollo_inputs <- apollo_validateInputs()
# Expect some warnings about NA or panel data – that’s fine.


### ==========================================================
### 4. DEFINE MODEL STRUCTURE (MNL)
### ==========================================================

apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality="estimate"){
  
  # Attach parameters and database for direct use
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  # 4.1 Define utility functions for each alternative
  b_time = -exp(l_time)
  b_cost = -exp(l_cost)
  
  V <- list()
  V[["car"]]  = 0 +
    b_time  * (time   * (alt=="car")) +
    b_cost  * (cost   * (alt=="car"))
  V[["bus"]]  =
    asc_bus +
    b_time  * (time   * (alt=="bus")) +
    b_cost  * (cost   * (alt=="bus")) +
    b_access* (access * (alt=="bus"))
  V[["air"]]  =
    asc_air +
    b_time  * (time   * (alt=="air")) +
    b_cost  * (cost   * (alt=="air")) +
    b_access* (access * (alt=="air")) +
    b_service*(service* (alt=="air"))
  V[["rail"]] =
    asc_rail +
    b_time  * (time   * (alt=="rail")) +
    b_cost  * (cost   * (alt=="rail")) +
    b_access* (access * (alt=="rail")) +
    b_service*(service* (alt=="rail"))
  
  # 4.2 Model settings
  avail_list <- list(car=1, bus=1, air=1, rail=1)
  
  mnl_settings <- list(
    alternatives = c(car=1, bus=2, air=3, rail=4),
    avail        = avail_list,
    choiceVar    = choice_num,
    V            = V
  )
  
  # 4.3 Compute probabilities
  P <- list()
  P[["model"]] <- apollo_mnl(mnl_settings, functionality)
  P <- apollo_panelProd(P, apollo_inputs, functionality)
  P <- apollo_prepareProb(P, apollo_inputs, functionality)
  
  return(P)
}
# --- validate only after database is final and probabilities are defined
if (exists("apollo_inputs")) rm(apollo_inputs)
apollo_inputs <- apollo_validateInputs()

# optional: quick check
invisible(apollo_probabilities(apollo_beta, apollo_inputs, functionality="validate"))
cat("Validation passed (choice_num visible to Apollo).\n")


### ==========================================================
### 5. ESTIMATE MODEL
### ==========================================================

mnl_model <- apollo_estimate(apollo_beta, apollo_fixed,
                             apollo_probabilities, apollo_inputs)


### ==========================================================
### 6. OUTPUT RESULTS
### ==========================================================

# Display results in console
apollo_modelOutput(mnl_model)

# Save full estimation output to file (in working directory)
apollo_saveOutput(mnl_model)


### ==========================================================
### 7. INTERPRETATION NOTES
### ==========================================================
# Parameterisation:
# - Travel time and cost are estimated in log-space:
#     b_time = -exp(l_time)  < 0   (ensured negative)
#     b_cost = -exp(l_cost)  < 0
# - Interpret substantive effects using b_time and b_cost (transformed values),
#   not the raw l_time/l_cost estimates.

# Signs and meanings:
# - b_time < 0 : longer travel time decreases utility (disutility of time)
# - b_cost < 0 : higher travel cost decreases utility (disutility of cost)
# - ASC_bus / ASC_air / ASC_rail are relative to car (ASC_car = 0)
# - b_access > 0 and b_service > 0 increase utility (improved access/service)

# Policy-relevant metrics (check units: time = minutes or hours; cost = £):
# - VTTS (£/min)  = b_time / b_cost
#   VTTS (£/hour) = 60 * (b_time / b_cost)
# - WTP_access  =  b_access / (-b_cost)   ( £ per 1 unit of access )
# - WTP_service =  b_service / (-b_cost)  ( £ per 1 unit of service )

# Goodness-of-fit:
# - LL(final) should be greater than LL(0)
# - Typical Rho-squared (ρ²) range for MNL: 0.1–0.3
# - |t-ratio| > 1.96 → statistically significant at 5% level

# Practical notes:
# - Coefficients and VTTS are sensitive to variable scaling (minutes vs hours; £ vs tens of £)
# - Extremely negative l_time (→ b_time ≈ 0) may indicate rescaling is needed
# - If signs and magnitudes are reasonable, proceed to R02_MixedLogit for random parameters
##############################################################


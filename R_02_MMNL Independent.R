#!/usr/bin/env Rscript
##############################################################
# R_02_MMNL_Independent.R
# ------------------------------------------------------------
# Mixed logit model with independent random coefficients.
# - Cleans and validates the processed long-format data
# - Scales continuous variables
# - Estimates an MMNL model with Halton/Sobol draws
# - Provides diagnostics (log-likelihood, gradient norm,
#   Hessian definiteness) and saves a concise summary
##############################################################

## STEP 0 – Load required packages and helper functions --------------------

suppressPackageStartupMessages({
  library(apollo)
  library(dplyr)
  library(readr)
  library(tidyr)
})

`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

mmnl_model <- function(processed_path = "DATA/processed/modechoice_long.csv",
                       output_dir = "output",
                       draws = 1000,
                       seed = 1234) {
  ## STEP 1 – Initialise Apollo session and output folders -----------------
  set.seed(seed)
  apollo_initialise()

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  ## STEP 2 – Load processed LONG data and validate structure --------------
  if (!file.exists(processed_path)) {
    stop(sprintf("Processed file not found at '%s'. Please run the data processing script first.",
                 processed_path))
  }

  raw_long <- read_csv(processed_path, show_col_types = FALSE)

  required_cols <- c("ID", "choice_id", "alt", "time", "cost", "avail",
                     "choice", "choice_num")
  missing_cols <- setdiff(required_cols, names(raw_long))
  if (length(missing_cols) > 0) {
    stop(sprintf("Processed data is missing required columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  ## STEP 3 – Clean and filter choice situations ---------------------------
  alt_levels <- c("car", "bus", "air", "rail")
  database <- raw_long %>%
    filter(.data$alt %in% alt_levels) %>%
    mutate(
      time = suppressWarnings(as.numeric(.data$time)),
      cost = suppressWarnings(as.numeric(.data$cost)),
      avail = if_else(is.na(.data$avail), 0L, as.integer(.data$avail)),
      choice_flag = as.integer(.data$choice),
      choice_num = as.integer(.data$choice_num)
    ) %>%
    filter(!is.na(.data$time), !is.na(.data$cost)) %>%
    group_by(.data$ID, .data$choice_id) %>%
    filter(sum(.data$choice_flag, na.rm = TRUE) == 1L) %>%
    ungroup()

  choice_lookup <- database %>%
    filter(.data$choice_flag == 1L) %>%
    transmute(ID = .data$ID, choice_id = .data$choice_id,
              chosen_alt = match(.data$alt, alt_levels))

  database <- database %>%
    left_join(choice_lookup, by = c("ID", "choice_id")) %>%
    mutate(choice = as.integer(.data$chosen_alt)) %>%
    filter(!is.na(.data$choice))

  ## STEP 4 – Build availability and scaling information -------------------
  availability_wide <- database %>%
    select(.data$ID, .data$choice_id, .data$alt, .data$avail) %>%
    distinct() %>%
    pivot_wider(names_from = .data$alt, values_from = .data$avail,
                names_prefix = "av_")

  database <- database %>%
    left_join(availability_wide, by = c("ID", "choice_id")) %>%
    mutate(across(starts_with("av_"), ~ replace_na(as.integer(.), 0L)))

  scale_stats <- database %>%
    filter(.data$avail == 1L) %>%
    summarise(
      time_mean = mean(.data$time),
      time_sd   = sd(.data$time),
      cost_mean = mean(.data$cost),
      cost_sd   = sd(.data$cost)
    )

  time_mean <- scale_stats$time_mean %||% 0
  time_sd   <- scale_stats$time_sd
  if (is.null(time_sd) || is.na(time_sd) || time_sd == 0) time_sd <- 1
  cost_mean <- scale_stats$cost_mean %||% 0
  cost_sd   <- scale_stats$cost_sd
  if (is.null(cost_sd) || is.na(cost_sd) || cost_sd == 0) cost_sd <- 1

  database <- database %>%
    mutate(
      time_sc = (.data$time - time_mean) / time_sd,
      cost_sc = (.data$cost - cost_mean) / cost_sd
    ) %>%
    arrange(.data$ID, .data$choice_id,
            factor(.data$alt, levels = alt_levels))

  scale_info <- list(
    time_mean = time_mean,
    time_sd = time_sd,
    cost_mean = cost_mean,
    cost_sd = cost_sd
  )

  ## STEP 5 – Configure apollo_control and simulation draws ----------------
  apollo_control <- list(
    modelName       = "MMNL_independent_scaled",
    modelDescr      = "MMNL with scaled time/cost and lognormal RC",
    indivID         = "ID",
    panelData       = TRUE,
    mixing          = TRUE,
    nCores          = {
      if (requireNamespace("parallel", quietly = TRUE)) {
        max(1L, min(4L, parallel::detectCores() - 1L))
      } else 1L
    },
    outputDirectory = output_dir
  )

  apollo_draws <- list(
    interDrawsType = "sobol",
    interNDraws    = draws,
    interUnifDraws = c(),
    interNormDraws = c("draws_time", "draws_cost")
  )

  ## STEP 6 – Define parameters, random coefficients, and utilities --------
  apollo_randCoeff <- function(apollo_beta, apollo_inputs) {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))

    sd_time <- exp(log_sigma_time)
    sd_cost <- exp(log_sigma_cost)

    list(
      b_time = -exp(mu_time + sd_time * draws_time),
      b_cost = -exp(mu_cost + sd_cost * draws_cost)
    )
  }

  apollo_beta <- c(
    asc_bus       = -0.3,
    asc_air       = -0.1,
    asc_rail      =  0.2,
    mu_time       = log(1.0),
    log_sigma_time= log(0.4),
    mu_cost       = log(1.2),
    log_sigma_cost= log(0.5)
  )
  apollo_fixed <- c()

  apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate") {
    apollo_attach(apollo_beta, apollo_inputs)
    on.exit(apollo_detach(apollo_beta, apollo_inputs))

    ## STEP 6 (continued) – Construct utilities for each alternative --------
    V <- list()
    V[["car"]]  <- 0 +        b_time * (time_sc * (alt == "car"))  + b_cost * (cost_sc * (alt == "car"))
    V[["bus"]]  <- asc_bus +  b_time * (time_sc * (alt == "bus"))  + b_cost * (cost_sc * (alt == "bus"))
    V[["air"]]  <- asc_air +  b_time * (time_sc * (alt == "air"))  + b_cost * (cost_sc * (alt == "air"))
    V[["rail"]] <- asc_rail + b_time * (time_sc * (alt == "rail")) + b_cost * (cost_sc * (alt == "rail"))

    mnl_settings <- list(
      alternatives = setNames(seq_along(alt_levels), alt_levels),
      avail        = list(
        car  = av_car,
        bus  = av_bus,
        air  = av_air,
        rail = av_rail
      ),
      choiceVar    = choice,
      V            = V
    )

    P <- list()
    P[["model"]] <- apollo_mnl(mnl_settings, functionality)
    P <- apollo_panelProd(P, apollo_inputs, functionality)
    P <- apollo_avgInterDraws(P, apollo_inputs, functionality)
    apollo_prepareProb(P, apollo_inputs, functionality)
  }

  ## STEP 7 – Validate inputs and estimate MMNL model ----------------------
  apollo_inputs <- apollo_validateInputs(
    database        = database,
    apollo_control  = apollo_control,
    apollo_draws    = apollo_draws,
    apollo_randCoeff= apollo_randCoeff
  )
  invisible(apollo_probabilities(apollo_beta, apollo_inputs, functionality = "validate"))

  model <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

  ## STEP 7 (continued) – Save Apollo output and custom diagnostics ---------
  apollo_modelOutput(model)
  apollo_saveOutput(model)

  diag_lines <- list()
  diag_lines[["timestamp"]] <- paste("Timestamp:", format(Sys.time(), tz = "UTC"))
  ll_start <- model$LLStart %||% model$LL0 %||% NA
  diag_lines[["ll_start"]] <- paste("Initial log-likelihood:", signif(ll_start, 6))
  diag_lines[["final_ll"]] <- paste("Final log-likelihood:", signif(model$maximum, 6))
  grad_norm <- sqrt(sum(model$grad^2))
  diag_lines[["grad_norm"]] <- paste("Gradient 2-norm:", signif(grad_norm, 6))

  hessian_pd <- NA
  eig_text <- "N/A"
  if (!is.null(model$hessian)) {
    eig_vals <- tryCatch(
      eigen(-model$hessian, symmetric = TRUE, only.values = TRUE)$values,
      error = function(cond) NA
    )
    if (all(is.finite(eig_vals))) {
      hessian_pd <- all(eig_vals > 0)
      eig_text <- paste(signif(head(eig_vals, 10), 6), collapse = ", ")
    }
  }
  diag_lines[["hessian_def"]] <- paste("-H positive definite:", hessian_pd)
  diag_lines[["eigenvalues"]] <- paste("Leading eigenvalues(-H):", eig_text)

  coef_table <- capture.output(apollo_modelOutput(model))

  scale_lines <- sprintf("%s (original): mean=%.4f sd=%.4f",
                         c("time", "cost"),
                         c(scale_info$time_mean, scale_info$cost_mean),
                         c(scale_info$time_sd, scale_info$cost_sd))

  summary_path <- file.path(output_dir, "MMNL_independent_summary.txt")
  writeLines(c(
    "MMNL Independent Model Summary",
    "================================",
    unlist(diag_lines, use.names = FALSE),
    "",
    "Scaling information (original units):",
    scale_lines,
    "",
    "Model output:",
    coef_table
  ), con = summary_path)

  invisible(list(model = model, scale_info = scale_info, summary_path = summary_path))
}

if (sys.nframe() == 0) {
  mmnl_model()
}

# mixedlogit-debug

Mixed logit (MMNL) debugging playground built around the Apollo R package. The
project contains data-preparation utilities, a baseline multinomial logit model,
and an MMNL specification with independent random coefficients. A convenience
runner script automates preprocessing and estimation from a single command.

The workflow is intentionally minimal: build the processed long-format data,
run the baseline MNL for a quick check, and estimate the mixed logit with
independent random coefficients for time and cost.

## Dependencies

Please install R (version 4.1 or newer is recommended) plus the following
packages before running the scripts:

- `apollo`
- `data.table`
- `dplyr`
- `readr`
- `tidyr`
- `utils`

Additional tidyverse packages may already be available on your system, but the
list above covers all packages explicitly used by the repository scripts.

## File structure

- `R00_Data Processing.R` – Reads the original wide-format data, cleans it, and
  produces the long-format dataset used for estimation.
- `R_01_MNL Base Line.R` – Estimates the baseline multinomial logit (MNL)
  model.
- `R_02_MMNL Independent.R` – Estimates the mixed logit model with independent
  random coefficients, including scaling and diagnostics.
- `run_model.R` – Orchestrator that checks packages, runs preprocessing, and
  launches the MMNL estimation from a single entry point.
- `R_03_MMNL Dependent.R` – **Planned** future extension for a mixed logit model
  with correlated random coefficients (script not yet implemented).
- `apollo_modeChoiceData.csv` – Raw source data in wide format.
- `modechoice_long.csv` / `DATA/processed/modechoice_long.csv` – Saved
  long-format dataset produced by the preprocessing script.
- `output/` – Folder created at runtime to store Apollo estimation logs and
  summaries.

## How to run

### Run the baseline MNL model directly

```bash
Rscript 'R_01_MNL Base Line.R'
```

### Run the MMNL (independent) model directly

Ensure the processed long-format data exists (either by running the data
processing script first or via the full pipeline below), then execute:

```bash
Rscript 'R_02_MMNL Independent.R'
```

### Run the full preprocessing + MMNL pipeline

From the repository root, execute:

```bash
Rscript run_model.R
```

The runner will copy the raw CSV into `DATA/`, generate the long-format dataset
with `R00_Data Processing.R`, and estimate the MMNL model with
`R_02_MMNL Independent.R`. Model outputs, diagnostic summaries, and Apollo log
files are written to the `output/` directory.

### Running from RStudio

1. Open the `mixedlogit-debug` project folder in RStudio.
2. Source `run_model.R` to build the processed data and estimate the MMNL model
   in one step; or open `R_02_MMNL Independent.R` and run `mmnl_model()`
   manually after ensuring `DATA/processed/modechoice_long.csv` exists.
3. Baseline checks: you can also source `R_01_MNL Base Line.R` directly to
   verify the non-mixed specification on the same processed dataset.


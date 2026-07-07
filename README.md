# An Analytical Framework for Public Debt Management

Replication package for Jean-Paul Renne, *An Analytical Framework for Public Debt Management*.

## Contents

* `main.R`: master replication script.
* `Data/data.Rda`: frozen input dataset used by the default replication run.
* `estimation/`: scripts that load data, construct targets, and estimate/calibrate the model.
* `procedures/`: shared R and C++ routines used by the model solution and simulations.
* `simulations/`: scripts for the demand/supply exercise and issuance-strategy analysis.
* `outputs/`: scripts that create paper figures and LaTeX table fragments.
* `resources/`: small static resources copied into the output folders.
* `results/`: saved model estimates and strategy-simulation results.
* `figures/`: output folder for replicated figures.
* `tables/`: output folder for replicated LaTeX tables.

## Quick Replication

Open R from the repository root and run:

```r
source("main.R")
```

Or run the same script from a terminal:

```sh
Rscript main.R
```

The default run loads the saved calibrated model in `results/res_26082024.Rdat` and regenerates:

* `figures/Figure_fit.pdf`
* `figures/Figure_avg_yc.pdf`
* `figures/formula.pdf`
* `tables/table_param.txt`
* `tables/table_moment_matching.txt`

This is the fastest way to check that the package is installed correctly.

## Optional Replication Blocks

The main switches are near the top of `main.R`.

Set:

```r
indic_DemSup <- TRUE
```

to regenerate the demand/supply exercise outputs:

* `figures/Figure_expected_returns_DemaSupp.pdf`
* `figures/Figure_nu_effect_DemaSupp.pdf`
* `tables/table_param_DemSup.txt`
* `tables/table_DemSup_elastsurplus*_nu*.txt`

Set:

```r
indic_run_performances <- TRUE
```

to regenerate the issuance-strategy simulation outputs:

* `results/results_strategies.Rda`
* `figures/Figure_strategies_perf.pdf`
* `tables/table_strategies.txt`

Set:

```r
indic_estim <- TRUE
```

to re-estimate the model. This is slower than the default run. If `indic_save_model <- TRUE`, the final model is saved to the file named by `file_with_saved_param`.

## Rebuilding The Data

The default replication uses the frozen dataset:

* `Data/data.Rda`

To rebuild the dataset from public sources, set:

```r
indic_load_data <- TRUE
```

near the top of `main.R`.

The data-refresh script downloads CPI and real GDP from FRED and nominal and real yield-curve data from the Federal Reserve Board. FRED access requires an API key. Set it before launching R:

```sh
export FRED_API_KEY="your_fred_key"
```

The rebuilt dataset is saved as `Data/data.Rda`. Public macroeconomic data are revised over time, so an updated-data run may not reproduce the frozen paper input exactly.

Large raw Federal Reserve yield-curve CSV files are not stored in the repository; `estimation/load_data.R` downloads them when the data are rebuilt. The default replication uses the compact frozen file `Data/data.Rda`.

The illustrative default-probability chart `figures/formula.pdf` is copied from `resources/formula.pdf`, which stores the publication-quality TikZ/LaTeX version used in the paper.

## R Package Dependencies

The replication code uses:

* `fredr`
* `Hmisc`
* `Rcpp`
* `RcppEigen`
* `optimx`
* `tidyverse`
* `zoo`

Install missing packages with, for example:

```r
install.packages(c("fredr", "Hmisc", "Rcpp", "RcppEigen", "optimx", "tidyverse", "zoo"))
```

The C++ routines in `procedures/library_cpp.cpp` are compiled automatically by `Rcpp::sourceCpp()` when `main.R` is run.

## Paper Integration

The LaTeX paper uses copies of the generated files in the Overleaf project folders `Figures/` and `Tables/`. After regenerating outputs in this repository, copy the updated files from `figures/` and `tables/` to the paper project before compiling the manuscript.

## License

The original code in this repository is released under the MIT License; see `LICENSE`. Data downloaded by the scripts remain subject to the terms of their original providers.

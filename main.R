# ==============================================================================
# AN ANALYTICAL FRAMEWORK FOR PUBLIC DEBT MANAGEMENT
# ==============================================================================
# Jean-Paul RENNE
# This version: July 2026.
# corresponding author: jean-paul.renne@unil.ch
# ==============================================================================

# This is the master replication script. Run it from the root of the replication
# package, e.g.:
#   setwd("/path/to/PDMAnalyt")
#   source("main.R")
# or, from a terminal:
#   Rscript main.R
#
# The default settings below reproduce the calibration tables and figures using
# the saved model parameters in results/. The estimation and strategy simulations
# are deliberately off by default because they are substantially slower.

# Clear the workspace to avoid accidental dependence on objects from an
# interactive R session.
rm(list = ls())

# --- Replication switches ------------------------------------------------------
# Set these flags to 1 to regenerate the corresponding groups of results.

# Demand/supply exercise: creates the stylized-economy tables and figures used
# in the demand-vs-supply discussion and appendix.
indic_DemSup <- FALSE

# Full model estimation: re-estimates the macro-finance model from the starting
# values specified below. Leave at 0 to load saved parameters from results/.
indic_estim <- TRUE

# Issuance-strategy simulations: computes performance statistics over the grid
# of maturities and indexation choices. Set to TRUE for a full replication.
indic_run_performances <- TRUE
# ------------------------------------------------------------------------------



# --- Estimation/data settings --------------------------------------------------
start_year  <- 1970 # first year of the estimation sample

# Saved model loaded when indic_estim is FALSE, and overwritten when
# indic_estim and indic_save_model are TRUE.
file_with_saved_param <- "res_26082024.Rdat"

# Starting point loaded when indic_estim and indic_use_last_res are TRUE.
file_with_input_param <- "res_24082024.Rdat"

# Refresh the raw dataset from FRED and the Federal Reserve Board. This requires
# internet access and a valid FRED API key in estimation/load_data.R.
indic_load_data    <- FALSE

# If re-estimating, start from file_with_input_param rather than from the generic
# initialization in estimation/set_ini_model.R.
indic_use_last_res <- TRUE

# Standard deviation of the random perturbation applied to starting parameters
# in each estimation attempt. The seed is fixed below for reproducibility.
random_factor      <- .1

# If re-estimating, save the final calibrated model to file_with_saved_param.
indic_save_model   <- TRUE
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Load packages and shared routines.
library(fredr)
library(Hmisc)
library(Rcpp)
library(RcppEigen)
library(optimx)

source("procedures/proc_model.R")
sourceCpp("procedures/library_cpp.cpp")
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# ---- Demand/Supply exercise --------------------------------------------------
if(isTRUE(indic_DemSup)){
  message("")
  message("------------------------------------------------")
  message(" Demand/Supply exercise")
  message("------------------------------------------------")
  
  # Numerical settings for the demand/supply exercise.
  maxH        <- 10 # maximum maturity considered for zero-coupon bonds
  nb_iter     <- 30 # iterations used to solve models
  nb_iter_sdf <- 10 # iterations used to solve SDF
  nb_grid     <- 25 # size of grids
  
  # min and max values of d and r:
  min_d  <- .4
  max_d  <- 1.6
  min_rr <- .0
  max_rr <- .15
  
  # Stylized macro states used to build the demand- and supply-driven economies.
  DemSup_low_pi <- .0
  DemSup_med_pi <- .03
  DemSup_hig_pi <- .06
  DemSup_low_y  <- .0
  DemSup_med_y  <- .02
  DemSup_hig_y  <- .04
  DemSup_rho    <- .8
  DemSup_d_bar  <- .8

  # Scenarios covered by the demand/supply tables and figures.
  DemSup_elasticities_of_surpluses <- c(0,1)
  DemSup_abs_nu_y_values          <- c(0,.1)
  DemSup_plot_nu_y_values         <- c(0,-.1)
  values_of_chi <- c(.2,.9)
  
  for(elasticity_of_surpluses in DemSup_elasticities_of_surpluses){
    for(abs_nu_y in DemSup_abs_nu_y_values){
      source("simulations/exercise_demand_supply.R")
    }
  }
  
  # Output files:
  #   tables/table_param_DemSup.txt
  #   tables/table_DemSup_elastsurplus*_nu*.txt
  #   figures/Figure_expected_returns_DemaSupp.pdf
  #   figures/Figure_nu_effect_DemaSupp.pdf
  source("outputs/make_table_figure_DemSup.R")
}
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# ---- Data loading and model calibration --------------------------------------

if(isTRUE(indic_load_data)){
  message(" --- Loading data ---")
  source("estimation/load_data.R")
  message(" --- Loading data: Done ---")
}else{
  load(file="Data/data.Rda")
}

# Keep the estimation sample starting in start_year and construct the moment
# targets used by the loss function and by the moment-matching table.
source("estimation/resize_and_compute_targets.R")

if(isTRUE(indic_estim)){
  
  set.seed(123) # for the exact replication of the results
  # Note: random numbers are used only to randomize starting values
  # in the loss function optimization.

  # Number of Markov regimes in the macro block.
  nb_m <- 5
  
  message("")
  message("------------------------------------------------")
  message(" Calibration of macro block")
  message("------------------------------------------------")
  
  # Parameter constraints:
  min_Pi <- -.02
  max_Pi <- +.20
  min_Dy <- -.10
  max_Dy <- +.06
  min_gamma <- 1
  max_gamma <- 10
  # Maximum value of log(param):
  max_abs_param <- 8
  
  # Optimization setup for the macro block.
  maxit.nlminb <- 50
  maxit.NlMd   <- 2000
  nb_loops     <- 2
  
  nb_attemps <- 40
  # Run a number nb_attemps of new estimations, starting from randomized
  # starting values:
  best <- 100000
  for(iii in 1:nb_attemps){
    source("estimation/run_estim.R")
    if(res.estim["value"]<best){
      message("--- new best model ---")
      best_param <- param
      best <- res.estim["value"]
    }
  }
  # Reload best param:
  param <- best_param
  Model <- make_model(param,Model_ini)
  
  message("")
  message("------------------------------------------------")
  message(" Calibration of nu_y, nu_pi, and mu_eta")
  message("------------------------------------------------")
  
  # Parameters controlling the macroeconomic effect of default.
  Model$nu_y  <- -.05
  Model$nu_pi <- -.021

  # Mean growth/inflation dynamics after default.
  Model$mu_eta <- .5 * Model$mu_y
  
  message("")
  message("------------------------------------------------")
  message(" Calibration of alpha, beta, d_star and s_star")
  message("------------------------------------------------")
  
  candidate_alpha_values  <- c(.1,.2)
  candidate_beta_values   <- c(.02,.05,.1,.2)
  candidate_d_star_values <- c(.9,1,1.1,1.2)
  
  nb_grid         <- 25 # number of values per state variable
  nb_iter         <- 30 # iterations to solve model
  nb_iter_sdf     <- 10 # iterations to solve SDF
  
  # min and max values of d and r:
  min_d  <- .4
  max_d  <- 1.6
  min_rr <- .0
  max_rr <- .15
  
  avgD <- 6 # targeted average debt maturity
  
  # Define targets:
  Targets <- list(spread_in_bps = 20,
                  mean_d_in_percent = 80,
                  stdv_d_in_percent = 15)
  
  # Grid search over fiscal/default parameters. The selected model is saved
# below when indic_save_model is TRUE.
  source("estimation/calibrate_alpha_beta.R")
  
  if(isTRUE(indic_save_model)){
    save(Model,file=paste("results/",file_with_saved_param,sep=""))
  }
  
}else{
  # Fast replication path: load the saved model parameters.
  load(file=paste("results/",file_with_saved_param,sep=""))
}

# Default outputs produced by every run:
#   figures/formula.pdf
#   figures/Figure_fit.pdf
#   figures/Figure_avg_yc.pdf
#   tables/table_moment_matching.txt
#   tables/table_param.txt
source("outputs/make_figure_illustrPD.R")
source("outputs/make_figures_fit.R")
source("outputs/make_table_moments.R")
source("outputs/make_table_param.R")
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# ---- Performance of issuance strategies --------------------------------------

# Strategy-performance settings. These are the main controls for the long run.
nb_grid         <- 30   # number of values per state variable
nb_iter         <- 30   # iterations used to solve the perpetuity price
nb_iter_sdf     <- 10   # iterations used to solve the SDF
nb_iter4probas  <- 1000 # iterations used to compute unconditional distributions
maxH            <- 10   # maximum maturity of zero-coupon bonds

# FAST ==================================
# FAST ==================================
# FAST ==================================
# FAST ==================================
nb_grid         <- 25   # number of values per state variable
nb_iter         <- 20   # iterations used to solve the perpetuity price
nb_iter_sdf     <- 10   # iterations used to solve the SDF
nb_iter4probas  <- 500 # iterations used to compute unconditional distributions
maxH            <- 10   # maximum maturity of zero-coupon bonds
# FAST ==================================
# FAST ==================================
# FAST ==================================
# FAST ==================================


# Parallelize independent strategy evaluations. This uses forked processes on
# macOS/Linux. Progress messages are shown, but may arrive out of order because
# workers finish at different times; set to FALSE for strictly sequential logs.
indic_parallel_strategies <- TRUE
nb_cores_strategies <- min(8, max(1, parallel::detectCores() - 1))

# Set to TRUE only when debugging a single strategy interactively. Keeping these
# plots off makes the batch run faster and avoids refreshing Rplots.pdf.
indic_strategy_diagnostic_plots <- FALSE

# min and max values of d and r:
min_d  <- .4
max_d  <- 1.6
min_rr <- .0
max_rr <- .15

# Determine strategies to explore:
values_of_chi      <- c(.1,.5,.9)
values_of_kappa_pi <- seq(0,1,by=.25)
values_of_kappa_y  <- seq(0,.3,by=.1)
values_of_chi      <- sort(unique(values_of_chi))
values_of_kappa_pi <- sort(unique(values_of_kappa_pi))
values_of_kappa_y  <- sort(unique(values_of_kappa_y))

# Add the strategy approximating the US government debt portfolio as one extra
# point, rather than adding each US parameter to the full Cartesian grid.
additional_strategies <- data.frame(chi = Model$chi,
                                    kappa_pi = 1935/25734,
                                    kappa_y = 0)

# Define variables to be plotted on scatter plots:
outputs4chart <- matrix(NaN,3,2)
outputs4chart[,1] <- "mean_d"
outputs4chart[,2] <- c("stdv_d","DaR95","avg_PD[maxH]")

if(isTRUE(indic_run_performances)){
  message("")
  message("------------------------------------------------")
  message(" Compute performances of issuance strategies")
  message("------------------------------------------------")
  
  # Output files:
  #   results/results_strategies.Rda
  #   figures/Figure_strategies_perf.pdf
  #   tables/table_strategies.txt
  source("simulations/run_strategies.R")  
  source("outputs/make_figure_strategies.R")  
  source("outputs/make_table_strategies.R")  
}
# ------------------------------------------------------------------------------

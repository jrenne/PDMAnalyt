# ==============================================================================
# AN ANALYTICAL FRAMEWORK FOR PUBLIC DEBT MANAGEMENT
# ==============================================================================
# Jean-Paul RENNE
# This version: August 2024.
# corresponding author: jean-paul.renne@unil.ch
# ==============================================================================

# clear workspace:
rm(list = ls())

# -- Run simulation in demand/supply economies? --------------------------------
indic_DemSup <- 0
# -- Run new estimation? -------------------------------------------------------
indic_estim <- 0 # binary variable -- if 1, run the whole estimation approach
# -- Run simulation in calibrated economy? -------------------------------------
indic_run_performances <- 0



# --- if estimation: -----------------------------------------------------------
start_year  <- 1970 # first year of the estimation sample
file_with_saved_param <- "res_26082024.Rdat"
file_with_input_param <- "res_24082024.Rdat"
indic_load_data    <- 0 # update dataset (from FRED and Board for yields)
indic_use_last_res <- 1 # binary variables -- start from last best param or not
random_factor      <- .1 # randomization to change initial values
indic_save_model   <- 1 # binary variables -- say if estimated model is saved
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Load packages:
library(fredr)
library(Hmisc)
library(Rcpp)
library(RcppEigen)
library(optimx)

# Load procedures:
source("procedures/proc_model.R")
sourceCpp("procedures/library_cpp.cpp") # load C++ functions
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# ---- Demand/Supply exercise --------------------------------------------------
if(indic_DemSup){
  print("")
  print("------------------------------------------------")
  print(" Demand/Supply exercise")
  print("------------------------------------------------")
  
  maxH        <- 10 # maximum maturity considered for zero-coupon bonds
  nb_iter     <- 30 # iterations used to solve models
  nb_iter_sdf <- 10 # iterations used to solve SDF
  nb_grid     <- 25 # size of grids
  
  # min and max values of d and r:
  min_d  <- .4
  max_d  <- 1.6
  min_rr <- .0
  max_rr <- .15
  
  elasticity_of_surpluses <- 1
  #abs_nu_y <- .1
  abs_nu_y <- 0
  
  values_of_chi <- c(.2,.9)
  
  for(elasticity_of_surpluses in c(0,1)){
    for(abs_nu_y in c(0,.1)){
      source("simulations/exercise_demand_supply.R")
    }
  }
  
  # Prepare table showing parameterization and figure showing yield curves:
  source("outputs/make_table_figure_DemSup.R")
}
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# ---- Model calibration -------------------------------------------------------

if(indic_load_data==1){# in that case refresh data (from the internet)
  print(" --- Loading data ---")
  source("estimation/load_data.R")
  print(" --- Loading data: Done ---")
}else{
  load(file="Data/data.Rda")
}
# Resize dataset (based on "start_year"), and compute moment targets:
source("estimation/resize_and_compute_targets.R")

if(indic_estim == 1){
  
  set.seed(123) # for the exact replication of the results
  # Note: random numbers are used only to randomize starting values
  # in the loss function optimization.

  # Number of regimes:
  nb_m <- 5
  
  print("")
  print("------------------------------------------------")
  print(" Calibration of macro block")
  print("------------------------------------------------")
  
  # Parameter constraints:
  min_Pi <- -.02
  max_Pi <- +.20
  min_Dy <- -.10
  max_Dy <- +.06
  min_gamma <- 1
  max_gamma <- 10
  # Maximum value of log(param):
  max_abs_param <- 8
  
  # Optimization set up:
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
      print("--- new best model ---")
      best_param <- param
      best <- res.estim["value"]
    }
  }
  # Reload best param:
  param <- best_param
  Model <- make_model(param,Model_ini)
  
  print("")
  print("------------------------------------------------")
  print(" Calibration of nu_y, nu_pi, and mu_eta")
  print("------------------------------------------------")
  
  # Calibration of nu_y and nu_pi
  Model$nu_y  <- -.05
  Model$nu_pi <- -.021
  # Modify mu_eta:
  Model$mu_eta <- .5 * Model$mu_y
  
  print("")
  print("------------------------------------------------")
  print(" Calibration of alpha, beta, d_star and s_star")
  print("------------------------------------------------")
  
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
  
  source("estimation/calibrate_alpha_beta.R")
  
  if(indic_save_model){
    save(Model,file=paste("results/",file_with_saved_param,sep=""))
  }
  
}else{
  load(file=paste("results/",file_with_saved_param,sep=""))
}

# Prepare figures illustrating the model fit:
source("outputs/make_figures_fit.R")
source("outputs/make_table_moments.R")
source("outputs/make_table_param.R")
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# ---- Analysis of the performances of issuance strategies ---------------------

# Specify solution approach:
nb_grid     <- 30 # number of values per state variable
nb_iter     <- 30 # iterations used to solve model
nb_iter_sdf <- 10 # iterations used to solve SDF
maxH        <- 10 # maximum maturity of zero-coupon bonds

# min and max values of d and r:
min_d  <- .4
max_d  <- 1.6
min_rr <- .0
max_rr <- .15

# Determine strategies to explore:
values_of_chi      <- c(.1,.5,.9)
values_of_kappa_pi <- seq(0,1,by=.2)
values_of_kappa_y  <- seq(0,.3,by=.1)

# Define variables to be plotted on scatter plots:
outputs4chart <- matrix(NaN,3,2)
outputs4chart[,1] <- "mean_d"
outputs4chart[,2] <- c("stdv_d","DaR95","avg_PD[maxH]")

if(indic_run_performances){
  print("")
  print("------------------------------------------------")
  print(" Compute performances of issuance strategies")
  print("------------------------------------------------")
  
  source("simulations/run_strategies.R")  
  source("outputs/make_figure_strategies.R")  
  source("outputs/make_table_strategies.R")  
}
# ------------------------------------------------------------------------------


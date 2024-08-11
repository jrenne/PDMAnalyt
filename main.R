# ==============================================================================
# ANALYTICAL FRAMEWORK FOR PUBLIC DEBT MANAGEMENT
# ==============================================================================
# Jean-Paul RENNE
# This version: July 2024.
# corresponding author: jean-paul.renne@unil.ch
# ==============================================================================

# clear workspace:
rm(list = ls())

# Number of cores to be used (relevant only for estimation):
number.of.cores <- 8

indic_estim     <- 1
indic_load_data <- 0

# # Load packages:
# library(tsDyn)
# library(vars)
# library(mnormt)
# library(Hmisc)
# library(expm)
# library(optimx)
# library(tikzDevice)
# library(stringr)
# library(doParallel)
# library(seasonal)
# library(mFilter)
# library(splines)
# library(zoo)
library(Rcpp)
library(RcppEigen)
library(optimx)
# library(pracma)

# Load procedures:
source("procedures/proc_toymodel.R")
sourceCpp("procedures/pricing_cpp.cpp") # load C++ functions

# ---- Model calibration -------------------------------------------------------
if(indic_estim == 1){
  if(indic_load_data==1){
    print(" --- Loading data ---")
    source("estimation/load_data.R")
    print(" --- Loading data: Done ---")
  }else{
    load(file="Data/data.Rda")
  }
  start_year <- "1968"
  source("estimation/run_estim.R")
}else{
  load(file="results/res_882024.Rdat")
}

# ---- Test PDM framework ------------------------------------------------------
nb_iter <- 10 # number of iterations used to solve for perpetuity yields
maxH    <- 20
source("make_grids.R")


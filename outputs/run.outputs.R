# ==============================================================================
# This script calls the scripts producing tables and figures
# ==============================================================================

parcex    <- 1.25 # not for prez
parcexleg <- 1.25 # size of legend font

# Forecasts:
print("=== Preparing forecast figures ===")
max.nb.rows.per.figure <- 8
max.nb.of.fcsts <- 5 # Maximum number of fcsts in DATA (across countries)
for(variable in c("debt","infl","growth")){
  source(paste("outputs/make.charts.fcsts.R", sep=""))
}

# Correlation between risk-free yields and CDS spreads
print("=== Preparing figure showing yields vs CDS spreads ===")
source(paste("outputs/make.charts.CorrelYDS_CDS.R", sep=""))

# Yields
print("=== Preparing bond yields figures ===")
max.nb.of.yds <- 3 # Maximum number of yields in DATA (across countries)
source(paste("outputs/make.charts.yds.R", sep=""))

# CDS Spreads
print("=== Preparing CDS figures ===")
nb.of.CDS.maturities <- 5 # maximum number of CDS maturities in DATA (across countries)
source(paste("outputs/make.charts.CDS.R", sep=""))

# Default Probabilities
print("=== Preparing Default Proba figures ===")
Horizons <- c(8,12,20)
Horizons <- c(8,20)
source(paste("outputs/make.charts.ProbaDef.R", sep=""))

#max.nb.rows.per.figure <- 4

# Fiscal limits:
print("=== Preparing Fiscal Limit figures ===")
#Proba.4.def.of.plotted.FL <- c(.01,.05)
#nb.Q <- 4
source("outputs/make.charts.FL.R")

# Table of parameter estimates (and stdv of parameters)
print("=== Preparing model parameterization table ===")
source("outputs/make.tables.param.paper.R")

print("=== Run simulations ===")
if(indic.comput.stdv.param==1){
  N.replic      <- 100
  max.diff.logl <- 300
  source("outputs/run.simulation.R")
}

print("=== Preparing Fiscal Space figures ===")
source("outputs/make.charts.FS.R")

print("=== Preparing S_star figures ===")
source("outputs/make.charts.Sstar.R")

print("=== Preparing w factors figures ===")
source("outputs/make.charts.factors.R")

# Exercise illustrating CDS sensitivity:
print("=== Preparing CDS sensitivity figures ===")
source("outputs/make.table.figure.sensitivity.R")


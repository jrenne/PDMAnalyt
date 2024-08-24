
# # Optimization set up:
# maxit.nlminb <- 100
# maxit.NlMd   <- 4000
# nb_loops     <- 2


start_year <- 1970
#start_year <- 2000

maxH = 10 # for charts

# Resize dataset:
indic_fst <- which(format(DATA$date,"%Y")==start_year)
DATA <- DATA[indic_fst:dim(DATA)[1],]
#DATA$surpl_1 <- c(DATA$surpl[2:length(DATA$surpl)],NaN)

# Number of regimes:
nb_m <- 5

min_Pi <- -.02
max_Pi <- +.20
min_Dy <- -.10
max_Dy <- +.06
min_gamma <- 1
max_gamma <- 10

max_abs_param <- 8

# Targetted moments:
targets <- list(
  target_10_nom = mean(DATA$SVENY10,na.rm=TRUE)/100,
  target_01_nom = mean(DATA$SVENY01,na.rm=TRUE)/100,
  target_10_rea = mean(DATA$TIPSY10,na.rm=TRUE)/100,
  target_02_rea = mean(DATA$TIPSY02,na.rm=TRUE)/100,
  target_slop_nom = mean(DATA$SVENY10,na.rm=TRUE)/100 - 
    mean(DATA$SVENY01,na.rm=TRUE)/100,
  target_slop_rea = mean(DATA$TIPSY10,na.rm=TRUE)/100 -
    mean(DATA$TIPSY02,na.rm=TRUE)/100,
  target_avg_Pi = mean(DATA$pi,na.rm=TRUE),
  target_avg_Dy = mean(DATA$dy,na.rm=TRUE),
  target_std_10_nom = sd(DATA$SVENY10,na.rm=TRUE)/100,
  target_std_10_rea = sd(DATA$TIPSY10,na.rm=TRUE)/100,
  IRP10 = 0 # 10-year inflation risk premium
  )

# targets <- list(
#   target_10_nom = .06,
#   target_01_nom = .045,
#   target_10_rea = .015,
#   target_01_rea = .01,
#   target_slop_nom = .015,
#   target_slop_rea = .005,
#   target_avg_Pi = .03,
#   target_avg_Dy = .015)

# Set initial Model ------------------------------------------------------------
gamma      <- 6      # risk aversion
delta      <- .99    # preference for present
beta       <- .3     # sensitivity of surplus to debt level
d_star     <- .6     # targeted debt level
sigma_eps  <- .02    # std dev of surplus shocks
alpha      <- 1      # elasticity of default proba
s_star     <- .06    # max surplus
std_Pi     <- .01    # standard deviation of inflation measurement errors (for Hamilton filter)
std_Dy     <- .005   # standard deviation of GDP growth measurement errors (for Hamilton filter)
std_nom_yd <- .008   # standard deviation of yield measurement errors (for Hamilton filter)
nu_y       <- -.2 * 0
nu_pi      <- .05 * 0
sigma_nu   <- .10    # uncertainty regarding debt threshold

# Perpetuity specification:
kappa_pi <- 0 # Inflation indexation (1 = TIPS)
kappa_y  <- 0 # Real GDP indexation
RR       <- .3 # Recovery rate
chi      <- .7 # coupon decay

rho <- .5
Omega <- rho * diag(nb_m)
Omega[1:(nb_m-1),2:nb_m] <- Omega[1:(nb_m-1),2:nb_m] + (1 - rho)/2 * diag(nb_m-1)
Omega[2:nb_m,1:(nb_m-1)] <- Omega[2:nb_m,1:(nb_m-1)] + (1 - rho)/2 * diag(nb_m-1)
Omega[1,2] <- 1 - Omega[1,1]
Omega[nb_m,nb_m-1] <- 1 - Omega[nb_m,nb_m]
Omega <- Omega + .01 * matrix(1,nb_m,nb_m)
Omega <- Omega / (matrix(apply(Omega,1,sum),ncol=1) %*% matrix(1,1,nb_m))

mu_y  <- matrix(c(min(DATA$dy,na.rm=TRUE),rep(mean(DATA$dy,na.rm=TRUE),nb_m-2),.05),nb_m,1)
mu_pi <- matrix(c(rep(mean(DATA$pi,na.rm=TRUE),nb_m-1),min(DATA$pi,na.rm = TRUE)),nb_m,1)

mu_eta <- .5 * mu_y

Model_ini <- list(gamma=1.2,
                  delta     = delta,
                  chi       = chi,
                  beta      = beta,
                  d_star    = d_star,
                  sigma_eps = sigma_eps,
                  sigma_nu  = sigma_nu,
                  alpha     = alpha,
                  s_star    = s_star,
                  RR        = RR,
                  Omega     = Omega,
                  mu_pi     = mu_pi,
                  mu_y      = mu_y,
                  mu_eta    = mu_eta,
                  nu_pi     = nu_pi,
                  nu_y      = nu_y,
                  kappa_pi  = kappa_pi,
                  kappa_y   = kappa_y,
                  std_Pi    = std_Pi,
                  std_Dy    = std_Dy,
                  std_nom_yd = std_nom_yd)

load(file="results/res_22082024.Rdat")
Model_ini <- Model


param <- model2param(Model_ini)
param <- param + .4*rnorm(param)

#param <- param + .1*rnorm(param)

algo <- "Nelder-Mead"
#algo <- "nlminb"
for(ii in 1:nb_loops){
  for(algo in c("Nelder-Mead","nlminb")){
    res.estim <- optimx(param,
                        #compute_distance,
                        #compute_loglik,
                        compute_total_distance,
                        targets = targets,
                        Model_ini = Model_ini,
                        method = algo,
                        control=
                          list(trace=1,
                               maxit=ifelse(algo=="Nelder-Mead",
                                            maxit.NlMd,maxit.nlminb),
                               kkt = FALSE))
    param     <- c(as.matrix(res.estim)[1:length(param)])
  }
}

compute_total_distance(param,targets,Model_ini)

Model <- make_model(param,Model_ini)

#save(Model,file="results/res_772024.Rdat")
#save(Model,file="results/res_882024.Rdat")
#save(Model,file="results/res_11082024.Rdat")
#save(Model,file="results/res_12082024.Rdat")
#save(Model,file="results/res_13082024.Rdat")
#save(Model,file="results/res_14082024.Rdat")
#save(Model,file="results/res_16082024.Rdat")
#save(Model,file="results/res_22082024.Rdat")

#load(file="results/res_22082024.Rdat")

source("make_chart_fit.R")



nb_m <- 5

maxit.nlminb <- 100
maxit.NlMd   <- 5000
nb_loops     <- 3

std_Pi <- .01
std_Dy <- .01
std_nom_yd <- .01

min_Pi <- -.02
max_Pi <- +.20
min_Dy <- -.10
max_Dy <- +.06
min_gamma <- 1
max_gamma <- 10

max_abs_param <- 10

targets <- list(
  target_10_nom = .06,
  target_01_nom = .045,
  target_10_rea = .015,
  target_01_rea = .01,
  target_slop_nom = .015,
  target_slop_rea = .005,
  target_avg_Pi = .03,
  target_avg_Dy = .015)



# Set initial Model ------------------------------------------------------------
gamma     <- 5       # risk aversion
delta     <- .99     # preference for present
beta      <- .3      # sensitivity of surplus to debt level
d_star    <- .6      # targeted debt level
sigma_eps <- .04     # std dev of surplus shocks
alpha     <- 1     # elasticity of default proba
s_star    <- .03     # max surplus
nu_y      <- -.2 * 0
nu_pi     <- .05 * 0

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

Model <- list(gamma=1.2,
              delta     = delta,
              chi       = chi,
              beta      = beta,
              d_star    = d_star,
              sigma_eps = sigma_eps,
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
              kappa_y   = kappa_y)

param <- model2param(Model)
param <- param + .6*rnorm(param)

library(optimx)
algo <- "Nelder-Mead"
#algo <- "nlminb"
for(ii in 1:nb_loops){
  for(algo in c("Nelder-Mead","nlminb")){
    res.estim <- optimx(param,
                        #compute_distance,
                        #compute_loglik,
                        compute_total_distance,
                        targets = targets,
                        method = algo,
                        control=
                          list(trace=1,
                               maxit=ifelse(algo=="Nelder-Mead",
                                            maxit.NlMd,maxit.nlminb),
                               kkt = FALSE))
    param     <- c(as.matrix(res.estim)[1:length(param)])
  }
}
#save(Model,file="results/res_772024.Rdat")


Model <- make_model(param)

# Compute (LT) bond prices:
Model_nominal <- Model
Model_nominal$kappa_pi <- 0
Model_nominal$kappa_y  <- 0
res_LTnominal_prices <- compute_LTRF_bond_prices(Model_nominal,maxH=10)

yds <- compute_LT_yds(Model)
yds_nom <- yds$yds_nom
yds_rea <- yds$yds_rea

stat_distri <- yds$stat_distri
avg_Pi <- sum(stat_distri * Model$mu_pi)
avg_Dy <- sum(stat_distri * Model$mu_Dy)

print(cbind(Model$mu_pi,Model$mu_y))
print(Model$Omega)

N <- rbind(rep(std_Pi,nb_m),
           rep(std_Dy,nb_m),
           rep(std_nom_yd,nb_m),
           rep(std_nom_yd,nb_m))
M <- rbind(c(Model$mu_pi),
           c(Model$mu_y),
           res_LTnominal_prices$all_LT_rth[,1],
           res_LTnominal_prices$all_LT_rth[,10])
F <- cbind(DATA$pi,
           DATA$dy,
           DATA$SVENY01/100,
           DATA$SVENY10/100)
F <- F[complete.cases(F),]
res_KH <- KH_filter(F,M,N,Model$Omega)

fitted <- res_KH$ksi_matrix %*% t(M)

par(mfrow=c(2,3))
for(i in 1:dim(M)[1]){
  plot(F[,i],type="l")
  lines(fitted[,i],col="red")
}

plot(yds_nom,type="l",ylim=c(0,.08))
lines(yds_rea,col="red")


# T <- 30
# res_simul <- simul_RS(Model,maxH=10,T)
# par(mfrow=c(2,1))
# plot(res_simul$Pi,type="l")
# plot(res_simul$Dy,type="l")
# 
# # Check KH filter formula:
# noise_Pi <- 3*c(.01,.02,.01)
# noise_Dy <- 3*c(.01,.01,.02)
# N <- rbind(noise_Pi,
#            noise_Dy)
# M <- rbind(Model$mu_pi,Model$mu_y)
# noisy_Pi <- res_simul$Pi + res_simul$states %*% noise_Pi * rnorm(T)
# noisy_Dy <- res_simul$Dy + res_simul$states %*% noise_Dy * rnorm(T)
# plot(noisy_Pi,type="l")
# plot(noisy_Dy,type="l")
# 
# F <- cbind(noisy_Pi,noisy_Dy)
# 
# res_KH <- KH_filter(F,M,N,Model$Omega)
# 
# par(mfrow = c(1,nb_m))
# for(i in 1:nb_m){
#   plot(res_simul$states[,i])
#   lines(res_KH$ksi_matrix[,i],col="red")
# }
# 







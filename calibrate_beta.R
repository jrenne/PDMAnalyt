# ==============================================================================
# Calibration of beta
# ==============================================================================

candidate_values <- c(.05,.1,.2)

# Issuance of nominal bonds:
Model$kappa_pi <- 0
Model$kappa_y  <- 0

nb_iter <- 30 # used to solve model

res_stat_distri_and_rbar <- compute_stat_distri_and_rbar(Model)
r_bar  <- res_stat_distri_and_rbar$r_bar
max_rr <- max((r_bar + .1) * Model$d_star,.15)
all_rr <- seq(0,max_rr,length.out = nb_grid)
all_rr <- matrix(all_rr,ncol=1)

# Determine Chi by fitting average debt maturity:
res_LTprices_nominal <- compute_LTRF_bond_prices(Model,maxH = 10)
stat_distri <- res_stat_distri_and_rbar$stat_distri
q <- c(t(stat_distri) %*% res_LTnominal_prices$all_LT_rth[,10])
D <- 6
Model$chi <- 1 + q - 1/D

# Values of debt and debt service at which we compute spread:
DD <- 1.0
indic_DD <- which.min(abs(all_d - DD))
rr <- .03
indic_rr <- which.min(abs(all_rr - rr))

# Consider a model without default risk (to determine debt2GDP distribuiton):
Model$beta <- .2
Model$d_star <- .4

Model_RF <- Model
Model_RF$alpha <- 0

res0_nominal <- solve_ToyModel(all_d,all_rr,all_eps,proba_eps,
                               Model_RF,nb_iter = nb_iter)
# Compute unconditional distribution of d:
p <- compute_uncond_distri(res0_nominal$indicators_x,res0_nominal$Probas,1000)
distri_d  <- compute_distri_x(all_d,res0_nominal$d,p)
plot(all_d,distri_d,type="l")
distri_rr  <- compute_distri_x(all_rr,res0_nominal$rr,p)
plot(all_rr,distri_rr,type="l")


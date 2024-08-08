

#Model$nu_y <- -.4

mean_pi <- .03
mean_y  <- .02

#Model$gamma <- .1
Model$RR    <- 0.3

chi <- .6

nb_iter <- 40 # used to solve the model

Model$mu_pi <- matrix(mean_pi,length(Model$mu_pi))
Model$mu_y  <- matrix(mean_y, length(Model$mu_pi))

Model$mu_eta <- 1 * Model$mu_y
Model$d_star <- .3
Model$alpha  <- .5
Model$beta   <- .1

# Issuance of nominal bonds:
Model$kappa_pi <- 0
Model$kappa_y  <- 0
res_stat_distri_and_rbar <- compute_stat_distri_and_rbar(Model)
r_bar  <- res_stat_distri_and_rbar$r_bar
max_rr <- (r_bar + .1) * Model$d_star
all_rr <- seq(0,max_rr,length.out = nb_grid)
all_rr <- matrix(all_rr,ncol=1)
Model$chi <- chi / exp((Model$kappa_pi - 1)*mean_pi + (Model$kappa_y - 1)*mean_y)
#Model$chi <- chi
res0_nominal <- solve_ToyModel(all_d,all_rr,all_eps,proba_eps,
                               Model,nb_iter = nb_iter)
res_LTprices_nominal <- compute_LTRF_bond_prices(Model,maxH)
avgLT_ExpReturns_nominal  <- c(t(res0_nominal$stat_distri) %*% res_LTprices_nominal$all_LT_ExpReturn_th)
PD_nominal <- compute_proba_def(maxH=maxH,
                                indicators_x  = res0_nominal$indicators_x,
                                all_proba_def = res0_nominal$all_proba_def,
                                Probas        = res0_nominal$Probas)

res_prices_nominal <- compute_bond_prices(Model, maxH,
                                          res0_nominal$indicators_x,
                                          res0_nominal$all_proba_def,
                                          res0_nominal$Probas)

Model_RF <- Model
Model_RF$RR <- 1
res_prices_nominal_RF <- compute_bond_prices(Model_RF, maxH,
                                             res0_nominal$indicators_x,
                                             res0_nominal$all_proba_def,
                                             res0_nominal$Probas)

plot(res_LTprices_nominal$all_LT_rth[1,],ylim=c(.05,.065))
lines(res_prices_nominal_RF$all_rth[1,])
lines(res_prices_nominal$all_rth[1,],col="red")

stop()

# Compute unconditional distribution of d:
p <- compute_uncond_distri(res0_nominal$indicators_x,res0_nominal$Probas,10000)
distri_d  <- compute_distri_x(all_d,res0_nominal$d,p)
plot(all_d,distri_d,type="l")
distri_Dy <- compute_distri_x(Model$mu_y,res0_nominal$Dy,p)
plot(Model$mu_y,distri_Dy)

# Issuance of ILBs:
Model$kappa_pi <- 1
Model$kappa_y  <- 0
res_stat_distri_and_rbar <- compute_stat_distri_and_rbar(Model)
r_bar  <- res_stat_distri_and_rbar$r_bar
max_rr <- (r_bar + .1) * Model$d_star
all_rr <- seq(0,max_rr,length.out = nb_grid)
all_rr <- matrix(all_rr,ncol=1)
Model$chi <- chi / exp((Model$kappa_pi - 1)*mean_pi + (Model$kappa_y - 1)*mean_y)
#Model$chi <- chi
res0_TIPS <- solve_ToyModel(all_d,all_rr,all_eps,proba_eps,
                            Model,nb_iter = nb_iter)
res_LTprices_TIPS <- compute_LTRF_bond_prices(Model,maxH)
avgLT_ExpReturns_TIPS  <- c(t(res0_TIPS$stat_distri) %*% res_LTprices_TIPS$all_LT_ExpReturn_th)
PD_TIPS <- compute_proba_def(maxH=maxH,
                             indicators_x  = res0_TIPS$indicators_x,
                             all_proba_def = res0_TIPS$all_proba_def,
                             Probas        = res0_TIPS$Probas)

# Issuance of GDP-linked bonds:
Model$kappa_pi <- 1
Model$kappa_y  <- 1
res_stat_distri_and_rbar <- compute_stat_distri_and_rbar(Model)
r_bar  <- res_stat_distri_and_rbar$r_bar
max_rr <- (r_bar + .1) * Model$d_star
all_rr <- seq(0,max_rr,length.out = nb_grid)
all_rr <- matrix(all_rr,ncol=1)
Model$chi <- chi / exp((Model$kappa_pi - 1)*mean_pi + (Model$kappa_y - 1)*mean_y)
#Model$chi <- chi
res0_GDPLB <- solve_ToyModel(all_d,all_rr,all_eps,proba_eps,
                             Model,nb_iter = nb_iter)
res_LTprices_GDPLB <- compute_LTRF_bond_prices(Model,maxH)
avgLT_ExpReturns_GDPLB  <- c(t(res0_GDPLB$stat_distri) %*% res_LTprices_GDPLB$all_LT_ExpReturn_th)
PD_GDPLB <- compute_proba_def(maxH=maxH,
                              indicators_x  = res0_GDPLB$indicators_x,
                              all_proba_def = res0_GDPLB$all_proba_def,
                              Probas        = res0_GDPLB$Probas)


par(mfrow=c(1,2))
v <- which(all_d==max(all_d))
v <- which(all_d==d_star)
v <- which(all_d==min(all_d))
v <- 16
Ts <- which((res0_nominal$d==all_d[v])&(res0_nominal$d_1==min(all_d))&(res0_nominal$rr==0))


plot(c(t(res0_nominal$stat_distri) %*% PD_nominal[Ts,]))
lines(c(t(res0_TIPS$stat_distri) %*% PD_TIPS[Ts,]),col="blue")
lines(c(t(res0_GDPLB$stat_distri) %*% PD_GDPLB[Ts,]),col="red")

# Ts <- Ts[3]
# plot(PD_nominal[Ts,])
# lines(PD_TIPS[Ts,],col="blue")
# lines(PD_GDPLB[Ts,],col="red")

plot(avgLT_ExpReturns_nominal)
lines(avgLT_ExpReturns_TIPS,col="blue")
lines(avgLT_ExpReturns_GDPLB,col="red")

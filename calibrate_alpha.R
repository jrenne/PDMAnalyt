# ==============================================================================
# Calibration of alpha
# ==============================================================================

candidate_values <- c(.05,.1,.2)

# Issuance of nominal bonds:
Model$kappa_pi <- 0
Model$kappa_y  <- 0

res_stat_distri_and_rbar <- compute_stat_distri_and_rbar(Model)
r_bar  <- res_stat_distri_and_rbar$r_bar
max_rr <- (r_bar + .1) * Model$d_star
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

# Average surplus:
Model$beta <- 0

# Debt threshold:
Model$d_star <- 1.1

for(alpha in candidate_values){
  
  Model$alpha <- alpha
  
  # Solve model:
  res0_nominal <- solve_ToyModel(all_d,all_rr,all_eps,proba_eps,
                                 Model,nb_iter = nb_iter)
  
  # Compute nominal yields:
  res_prices_nominal <- compute_bond_prices(Model, maxH,
                                            res0_nominal$indicators_x,
                                            res0_nominal$all_proba_def,
                                            res0_nominal$Probas)
  # Compute risk-free nominal yields:
  Model_RF <- Model
  Model_RF$RR <- 1
  res_prices_nominal_RF <- compute_bond_prices(Model_RF, maxH,
                                               res0_nominal$indicators_x,
                                               res0_nominal$all_proba_def,
                                               res0_nominal$Probas)
  
  spreads <- res_prices_nominal$all_rth - res_prices_nominal_RF$all_rth
  
  Ts <- which((res0_nominal$d==all_d[indic_DD])&
                (res0_nominal$d_1==all_d[indic_DD])&(res0_nominal$rr==all_rr[indic_rr]))
  
  avg_spreads <- c(t(stat_distri) %*% spreads[Ts,])
  
  print(avg_spreads[10])
}




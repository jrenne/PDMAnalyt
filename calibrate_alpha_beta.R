# ==============================================================================
# Calibration of beta
# ==============================================================================

Model$s_star <- .05
candidate_alpha_values  <- c(.1,.3,.5)
candidate_beta_values   <- c(.1,.2,.3)
candidate_d_star_values <- c(.4,.6,.8)

Model$RR <- .5
Model$mu_eta <- 1 * Model$mu_y

# Issuance of nominal bonds:
Model$kappa_pi <- 0
Model$kappa_y  <- 0

nb_iter <- 20 # used to solve model

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

# # Values of debt and debt service at which we compute spread:
# DD <- 1.0
# indic_DD <- which.min(abs(all_d - DD))
# rr <- .04
# indic_rr <- which.min(abs(all_rr - rr))

# Define targets:
Targets <- list(spread_in_bps = 20,
                mean_d_in_percent = 80)

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

# res0_nominal <- solve_ToyModel(all_d,all_rr,all_eps,proba_eps,
#                                Model,nb_iter = nb_iter)
# # Compute unconditional distribution of d:
# p <- compute_uncond_distri(res0_nominal$indicators_x,res0_nominal$Probas,1000)
# distri_d  <- compute_distri_x(all_d,res0_nominal$d,p)
# lines(all_d,distri_d,col="red")

# distri_rr  <- compute_distri_x(all_rr,res0_nominal$rr,p)
# plot(all_rr,distri_rr,type="l")

best_distance <- 1000000

for(alpha in candidate_alpha_values){
  for(beta in candidate_beta_values){
    for(d_star in candidate_d_star_values){
      
      Model$alpha  <- alpha
      Model$beta   <- beta
      Model$d_star <- d_star
      
      # Solve model:
      res0_nominal <- solve_ToyModel(all_d,all_rr,all_eps,proba_eps,
                                     Model,nb_iter = nb_iter)
      
      # Compute average of debt:
      p <- compute_uncond_distri(res0_nominal$indicators_x,res0_nominal$Probas,1000)
      distri_d  <- compute_distri_x(all_d,res0_nominal$d,p)
      mean_d <- sum(distri_d * all_d)
      
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
      
      avg_spreads <- c(t(p) %*% spreads)
      
      spread_in_bps    <- 10000*avg_spreads[10]
      mean_d_in_percent <- 100*mean_d
      
      distance <- (spread_in_bps - Targets$spread_in_bps)^2 +
        (mean_d_in_percent - Targets$mean_d_in_percent)^2
      
      print("------------------------------")
      print(paste("Value of alpha: ",alpha,"; value of beta: ",beta,
                  "; value of d_star: ",d_star,sep=""))
      print(paste("Average spread (in bps): ",round(spread_in_bps,1),sep=""))
      print(paste("Average debt (in percent): ",round(mean_d_in_percent,1),sep=""))
      print(paste("Distance to target: ",distance,sep=""))
      
      if(distance < best_distance){
        best <- list(alpha = alpha,
                     beta = beta,
                     d_star = d_star,
                     spread_in_bps = spread_in_bps,
                     mean_d_in_percent = mean_d_in_percent,
                     distance = distance,
                     p = p,
                     distri_d = distri_d)
        best_distance <- distance
      }
    }
  }
}

Model$alpha  <- best$alpha
Model$beta   <- best$beta
Model$d_star <- best$d_star

plot(all_d,best$distri_d,type="l")


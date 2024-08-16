# ==============================================================================
# Calibration of beta
# ==============================================================================

nb_grid <- 25 # number of values per state variable
Model$sigma_eps <- .03

grids <- make_grid(nb_grid = 23,
                   min_d = .4,
                   max_d = 1.5,
                   min_rr = 0,
                   max_rr = .15,
                   sigma_eps = .03,
                   all_quantiles_eps = c(-2,-1,1,2))

candidate_alpha_values  <- c(.1,.2)
candidate_beta_values   <- c(.02,.05,.1,.2)
candidate_d_star_values <- c(.9,1)

Model$RR <- .5
Model$mu_eta <- 0 * Model$mu_y

maxH <- 10

# Issuance of nominal bonds:
Model$kappa_pi <- 0
Model$kappa_y  <- 0

nb_iter <- 20 # used to solve model

# Determine Chi by fitting average debt maturity:
res_LTnominal_prices <- compute_LTRF_bond_prices(Model,maxH)
stat_distri <- res_LTnominal_prices$stat_distri
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
                stdv_d_in_percent = 15)

d_bar <- .8 # used tom compute s_star based on deterministic steady state

best_distance <- 1000000

for(alpha in candidate_alpha_values){
  for(beta in candidate_beta_values){
    for(d_star in candidate_d_star_values){
      
      Model$alpha  <- alpha
      Model$beta   <- beta
      Model$d_star <- d_star
      
      # Determine s_star to have 
      res_aux <- compute_determ_steady_state(Model,
                                             indic_d_bar_from_s_star = 1,
                                             d_bar = d_bar)
      Model$s_star <- res_aux$s_star
      
      # Solve model:
      Model_solved <- solve_ToyModel(Model,grids,nb_iter = nb_iter)
      
      # Compute average of debt:
      p <- compute_uncond_distri(Model_solved$indicators_x,Model_solved$Probas,1000)
      distri_d <- compute_distri_x(grids$all_d,Model_solved$d,p)
      mean_d   <- sum(distri_d * grids$all_d)
      
      # Compute nominal yields:
      res_prices_nominal <- compute_bond_prices(Model_solved, maxH)
      # Compute risk-free nominal yields:
      Model_solved_RF          <- Model_solved
      Model_solved_RF$Model$RR <- 1
      res_prices_nominal_RF <- compute_bond_prices(Model_solved_RF, maxH)
      
      spreads <- res_prices_nominal$all_rth - res_prices_nominal_RF$all_rth
      
      avg_spreads <- c(t(p) %*% spreads)
      
      spread_in_bps     <- 10000 * avg_spreads[10]
      mean_d_in_percent <- 100 * mean_d
      var_d <- sum(distri_d * grids$all_d^2) - mean_d^2
      stdv_d_in_percent <- sqrt(var_d) * 100
      
      distance <- (spread_in_bps - Targets$spread_in_bps)^2 +
        (stdv_d_in_percent - Targets$stdv_d_in_percent)^2
      
      print("------------------------------")
      print(paste("Value of alpha: ",alpha,"; value of beta: ",beta,
                  "; value of d_star: ",d_star,sep=""))
      print(paste("Average spread (in bps): ",round(spread_in_bps,1),sep=""))
      print(paste("Average debt (in percent): ",round(mean_d_in_percent,1),sep=""))
      print(paste("Stdv debt (in percent): ",round(stdv_d_in_percent,1),sep=""))
      print(paste("Distance to target: ",distance,sep=""))
      
      if(distance < best_distance){
        best <- list(alpha = alpha,
                     beta = beta,
                     d_star = d_star,
                     s_star = res_aux$s_star,
                     spread_in_bps = spread_in_bps,
                     mean_d_in_percent = mean_d_in_percent,
                     stdv_d_in_percent = stdv_d_in_percent,
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
Model$s_star <- best$s_star

plot(grids$all_d,best$distri_d,type="l")

distri_rr  <- compute_distri_x(grids$all_rr,Model_solved$rr,best$p)
plot(grids$all_r,distri_rr,type="l")


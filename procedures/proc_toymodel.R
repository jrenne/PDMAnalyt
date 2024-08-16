
solve_ToyModel_notRcpp <- function(Model,grids,nb_iter){
  
  all_d     <- grids$all_d
  all_rr    <- grids$all_rr
  all_eps   <- grids$all_eps
  proba_eps <- grids$proba_eps
  
  nb_grid_d  <- dim(all_d)[1]
  nb_grid_rr <- dim(all_rr)[1]
  nb_m       <- dim(Model$Omega)[1]
  nb_states  <- nb_grid_d^2 * nb_grid_rr * nb_m
  nb_eps     <- dim(all_eps)[1]
  
  vec_ones_d   <- matrix(1,nb_grid_d)
  vec_ones_rr  <- matrix(1,nb_grid_rr)
  vec_ones_m   <- matrix(1,nb_m)
  vec_ones_eps <- matrix(1,nb_eps)
  
  kappa_pi <- Model$kappa_pi
  kappa_y  <- Model$kappa_y
  nu_pi    <- Model$nu_pi
  nu_y     <- Model$nu_y
  mu_pi    <- Model$mu_pi
  mu_y     <- Model$mu_y
  mu_eta   <- Model$mu_eta
  Gamma    <- Model$gamma
  delta    <- Model$delta
  Dy_bar   <- Model$Dy_bar
  Pi_bar   <- Model$Pi_bar
  beta     <- Model$beta
  chi      <- Model$chi
  Omega    <- Model$Omega
  d_star   <- Model$d_star
  alpha    <- Model$alpha
  s_star   <- Model$s_star
  RR       <- Model$RR
  sigma_eps <- Model$sigma_eps
  sigma_nu  <- Model$sigma_nu
  
  # Compute stationary macro distribution:
  Omega_2h <- t(Omega)
  for(i in 1:10){
    Omega_2h <- Omega_2h %*% Omega_2h
  }
  stat_distri <- Omega_2h[,1]
  
  # s.d.f. specification (associated with composite index)
  res_SDF <- compute_SDF(Model)
  mu_f0   <- res_SDF$mu_f0
  mu_f1   <- res_SDF$mu_f1
  
  nu <- - Gamma * nu_y + (kappa_pi - 1) * nu_pi + kappa_y * nu_y
  
  d     <- vec_ones_m  %x% vec_ones_rr %x% vec_ones_d %x% all_d
  d_1   <- vec_ones_m  %x% vec_ones_rr %x% all_d      %x% vec_ones_d
  rr    <- vec_ones_m  %x% all_rr      %x% vec_ones_d %x% vec_ones_d
  Pi    <- mu_pi       %x% vec_ones_rr %x% vec_ones_d %x% vec_ones_d
  Dy    <- mu_y        %x% vec_ones_rr %x% vec_ones_d %x% vec_ones_d
  s_m   <- mu_eta      %x% vec_ones_rr %x% vec_ones_d %x% vec_ones_d
  
  all_d_t     <- matrix(d,   nb_states,nb_eps*nb_m)
  all_d_t_1   <- matrix(d_1, nb_states,nb_eps*nb_m)
  all_rr_t    <- matrix(rr,  nb_states,nb_eps*nb_m)
  all_Pi_t    <- matrix(Pi,  nb_states,nb_eps*nb_m)
  all_Dy_t    <- matrix(Dy,  nb_states,nb_eps*nb_m)
  all_s_m_t   <- matrix(s_m, nb_states,nb_eps*nb_m)
  
  Dlast <- diag(c(exp(mu_f1)))
  Dbetw <- diag(c(exp(mu_f0 + mu_f1)))
  D1rst <- diag(c(exp(mu_f0)))
  Mlast <- Dlast %*% t(Omega)
  Mbetw <- Dbetw %*% t(Omega)
  M1rst <- D1rst
  Pstar = t(vec_ones_m) %*% Mlast %*% solve(diag(nb_m) - chi * Mbetw) %*% M1rst
  OnepChiPstar <- 1 + chi * Pstar
  avg_LT_price <- c(OnepChiPstar %*% stat_distri)
  rstar <- t(1/Pstar - 1 + chi) # will be used to initialize q
  
  #r_bar <- 1/avg_LT_price - 1 + chi
  r_bar <- ((1 - chi) * sum(OnepChiPstar * stat_distri) - 1) / (1 - sum(OnepChiPstar * stat_distri)) ;
  
  all_eps_tp1      <- t(matrix(vec_ones_m   %x% all_eps, nb_eps * nb_m,nb_states))
  all_Pi_tp1       <- t(matrix(mu_pi        %x% vec_ones_eps, nb_eps * nb_m,nb_states))
  all_Dy_tp1       <- t(matrix(mu_y         %x% vec_ones_eps, nb_eps * nb_m,nb_states))
  all_f_tp1        <- t(matrix(mu_f1        %x% vec_ones_eps, nb_eps * nb_m,nb_states)) +
    matrix(mu_f0 %x% vec_ones_rr %x% vec_ones_d %x% vec_ones_d,nb_states,nb_eps * nb_m)
  all_OnepChiPstar <- t(matrix(OnepChiPstar %x% vec_ones_eps, nb_eps * nb_m,nb_states))
  all_eta_tp1      <- all_eps_tp1 +
    t(matrix(mu_eta %x% vec_ones_eps, nb_eps * nb_m,nb_states)) - all_s_m_t
  
  # Probas <- (Omega %x% matrix(1/nb_eps,1,nb_eps)) %x%
  #   matrix(1,nb_grid_d^2*nb_grid_rr,1)
  Probas <- (Omega %x% t(proba_eps)) %x% matrix(1,nb_grid_d^2*nb_grid_rr,1)
  
  q   <- (rstar + .001) %x% matrix(1,nb_grid_d^2*nb_grid_rr,1)
  q0  <- (rstar + .00 ) %x% matrix(1,nb_grid_d^2*nb_grid_rr,1) # risk-free rate
  
  q_1 <- q # this is used to compute chges in the recursions
  
  all_zeta_t   <- exp((kappa_pi - 1) * all_Pi_t   + (kappa_y - 1) * all_Dy_t)
  all_zeta_tp1 <- exp((kappa_pi - 1) * all_Pi_tp1 + (kappa_y - 1) * all_Dy_tp1)
  
  indicators_d_t   <- matrix(1:nb_grid_d,nb_states,nb_eps*nb_m)
  indicators_m_tp1 <- t(matrix((1:nb_m) %x% vec_ones_eps,nb_m*nb_eps,nb_states))
  
  
  indicators_x <- NULL
  
  # all_lambdas   <- pmax(beta * (all_d_t - d_star) + all_eta_tp1 - s_star,0)
  all_lambdas   <- pmax(all_d_t - d_star,0)
  all_proba_def <- 1 - exp(- alpha * all_lambdas)
  
  if(nb_iter>0){
    for(i in 1:nb_iter){
      print(i)
      
      # update all_q matrix:
      all_q_t  <- matrix(q,nb_states,nb_eps*nb_m)
      
      all_rr_tp1 <- all_q_t * all_zeta_tp1 * (all_d_t - chi * all_zeta_t * all_d_t_1) + 
        chi * all_zeta_tp1 * all_rr_t
      # all_d_tp1  <- all_zeta_tp1 * all_d_t - beta * (all_d_t - d_star) -
      #   all_eta_tp1 + all_rr_tp1
      all_d_tp1  <- all_zeta_tp1 * all_d_t - s_star - beta * all_d_t -
        all_eta_tp1 + all_rr_tp1
      
      # Match the previous states to the closest ones in the grid
      indicators_d_tp1   <- apply(all_d_tp1,   c(1,2),function(x){which((x-all_d )^2==min((x-all_d )^2))[1]})
      indicators_rr_tp1  <- apply(all_rr_tp1,  c(1,2),function(x){which((x-all_rr)^2==min((x-all_rr)^2))[1]})
      
      indicators_x <- indicators_d_tp1 + (indicators_d_t - 1) * nb_grid_d + 
        (indicators_rr_tp1 - 1) * nb_grid_d^2 +
        (indicators_m_tp1 - 1) * nb_grid_d^2 * nb_grid_rr
      
      # Compute right-hand side of the equation, conditional on varepsilon:
      all_q_tp1 <- matrix(q[c(indicators_x)],nb_states,nb_eps*nb_m)
      E <- exp(all_f_tp1) * ((1 + all_q_tp1)/(1 + all_q_tp1 - chi) +
                               all_proba_def * (exp(nu)*RR*all_OnepChiPstar -
                                                  (1 + all_q_tp1)/(1 + all_q_tp1 - chi)))
      q <- chi - 1 + 1/apply(E * Probas,1,sum)
      
      # Computation of risk free rate (RR=1):
      all_q0_tp1 <- matrix(q0[c(indicators_x)],nb_states,nb_eps*nb_m)
      E0 <- exp(all_f_tp1) * ((1 + all_q0_tp1)/(1 + all_q0_tp1 - chi) +
                                all_proba_def * (exp(nu)*all_OnepChiPstar -
                                                   (1 + all_q0_tp1)/(1 + all_q0_tp1 - chi)))
      Q0 <- (chi - 1 + 1/E0) * Probas
      q0 <- apply(Q0,1,sum)
      
      # Compute change in q:
      chges_in_q <- q - q_1
      print(
        paste("Max percent chge in q (over states): ",max(abs(chges_in_q)),sep="")
      )
      q_1 <- q
    }
  }
  
  return(list(
    q = q,
    q0 = q0,
    d = d,
    d_1 = d_1,
    rr = rr,
    Pi = Pi,
    Dy = Dy,
    r_bar = r_bar,
    #q_chge = chges_in_q,
    indicators_x = indicators_x,
    all_proba_def = all_proba_def,
    Probas = Probas,
    mu_f0 = mu_f0,
    mu_f1 = mu_f1,
    nu = nu,
    stat_distri = stat_distri,
    M1rst = M1rst,
    Mbetw = Mbetw,
    Mlast = Mlast,
    Pstar = Pstar,
    rstar = rstar,
    Model = Model,
    all_d = all_d,
    all_rr = all_rr))
}

compute_LTRF_bond_prices_notRcpp <- function(Model,
                                             maxH){
  # This procedure computes bond prices (after default, i.e. D_t = 1).
  
  Omega <- Model$Omega
  
  nb_m <- dim(Model$Omega)[1]
  vec_1_m <- matrix(1,nb_m,1)
  
  # Get SDF specification:
  res_SDF <- compute_SDF(Model)
  mu_f0         <- res_SDF$mu_f0
  mu_f1         <- res_SDF$mu_f1
  mu_f1_nominal <- res_SDF$mu_f1_nominal +
    Model$kappa_pi * Model$mu_pi + Model$kappa_y * Model$mu_y
  
  Mlast <- diag(c(exp(mu_f1_nominal))) %*% t(Omega)
  M1rst <- diag(c(exp(mu_f0)))
  Mbetw <- diag(c(exp(mu_f1_nominal + mu_f0))) %*% t(Omega)
  
  # Compute stationary distribution:
  Omega_2h <- t(Model$Omega)
  for(i in 1:10){
    Omega_2h <- Omega_2h %*% Omega_2h
  }
  stat_distri <- Omega_2h[,1]
  
  # To compute expected returns:
  D_exp_indextOmega = diag(c(exp(Model$kappa_pi * Model$mu_pi +
                                   Model$kappa_y * Model$mu_y))) %*% t(Omega)
  D_exp_indextOmega_h <- diag(nb_m)
  
  all_LT_Bth <- NULL
  all_LT_rth <- NULL
  all_LT_ExpReturn_th <- NULL
  
  # Compute bond prices:
  Mbetw_h <- diag(nb_m)
  for(h in 1:maxH){
    Bth <- t(vec_1_m) %*% Mlast %*% Mbetw_h %*% M1rst
    Mbetw_h <- Mbetw_h %*% Mbetw
    
    rth <- - log(Bth) / h
    
    D_exp_indextOmega_h = D_exp_indextOmega_h %*% D_exp_indextOmega
    ExpReturn_th <- 1/Bth * (t(vec_1_m) %*% D_exp_indextOmega_h)
    
    all_LT_Bth          <- cbind(all_LT_Bth,c(Bth))
    all_LT_rth          <- cbind(all_LT_rth,c(rth))
    all_LT_ExpReturn_th <- cbind(all_LT_ExpReturn_th,c(ExpReturn_th))
  }
  
  return(list(all_LT_Bth = all_LT_Bth,
              all_LT_rth = all_LT_rth,
              all_LT_ExpReturn_th = all_LT_ExpReturn_th,
              stat_distri = stat_distri))
}


compute_proba_def_notRcpp <- function(Model_solved,maxH){
  
  indicators_x  <- Model_solved$indicators_x
  all_proba_def <- Model_solved$all_proba_def
  Probas        <- Model_solved$Probas
  
  # Compute probabilities of default:
  nb_states <- dim(indicators_x)[1]
  nb_eps    <- dim(indicators_x)[2]
  all.prob.def <- matrix(NaN,nb_states,maxH)
  p_h_1 <- matrix(0,nb_states,1)
  for(h in 1:maxH){
    all_p_h_1_tp1    <- matrix(p_h_1[c(indicators_x)],nb_states,nb_eps)
    p_h_aux          <-  all_p_h_1_tp1 + (1 - all_p_h_1_tp1)*all_proba_def
    p_h              <-  apply(p_h_aux*Probas,1,sum)
    all.prob.def[,h] <- p_h
    p_h_1 <- p_h
  }
  return(all.prob.def)
}


simul_RS <- function(Model,maxH,T){
  nb_regimes = dim(Model$Omega)[1]
  
  res_LTprices <- compute_LTRF_bond_prices(Model,maxH)
  stat_distri  <- res_LTprices$stat_distri
  state <- which(stat_distri==max(stat_distri))[1]
  states <- matrix(0,T,nb_regimes)
  Pi <- NULL
  Dy <- NULL
  for(t in 1:T){
    p    <- Model$Omega[state,]
    cump <- cumsum(p)
    u    <- runif(1)
    state <- which(cump > u)[1]
    
    # Store values:
    Pi <- c(Pi,Model$mu_pi[state])
    Dy <- c(Dy,Model$mu_y[state])
    states[t,state] <- 1
  }
  
  return(list(
    Pi = Pi,
    Dy = Dy,
    states = states
  ))
}



logist <- function(x){
  return(exp(x)/(1+exp(x)))
}

inv_logist <- function(x){
  return(log(x/(1-x)))
}

# compute_LT_yds <- function(Model){
#   Model_nominal <- Model
#   Model_nominal$kappa_pi <- 0
#   Model_nominal$kappa_y  <- 0
#   res_LTprices <- compute_LTRF_bond_prices(Model_nominal,maxH=10)
#   avgLT_nom_yields  <- c(t(res_LTprices$stat_distri) %*% res_LTprices$all_LT_rth)
#   Model_real <- Model
#   Model_real$kappa_pi <- 1
#   Model_real$kappa_y  <- 0
#   res_LTprices <- compute_LTRF_bond_prices(Model_real,maxH=10)
#   avgLT_rea_yields  <- c(t(res_LTprices$stat_distri) %*% res_LTprices$all_LT_rth)
#   return(list(yds_nom = avgLT_nom_yields,
#               yds_rea = avgLT_rea_yields,
#               stat_distri = res_LTprices$stat_distri))
# }


make_model <- function(param,Model_ini){
  Model0 <- Model_ini
  nb_m <- dim(Model$Omega)[1]
  Model0$mu_pi <- matrix(min_Pi + (max_Pi - min_Pi) * logist(param[1:nb_m]),ncol=1)
  Model0$mu_y  <- matrix(min_Dy + (max_Dy - min_Dy) * logist(param[(nb_m+1):(2*nb_m)]),ncol=1)
  Omega <- matrix(0,nb_m,nb_m)
  count <- 2*nb_m
  for(i in 1:nb_m){
    sum_omegas <- 0
    for(j in 1:(nb_m-1)){
      count <- count + 1
      Omega[i,j] <- logist(param[count]) * (1 - sum_omegas)
      sum_omegas <- sum_omegas + Omega[i,j]
    }
    Omega[i,nb_m] <- 1 - sum_omegas
  }
  Model0$Omega <- Omega
  count <- count + 1
  Model0$gamma <- min_gamma + (max_gamma - min_gamma) * logist(param[count])
  return(Model0)
}

model2param <- function(Model){
  param <- NULL
  param[1:nb_m] <- inv_logist((Model$mu_pi - min_Pi)/(max_Pi - min_Pi))
  param[(nb_m+1):(2*nb_m)] <- inv_logist((Model$mu_y - min_Dy)/(max_Dy - min_Dy))
  count <- 2*nb_m
  Omega <- Model$Omega
  for(i in 1:nb_m){
    sum_omegas <- 0
    for(j in 1:(nb_m-1)){
      count <- count + 1
      param[count] <- inv_logist(Omega[i,j]/(1 - sum_omegas))
      sum_omegas <- sum_omegas + Omega[i,j]
    }
  }
  count <- count + 1
  param[count] <- inv_logist((Model$gamma - min_gamma)/(max_gamma - min_gamma))
  return(param)
}

compute_distance <- function(param,targets,Model_ini){
  Model   <- make_model(param,Model_ini)
  
  # Nominal yields:
  Model_nominal <- Model
  Model_nominal$kappa_pi <- 0
  Model_nominal$kappa_y  <- 0
  res_LTnominal_prices <- compute_LTRF_bond_prices(Model_nominal,maxH=10)
  stat_distri <- res_LTnominal_prices$stat_distri
  avg_nom_yds <- c(t(stat_distri) %*% res_LTnominal_prices$all_LT_rth)
  
  # Real yields:
  Model_real <- Model
  Model_real$kappa_pi <- 1
  Model_real$kappa_y  <- 0
  res_LTreal_prices <- compute_LTRF_bond_prices(Model_real,maxH=10)
  avg_real_yds <- c(t(stat_distri) %*% res_LTreal_prices$all_LT_rth)
  
  # Average inflation and GDP growth:
  avg_Pi <- sum(stat_distri * Model$mu_pi)
  avg_Dy <- sum(stat_distri * Model$mu_Dy)
  
  Var_nom_yds  <- t(stat_distri) %*% res_LTnominal_prices$all_LT_rth^2 - avg_nom_yds^2
  Var_real_yds <- t(stat_distri) %*% res_LTreal_prices$all_LT_rth^2    - avg_real_yds^2
  
  Std_nom_yds  <- sqrt(Var_nom_yds)
  Std_real_yds <- sqrt(Var_real_yds)
  
  distance <- 5000000 * ((avg_nom_yds[10] - avg_nom_yds[1] - targets$target_slop_nom)^2 +
                           .25*(avg_nom_yds[10] - targets$target_10_nom)^2 +
                           (avg_real_yds[10] - avg_real_yds[2] - targets$target_slop_rea)^2 +
                           .25*(avg_real_yds[10] - targets$target_10_rea)^2 + 
                           .1*(avg_Pi - targets$target_avg_Pi)^2 +
                           .1*(avg_Dy - targets$target_avg_Dy)^2 +
                           .0 * (Std_nom_yds[10]  - targets$target_std_10_nom)^2 +
                           .2 * (Std_real_yds[10] - targets$target_std_10_rea)^2 +
                           1000 * (avg_nom_yds[10] - avg_real_yds[10] -
                                   avg_Pi - targets$IRP10)^2)
  
  # Add penalty when yield curves are not monotonously increasing:
  distance <- distance + 10000*(avg_nom_yds[2]<avg_nom_yds[1])*(avg_nom_yds[1]-avg_nom_yds[2])
  #distance <- distance + 10000*(avg_real_yds[2]<avg_real_yds[1])*(avg_real_yds[1]-avg_real_yds[2])
  
  return(distance)
}


make_StateSpace <- function(Model){
  # Compute (LT) bond prices:
  Model_nominal <- Model
  Model_nominal$kappa_pi <- 0
  Model_nominal$kappa_y  <- 0
  res_LTnominal_prices <- compute_LTRF_bond_prices(Model_nominal,maxH=10)
  
  nb_m <- dim(Model$Omega)[1]
  
  N <- rbind(rep(Model$std_Pi,nb_m),
             rep(Model$std_Dy,nb_m),
             rep(Model$std_nom_yd,nb_m),
             rep(Model$std_nom_yd,nb_m))
  M <- rbind(c(Model$mu_pi),
             c(Model$mu_y),
             res_LTnominal_prices$all_LT_rth[,1],
             res_LTnominal_prices$all_LT_rth[,10])
  F <- cbind(DATA$pi,
             DATA$dy,
             DATA$SVENY01/100,
             DATA$SVENY10/100)
  dates <- DATA$date[complete.cases(F)]
  F <- F[complete.cases(F),]
  
  return(list(M=M,N=N,F=F,dates=dates))
}

compute_loglik <- function(param,Model_ini){
  
  Model <- make_model(param,Model_ini)
  
  res_SS <- make_StateSpace(Model)
  M <- res_SS$M
  N <- res_SS$N
  F <- res_SS$F
  
  res_KH <- KH_filter(F,M,N,Model$Omega)
  # res_KH_noRcpp <- KH_filter_notRcpp(F,M,N,Model$Omega)
  # print(c(res_KH$loglik,res_KH_noRcpp$loglik))
  
  penalty <- 100 * sum((abs(param) > max_abs_param) * 
                         (abs(param) - max_abs_param)^2)
  
  return(-res_KH$loglik + penalty)
}


compute_total_distance <- function(param,targets,Model_ini){
  
  loss1 <- compute_distance(param,targets,Model_ini)
  loss2 <- compute_loglik(param,Model_ini)
  
  return(loss1 + loss2)
}


compute_multivariate_normal_notRcpp <- function(epsilon_matrix,
                                                Covariance){
  n <- dim(epsilon_matrix)[2]
  T <- dim(epsilon_matrix)[1]
  
  vec_Covariance_1 <- matrix(solve(Covariance),nrow=1)
  
  aux <- - (epsilon_matrix %x% matrix(1,1,n)) * (matrix(1,1,n) %x% epsilon_matrix) *
    (matrix(1,T,1) %*% vec_Covariance_1)/ 2
  
  loglik_vec <- - n/2 * log(2*pi) - 1/2*log(det(Covariance)) + apply(aux,1,sum)
  
  return(loglik_vec)
}

# Sigma <- diag(c(.2,2^2,1))
# Sigma[2,1] <- 2
# epsilon_matrix <- matrix(rnorm(12),4,3) %*% t(Sigma)
# Covariance <- Sigma %*% t(Sigma)
# compute_multivariate_normal_notRcpp(epsilon_matrix,Covariance)
# compute_multivariate_normal(epsilon_matrix,Covariance)


KH_filter_notRcpp <- function(F, M, N, Omega){
  
  J         <- ncol(M)  # number of regimes
  nb_dates  <- nrow(F)
  nb_obsvar <- ncol(F)
  
  ksi_matrix   <- matrix(0, nrow = nb_dates, ncol = J)
  ksi_1_matrix <- matrix(0, nrow = nb_dates, ncol = J)
  eta_matrix   <- matrix(0, nrow = nb_dates, ncol = J)
  ksi_t        <- matrix(0, nrow = J, ncol = 1)
  eta_t        <- matrix(0, nrow = J, ncol = 1)
  
  vec_1_dates <- matrix(1, nrow = nb_dates, ncol = 1)
  vec_1_J     <- matrix(1, nrow = J, ncol = 1)
  
  # Compute log-likelihood conditional on each regime:
  variances4regime <- matrix(0, nrow = nb_obsvar, ncol = 1)
  Covariance       <- matrix(0, nrow = nb_obsvar, ncol = nb_obsvar)
  for (j in 1:J){
    variances4regime <- (N[, j])^2
    Covariance <- diag(variances4regime)
    eta_matrix[, j] <- compute_multivariate_normal(F - vec_1_dates %*% t(M[, j]), Covariance)
  }
  
  # Compute stationary distribution:
  Omega_2h <- matrix(0, nrow = J, ncol = J)
  stat_distri <- matrix(0, nrow = J, ncol = 1)
  Omega_2h <- t(Omega)
  for (i in 1:10){
    Omega_2h <- Omega_2h %*% Omega_2h
  }
  stat_distri[] <- Omega_2h[, 1]
  ksi_t[] <- stat_distri
  
  for (t in 1:nb_dates){
    ksi_1_matrix[t, ] <- t(ksi_t)  # will be used to compute log-likelihood
    
    eta_t <- matrix(eta_matrix[t, ],ncol=1)  # eta_matrix is LOG-likelihood
    eta_t <- exp(eta_t)
    # Use update formula:
    ksi_t <- (t(Omega) %*% ksi_t) * eta_t
    normalisation_factor <- sum(ksi_t)
    ksi_t <- ksi_t / normalisation_factor
    
    ksi_matrix[t, ] <- t(ksi_t)
    eta_matrix[t, ] <- t(eta_t)
  }
  
  # Compute log-likelihood:
  loglik_mat <- (ksi_1_matrix %*% Omega) * eta_matrix
  loglik_vec <- log(loglik_mat %*% vec_1_J)
  loglik <- sum(loglik_vec)
  
  return(list(ksi_matrix = ksi_matrix,
              eta_matrix = eta_matrix,
              loglik_vec = loglik_vec,
              loglik = loglik))
}


run_strategy <- function(Model_solved,maxH,
                         nb_iter4probas = 1000){
  
  # Compute bond prices (no default):
  res_LTprices <- compute_LTRF_bond_prices(Model_solved$Model,maxH)
  
  # Compute average expected returns (no default):
  avgLT_ExpReturns  <- c(t(Model_solved$stat_distri) %*% res_LTprices$all_LT_ExpReturn_th)
  
  # Compute probabilities of default:
  PD <- compute_proba_def(Model_solved,maxH=maxH)
  
  # Compute bond prices:
  res_prices <- compute_bond_prices(Model_solved, maxH)
  
  # Compute risk-free yields:
  Model_solved_RF    <- Model_solved
  Model_solved_RF$Model$RR <- 1
  res_prices_RF <- compute_bond_prices(Model_solved_RF, maxH)
  
  # Compute unconditional distribution:
  p <- compute_uncond_distri(Model_solved$indicators_x,Model_solved$Probas,nb_iter4probas)
  
  # Compute cost and risk measures: --------------------------------------------
  
  # Debt-to-GDP:
  distri_d <- compute_distri_x(Model_solved$all_d,Model_solved$d,p)
  mean_d   <- sum(Model_solved$all_d * distri_d)
  stdv_d   <- sqrt(sum(Model_solved$all_d^2 * distri_d) - mean_d^2)
  
  # Changes in debt-to-GDP:
  all_Delta_d <- as.numeric(levels(as.factor(Model_solved$d - Model_solved$d_1)))
  all_Delta_d <- matrix(all_Delta_d,ncol=1)
  distri_Delta_d <- compute_distri_x(all_Delta_d,
                                     Model_solved$d - Model_solved$d_1,p)
  mean_Delta_d <- sum(all_Delta_d * distri_Delta_d)
  stdv_Delta_d <- sqrt(sum(all_Delta_d^2 * distri_Delta_d) - mean_Delta_d^2)
  
  # Debt service:
  distri_rr <- compute_distri_x(Model_solved$all_rr,Model_solved$rr,p)
  mean_rr   <- sum(Model_solved$all_rr * distri_rr)
  stdv_rr   <- sqrt(sum(Model_solved$all_rr^2 * distri_rr) - mean_rr^2)
  
  # average PDs:
  avg_PD <- c(t(p) %*% PD)
  
  # average spreads:
  spreads <- res_prices$all_rth - res_prices_RF$all_rth
  avg_spreads   <- c(t(p) %*% spreads)
  
  # Debt-at-Risk:
  DaR95 <- interpolate_quantile(Model_solved$all_d,distri_d,.95)
  DaR99 <- interpolate_quantile(Model_solved$all_d,distri_d,.99)
  
  return(list(p=p,
              DaR95 = DaR95,
              DaR99 = DaR99,
              res_prices = res_prices,
              res_prices_RF = res_prices_RF,
              PD = 100*PD, avg_PD = 100*avg_PD,
              spreads = 10000*spreads, avg_spreads = 10000*avg_spreads,
              res_LTprices = res_LTprices, avgLT_ExpReturns = avgLT_ExpReturns,
              distri_d=distri_d,mean_d=100*mean_d,stdv_d=100*stdv_d,
              distri_Delta_d=distri_Delta_d,mean_Delta_d=100*mean_Delta_d,stdv_Delta_d=100*stdv_Delta_d,
              distri_rr=distri_rr,mean_rr=100*mean_rr,stdv_rr=100*stdv_rr))
}


prepare_returns_yds <- function(Model,maxH){
  
  nb_m <- dim(Model$Omega)[1]
  
  # Compute (LT) nominal bond prices:
  Model_nominal <- Model
  Model_nominal$kappa_pi <- 0
  Model_nominal$kappa_y  <- 0
  res_LTnominal_prices <- compute_LTRF_bond_prices(Model_nominal,maxH=maxH)
  
  # Compute (LT) real bond prices:
  Model_TIPS <- Model
  Model_TIPS$kappa_pi <- 1
  Model_TIPS$kappa_y  <- 0
  res_LTTIPS_prices <- compute_LTRF_bond_prices(Model_TIPS,maxH=maxH)
  
  # Compute (LT) TIPS bond prices:
  Model_GDPLB <- Model
  Model_GDPLB$kappa_pi <- 1
  Model_GDPLB$kappa_y  <- 1
  res_LTGDPLB_prices   <- compute_LTRF_bond_prices(Model_GDPLB,maxH=maxH)
  
  stat_distri <- compute_stat_distri(Model)
  
  nominal_yds <- res_LTnominal_prices$all_LT_rth
  TIPS_yds    <- res_LTTIPS_prices$all_LT_rth
  GDPLB_yds   <- res_LTGDPLB_prices$all_LT_rth
  
  avg_nominal_yds <- t(stat_distri) %*% nominal_yds
  avg_TIPS_yds    <- t(stat_distri) %*% TIPS_yds
  avg_GDPLB_yds   <- t(stat_distri) %*% GDPLB_yds
  
  avg_nominal_returns <- t(stat_distri) %*% res_LTnominal_prices$all_LT_ExpReturn_th
  avg_TIPS_returns    <- t(stat_distri) %*% res_LTTIPS_prices$all_LT_ExpReturn_th
  avg_GDPLB_returns   <- t(stat_distri) %*% res_LTGDPLB_prices$all_LT_ExpReturn_th
  
  Var_nominal_yds <- t(stat_distri) %*% nominal_yds^2 - avg_nominal_yds^2
  Var_TIPS_yds    <- t(stat_distri) %*% TIPS_yds^2    - avg_TIPS_yds^2
  Var_GDPLB_yds   <- t(stat_distri) %*% GDPLB_yds^2   - avg_GDPLB_yds^2
  Std_nominal_yds <- sqrt(Var_nominal_yds)
  Std_TIPS_yds    <- sqrt(Var_TIPS_yds)
  Std_GDPLB_yds   <- sqrt(Var_GDPLB_yds)
  
  exp_log_nominal_returns <- log(res_LTnominal_prices$all_LT_ExpReturn_th)
  exp_log_TIPS_returns    <- log(res_LTTIPS_prices$all_LT_ExpReturn_th)
  exp_log_GDPLB_returns   <- log(res_LTGDPLB_prices$all_LT_ExpReturn_th)
  
  exp_annual_nominal_returns <- exp_log_nominal_returns/t(matrix(1:maxH,maxH,nb_m))
  exp_annual_TIPS_returns    <- exp_log_TIPS_returns/t(matrix(1:maxH,maxH,nb_m))
  exp_annual_GDPLB_returns   <- exp_log_GDPLB_returns/t(matrix(1:maxH,maxH,nb_m))
  
  avg_annual_nominal_returns <- c(t(stat_distri) %*% exp_annual_nominal_returns)
  avg_annual_TIPS_returns    <- c(t(stat_distri) %*% exp_annual_TIPS_returns)
  avg_annual_GDPLB_returns   <- c(t(stat_distri) %*% exp_annual_GDPLB_returns)
  
  Var_annual_nominal_returns <- t(stat_distri) %*% exp_annual_nominal_returns^2 - avg_annual_nominal_returns^2
  Var_annual_TIPS_returns    <- t(stat_distri) %*% exp_annual_TIPS_returns^2    - avg_annual_TIPS_returns^2
  Var_annual_GDPLB_returns   <- t(stat_distri) %*% exp_annual_GDPLB_returns^2   - avg_annual_GDPLB_returns^2
  Std_annual_nominal_returns <- sqrt(Var_annual_nominal_returns)
  Std_annual_TIPS_returns    <- sqrt(Var_annual_TIPS_returns)
  Std_annual_GDPLB_returns   <- sqrt(Var_annual_GDPLB_returns)
  
  return(list(res_LTnominal_prices = res_LTnominal_prices,
              res_LTTIPS_prices = res_LTTIPS_prices,
              res_LTGDPLB_prices = res_LTGDPLB_prices,
              nominal_yds = nominal_yds,
              TIPS_yds = TIPS_yds,
              GDPLB_yds = GDPLB_yds,
              avg_nominal_yds = avg_nominal_yds,
              avg_TIPS_yds = avg_TIPS_yds,
              avg_GDPLB_yds = avg_GDPLB_yds,
              avg_nominal_returns = avg_nominal_returns,
              avg_TIPS_returns = avg_TIPS_returns,
              avg_GDPLB_returns = avg_GDPLB_returns,
              Std_nominal_yds = Std_nominal_yds,
              Std_TIPS_yds = Std_TIPS_yds,
              Std_GDPLB_yds = Std_GDPLB_yds,
              exp_annual_nominal_returns = exp_annual_nominal_returns,
              exp_annual_TIPS_returns = exp_annual_TIPS_returns,
              exp_annual_GDPLB_returns = exp_annual_GDPLB_returns,
              avg_annual_nominal_returns = avg_annual_nominal_returns,
              avg_annual_TIPS_returns = avg_annual_TIPS_returns,
              avg_annual_GDPLB_returns = avg_annual_GDPLB_returns,
              Std_annual_nominal_returns = Std_annual_nominal_returns,
              Std_annual_TIPS_returns = Std_annual_TIPS_returns,
              Std_annual_GDPLB_returns = Std_annual_GDPLB_returns
  ))
}



prepare_and_solve_3 <- function(Model,grids,nb_iter){
  
  # Compute average growth and inflation:
  stat_distri <- compute_stat_distri(Model)
  mean_pi <- c(t(stat_distri) %*% Model$mu_pi)
  mean_y  <- c(t(stat_distri) %*% Model$mu_y)
  
  # Issuance of nominal bonds: ---------------------------------------------------
  Model_nominal <- Model
  Model_nominal$kappa_pi <- 0
  Model_nominal$kappa_y  <- 0
  # Correct chi:
  Model_nominal$chi <- Model$chi / exp(Model_nominal$kappa_pi*mean_pi + Model_nominal$kappa_y*mean_y)
  print("--- Solving nominal-bond model ---")
  Model_solved_nominal <- solve_ToyModel(Model_nominal,grids,nb_iter = nb_iter)
  
  # Issuance of TIPS: ------------------------------------------------------------
  Model_TIPS <- Model
  Model_TIPS$kappa_pi <- 1
  Model_TIPS$kappa_y  <- 0
  # Correct chi:
  Model_TIPS$chi <- Model$chi / exp(Model_TIPS$kappa_pi*mean_pi + Model_TIPS$kappa_y*mean_y)
  print("--- Solving TIPS model ---")
  Model_solved_TIPS <- solve_ToyModel(Model_TIPS,grids,nb_iter = nb_iter)
  
  # Issuance of GDP-LBs: ---------------------------------------------------------
  Model_GDPLB <- Model
  Model_GDPLB$kappa_pi <- 1
  Model_GDPLB$kappa_y  <- 1
  # Correct chi:
  Model_GDPLB$chi <- Model$chi / exp(Model_GDPLB$kappa_pi*mean_pi + Model_GDPLB$kappa_y*mean_y)
  print("--- Solving GDPLB model ---")
  Model_solved_GDPLB <- solve_ToyModel(Model_GDPLB,grids,nb_iter = nb_iter)
  
  return(list(Model_solved_nominal = Model_solved_nominal,
              Model_solved_TIPS    = Model_solved_TIPS,
              Model_solved_GDPLB   = Model_solved_GDPLB))
}


make_grid <- function(nb_grid,min_d=0,max_d,min_rr=0,max_rr,
                      sigma_eps,all_quantiles_eps){
  
  # Determine debt's grid:
  all_d <- matrix(seq(min_d,max_d,length.out=nb_grid),ncol=1)
  
  # Determine debt service's grid:
  all_rr <- matrix(seq(min_rr,max_rr,length.out = nb_grid),ncol=1)
  
  nb_states <- length(all_d)^2*length(all_rr)
  
  # Determine considered values of epsilon:
  # (Approach: based on quantiles of N(0,1), deduce probas):
  aa <- c(-Inf,all_quantiles_eps)
  bb <- c(all_quantiles_eps,Inf)
  phi_a <- dnorm(aa)
  phi_b <- dnorm(bb)
  Phi_a <- pnorm(aa)
  Phi_b <- pnorm(bb)
  all_eps <- sigma_eps * (phi_a - phi_b)/(Phi_b - Phi_a)
  all_eps <- matrix(all_eps,ncol=1)
  proba_eps <- matrix(pnorm(bb) - pnorm(aa),ncol=1)
  
  return(list(all_d = all_d,
              all_rr = all_rr,
              all_eps = all_eps,
              proba_eps = proba_eps))
}


compute_determ_steady_state <- function(Model,
                                        indic_d_bar_from_s_star = 1,
                                        d_bar = .8){
  
  # Compute average growth and inflation:
  stat_distri <- compute_stat_distri(Model)
  mean_pi <- c(t(stat_distri) %*% Model$mu_pi)
  mean_y  <- c(t(stat_distri) %*% Model$mu_y)
  
  # Approximate q mean:
  res_LT_bond_prices <- compute_LTRF_bond_prices(Model,maxH = 20)
  # avg_yields <- t(stat_distri) %*% res_LT_bond_prices$all_LT_rth
  avg_Bondprices <- t(stat_distri) %*% res_LT_bond_prices$all_LT_Bth
  avg_yields <- -1/1:maxH * log(avg_Bondprices)
  duration <- 10 # initial (guess) value
  for(i in 1:10){
    q_bar <- avg_yields[duration]
    D <- 1/(1 + q_bar - Model$chi)
    duration <- round(D)
  }
  
  zeta_bar     <- c(t(stat_distri) %*% exp((Model$kappa_pi-1)*Model$mu_pi + (Model$kappa_y-1)*Model$mu_y))
  zeta_bar_bar <- c(t(stat_distri) %*% exp(- Model$mu_pi - Model$mu_y))
  
  s_star <- Model$s_star
  beta   <- Model$beta
  
  if(indic_d_bar_from_s_star == 1){
    d_bar <- - s_star/(1 + beta - zeta_bar * (1 + q_bar))
  }else{
    s_star <- - d_bar * (1 + beta - zeta_bar * (1 + q_bar))
  }
  r_bar <- ((zeta_bar - zeta_bar_bar) + q_bar * zeta_bar) * d_bar
  s_bar <- s_star + beta * d_bar
  
  return(list(d_bar = d_bar,
              r_bar = r_bar,
              q_bar = q_bar,
              s_bar = s_bar,
              s_star = s_star))
}


interpolate_quantile <- function(all_x,distri_x,proba){
  cumdf <- cumsum(distri_x)
  
  if(sum(cumdf == proba)>0){
    index <- which(cumdf==proba)
    quantil <- all_x[index]
  }
  
  below <- which(cumdf<proba)
  below <- below[length(below)]
  value_below <- all_x[below]
  distance2below <- proba - cumdf[below]
  
  above <- which(cumdf>proba)[1]
  value_above <- all_x[above]
  distance2above <- cumdf[above] - proba
  
  Dtot <- distance2below + distance2above
  quantil <- value_below * distance2above/Dtot +
    value_above * distance2below/Dtot
  
  return(quantil)
}

ShadowInt_PD_not <- function(fs,alpha,sigma){
  # values.fs are the values of the (log) fiscal space
  P <- pnorm(fs/sigma) + 
    exp(alpha*fs + alpha^2 * sigma^2/2) * (1 - pnorm(fs/sigma + sigma * alpha))
  P <- 1 - P
  return(P)
}




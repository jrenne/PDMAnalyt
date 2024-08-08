
solve_ToyModel_notRcpp <- function(all_d,all_rr,all_eps,proba_def,
                                   Model,nb_iter){
  
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
  OnepChiPstar <- 1 +
    chi * t(vec_ones_m) %*% Mlast %*% solve(diag(nb_m) - chi * Mbetw) %*% M1rst
  avg_LT_price <- c(OnepChiPstar %*% stat_distri)
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
  
  q   <- matrix(r_bar + .01, nb_states, 1)
  q0  <- matrix(r_bar, nb_states, 1) # risk-free rate
  
  q_1 <- q # this is used to compute chges in the recursions
  
  all_zeta_t   <- exp((kappa_pi - 1) * all_Pi_t   + (kappa_y - 1) * all_Dy_t)
  all_zeta_tp1 <- exp((kappa_pi - 1) * all_Pi_tp1 + (kappa_y - 1) * all_Dy_tp1)
  
  indicators_d_t   <- matrix(1:nb_grid_d,nb_states,nb_eps*nb_m)
  indicators_m_tp1 <- t(matrix((1:nb_m) %x% vec_ones_eps,nb_m*nb_eps,nb_states))
  
  indicators_x <- NULL
  all_proba_def = NULL
  
  if(nb_iter>0){
    for(i in 1:nb_iter){
      print(i)
      
      # update all_q matrix:
      all_q_t  <- matrix(q,nb_states,nb_eps*nb_m)
      
      all_rr_tp1 <- all_q_t * all_zeta_tp1 * (all_d_t - chi * all_zeta_t * all_d_t_1) + 
        chi * all_zeta_tp1 * all_rr_t
      all_d_tp1  <- all_zeta_tp1 * all_d_t - beta * (all_d_t - d_star) -
        all_eta_tp1 + all_rr_tp1
      
      # match the previous states to the closest ones in the grid
      indicators_d_tp1   <- apply(all_d_tp1,   c(1,2),function(x){which((x-all_d )^2==min((x-all_d )^2))[1]})
      indicators_rr_tp1  <- apply(all_rr_tp1,  c(1,2),function(x){which((x-all_rr)^2==min((x-all_rr)^2))[1]})
      
      indicators_x <- indicators_d_tp1 + (indicators_d_t - 1) * nb_grid_d + 
        (indicators_rr_tp1 - 1) * nb_grid_d^2 +
        (indicators_m_tp1 - 1) * nb_grid_d^2 * nb_grid_rr
      
      all_lambdas   <- pmax(beta * (all_d_t - d_star) + all_eta_tp1 - s_star,0)
      all_proba_def <- 1 - exp(- alpha * all_lambdas)
      
      # Compute right-hand side of the equation, conditional on varepsilon:
      all_q_tp1 <- matrix(q[c(indicators_x)],nb_states,nb_eps*nb_m)
      E <- exp(all_f_tp1) * ((1 + all_q_tp1)/(1 + all_q_tp1 - chi) +
                               all_proba_def * (exp(nu)*RR*all_OnepChiPstar -
                                                  (1 + all_q_tp1)/(1 + all_q_tp1 - chi)))
      Q <- (chi - 1 + 1/E) * Probas
      q <- apply(Q,1,sum) # update q
      
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
    OnepChiPstar = OnepChiPstar,
    stat_distri = stat_distri,
    M1rst = M1rst,
    Mbetw = Mbetw,
    Mlast = Mlast))
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


compute_proba_def_notRcpp <- function(maxH,
                                      indicators_x,
                                      all_proba_def,
                                      Probas){
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

compute_LT_yds <- function(Model){
  Model_nominal <- Model
  Model_nominal$kappa_pi <- 0
  Model_nominal$kappa_y  <- 0
  res_LTprices <- compute_LTRF_bond_prices(Model_nominal,maxH=10)
  avgLT_nom_yields  <- c(t(res_LTprices$stat_distri) %*% res_LTprices$all_LT_rth)
  Model_real <- Model
  Model_real$kappa_pi <- 1
  Model_real$kappa_y  <- 0
  res_LTprices <- compute_LTRF_bond_prices(Model_real,maxH=10)
  avgLT_rea_yields  <- c(t(res_LTprices$stat_distri) %*% res_LTprices$all_LT_rth)
  return(list(yds_nom = avgLT_nom_yields,
              yds_rea = avgLT_rea_yields,
              stat_distri = res_LTprices$stat_distri))
}


make_model <- function(param){
  Model0 <- Model
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

compute_distance <- function(param,targets){
  Model <- make_model(param)
  yds <- compute_LT_yds(Model)
  yds_nom <- yds$yds_nom
  yds_rea <- yds$yds_rea
  
  stat_distri <- yds$stat_distri
  avg_Pi <- sum(stat_distri * Model$mu_pi)
  avg_Dy <- sum(stat_distri * Model$mu_Dy)
  
  distance <- 10000000 * ((yds_nom[10] - yds_nom[1] - targets$target_slop_nom)^2 +
                            .25*(yds_nom[10] - targets$target_10_nom)^2 +
                            (yds_rea[10] - yds_rea[1] - targets$target_slop_rea)^2 +
                            .25*(yds_rea[10] - targets$target_10_rea)^2 + 
                            0*(avg_Pi - targets$target_avg_Pi)^2 +
                            0*(avg_Dy - targets$target_avg_Dy)^2)
  return(distance)
}


compute_loglik <- function(param){
  
  Model <- make_model(param)
  
  # Compute (LT) bond prices:
  Model_nominal <- Model
  Model_nominal$kappa_pi <- 0
  Model_nominal$kappa_y  <- 0
  res_LTnominal_prices <- compute_LTRF_bond_prices(Model_nominal,maxH=10)
  
  nb_m <- dim(Model$Omega)[1]
  
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
  # res_KH_noRcpp <- KH_filter_notRcpp(F,M,N,Model$Omega)
  # print(c(res_KH$loglik,res_KH_noRcpp$loglik))
  
  penalty <- 100 * sum((abs(param) > max_abs_param) * 
                         (abs(param) - max_abs_param)^2)
  
  return(-res_KH$loglik + penalty)
}


compute_total_distance <- function(param,targets){
  
  loss1 <- compute_distance(param,targets)
  loss2 <- compute_loglik(param)
  
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


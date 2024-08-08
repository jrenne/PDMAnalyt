
mod.sol.ext <- function(model){
  # If indic.compute.l = 1, then computation of a.l and b.l using compute.AB.ell
  
  delta      <- model$delta
  gamma      <- model$gamma
  Phi        <- model$Phi
  Sigma.w    <- model$Sigma.w
  mu_c       <- model$mu_c
  sigma_c    <- model$sigma_c
  mu_y       <- model$mu_y
  sigma_y    <- model$sigma_y
  mu_pi      <- model$mu_pi
  sigma_pi   <- model$sigma_pi
  sigma_s    <- model$sigma_s
  mu_star    <- model$mu_star
  sigma_star <- model$sigma_star
  beta       <- model$beta
  d_star     <- model$d_star
  max.h      <- model$max.h
  Chi        <- model$Chi
  
  n_w   <- ncol(Phi)
  n_x   <- n_w + 9
  
  model.sol <- model
  
  # Compute unconditional variance of w:
  I <- diag(n_w)
  I_Phi <- solve(I - model.sol$Phi)
  Omega.w <- matrix(solve(diag(n_w^2) - Phi %x% Phi) %*% 
                      c(Sigma.w %*% t(Sigma.w)),n_w,n_w)
  model.sol$Omega.w <- Omega.w
  
  # Get parameterizations of zero-coupon sovereign bond yields:
  AB.zc.yds <- compute_AB_nom_yields(model.sol)
  # store results:
  model.sol$a.nom <- AB.zc.yds$all.a.h
  model.sol$b.nom <- AB.zc.yds$all.b.h
  
  res_determine_h_star <- determine_h_star(model.sol)
  frac.h.star  <- res_determine_h_star$frac.h.star
  h.star.lower <- res_determine_h_star$h.star.lower
  h.star.upper <- res_determine_h_star$h.star.upper
  h_star       <- (1 - frac.h.star) * h.star.lower +
    frac.h.star * h.star.upper
  q_bar <- (1 - frac.h.star) * model.sol$b.nom[h.star.lower] +
    frac.h.star * model.sol$b.nom[h.star.upper]
  b_h_star <- q_bar
  a_h_star <- matrix(
    (1 - frac.h.star) * model.sol$a.nom[,h.star.lower] +
      frac.h.star * model.sol$a.nom[,h.star.upper],
    nrow = 1)

  ### Solve for d_bar and r_bar:
  g_bar   <- mu_y + mu_pi
  M <- matrix(0,2,2)
  M[1,1] <- g_bar + beta
  M[1,2] <- -1
  M[2,1] <- q_bar * (1 - Chi*(1 - g_bar))
  M[2,2] <- - (1 + g_bar - Chi)
  aux      <- matrix(0,2,1)
  aux[1,1] <- beta * d_star
  solut    <- solve(M) %*% aux
  d_bar    <- solut[1]
  r_bar    <- solut[2]
  s_bar    <- beta * (d_bar - d_star)
  
  model.sol$q_bar    <- q_bar
  model.sol$b_h_star <- b_h_star
  model.sol$a_h_star <- a_h_star
  model.sol$h_star   <- h_star
  model.sol$d_bar    <- d_bar
  model.sol$r_bar    <- r_bar
  model.sol$g_bar    <- g_bar
  model.sol$s_bar    <- s_bar
  
  ### Construct the parameters equation by equation for x_t
  
  A0     <- diag(n_x)
  A1     <- matrix(0,n_x,n_x)
  Sigma0 <- matrix(0,n_x,n_w)
  
  # w equation:
  A1[1:n_w,1:n_w]     <- Phi
  Sigma0[1:n_w,1:n_w] <- Sigma.w
  
  # Delta(c) equation:
  A0[n_w+1,1:n_w] <- -sigma_c

  # Delta(y) equation:
  A0[n_w+2,1:n_w] <- -sigma_y

  # pi equation:
  A0[n_w+3,1:n_w] <- -sigma_pi

  # s equation:
  A0[n_w+4,1:n_w] <- -sigma_s
  A1[n_w+4,n_w+8] <- beta

  # s* equation:
  A0[n_w+5,1:n_w] <- -sigma_star

  # r equation:
  A0[n_w+6,n_w+6] <- 1 + g_bar
  A0[n_w+6,n_w+2] <- r_bar
  A0[n_w+6,n_w+3] <- r_bar
  A1[n_w+6,n_w+2] <- Chi * q_bar * d_bar
  A1[n_w+6,n_w+3] <- Chi * q_bar * d_bar
  A1[n_w+6,n_w+6] <- Chi
  A1[n_w+6,n_w+7] <- d_bar * ( 1 - Chi * (1 - g_bar))
  A1[n_w+6,n_w+8] <- q_bar
  A1[n_w+6,n_w+9] <- - q_bar * Chi * (1 - g_bar)
  
  # q equation:
  A0[n_w+7,1:n_w] <- - a_h_star
  
  # d equation:
  A0[n_w+8,n_w+2] <- d_bar
  A0[n_w+8,n_w+3] <- d_bar
  A0[n_w+8,n_w+4] <- 1
  A0[n_w+8,n_w+6] <- - 1
  A0[n_w+8,n_w+8] <- 1
  A1[n_w+8,n_w+8] <- 1
  
  # d_1 equation:
  A0[n_w+9,n_w+9] <- 1
  A1[n_w+9,n_w+8] <- 1
  
  x_bar <- matrix(0,n_x,1)
  x_bar[n_w+1] <- mu_c
  x_bar[n_w+2] <- mu_y
  x_bar[n_w+3] <- mu_pi
  x_bar[n_w+4] <- s_bar
  x_bar[n_w+5] <- mu_star
  x_bar[n_w+6] <- r_bar
  x_bar[n_w+7] <- q_bar
  x_bar[n_w+8] <- d_bar
  x_bar[n_w+9] <- d_bar
  
  A0_1 <- solve(A0)
  Phi_x   <- A0_1 %*% A1
  Sigma_x <- A0_1 %*% Sigma0
  Mu_x    <- (diag(n_x)-Phi_x) %*% x_bar
  
  model.sol$Mu_x     <- Mu_x
  model.sol$Phi_x    <- Phi_x
  model.sol$Sigma_x  <- Sigma_x
  model.sol$E_x      <- x_bar
  SS_x <- Sigma_x %*% t(Sigma_x)
  model.sol$V_x      <- matrix(solve(diag(n_x^2) - Phi_x %x% Phi_x) %*%
                                 c(SS_x),n_x,n_x)

  return(model.sol)
}











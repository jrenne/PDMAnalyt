
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
RR       <- .5 # Recovery rate
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

# ==============================================================================
# This script XXX
# ==============================================================================


gamma     <- 5       # risk aversion
delta     <- .99     # preference for present
chi       <- .7      # coupon decay
beta      <- .3      # sensitivity of surplus to debt level
d_star    <- .6      # targeted debt level
sigma_eps <- .04     # std dev of surplus shocks
alpha     <- 1     # elasticity of default proba
#alpha    <- 0     # elasticity of default proba
s_star    <- .03     # max surplus

n_m <- 3
rho <- .9
Omega <- rho * diag(n_m)
Omega[2:n_m,1:(n_m-1)] <- Omega[2:n_m,1:(n_m-1)] + (1-Omega[1,1])/2 * diag(n_m-1)
Omega[1:(n_m-1),2:n_m] <- Omega[1:(n_m-1),2:n_m] + (1-Omega[1,1])/2 * diag(n_m-1)
Omega[1,] <- c(rep(0,n_m-1),1)
Omega[n_m,] <- c(rep(0,n_m-2),1,0)
Omega[2,] <- c(1-rho,rho,rep(0,n_m-2))
#Omega[n_m,n_m-1] <- 1 - rho



# Omega[1,1] <- .6
# Omega[1,3] <- .4
# 
# Omega[3,3] <- .6
# Omega[3,2] <- .4


mu_star <- .02

# Post-default dynamics:
Dy_bar <- mu_star
Pi_bar <- .03

mult_pi <- gamma
mult_y  <- 2
mu_pi <- matrix(Pi_bar,n_m,1)
mu_y  <- matrix(Dy_bar,n_m,1)
mu_pi[1]   <- (1 + mult_pi)*Pi_bar
#mu_pi[n_m] <- (1 - mult_pi)*Pi_bar
mu_y[1]    <- (1 - mult_y) *Dy_bar
mu_y[n_m]  <- (1 + 0*mult_y) *Dy_bar

mu_pi[1] <- -.03
mu_y[1]  <- -.05



# 0+ ++ +0 0 0- -- -0 
Omega <- matrix(0,7,7)
Omega[1,1] <- 0
Omega[1,2] <- 1
Omega[2,2] <- .8
Omega[2,3] <- .2
Omega[3,3] <- 0
Omega[3,4] <- 1
Omega[4,1] <- .1
Omega[4,4] <- .8
Omega[4,5] <- .1
Omega[5,5] <- 0
Omega[5,6] <- 1
Omega[6,6] <- .8
Omega[6,7] <- .2
Omega[7,4] <- 1
Omega[7,7] <- 0
library(expm)
(t(Omega)) %^% 100

mu_y  <- matrix(c(2*mu_star,mu_star,-mu_star,mu_star,-mu_star,mu_star,2*mu_star))
# Supply:
mu_pi <- matrix(c(0, Pi_bar, 2*Pi_bar, Pi_bar, 2*Pi_bar, Pi_bar, 0))
# Demand:
mu_pi <- matrix(c(2*Pi_bar, Pi_bar, 0, Pi_bar, 0, Pi_bar, 2*Pi_bar))
# Deterministic inflation:
#mu_pi <- matrix(rep(Pi_bar,dim(Omega)[1]),ncol=1)


# # 0+ ++ +0 0 0- -- -0 
# Omega <- matrix(0,3,3)
# Omega[1,1] <- .5
# Omega[1,2] <- .5
# Omega[2,1] <- .2
# Omega[2,2] <- .6
# Omega[2,3] <- .2
# Omega[3,2] <- .5
# Omega[3,3] <- .5
# 
# mu_y  <- matrix(c(2*mu_star,mu_star,-2*mu_star))
# # Supply:
# mu_pi <- matrix(c(-2*Pi_bar, Pi_bar, 2*Pi_bar))
# # Demand:
# #mu_pi <- matrix(c(2*Pi_bar, Pi_bar, -2*Pi_bar))
# # Deterministic inflation:
# #mu_pi <- matrix(rep(Pi_bar,dim(Omega)[1]),ncol=1)


mu_eta <- .5 * mu_y

nu_y  <- -.2 * 0
nu_pi <- .05 * 0

kappa_pi <- 0
kappa_y  <- 0

# RR    <- 0*exp(-b*gamma)  # recovery rate
# RR    <- exp(-b*gamma)    # recovery rate
RR    <- .2

Model <- list(gamma=gamma,
              delta = delta,
              chi = chi,
              beta = beta,
              d_star = d_star,
              sigma_eps = sigma_eps,
              alpha = alpha,
              s_star = s_star,
              RR = RR,
              Omega = Omega,
              mu_pi = mu_pi,
              mu_y = mu_y,
              mu_eta = mu_eta,
              nu_pi = nu_pi,
              nu_y = nu_y,
              kappa_pi = kappa_pi,
              kappa_y = kappa_y)

nb_grid <- 17 # number of values per state variable

# Determine debt's grid:
center_d <- d_star
disp_d   <- d_star/3
quantiles <- seq(.1/nb_grid,1-.1/nb_grid,length.out = nb_grid)
all_d <- d_star + disp_d * qnorm(quantiles)
all_d <- matrix(all_d,ncol=1)
par(mfrow=c(1,2))
plot(all_d)

# Determine debt service's grid:
r      <- - log(delta) + mu_star
max_rr <- .15
all_rr <- seq(r - .03,max_rr,length.out = nb_grid)
all_rr <- matrix(all_rr,ncol=1)
plot(all_rr)

nb_states <- length(all_d)^2*length(all_rr)

# Determine considered values of epsilon:
# # ==> Note: they will be treated as equi-probable <==
# nb_eps <- 5
# smallest_proba <- .01
# proba_bin_eps <- 1/nb_eps
# proba_eps <- seq(proba_bin_eps,1-proba_bin_eps,by=proba_bin_eps)
# all_quantiles_eps <- qnorm(proba_eps)
# # Compute expectations of eps given that eps is in a given bin:
# # (using truncated normal distribution)
# aa <- c(-Inf,all_quantiles_eps)
# bb <- c(all_quantiles_eps,Inf)
# phi_a <- dnorm(aa)
# phi_b <- dnorm(bb)
# Phi_a <- pnorm(aa)
# Phi_b <- pnorm(bb)
# all_eps <- sigma_eps * (phi_a - phi_b)/(Phi_b - Phi_a)
# all_eps <- matrix(all_eps,ncol=1)
# proba_eps <- matrix(1/nb_eps,nb_eps,1)

# Inverse option to get all_eps and proba_eps: define quantiles and deduce probas:
all_quantiles_eps <- c(-2,-1,1,2)
aa <- c(-Inf,all_quantiles_eps)
bb <- c(all_quantiles_eps,Inf)
phi_a <- dnorm(aa)
phi_b <- dnorm(bb)
Phi_a <- pnorm(aa)
Phi_b <- pnorm(bb)
all_eps <- sigma_eps * (phi_a - phi_b)/(Phi_b - Phi_a)
all_eps <- matrix(all_eps,ncol=1)
proba_eps <- matrix(pnorm(bb) - pnorm(aa),ncol=1)


# # ==============================================================================
# # Prepare table with parametrization
# # ==============================================================================
# 
# Top.of.Table <- rbind(
#   "\\begin{table}[!h]",
#   "\\caption{Parameterization of stylized model}",
#   "\\label{tab:stylized}",
#   "\\begin{footnotesize}",
#   "\\begin{tabular*}{\\textwidth}{l@{\\extracolsep{\\fill}}cc}",
#   "")
# 
# latex.table <- NULL
# 
# latex.table <- rbind(Top.of.Table,
#                      "\\hline",
#                      "\\hline",
#                      "{\\bf Parameter}&{\\bf Notation}&{\\bf Value}\\\\",
#                      "\\hline",
#                      paste("Risk aversion&$\\gamma$&",
#                            round.fixed.length(gamma,0),"\\\\",sep=""),
#                      paste("GDP growth&$\\mu_y$&",
#                            round.fixed.length(mu,3),"\\\\",sep=""),
#                      paste("Preference for present&$\\delta$&",
#                            round.fixed.length(delta,3),"\\\\",sep=""),
#                      paste("GDP fall upon default&$b_{c}$&",
#                            round.fixed.length(b,2),"\\\\",sep=""),
#                      paste("Coupon decay rate&$\\chi$&",
#                            round.fixed.length(chi,2),"\\\\",sep=""),
#                      paste("Sensitivity of surplus to debt&$\\beta$&",
#                            round.fixed.length(beta,2),"\\\\",sep=""),
#                      paste("Targetted debt level&$d^*$&",
#                            round.fixed.length(d_star,2),"\\\\",sep=""),
#                      paste("Standard deviation of surplus shocks&$\\sigma_\\varepsilon$&",
#                            round.fixed.length(sigma_eps,3),"\\\\",sep=""),
#                      paste("Semi-elasticity of default proba. to fiscal space&$\\alpha$&",
#                            round.fixed.length(alpha,3),"\\\\",sep=""),
#                      paste("Maximum surplus&$s^*$&",
#                            round.fixed.length(s_star,3),"\\\\",sep="")
# )
# 
# latex.table <- rbind(latex.table,
#                      "\\hline",
#                      "\\hline",
#                      "\\end{tabular*}",
#                      "\\end{footnotesize}",
#                      "\\begin{footnotesize}",
#                      "\\parbox{\\linewidth}{",
#                      "Note: This table presents the parameterization of the stylized model.",
#                      "}",
#                      "\\end{footnotesize}",
#                      "\\end{table}")
# 
# latex.file <- "tables/table_toymodel_param.txt"
# write(latex.table, latex.file)
# 
# 
# 
# 
# 

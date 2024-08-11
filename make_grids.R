
nb_grid <- 27 # number of values per state variable

# Determine debt's grid:
center_d <- d_star
disp_d   <- d_star/3
quantiles <- seq(.1/nb_grid,1-.1/nb_grid,length.out = nb_grid)
all_d <- d_star + disp_d * qnorm(quantiles)
all_d <- matrix(all_d,ncol=1)

all_d <- matrix(seq(0,1.5,length.out=nb_grid),ncol=1)

par(mfrow=c(1,2))
plot(all_d)


# Determine debt service's grid:
r      <- d_star * .06
max_rr <- .15
all_rr <- seq(0,max_rr,length.out = nb_grid)
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


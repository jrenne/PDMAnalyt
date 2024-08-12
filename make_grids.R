
# Determine debt's grid:

# center_d <- Model$d_star
# disp_d   <- Model$d_star/3
# quantiles <- seq(.1/nb_grid,1-.1/nb_grid,length.out = nb_grid)
# all_d <- d_star + disp_d * qnorm(quantiles)
# all_d <- matrix(all_d,ncol=1)

all_d <- matrix(seq(0,1.5,length.out=nb_grid),ncol=1)

par(mfrow=c(1,2))
plot(all_d)

# Determine debt service's grid:
max_rr <- .12
all_rr <- seq(0,max_rr,length.out = nb_grid)
all_rr <- matrix(all_rr,ncol=1)
plot(all_rr)

nb_states <- length(all_d)^2*length(all_rr)

# Determine considered values of epsilon (define quantiles and deduce probas):
all_quantiles_eps <- c(-2,-1,1,2)
aa <- c(-Inf,all_quantiles_eps)
bb <- c(all_quantiles_eps,Inf)
phi_a <- dnorm(aa)
phi_b <- dnorm(bb)
Phi_a <- pnorm(aa)
Phi_b <- pnorm(bb)
all_eps <- Model$sigma_eps * (phi_a - phi_b)/(Phi_b - Phi_a)
all_eps <- matrix(all_eps,ncol=1)
proba_eps <- matrix(pnorm(bb) - pnorm(aa),ncol=1)


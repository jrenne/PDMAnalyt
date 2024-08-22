
# Figure illustrating credit-event pricing

Model_nominal <- Model
Model_nominal$kappa_pi <- 0
Model_nominal$kappa_y  <- 0


value_4_nu_y <- c(0,-.05)

nb_grid <- 17
nb_iter <- 30
nb_iter_sdf <- 10

grids <- make_grid(nb_grid = nb_grid,
                   min_d = .5,
                   max_d = 1.6,
                   min_rr=.0,
                   max_rr=.15,
                   sigma_eps = Model$sigma_eps,
                   all_quantiles_eps = c(-2,-1,1,2))


# Model_solved <- solve_ToyModel(Model,
#                                grids,
#                                nb_iter = nb_iter,
#                                nb_iter_sdf = nb_iter_sdf)
# a <- compute_SDF(Model)
# b <- compute_SDF_D(Model,a$mu_u,
#                    Model_solved$indicators_x,
#                    Model_solved$all_proba_def,
#                    Model_solved$Probas,
#                    nb_iter_sdf = 10)
# n_states <- dim(Model_solved$indicators_x)
# n_m <- dim(Model$Omega)
# y <- b$Mu_f1_real + (Model$kappa_pi - 1) * (Model$mu_pi %x% matrix(1,n_states/n_m,1)) +
#   (Model$kappa_y) * (Model$mu_y %x% matrix(1,n_states/n_m,1))
# y <- b$Mu_nu_real + (Model$kappa_pi - 1) * Model$nu_pi + 
#   (Model$kappa_y) * Model$nu_y


par(mfrow=c(1,2))

for(nu_y in value_4_nu_y){
  Model_nominal$nu_y <- nu_y
  
  Model_solved_nominal <- solve_ToyModel(Model_nominal,
                                         grids,
                                         nb_iter = nb_iter,
                                         nb_iter_sdf = nb_iter_sdf)
  # Compute uncond. distri:
  p <- compute_uncond_distri(Model_solved_nominal$indicators_x,
                             Model_solved_nominal$Probas,
                             nb_iter = 1000)
  
  # Compute nominal yields:
  res_prices <- compute_bond_prices(Model_solved_nominal,maxH = 10,
                                    nb_iter_sdf = nb_iter_sdf)
  Model_solved_nominal_RF <- Model_solved_nominal
  Model_solved_nominal_RF$Model$RR <- 1
  res_prices_RF <- compute_bond_prices(Model_solved_nominal_RF,maxH = 10,
                                       nb_iter_sdf = nb_iter_sdf)
  avg_yds    <- c(t(p) %*% res_prices$all_rth)
  avg_yds_RF <- c(t(p) %*% res_prices_RF$all_rth)
  
  plot(avg_yds,type="l",ylim=c(min(avg_yds,avg_yds_RF),
                               max(avg_yds,avg_yds_RF)))
  lines(avg_yds_RF,col="red")
}

#plot(Model_solved_nominal$q)


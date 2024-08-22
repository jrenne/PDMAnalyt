

grids <- make_grid(nb_grid = 27,
                   min_d = 0.2,
                   max_d = 1.6,
                   min_rr = 0,
                   max_rr=.15,
                   sigma_eps = Model$sigma_eps,
                   all_quantiles_eps = c(-2,-1,1,2))

RES <- prepare_returns_yds(Model_Demand,maxH)

plot(RES$avg_annual_nominal_returns,type="l")

Model_Demand$kappa_pi <- 0
Model_Demand$kappa_y  <- 0

Model_Demand$alpha <- 0
#Model_Demand$nu_y <- -.001

Model_Demand_solved <- solve_ToyModel(Model_Demand,
                                      grids,nb_iter,
                                      nb_iter_sdf)
p <- compute_uncond_distri(Model_Demand_solved$indicators_x,
                           Model_Demand_solved$Probas,nb_iter = 2000)

res_bond_prices <- compute_bond_prices(Model_Demand_solved,maxH=10,nb_iter_sdf)

avg_yds <- t(p) %*% res_bond_prices$all_rth
lines(c(avg_yds),col="red")



# Figure illustrating credit-event pricing

Model_nominal <- Model
Model_nominal$kappa_pi <- 0
Model_nominal$kappa_y  <- 0

nb_grid <- 27
source("make_grids.R")

nb_iter <- 30

Model_solved_nominal <- solve_ToyModel(all_d,all_rr,all_eps,proba_eps,
                                       Model_nominal,nb_iter = nb_iter)

# Compute nominal yields:
a <- compute_bond_prices(Model_solved_nominal,maxH = 10)



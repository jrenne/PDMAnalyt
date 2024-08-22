

# Test new compute_sdf (with default events)

Model$nu_y <- 0

Model$kappa_pi <- 0
Model$kappa_y  <- 0

nb_grid <- 20
nb_iter <- 20
nb_iter_sdf <- 10

grids <- make_grid(nb_grid = nb_grid,
                   min_d = .5,
                   max_d = 1.6,
                   min_rr=.0,
                   max_rr=.15,
                   sigma_eps = Model$sigma_eps,
                   all_quantiles_eps = c(-2,-1,1,2))

Model_solved <- solve_ToyModel(Model,grids,
                               nb_iter = nb_iter,
                               nb_iter_sdf = nb_iter_sdf)

Model_solved_notRcpp <- solve_ToyModel_notRcpp(Model,grids,
                                               nb_iter = nb_iter,
                                               nb_iter_sdf = nb_iter_sdf)
Model_solved_notRcpp$q0[1:5]
Model_solved$q0[1:5]

max(abs(Model_solved_notRcpp$all_proba_def - Model_solved$all_proba_def))

res_sdf_noEvent <- compute_SDF(Model)

res_sdf <- compute_SDF_D(Model,
                         mu_u_bar = res_sdf_noEvent$mu_u,
                         Model_solved$indicators_x,
                         Model_solved$all_proba_def,
                         Model_solved$Probas,
                         nb_iter_sdf = 0)

cbind(as.numeric(levels(as.factor(res_sdf$Mu_u))),sort(res_sdf_noEvent$mu_u))
cbind(as.numeric(levels(as.factor(res_sdf$Mu_f0))),sort(res_sdf_noEvent$mu_f0))
cbind(as.numeric(levels(as.factor(res_sdf$Mu_f1_nominal))),
      sort(res_sdf_noEvent$mu_f1_nominal))
cbind(as.numeric(levels(as.factor(res_sdf$Mu_f1))),sort(res_sdf_noEvent$mu_f1))

res_prices <- compute_bond_prices(Model_solved,maxH=10,nb_iter_sdf=10)





# Test new compute_sdf (with default events)

Model$kappa_pi <- 0
Model$kappa_y  <- 0

grids <- make_grid(nb_grid = 21,
                   min_d = .5,
                   max_d = 1.6,
                   min_rr=.0,
                   max_rr=.15,
                   sigma_eps = Model$sigma_eps,
                   all_quantiles_eps = c(-2,-1,1,2))

Model_solved <- solve_ToyModel(Model,grids,nb_iter = 30)

res_sdf_noEvent <- compute_SDF(Model)

res_sdf <- compute_SDF_D(Model,
              mu_u_bar = res_sdf_noEvent$mu_u,
              Model_solved$indicators_x,
              Model_solved$all_proba_def,
              Model_solved$Probas,
              nb_iter_sdf = 0)

cbind(as.numeric(levels(as.factor(res_sdf$mu_u))),sort(res_sdf_noEvent$mu_u))
cbind(as.numeric(levels(as.factor(res_sdf$mu_f0))),sort(res_sdf_noEvent$mu_f0))
cbind(as.numeric(levels(as.factor(res_sdf$mu_f1_nominal))),
      sort(res_sdf_noEvent$mu_f1_nominal))
cbind(as.numeric(levels(as.factor(res_sdf$mu_f1))),sort(res_sdf_noEvent$mu_f1))


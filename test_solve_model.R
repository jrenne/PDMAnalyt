
#. Test solve model

grids <- make_grid(nb_grid = 4,
                   min_d = .4,
                   max_d = 1.5,
                   min_rr = 0,
                   max_rr = .15,
                   sigma_eps = .03,
                   all_quantiles_eps = c(-2,-1,1,2))


Model_solved_notRcpp <- solve_ToyModel_notRcpp(Model,grids,nb_iter = 2)
Model_solved <- solve_ToyModel(Model,grids,nb_iter = 2)






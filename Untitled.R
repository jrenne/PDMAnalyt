
res0_nominal <- solve_ToyModel(all_d,all_rr,all_eps,proba_eps,
                               Model,nb_iter = 10)
res0_nominal_notRcpp <- solve_ToyModel_notRcpp(all_d,all_rr,all_eps,proba_eps,
                                               Model,nb_iter = 10)
# res0_nominal$q[c(1,400,1000,8000)]
# res0_nominal_notRcpp$q[c(1,400,1000,8000)]

maxH <- 10
res_LTprices_nominal <- compute_LTRF_bond_prices(Model,maxH)
avgLT_ExpReturns_nominal  <- c(t(res0_nominal$stat_distri) %*% res_LTprices_nominal$all_LT_ExpReturn_th)
PD_nominal <- compute_proba_def(maxH=maxH,
                                indicators_x  = res0_nominal$indicators_x,
                                all_proba_def = res0_nominal$all_proba_def,
                                Probas        = res0_nominal$Probas)

H <- 200
res_prices_nominal <- compute_bond_prices(Model, maxH = H,
                                          res0_nominal$indicators_x,
                                          res0_nominal$all_proba_def,
                                          res0_nominal$Probas)

res_prices_nominal$all_Bth %*% matrix(Model$chi^(0:(H-1)),ncol=1)
res0_nominal$Pstar


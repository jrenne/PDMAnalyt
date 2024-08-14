
nb_grid <- 27 # number of values per state variable
nb_iter <- 30 # to solve model

grids <- make_grid(nb_grid = 27,
                   max_d = 1.5,
                   max_rr=.12,
                   sigma_eps = .03,
                   all_quantiles_eps = c(-2,-1,1,2))

nb_iter <- 30 # for model solution

chi <- Model$chi
#chi <- .5

Models <- prepare_and_solve_3(Model,grids,nb_iter)

strat_nominal <- run_strategy(Models$Model_solved_nominal,maxH=10)
strat_TIPS    <- run_strategy(Models$Model_solved_TIPS,maxH=10)
strat_GDPLB   <- run_strategy(Models$Model_solved_GDPLB,maxH=10)

M <- rbind(c(strat_nominal$mean_d,strat_TIPS$mean_d,strat_GDPLB$mean_d),
           c(strat_nominal$stdv_d,strat_TIPS$stdv_d,strat_GDPLB$stdv_d),
           c(strat_nominal$mean_rr,strat_TIPS$mean_rr,strat_GDPLB$mean_rr),
           c(strat_nominal$stdv_rr,strat_TIPS$stdv_rr,strat_GDPLB$stdv_rr),
           c(strat_nominal$stdv_Delta_d,strat_TIPS$stdv_Delta_d,strat_GDPLB$stdv_Delta_d),
           c(strat_nominal$avg_PD[maxH],strat_TIPS$avg_PD[maxH],strat_GDPLB$avg_PD[maxH]),
           c(strat_nominal$avg_spreads[maxH],strat_TIPS$avg_spreads[maxH],strat_GDPLB$avg_spreads[maxH]))

colnames(M) <- c("Nominal","TIPS","GDPLB")
rownames(M) <- c("Mean d","Stdv d",
                 "Mean r","Stdv r",
                 "Stdv D(d)","avg PD","avg spds")

print(M)






nb_grid <- 27 # number of values per state variable
source("make_grids.R")

nb_iter <- 30 # for model solution

chi <- Model$chi
#chi <- .5

# Determine grid of debt service:
max_rr <- .10
all_rr <- seq(0,max_rr,length.out = nb_grid)
all_rr <- matrix(all_rr,ncol=1)

Models <- prepare_and_solve_3(Model,
                              all_d,all_rr,all_eps,proba_eps)

strat_nominal <- run_strategy(Models$Model_solved_nominal,maxH=10)
strat_TIPS    <- run_strategy(Models$Model_solved_TIPS,maxH=10)
strat_GDPLB   <- run_strategy(Models$Model_solved_GDPLB,maxH=10)

M <- rbind(c(strat_nominal$mean_d,strat_TIPS$mean_d,strat_GDPLB$mean_d),
           c(strat_nominal$stdv_d,strat_TIPS$stdv_d,strat_GDPLB$stdv_d),
           c(strat_nominal$mean_rr,strat_TIPS$mean_rr,strat_GDPLB$mean_rr),
           c(strat_nominal$stdv_rr,strat_TIPS$stdv_rr,strat_GDPLB$stdv_rr),
           c(strat_nominal$mean_Delta_d,strat_TIPS$mean_Delta_d,strat_GDPLB$mean_Delta_d),
           c(strat_nominal$stdv_Delta_d,strat_TIPS$stdv_Delta_d,strat_GDPLB$stdv_Delta_d))

colnames(M) <- c("Nominal","TIPS","GDPLB")
rownames(M) <- c("Mean d","Stdv d",
                 "Mean r","Stdv r",
                 "Mean D(d)","Stdv D(d)")

print(M)





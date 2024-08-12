
nb_grid <- 27 # number of values per state variable
source("make_grids.R")

nb_iter <- 30 # for model solution

chi <- Model$chi
chi <- .5

# Compute average growth and inflation:
res_stat_distri_and_rbar <- compute_stat_distri_and_rbar(Model)
stat_distri <- res_stat_distri_and_rbar$stat_distri
mean_pi <- c(t(stat_distri) %*% Model$mu_pi)
mean_y  <- c(t(stat_distri) %*% Model$mu_y)

# Determine grid of debt service:
max_rr <- .10
all_rr <- seq(0,max_rr,length.out = nb_grid)
all_rr <- matrix(all_rr,ncol=1)

# Issuance of nominal bonds: ---------------------------------------------------
Model_nominal <- Model
Model_nominal$kappa_pi <- 0
Model_nominal$kappa_y  <- 0
# Correct chi:
Model_nominal$chi <- chi / exp(Model_nominal$kappa_pi*mean_pi + Model_nominal$kappa_y*mean_y)
Model_solved_nominal <- solve_ToyModel(all_d,all_rr,all_eps,proba_eps,
                                       Model_nominal,nb_iter = nb_iter)
strat_nominal <- run_strategy(Model_solved_nominal,maxH=10)


# Issuance of TIPS: ------------------------------------------------------------
Model_TIPS <- Model
Model_TIPS$kappa_pi <- 1
Model_TIPS$kappa_y  <- 0
# Correct chi:
Model_TIPS$chi <- chi / exp(Model_TIPS$kappa_pi*mean_pi + Model_TIPS$kappa_y*mean_y)
Model_solved_TIPS <- solve_ToyModel(all_d,all_rr,all_eps,proba_eps,
                                    Model_TIPS,nb_iter = nb_iter)
strat_TIPS <- run_strategy(Model_solved_TIPS,maxH=10)


# Issuance of GDP-LBs: ---------------------------------------------------------
Model_GDPLB <- Model
Model_GDPLB$kappa_pi <- 1
Model_GDPLB$kappa_y  <- 1
# Correct chi:
Model_GDPLB$chi <- chi / exp(Model_GDPLB$kappa_pi*mean_pi + Model_GDPLB$kappa_y*mean_y)
Model_solved_GDPLB <- solve_ToyModel(all_d,all_rr,all_eps,proba_eps,
                                     Model_GDPLB,nb_iter = nb_iter)
strat_GDPLB <- run_strategy(Model_solved_GDPLB,maxH=10)

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




